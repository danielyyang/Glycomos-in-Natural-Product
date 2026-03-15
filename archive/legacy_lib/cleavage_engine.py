import copy
from rdkit import Chem

def get_branch_atoms(mol, start_idx, avoid_set):
    """获取一个侧链的所有原子索引，避开 avoid_set 中的原子"""
    branch = set()
    queue = [start_idx]
    while queue:
        curr = queue.pop(0)
        if curr in branch or curr in avoid_set:
            continue
        branch.add(curr)
        for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
            if nbr.GetIdx() not in avoid_set and nbr.GetIdx() not in branch:
                queue.append(nbr.GetIdx())
    return branch

def cleave_glycan_aglycan(mol, mapped_sugar_units):
    """
    外科手术级切割：安全分离糖链与苷元，并注入同位素 Dummy Atoms 保护拓扑连接点。
    :param mol: RDKit Mol 对象 (原始完整分子)
    :param mapped_sugar_units: 第一阶段引擎输出的真糖单元列表
    :return: tuple (Glycan_SMILES, Aglycan_SMILES)
    """
    if not mapped_sugar_units:
        # 拓扑铁律：未识别到任何真糖，直接作为假阳性剔除 (Phase 2 规则)
        return None, None

    # 1. 收集所有的 Glycan Zone 原子索引
    core_sugar_atoms = set()
    for u in mapped_sugar_units:
        ring = set(u['ring_atoms'])
        core_sugar_atoms.update(ring)
        # 加入直接相连的碳原子（如 C6 等各种取代碳）
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == 'C' and nbr.GetIdx() not in ring:
                    core_sugar_atoms.add(nbr.GetIdx())
                    
    # 加入直接相连的杂原子（O, N, S, P），确保羟基、氨基以及连接氧全部保留在糖区
    sugar_zone = set(core_sugar_atoms)
    for idx in list(sugar_zone):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetSymbol() in ('O', 'N', 'S', 'P'):
                sugar_zone.add(nbr.GetIdx())

    aglycan_zone = set()
    non_sugar_zone = set(range(mol.GetNumAtoms())) - sugar_zone
    
    # 将剩余分支按特征划入 Aglycan Zone
    for idx in non_sugar_zone:
        if idx in aglycan_zone:
            continue
        branch = get_branch_atoms(mol, idx, avoid_set=sugar_zone)
        
        # Aglycan 判断标准：包含非糖环，或者碳链>=8
        has_ring = False
        carbon_count = 0
        for b_idx in branch:
            b_atom = mol.GetAtomWithIdx(b_idx)
            if b_atom.GetSymbol() == 'C':
                carbon_count += 1
            if b_atom.IsInRing() and b_idx not in core_sugar_atoms:
                has_ring = True
                
        if has_ring or carbon_count >= 8:
            aglycan_zone.update(branch)
        else:
            sugar_zone.update(branch)
            
    if not aglycan_zone:
        return Chem.MolToSmiles(mol), "NULL"

    # 2. 扫描分子，寻找横跨 Glycan 和 Aglycan 的边界键 (Boundary Bonds)
    bonds_to_cut = []
    dummy_labels = [] 
    
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        
        if a1 in sugar_zone and a2 in aglycan_zone:
            bonds_to_cut.append(bond.GetIdx())
            # 留在 Glycan 侧的断点标记为 [14*], Aglycan 侧为 [15*]
            dummy_labels.append((14, 15))
        elif a2 in sugar_zone and a1 in aglycan_zone:
            bonds_to_cut.append(bond.GetIdx())
            dummy_labels.append((15, 14))

    if not bonds_to_cut:
        return Chem.MolToSmiles(mol), "NULL"

    # 3. 执行安全切割 (使用 RDKit 原生碎片化工具)
    fragmented_mol = Chem.FragmentOnBonds(
        mol, 
        bondIndices=bonds_to_cut, 
        dummyLabels=dummy_labels
    )

    # 4. 分离碎片 (Chem.GetMolFrags)
    frags_indices = Chem.GetMolFrags(fragmented_mol)
    glycan_smi_list = []
    aglycan_smi_list = []
    
    true_ring_atoms = set()
    for u in mapped_sugar_units:
        true_ring_atoms.update(u['ring_atoms'])
        
    for indices in frags_indices:
        is_glycan = False
        for idx in indices:
            if idx < mol.GetNumAtoms() and idx in true_ring_atoms:
                is_glycan = True
                break
                
        frag_smi = Chem.MolFragmentToSmiles(fragmented_mol, atomsToUse=list(indices), isomericSmiles=True)
        
        if is_glycan:
            glycan_smi_list.append(frag_smi)
        else:
            aglycan_smi_list.append(frag_smi)

    glycan_final = ".".join(glycan_smi_list) if glycan_smi_list else "NULL"
    aglycan_final = ".".join(aglycan_smi_list) if aglycan_smi_list else "NULL"
    
    return glycan_final, aglycan_final

if __name__ == "__main__":
    import os
    import sys
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    
    # 模拟一次主流程
    rg1_smiles = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(O)C[C@@H](C)O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)C[C@H](O)[C@@H]5[C@@]3(C)CC[C@H](C5(C)C)O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O"
    mol = Chem.MolFromSmiles(rg1_smiles)
    
    from lib.sugar_utils import find_mapped_sugar_units
    
    mapped_units = find_mapped_sugar_units(mol)
    
    print("Testing Ginsenoside Rg1 with Cleavage Engine...")
    glycan, aglycan = cleave_glycan_aglycan(mol, mapped_units)
    print("Glycan Part:", glycan)
    print("Aglycan Part:", aglycan)
