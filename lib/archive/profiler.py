import pubchempy as pcp
from chembl_webresource_client.new_client import new_client
import logging
import time
from typing import Dict, Any, Optional

# 配置日志 (Configure logging)
logger = logging.getLogger(__name__)

class PropertyProfiler:
    """
    负责从 PubChem 和 ChEMBL 获取化合物的物理化学性质与生物活性数据
    (Responsible for retrieving physicochemical properties and bioactivity data from PubChem and ChEMBL)
    """

    def __init__(self):
        self.target_client = new_client.target
        self.activity_client = new_client.activity

    def get_pubchem_props(self, inchikey: str) -> Optional[Dict]:
        """
        从 PubChem 获取详细属性 (Fetch detailed properties from PubChem)
        """
        try:
            # 使用 InChIKey 搜索
            compounds = pcp.get_compounds(inchikey, 'inchikey')
            if not compounds:
                return None
            
            comp = compounds[0]
            
            # 获取所有属性并转换为字典 (Get all properties and convert to dict)
            # pcp.Compound.to_dict() returns many fields.
            full_data = comp.to_dict(properties=[
                'molecular_weight', 'xlogp', 'tpsa', 'complexity', 
                'heavy_atom_count', 'h_bond_donor_count', 'h_bond_acceptor_count',
                'rotatable_bond_count', 'exact_mass', 'monoisotopic_mass',
                'charge', 'isotope_atom_count', 'defined_atom_stereocenter_count',
                'undefined_atom_stereocenter_count', 'defined_bond_stereocenter_count',
                'undefined_bond_stereocenter_count', 'covalent_unit_count',
                'volume_3d', 'xsheric_3d', 'y3d', 'z3d'
            ])
            
            # Clean up keys (remove None values if preferred)
            cleaned_data = {k: v for k, v in full_data.items() if v is not None}
            # Add InChIKey explicit
            cleaned_data['InChIKey'] = inchikey
            
            logger.info(f"PubChem data found for {inchikey}")
            
            # Rate limiting prevention
            time.sleep(1) # 遵守 API 频率限制
            
            return cleaned_data
            
        except Exception as e:
            logger.error(f"Error fetching PubChem props for {inchikey}: {e}")
            return None

    def get_chembl_data(self, inchikey: str) -> Dict[str, Any]:
        """
        从 ChEMBL 获取活性数据 (Retrieve activity data from ChEMBL)
        """
        data = {}
        try:
            # 搜索分子 (Search for molecule similar functionality if necessary, but here exact match)
            # ChEMBL API is complex, simplified example:
            molecule = new_client.molecule
            res = molecule.filter(molecule_structures__standard_inchikey=inchikey).only('molecule_chembl_id')
            
            if res:
                chembl_id = res[0]['molecule_chembl_id']
                data['ChEMBL_ID'] = chembl_id
                logger.info(f"ChEMBL ID found: {chembl_id}")
                
                # Fetch activities (IC50, etc.)
                activities = self.activity_client.filter(molecule_chembl_id=chembl_id, type='IC50').only(['target_chembl_id', 'value', 'units'])
                if activities:
                    data['Activities'] = list(activities)[:5] # Limit to top 5
            else:
                logger.warning(f"No ChEMBL entry for {inchikey}")
                
            time.sleep(1)
            
        except Exception as e:
            logger.error(f"Error fetching ChEMBL data: {e}")
            
        return data

if __name__ == "__main__":
    profiler = PropertyProfiler()
    # Test with Solamargine InChIKey
    test_key = "KUXHORMIBTWVTM-SNAWJCMRSA-N"
    print(profiler.get_pubchem_props(test_key))
