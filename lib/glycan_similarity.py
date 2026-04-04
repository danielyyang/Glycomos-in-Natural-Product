#!/usr/bin/env python3
"""
==========================================================================
  [EN] Glycan Tree Edit Distance (TED) — Glycan Sequence Similarity Engine
       A graph-similarity algorithm that treats sugar sequences as rooted,
       ordered trees and computes a weighted edit distance across three
       scoring dimensions: node (monosaccharide), edge (linkage), and
       topology (branching pattern).

  [CN] 糖链树编辑距离 (TED) — 糖序列相似性匹配引擎
       将糖链序列视为有根有序树, 在三个打分维度 (节点/边/拓扑) 上
       计算加权编辑距离, 用于在数据库中检索"最相似的糖链"。

  核心算法: Zhang-Shasha Tree Edit Distance (Zhang & Shasha, 1989)
  时间复杂度: O(n1 * n2 * min(d1, l1) * min(d2, l2))
  其中 n = 节点数, d = 深度, l = 叶子数

  Reference: K. Zhang and D. Shasha, "Simple Fast Algorithms for the
  Editing Distance between Trees and Related Problems," SIAM J. Comput.,
  vol. 18, no. 6, pp. 1245-1262, 1989.

  [TEST DATA ONLY]
==========================================================================
"""
import re
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple

# =========================================================================
# 一、糖链树节点数据结构 (Glycan Tree Node Data Structure)
#
# 设计意图: 每个节点表示一个单糖残基, 边表示糖苷键连接。
# 节点携带单糖名称 (如 D-Glc)、修饰标签。
# 边携带异头碳构型 (α/β) 和连接位置 (如 1→4)。
#
# Design: Each node = one monosaccharide residue; edges = glycosidic bonds.
# Node stores sugar name + modifications; edge stores anomer + linkage position.
# =========================================================================


@dataclass
class GlycanNode:
    """树中的一个单糖节点 (One monosaccharide node in the tree).

    Attributes:
        sugarName: 单糖名称 (如 "D-Glc", "L-Rha")
        anomerConfig: 与父节点连接的异头碳构型 ("α" / "β" / "?")
        linkagePosition: 与父节点连接的位置 (如 "1→4" / "1→2")
        modifications: 修饰标签列表 (如 ["O-Ac", "NAc"])
        children: 子节点列表 (Child nodes)
    """
    sugarName: str = ""
    anomerConfig: str = "?"
    linkagePosition: str = ""
    modifications: List[str] = field(default_factory=list)
    children: List["GlycanNode"] = field(default_factory=list)

    def nodeCount(self) -> int:
        """递归统计树的节点总数 (Recursively count total nodes)."""
        return 1 + sum(c.nodeCount() for c in self.children)

    def depth(self) -> int:
        """计算树的最大深度 (Compute maximum tree depth)."""
        if not self.children:
            return 1
        return 1 + max(c.depth() for c in self.children)

    def leafCount(self) -> int:
        """统计叶子节点数 (Count leaf nodes)."""
        if not self.children:
            return 1
        return sum(c.leafCount() for c in self.children)

    def toSequenceString(self) -> str:
        """将树还原为序列字符串 (Reconstruct sequence string from tree)."""
        if not self.children:
            return self.sugarName
        parts = []
        for child in self.children:
            childStr = child.toSequenceString()
            linkage = f"({child.anomerConfig}{child.linkagePosition})" if child.linkagePosition else ""
            parts.append(f"{childStr}-{linkage}-" if linkage else f"{childStr}-")
        return "".join(parts) + self.sugarName

    def getAllNodes(self) -> List["GlycanNode"]:
        """后序遍历获取所有节点 (Post-order traversal to get all nodes)."""
        result = []
        for child in self.children:
            result.extend(child.getAllNodes())
        result.append(self)
        return result


# =========================================================================
# 二、糖链序列解析器 (Glycan Sequence Parser)
#
# 设计意图: 将 GlycoNP 格式的糖链序列字符串解析为树结构。
# 支持线性链和分支结构 (用 [...] 表示分支)。
#
# Design: Parses GlycoNP-format sugar sequence strings into tree structure.
# Supports both linear chains and branched structures ([...] = branch).
# =========================================================================

def parseGlycanSequence(sequence: str) -> Optional[GlycanNode]:
    """将 GlycoNP 格式的糖链序列解析为树结构。
    Parse GlycoNP-format glycan sequence into a rooted tree.

    支持的格式 (Supported formats):
        线性 (Linear): "D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA"
        分支 (Branched): "[D-Qui-(b1-2)]-D-Gal-(b1-4)-D-Xyl"
        多链 (Multi-chain): "L-Rha ; L-Ara-(a1-2)-D-Glc"

    Args:
        sequence: GlycoNP 格式的糖链序列字符串

    Returns:
        GlycanNode: 解析后的树根节点, 或 None (解析失败时)
    """
    if not sequence or not isinstance(sequence, str):
        return None

    sequence = sequence.strip()
    if not sequence:
        return None

    # 处理多链 (Handle multi-chain with ";")
    # 多链皂苷中, 每条链独立连接到苷元 → 创建虚拟根节点
    # Multi-chain: each chain connects independently → create virtual root
    chains = [c.strip() for c in sequence.split(";") if c.strip()]
    if len(chains) > 1:
        virtualRoot = GlycanNode(sugarName="[ROOT]")
        for chain in chains:
            subtree = _parseSingleChain(chain)
            if subtree:
                virtualRoot.children.append(subtree)
        return virtualRoot if virtualRoot.children else None

    return _parseSingleChain(chains[0])


def _parseSingleChain(chain: str) -> Optional[GlycanNode]:
    """解析单条糖链 (可含分支) 为树。
    Parse a single chain (may contain branches) into a tree.

    解析策略 (Parsing strategy):
    1. 从右到左读取 (根糖在最右边)
    2. 遇到 [...] 则递归解析为分支子树
    3. 遇到 (a1-4) 则记录为连接信息

    Args:
        chain: 单条糖链字符串 (如 "[D-Glc-(b1-3)]-D-Gal-(b1-4)-D-Xyl")

    Returns:
        GlycanNode: 树根节点
    """
    if not chain:
        return None

    # 提取分支: 找到所有 [...] 块 (Extract branches: find all [...] blocks)
    branches = []
    cleanChain = chain
    branchPattern = re.compile(r"\[([^\[\]]+)\]")
    while branchPattern.search(cleanChain):
        match = branchPattern.search(cleanChain)
        branches.append(match.group(1))
        cleanChain = cleanChain[:match.start()] + cleanChain[match.end():]

    # 处理嵌套分支 [[...]] — 展平 (Flatten nested branches)
    nestedPattern = re.compile(r"\[\[([^\]]+)\]\]")
    while nestedPattern.search(cleanChain):
        match = nestedPattern.search(cleanChain)
        branches.append(match.group(1))
        cleanChain = cleanChain[:match.start()] + cleanChain[match.end():]

    # 清理残留的连接符 (Clean up residual separators)
    cleanChain = re.sub(r"^[\s\-]+|[\s\-]+$", "", cleanChain)
    cleanChain = re.sub(r"\-\-+", "-", cleanChain)

    # 解析主链 (Parse main chain)
    mainNode = _parseLinearChain(cleanChain)
    if mainNode is None:
        return None

    # 将分支附加到主链的适当位置 (Attach branches to main chain)
    # 简化策略: 分支总是附加到主链从根数第一个糖上
    # Simplified: branches attach to the root of the main chain
    for branchStr in branches:
        branchNode = _parseLinearChain(branchStr)
        if branchNode:
            mainNode.children.append(branchNode)

    return mainNode


def _parseLinearChain(chain: str) -> Optional[GlycanNode]:
    """解析线性 (无分支) 糖链为链表结构的树。
    Parse a linear chain into a linked-list shaped tree.

    格式: "D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA"
    根节点 = 最右边的糖 (D-GlcA); 最左边 = 叶子 (非还原端)

    Format: "D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA"
    Root = rightmost sugar (D-GlcA); leftmost = leaf (non-reducing end)
    """
    if not chain:
        return None

    # 分割: sugar-(link)-sugar-(link)-sugar
    # 模式: 糖名 + 可选的连接
    tokens = re.split(r"-(\([ab]\d-\d\))-", chain)
    # tokens 结构: [sugar1, link1, sugar2, link2, sugar3, ...]
    # 当只有一个糖时: [sugar1]

    if not tokens:
        return None

    # 从右(根)到左(叶)构建 (Build from right/root to left/leaf)
    sugarNames = tokens[0::2]  # 索引 0, 2, 4, ...
    linkages = tokens[1::2]     # 索引 1, 3, 5, ...

    # 清理糖名 (Clean sugar names)
    sugarNames = [s.strip().strip("-") for s in sugarNames if s.strip().strip("-")]

    if not sugarNames:
        return None

    # 根节点 = 最右边 (还原端连接苷元的糖)
    # Root = rightmost (reducing-end sugar connecting to aglycone)
    root = GlycanNode(sugarName=sugarNames[-1])
    current = root

    # 从倒数第二个糖开始向左构建 (Build leftward from second-to-last)
    for i in range(len(sugarNames) - 2, -1, -1):
        linkIdx = i  # linkages[i] 连接 sugarNames[i] 和 sugarNames[i+1]
        anomer = "?"
        position = ""
        if linkIdx < len(linkages):
            linkStr = linkages[linkIdx]  # 如 "(b1-4)"
            linkMatch = re.match(r"\(([ab])(\d)-(\d)\)", linkStr)
            if linkMatch:
                anomer = "α" if linkMatch.group(1) == "a" else "β"
                position = f"{linkMatch.group(2)}→{linkMatch.group(3)}"

        childNode = GlycanNode(
            sugarName=sugarNames[i],
            anomerConfig=anomer,
            linkagePosition=position,
        )
        current.children.insert(0, childNode)
        current = childNode

    return root


# =========================================================================
# 三、打分矩阵 (Scoring Matrices)
#
# 设计意图: 三维打分系统, 分别评估节点/边/拓扑的差异程度。
# 权重可调, 允许用户根据研究目的侧重不同维度。
#
# Design: Three-dimensional scoring system evaluating differences across
# node/edge/topology dimensions. Weights are adjustable for different
# research objectives.
# =========================================================================

# 单糖差向异构体组 (Monosaccharide epimer groups)
# 设计意图: 同一组内的糖仅在某个碳位手性不同, 结构相似度高
# Design: Sugars in same group differ only at one chiral center → high similarity
EPIMER_GROUPS: Dict[str, int] = {
    # C4 差向异构体 (C4 epimers)
    "D-Glc": 1, "D-Gal": 1,  # Glc ↔ Gal (C4 epimerization)
    "D-Man": 2,                # Man 独立组 (独立 C2 差向异构体)
    "D-GlcNAc": 3, "D-GalNAc": 3,  # GlcNAc ↔ GalNAc (C4)
    "D-GlcA": 4, "D-GalA": 4,      # GlcA ↔ GalA (C4)
    "L-Rha": 5, "L-Fuc": 5,  # Rha ↔ Fuc (6-deoxy 系列 C4)
    "D-Xyl": 6, "L-Ara": 6,  # Xyl ↔ Ara (戊糖 C4)
}

# 氧化态关系 (Oxidation state relationships)
# 设计意图: Glc → GlcA 是 C6 氧化, 结构变化较大但仍有骨架关联
# Design: Glc → GlcA is C6 oxidation; significant but related
OXIDATION_PAIRS = {
    frozenset({"D-Glc", "D-GlcA"}),
    frozenset({"D-Gal", "D-GalA"}),
    frozenset({"D-Man", "D-ManA"}),
}


@dataclass
class ScoringWeights:
    """三维打分权重 (Three-dimensional scoring weights).

    Attributes:
        nodeWeight: 节点 (单糖) 匹配权重 — 影响"原料替换代价"评估
        edgeWeight: 边 (连接键) 匹配权重 — 影响"构象相似性"评估
        topologyWeight: 拓扑 (分支) 匹配权重 — 影响"空间排布"评估
    """
    nodeWeight: float = 1.0
    edgeWeight: float = 1.0
    topologyWeight: float = 0.5


def computeNodeCost(nodeA: GlycanNode, nodeB: GlycanNode) -> float:
    """计算两个单糖节点之间的替换代价。
    Compute substitution cost between two monosaccharide nodes.

    打分规则 (Scoring rules):
        完全匹配 (Exact match): 0.0
        差向异构体 (C4 epimer): 0.3 — 结构非常相似
        氧化态变化 (Oxidation): 0.6 — 骨架关联但化学性质改变
        同为该泛指糖 (Same generic class): 0.5
        完全不同 (Different): 1.0

    Args:
        nodeA: 第一个节点
        nodeB: 第二个节点

    Returns:
        float: 替换代价 [0.0, 1.0]
    """
    nameA = nodeA.sugarName
    nameB = nodeB.sugarName

    # 虚拟根节点匹配 (Virtual root match)
    if nameA == "[ROOT]" or nameB == "[ROOT]":
        return 0.0 if nameA == nameB else 0.5

    # 完全匹配 (Exact match)
    if nameA == nameB:
        return 0.0

    # 差向异构体 (Epimer check — same group)
    groupA = EPIMER_GROUPS.get(nameA)
    groupB = EPIMER_GROUPS.get(nameB)
    if groupA is not None and groupA == groupB:
        return 0.3

    # 氧化态关系 (Oxidation relationship)
    if frozenset({nameA, nameB}) in OXIDATION_PAIRS:
        return 0.6

    # D/L 互换 (D/L enantiomer swap)
    # 设计意图: D-Glc ↔ L-Glc 是镜像异构, 空间完全不同
    # Design: D-Glc ↔ L-Glc are mirror images, spatially very different
    baseA = re.sub(r"^[DL]-", "", nameA)
    baseB = re.sub(r"^[DL]-", "", nameB)
    if baseA == baseB:
        return 0.4  # 同骨架不同手性 (Same backbone, different chirality)

    # 泛指糖匹配 (Generic sugar matching)
    # Hex 匹配任何六碳糖, Pen 匹配任何五碳糖
    GENERIC_MAP = {
        "Hex": {"D-Glc", "D-Gal", "D-Man", "L-Glc", "L-Gal", "L-Man"},
        "HexA": {"D-GlcA", "D-GalA", "D-ManA", "L-GlcA", "L-GalA"},
        "Pen": {"D-Xyl", "L-Ara", "D-Rib", "D-Lyx", "L-Lyx"},
        "dHex": {"L-Rha", "D-Rha", "L-Fuc", "D-Fuc", "D-Qui", "L-Qui"},
    }
    for generic, members in GENERIC_MAP.items():
        if (nameA == generic and nameB in members) or (nameB == generic and nameA in members):
            return 0.5

    # 完全不同 (Completely different)
    return 1.0


def computeEdgeCost(nodeA: GlycanNode, nodeB: GlycanNode) -> float:
    """计算两条糖苷键 (边) 之间的差异代价。
    Compute difference cost between two glycosidic bonds (edges).

    打分规则 (Scoring rules):
        完全匹配: 0.0 (如 β(1→4) == β(1→4))
        仅位置不同 (区域异构): 0.6 (如 1→4 变为 1→3)
        仅异头碳不同: 0.8 (如 α 变为 β — 空间结构完全翻转)
        位置和异头碳都不同: 1.0

    Args:
        nodeA: 第一个节点 (其连接信息 = 与父节点的边)
        nodeB: 第二个节点

    Returns:
        float: 边差异代价 [0.0, 1.0]
    """
    # 如果两个节点都没有连接信息 (都是根节点), 边代价为 0
    # If both nodes lack linkage info (both are roots), edge cost is 0
    if not nodeA.linkagePosition and not nodeB.linkagePosition:
        if nodeA.anomerConfig == nodeB.anomerConfig or nodeA.anomerConfig == "?" or nodeB.anomerConfig == "?":
            return 0.0

    # 如果任一边缺失连接信息, 给中等惩罚
    # If either edge has missing linkage info, give moderate penalty
    if not nodeA.linkagePosition or not nodeB.linkagePosition:
        if nodeA.anomerConfig == nodeB.anomerConfig:
            return 0.2
        return 0.5

    anomerMatch = (nodeA.anomerConfig == nodeB.anomerConfig)
    posMatch = (nodeA.linkagePosition == nodeB.linkagePosition)

    if anomerMatch and posMatch:
        return 0.0
    elif anomerMatch and not posMatch:
        return 0.6  # 区域异构 (Regioisomer — position changed)
    elif not anomerMatch and posMatch:
        return 0.8  # 异头碳翻转 (Anomeric flip — major structural change)
    else:
        return 1.0  # 完全不同 (Completely different)


# =========================================================================
# 四、Zhang-Shasha 树编辑距离算法 (Zhang-Shasha TED Algorithm)
#
# 设计意图: 经典的树编辑距离算法, 支持三种操作:
#   1. 删除 (Delete): 移除一个节点, 其子节点上移
#   2. 插入 (Insert): 在某位置插入一个新节点
#   3. 替换 (Substitute): 将一个节点替换为另一个节点
# 每种操作都分配一个由打分矩阵决定的代价。
#
# Design: Classic tree edit distance algorithm supporting:
#   1. Delete: remove a node, children move up
#   2. Insert: add a new node at some position
#   3. Substitute: replace one node with another
# Each operation has costs determined by the scoring matrices.
# =========================================================================

def _postOrderIndex(root: GlycanNode) -> Tuple[List[GlycanNode], Dict[int, int]]:
    """后序遍历并返回节点列表 + 最左叶子索引映射。
    Post-order traversal returning node list + leftmost leaf index mapping.

    Zhang-Shasha 算法需要后序遍历编号和每个节点的"最左叶子后代"索引。
    The ZS algorithm requires post-order numbering and each node's
    leftmost leaf descendant index.

    Returns:
        nodes: 后序排列的节点列表 (Post-order indexed node list)
        leftmost: 每个后序索引 → 其最左叶子的后序索引
    """
    nodes: List[GlycanNode] = []
    leftmost: Dict[int, int] = {}

    def _traverse(node: GlycanNode) -> int:
        """递归后序遍历, 返回该子树最左叶子的索引 (Recursive post-order)."""
        myLeftmost = -1
        for child in node.children:
            childLeftmost = _traverse(child)
            if myLeftmost == -1:
                myLeftmost = childLeftmost
        idx = len(nodes)
        nodes.append(node)
        if myLeftmost == -1:
            myLeftmost = idx  # 叶子节点的最左叶子就是自己 (Leaf's leftmost is itself)
        leftmost[idx] = myLeftmost
        return myLeftmost

    _traverse(root)
    return nodes, leftmost


def _computeKeyRoots(leftmost: Dict[int, int], nodeCount: int) -> List[int]:
    """计算关键根集合 (Compute key-root set for Zhang-Shasha).

    关键根 = 最左叶子索引唯一的节点。这是 ZS 算法的核心优化:
    只需在关键根对之间计算子问题。
    Key roots = nodes with unique leftmost leaf indices. This is the
    core optimization: only compute sub-problems between key root pairs.
    """
    seen = {}
    for i in range(nodeCount):
        seen[leftmost[i]] = i
    return sorted(seen.values())


def computeTreeEditDistance(
    treeA: GlycanNode,
    treeB: GlycanNode,
    weights: Optional[ScoringWeights] = None,
) -> float:
    """计算两棵糖链树之间的加权编辑距离。
    Compute weighted edit distance between two glycan trees.

    实现 Zhang-Shasha 1989 算法, 使用自定义的三维代价函数。
    Implements Zhang-Shasha 1989 with custom three-dimensional cost functions.

    Args:
        treeA: 第一棵树 (查询树)
        treeB: 第二棵树 (候选树)
        weights: 打分权重 (Scoring weights, defaults to equal)

    Returns:
        float: 加权编辑距离 (Weighted edit distance, lower = more similar)
    """
    if weights is None:
        weights = ScoringWeights()

    nodesA, leftmostA = _postOrderIndex(treeA)
    nodesB, leftmostB = _postOrderIndex(treeB)
    sizeA = len(nodesA)
    sizeB = len(nodesB)

    keyRootsA = _computeKeyRoots(leftmostA, sizeA)
    keyRootsB = _computeKeyRoots(leftmostB, sizeB)

    # 删除/插入代价 (Deletion/insertion costs)
    # 设计意图: 删除/插入一个糖 = 合成路线增加/减少一步
    # Design: Delete/insert a sugar = one more/fewer synthetic step
    DELETE_COST = 1.0
    INSERT_COST = 1.0

    def substitutionCost(idxA: int, idxB: int) -> float:
        """加权替换代价 (Weighted substitution cost)."""
        nA, nB = nodesA[idxA], nodesB[idxB]
        nodeCost = computeNodeCost(nA, nB)
        edgeCost = computeEdgeCost(nA, nB)
        return (weights.nodeWeight * nodeCost +
                weights.edgeWeight * edgeCost)

    # 初始化全局距离矩阵 (Initialize global distance matrix)
    td = [[0.0] * (sizeB + 1) for _ in range(sizeA + 1)]

    # Zhang-Shasha 核心循环 (Core ZS loop)
    for krA in keyRootsA:
        for krB in keyRootsB:
            lA = leftmostA[krA]
            lB = leftmostB[krB]

            # 临时子树距离矩阵 (Temporary forest distance matrix)
            m = krA - lA + 2
            n = krB - lB + 2
            fd = [[0.0] * n for _ in range(m)]

            # 初始化边界 (Initialize boundaries)
            fd[0][0] = 0.0
            for i in range(1, m):
                fd[i][0] = fd[i - 1][0] + DELETE_COST
            for j in range(1, n):
                fd[0][j] = fd[0][j - 1] + INSERT_COST

            for i in range(1, m):
                for j in range(1, n):
                    idxA = lA + i - 1
                    idxB = lB + j - 1

                    costDel = fd[i - 1][j] + DELETE_COST
                    costIns = fd[i][j - 1] + INSERT_COST

                    if leftmostA[idxA] == lA and leftmostB[idxB] == lB:
                        # 两个子树的根: 直接替换 (Both are subtree roots: direct substitution)
                        costSub = fd[i - 1][j - 1] + substitutionCost(idxA, idxB)
                        fd[i][j] = min(costDel, costIns, costSub)
                        td[idxA + 1][idxB + 1] = fd[i][j]
                    else:
                        # 使用已计算的子问题 (Use precomputed sub-problem)
                        p = leftmostA[idxA] - lA
                        q = leftmostB[idxB] - lB
                        costSub = fd[p][q] + td[idxA + 1][idxB + 1]
                        fd[i][j] = min(costDel, costIns, costSub)

    return td[sizeA][sizeB]


# =========================================================================
# 五、归一化相似度 (Normalized Similarity Score)
#
# 设计意图: 原始编辑距离与树大小成正比, 不同大小的树之间不可比。
# 归一化到 [0, 1] 区间, 1 = 完全相同, 0 = 完全不同。
#
# Design: Raw edit distance scales with tree size and isn't comparable
# across different-sized trees. Normalize to [0, 1] where 1 = identical.
# =========================================================================

def computeSimilarity(
    treeA: GlycanNode,
    treeB: GlycanNode,
    weights: Optional[ScoringWeights] = None,
) -> float:
    """计算两棵糖链树的归一化相似度。
    Compute normalized similarity between two glycan trees.

    归一化公式: similarity = 1 - (distance / maxPossibleDistance)
    maxPossibleDistance = max(sizeA, sizeB) * maxCostPerOperation

    Args:
        treeA: 第一棵树
        treeB: 第二棵树
        weights: 打分权重

    Returns:
        float: 归一化相似度 [0.0, 1.0], 1.0 = 完全相同
    """
    if weights is None:
        weights = ScoringWeights()

    sizeA = treeA.nodeCount()
    sizeB = treeB.nodeCount()

    if sizeA == 0 and sizeB == 0:
        return 1.0
    if sizeA == 0 or sizeB == 0:
        return 0.0

    distance = computeTreeEditDistance(treeA, treeB, weights)

    # 最大可能距离: 删除所有 A 的节点 + 插入所有 B 的节点
    # Max possible distance: delete all of A + insert all of B
    maxDistance = float(sizeA + sizeB)

    similarity = 1.0 - (distance / maxDistance) if maxDistance > 0 else 1.0
    return max(0.0, min(1.0, similarity))


# =========================================================================
# 六、批量相似性检索 (Batch Similarity Search)
#
# 设计意图: 给定一条查询糖链, 在整个数据库中找出最相似的 Top-N。
# 这是"糖序列相似性匹配系统"的核心功能入口。
#
# Design: Given a query glycan, find Top-N most similar sequences in
# the entire database. Core entry point for the similarity matching system.
# =========================================================================

def searchSimilar(
    querySequence: str,
    candidateSequences: List[str],
    topN: int = 10,
    weights: Optional[ScoringWeights] = None,
) -> List[Dict]:
    """在候选序列中检索与查询序列最相似的 Top-N。
    Search for Top-N most similar sequences to the query.

    Args:
        querySequence: 查询糖链序列 (如 "D-Glc-(b1-4)-L-Rha")
        candidateSequences: 候选序列列表
        topN: 返回的最相似结果数
        weights: 打分权重

    Returns:
        List[Dict]: 按相似度降序排列的结果列表, 每项包含:
            - sequence: 候选序列字符串
            - similarity: 相似度分数 [0, 1]
            - distance: 原始编辑距离
    """
    queryTree = parseGlycanSequence(querySequence)
    if queryTree is None:
        return []

    results = []
    for i, candidateSeq in enumerate(candidateSequences):
        if not candidateSeq or not isinstance(candidateSeq, str):
            continue
        candidateTree = parseGlycanSequence(candidateSeq)
        if candidateTree is None:
            continue

        distance = computeTreeEditDistance(queryTree, candidateTree, weights)
        similarity = computeSimilarity(queryTree, candidateTree, weights)
        results.append({
            "index": i,
            "sequence": candidateSeq,
            "similarity": round(similarity, 4),
            "distance": round(distance, 4),
        })

    # 按相似度降序排列 (Sort by similarity descending)
    results.sort(key=lambda x: -x["similarity"])
    return results[:topN]
