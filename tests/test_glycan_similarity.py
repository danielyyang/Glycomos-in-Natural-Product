"""
test_glycan_similarity.py
=========================
[EN] Pytest suite for the Glycan Tree Edit Distance (TED) module.
     Tests the sequence parser, scoring matrices, TED algorithm,
     and batch similarity search across diverse glycan structures.

[CN] 糖链树编辑距离 (TED) 模块的 Pytest 测试套件。
     覆盖序列解析器、打分矩阵、TED 算法和批量相似性检索。

[TEST DATA ONLY]
"""
import pytest
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.glycan_similarity import (
    GlycanNode,
    parseGlycanSequence,
    computeNodeCost,
    computeEdgeCost,
    computeTreeEditDistance,
    computeSimilarity,
    searchSimilar,
    ScoringWeights,
)


# =========================================================================
# 序列解析器测试 (Sequence Parser Tests)
# =========================================================================

class TestGlycanParser:
    """测试糖链序列解析器的各种输入格式 (Test parser for various input formats)."""

    def test_parseSingleSugar(self):
        """单个糖: 应解析为单节点树 (Single sugar → single-node tree)."""
        tree = parseGlycanSequence("D-Glc")
        assert tree is not None
        assert tree.sugarName == "D-Glc"
        assert tree.nodeCount() == 1
        assert tree.children == []

    def test_parseLinearDisaccharide(self):
        """线性二糖: 根 = 还原端, 叶 = 非还原端
        Linear disaccharide: root = reducing end, leaf = non-reducing end."""
        tree = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        assert tree is not None
        # 根节点 = L-Rha (还原端, 连接苷元)
        assert tree.sugarName == "L-Rha"
        assert tree.nodeCount() == 2
        assert len(tree.children) == 1
        # 子节点 = D-Glc (非还原端)
        child = tree.children[0]
        assert child.sugarName == "D-Glc"
        assert child.anomerConfig == "β"
        assert child.linkagePosition == "1→4"

    def test_parseLinearTrisaccharide(self):
        """线性三糖链 (Linear trisaccharide)."""
        tree = parseGlycanSequence("L-Rha-(a1-2)-D-GlcA-(b1-3)-D-Glc")
        assert tree is not None
        assert tree.sugarName == "D-Glc"  # 根 = 还原端
        assert tree.nodeCount() == 3
        assert tree.depth() == 3

    def test_parseBranchedSequence(self):
        """分支糖链: [...] 表示分支 (Branched: [...] notation)."""
        tree = parseGlycanSequence("[D-Qui-(b1-2)]-D-Gal-(b1-4)-D-Xyl")
        assert tree is not None
        assert tree.sugarName == "D-Xyl"  # 根
        assert tree.nodeCount() >= 3  # 至少 3 个节点

    def test_parseMultiChain(self):
        """多链 (用分号分隔): 创建虚拟根节点
        Multi-chain (semicolon-separated): creates virtual root."""
        tree = parseGlycanSequence("L-Rha ; D-Glc-(a1-4)-D-Gal")
        assert tree is not None
        assert tree.sugarName == "[ROOT]"  # 虚拟根
        assert len(tree.children) == 2

    def test_parseEmptyAndInvalid(self):
        """空字符串和无效输入 → None (Empty/invalid → None)."""
        assert parseGlycanSequence("") is None
        assert parseGlycanSequence(None) is None
        assert parseGlycanSequence("   ") is None

    def test_parseGenericSugar(self):
        """泛指糖 (如 Hex, dHex) 也应可解析 (Generic sugars should parse)."""
        tree = parseGlycanSequence("Hex-(b1-4)-dHex")
        assert tree is not None
        assert tree.nodeCount() == 2


# =========================================================================
# 打分矩阵测试 (Scoring Matrix Tests)
# =========================================================================

class TestScoringMatrices:
    """测试节点和边的打分函数 (Test node and edge cost functions)."""

    def test_exactNodeMatch(self):
        """完全相同的糖 → 代价 0 (Exact match → cost 0)."""
        nodeA = GlycanNode(sugarName="D-Glc")
        nodeB = GlycanNode(sugarName="D-Glc")
        assert computeNodeCost(nodeA, nodeB) == 0.0

    def test_epimerNodeCost(self):
        """差向异构体 → 低代价 (Epimers → low cost)."""
        nodeA = GlycanNode(sugarName="D-Glc")
        nodeB = GlycanNode(sugarName="D-Gal")
        cost = computeNodeCost(nodeA, nodeB)
        assert 0 < cost < 0.5  # 应该是低代价

    def test_oxidationNodeCost(self):
        """氧化态变化 → 中等代价 (Oxidation change → moderate cost)."""
        nodeA = GlycanNode(sugarName="D-Glc")
        nodeB = GlycanNode(sugarName="D-GlcA")
        cost = computeNodeCost(nodeA, nodeB)
        assert 0.3 < cost < 0.8

    def test_completelyDifferentNodeCost(self):
        """完全不同的糖 → 高代价 (Completely different → high cost)."""
        nodeA = GlycanNode(sugarName="D-Glc")
        nodeB = GlycanNode(sugarName="L-Fuc")
        cost = computeNodeCost(nodeA, nodeB)
        assert cost >= 0.8

    def test_exactEdgeMatch(self):
        """完全相同的连接 → 0 (Exact edge match → 0)."""
        nodeA = GlycanNode(anomerConfig="β", linkagePosition="1→4")
        nodeB = GlycanNode(anomerConfig="β", linkagePosition="1→4")
        assert computeEdgeCost(nodeA, nodeB) == 0.0

    def test_positionIsomerEdgeCost(self):
        """仅位置不同 → 中等代价 (Position isomer → moderate)."""
        nodeA = GlycanNode(anomerConfig="β", linkagePosition="1→4")
        nodeB = GlycanNode(anomerConfig="β", linkagePosition="1→3")
        cost = computeEdgeCost(nodeA, nodeB)
        assert 0.3 < cost < 0.8

    def test_anomerFlipEdgeCost(self):
        """异头碳翻转 → 高代价 (Anomeric flip → high cost)."""
        nodeA = GlycanNode(anomerConfig="α", linkagePosition="1→4")
        nodeB = GlycanNode(anomerConfig="β", linkagePosition="1→4")
        cost = computeEdgeCost(nodeA, nodeB)
        assert cost >= 0.7


# =========================================================================
# 树编辑距离算法测试 (TED Algorithm Tests)
# =========================================================================

class TestTreeEditDistance:
    """测试 Zhang-Shasha 树编辑距离算法 (Test ZS TED algorithm)."""

    def test_identicalTrees(self):
        """完全相同的树 → 距离 0 (Identical trees → distance 0)."""
        treeA = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        treeB = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        distance = computeTreeEditDistance(treeA, treeB)
        assert distance == 0.0

    def test_singleDeletion(self):
        """删除一个糖 → 距离 > 0 (Delete one sugar → distance > 0)."""
        treeA = parseGlycanSequence("D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA")
        treeB = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        distance = computeTreeEditDistance(treeA, treeB)
        assert distance > 0

    def test_singleSubstitution(self):
        """替换一个糖 → 距离反映代价 (Substitute one sugar → distance reflects cost)."""
        treeA = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        treeB = parseGlycanSequence("D-Gal-(b1-4)-L-Rha")  # Glc → Gal (差向异构体)
        distance = computeTreeEditDistance(treeA, treeB)
        # 应该是低距离 (差向异构体代价低)
        assert 0 < distance < 1.0

    def test_similarityIdentical(self):
        """相同序列 → 相似度 1.0 (Identical → similarity 1.0)."""
        treeA = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        treeB = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        sim = computeSimilarity(treeA, treeB)
        assert sim == 1.0

    def test_similarityRange(self):
        """相似度应在 [0, 1] 范围内 (Similarity must be in [0, 1])."""
        treeA = parseGlycanSequence("D-Glc-(b1-4)-D-Gal-(a1-6)-L-Rha")
        treeB = parseGlycanSequence("L-Fuc")
        sim = computeSimilarity(treeA, treeB)
        assert 0.0 <= sim <= 1.0

    def test_symmetry(self):
        """距离应该是对称的: d(A,B) == d(B,A) (Distance is symmetric)."""
        treeA = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        treeB = parseGlycanSequence("D-Gal-(a1-2)-D-GlcA")
        distAB = computeTreeEditDistance(treeA, treeB)
        distBA = computeTreeEditDistance(treeB, treeA)
        assert abs(distAB - distBA) < 1e-6

    def test_weightedScoring(self):
        """自定义权重应改变距离 (Custom weights should change distance)."""
        treeA = parseGlycanSequence("D-Glc-(b1-4)-L-Rha")
        treeB = parseGlycanSequence("D-Glc-(a1-4)-L-Rha")  # 仅 β→α 差异

        # 默认权重 (Default weights)
        distDefault = computeTreeEditDistance(treeA, treeB, ScoringWeights())
        # 高边权重 (High edge weight — anomeric changes matter more)
        distHighEdge = computeTreeEditDistance(treeA, treeB, ScoringWeights(edgeWeight=3.0))
        # 低边权重 (Low edge weight — anomeric changes matter less)
        distLowEdge = computeTreeEditDistance(treeA, treeB, ScoringWeights(edgeWeight=0.1))

        assert distHighEdge > distDefault > distLowEdge


# =========================================================================
# 批量检索测试 (Batch Search Tests)
# =========================================================================

class TestBatchSearch:
    """测试批量相似性检索 (Test batch similarity search)."""

    def test_selfIsTopResult(self):
        """查询自己应排第一 (Query itself should rank #1)."""
        candidates = [
            "D-Glc-(b1-4)-L-Rha",
            "D-Gal-(a1-2)-D-GlcA",
            "L-Fuc",
        ]
        results = searchSimilar("D-Glc-(b1-4)-L-Rha", candidates, topN=3)
        assert len(results) > 0
        assert results[0]["similarity"] == 1.0
        assert results[0]["sequence"] == "D-Glc-(b1-4)-L-Rha"

    def test_epimerRanksHigher(self):
        """差向异构体应比完全不同的糖排名更高
        Epimer should rank higher than completely different sugar."""
        candidates = [
            "D-Gal-(b1-4)-L-Rha",   # Glc→Gal 差向异构体 (epimer)
            "L-Fuc-(a1-2)-D-Man",    # 完全不同 (very different)
        ]
        results = searchSimilar("D-Glc-(b1-4)-L-Rha", candidates, topN=2)
        assert len(results) == 2
        # 差向异构体应排在前面 (Epimer should rank first)
        assert "Gal" in results[0]["sequence"]

    def test_emptyQuery(self):
        """空查询 → 空结果 (Empty query → empty results)."""
        results = searchSimilar("", ["D-Glc"], topN=5)
        assert results == []

    def test_topNRespected(self):
        """topN 限制应被遵守 (topN limit should be respected)."""
        candidates = [f"D-Glc" for _ in range(20)]
        results = searchSimilar("D-Glc", candidates, topN=5)
        assert len(results) <= 5
