# GlycoNP-Pipeline

## 这个项目在做什么？ (What Does This Project Do?)

**GlycoNP-Pipeline** 是一个从天然产物数据库 (COCONUT) 中系统性挖掘**糖缀合物** (Glycoconjugates) 规律的计算化学管线。

在天然产物化学中，糖链不是装饰品 —— 它决定了化合物的水溶性、细胞识别、生物利用度甚至药理活性。但传统的天然产物研究几乎都聚焦在苷元骨架上，**糖链被当作噪声丢弃了**。

我们的目标是：把这些被忽视的糖链"翻出来"，看看大自然在不同物种、不同骨架、不同生物功能的化合物上，**到底选择了哪些糖、做了哪些修饰、为什么**。

## What question are we answering?

> **"在 94,242 个天然含糖化合物中，特定的糖链序列和修饰模式，是否与物种来源、苷元骨架类型存在系统性的关联？"**

---

## 七阶段处理流程 (Seven-Phase Pipeline)

```
原始数据 ──► Phase 1 ──► Phase 2 ──► Phase 3 ──► Phase 4
                                                     │
  报告  ◄── Phase 7 ◄── Phase 6 ◄── Phase 5 ◄───────┘
```

### Phase 1: 含糖筛选 (Sugar Filtering)

**解决的化学问题:** COCONUT 数据库有 40 万+化合物，大部分不含糖。我们需要精准筛出所有含糖苷的分子。

**实现:** 基于 RDKit 的糖环检测 (五元/六元含氧杂环 + 多羟基特征)，从原始数据中筛选出 94,242 条含糖化合物，按 `standard_inchi_key` 严格去重。

**产出:** `Coconut_Sugar_Check.csv` — 全部后续分析的起点。

---

### Phase 2: 糖-苷元精确切分 (Glycan–Aglycon Cleavage)

**解决的化学问题:** 要分别研究糖链和苷元，必须在**正确的位置**切开分子。切错了（比如切断了 C5-C6 键），葡萄糖就变成了木糖，整个后续分析全废。

**实现:**
1. 定位异头碳 (Anomeric Carbon, C1)
2. 识别 C1 发出的糖苷键 (C-O/C-N/C-S)
3. 区分"糖-苷元键"和"糖-糖键"（保留糖链完整性）
4. 碳原子守恒断言: C(糖) + C(苷元) = C(总)

**产出:** `Glycan_SMILES` + `Aglycon_SMILES` — 保证原子不丢失的精确碎片。

---

### Phase 3: 核苷酸糖 & 糖肽识别 (Nucleotide Sugar & Glycopeptide Detection)

**解决的化学问题:** 某些含糖化合物不是普通的苷，而是**核苷酸糖** (如 UDP-Glc, GDP-Man) 或**糖肽** (如万古霉素类)。它们的生物学意义完全不同，需要独立标记。

**实现:**
- 核苷酸: 碱基 (嘌呤/嘧啶 SMARTS) + 磷酸基团 **必须共存**
- 肽键: `[NX3]-C(=O)-[CX4;!CH3]` — 明确排除 NAc (N-乙酰基) 的甲基碳以防误判

**产出:** `Has_Nucleotide`, `Has_Peptide` 布尔标志 + 详细信息列。

---

### Phase 4: 分类学填补 (Taxonomy Enrichment)

**解决的化学问题:** 要研究"为什么这种植物选择了这种糖"，我们需要知道**化合物来自哪个物种、哪个科属**。COCONUT 原始数据中有约 11% 的分类学字段缺失。

**实现:** 使用 LOTUS 数据库的 NP Classifier 超类信息，通过 InChIKey Block-1 映射精确填补空缺。

**产出:** 补全的 `np_classifier_superclass` + 物种来源 `organisms` 字段。

---

### Phase 5: 糖序列 & 苷元骨架提取 (Sugar Sequence & Scaffold Extraction)

**解决的化学问题:** 把 SMILES 转化为化学家能看懂的语言。`OC[C@H]1OC(O)...` 需要变成 **"D-Glc-(α1→4)-D-Gal"**，苷元需要提取出 **Murcko 骨架**来做同族对比。

**实现:**
- 糖序列: 分层匹配 (精确单糖库 → 通用 Hex/Pen 兜底) + 连接键推断
- 苷元骨架: Bemis-Murcko 分解, 直链脂肪烃标记为 `Aliphatic Chain`
- Morgan 指纹: 2048-bit (radius=2), 用于 Phase 6 的 Tanimoto 回退

**产出:** `Sugar_Sequence`, `Murcko_Scaffold`, Morgan 指纹向量。

---

### Phase 6: 智能化学分类 (Intelligent Chemical Classification)

**解决的化学问题:** 化合物的大类归属 (黄酮、甾体、生物碱...)。这决定了后续分析中"哪些化合物比较"的分组维度。

**实现:**
1. SMARTS 规则匹配 (黄酮、香豆素、生物碱、甾体)
2. 糖脂特异捕获 (长链 C≥10 + 鞘氨醇/甘油骨架)
3. Tanimoto 回退 (对 SMARTS 未命中的，用 Morgan 指纹找最近邻)
4. Phase 3 覆盖: 核苷酸糖/糖肽强制归类

**产出:** `Superclass`, `Classification_Method`, `Glycolipid_Flag`。

---

### Phase 7: 可视化与报告 (Visualization & Reports)

**解决的化学问题:** 让化学家 **用眼睛** 验证计算结果是否合理。一张错误的高亮图能立刻暴露切分 Bug。

**实现:**
- 三色高亮: 🔴 糖核心 / 🟡 修饰基团 / 🔵 苷元
- BFS 洪水填充归属 (以糖苷键为边界，不是"环+1邻居")
- HTML/Excel 报告嵌入 Base64 分子图

**产出:** `Sample_Visual_Report.html`, 可内嵌的 PNG 分子图。

---

## 糖链修饰扫描 (Glycan Modification Scanning)

贯穿 Phase 2→7 的独立模块。仅扫描 **Glycan_SMILES** (不碰苷元), 识别 7 类修饰:

| 修饰 | SMARTS Pattern | 化学意义 |
|------|----------------|----------|
| **NAc** | N-C(=O)-CH₃ | GlcNAc/GalNAc 标志 |
| **O-Ac** | O-C(=O)-CH₃ | 乙酰化酯 |
| **Sulfate** | O-SO₃⁻ | 硫酸化 (肝素等) |
| **Phosphate** | O-PO₃²⁻ | 磷酸化 (Man-6-P) |
| **O-Me** | O-CH₃ | 甲醚化 |
| **COOH** | C(=O)OH | 羧基 (醛糖酸) |
| **NH₂** | 游离氨基 | 氨基糖 |

---

## 项目目录结构 (Project Structure)

```
D:\Glycan_Database\
├── data/                      # 原始和中间数据
├── lib/                       # 核心计算库
│   ├── sugar_utils.py         # 糖环识别
│   ├── glycosidic_cleavage.py # Phase 2: 精确切分
│   ├── phase3_secondary_scan.py # Phase 3: 核苷酸/肽键
│   ├── phase5_features.py     # Phase 5: 糖序列+骨架
│   ├── phase6_classifier.py   # Phase 6: 智能分类
│   ├── phase7_visualizer.py   # Phase 7: 三色可视化
│   ├── glycan_modifications.py# 修饰基团 SMARTS 扫描
│   └── chemical_dictionaries.py # 单糖/SMARTS 字典
├── scripts/                   # 可执行脚本
│   ├── run_glyconp_pipeline.py # 一键全量管线
│   └── analyze_glycan_frequencies.py # 科学分析
├── reports/                   # 输出报告
└── tests/                     # 单元测试
```

---

## 依赖 (Dependencies)

| 包 | 版本 | 用途 |
|----|------|------|
| `rdkit` | ≥2023.03 | 分子操作、SMARTS 匹配、指纹、Murcko 骨架 |
| `pandas` | ≥1.5 | 数据表操作 |
| `tqdm` | ≥4.65 | 进度条 |
| `xlsxwriter` | ≥3.0 | Excel 图片嵌入 |
| `numpy` | ≥1.24 | 数值计算 |

---

## 快速开始 (Quick Start)

```bash
# 全量运行 (约 30 分钟)
python scripts/run_glyconp_pipeline.py

# 限量测试 (1000 条, 约 20 秒)
python scripts/run_glyconp_pipeline.py --limit 1000

# 生成科学分析报告
python scripts/analyze_glycan_frequencies.py
```

---

## 许可 (License)

本项目仅供学术研究使用。COCONUT 数据遵循其原始许可协议。
