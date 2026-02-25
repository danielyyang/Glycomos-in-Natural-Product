# Glycan Database 升级计划 (Glycan Database Upgrade Plan)

## 0. 项目初始化与架构优化 (Project Initialization & Architecture)
经梳理，根目录含有较多测试与排错日志文件。为了使项目更清晰，我们建立以下标准文件夹规则并归档，以便于上传 GitHub：
- `data/`: 存放唯一的原始输入数据集 (如 `Coconut.csv`)。删除了冗余的 raw 嵌套文件夹。
- `reports/`: 所有经过程序处理的中间产物、清洗后的最终数据集 (如 `Coconut_Sugar_Final.csv` / `Excel`) 统一输出至此。
- `lib/`: 存放核心支撑库 (如 `sugar_utils.py`, `visualizer.py`)。
- `scripts/`: 存放所有的核心流程和解析脚本 (如 `unified_pipeline.py`)。过往测试用的冗余脚本均已移入 `scripts/archive/`。
- `docs/` 和 `images/`: 存放说明文档与可视化输出图形。
- `log/`: 集中收集所有操作记录。

## 1. 核心库注入 (Injected Libraries)
- **核碱基库 (Base Library)**: 内置 A, T, C, G, U 的碱基 SMARTS，用于检测五元糖的特异性 N-糖苷。
- **氨基酸库 (Amino Acid Library)**: 包含天然氨基酸侧链的明确 SMARTS，能够抓取到精确的氨基酸肽链连接结构。

## 2. 核心功能与 BUG 修复逻辑 (Core Functions & Bug Fix Logic)

### A. [BUG FIX] 糖链 vs 稠环/螺环判定 (Glycan vs Polycyclic Framework)
- **识别规则**: 两个含氧环必须通过 **非环内桥接原子** (`COC`, `CNC`, `CSC`, `CCC`) 连接。
- **排除规则**: 若两个环共用原子 (连接键对于 `Bond.IsInRing()` 判定为 `True`)，则整个框架标记为 `Polycyclic_Organic_Framework`，严禁按多糖序列粗暴拆分。

### B. [UPGRADE] 核苷酸 N-糖苷识别 (Nucleotide Recognition)
- **判定逻辑**: 检测并追踪五元糖环 (呋喃糖) 的异头碳 `C1` 位置。
- **条件**: 如果 `C1` 直接连接到一个 `N` 原子，且该含 `N` 环系符合 Base Library，则视为核苷类结构。
- **分类**: 若结构中带有磷酸基团，标记为 `Nucleotide`；否则标记为 `Nucleoside`。

### C. [NEW] 氨基酸绿色标注系统 (Amino Acid Green Tagging System)
- 在终端处理日志中，只要捕获到氨基酸/肽类片段，立即使用绿色字体 (`\033[92m`) 高亮打印。
- 自动写入持久化日志至 `log/green_tags.txt`，存储格式为：`Compound_ID | AA_Type | Connection_Point`。

### D. [SCANNED] 全局 Aglycan 骨架分析 (Global Aglycan Scaffold Analysis)
- 精准定位并仅提取连在 `C1` 上的非糖与非氨基酸部分 (即 `Primary_Aglycan`)。
- 生成其特有的 `Murcko Scaffold` 并求解 `InChIKey`。取该键的首段 (First Segment) 建立该分子的“家族 ID”。

## 3. 预期产出 (Expected Outputs)
- 输出升级版数据集：`output/Sugar_Sort_Refined.CSV`，该数据集增加或重构以下列：
  - `Refined_Sequence`: 优化后的含结构修正的修饰序列。
  - `Aglycan_Scaffold_ID`: 本配糖体所属的结构大类 InChIKey。
  - `Linker_Type`: 标注是否带有肽键与侧链特征。
  - `Is_Nucleotide`: 布尔值或标注。
- 对未能分类 (`No Classified`) 的记录，通过大样本分析同一骨架 ID 的聚集性，给出智能分类建议。

## 4. 深度重分类与图像归档 (Phase 4: Deep Reclassification & Image Sorting)
- **子集筛查 (Subset Filtering)**: 筛选来源表中 `Sheet_Source` 为 `No Classified` 的所有目标分子。
- **降维比对 (Dimensionality Reduction & Mapping)**: 将已识别且带有智能标注提示的 `Scaffold_Class_Hint`，或是大样本中与已知其他分类 (`Steroid`, `Flavonoid` 等) 完全重叠的 `Aglycan_Scaffold_ID` 用作其新归属大类。
- **表格改写 (Excel Modification)**: 开启原本的 `Coconut_Sugars_Unique.xlsx`，建立一个全新的 `New_Classfied` Sheet 并填入这批成功被溯源分类的数据，新增 `Classification` 列记载匹配结果。
- **图像物理分类 (Physical Image Sorting)**: 对被闲置在 `images/No_Classificed/` 这座“孤岛”中的数百/数千张结构分布图进行扫描，并根据分配的新类别创建对应的物理子级文件夹，一键归档。
