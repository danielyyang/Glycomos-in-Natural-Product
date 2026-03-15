# 项目状态报告 — 供专家咨询用
# Project Status Report — For Expert Consultation

> 更新日期 / Updated: 2026-03-14

---

## 一、项目定位 (Project Positioning)

**GlycoNP-Pipeline** 是一个面向糖化学与天然产物化学的自动化计算管线，核心目标是：
- 从大型天然产物数据库（COCONUT, ~94,000 条含糖化合物）中自动化地**识别、切分和注释糖苷结构**
- 为糖基化天然产物提供**系统化的骨架分析、序列标注和分类体系**

---

## 二、已实现功能 (Implemented Features) ✅

### Phase 1: 数据初始化与去重
- [x] 读取 COCONUT 全量 CSV (~660MB, 94,242 条含糖化合物)
- [x] 基于 `contains_sugar` 字段筛选含糖子集
- [x] 基于 `standard_inchi_key` 精确去重

### Phase 2: 核心结构切分 (Main Innovation)
- [x] **糖环识别**：自动检测五元/六元含氧环并验证 (sugar ring detection)
- [x] **不饱和环拒止**：仅检查环内键，环外 C=O 不影响判定
- [x] **多环聚醚拓扑铁律**：共享边的含氧环视为聚醚（非糖），自动拒绝
- [x] **糖键连方式约束**：仅允许 O/N/S/C/P 桥原子
- [x] **核苷排除**：呋喃糖环邻接嘌呤/嘧啶碱基时自动排除
- [x] **宏环糖豁免**：大环（≥8 原子环）交联不影响糖识别
- [x] **缩醛桥豁免**：C/O 组成的小环（1,6-脱水桥等）不影响糖识别
- [x] **糖苷/苷元分离**：基于分枝大小启发式+虚拟原子断键
- [x] **假阳性标记**：无法通过拓扑验证的化合物标记为 FALSE_POSITIVE

### Phase 3: 二级片段扫描
- [x] **核苷酸识别**：在糖基片段中扫描嘌呤/嘧啶碱基
- [x] **氨基酸/肽键识别**：在苷元中扫描肽键 (N-C=O-C) 模式

### Phase 4: 生物分类学填补
- [x] **多层级 API 搜索**：COCONUT → LOTUS/Wikidata → PubChem → Wikipedia
- [x] **GBIF 物种匹配**：organism → Family 自动映射
- [x] **本地缓存**：避免重复 API 调用

### Phase 5: 骨架提取与糖序列
- [x] **Morgan 指纹 (radius=2, 2048-bit)** 生成
- [x] **Murcko 骨架** 提取
- [x] **IUPAC-like 糖序列** 生成：`α-D-Glc-(1→4)-β-D-Gal` 格式
- [x] **修饰基团标注**：`*NAc, *S, *P, *deoxy` 等 `*` 前缀格式

### Phase 6: 智能化学分类
- [x] **chemical_super_class 推断**：InChIKey Block-1 匹配 + Tanimoto ≥ 0.8
- [x] **np_classifier_superclass 推断**：同策略
- [x] **Excel 分 Sheet 导出**：按 super_class 分组

### Phase 7: 可视化
- [x] **彩色 2D 分子图**：糖环=红、苷元=蓝、取代基=黄、核苷/氨基酸=浅绿
- [x] **Excel 图片嵌入**：Structure_Image 列
- [x] **按大分类/子分类分层文件夹** 输出图片

### 单糖库 (Monosaccharide Library)
- [x] **85+ 条** SMILES 条目（全立体化学闭环形式）
- [x] **21 条** SMARTS 子结构匹配模式
- [x] 覆盖己醛糖 D/L 全构型、脱氧糖、二/三脱氧糖、氨基糖、糖醛酸、戊糖、呋喃糖、酮糖、支链糖、庚糖

---

## 三、已知问题与待完善功能 (Known Issues & TODO) ⚠️

### 3.1 API 连通性问题
| API | 状态 | 原因 |
|-----|------|------|
| COCONUT API | ❌ 404 | API 端点已失效，需更新 URL |
| Wikidata SPARQL | ❌ 403 | 被网络防火墙拦截 |
| Wikipedia API | ❌ 403 | 被网络防火墙拦截 |
| PubChem REST | ✅ 200 | 可连通，但 Taxonomy 解析逻辑需优化 |
| GBIF | ✅ 200 | 正常工作 |

**影响**：约 49.4% 的含糖化合物 (46,574 条) 缺少 organism 信息，当前 API 搜索补全率 ≈ 0%。

### 3.2 功能缺口 (Feature Gaps)
| 功能 | 优先级 | 说明 |
|------|--------|------|
| **COCONUT v2 API 适配** | 🔴 高 | 更新失效的 COCONUT API 端点 |
| **PubChem Taxonomy 解析修复** | 🔴 高 | PUG View 返回了 Taxonomy Section 但解析逻辑不匹配 |
| **唾液酸 (Sialic Acid) 支持** | 🟡 中 | KDO, Neu5Ac 等九碳糖尚未加入库 |
| **糖脂自动分类** | 🟡 中 | 基于苷元骨架自动判定糖脂类型 (鞘糖脂/甾体皂苷/三萜皂苷等) |
| **异头碳 C1 还原端精确判定** | 🟡 中 | 当前使用 out-degree 启发式，可改为 BFS 最远端 |
| **开环糖/氧化糖** 处理 | 🟡 中 | 线形醛糖/酮糖和氧化修饰糖的识别 |
| **废弃测试文件清理** | 🟢 低 | tests/ 中有 ~40 个 debug_* 脚本和 ~5 个引用废弃函数的测试 |
| **CI/CD 集成** | 🟢 低 | GitHub Actions 自动测试 |

### 3.3 废弃文件 (Deprecated Files)
以下文件已被标记为 DEPRECATED，保留仅供历史参考：
- `lib/secondary_fragments.py` → 功能已由 `sugar_utils.py` 替代
- `scripts/utils/run_pipeline.py` → 引用不存在的模块 (aggregator, profiler)
- `lib/archive/` → 早期版本的 aggregator, profiler, researcher

---

## 四、项目总流程图 (Overall Workflow)

```
┌─────────────────────────────────────────────────────────────┐
│                    COCONUT Database (660MB)                  │
│              ~400,000 天然产物, ~94,000 含糖                │
└───────────────────────┬─────────────────────────────────────┘
                        │ Phase 1: 含糖筛选 + InChIKey 去重
                        ▼
┌─────────────────────────────────────────────────────────────┐
│              Coconut_Sugar_Check.csv (94,242 条)            │
└───────────────────────┬─────────────────────────────────────┘
                        │ Phase 2: 核心结构切分
                        │ ┌─────────────────────────────────┐
                        │ │ 糖环识别 (五/六元含氧饱和环)    │
                        │ │ ↓                               │
                        │ │ 拓扑验证 (聚醚铁律/键连约束...)  │
                        │ │ ↓                               │
                        │ │ 糖苷/苷元分离 (Glycan/Aglycan)  │
                        │ └─────────────────────────────────┘
                        ▼
┌──────────────────┐    ┌──────────────────┐
│  Glycan_SMILES   │    │ Aglycan_SMILES   │
│  (糖基部分)      │    │ (苷元部分)        │
└────────┬─────────┘    └────────┬─────────┘
         │ P3: 核苷酸扫描        │ P3: 氨基酸/肽键扫描
         ▼                       ▼
         │                       │ P5: Morgan FP + Murcko 骨架
         │ P5: 单糖鉴定+序列     │
         │  (85+条库匹配)        │
         ▼                       ▼
┌──────────────────────────────────────────────┐
│                 P4: API 分类学填补            │
│  COCONUT → LOTUS → PubChem → Wikipedia → GBIF│
├──────────────────────────────────────────────┤
│                 P6: 智能分类推断              │
│  InChIKey Block-1 匹配 + Tanimoto ≥ 0.8      │
├──────────────────────────────────────────────┤
│                 P7: 可视化导出                │
│  彩色 2D 图 + Excel 分 Sheet + 分层文件夹     │
└──────────────────────────────────────────────┘
```

---

## 五、未来可拓展方向 (Future Directions) 🔮

### 近期 (Short-term)
1. **修复 API 连通性** — 更新 COCONUT v2 端点 + 优化 PubChem Taxonomy 解析
2. **唾液酸&特殊糖扩展** — Neu5Ac, KDO, KDN 等九碳糖加入单糖库
3. **开环糖/醛糖酸处理** — 线形醛糖 (aldose) 和醛糖酸 (aldonic acid) 的识别

### 中期 (Mid-term)
4. **糖脂自动分类器** — 基于苷元骨架自动判定: 鞘糖脂 (Sphingolipid)、甾体皂苷 (Steroidal Saponin)、三萜皂苷 (Triterpenoid Saponin)、黄酮苷 (Flavonoid Glycoside) 等
5. **糖苷键类型精确标注** — α/β-1→2, α-1→3, β-1→4, α-1→6 等
6. **与 GlyTouCan 数据库集成** — 自动比对标准糖链序列编号
7. **图神经网络 (GNN) 糖苷分类** — 利用分子图特征进行机器学习分类

### 长期 (Long-term)
8. **糖苷差异组学 (Glycomics)** — 不同物种/科属间的糖基化修饰模式比较
9. **活性预测** — 基于糖链结构特征预测生物活性 (抗肿瘤/免疫调节等)
10. **Web 界面** — 在线版管线，支持单分子查询和批量分析
11. **知识图谱** — 糖苷-物种-活性-靶点关系网络
