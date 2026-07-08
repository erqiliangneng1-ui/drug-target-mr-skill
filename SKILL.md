---
name: drug-target-mr
description: >
  药靶孟德尔随机化 (Drug-target MR) 完整工作流。涵盖 GWAS 人群匹配评估、
  GTEx/eQTLGen 数据选择、单SNP vs 多SNP工具变量策略、结局 GWAS 数据源评估
  (UK Biobank, FinnGen, BBJ)、MR 方法选择 (IVW, MR-Egger, Weighted median)。
  专为生物医学研究者设计。
category: research
---

# Drug-target MR (药靶孟德尔随机化)

## 使用场景

- 用户问孟德尔随机化、药靶 MR 或 cis-eQTL MR
- 需要评估暴露和结局之间的 GWAS 祖源兼容性
- 用 GTEx eQTL 数据作为暴露做 MR
- 需要判断结局 GWAS ID 是否适合当前 MR
- 需要决定工具变量策略（lead SNP vs 多SNP + LD clumping）

---

## 一、人群匹配检查

### 1. 暴露数据祖源

**GTEx v8/v10/v11**：
- 约 85% 欧洲裔、13% 非裔、1% 亚裔
- 没有分层分析的 eQTL 数据，整个队列一起分析
- eQTL 效应反映欧洲人群的 LD 结构和等位基因频率

**eQTLGen**：
- 约 31,000 样本，主要是欧洲裔
- 使用 Illumina 芯片（杂交探针），不如 RNA-seq 灵敏
- 部分神经肽基因因表达量低未被覆盖

### 2. 结局 GWAS 祖源

**与 GTEx 兼容的欧洲裔来源**：

| 来源 | 祖源 | 说明 |
|------|------|------|
| UK Biobank (ieu-b-* / ukb-b-*) | >99% 英国白人 | 自报表型 |
| FinnGen (finn-b-*) | 芬兰人（欧洲裔） | ICD 编码，R11 ~37.7万样本 |
| 欧洲荟萃分析 | 欧洲多国 | 需查原文献 |

**跨人群研究中需要注意**：
- Sakaue 2021 (Nat Genet) 混合了 FinnGen + BBJ（日本生物银行）
- 但在 **IEU OpenGWAS 中可能只托管了欧洲子集**
- 验证方法：查看 IEU 页面 `population` 字段是否为 "European"，`ncase` 是否与欧洲子集一致
- 示例：`ebi-a-GCST90018883` 显示 `population: European, ncase: 5093` — 这是纯 FinnGen 欧洲子集 ✅

**与 GTEx 不兼容的来源**：
- Biobank Japan (BBJ) — 纯日本人
- China Kadoorie Biobank (CKB) — 中国人
- 需要对应的东亚 eQTL 数据

### 3. 决策树

```
暴露数据是什么祖源？
├── 欧洲裔 (GTEx, eQTLGen) ─────────────────────
│   结局 GWAS 是什么祖源？
│   ├── 欧洲裔 → ✅ 兼容，直接做 MR
│   ├── 跨人群 (EUR+EAS混合) → ⚠️ 有条件兼容
│   │   ├── >95% 欧洲裔 → 可用，讨论中注明
│   │   └── <95% 欧洲裔 → ❌ 找纯欧洲替代
│   └── 东亚裔 → ❌ 不兼容，需东亚 eQTL 数据
│
└── 东亚裔 → 需东亚结局 GWAS
```

### 4. FinnGen 查找指南

**FinnGen R11**（旧版，样本量较小）：
- 访问 https://r11.finngen.fi/
- 按表型名搜索
- IEU OpenGWAS ID 格式：`finn-b-<PHENOCODE>`

**FinnGen R13**（最新版，推荐直接下载原始数据）：
- Manifest 文件 `finngen_R13_manifest.tsv` 列出所有表型
- 直接下载格式：
  ```
  https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_<PHENOCODE>.gz
  ```
- R13 样本量通常比 R11 大 10-40%
- 文件列：`#chrom`, `pos`, `ref`, `alt`, `rsids`, `nearest_genes`, `pval`, `mlogp`, `beta`, `sebeta`, `af_alt`, `af_alt_cases`, `af_alt_controls`
- 参考基因组：**GRCh38**

### 5. IEU OpenGWAS VCF 下载文件类型

| 文件 | 用途 | 需要？ |
|------|------|:------:|
| `<id>.vcf.gz` | **插补版 GWAS（100万-1000万 SNP，MR 标准用）** | **是** |
| `<id>.vcf.gz.tbi` | Tabix 索引 | **是** |
| `<id>_raw.vcf.gz` | 只有芯片上直接分型的 SNP（约 50万） | 否 |

**一定要下载不带 `_raw` 的版本**。`_raw` 版会漏掉大多数 lead SNP。

---

## 二、工具变量 (IV) 策略

### 标准做法：lead cis-eQTL SNP

- 每个靶基因选 P 值最小的 cis 区域 SNP
- 优点：简单、可解释，药靶 MR 的标准做法
- F 统计量 > 10（弱工具变量阈值）
- 结果：**N 个基因 = N 个 IV**

### 替代做法：多 SNP + LD clumping（有条件时）

- 保留所有显著 cis-eQTL SNP（FDR < 0.05）
- LD clumping：r² < 0.001，窗口 10,000kb，P1 < 5e-8
- 需要：PLINK + 1000 Genomes 参考面板（欧洲裔）
- **典型产出**：1000 个候选 → 几十到上百个独立 SNP
- PLINK 安装参考：`references/install-plink.md`

#### PLINK 2.0 LD clumping 技巧

当参考面板用 rsID 但 GTEx 用 `variant_id` 时，用 `--set-all-var-ids` 统一为 `chr:pos`：

```bash
plink2 --bfile <面板路径> \
  --set-all-var-ids @:# \
  --clump <汇总统计文件> \
  --clump-p1 5e-8 --clump-p2 5e-8 \
  --clump-r2 0.001 --clump-kb 10000 \
  --out <输出前缀>
```

#### ⚠️ 面板基因组版本不匹配

`1kg.v3` 面板（来自 IEU OpenGWAS）是 **GRCh37**。GTEx v11 和 FinnGen R13 是 **GRCh38**。用 GRCh37 面板对 GRCh38 数据做 clumping，结果完全错误。

**检测方法**：查一个已知 SNP 的坐标。`rs429358`（APOE ε4）：
- GRCh37：chr19:45411941
- GRCh38：chr19:44908684

**解决方法（优先顺序）**：
1. 下载 GRCh38 面板：`wget https://ctg.cncr.nl/software/magma/ref_data/g1000_eur.zip`
2. LiftOver：将暴露 SNP 坐标从 GRCh38 转 GRCh37，clumping 后再转回
3. 降级方案：距离剪枝（不基于 r²，只是按距离过滤，精度较低）

---

## 三、网络受限环境的替代方案

当环境无法访问 OpenGWAS API 时：

1. 在联网机器上用 `ieugwasr::extract_outcome()` 提取 lead SNP
2. 导出 CSV
3. 拷回本机
4. 用 `data.table::fread()` 加载
5. 用 `TwoSampleMR::harmonise_data()` 继续处理

---

## 四、流水线结构与职责分离

脚本按顺序执行，每个脚本只做一件事：

```
0_0folder_initi.R          →  配置：路径、基因、IV策略、阈值
        ↓
0_1data_formating.R        →  数据进去，标准化列名出来
                               (不做LD clumping，不做F-statistic)
        ↓
1_correlation_analysis.R   →  P < 5e-8 过滤 + 3个输出：
                               曼哈顿图、基因散点图、基因统计CSV
        ↓
2_linkage_disequilibrium.R →  LD clumping (PLINK 2.0 + 1000G面板)
        ↓
3_remove_weak_IV.R         →  F 统计量过滤（F > 10）
        ↓
4_remove_confounder.R      →  混杂因素剔除（FastTraitR）
        ↓
5_do_MR.R                  →  MR 分析（harmonise + mr()）
```

**核心原则**：一个脚本做一件事。不要在上游脚本里做下游脚本的事（比如在 `0_1` 里做 LD clumping 或 F 检验）。

**修改顺序**：按流水线顺序逐一修改。改下游脚本前先确认不会与上游职责重叠。

### 结果报告规范

每个脚本跑通后报告：
1. **输出目录** — 完整路径
2. **生成的文件** — 文件名 + 说明
3. **关键数字** — SNP 数量变化、F 值范围等

---

## 五、GTEx eQTL + FinnGen 数据格式化详解

### 暴露数据：GTEx parquet → TwoSampleMR exposure 格式

1. `arrow::read_parquet()` → 按靶基因过滤
2. `parse_gtex_variant_id()`：`chr1_665098_G_A_b38` → CHR、BP、ref、alt、SNP_ID
3. IV 策略选择
4. 列名标准化为 TwoSampleMR 格式：`SNP, CHR, BP, A1, A2, BETA, SE, FRQ, P, N, gene`
5. `data.table::fwrite()` → `FMT_DIR/` 下面

#### ⚠️ phenotype_id 版本号后缀匹配

GTEx 的 `phenotype_id` 带版本号后缀（如 `ENSG00000006128.12`），但配置文件中是裸 ENSG ID（如 `ENSG00000006128`）。直接 `%in%` 匹配结果为 0 行。

**修复**：先去掉后缀再匹配：
```r
gtex_data <- gtex_data %>%
  mutate(phenotype_id_clean = gsub("\\.[0-9]+$", "", phenotype_id)) %>%
  filter(phenotype_id_clean %in% target_ensg)
```

#### ⚠️ file.exists 检查时需返回路径而非字符串

`format_gwas_data()` 中文件已存在时的 `return()` 需要返回实际路径，否则下游 `file.exists()` 检查失败。

### 结局数据：FinnGen TSV.gz → TwoSampleMR outcome 格式

1. `data.table::fread(cmd = "gzcat <文件>")` + `select=` 只取需要的列
2. chr:pos 键值匹配
3. 等位基因对齐（正向匹配、互换匹配翻转 β 符号）
4. 输出列：`SNP, chr, pos, effect_allele.outcome, beta.outcome, ...`
5. 写入 `OUTCOME_DIR/`

#### ⚠️ ~ 路径在 shell 中不展开

`fread(cmd = "gzcat ~/Desktop/MR/...")` 中 `~` 不会被 shell 展开。**一定要用 `path.expand()`**。最佳实践：在配置文件中直接用绝对路径或 `path.expand('~/Desktop/...')`。

### 相关性分析（1_correlation）

三个输出：

| 输出 | 说明 |
|------|------|
| `exposure.pvalue.csv` | P < 5e-8 过滤后的 SNP 列表 |
| 曼哈顿图 (PNG) | 全基因组概览 |
| 基因散点图 (PNG) | x=基因, y=-log10(P)，对小基因集更有用 |
| `exposure_gene_summary.csv` | 每基因统计：n_SNPs, min_P, lead_SNP, mean_BETA, mean_FRQ |

#### ⚠️ SNP_ID vs SNP 列名

CSV 中的列名叫 `SNP`（chr:pos），不是 `SNP_ID`。做基因汇总时用 `SNP[which.min(P)]`，不要用 `SNP_ID`。

### LD clumping（2_linkage）

常见陷阱一览：

| 问题 | 修复 |
|------|------|
| **CLUMP_KB/R2 被硬编码覆盖** — 原模板写死了 CLUMP_KB=500, R2=0.1 | 删掉这些行，用配置文件的全局变量 |
| **bfile 路径双后缀** — `file.path(LD_REF, EXPOSURE_POP)` 拼出 `.../EUR/EUR` | 直接用 `LD_REF` |
| **缺 library()** — 没有导入 dplyr, ieugwasr, plinkbinr | 补上 |
| **GRCh38 非标染色体** — PLINK 不认识 PAR1, chr1_KI270706v1_random 等 | **绕过 ieugwasr，直接 system() + --allow-extra-chr** |
| **bim 中重复 rsID** — 同一 rsID 出现在不同位置 | `make.unique(paste0(chr, ":", bp))` 生成唯一 ID |
| **--extract range 格式** — 必须是空格分隔的 4 列 | `fwrite(..., sep = " ")` + label 列 |
| **--clump 输入需表头** — PLINK 要求第一行是 `SNP P` | `writeLines("SNP P", file)` 写表头 |
| **大面板太慢** — 7493 万位点 × 单线程跑几个小时 | 切迷你面板 |

**迷你面板提取**：
1. 读 SNP 位置
2. 按染色体算 `min(pos)-WINDOW` 到 `max(pos)+WINDOW`
3. 合并重叠区间
4. 写 range 文件 → `plink --extract range` → 得到 ~12K 变体的小面板
5. 在这个小面板上做 clumping，秒级完成

---

## 六、弱 IV 剔除（3_remove_weak_IV）

#### ⚠️ Fval vs F 列名

原模板计算了 `exposure_data$Fval` 但过滤时用了 `exposure_data$F`（无此列），结果全被删除。**用 `$Fval` 过滤**。

#### ⚠️ 硬编码 F_THRESHOLD

原模板第 1 行写死了 `F_THRESHOLD <- 10`，覆盖了配置文件的值。**删掉这行**，用配置中的全局变量。

---

## 七、混杂剔除（4_remove_confounder）

通过 `FastTraitR::look_trait()` 查询各 SNP 是否关联到混杂表型（如 BMI、吸烟、免疫指标）。

### ⚠️ FastTraitR 需要 rsID，不是 chr:pos

`look_trait()` 查询 IEU OpenGWAS API，用 rsID 作为查找键。chr:pos 格式会返回"输入的rs号没有检索到数据"。

**修复**：从结局文件提取 rsID ↔ chr:pos 映射，用 named vector 做转换：
```r
rsid_lookup <- setNames(mapping$SNP, mapping$SNP_ID)
snp_ids_rs <- rsid_lookup[exposure_data$SNP]
```

### ⚠️ eQTL SNP 返回 NULL 时的优雅降级

GTEx eQTL SNP（调节性变异）往往不在 GWAS 数据库中，`look_trait()` 返回 NULL 是正常现象。此时应：
1. 创建空的 `#confounder_SNPs.txt`
2. 全部 SNP 保留不变
3. 日志中说明跳过原因
4. **不要**尝试 PhenoScanner 或其他 API 备用（有同样的 GWAS 数据库偏差）

---

## 八、MR 分析（5_do_MR）

### 常见陷阱

| 问题 | 修复 |
|------|------|
| **结局路径错** — 模板读 `FMT_DIR/` 但结局在 `OUTCOME_DIR/` | 用 `OUTCOME_DIR` |
| **结局列名错** — 模板用 `A1, A2, BETA` 但结局是 outcome 格式 | 用 `effect_allele.outcome`, `beta.outcome` 等 |
| **SNP 匹配错** — 用 chr:pos 匹配 rsID（100% 失败） | 用 `SNP_ID`（chr:pos）列 |
| **单结局模板** — 原模板只处理一个结局 | 改造为函数 `run_mr(outcome_file, outcome_name, subdir)` |
| **MR-PRESSO 易崩溃** — 参数不匹配会卡死整个 MR | 用 `tryCatch()` 包裹 |

### 每个结局的输出

保存到 `MR_DIR/<子目录>/`：

| 文件 | 内容 |
|------|------|
| `MR-Result.csv` | IVW、MR-Egger、Weighted median 等 + OR |
| `heterogeneity.csv` | Cochran's Q 检验 |
| `pleiotropy.csv` | MR-Egger 截距检验 |
| `MR-PRESSO.txt` | 异常值检测（可选） |
| `scatter_plot.pdf` | 散点图 |
| `forest.pdf` | 森林图 |
| `funnel_plot.pdf` | 漏斗图 |
| `leaveoneout.pdf` | 留一法分析 |

---

## 九、错误报告格式（4 步法）

用户要求每次报错都按这个格式：

1. **分析原因** — 错误信息实际含义、触发的根因，不用抽象术语
2. **修改方案** — 改前/改后代码对比，说明为什么修好了
3. **背景补充** — 相关工具/格式/概念的知识，用具体例子而非抽象术语
4. **完成修改** — 应用修复

⚠️ **永远不要在脚本中用 `quit(save="no")`** — 这会杀死整个 R 进程，在 GUI 环境（RStudio、Hermes Desktop、Jupyter）中导致 fatal error 崩溃。只在批处理模式（`Rscript`）中安全。用 `if/else` 分支或 `return()` 替代。

---

## 十、MR 结果解读

### 多重校正在药靶 MR 中的应用

| 场景 | 结局数 | Bonferroni α | 判断 |
|------|:------:|:------------:|------|
| 主分析 + 敏感性 | 2（AR、CRSwNP） | 0.05/2 = **0.025** | P=0.047 不显著 |
| 主分析 + 2个敏感性 | 3 | 0.05/3 = **0.017** | |
| 基因 × 结局 | 5×2 = 10 | 0.05/10 = **0.005** | 过于保守，有争议 |

**措辞约定**：
- P < 校正后 α → "有统计学显著性"
- P < 0.05 但 > 校正后 α → "名义显著，提示性关联"
- P > 0.05 → "未观察到关联"（避免绝对表述）

### 讨论中需说明的已知弱点

1. **多重校正未做** — 必须做。P=0.047 在 Bonferroni 后不显著
2. **基因打包混杂** — 抗炎/促炎基因合并相互抵消。**必须拆分做基因特异性 MR**
3. **组织特异性局限** — 肺 eQTL ≠ 鼻黏膜 eQTL
4. **缺 COLOC** — 无法排除 LD 污染
5. **阴性结论勿过于绝对** — "未观察到"而非"不存在"

---

## 十一、基因特异性 MR（Gene-specific MR）

### 何时使用

- 打包版 MR 得到边界性结果（如 P=0.047）— 想知道哪个基因在驱动信号
- 不同基因的 β 方向不一致（如 NMU 正向、VIP 负向）— 互相抵消
- 审稿人要求展示每个基因的结果
- 每个基因有 1-4 个 SNP

### 输出结构

```
MR_DIR/gene_specific/
├── gene_MR_summary.csv          ← 所有基因 × 结局汇总表
├── ADCYAP1/
│   ├── AR/ (MR-Result.csv, 异质性, 多效性, MR-PRESSO)
│   └── CRSwNP/ (同上)
├── GRP/...
├── NMU/...
└── NGF/ (1 SNP → 仅 Wald ratio)
```

### SNP 数量解读

| SNP 数 | 能做的 | 限制 |
|:------:|--------|------|
| 1 | 仅 Wald ratio | 无法做异质性/多效性检验，CI 通常宽 |
| 2-3 | 全部 MR 方法 + Q 检验 | Egger 截距统计效力不足 |
| 4+ | 完整流程 | IVW 可靠，MR-PRESSO 可行 |

---

## 十二、三阶段升级计划

### Phase 1 🔴 基因拆分（立刻做，无需新数据）
- 脚本：`5_do_MR1_gene_specific.R`
- 几小时跑完

### Phase 2 🟡 组织验证（需下载新数据）
- 第二个组织 eQTL（如 GTEx Skin — 皮肤神经末梢丰富）
- 新增结局（如特应性皮炎 AD）
- 逻辑：如果神经肽→过敏通路是真实的，皮肤 eQTL → AD 应为阳性

### Phase 3 🟢 COLOC + 阳性对照
- 共定位分析排除 LD 污染
- UK Biobank 嗜酸性粒细胞计数 / IgE 水平
- 如果 IV 能预测过敏生物标志物，通路就得到了生物学验证

---

## 十三、论文与 PPT 生成

### 论文结构（manuscript_summary.md）

1. **背景** — 神经肽在神经源性炎症中的作用、靶基因选择依据、药靶 MR 原理
2. **方法** — eQTL 来源 (GTEx v11 Lung)、IV 策略、结局来源 (FinnGen R13)、MR 方法、软件
3. **结果** — 按分析阶段列表、方法级森林图、敏感性检验汇总
4. **讨论** — 优势（本地数据、多结局、药靶 MR 框架）、局限（组织特异性、无 COLOC、样本量）
5. **结论** — 基因优先级排序、治疗意义、下一步建议

### PPT 结构（8 页）

1. 标题 + 背景
2. 研究设计图
3. 靶基因选择（文献证据）
4. eQTL 结果（基因散点图 + 汇总表）
5. MR 结果（方法级森林图 — 核心图）
6. 敏感性分析
7. 基因特异性 MR 热图
8. 局限 + 下一步

### 画图原则

- **曼哈顿图和基因散点图都要生成**。曼哈顿图符合期刊格式，基因散点图（x=基因, y=-log10(P)）对分析者更直观
- 用户习惯："图多我后期写作时可以删减"
- 输出 `exposure_gene_summary.csv`：基因、SNP数、最小P、lead_SNP、平均BETA、平均FRQ。0 SNP 的基因在表格中展示为 `n_SNPs = 0`

---

## 参考资料

- `references/gwas-ancestry-compatibility.md` — GTEx 祖源与 GWAS 人群数据详细笔记
- `references/install-plink.md` — macOS 安装 PLINK 1.9/2.0 指南
- `references/ld-debugging-narrative.md` — LD clumping 7 次报错调试全程记录
- `references/plink-pgen-to-bed.md` — 1000 Genomes pgen → PLINK 1.9 bed 转换
- `references/plink-reference-panel-errors.md` — hg38 面板常见报错汇总
- `references/project-structure-neuropeptide.md` — 神经肽 MR 项目专属结构说明
- `references/publishing-workflow.md` — 将 skill 发布到 GitHub 的工作流
