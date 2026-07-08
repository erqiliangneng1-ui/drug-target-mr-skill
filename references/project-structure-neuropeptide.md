# Project structure: neuropeptide eQTL → AR/CRSwNP

## Pipeline scripts (ordered by execution)

| Script | Role | Key output |
|--------|------|------------|
| `~/Desktop/MR/0_0folder_initi.R` | Config: IV_STRATEGY, target genes, GWAS paths, directories | Global variables |
| `~/Desktop/MR/0_1data_formating_zhou_v2.R` | **Data formatting only** — GTEx parquet + FinnGen TSV.gz → TwoSampleMR format. NO LD clumping, NO F-statistic | `_formatted.csv.gz` (FMT_DIR), `_formatted.csv.gz` (OUTCOME_DIR) |
| `~/Desktop/MR/1_correlation_analysis_zhou.R` | P < 5e-8 filter + Manhattan plot + gene scatter plot + gene summary CSV | `exposure.pvalue.csv`, `exposure_gene_summary.csv`, 2 PNGs |
| `~/Desktop/MR/2a_extract_mini_ld_panel.R` | Extract ~1.48M variant mini-panel from 11GB hg38 panel | `EUR_hg38_v19_mini_10mb.{bed,bim,fam}` |
| `~/Desktop/MR/2_linkage_disequilibrium_analysis_v3.R` | LD clumping (PLINK + mini panel, --allow-extra-chr) | `exposure.LD.csv` (357→11 SNPs) |
| `~/Desktop/MR/3_remove_weak_IV.R` | F-statistic filter (F > 10) | `exposure.F.csv` (11/11 passed, F: 31.9–148.3) |
| `~/Desktop/MR/4_remove_confounder.R` | Confounder removal via FastTraitR | `exposure.confounder.csv` (11/11 retained, API no data) |
| `~/Desktop/MR/5_do_MR1.R` | MR analysis (pooled, all methods + sensitivity) | AR P=0.047 / CRSwNP P=0.78 |
| `~/Desktop/MR/5_do_MR1_gene_specific.R` | MR by individual gene (★ **next to run**) | gene_specific/ subdir with per-gene results |
| `~/Desktop/MR/5_do_MR2.R` | Method-level forest plot | PDF per outcome |
| `~/Desktop/MR/manuscript_summary.md` | Paper draft | Full article structure |
| `~/Desktop/MR/神经肽eQTL_MR研究.pptx` | 8-slide academic PPT | Presentation |
| `~/Desktop/MR/progress.md` | Project progress log | Updated 2026-07-07 |

## IV strategy determination

The **single source of truth** for which IV strategy is active is the `IV_STRATEGY` variable in the config script:

```r
# In 0_0folder_init_zhou.R:
IV_STRATEGY <- 'all_signif'  # dozens to ~100 IVs after LD clumping (in 2_ld.R)
# IV_STRATEGY <- 'lead_snp'  # exactly 9 IVs (one per gene)
```

Do NOT infer the strategy from file names alone.

## Key data files

| File | Description |
|------|-------------|
| `~/Desktop/GTEx_Analysis_v11_eQTL/Lung.v11.eQTLs.signif_pairs.parquet` | Source GTEx eQTL (91 MB, 3.1M rows) |
| `~/Desktop/MR/finngen_R13_ALLERG_RHINITIS.gz` | Outcome 1: Allergic rhinitis (15,958 cases) |
| `~/Desktop/MR/finngen_R13_J10_NASALPOLYP.gz` | Outcome 2: CRSwNP (8,793 cases) |
| `~/Desktop/MR/1kg.v3/EUR/` | LD reference panel (GRCh37 — **build mismatch** with GRCh38 data) |

## 1_correlation_analysis.R outputs

| Output | Path (relative to COR_DIR) | Purpose |
|--------|---------------------------|---------|
| `exposure.pvalue.csv` | `COR_DIR/exposure.pvalue.csv` | Input to 2_ld.R |
| Manhattan plot | `COR_DIR/*.png` | Genome-wide overview (CMplot) |
| Gene scatter plot | `COR_DIR/gene_scatter.png` | x=gene, y=-log10(P), colored by gene |
| Gene summary CSV | `COR_DIR/exposure_gene_summary.csv` | Per-gene: n_SNPs, min_P, lead_SNP, mean_BETA |

## Common mistakes

- ❌ "Exposure only has 9 SNPs" — wrong when IV_STRATEGY is 'all_signif'
- ❌ Merging LD clumping or F-statistic into `0_1data_formating.R` — violates pipeline separation and duplicates downstream scripts
- ❌ Using a GRCh37 LD reference panel (1kg.v3) with GRCh38 data without LiftOver
- ❌ GTEx `phenotype_id` version suffix mismatch — GTEx v11 uses `ENSG00000006128.12` format but config lists `ENSG00000006128`; filter must strip `\.[0-9]+$` before matching
- ❌ `~` in file paths not expanded in `fread(cmd = "gzcat ...")` — shell receives literal `~`; use `path.expand()` in config definitions
- ✅ Each script has one responsibility; respect the serial pipeline
- ✅ The gene summary CSV flags genes with zero significant cis-eQTLs (e.g., GRP with n_SNPs=0)
- ✅ Run order: `source("#_env_initi.R")` → `source("0_0folder_initi.R")` → `source("0_1data_formating_zhou_v2.R")`. The env init script ends with `rm(list=ls())`, so the config must be sourced AFTER it.

## Completed MR results (pooled analysis)

| Outcome | Method | OR | 95% CI | P | Sensitivity |
|---------|--------|-----|--------|---|-------------|
| AR | IVW | 1.023 | 1.000–1.046 | **0.047** | ✅ Q=0.57, Egger=0.78, PRESSO=no outliers |
| CRSwNP | IVW | 1.004 | 0.975–1.035 | 0.78 | — |

**Bonferroni**: α=0.025 → AR P=0.047 **not significant**. Write as "提示性关联" (suggestive).

## Run order reminders

- **New session recovery**: `source("~/Desktop/MR/0_0folder_initi.R")` + `readLines("~/Desktop/MR/progress.md")`
- **env init → config → analysis**: The env init script (`#_env_initi.R`) ends with `rm(list=ls())`, so config must be sourced AFTER it.
- **Gene-specific MR next**: `source("~/Desktop/MR/5_do_MR1_gene_specific.R")` — before running this, check that confounder-removed exposure data exists in MR_CON_DIR.

## Post-MR planned phases

1. 🔴 **Gene-specific MR** — script ready, produce `gene_MR_summary.csv`
2. 🟡 **Skin eQTL + AD** — download GTEx Skin parquet + FinnGen AD; replicate pipeline
3. 🟢 **COLOC + IgE/eosinophil** — co-localization + positive control biomarkers

## Three-phase upgrade plan

- Phase 1 (gene-specific MR): script exists at `~/Desktop/MR/5_do_MR1_gene_specific.R`
- Phase 2 (tissue validation): GTEx Skin eQTL → AD (atopic dermatitis); FinnGen R13 has AD endpoint
- Phase 3 (COLOC + positive control): requires coloc R package; UK Biobank IgE/eosinophil via IEU OpenGWAS

- **Script**: 2_linkage_disequilibrium_analysis_annotated.R
- **Technique**: `--set-all-var-ids @:#` renames reference panel variant IDs to `chr:pos`
- **Input format**: TSV with `SNP` and `P` columns, where SNP = `chr:pos` (e.g., `7:97270465`)
- **Output**: `.clumps` file — one row per clump, `ID` column = index SNP
- **Build mismatch**: 1kg.v3 from IEU is GRCh37. GTEx v11 and FinnGen R13 are GRCh38. Detection via `rs429358` position in .bim.
