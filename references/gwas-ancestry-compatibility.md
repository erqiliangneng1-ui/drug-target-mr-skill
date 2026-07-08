# GWAS Ancestry Compatibility Reference

## GTEx donor demographics

Based on GTEx v8 (v10/v11 proportions are similar):

| Population | Proportion | Notes |
|------------|:----------:|-------|
| White / European | ~84.6% | Primary ancestry driving eQTL estimates |
| African American | ~12.9% | Second largest group |
| Asian American | ~1.3% | Too few for stratified analysis |
| Other / Unknown | ~1.2% | Includes Hispanic/Latino |

**Key implication**: GTEx eQTL effects and LD structure are European-dominated. For MR, outcome GWAS should also be European-ancestry to avoid population stratification.

## GTEx v11 sample counts (Lung tissue)

- GTEx v10 (most recent public release): 943 donors with genotype, 50 tissues with eQTL analysis
- Lung sample size: ~578 samples (one of the best-powered tissues in GTEx)
- GTEx v11 is an incremental release, sample counts similar or slightly larger

## Common MR outcome GWAS by phenotype

### Allergic Rhinitis (AR)

| IEU ID / Source | Ancestry | Cases | Controls | Notes |
|----------------|----------|:-----:|:--------:|-------|
| `ukb-b-7178` | European (UKB) | NA (binary) | ~385K | Self-reported "Hay fever / allergic rhinitis" |
| FinnGen R11 "Allergic rhinitis" | European (Finnish) | 13,846 | 430,927 | ICD-based, 17 genome-wide sig loci |
| Ferreira et al. 2017 | European meta | Various | Various | Classic AR GWAS from meta-analysis |

### Chronic Rhinosinusitis with Nasal Polyps (CRSwNP / Nasal polyps)

| IEU ID / Source | Ancestry | Cases | Controls | Notes |
|----------------|----------|:-----:|:--------:|-------|
| `ebi-a-GCST90018883` | **European** (IEU subset) | 5,093 | 444,966 | Sakaue 2021 Nat Genet; **IEU hosts the pure FinnGen EUR subset**, not the trans-ancestry meta |
| `finn-b-J10_NASALPOLYP` | European (FinnGen R11) | 3,236 | 167,849 | IEU OpenGWAS FinnGen R11 endpoint |
| FinnGen R13 `J10_NASALPOLYP` | European (FinnGen R13) | **8,793** | **369,050** | **Recommended**: download raw TSV from GCS; cases ~2.7× R11 |

**Critical correction on `ebi-a-GCST90018883`**:
The original Sakaue 2021 publication performed a **cross-population meta-analysis** (EUR + EAS). However, the version deposited in **IEU OpenGWAS** (`ebi-a-GCST90018883`) is the **European-only subset**:
- IEU page shows: `population: European`, `ncase: 5093`, `ncontrol: 444966`
- These numbers exactly match the FinnGen European subset (5,093 cases + 444,966 controls)
- Therefore, this ID is ✅ **directly compatible with GTEx** without ancestry-mismatch caveats
- If you need even larger sample size, use **FinnGen R13** directly (8,793 cases)

**FinnGen R13 data access**:
- Manifest: `finngen_R13_manifest.tsv` (maps phenocode → case/control counts → GCS URL)
- Download pattern: `https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_<PHENOCODE>.gz`
- File format: bgzip-compressed TSV with columns: `#chrom`, `pos`, `ref`, `alt`, `rsids`, `nearest_genes`, `pval`, `mlogp`, `beta`, `sebeta`, `af_alt`, `af_alt_cases`, `af_alt_controls`
- Build: **GRCh38** (same as GTEx v11)

## How to check ancestry in GWAS Catalog

URL pattern: `https://www.ebi.ac.uk/gwas/studies/GCST90018883`

Key fields to inspect:
- **Disclosure sample description**: explicitly lists ancestry and counts
- **Discovery ancestry label (country of recruitment)**: tells you which populations
- **Full Summary Statistics**: FTP link — may have population-stratified files

**Important**: GWAS Catalog shows the **original study design** (which may be trans-ancestry). IEU OpenGWAS may host a **population-specific subset**. Always verify both sources. If IEU shows `population: European` and the case count matches the European subset from the original paper, the IEU version is already filtered and safe to use.

## IEU OpenGWAS: VCF download practical notes

When downloading directly from IEU OpenGWAS (presigned URLs):

**File types**:
| File | Purpose | Required? |
|------|---------|:---------:|
| `<id>.vcf.gz` | Imputed GWAS (~1–10M SNPs) | **Yes** |
| `<id>.vcf.gz.tbi` | Tabix index | **Yes** |
| `<id>_raw.vcf.gz` | Directly genotyped only (~500K) | No |
| `<id>_raw.vcf.gz.tbi` | Raw index | No |

**Always use the non-`_raw` version** — lead SNPs are often imputed and absent from raw files.

**Post-download SNP extraction**:
```bash
# Query a specific SNP by position
tabix ukb-b-7178.vcf.gz chr7:97270465-97270465

# Extract multiple SNPs from a BED file
tabix -B ukb-b-7178.vcf.gz snps.bed
```

## FinnGen R13: direct download workflow

When IEU OpenGWAS API is unavailable or sample sizes are too small:

1. Download manifest: `finngen_R13_manifest.tsv`
2. Search for phenocode: `grep NASALPOLYP finngen_R13_manifest.tsv`
3. Columns: `phenocode`, `phenotype`, `category`, `num_cases`, `num_controls`, `path_bucket`, `path_https`
4. Download with `curl -L -o <file> "<path_https>"`
5. The file is a **bgzip-compressed TSV** — use `zcat` or `data.table::fread(cmd = "zcat ...")`
6. **Build is GRCh38** — verify SNP positions match your exposure data

## MR ancestry mismatch: practical impact

Using exposure from European eQTL with outcome from East Asian GWAS introduces:
1. LD structure differences → SNP may not proxy the same gene
2. Allele frequency differences → harmonisation may fail or be unreliable
3. Effect heterogeneity → real biological differences in gene regulation across populations
4. Population stratification → spurious associations

**Exception**: If the variant is a well-validated functional SNP (e.g., missense variant), cross-population MR may be justified with strong biological rationale. Always address in limitations.
