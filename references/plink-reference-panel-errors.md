# PLINK Reference Panel Errors — Common Patterns & Fixes

This documents PLINK errors encountered when using hg38 reference panels for LD clumping/extraction in MR pipelines.

---

## Error 1: Non-standard chromosome names

```
Error: Invalid chromosome code 'PAR1' on line 70692016 of .bim file.
Error: Invalid chromosome code 'chr1_KI270706v1_random' on line 73617919 of .bim file.
```

**Cause**: hg38 panels name pseudoautosomal regions as `PAR1`, `PAR2` and include alternate contigs (`chr1_KI270706v1_random`, `chrUn_*`, etc.). PLINK 1.9 only accepts `1-22`, `X`, `Y`, `MT`.

**Fix**: `--allow-extra-chr` flag. `ieugwasr::ld_clump()` doesn't expose this, so call PLINK directly.

---

## Error 2: Duplicate variant IDs (rsID duplicates)

```
Error: Duplicate ID 'rs1557426849'.
```

**Cause**: The `.bim` file contains duplicate rsIDs (same rsID at different positions). PLINK requires unique IDs.

**Fix**: Create temp copy where `.bim` column 2 is `chr:pos` (naturally unique):

```r
bim <- data.table::fread(paste0(bfile_path, ".bim"), header = FALSE)
bim[[2]] <- paste0(bim[[1]], ":", bim[[4]])  # chr:pos
```

Also improves matching: GTEx SNP IDs use `chr:pos` format.

---

## Error 3: Duplicate variant IDs (chr:pos duplicates within mini panel)

```
Error: Duplicate ID '1:114820381'.
```

**Cause**: Even after converting rsIDs → `chr:pos`, duplicates still exist. Multiple variants (different alleles, indels) can share the same chr:pos in the reference panel. `paste0(chr, ":", bp)` is NOT guaranteed unique.

**Fix**: Use `make.unique()` — adds `.1`, `.2` suffix only to duplicates, leaves first occurrence unchanged (so our input SNPs still match):

```r
bim[[2]] <- make.unique(paste0(bim[[1]], ":", bim[[4]]))  # "1:114820381" + "1:114820381.1"
```

---

## Error 4: `--extract range` file format

```
Error: Line 1 of --extract range file has fewer tokens than expected.
```

**Cause 1 (tab separators)**: PLINK's `--extract range` parser requires **space-separated** columns, not tabs.

**Cause 2 (3 columns)**: PLINK 1.9 expects **4 columns**: `CHR BP1 BP2 LABEL`. Three columns (no label) triggers "fewer tokens".

**Fix**: Write 4 space-separated columns:

```r
# WRONG:
fwrite(regions, range_file, sep = "\t", col.names = FALSE)

# CORRECT:
fwrite(regions %>% mutate(label = paste0("r", row_number())),
       range_file, sep = " ", col.names = FALSE)
```

---

## Error 5: `--clump` input file header format

```
Error: No variant ID field found in clump_input.txt.
```

**Cause**: When calling PLINK directly (not via `ieugwasr`), the `--clump` input file must start with a header line `SNP P` (space-separated, no quotes). `write.table(col.names = TRUE)` writes quoted headers `"rsid" "pval"` that PLINK can't parse, causing "No variant ID field found". `write.table(col.names = FALSE)` omits headers entirely, also causing the same error.

**Fix**: Write header manually with `writeLines()`, then append data:

```r
writeLines("SNP P", tmp_in)
write.table(clump_input[, c("rsid", "pval")], file = tmp_in,
            row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
```

---

## Oversized Panel Performance

`EUR_hg38_v19` has **75M variants** (11GB). Single-threaded PLINK clumping of 357 SNPs takes hours.

**Solutions** (ordered by preference):
1. **Mini-panel extraction**: Extract only relevant regions (±500kb around SNPs) with `--extract range`, then clump on the mini panel (~1.5M variants, seconds).
2. **API clumping**: Use `ieugwasr::ld_clump(bfile = NULL)` for IEU OpenGWAS API.
3. **Smaller panel**: Use standard 1000 Genomes Phase 3 panel (~10M variants).

---

## Mini-panel extraction workflow (complete)

```r
# 1. Read SNP positions
snp_df <- fread("exposure.pvalue.csv", select = c("SNP", "chr.exposure", "pos.exposure"))

# 2. Compute regions (±500kb)
WINDOW <- 500000
regions <- snp_df %>%
  group_by(chr = chr.exposure) %>%
  summarise(min_pos = max(1, min(pos.exposure) - WINDOW),
            max_pos = max(pos.exposure) + WINDOW,
            .groups = "drop")

# 3. Write range file (4 space-separated columns)
fwrite(regions %>% mutate(label = paste0("r", row_number())),
       "ld_mini_ranges.txt", sep = " ", col.names = FALSE)

# 4. Extract
system("plink --bfile EUR_hg38_v19 --extract range ld_mini_ranges.txt --allow-extra-chr --make-bed --out ld_mini/mini")
```

---

## Final Working Solution (all errors combined)

```r
# 1. Write clump input with proper header
writeLines("SNP P", tmp_in)
write.table(clump_input[, c("rsid", "pval")], file = tmp_in,
            row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)

# 2. Fix bim: duplicate IDs → make.unique(chr:pos)
bim <- fread(paste0(mini_bfile, ".bim"), header = FALSE)
bim[[2]] <- make.unique(paste0(bim[[1]], ":", bim[[4]]))
fwrite(bim, file.path(tmp_dir, "ref.bim"), sep = "\t", col.names = FALSE)
file.copy(paste0(mini_bfile, ".bed"), file.path(tmp_dir, "ref.bed"))
file.copy(paste0(mini_bfile, ".fam"), file.path(tmp_dir, "ref.fam"))

# 3. Call PLINK directly
plink_cmd <- paste(
  shQuote(plink_path),
  "--bfile", shQuote(file.path(tmp_dir, "ref")),
  "--clump", shQuote(tmp_in),
  "--clump-kb", CLUMP_KB,
  "--clump-r2", CLUMP_R2,
  "--clump-p1", 0.99,
  "--out", shQuote(tmp_out)
)
system(plink_cmd)

# 4. Parse .clumped output
clumped <- fread(paste0(tmp_out, ".clumped"), header = TRUE, fill = TRUE)
kept_snps <- clumped[!is.na(SNP) & SNP != "" & SNP != "NONE", SNP]
```
