# LD Clumping Debugging Narrative — 2026-07-06 Session

## Context
- Reference panel: `EUR_hg38_v19` (GRCh38, 75M variants, 11GB, from `~/Desktop/GRCh38/`)
- 357 input SNPs (5 genes: NMU, ADCYAP1, GRP, VIP, NGF), chr:pos format
- Goal: LD clumping with kb=10000, r²=0.001
- PLINK 1.9 via Homebrew (`/usr/local/bin/plink`)

## Error Sequence

### Error 1: `Invalid chromosome code 'PAR1'`
- **Cause**: GRCh38 reference panel uses `PAR1`, `PAR2` for pseudoautosomal regions. PLINK 1.9 only recognizes `1-22, X, Y, MT, XY`.
- **Attempted fix**: Filter `.bim` rows starting with `^PAR`
- **Result**: Failed — missed other non-standard chromosomes

### Error 2: `Invalid chromosome code 'chr1_KI270706v1_random'`
- **Cause**: Same root cause — GRCh38 includes many alternate contig names
- **Fix**: Bypass `ieugwasr::ld_clump()` and call PLINK directly with `--allow-extra-chr`
- **Result**: Fixed non-standard chromosome issue

### Error 3: `Duplicate ID 'rs1557426849'`
- **Cause**: `.bim` column 2 (variant ID) has duplicate rsIDs at different positions
- **Fix**: Replace rsIDs with `chr:pos` format in a temp copy of `.bim`
- **Result**: Fixed duplicate IDs, but...

### Error 4: Big panel clumping extremely slow
- **Cause**: 75M variants × 357 SNPs, single-threaded LD calculation
- **Fix**: Extract mini panel first, then clump on mini panel
- **Result**: 2a_extract_mini_ld_panel.R created

### Error 5: `--extract range` file format rejected
- **Cause**: `--extract range` requires 4 **space-separated** columns with a label. Tab-separated 3-column input triggers "fewer tokens than expected"
- **Fix**: `fwrite(regions %>% mutate(label = paste0("r", row_number())), file, sep = " ")`
- **Result**: Mini panel extracted successfully (148万 variants, 267MB vs 75M/11GB)

### Error 6: `Duplicate ID '1:114820381'` (mini panel)
- **Cause**: Same chr:pos can have multiple variants (different alleles/indels) → chr:pos alone is not unique
- **Fix**: `make.unique(paste0(chr, ":", pos))` — adds .1/.2 suffixes only to duplicates, keeping first occurrence clean for matching our input SNPs
- **Result**: Duplicate IDs resolved

### Error 7: `No variant ID field found in clump_input.txt`
- **Cause**: Direct PLINK `--clump` call requires header `SNP P` (space-separated, no quotes). `write.table(..., col.names = TRUE)` writes `"rsid" "pval"` with quotes
- **Fix**: `writeLines("SNP P", file)` then `write.table(..., append = TRUE, col.names = FALSE)`
- **Result**: ✅ 357 → 11 independent SNPs

## Final Working Pipeline

```
EUR_hg38_v19 (75M variants, 11GB)
    ↓ plink --extract range (4 space-separated cols, --allow-extra-chr)
ld_mini/mini (1.48M variants, 267MB)
    ↓ make.unique(chr:pos) bim fix
    ↓ plink --clump (header "SNP P", space-separated)
11 independent SNPs → exposure.LD.csv
```

## Key Lessons

1. **Never use `ieugwasr::ld_clump()` with GRCh38 panels** — it can't pass `--allow-extra-chr`
2. **Always use `make.unique()` not plain `paste0()` for variant IDs** — even chr:pos can have duplicates
3. **Mini-panel strategy**: extract regions first, clump on mini = seconds vs hours
4. **`--extract range`**: MUST be 4 space-separated columns with a label
5. **`--clump`**: MUST start with header `SNP P` (no quotes, space-separated)
6. **`fwrite()` default separator is tab** — must explicitly set `sep = " "` for PLINK range files
7. **API clumping (`bfile=NULL`)** is a viable fallback when local panel issues become too complex
