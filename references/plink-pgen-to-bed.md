# PLINK 2.0: Filter 1000 Genomes by population & convert to PLINK 1.9

## Use case

You have 1000 Genomes GRCh38 data in PLINK 2.0 `.pgen/.pvar/.psam` format and need:
1. Filter to a specific super-population (e.g., EUR for LD reference panel)
2. Convert to PLINK 1.9 `.bed/.bim/.fam` format (for `ld_clump_local`, `bigsnpr`, etc.)

## 🔴 Critical pitfall: FID = 0 in 1000G .psam files

**Symptom:** `--keep: 0 samples remaining` or `--keep-if` returns no matches.

**Root cause:** The 1000 Genomes `.psam` file uses `#IID` as its first column header (no `#FID`). PLINK 2.0 responds by setting **FID = 0** for every sample. If your keep file writes FID = IID (e.g., `HG00096  HG00096`), nothing matches — PLINK expects `0  HG00096`.

**Detection:** Run `plink2 --pfile <prefix> --make-just-fam --out tmp` and check the first column:
```
# Correct (FID=0, 1000G-format psam):
0   HG00096   0   0   1   -9

# If you see HG00096 in col 1, the psam has #FID:
HG00096   HG00096   0   0   1   -9
```

**Fix:** In the keep file, use `0` for FID:
```
0   HG00096
0   HG00097
0   HG00099
```

## R script: filter EUR + convert to PLINK 1.9

```r
library(data.table)

work_dir <- "/Users/zhouyumin/Desktop/GRCh38"
setwd(work_dir)

# Step 1: Inspect population breakdown
psam <- fread("all_hg38.psam")
eur_n <- psam[SuperPop == "EUR", .N]
cat(sprintf("EUR samples: %d / %d\n", eur_n, nrow(psam)))
print(psam[SuperPop == "EUR", .N, by = Population])

# Step 2: Write keep file with FID=0
eur_ids <- psam[SuperPop == "EUR", .(FID = "0", IID = `#IID`)]
fwrite(eur_ids, "eur_keep.txt", sep = "\t", col.names = FALSE)

# Step 3: Filter + convert
system2("plink2", c(
  "--pfile",  "all_hg38",
  "--keep",   "eur_keep.txt",
  "--make-bed",
  "--out",    "EUR_hg38_v19"
), stdout = "", stderr = "")

# Step 4: Verify
for (f in c("EUR_hg38_v19.bed", "EUR_hg38_v19.bim", "EUR_hg38_v19.fam")) {
  if (file.exists(f)) {
    cat(sprintf("✅ %s  %s bytes\n", f, format(file.size(f), big.mark = ",")))
  } else {
    cat(sprintf("❌ %s  NOT FOUND\n", f))
  }
}

file.remove("eur_keep.txt")
```

## Typical output

| Population | N |
|------------|---|
| GBR | 91 |
| FIN | 99 |
| IBS | 157 |
| CEU | 179 |
| TSI | 107 |
| **Total** | **633** |

Output: `EUR_hg38_v19.bed/.bim/.fam` — ~633 samples, ~75M variants, PLINK 1.9 binary format.

## Alternative: `--keep-if` (does NOT work for categorical phenotypes)

PLINK 2.0 treats string columns (`SuperPop`, `Population`) as **categorical phenotypes**, not text filters. `--keep-if "SuperPop == 'EUR'"` fails with:

```
Warning: --keep-if categorical phenotype/covariate 'SuperPop' does not have a
category named ''EUR''.
```

Use the `--keep` file approach above — it's reliable and avoids quote-escaping issues.
