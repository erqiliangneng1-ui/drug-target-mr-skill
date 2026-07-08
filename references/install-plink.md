# PLINK Installation Guide (macOS)

PLINK is required for multi-SNP LD clumping in drug-target MR workflows. Two versions are used:

| Command | Version | Use case |
|---------|---------|----------|
| `plink`  | **1.9** (stable) | LD clumping, data management, legacy formats |
| `plink2` | **2.0** (alpha) | Faster I/O, pgen format, newer features |

## Installation on macOS

### Apple Silicon (M1/M2/M3/M4)

```bash
# PLINK 1.9 — single binary, x86_64 (works via Rosetta)
curl -fsSL -o /tmp/plink_mac.zip https://s3.amazonaws.com/plink1-assets/plink_mac_20241022.zip
cd /tmp && unzip -o plink_mac.zip
sudo cp plink /usr/local/bin/plink
sudo chmod +x /usr/local/bin/plink

# PLINK 2.0 — Apple Silicon native (arm64), check for latest alpha version
curl -fsSL -o /tmp/plink2_mac_arm64.zip https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260504.zip
cd /tmp && unzip -o plink2_mac_arm64.zip
sudo cp plink2 /usr/local/bin/plink2
sudo chmod +x /usr/local/bin/plink2
```

> **Finding the latest alpha**: Go to https://www.cog-genomics.org/plink/2.0/, scroll to "Introduction, downloads", and look for the `plink2_mac_arm64_*.zip` link. The alpha version number increments over time (alpha5, alpha6, alpha7...).

### Intel Mac (x86_64)

Replace `arm64` with `amd_avx2` or `x86_64` in the PLINK 2.0 URL. PLINK 1.9 uses the same `plink_mac_20241022.zip`.

## Verification

```bash
plink --version    # PLINK v1.9.0-b.7.7
plink2 --version   # PLINK v2.0.0-a.7.1 M1
```

## LD clumping usage

### PLINK 1.9 (preferred, universally supported in MR tools)
```bash
plink \
  --bfile 1000G_EUR_Phase3 \
  --clump significant_snps.txt \
  --clump-p1 1e-5 \
  --clump-p2 1e-5 \
  --clump-r2 0.001 \
  --clump-kb 10000 \
  --out ld_clump_output
```

### PLINK 2.0
```bash
plink2 \
  --bfile 1000G_EUR_Phase3 \
  --clump significant_snps.txt \
  --clump-p1 1e-5 \
  --clump-p2 1e-5 \
  --clump-r2 0.001 \
  --clump-kb 10000 \
  --out ld_clump_output
```

## Notes

- PLINK is **not available via Homebrew** — must be downloaded from the official S3 buckets.
- PLINK 1.9 is the safer choice for MR workflows since `TwoSampleMR` and `ieugwasr` were tested against it.
- PLINK 2.0 is faster on Apple Silicon (native arm64 binary vs x86_64 through Rosetta).
- `vcf_subset` (from PLINK 2.0 zip) is also installed — a helper for extracting VCF regions.
