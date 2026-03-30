# fossilized-ark

Benchmarking ensilicated (silica-cached) vs frozen (-80C) DNA preservation for Oxford Nanopore whole-genome sequencing using the GIAB Ashkenazi trio (HG002, HG003, HG004).

**Manuscript**: [medRxiv 2025.10.26.25338579v1](https://www.medrxiv.org/content/10.1101/2025.10.26.25338579v1)

## Overview

This repository contains analysis scripts, summary results, and figures for a systematic comparison of ensilicated and frozen DNA preservation methods. Six whole-genome ONT libraries (3 frozen, 3 ensilicated) were sequenced, aligned to GRCh38, and benchmarked across multiple axes:

- **Read quality**: read length distributions, Q-scores, per-read identity, mismatch rates
- **Small variant calling**: DeepVariant SNV/indel calls benchmarked against GIAB v4.2.1 truth sets with hap.py
- **Structural variant calling**: Sniffles2 SV calls benchmarked against GIAB Tier 1 SV v0.6 with Truvari
- **Mendelian concordance**: all 8 frozen/ensilicated trio combinations tested for inheritance consistency
- **CpG methylation**: modkit pileup calls compared between conditions and against EM-seq ground truth (Pearson correlation, discordance analysis, binomial null model, single-read calibration, pairwise AUROC)
- **Titration experiment**: 7 frozen/ensilicated mixture ratios (0-100%) to test variant calling robustness

## Samples

| Sample | Role | Frozen ID | Ensilicated ID |
|--------|------|-----------|----------------|
| HG002 | Son (proband) | CD_3033_GIAB | CD_3032_Cache |
| HG003 | Father | CD_3031_GIAB | CD_3030_Cache |
| HG004 | Mother | CD_3029_GIAB | CD_3028_Cache |

Raw sequencing data are available at NCBI SRA (accession pending).

## Repository structure

```
scripts/          analysis and pipeline scripts
  *.sh            shell pipeline (alignment, variant calling, methylation extraction)
  *.py            python analysis and plotting scripts
results/          summary tables
  error_profiles/ mismatch rates, per-read identity histograms
  methylation_discordance/  discordance, calibration, AUROC, binomial null
  mendelian/      Mendelian concordance summaries
```

## Pipeline

### Alignment and QC
Reads were basecalled with Dorado v1.4.0 and aligned to GRCh38 (GIABv3 analysis set) with minimap2 v2.28 (`-ax map-ont -K5g`). Alignments were sorted and indexed with samtools v1.19.2/v1.21.

### Variant calling
- **Small variants**: DeepVariant v1.6.1 (ONT_R104 model, GPU, 32 shards)
- **Structural variants**: Sniffles2 v2.2

### Benchmarking
- **Small variants**: hap.py v0.3.12 with vcfeval engine against GIAB v4.2.1 truth sets
- **Structural variants**: Truvari v4.1.0 against GIAB Tier 1 SV v0.6
- **Additional benchmarks** (HG002): CMRG small variant v1.00, CMRG SV v1.00, tandem repeat v1.0

### Methylation
5mC calls were extracted with modkit v0.4.1 (`--combine-strands --cpg --bedgraph`). Analyses include:
- Pairwise Pearson correlation matrix (9 datasets: 3 samples x {frozen, ensilicated, EM-seq})
- Discordance analysis with binomial null model
- Single-read concordance against EM-seq (accuracy, sensitivity, specificity, ROC)
- Calibration plots and pairwise AUROC heatmap
- Coverage-stratified discordance
- CpG island context and GC content characterization

## Software versions

| Tool | Version |
|------|---------|
| Dorado | 1.4.0 |
| minimap2 | 2.28 |
| samtools | 1.19.2 / 1.21 |
| DeepVariant | 1.6.1 |
| Sniffles2 | 2.2 |
| modkit | 0.4.1 |
| hap.py | 0.3.12 |
| Truvari | 4.1.0 |
| bcftools | 1.21 |
| bedtools | 2.31.1 |
| Python | 3.9.19 |
| pysam | 0.21.0 |
| numpy | 1.26.4 |
| pandas | 2.2.2 |
| matplotlib | 3.8.4 |
| scipy | 1.13.1 |
| scikit-learn | 1.5.1 |
| seaborn | 0.13.2 |

## Reference data

- **Genome**: GRCh38 GIABv3 no-alt analysis set with masked GRC decoys (MAP2K3/KMT2C/KCNJ18)
- **Small variant benchmark**: GIAB v4.2.1 (HG002, HG003, HG004)
- **SV benchmark**: GIAB Tier 1 v0.6 (HG002)
- **CpG islands**: UCSC cpgIslandExt (GRCh38)
- **EM-seq truth**: ENCODE/GIAB EMSeq LAB01 REP01 (HG002, HG003, HG004)
