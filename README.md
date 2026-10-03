# Tissue composition in gastric peritoneal metastasis

This repository contains the analysis code and reproducibility materials for:

> Chao Huang, Qin Yu, Zheng-Ming Zhu, and Yinfei Kong. **Tissue composition quantitatively reproduces bulk transcriptomic differences between paired primary gastric tumours and solid peritoneal metastases.**

## Purpose

The code tests how much of a primary-versus-peritoneal bulk transcriptomic contrast can be reproduced by tissue composition. It implements patient-paired analyses, independent scoring checks, pseudo-bulk calibration and an external histological-purity analysis. Scientific results are reported in the associated manuscript.

## Repository contents

- `03_scripts/`: analysis, figure, verification and submission-asset scripts.
- `00_admin/protocol_v1.md`: analysis protocol frozen before the primary site contrast.
- `00_admin/protocol_lock.json`: hashes for the frozen scientific configuration and programme definitions.
- `00_admin/data_manifest.csv`: public input sources, roles and checksums.
- `00_admin/deviations.md`: transparent record of deviations and post-hoc additions.
- `renv.lock`: exact R package dependency record.
- `MANIFEST.csv`: SHA-256 hash and byte size for every released file.

Raw data, processed matrices, result tables, figures, manuscripts, submission files and internal review materials are excluded. Public inputs are rebuilt by the pipeline; their accessions, URLs and checksums are recorded in `config.yml` and `00_admin/data_manifest.csv`.

## Reproduce the analysis

The archived run used R 4.5.2 on macOS. From the repository root:

```bash
Rscript -e 'install.packages("renv")'
Rscript -e 'renv::restore(prompt = FALSE)'
Rscript --vanilla run_all.R --stage all
```

The `all` stage downloads the public inputs, prepares bulk and single-cell objects, rebuilds the frozen programmes, runs the primary and secondary analyses, generates the figures and tables, performs verification, and runs the post-hoc scoring and histological-purity checks used in the submitted manuscript.

The `tools` entries in `config.yml` record the macOS paths used for Pandoc and LibreOffice. These tools are reported by preflight and used for manuscript rendering; the statistical analyses do not depend on those exact executable paths.

The full rebuild requires substantial disk space and memory because the public single-cell inputs are large. Individual stages can be run with `--stage`; completed stages can be resumed with `--resume` when their recorded hashes match.

## Reproducibility safeguards

The primary endpoint, programme membership, eligibility rules and scoring procedure were locked before any site contrast. `00_admin/protocol_lock.json` records the lock and hashes. Each stage writes a success marker with input, output, script, configuration and protocol hashes. The completed frozen pipeline passed 22 verification checks, and an isolated rebuild reproduced the audited result files exactly. Analyses added after inspection of the primary result are labelled post hoc in code and outputs.

## Data availability

All analysed data were previously published and publicly deposited. The principal accessions are GSE314812, GSE237876, GSE183904, GSE163558, GSE308231, GSE239676, GSE234129 and GSE228598. See `00_admin/data_manifest.csv` for source roles and checksums.

## Citation and contact

Please cite the associated article and this repository. Machine-readable citation metadata are in `CITATION.cff`.

- Zheng-Ming Zhu: zzm8654@163.com
- Yinfei Kong: yikong@fullerton.edu

Repository: https://github.com/kingterkong/gastric-pm-tissue-composition
