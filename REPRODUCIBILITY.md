# Reproducibility notes

## Analysis status

The principal within-patient endpoint and its analysis rules were pre-specified. The scoring-method comparison, independent epithelial-marker panel and external within-patient histological-purity model were added after the primary result and are explicitly labelled post hoc.

## Expected generated directories

The complete run creates `01_data_raw/`, additional subdirectories under `02_data_processed/`, `04_results/`, `05_figures/`, `06_manuscript/`, and execution markers under `00_admin/execution_logs/markers/`. Large raw and intermediate matrices are deliberately absent from Git.

## Output-root isolation

`run_all.R` accepts `--output-root <path>` for an isolated rebuild. The path must be inside the project directory or its `_rebuilds` directory. Hash-based success markers prevent `--resume` from silently accepting changed inputs, scripts, configuration or protocol files.

## Generated audit trail

During a complete run, the pipeline creates stage markers with hashes and a number-tracing table that maps manuscript quantities to result files and rows. These generated products are intentionally excluded from the code-only public repository. `MANIFEST.csv` provides checksums for the released source files.
