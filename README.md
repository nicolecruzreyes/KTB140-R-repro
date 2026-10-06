# KTB140 analysis code

R code for the statistical analyses, figures and supplementary tables in:

> (pending) [AUTHORS]. [TITLE]. [JOURNAL, YEAR]. doi: [ARTICLE DOI]

Archived version: (pending) [ZENODO DOI]

This repository contains code only. No study data are included (see **Data availability**).

## What the code does

| Script | Produces |
|---|---|
| `RUN_MODULE.R --module figure1` | Figure 1 and Figures S1–S4 (lobule-level mixed models, donor bootstrap decompositions) |
| `RUN_MODULE.R --module figure2` | Figure 2 and Figure S5 (GeoMx DSP parity models, tissue states) |
| `RUN_MODULE.R --module figure3` | Figure 3 and Figures S6–S7 (time since last birth, spline and interaction models) |
| `RUN_MODULE.R --module figure4` | Figure 4 and Figures S8–S9 (serum screens, cross-validated prediction) |
| `RUN_TABLE_ONLY.R` | Models reported only in the supplementary tables (S3B–S3E, S4, S7 slopes) |
| `RUN_ALL.R` | All of the above, in order |

Each run writes result tables (CSV and Excel), plot data, fitted model objects, figures
(PDF, SVG, 600-dpi PNG), individual panels and, optionally, PowerPoint assemblies to
`outputs/<run-id>/`.

Code layout:

- `R/core/` – shared helpers: data assembly (`02_data.R`), model fitting and multiple-testing
  correction (`01_numeric.R`), and graphics/export (`03_graphics.R`, `05_portable_exports.R`).
- `R/modules/` – one model file and one graphics file per figure, plus `tableonly_models.R`.
- `config/` – analysis settings (`analysis_config.R`, `figure3_config.R`, `figure4_config.R`).

## Requirements

- R (analyses were run with R 4.5.3 on macOS and R 4.3.3 on Windows)
- Packages: ggplot2 (≥ 3.4.0), systemfonts, ragg, svglite, jsonlite, dplyr, sandwich, lme4,
  lmerTest, openxlsx, Cairo, pdftools, officer, zip, xml2
- The Arial font installed locally (not bundled)

Check for missing packages, then install them if needed:

```
Rscript 00_INSTALL_PACKAGES.R
Rscript 00_INSTALL_PACKAGES.R --yes
```

## Running

Run from the project root with `Rscript` (the scripts read command-line options, so they are
not intended to be run line by line in the RStudio console):

```
Rscript RUN_ALL.R --run-id my_run
Rscript RUN_MODULE.R --module figure2 --run-id my_fig2
Rscript RUN_TABLE_ONLY.R --run-id my_tables
```

On Windows, use the full path to `Rscript.exe`, for example
`& "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" RUN_ALL.R --run-id my_run` in PowerShell.
Existing output folders are never overwritten; use a new `--run-id` for each run.

## Statistical notes

- Lobule-level mixed models (lme4/lmerTest, maximum likelihood, donor random intercept) report
  Satterthwaite degrees of freedom; normal-Wald results are reported separately as a sensitivity.
- Regression models use HC3 or donor-clustered HC1 robust standard errors (sandwich).
- Multiple testing is corrected with Holm or Benjamini–Hochberg within the families stated in the paper.
- The donor bootstrap (Figure 1) and cross-validation (Figure 4) use fixed, pre-generated
  resampling plans, so repeated runs give identical results.

## Data availability

Data can be provided upon reasonable request to corresponding author.

The code expects the following data and supporting files, which are **not** in this repository:

- `data_common/`, `data_cells/` – analysis input tables
- `resampling/` – the fixed bootstrap and cross-validation plans
- `config/visual_jitter.csv`, `config/figure2_visual_jitter.csv`, `config/figure3_visual_jitter.csv`,
  `config/figure4_visual_jitter.csv` – per-donor point offsets for plotting (keyed by donor ID)
- `templates/blank_*.pptx` – blank PowerPoint templates at figure dimensions (PowerPoint export only)

Steps upstream of these inputs (whole-cell reconstruction from the MxIF images) and the conditional
spatial randomization reported in Table S3I were performed separately in Python [and are available at …].

## Citation

See `CITATION.cff`. Please cite both the article and the archived code DOI.

## License

MIT
