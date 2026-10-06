# config/

- `analysis_config.R` – shared settings (input locations, font and resolution, bootstrap sizes,
  inference conventions, tissue-state labels)
- `figure3_config.R` – Figure 3 / S6–S7 settings (spline knots, endpoints, curve grid)
- `figure4_config.R` – Figure 4 / S8–S9 settings (covariates, features, cross-validation design)

The per-donor plotting offset files (`visual_jitter.csv`, `figure2_visual_jitter.csv`,
`figure3_visual_jitter.csv`, `figure4_visual_jitter.csv`) must also be placed in this folder to
run the code. They contain donor IDs, are distributed with the data, and are excluded from the
repository by `.gitignore`.
