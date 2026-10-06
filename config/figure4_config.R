# Figure 4 / S8-S9 settings (shared settings are in analysis_config.R).
figure4_config <- list(
  serum_covariates=c('Age','BMI','ParousFlag'),
  sensitivity_covariates=c('Age','BMI'),
  lobule_covariates=c('Age','BMI'),
  states=c('epithelial_balance','neighborhood','myeloid','stromal'),
  targeted_features=c('CTAG1B.1','CTSD','FOXP3'),
  n_common=109L,n_serum_parity=131L,n_postpartum=43L,n_postpartum_training=40L,
  features=209L,cv_v=10L,cv_comparison_repeats=50L,cv_fixed_repeats=100L,
  cv_outer_repeats=50L,cv_inner_v=5L,
  folds_dir='resampling/figure4',curve_points=161L,
  endpoints=c('p63_fraction','ER_fraction','ERp63_enrichment',
    'all_lobules_p63_high','all_lobules_ER_low','all_lobules_ERp63_high','all_lobules_composite',
    'count_rich_p63_high','count_rich_ER_low','count_rich_ERp63_high','count_rich_composite',
    'coupling_er','coupling_p63')
)
