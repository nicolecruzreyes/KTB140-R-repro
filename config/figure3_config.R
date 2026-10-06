# Figure 3 / S6-S7 settings (shared settings are in analysis_config.R).
figure3_config <- list(knots=c(2.7,4,9),curve_years=seq(2,10,length.out=161L),
  composition=c('er_pooled','p63_pooled','ki67_pooled'),
  equal_lobule=c('er_mean_lobule','p63_mean_lobule','ki67_mean_lobule'),
  denominator=c('er_p63neg_pooled','ki67_p63neg_pooled'),
  neighborhood_dsp=c('erp63_mean','DSP_immune','DSP_myeloid','DSP_ecm'),
  spline=c('er_pooled','erp63_mean','DSP_ecm'),
  interactions=c('er_fraction','p63_fraction'),
  figure_mm=c(183,212),panel_d=c(30.5,96,122,114),S6_mm=c(183,193),S7_mm=c(183,200))
