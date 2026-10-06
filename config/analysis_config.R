config <- list(
 version='KTB140_analysis_v1', input_dir='data_common',
 font='Arial', base_size=9, dpi=600,
 neighborhood_min=100L, neighborhood_sensitivity_min=50L,
 ki67_min=10L, radius_um=15, epsilon=1e-6,
 dsp_status='source_reconciled', dsp_file='data_common/dsp_verified.csv',
 dsp_provenance_file='data_common/dsp_provenance.csv', dsp_instrument_column='instrument',
 cluster_ci='normal', cluster_p='residual_t', serum_reference='normal',
 bootstrap_simple=5000L, bootstrap_serial=20000L,
 deterministic_atol=1e-8, deterministic_rtol=1e-7)
state_labels <- c(epithelial_balance='ER–p63 epithelial balance',
 neighborhood='ER–p63 neighborhood enrichment',
 myeloid='Myeloid/antigen-presentation microenvironment',
 stromal='Stromal/fibrovascular state')
