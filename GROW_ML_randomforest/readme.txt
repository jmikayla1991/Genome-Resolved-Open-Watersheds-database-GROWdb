Random forest analysis for GROW paper
Ikaia Leleiwi, May 31st 2023

Directory Structure
wd = GROW_ML_randomforest
├── data
│   ├── SOMfiles
│   │   ├── GROW_SuppFile1.xlsx
│   │   └── GROW_SuppFile2.xlsx
│   ├── SOMfiles.zip
│   ├── metag
│   │   ├── GROWdb_with_vars_20220715_115_metaG.csv
│   │   ├── GROWdb_with_vars_20220715_142_metaG.csv
│   │   ├── GROWdb_with_vars_20220715_50_metaG.csv
│   │   ├── strict_mapping_table_95id_3x_60%cov_relabundance_115.csv
│   │   ├── strict_mapping_table_95id_3x_60%cov_relabundance_142.csv
│   │   └── strict_mapping_table_95id_3x_60%cov_relabundance_50.csv
│   ├── metat
│   │   ├── corrplot_mergedtable_43.csv
│   │   ├── corrplot_mergedtable_52.csv
│   │   ├── geTMM_norm.counts.rpk_edger_genes.csv
│   │   ├── geTMM_norm.counts.rpk_edger_genome_relabund_43.csv
│   │   └── geTMM_norm.counts.rpk_edger_genome_relabund_52.csv
│   ├── metat_metag.csv
│   ├── ml_input_scaled.csv
│   ├── ml_input_taxa.csv
│   ├── rawfiles
│   │   ├── geTMM_norm.counts.rpk_edger_genes.csv
│   │   ├── geTMM_norm.counts.rpk_edger_genome_relabund.csv
│   │   └── strict_mapping_table_95id_3x_60%cov_relabundance.csv
│   ├── rf_taxa_byClass_df.csv
│   ├── rf_taxa_overall_stats_df.csv
│   └── rf_taxa_varimp_df.csv
├── figures  #####only produced figures for successful models
│   ├── rf_MAGrichness_varimp.svg
│   ├── rf_generichness_varimp.svg
│   ├── rf_taxa_model_accuracy_by_most_imp_var.svg
│   └── rf_taxa_varimp_heatmap.svg
├── readme.txt
├── references
│   ├── Figure4_draft1.pdf
│   ├── ForKai.zip
│   ├── StreamCat_variable_info.csv
│   └── lakecatvariablelist-quickreference.xlsx
└── scripts
    ├── GROW_ML_explore_data.R #####initial data exploration
    ├── GROW_ML_input_data.R  #####produces machine learning input dataframes
    ├── GROW_ML_randomforest_diversity.R  #####runs random forest regression to predict diversity metrics from landuse
    ├── GROW_ML_randomforest_fctcat.R  #####runs random forest classification to predict function class from landuse
    ├── GROW_ML_randomforest_taxa.R  #####runs random forest classification to predict Class (taxonomy) from landuse
    ├── GROW_ML_randomforest_taxa_figs.R  #####produces basic figures from bootstrapped Class (taxonomy) classification models
    ├── GRWO_ML_randomforest_fctreg.R  #####runs random forest regression with each function class as target
    └── GRWO_ML_randomforest_landusereg.R  #####runs random forest regression and classification to predict landuse from function class values


#notes
The only successful predictions using landuse were made for Class taxonmy (GROW_ML_randomforest_taxa.R) and either gene or MAG richness (GROW_ML_randomforest_diversity.R)
