## GO_terms_analysis.ipynb

*de_gene_list* selects differentially regulated proteins (DRP) based on selected method (for more detailed description see GO_terms_analysis.ipynb)

*show_string_picture* downloads GO terms enrichment results network from https://string-db.org

*load_go_entrichment* downloads tables with GO terms enrichment analysis results from https://string-db.org

*enrichment_calculation* calculates enrichment and GO_score (formulas can be found at GO_terms_analysis.ipynb)

*top_processes* selects top N (default 10) processes from *load_go_enrichment* output files based on GO_score

*processes_scatterplot* creates scatterplot based on *top_processes* output 

## metrics.ipynb
*metrics* selects differentially regulated proteins (DRP) based on selected method and calculate proteomic metrics (for more detailed description see metrics.ipynb)

## stat_test_kNN.ipynb
*stat_test_kNN* performs missing value imputation with KNNImputer and performs statistical testing

## stat_test_MinDet.ipynb
*stat_test_MinDet* performs missing value imputation with minimum detected NSAF value and performs statistical testing

## volcano.ipynb
*volcano* creates volcano plot based on *stat_test* results
