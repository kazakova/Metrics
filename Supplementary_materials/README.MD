# Supplemetary tables

Table_S1a. Results of statistical analysis of IFN-induced proteome changes: imputation with minimal NSAF  [DOI: 10.1021/pr060161n].

Table_S1b. Results of statistical analysis of IFN-induced proteome changes: imputation with k-nearest neighbors machine learning.

Table_S2. Summary of virus assay-based and proteomics scores calculated for ranking the functionality of IFN-dependent antiviral
  mechanisms in cancer and normal cells. Table S2 contains descriptions of proteomic datasets used in this study. 

Table_S3. Results of GO analyses for differentially regulated proteins selected using different workflows, sorted by GO_score: 
  a) imputation with minimal NSAF and satisfying fdr < 0.05 and |log2FC| > 0.585; 
  b) imputation with k-nearest neighbors machine learning and satisfying fdr < 0.05 and |log2FC| > 0.585; 
  c) imputation by minimal detected NSAFs with quartile-based dynamic selection; 
  d) imputation with k-nearest neighbors machine learning and quartile-based dynamic selection.

Table_S4. Summary on missing value percentages across datasets considered in this study:
  label-free quantification at protein level with NSAF [DOI: 10.1021/pr060161n].

Table_S5. Results of Shapiro-Wilk testing for normal data distribution using two imputation strategies: 
  a) the minimal detected NSAFs, b) the kNN-assisted imputation.

# Data processing pipeline:
1. Download raw files from datasets PXD022836, PXD022868, PXD006442 (MS2) and PXD035055 (MS1)
2. Convert raw files to mzML with *msConvert*

        input: raw files (step 1)
        
        output: mzML files
        
3. Perform peptide spectrum matching of LC-MS/MS data with search engines *IdentiPy*, *MSFragger*, or *X!Tandem*. For LC-MS1 data *ms1searchpy* can be used

        input: mzML file (step 2)
        
        output: pepxml file
        
4. Perform search post-processing and calculation of NSAF values with *Scavager*$

        input: pepxml file (step 3)
        
        output: PSMs, peptides, proteins and protein_groups files
        
5. Perform quantitative analysis with *diffacto* or with *stat_test.ipynb* from https://github.com/kazakova/Metrics/Notebooks

        input: Scavager output files for Diffacto or protein_groups file for stat_test (step 4)
        
        ouput: tab-separated file with quantitative analysis results (for diffacto fold_change and fdr_BH should
        be calculated manually):
            NSAF_t-test_*.tsv
            
6. Calculate metrics with *metrics.ipynb* from https://github.com/kazakova/Metrics/Notebooks

        input: tab-separated file(-s) with quantitative analysis results (step 5)
        
        output: tab-separated file and scatterplots with calculated metrics for all the files provided:
            metrics_*.tsv
            
7. Select DE expressed proteins (*de_gene_list* from https://github.com/kazakova/Metrics/Notebooks/GO_terms_analysis.ipynb)

        input: tab-separated file with quantitative analysis results (step 5)
        
        output: 3 comma-separated files with positively and negatively regulated DE proteins and their aggregate (UP, DOWN,
        DE):
            genes_*_UP.csv, genes_*_DOWN.csv, genes_*_DE.csv
            
8. Perform GO terms enrichment analysis (*show_string_picture* and *load_go_enrichment* from
    https://github.com/kazakova/Metrics/Notebooks/GO_terms_analysis.ipynb)
    
        input: comma-separated file with DE proteins (step 7)
        
        output: tab-separated file with identified GO terms and protein network plot
        
9. Select processes from GO analysis results and calculate Enrichment and GO score (formulas can be found below;
    *enrichment_calculation* from https://github.com/kazakova/Metrics/Notebooks/GO_terms_analysis.ipynb)
    
        input: tap-separated file with identified GO terms (step 8)
        
        output: tab-separated file with identified processes, their enrichments and GO scores (sorted by GO score):
            *.string.enrichment_processes.tsv
            
10. Сombine topN processes for different data in one dataframe (*top_processes* from
    https://github.com/kazakova/Metrics/blob/main/GO_terms_analysis.ipynb)
    
        input: tab-separated files with identified processes, their enrichments and GO scores (sorted by GO score) (step 9)
        
        ouput: tab-separated file with processes and their GO scores calculated on different data:
            topN_processes_*.tsv
