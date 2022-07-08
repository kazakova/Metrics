# Metrics
**Supplemetary tables**

TabS1a. Results of statistical analysis of IFN-induced proteome changes: imputation with minimal NSAF  [DOI: 10.1021/pr060161n].

TabS1b. Results of statistical analysis of IFN-induced proteome changes: imputation with k-nearest neighbors machine learning.

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

**Jupyter notebooks**

https://github.com/kazakova/Metrics/blob/main/metrics.ipynb contains function that selects DE proteins based on selected method, calculate metrics on fold change and fdr values of these proteins with the usage example

Example of output plot ![metric_pi2_dynamic_UP+DOWN](https://user-images.githubusercontent.com/107166264/178030413-05b7cfda-6b16-4e82-ab37-81197aee3126.png)

https://github.com/kazakova/Metrics/blob/main/volcano.ipynb contains function that creates volcano plot and counts DE proteins baised on selected methods with the usage example

Example of output plot ![volcano_plants_US_Leaves_US_Roots_US_static](https://user-images.githubusercontent.com/107166264/178030901-37c42ac7-3ebf-419f-b35a-55ae0bd91ee6.png)


https://github.com/kazakova/Metrics/blob/main/GO_terms_analysis.ipynb contains functions for DE proteins selection, GO terms enrichment analysis, GO score calculation and GO processes scatterplot with the usage example 

Processes scatterplot example ![processes_scatterplot](https://user-images.githubusercontent.com/107166264/178031440-332343ea-aedd-4a52-81ea-cc857e8ab191.png)

