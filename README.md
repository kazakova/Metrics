![Logo](/QRePS_logo.png)

## Overview
**QRePS** is a tool for shotgun proteomics that performs statistical analysis on NSAF values in the results of proteome analysis.
QRePS visualizes results of statistical testing with volcano plot and selects differentially regulated proteins (DRP) with different methods (listed below). 
Based on selected features QRePS calculates set of proteomic metrics and performs GO terms enrichment analysis with the use of [STRING](https://string-db.org).
New version of QRePS calculates proteomic metrics and performs GO analysis for [DirectMS1Quant](https://github.com/markmipt/ms1searchpy) results.

### DRP selection

QRePS proviDes three methods to select DRP:

1. *static* - fold change and fdr thresholds are given by user
2. *semi-dynamic* - fold change threshold is given by user and fdr threshold is calculated according to outliers rule: Q3 + 1.5 IQR
3. *dynamic* - lower and upper fold change thresholds are calculated as Q1 - 1.5 IQR and Q3 + 1.5 IQR respectively, fdr threhold is calculated as in *semi-dynamic*

### Metrics calculation
QRePS calculates following metrics

<img src="https://latex.codecogs.com/png.image?\dpi{110}&space;\bg_white&space;\\&space;E=\sqrt{\left(\overline{log_{2}FC}\right)^{2}&plus;\left(\overline{-log_{10}FDR}\right)^{2}}&space;\\&space;E_{m}=\sqrt{\left(\overline{log_{2}FC}-T_{FC}\right)^{2}&plus;\left(\overline{-log_{10}FDR}-T_{FDR}\right)^{2}}&space;\\&space;\pi_{1}&space;=&space;\sum_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;\left(-log_{10}FDR_{i}\right)\right|&space;\\&space;\pi_{2}&space;=&space;log_{10}&space;\left(\prod_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;(-log_{10}FDR_{i})\right|&space;\right)\\\\&space;T_{FC},&space;T_{FDR}&space;\&space;stand&space;\&space;for&space;\&space;fdr&space;\&space;and&space;\&space;fc&space;\&space;thresholds&space;\&space;respectively" title="\bg_white \\ E=\sqrt{\left(\overline{log_{2}FC}\right)^{2}&plus;\left(\overline{-log_{10}FDR}\right)^{2}} \\ E_{m}=\sqrt{\left(\overline{log_{2}FC}-T_{FC}\right)^{2}&plus;\left(\overline{-log_{10}FDR}-T_{FDR}\right)^{2}} \\ \pi_{1}&space;=&space;\sum_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;\left(-log_{10}FDR_{i}\right)\right| \\ \pi_{2}&space;=&space;log_{10}&space;\left(\prod_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;(-log_{10}FDR_{i})\right|&space;\right)\\\\ T_{FC},&space;T_{FDR} \ stand \ for \ fdr \ and \ fc \ thresholds \ respectively" />


## Installation
Install from PyPI
```
pip install QRePS
```

Alternatively, you can install directly from GitHub::
```
pip install git+https://github.com/kazakova/Metrics
```
## Usage
```
qreps [-h]
             (--sample-file SAMPLE_FILE | --quantitation-file QUANTITATION_FILE | --ms1-file MS1_FILE)
             [--pattern PATTERN] [--labels LABELS [LABELS ...]]
             [--input-dir INPUT_DIR] [--output-dir OUTPUT_DIR]
             [--imputation {kNN,MinDet}] [--max-mv MAX_MV]
             [--thresholds {static,semi-dynamic,dynamic,ms1}]
             [--regulation {UP,DOWN,all}] [--species SPECIES]
             [--fold-change FOLD_CHANGE] [--alpha ALPHA]
             [--fasta-size FASTA_SIZE] [--report REPORT]

options:
  -h, --help            show this help message and exit
  --sample-file SAMPLE_FILE
                        Path to sample file.
  --quantitation-file QUANTITATION_FILE
                        Path to quantitative analysis results file.
  --ms1-file MS1_FILE   Path to DirectMS1Quant results file.
  --pattern PATTERN     Input files common endpattern. Default is "_protein_groups.tsv".
  --labels LABELS [LABELS ...]
                        Groups to compare.
  --input-dir INPUT_DIR
  --output-dir OUTPUT_DIR
                        Directory to store the results. Default value is current directory.
  --imputation {kNN,MinDet}
                        Missing value imputation method.
  --max-mv MAX_MV       Maximum ratio of missing values. Default value is 0.5.
  --thresholds {static,semi-dynamic,dynamic,ms1}
                        DE thresholds method.
  --regulation {UP,DOWN,all}
                        Target group of DE proteins.
  --species SPECIES     NCBI species identifier. Default value is 9606 (H. sapiens).
  --fold-change FOLD_CHANGE
                        Fold change threshold.
  --alpha ALPHA         False discovery rate threshold.
  --fasta-size FASTA_SIZE
                        Number of proteins in database for enrichment calculation.
  --report REPORT       Generate report.txt file, default False.
  ```
### Input files
QRePS can be used in three different ways:
1. Perform quantitative analysis (--input-dir, --pattern, --imputation, --sample-file parameters)
2. Use external quantitative analysis results (--quantitation-file)
3. Use results of MS1-based quantitative analysis (--ms1-file) 

Input files for **quantitative analysis** should contain following columns: 
1. 'dbname' (i.e. *sp|P14866|HNRPL_HUMAN*) 
2. 'description' (i.e. *Heterogeneous nuclear ribonucleoprotein L OS=Homo sapiens OX=9606 GN=HNRNPL PE=1 SV=2*) 
3. 'NSAF'
We suggest using [Scavager](https://github.com/markmipt/scavager) *protein_groups* result files. If you use something else, you should specify what files are to be taken from *--input-dir* with common endpattern *--pattern*.

**Quantitation file** should contain 'log2(fold_change)', '-log10(fdr_BH)', 'Gene', 'Protein' columns

**MS1 file** should be *quant_full.tsv file from [DirectMS1Quant](https://github.com/markmipt/ms1searchpy) results.

### Sample file
QRePS tool needs a **sample** file and at least one **data** file for each of groups to perform quantitative analysis.
Sample file should be comma-separated and contain columns 'File Name' and 'SampleID'. 

Input directory can be given either with *--input_dir* or with 'File Name' in sample file.
If both *--input-dir* and path with sample file are given, directory given with *--input-dir* will be used.

Pattern may or may not be included in 'File Name' (see example).
  
SampleID contain labels of groups to be compared and should match those given by *--labels*.
 
### Output files
QRePS produces the following files:
1. volcano plot (volcano.png)
2. missing value ration distribution plot (NaN_distribution.png) (*only if quantitative analysis is performed*)
3. summary table with the results of statistical testing (Quant_res.tsv)
4. summary table with differentially regulated genens (DRF.tsv)
5. symmary table with calculated proteomic metrics (metrics.tsv)
6. summary table with the results of GO terms enrichment analysis (GO_res.tsv)
7. STRING network plot (GO_network.svg)
8. report file (report.txt *if --report True*)

## Example
Input and output files can be found in /example
1. Quantiative analysis
```
qreps --sample-file example_1/a172_dbtrg_sample.csv --labels DBTRG_I,DBTRG_K A172_I,A172_K --input-dir example_1 --output-dir example_1 --imputation kNN --thresholds dynamic --regulation UP 
```

2. External quantitative analysis results
```
qreps --quantitation-file example_2/ms1diffacto_out_DE_A2780_0.5_sum_each_run.txt --labels Chemprot_0.5,Chemprot_K --output-dir example_2 --thresholds semi-dynamic --fold-change 1.5 --regulation all --report True
```
3. DirectMS1Quant resuls analysis
```
qreps --ms1-file ms1quant_out_DE_A2780_5_DE_A2780_K1_quant_full.tsv --labels Chemprot_0.5,Chemprot_K --output-dir output --thresholds ms1 --regulation all --report True
```

## Extra Materials for Publication
Apart from the **QRePS** tool, this repository contains additional resources referenced in the article ["PROTEOMICS-BASED SCORING OF CELLULAR RESPONSE TO STIMULI FOR IMPROVED CHARACTERIZATION OF SIGNALING PATHWAY ACTIVITY"](https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/pmic.202200275):

• [Supplementary Tables](https://github.com/kazakova/Metrics/tree/main/Supplementary_materials)

• [Jupyter Notebooks](https://github.com/kazakova/Metrics/tree/main/Notebooks) with original calculations (you can just use QRePS on your data now) 
