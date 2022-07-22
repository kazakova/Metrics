# QMetrics

## Overview
**QMetrics** is a tool for shotgun proteomics that performs statistical analysis on NSAF values in the results of proteome analysis. 
QMetrics visualizes results of statistical testing with volcano plot and selects differentially regulated proteins (DRP) with different methods (listed below). 
Based on selected features QMetrics calculates set of proteomic metrics and performs GO terms enrichment analysis with the use of https://string-db.org.

### DRP selection

QMetrics provides three methods to select DRP:

1. *static* - fold change and fdr thresholds are given by user
2. *semi-dynamic* - fold change threshold is given by user and fdr threshold is calculated according to outliers rule: Q3 + 1.5 IQR
3. *dynamic* - lower and upper fold change thresholds are calculated as Q1 - 1.5 IQR and Q3 + 1.5 IQR respectively, fdr threhold is calculated as in *semi-dynamic*

### Metrics calculation
QMetrics calculated following metrics

<img src="https://latex.codecogs.com/svg.image?E&space;=&space;\sqrt{(\overline{log_{2}FC})^{2}&space;&plus;&space;(\overline{-log_{10}FDR})^{2}}" title="E = \sqrt{(\overline{log_{2}FC})^{2} + (\overline{-log_{10}FDR})^{2}}" />
<img src="https://latex.codecogs.com/svg.image?E_{m}&space;=&space;\sqrt{(\overline{log_{2}FC}&space;-&space;T_{FC})^{2}&space;&plus;&space;(\overline{-log_{10}FDR}&space;-&space;T_{FDR})^{2}}" title="E_{m} = \sqrt{(\overline{log_{2}FC} - T_{FC})^{2} + (\overline{-log_{10}FDR} - T_{FDR})^{2}}" />
<img src="https://latex.codecogs.com/svg.image?\pi_{1}&space;=&space;\sum_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;(-log_{10}FDR_{i})&space;\right|" title="\pi_{1} = \sum_{i = 1} ^{n} \left| log_{2}FC_{i} \cdot (-log_{10}FDR_{i}) \right|" />
<img src="https://latex.codecogs.com/svg.image?\pi_{2}&space;=&space;log_{10}&space;\left(&space;\prod_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;(-log_{10}FDR_{i})\right|&space;\right)" title="\pi_{2} = log_{10} \left( \prod_{i = 1} ^{n} \left| log_{2}FC_{i} \cdot (-log_{10}FDR_{i})\right| \right)" />

<img src="https://latex.codecogs.com/svg.image?T_{FC},&space;T_{FDR}" title="T_{FC},\  T_{FDR}" /> stand for fdr and fc thresholds respectively

## Installation
Download from Github repository https://github.com/kazakova/Metrics. In the directory containing setup.py file run the following command:

```

pip install .
```

## Usage
```
QMetrics [-h] [--sample-file SAMPLE_FILE] [--labels LABELS [LABELS ...]] [--input-dir INPUT_DIR] 
[--output-dir OUTPUT_DIR] [--imputtation {kNN,MinDet}] [--thresholds {static,semi-dynamic,dynamic}]
[--regulation {UP,DOWN,all}] [--species SPECIES] [--fold-change FOLD_CHANGE] [--alpha ALPHA]

options:
  -h, --help            show this help message and exit
  --sample-file SAMPLE_FILE
                        Path to sample file.
  --labels LABELS [LABELS ...]
                        Groups to compare.
  --input-dir INPUT_DIR
  --output-dir OUTPUT_DIR
                        Directory to store the results. Default value is
                        current directory.
  --imputation {kNN,MinDet}
                        Missing value imputation method.
  --thresholds {static,semi-dynamic,dynamic}
                        DRP thresholds selection method.
  --regulation {UP,DOWN,all}
                        Target group of DRP
  --species SPECIES     NCBI species identifier. Default value 9606 (H.
                        sapiens).
  --fold-change FOLD_CHANGE
                        Fold change threshold.
  --alpha ALPHA         False discovery rate threshold.
  ```
### Sample file
The QMetrics tool needs a **sample** file and at least one **data** file for each of groups to compare.
Sample file should be comma-separated and contain columns 'File Name' and 'SampleID'. 

Input directory can be given either with *--input_dir* or with 'File Name' in sample file.
If both *--input-dir* and path with sample file are given, directory given with *--input-dir* will be used. 
  
SampleID contain labels of groups to compare and should match those given by *--labels*.

Sample file example

| File Name                     | SampleID |
|-------------------------------|----------|
| LUM00925_VG                   | DBTRG_K  |
| LUM00931_VG                   | DBTRG_I  |
| LUM00968_VG                   | A172_K   |
| LUM00973_VG                   | A172_I   |
 
### Output files
QMetrics produces the following files:
1. volcano plot (volcano.png)
2. missing value ration distribution plot (NaN_distribution.png)
3. summary table with the results of statistical testing (Quant_res.tsv)
4. summary table with the results of GO terms enrichment analysis (GO_res.tsv)
5. STRING network plot (GO_network.png) 



