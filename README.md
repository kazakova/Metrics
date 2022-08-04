# QRePS
An open source software to quantify cellular response in proteomics.

## Overview
**QRePS** is a tool for shotgun proteomics that performs statistical analysis on NSAF values in the results of proteome analysis. 
QRePS visualizes results of statistical testing with volcano plot and selects differentially regulated proteins (DRP) with different methods (listed below). 
Based on selected features QRePS calculates set of proteomic metrics and performs GO terms enrichment analysis with the use of [STRING](https://string-db.org).

### DRP selection

QRePS provides three methods to select DRP:

1. *static* - fold change and fdr thresholds are given by user
2. *semi-dynamic* - fold change threshold is given by user and fdr threshold is calculated according to outliers rule: Q3 + 1.5 IQR
3. *dynamic* - lower and upper fold change thresholds are calculated as Q1 - 1.5 IQR and Q3 + 1.5 IQR respectively, fdr threhold is calculated as in *semi-dynamic*

### Metrics calculation
QRePS calculates following metrics

<img src="https://latex.codecogs.com/svg.image?\bg_white&space;E&space;=&space;\sqrt{\left(\overline{log_{2}FC}\right)^{2}&space;&plus;&space;\left(\overline{-log_{10}FDR}\right)^{2}}" title="\bg_white E = \sqrt{\left(\overline{log_{2}FC}\right)^{2} + \left(\overline{-log_{10}FDR}\right)^{2}}" />
<img src="https://latex.codecogs.com/svg.image?\bg_white&space;E_{m}&space;=&space;\sqrt{\left(\overline{log_{2}FC}&space;-&space;T_{FC}\right)^{2}&space;&plus;&space;\left(\overline{-log_{10}FDR}&space;-&space;T_{FDR}\right)^{2}}" title="\bg_white E_{m} = \sqrt{\left(\overline{log_{2}FC} - T_{FC}\right)^{2} + \left(\overline{-log_{10}FDR} - T_{FDR}\right)^{2}}" />
<img src="https://latex.codecogs.com/svg.image?\bg_white&space;\pi_{1}&space;=&space;\sum_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;\left(-log_{10}FDR_{i}\right)\right|" title="\bg_white \pi_{1} = \sum_{i = 1} ^{n} \left| log_{2}FC_{i} \cdot \left(-log_{10}FDR_{i}\right)\right|" />
<img src="https://latex.codecogs.com/svg.image?\bg_white&space;\pi_{2}&space;=&space;log_{10}&space;\left(\prod_{i&space;=&space;1}&space;^{n}&space;\left|&space;log_{2}FC_{i}&space;\cdot&space;(-log_{10}FDR_{i})\right|&space;\right)" title="\bg_white \pi_{2} = log_{10} \left(\prod_{i = 1} ^{n} \left| log_{2}FC_{i} \cdot (-log_{10}FDR_{i})\right| \right)" />

<img src="https://latex.codecogs.com/svg.image?\bg_white&space;T_{FC},&space;T_{FDR}" title="\bg_white T_{FC}, T_{FDR}" /> stand for fdr and fc thresholds respectively

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
qreps [-h] [--sample-file SAMPLE_FILE] [--labels LABELS [LABELS ...]] [--input-dir INPUT_DIR] 
[--output-dir OUTPUT_DIR] [--imputation {kNN,MinDet}] [--thresholds {static,semi-dynamic,dynamic}]
[--regulation {UP,DOWN,all}] [--species SPECIES] [--fold-change FOLD_CHANGE] [--alpha ALPHA]

options:
  -h, --help            show this help message and exit
  --pattern PATTERN     Input files common endpattern. Default value is "_protein_groups.tsv".
  --sample-file SAMPLE_FILE
                        Path to sample file.
  --labels LABELS [LABELS ...] 
                        Groups to compare.
  --input-dir INPUT_DIR
  --output-dir OUTPUT_DIR
                        Directory to store the results. Default value is current directory.
  --imputation {kNN,MinDet}
                        Missing value imputation method.
  --thresholds {static,semi-dynamic,dynamic}
                        DRP thresholds selection method.
  --regulation {UP,DOWN,all}
                        Target group of DRP
  --species SPECIES     
                        NCBI species identifier. Default value is 9606 (H.sapiens).
  --fold-change FOLD_CHANGE
                        Fold change threshold.
  --alpha ALPHA         
                        False discovery rate threshold.
  ```
### Input files
Input file should contain following columns: 
1. 'dbname' (i.e. *sp|P14866|HNRPL_HUMAN*) 
2. 'description' (i.e. *Heterogeneous nuclear ribonucleoprotein L OS=Homo sapiens OX=9606 GN=HNRNPL PE=1 SV=2*) 
3. 'NSAF'

We suggest using [Scavager](https://github.com/markmipt/scavager) *protein_groups* result files. If you use something else, you should specify what files to use with *--pattern*.

### Sample file
QRePS tool needs a **sample** file and at least one **data** file for each of groups to compare.
Sample file should be comma-separated and contain columns 'File Name' and 'SampleID'. 

Input directory can be given either with *--input_dir* or with 'File Name' in sample file.
If both *--input-dir* and path with sample file are given, directory given with *--input-dir* will be used.

Pattern may or may not be included in 'File Name' (see example).
  
SampleID contain labels of groups to compare and should match those given by *--labels*.
 
### Output files
QRePS produces the following files:
1. volcano plot (volcano.png)
2. missing value ration distribution plot (NaN_distribution.png)
3. summary table with the results of statistical testing (Quant_res.tsv)
4. summary table with the results of GO terms enrichment analysis (GO_res.tsv)
5. STRING network plot (GO_network.png) 

## Example
Input and output files can be found in /example

```
qreps --sample-file example/a172_dbtrg_sample.csv --labels DBTRG_I,DBTRG_K A172_I,A172_K --input-dir example --output-dir example --imputation kNN --thresholds dynamic --regulation UP 
```
You will get following command line output 
```
Running DBTRG_I,DBTRG_K

Euclidean distance = 9.574977793545823
Modified euclidean distance = 6.9883518660824535
pi1 = 2660.1372095391735
pi2 = 127.02282548860697

Running A172_I,A172_K

Euclidean distance = 6.189393688035595
Modified euclidean distance = 4.422287959089018
pi1 = 1286.229911664442
pi2 = 82.68135920680737

```

## Extra Materials for Publication
Apart from the **QRePS** tool, this repository contains additional resources referenced in the article "PROTEOMICS-BASED SCORING OF CELLULAR RESPONSE TO STIMULI FOR IMPROVED CHARACTERIZATION OF SIGNALING PATHWAY ACTIVITY" (PROTEOMICS, submitted):

• [Supplementary Tables](https://github.com/kazakova/Metrics/tree/main/Supplementary_materials)

• [Jupyter Notebooks](https://github.com/kazakova/Metrics/tree/main/Notebooks) with original calculations (you can just use QRePS on your data now) 
