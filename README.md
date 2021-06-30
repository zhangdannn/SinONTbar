# SinONTbar
A R package designed for singleron and ONT data barcode assignment

# Requirements
B2B (https://github.com/tangchao7498/CBUC) \
Biostrings \
IRanges \
rDNAse \
data.table \
ShortRead \
DECIPHER

# Installation
Your can install B2B from GitHub by:

```
library(devtools)
devtools::install_github("zhangdannn/SinONTabr")
```

# Run example
```
library(SinONTbar)
library(B2B)
library(data.table)
library(ShortRead)
library(Biostrings)
library(DECIPHER)

table = SinONTbar::BarcodeAssign(ONTfastq = system.file("data", "nanopore_small.fq", package = "SinONTbar"),
              Sinbarcode = system.file("data", "Singleron.barcodes.tsv.gz", package = "SinONTbar"),
              MaxMisMatchvalue = 10)
```
![image](https://github.com/zhangdannn/SinONTbar/result.png)
