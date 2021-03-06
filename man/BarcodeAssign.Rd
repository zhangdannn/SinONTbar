% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SinONTbar.R
\name{BarcodeAssign}
\alias{BarcodeAssign}
\title{Match the ONT barcodes with Illunima NGS barcodes}
\usage{
BarcodeAssign(ONTfastq, Sinbarcode, MaxMisMatchvalue = 10)
}
\arguments{
\item{ONTfastq}{ONTfastq site, eg. /mnt/data1/zhangdan/Proj/Mouse_iso_atlas/scnapbar/sc_NGScDNA_C57_F_2m_c-kit+/data/nanopore_small.fq}

\item{Sinbarcode}{Singleron barcode site, eg. /mnt/data1/zhangdan/Proj/Mouse_iso_atlas/scnapbar/sc_NGScDNA_C57_F_2m_c-kit+/data/barcodes.tsv.gz}

\item{MaxMisMatch}{barcode assignment MaxMisMatch, we got 57 base in total barcode, and 10 MaxMisMatch is set as default.}
}
\value{
Singleron NGS barcode and ONT read id assignment table
}
\description{
Match the ONT barcodes with Illunima NGS barcodes
}
