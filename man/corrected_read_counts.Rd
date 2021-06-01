\name{corrected_read_counts}
\alias{corrected_read_counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate the Corrected Read Count Matrix within the Output of 
  the \code{BUSseq_MCMC}
}
\description{
  The function generates a version of count data, for which the batch
  effects are removed and the biological variabilities are retained. 
  We develop a quantile match approach based on the idea of inverse
  transform sampling. The users can perform downstream analysis on the
  corrected read count matrix, such as clustering, differentially
  expressed gene identification and so on, as if all the data were
  measured in a single batch.
}
\usage{
  corrected_read_counts(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
  \item{corrected_data}{The corrected read count matrix, which is added to 
  the assays of \code{sce_BUSseqfit} and can be loaded 
  by \code{assay(sce_BUSseqfit, "corrected_data")}.}
}
\references{
Song, Fangda, Ga Ming Angus Chan, and Yingying Wei. 
Flexible experimental designs for valid single-cell RNA-sequencing 
experiments allowing batch effects correction. 
Nature communications 11, no. 1 (2020): 1-15.
}
\author{
  Fangda Song
}
\examples{
# "BUSseqfits_example" is an example output
library(SingleCellExperiment)
BUSseqfits_example
BUSseqfits_corrected <- corrected_read_counts(BUSseqfits_example)
corrected_data <- assay(BUSseqfits_corrected, "corrected_data")
}