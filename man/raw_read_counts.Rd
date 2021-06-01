\name{raw_read_counts}
\alias{raw_read_counts}
\title{
  Obtain the Raw Read Count Matrix from the Output of the \code{BUSseq_MCMC} 
  Function
}
\description{
  The function gives the original read count matrix. 
  This function is actually equivalent to 
\code{assay(sce_BUSseqfit, "counts")}.
}
\usage{
  raw_read_counts(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
  \item{CountData_raw}{The raw read counts, a matrix  with each row for a gene
  and each column for a cell.}
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
Example_CountData_raw <- raw_read_counts(BUSseqfits_example)
}