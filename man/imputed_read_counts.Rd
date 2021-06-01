\name{imputed_read_counts}
\alias{imputed_read_counts}
\title{
  Obtain the Imputed Read Count Matrix from the Output of the \code{BUSseq_MCMC} 
  Function
}
\description{
  The function gives the imputed read counts in the output 
\code{SingleCellExperiment} object.
}
\usage{
  imputed_read_counts(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
  \item{CountData_imputed}{The imputed read counts, a \code{CountData}
  object with length equal to the batch number. Each element is a read
  count matrix for a specific batch, where each row corresponds to 
  a gene and each column represents a cell.}
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
Example_CountData_imputed <- imputed_read_counts(BUSseqfits_example)
}