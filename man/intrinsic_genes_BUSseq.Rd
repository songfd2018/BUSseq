\name{intrinsic_genes_BUSseq}
\alias{intrinsic_genes_BUSseq}
\title{
Obtain the Intrinsic Gene Indicators from the Output of 
the \code{BUSseq_MCMC} Function
}
\description{
The function gives the estimated intrinsic gene indicators in the output 
\code{SingleCellExperiment} object.
}
\usage{
intrinsic_genes_BUSseq(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
\item{intrinsic_genes}{A vector indicating whether a gene is 
intrinsic or not.}
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
BUSseqfits_example
intri_genes <- intrinsic_genes_BUSseq(BUSseqfits_example)
}