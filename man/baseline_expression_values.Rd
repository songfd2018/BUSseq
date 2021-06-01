\name{baseline_expression_values}
\alias{baseline_expression_values}
\title{
Obtain the Log-Scale Baseline Expression Levels from the Output of the 
\code{BUSseq_MCMC} Function
}
\description{
The function gives the estimated log-scale baseline expression levels from the 
output \code{SingleCellExperiment} object.
}
\usage{
baseline_expression_values(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object by 
  the function \code{BUSseq_MCMC}.}
}
\value{
\item{alpha.est}{The estimated log-scale baseline expression levels, a 
G-dimensional vector whose g-th element is the estimated log-scale mean gene 
expression level of gene g in the first cell type.}
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
alpha.est <- baseline_expression_values(BUSseqfits_example)
}
