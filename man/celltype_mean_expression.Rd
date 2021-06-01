\name{celltype_mean_expression}
\alias{celltype_mean_expression}
\title{
  Obtain the Cell-Type-Specific Mean Expression Levels from the Output of 
  the \code{BUSseq_MCMC} Function
}
\description{
  The function gives the estimated cell-type-specific mean expression
  levels in the output \code{SingleCellExperiment} object.
}
\usage{
  celltype_mean_expression(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object by 
  the function \code{BUSseq_MCMC}.}
}
\value{
  \item{mu.est}{The estimated cell-type-specific mean expression levels,
  a G by K matrix, whose [g,k] element is the mean expression levels of
  cell type k on gene g.}
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
mu.est <- celltype_mean_expression(BUSseqfits_example)
}
