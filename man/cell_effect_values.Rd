\name{cell_effect_values}
\alias{cell_effect_values}
\title{
Obtain the cell-specific size effects from the Output of 
the \code{BUSseq_MCMC} Function
}
\description{
The function gives the estimated cell-specific size effects in the output 
\code{SingleCellExperiment} object.
}
\usage{
cell_effect_values(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
\item{delta.est}{The estimated cell-specific global effects, an N-dimensional 
vector. Note that the first element in each vector is zero as the first cell 
in each batch is taken as the reference cell.}
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
delta.est <- cell_effect_values(BUSseqfits_example)
}
