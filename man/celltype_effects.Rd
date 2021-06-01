\name{celltype_effects}
\alias{celltype_effects}
\title{
Obtain the Cell-type Effects from the Output of the \code{BUSseq_MCMC} Function
}
\description{
The function gives the estimated cell-type effects in the output 
\code{SingleCellExperiment} object.
}
\usage{
celltype_effects(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object by 
  the function \code{BUSseq_MCMC}.}
}
\value{
\item{beta.est}{The estimated cell-type effects, a G by K matrix, 
whose [g,k] element is the effects of cell type k on gene g 
compared with the first cell type. Note that the first column 
is zero as the first cell type is taken as the baseline cell type.}
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
beta.est <- celltype_effects(BUSseqfits_example)
}
