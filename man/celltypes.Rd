\name{celltypes}
\alias{celltypes}
\title{
Obtain the Cell-type Indicators from the Output of 
the \code{BUSseq_MCMC} Function
}
\description{
The function gives the cell-type indicators of the output
\code{SingleCellExperiment} object.
}
\usage{
celltypes(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object by 
  the function \code{BUSseq_MCMC}.}
}
\value{
\item{w.est}{The estimated cell-type indicators, an N-dimensional vector, 
where N is the total number of cells.}
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
celltypes_est <- celltypes(BUSseqfits_example)
}
