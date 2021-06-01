\name{BIC_BUSseq}
\alias{BIC_BUSseq}
\title{
Obtain BIC from the Output of the \code{BUSseq_MCMC} Function
}
\description{
The function gives the Bayesian Informtion Criterion (BIC) value in the output 
\code{SingleCellExperiment} object.
}
\usage{
BIC_BUSseq(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
\item{BIC_val}{The BIC value.}
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
Example_BIC <- BIC_BUSseq(BUSseqfits_example)
}