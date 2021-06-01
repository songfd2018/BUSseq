\name{overdispersions}
\alias{overdispersions}
\title{
Obtain the Overdispersion Parameters from the Output of 
the \code{BUSseq_MCMC} Function
}
\description{
The function gives the estimated overdispersion parameters in the output 
\code{SingleCellExperiment} object.
}
\usage{
overdispersions(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
\item{phi.est}{The estimated overdispersion parameters, a G by B matrix,
where [g,b] element is the overdispersion parameter of gene g in batch b.}
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
phi.est <- overdispersions(BUSseqfits_example)
}