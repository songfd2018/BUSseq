\name{location_batch_effects}
\alias{location_batch_effects}
\title{
Obtain the Location Batch Effects from the Output of 
the \code{BUSseq_MCMC} Function
}
\description{
The function gives the estimated location batch effects in the output 
\code{SingleCellExperiment} object.
}
\usage{
location_batch_effects(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
}
\value{
\item{nu.est}{The estimated location batch effects, a G by B matrix,
where [g,b] element is the location batch effect on gene g in the 
batch b compared with the first batch. Note that the first column is 
zero as the first batch is taken as the reference batch without batch 
effects.}
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
nu.est <- location_batch_effects(BUSseqfits_example)
}
