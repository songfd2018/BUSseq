\name{overdispersions}
\alias{overdispersions}
\title{
Obtain the Overdispersion Parameters from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated overdispersion parameters of a
\code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
overdispersions(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained from the 
  function \code{BUSseq_MCMC}.}
}
\value{
\item{phi.est}{The estimated overdispersion parameters, a B by G matrix,
where phi_{bg} is the overdispersion parameter of gene g in batch b.}
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
#"BUSseqfits_example" is a sample BUSseqfits object
class(BUSseqfits_example)
phi.est <- overdispersions(BUSseqfits_example)
}