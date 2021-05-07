\name{BIC_BUSseq}
\alias{BIC_BUSseq}
\title{
Obtain BIC from a \code{BUSseqfits} Object
}
\description{
The function gives the Bayesian Informtion Criterion (BIC) value for a 
\code{BUSseqfits} object.
}
\usage{
BIC_BUSseq(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object obtained from the function 
  \code{BUSseq_MCMC}.}
}
\value{
\item{BIC_val}{The BIC value for a BUSseqfits object.}
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
Example_BIC <- BIC_BUSseq(BUSseqfits_example)
}