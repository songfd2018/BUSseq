\name{BIC_BUSseq}
\alias{BIC_BUSseq}
\title{
Obtain BIC from a \code{BUSseqfits} Object
}
\description{
The function gives the Bayesian Informtion Criterion (BIC) value for a \code{BUSseqfits} object.
}
\usage{
BIC_BUSseq(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
\value{
\item{BIC_val}{The BIC value for a BUSseqfits object.}
}
\references{
Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effects Correction, Manuscript.
}
\author{
Fangda Song
}
\examples{
#"BUSseqfits_example" is a sample BUSseqfits object
class(BUSseqfits_example)
Example_BIC <- BIC_BUSseq(BUSseqfits_example)
}