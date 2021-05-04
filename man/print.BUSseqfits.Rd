\name{print.BUSseqfits}
\alias{print.BUSseqfits}
\title{
Print the Output by \code{BUSseq_MCMC}
}
\description{
The function prints a \code{BUSseqfits} Object.
}
\usage{
\method{print}{BUSseqfits}(x, ...)
}
\arguments{
  \item{x}{
A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.
}
  \item{\dots}{not used.}
}
\value{
Print the results from the output by the function \code{BUSseq_MCMC}.
}
\references{
Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effect Correction, Manuscript.
}
\author{
Fangda Song
}
\examples{
#"BUSseqfits_example" is a sample BUSseqfits object
class(BUSseqfits_example)
print(BUSseqfits_example)
}
