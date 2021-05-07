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
A \code{BUSseqfits} object obtained from the function 
\code{BUSseq_MCMC}.
}
  \item{\dots}{not used.}
}
\value{
Print the results from the output by the function \code{BUSseq_MCMC}.
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
print(BUSseqfits_example)
}
