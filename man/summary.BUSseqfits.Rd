\name{summary.BUSseqfits}
\alias{summary.BUSseqfits}
\title{
Summarize a \code{BUSseqfits} Object
}
\description{
Call the function to summarize a \code{BUSseqfits} object, 
the output of \code{BUSseq_MCMC}.
}
\usage{
\method{summary}{BUSseqfits}(object, ...)
}
\arguments{
  \item{object}{
A \code{BUSseqfits} object obtained from the function 
\code{BUSseq_MCMC}.
}
\item{...}{not used.}
}
\value{
Summarize the results from the output of \code{BUSseq_MCMC}.
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
summary(BUSseqfits_example)
}
