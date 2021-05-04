\name{summary.BUSseqfits}
\alias{summary.BUSseqfits}
\title{
Summarize a \code{BUSseqfits} Object
}
\description{
Call the function to summarize a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
\method{summary}{BUSseqfits}(object, ...)
}
\arguments{
  \item{object}{
A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.
}
\item{...}{not used.}
}
\value{
Summarize the results from the output of \code{BUSseq_MCMC}.
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
summary(BUSseqfits_example)
}
