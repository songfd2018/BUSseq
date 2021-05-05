\name{baseline_expression_values}
\alias{baseline_expression_values}
\title{
Obtain the Log-Scale Baseline Expression levels from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated log-scale baseline expression levels from a
\code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
baseline_expression_values(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{A BUSseqfits object obtained from the function 
  \code{BUSseq_MCMC}.}
}
\value{
\item{alpha.est}{The estimated log-scale baseline expression levels, a 
G-dimensional vector whose g-th element is the estimated log-scale mean gene 
expression level of gene g in the first cell type.}
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
alpha.est <- baseline_expression_values(BUSseqfits_example)
}
