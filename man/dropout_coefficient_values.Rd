\name{dropout_coefficient_values}
\alias{dropout_coefficient_values}
\title{
Obtain the Coefficients of the Logistic Regression for the Dropout Events from a \code{BUSseqfits} Object
}
\description{
The function gives the intercept and odds ratio of the logistic regression for dropout events in a \code{BUSseqfits} object.
}
\usage{
dropout_coefficient_values(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
\value{
\item{gamma.est}{The estimated intercept and log ratio of the logistic regression for dropout events, a 2-diminsional vector.}
}
\references{
Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effect Correction, Manuscript.
}
\author{
Fangda Song
}
\examples{
#"BUSseqfits_example" is a sample of BUSseqfits object
class(BUSseqfits_example)
gamma.est <- dropout_coefficient_values(BUSseqfits_example)
}