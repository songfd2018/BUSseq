\name{dropout_coefficient_values}
\alias{dropout_coefficient_values}
\title{
Obtain the Coefficients of the Logistic Regression for the Dropout
Events from the Output of the \code{BUSseq_MCMC} Function
}
\description{
The function gives the intercept and odds ratio of the logistic regression 
for dropout events in the output \code{SingleCellExperiment} object.
}
\usage{
dropout_coefficient_values(sce_BUSseqfit)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object by 
  the function \code{BUSseq_MCMC}.}
}
\value{
\item{gamma.est}{The estimated intercept and log ratio of 
the logistic regression for dropout events, a 2-diminsional vector.}
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
gamma.est <- dropout_coefficient_values(BUSseqfits_example)
}