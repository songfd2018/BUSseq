\name{summary.CountData}
\alias{summary.CountData}
\title{
Summarize a \code{CountData} Object
}
\description{
Call the function to summarize a \code{CountData} object, the output of \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
}
\usage{
\method{summary}{CountData}(object, ...)
}
\arguments{
  \item{object}{
A \code{CountData} object obtained from the function \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
}
\item{...}{not used.}
}
\value{
Summarize the results from the output of \code{raw_read_counts}, \code{imputed_read_counts} or \code{corrected_read_counts}.
}
\references{
Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effect Correction, Manuscript.
}
\author{
Fangda Song
}
\examples{
#"BUSseqfits_example" is a sample BUSseqfits object
Example_CountData_raw <- raw_read_counts(BUSseqfits_example)
summary(Example_CountData_raw)
}
