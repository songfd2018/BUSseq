\name{print.CountData}
\alias{print.CountData}
\title{
  Print the Output by the Function \code{raw_read_counts}, 
  \code{imputed_read_counts} or \code{corrected_read_counts}
}
\description{
  The function prints a \code{CountData} Object.
}
\usage{
  \method{print}{CountData}(x, ...)
}
\arguments{
  \item{x}{
    A \code{CountData} object obtained from the function 
    \code{raw_read_counts}, \code{imputed_read_counts} or 
    \code{corrected_read_counts}.
  }
  \item{\dots}{not used.}
}
\value{
  Print the results from the output by the function 
  \code{raw_read_counts}, \code{imputed_read_counts} or 
  \code{corrected_read_counts}.
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
  Example_CountData_raw <- raw_read_counts(BUSseqfits_example)
  print(Example_CountData_raw)
}