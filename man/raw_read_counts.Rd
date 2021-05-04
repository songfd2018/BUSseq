\name{raw_read_counts}
\alias{raw_read_counts}
\title{
  Obtain the Raw Read Count Matrices from a \code{BUSseqfits} Object
}
\description{
  The function gives the raw read counts of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
  raw_read_counts(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
\value{
  \item{CountData_raw}{The raw read counts, a \code{CountData} object with length equal to the batch number. Each element is a read count matrix for a specific batch, where each row corresponds to a gene and each column represents a cell.}
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
  Example_CountData_raw <- raw_read_counts(BUSseqfits_example)
  class(Example_CountData_raw)
}