\name{corrected_read_counts}
\alias{corrected_read_counts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Obtain the Corrected Read Count Matrices from 
  a \code{BUSseqfits} Object
}
\description{
  The function generates a version of count data, for which the batch
  effects are removed and the biological variabilities are retained. 
  We develop a quantile match approach based on the idea of inverse
  transform sampling. The users can perform downstream analysis on the
  corrected read count matrix, such as clustering, differentially
  expressed gene identification and so on, as if all the data were
  measured in a single batch.
}
\usage{
  corrected_read_counts(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object output by 
  the function \code{BUSseq_MCMC}.}
}
\value{
  \item{CountData_corrected}{The corrected read counts, 
  a \code{CountData} with length equal to the batch number. 
  Each element is a read count matrix for a specific batch, 
  where each row corresponds to a gene and each column 
  represents a cell.}
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
  Example_CountData_corrected <- corrected_read_counts(BUSseqfits_example)
  class(Example_CountData_corrected)
}