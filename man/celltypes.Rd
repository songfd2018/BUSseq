\name{celltypes}
\alias{celltypes}
\title{
Obtain the Cell-type Indicators from a \code{BUSseqfits} Object
}
\description{
The function gives the cell-type indicators of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
celltypes(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object output by the function \code{BUSseq_MCMC}.}
}
\value{
\item{w.est}{The estimated cell-type indicators, an R list with length equal to the batch number. The b-th element is the estimated cell-type indicator vector for cells in batch b.}
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
celltypes_est <- celltypes(BUSseqfits_example)
}
