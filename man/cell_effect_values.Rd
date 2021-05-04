\name{cell_effect_values}
\alias{cell_effect_values}
\title{
Obtain the cell-specific size effects from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated cell-specific size effects of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
cell_effect_values(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained by the function \code{BUSseq_MCMC}.}
}
\value{
\item{delta.est}{The estimated cell-specific size effects, a list with length equal to the number of batch, whose b-th element is a vector of cell-specific effect in the b-th batch. Note that the first element in each vector is zero as the first cell in each batch is taken as the reference cell.}
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
delta.est <- cell_effect_values(BUSseqfits_example)
}
