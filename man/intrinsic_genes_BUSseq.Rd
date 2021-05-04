\name{intrinsic_genes_BUSseq}
\alias{intrinsic_genes_BUSseq}
\title{
Obtain the Intrinsic Gene Indices from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated intrinsic gene indices of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
intrinsic_genes_BUSseq(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{ A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
\value{
\item{intrinsic_genes}{The intrinsic gene indices. The return is a vector.}
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
D.est <- intrinsic_genes_BUSseq(BUSseqfits_example)
}