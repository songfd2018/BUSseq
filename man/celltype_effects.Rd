\name{celltype_effects}
\alias{celltype_effects}
\title{
Obtain the Cell-type Effects from a \code{BUSseqfits} Object
}
\description{
The function gives the estimated cell-type effects of a \code{BUSseqfits} object, the output of \code{BUSseq_MCMC}.
}
\usage{
celltype_effects(BUSseqfits_obj)
}
\arguments{
  \item{BUSseqfits_obj}{A \code{BUSseqfits} object obtained from the function \code{BUSseq_MCMC}.}
}
\value{
\item{beta.est}{The estimated cell-type effects, a G by K matrix, whose [g,k] element is the effects of cell type k on gene g compared with the first cell type. Note that the first column is zero as the first cell type is taken as the baseline cell type.}
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
beta.est <- celltype_effects(BUSseqfits_example)
}
