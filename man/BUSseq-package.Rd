\name{BUSseq-package}
\alias{BUSseq-package}
\alias{BUSseq}
\docType{package}
\title{
\packageTitle{BUSseq}
}
\description{
\packageDescription{BUSseq}
}
%\details{
%The DESCRIPTION file:
%\packageDESCRIPTION{BUSseq}
%}
%\packageIndices{BUSseq}{
%~ An overview of how to use the package, 
% including the most important functions ~~
%}
\author{
\packageAuthor{BUSseq}

Maintainer: \packageMaintainer{BUSseq}
}
\references{
Song, Fangda, Ga Ming Angus Chan, and Yingying Wei. 
Flexible experimental designs for valid single-cell RNA-sequencing 
experiments allowing batch effects correction. 
Nature communications 11, no. 1 (2020): 1-15.
}
%\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
\donttest{
#######################################
# Apply BUSseq to the Simulation Data #
#######################################
library(SingleCellExperiment)
RawCountData <- assay(BUSseqfits_example, "counts")
batch_ind <- colData(BUSseqfits_example)
sce <- SingleCellExperiment(assays = list(counts = RawCountData),
                            colData = DataFrame(Batch_ind = batch_ind))
BUSseqfits_res <- BUSseq_MCMC(ObservedData = sce, 
                              seed = 1234, n.cores = 2,
                              n.celltypes = 4, n.iterations = 500)

################################################
# Extract Estimates from the BUSseqfits Object #
################################################

#return cell type indicators
w.est <- celltypes(BUSseqfits_res)

#return the intercept and odds ratio of the logistic regression
#for dropout events
gamma.est <- dropout_coefficient_values(BUSseqfits_res)

#return the log-scale baseline expression values
alpha.est <-  baseline_expression_values(BUSseqfits_res)

#return the cell-type effects
beta.est <- celltype_effects(BUSseqfits_res)

#return the mean expression levels
mu.est <- celltype_mean_expression(BUSseqfits_res)

#return the cell-specific global effects
delta.est <- cell_effect_values(BUSseqfits_res)

#return the location batch effects
nu.est <- location_batch_effects(BUSseqfits_res)

#return the overdispersion parameters
phi.est <- overdispersions(BUSseqfits_res)

#return the intrinsic gene indices
D.est <- intrinsic_genes_BUSseq(BUSseqfits_res)

#return the BIC value
BIC <- BIC_BUSseq(BUSseqfits_res)

#return the raw read count matrix
CountData_raw <- raw_read_counts(BUSseqfits_res)

#return the imputed read count matrix
CountData_imputed <- imputed_read_counts(BUSseqfits_res)

#return the corrected read count matrix
CountData_corrected <- corrected_read_counts(BUSseqfits_res)

#################
# Visualization #
#################
#generate the heatmap of raw read count data
heatmap_data_BUSseq(BUSseqfits_res, project_name="Heatmap_raw")

#generate the heatmap of imputed read count data
heatmap_data_BUSseq(BUSseqfits_res, data_type = "Imputed",
                    project_name="Heatmap_imputed")

#generate the heatmap of corrected read count data
heatmap_data_BUSseq(BUSseqfits_res, data_type = "Corrected", 
                    project_name="Heatmap_corrected")
}
}
