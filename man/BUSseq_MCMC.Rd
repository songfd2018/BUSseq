\name{BUSseq_MCMC}
\alias{BUSseq_MCMC}
\title{
Condcut MCMC sampling and posterior inference for the BUSseq Model
}
\description{
The function \code{BUSseq_MCMC} runs a Markov chain Monte Carlo (MCMC) algorithm to fit the Batch effects correction with Unknown Subtypes for scRNA-seq data (BUSseq) model. BUSseq is an interpretable Bayesian hierarchical model that closely follows the data-generating mechanism of scRNA-seq experiments. BUSseq can simultaneously correct batch effects, cluster cell types, impute missing data caused by dropout events and detect differentially expressed genes without requiring a preliminary normalization step. We adopt MCMC algorithm to conduct posterior inference for the BUSseq model. Here, we denote the batch number as B, the gene number as G and the user-specific cell-type number as K.
}
\usage{
BUSseq_MCMC(ObservedData, n.celltypes, seed = round(runif(1,1,10000)), n.cores = 8,
            n.iterations = 2000, n.burnin = floor(n.iterations/2),
            n.unchanged = min(floor(n.iterations * 0.3), 500),
            n.output = n.iterations/10, working_dir = getwd(),
            Drop_ind = rep(TRUE, length(ObservedData)), fdr = 0.05, 
            hyper_pi = 2, hyper_gamma0 = 3, hyper_gamma1 = c(0.001, 0.01), 
            hyper_alpha = 5, tau1sq = 50, hyper_p = c(1, 3), 
            hyper_tau0sq = c(2, 0.01), hyper_nu = 5,
            hyper_delta = 5, hyper_phi = c(1, 0.1))
}
\arguments{
  \item{ObservedData}{
  \code{ObservedData} is an R list with the length equal to the batch number. The \code{b}-th element of \code{ObservedData} is the read count matrix of batch \code{b}, where each row corresponds to a gene, and each column corresponds to a cell.
}
  \item{n.celltypes}{
  \code{n.celltypes} is an integer and represents the number of cell types, which needs to be specified by the user.
}
  \item{seed}{
  \code{seed} is an integer for the Random Number Generator.
}
  \item{n.cores}{
  \code{n.cores} is an integer and denotes the number of cores used for the parallel MCMC algorithm.
}
  \item{n.iterations}{
  \code{n.iterations} is the total number  of iterations of the MCMC algorithm. The default is 2000.
}
  \item{n.burnin}{
  \code{n.burnin} is the number of burn-in iterations of the MCMC sampling. The default is a half of \code{n.iterations} 
}
  \item{n.unchanged}{
  \code{n.unchanged} is the number of iterations where the MCMC does not update the hyperparameter \code{p} and \code{tau0} of the slab and spike prior of cell-type effects. The default is the minimum of one third of \code{n.iterations} and 500.
}
  \item{n.output}{
  \code{n.output} is the number of iterations per hard-disk writing of the posterior sampling. The default is one tenth of \code{n.iterations}.
}
  \item{working_dir}{
  \code{working_dir} is the directory to store the posterior samples. The default is the current directory.
}
  \item{Drop_ind}{
  \code{Drop_ind} is a Boolean vector with length equal to the number of batches and indicates which batch suffers from dropout events. 
}
  \item{fdr}{The false discovery rate level we want to control in order to identify intrinsic genes. The default is 0.05.}
  \item{hyper_pi}{
  \code{hyper_pi} is a scalar representing the hyperparameter of the Dirichlet prior for cell type proportions. The default is \code{2}.
}
  \item{hyper_gamma0}{
  \code{hyper_gamma0} is a scalar representing the variance of the normal prior for the intercept of the logistic regression for dropout events. The default is \code{3}
}
  \item{hyper_gamma1}{
  \code{hyper_gamma1} is a vector representing the hyperparameter of the Gamma prior for the log-odds ratio of dropout events. The default is \code{(0.001, 0.01)} such that the prior mean of \code{gamma1} is \code{0.1}.
}
  \item{hyper_alpha}{
  \code{hyper_alpha} is a scalar representing the variance of the normal prior for the log-scale baseline expression levels. The default is \code{5}.
}
  \item{tau1sq}{
  \code{tau1sq} is the slab variance of the spike-and-slab prior for the cell type effects. The default is \code{50}.
}
  \item{hyper_p}{
  \code{hyper_p} is a two-dimensional vector representing the two shape parameters of the beta prior for the proportion of intrinsic genes. The default is \code{c(1,3)}.
}
  \item{hyper_tau0sq}{
  \code{hyper_tau0sq} is a two-dimensional vector representing the shape and scale of the inverse gamma prior for the variance of the spike normal prior. The default is \code{c(2,0.01)}.
}
  \item{hyper_nu}{
  \code{hyper_nu} is a scalar representing the variance of the normal prior for the batch effects. The default is \code{5}.
}
  \item{hyper_delta}{
  \code{hyper_delta} is a scalar representing the variance of the normal prior for the cell-specific size effects. The default is \code{sqrt(5)}.
}
  \item{hyper_phi}{
  \code{hyper_phi} is a two-dimensional vector representing the shape and rate of the gamma prior for the overdispersion parameters. The default is \code{c(1,0.1)}.
}
}
\value{

List with the following elements:
\item{n.cell}{The total number of cells in all batches, a scalar.}
\item{n.gene}{The number of genes, a scalar.}
\item{n.batch}{The number of batches, a scalar.}
\item{n.perbatch}{The number of cells in each batch, a B-dimensional vector.}
\item{n.celltype}{The number of cell types specified by user, a scalar.}
\item{n.iter}{The total number of iterations applied in the MCMC algorithm, a scalar.}
\item{seed}{The seed for the MCMC algorithm.}
\item{n.burnin}{The number of iterations as burnin, a scalar.}
\item{CountData_raw}{The observed raw read count data. The return is a list with length equal to the batch number B. Each element is a read count matrix where each row represents a gene and each column represents a cell.}
\item{CountData_imputed}{The imputed true read count data. The return is a list with length equal to the batch number B. Each element is a read count matrix where each row represents a gene and each column represents a cell.}
\item{gamma.est}{The estimated intercept and odds ratio of the logistic regression for dropout events, a B-by-2 matrix.}
\item{alpha.est}{The estimated log-scale baseline expression levels, a G-dimensional vector whose g-th element is the estimated log-scale mean gene expression level of gene g in the first cell type.}
\item{beta.est}{The estimated cell-type effects, a G by K matrix, whose [g,k] element is the effects of cell type k on gene g compared with the first cell type. Note that the first column is zero as the first cell type is taken as the baseline cell type.}
\item{nu.est}{The estimated location batch effects, a B by G matrix, where [b,g] element is the location batch effect on gene g in the batch b compared with the first batch. Note that the first column is zero as the first batch is taken as the reference batch without batch effects.}
\item{delta.est}{The estimated cell-specific global effects, a list with length equal to the batch number B, whose b-th element is a vector of cell-specific effect in batch b. Note that the first element in each vector is zero as the first cell in each batch is taken as the reference cell.}
\item{phi.est}{The estimated overdispersion parameters, a B by G matrix, where phi_{bg} is the overdispersion parameter on gene g in batch b.}
\item{pi.est}{The estimated cell-type proportions across batches, a K by B matrix, whose [k,b] element is the estimated proportion of cell type k in batch b.}
\item{w.est}{The estimated cell-type indicators of each cell, a list with length equal to the batch number B whose b-th element is the estimated cell type in the b-th batch.}
\item{p.est}{The estimated proportion of differentially expressed genes compared with the first cell type, a scalar.}
\item{PPI.est}{The estimated posterior marginal probability of being differentially expressed genes compared with the first cell type. The return is G by K matrix. Noted that the first column consist of zeros as there is no differentially expressed gene compared with the cell type of itself.}
\item{D.est}{The intrinsic gene indicators. The return is an N-dimensional vector.}
\item{BIC}{The BIC value when K = \code{n.celltypes}, which is used to determine the number of cell types by varying the value of K.}
}
\references{
Fangda Song, Ga Ming Chan and Yingying Wei. Flexible Experimental Designs for Valid Single-cell RNA-sequencing Experiments Allowing Batch Effect Correction, Manuscript.
}
\author{
Fangda Song
}

\examples{

#######################################
# Apply BUSseq to the Simulation Data #
#######################################
library(BUSseq)
\dontrun{
ObservedCounts <- BUSseqfits_example$CountData_raw
BUSseqfits_res <- BUSseq_MCMC(ObservedCounts, n.celltypes = 4, seed = 123, n.iterations = 1000)
}
BUSseqfits_res <- BUSseqfits_example

################################################
# Extract Estimates from the BUSseqfits Object #
################################################

#return cell type indicators
w.est <- celltypes(BUSseqfits_res)
table(w.est)

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
heatmap_data_BUSseq(CountData_raw, project_name="Heatmap_raw")

#generate the heatmap of imputed read count data
heatmap_data_BUSseq(CountData_imputed, project_name="Heatmap_imputed")

#generate the heatmap of corrected read count data
heatmap_data_BUSseq(CountData_corrected, project_name="Heatmap_corrected")

}
