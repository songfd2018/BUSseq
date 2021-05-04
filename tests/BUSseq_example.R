#######################################
# Apply BUSseq to the simulation data #
#######################################
library(BUSseq)
ObservedCounts <- BUSseqfits_example$CountData_raw
BUSseqfits_res <- BUSseq_MCMC(ObservedCounts, n.celltypes = 4, seed = 2353, n.iterations = 1000)
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