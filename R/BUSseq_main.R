##############################
# The MCMC Sampling Function #
##############################
BUSseq_MCMC <- function(ObservedData, n.celltypes,
                        seed = round(runif(1,1,10000)), n.cores = 8, 
                        # set iterations
                        n.iterations=2000, n.burnin=floor(n.iterations / 2), 
                        n.unchanged=min(floor(n.iterations * 0.3), 500),
                        n.output = n.iterations / 10,  
                        working_dir = getwd(), 
                        # batch indicator of dropout events
                        Drop_ind = rep(TRUE, length(ObservedData)),
                        # false discovery rate
                        fdr = 0.05,
                        # hyperparameters
                        hyper_pi = 2, hyper_gamma0 = 3, 
                        hyper_gamma1 = c(0.001,0.01),
                        hyper_alpha = 5, tau1sq = 50,
                        hyper_p = c(1,3), hyper_tau0sq = c(2,0.01),
                        hyper_nu = 5, hyper_delta = 5,
                        hyper_phi = c(1,0.1)){
  
  Read <- NULL #input data Y
  nb <- NULL
  batch_name_vec <- NULL
  if(is(ObservedData, "list")){       #The input data format is list
     #Each row represents a gene and each column is a cell
     B <- length(ObservedData)# batch number   
     for(b in seq_len(B)){#input data
        if(is(ObservedData[[b]],"data.frame") | is(ObservedData[[b]],"matrix")){
           ObservedData[[b]] <- as.matrix(ObservedData[[b]])
        }else{
           stop("Each element of ObservedData must be a",
                " \"matrix\" or \"data.frame\" object!\n")
        }
        Read <- cbind(Read,ObservedData[[b]])
        nb <- c(nb,ncol(ObservedData[[b]]))
        batch_name_vec <- c(batch_name_vec, paste("Batch",b))
     }
     
     G <- nrow(Read)
     N <- sum(nb)
     
     res <- SingleCellExperiment(assays = list(counts = Read),
            colData = DataFrame(Batch_ind = factor(rep(batch_name_vec,nb))))
     ###################################
     ###Test the consistency of genes###
  }else if(is(ObservedData, "SingleCellExperiment")){
     
     # check whether batch indicators exist in the colData 
     # of SingleCellExperiment
     if(sum(names(colData(ObservedData)) == "Batch_ind") == 0){
        message("Please set the batch indicators of each cell by",
                " \"colData(ObservedData) = DataFrame(Batch_ind = ...)\"")
        stop("The batch indicators does not exist!\n")
     }
     
     Read <- assay(ObservedData, "counts")
     if(is.null(Read)){
        stop("Fail to load read count data from the input!\n")
     }
     
     if(is(Read,"data.frame") | is(Read,"matrix")){
        Read <- as.matrix(Read)
     }else{
        stop("The \"counts\" assay must be a",
             " \"matrix\" or \"data.frame\" object!\n")
     }
     
     N <- ncol(Read)
     G <- nrow(Read)
     
     batch_ind <- unlist(colData(ObservedData))
     
     # search the batch indicators to get batch names
     B <- length(unique(batch_ind))# batch number 
     num_change <- 1
     batch_name_vec <- batch_ind[1]
     nb <- rep(NA, B)
     
     end_last_batch <- 0
     for(i in seq_len(N-1)){
        if(batch_ind[i] != batch_ind[i+1]){
           num_change <- num_change + 1
           if(num_change > B){
              stop("Cells from the same batch should be arranged together.\n")
           }
           batch_name_vec <- c(batch_name_vec, batch_ind[i+1])
           nb[num_change - 1] <- i - end_last_batch
           end_last_batch <- i
        }
     }
     nb[B] <- N - end_last_batch
     
     res <- ObservedData
     
  }else{
     stop("ObservedData must be a \"list\" or \"SingleCellExperiment\" 
          object!\n")
  }
  
  if(B < 2){
    stop("The batch number must be greater than one.\n")
  }
  
  K <- n.celltypes

  if(sum(K > nb) > 0){
    stop("The sample size in any batch must be greater",
                " than the assumed cell type number.\n")
  }
  
  ####Record the posterior sampling on the hard disk
  if(!dir.exists(working_dir)){
    dir.create(working_dir)
  }
  sampling_dir <- paste0(working_dir,"/MCMC_sampling_K",K)
  if(dir.exists(sampling_dir)){
    unlink(sampling_dir,recursive = TRUE)
  }
  dir.create(sampling_dir, showWarnings=FALSE)
  
  #################
  # MCMC sampling #
  #################
  t.start <- Sys.time()
  message("   conducting the posterior sampling...\n")
  
  # prepare the input to C++ program
  dim <- c(N, G, B, K, nb, Drop_ind)
  iter_infor <- c(n.iterations, n.output, n.unchanged, n.burnin)
  hyper <- c(hyper_pi, hyper_gamma0, hyper_gamma1,
             hyper_alpha, tau1sq,
             hyper_p, hyper_tau0sq,
             hyper_nu, hyper_delta,
             hyper_phi)
  
  mcmc_sample<-.C("BUSseq_MCMC", 
                  # count data
                  y_obs = as.integer(t(Read)), 
                  # dimension information
                  dim = as.integer(dim),
                  # seed and num of cores
                  seed = as.integer(seed),
                  n.cores = as.integer(n.cores),
                  # iteration setting
                  iter_infor = as.integer(iter_infor),
                  # output directory
                  dir_output = as.character(sampling_dir),
                  # hyperparameter
                  hyper = as.double(hyper),
                  # x_imputed
                  x_imputed = as.integer(rep(0,N * G))
  )
  
  x_imputed <- matrix(mcmc_sample$x_imputed, G, N, byrow = TRUE)
  
  t.end <- Sys.time()
  message("   The MCMC sampling takes: ", 
                 round(difftime(t.end, t.start,units="mins"), 3), 
                 " mins", "\n")
  
  #######################
  # Posterior inference #
  #######################
  message("   conducting the posterior inference...\n")
  t.start<-Sys.time()
  post_inference <- .C("BUSseq_inference",
                       # count data
                       y_obs = as.integer(t(Read)), 
                       # dimension information
                       dim = as.integer(dim),
                       # num of cores
                       n.cores = as.integer(n.cores),
                       # iteration setting
                       iter_infor = as.integer(iter_infor),
                       # output directory
                       dir_output = as.character(sampling_dir),
                       # false discovery rate
                       fdr = as.double(fdr),
                       # posterior mean, mode or standard deviation
                       alpha_est = as.double(rep(0,G)), 
                       alpha_sd = as.double(rep(0,G)),
                       beta_est = as.double(matrix(0,K,G)), 
                       beta_sd = as.double(matrix(0,K,G)),
                       nu_est = as.double(matrix(0,B,G)), 
                       nu_sd = as.double(matrix(0,B,G)),
                       delta_est = as.double(rep(0,N)), 
                       delta_sd = as.double(rep(0,N)),
                       gamma_est = as.double(matrix(0, 2, B)),
                       gamma_sd = as.double(matrix(0,2,B)),
                       phi_est = as.double(matrix(0,B,G)), 
                       phi_sd = as.double(matrix(0,B,G)),
                       pi_est = as.double(matrix(0, K, B)), 
                       pi_sd = as.double(matrix(0, K, B)),
                       tau0_est = as.double(0), tau0_sd = as.double(0),
                       p_est = as.double(0), p_sd = as.double(0),
                       w_est = as.integer(rep(0, N)), 
                       PPI_est = as.double(matrix(0, K, G)),
                       D_est = as.integer(rep(0,G)), BIC = as.double(0))
  
  alpha.est <- post_inference$alpha_est
  alpha.sd <- post_inference$alpha_sd
  beta.est <- matrix(post_inference$beta_est, G, K, byrow = TRUE)
  beta.sd <- matrix(post_inference$beta_sd, G, K, byrow = TRUE)
  nu.est <- matrix(post_inference$nu_est, G, B, byrow = TRUE)
  nu.sd <- matrix(post_inference$nu_sd, G, B, byrow = TRUE)
  delta.est <- post_inference$delta_est
  delta.sd <- post_inference$delta_sd
  gamma.est <- matrix(post_inference$gamma_est, B, 2, byrow = TRUE)
  gamma.sd <- matrix(post_inference$gamma_sd, B, 2, byrow = TRUE)
  phi.est <- matrix(post_inference$phi_est, G, B, byrow = TRUE)
  phi.sd <- matrix(post_inference$phi_sd, G, B, byrow = TRUE)
  pi.est <- matrix(post_inference$pi_est, B, K, byrow = TRUE)
  pi.sd <- matrix(post_inference$pi_sd, B, K, byrow = TRUE)
  tau0.est <- post_inference$tau0_est
  tau0.sd <- post_inference$tau0_sd
  p.est <- post_inference$p_est
  p.sd <- post_inference$p_sd
  w.est <- post_inference$w_est + 1
  PPI.est <- matrix(post_inference$PPI_est, G, K, byrow = TRUE)
  D.est <- post_inference$D_est
  BIC <- post_inference$BIC
  
  t.end<-Sys.time()
  message("   Posterior inference takes: ", 
                 round(difftime(t.end, t.start,units="mins"), 3),
                 " mins", "\n")
  
  # Construct the output
  assay(res, withDimnames = FALSE, "imputed_data") <- x_imputed
  
  # add estimated cell type labels to column-level metadata
  int_colData(res)$BUSseq <- DataFrame(CellLabels = w.est,
                                       RelativeSizeFactors = exp(delta.est))
  
  # add intrinsic gene indicators to row-level metadata
  intri_ind <- rep("No", G)
  intri_ind[which(D.est == 1)] <- "Yes"
  int_elementMetadata(res)$BUSseq <- DataFrame(IntrinsicGene = intri_ind)
  
  int_metadata(res)$BUSseq <- list(
     #dimensions
     n.cell=N, n.gene=G, n.batch=B, n.perbatch=nb, n.celltype=K,
     n.iter=n.iterations, n.burnin = n.burnin, seed=seed,
     #posterior mean or mode of parameters
     gamma.est=gamma.est, alpha.est=alpha.est, 
     beta.est=beta.est, nu.est=nu.est,
     delta.est=delta.est, phi.est=phi.est, pi.est=pi.est,
     w.est=w.est, p.est=p.est, tau0.est=tau0.est,
     PPI.est=PPI.est, D.est = D.est,
     #posterior sd of pararmeters
     gamma.sd = gamma.sd, alpha.sd=alpha.sd, 
     beta.sd=beta.sd, nu.sd=nu.sd,
     delta.sd=delta.sd, phi.sd=phi.sd, pi.sd=pi.sd,
     p.sd=p.sd, tau0.sd=tau0.sd,
     BIC = BIC)
  
  return(res)
  
}

##################################
# Useful Outputs from BUSseqfits #
##################################
#obtain the cell type indicators for samples
celltypes <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         
         .B <- int_metadata(sce_BUSseqfit)$BUSseq$n.batch
         .nb <-int_metadata(sce_BUSseqfit)$BUSseq$n.perbatch
         .w <- int_metadata(sce_BUSseqfit)$BUSseq$w.est
         
         index_cell <- 0
         for(b in seq_len(.B)){
            
            message("Batch ", b, " cells' cell type indicators: ",
                    .w[index_cell + 1],",",.w[index_cell + 2],",",
                    .w[index_cell + 3], "... ...\n")
            
            index_cell <- index_cell + .nb[b]
         }
         message("The output format is an N-dimensional verctor.\n")
         return(.w)
         
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the dropout intercept and odds ratio
dropout_coefficient_values <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .gamma<-int_metadata(sce_BUSseqfit)$BUSseq$gamma.est
         message("The output format is a matrix.\n")
         message("Each row represents a batch, the first column corresponds",
                 " to intercept and the second column is the odd ratio.\n")
         return(.gamma)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the log-scale baseline expression values
baseline_expression_values <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .alpha <- int_metadata(sce_BUSseqfit)$BUSseq$alpha.est
         message("The output format is a vector.\n")
         return(.alpha)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}
   
#obtain the cell type effects
celltype_effects <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .beta <- int_metadata(sce_BUSseqfit)$BUSseq$beta.est
         message("The output format is a matrix.\n")
         message("Each row represents a gene, and each column corresponds",
                 " to a cell type.\n")
         return(.beta)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the cell-tpye-specific mean expression levels
celltype_mean_expression <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .alpha <- int_metadata(sce_BUSseqfit)$BUSseq$alpha.est
         .beta <- int_metadata(sce_BUSseqfit)$BUSseq$beta.est
         mu <- exp(.alpha+.beta)
         message("The output format is a matrix.\n")
         message("Each row represents a gene, and each column corresponds",
                 " to a cell type.\n")
         return(mu)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the location batch effects
location_batch_effects <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .nu <- int_metadata(sce_BUSseqfit)$BUSseq$nu.est
         message("The output format is a matrix.\n")
         message("Each row represents a gene, and each column",
                 " corresponds to a batch.\n")
         return(.nu)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the scale batch effects
overdispersions <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .phi <- int_metadata(sce_BUSseqfit)$BUSseq$phi.est
         message("The output format is a matrix.\n")
         message("Each row represents a gene, and each column",
                 " corresponds to a batch.\n")
         return(.phi)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the cell-specific size effects
cell_effect_values <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .delta <- int_metadata(sce_BUSseqfit)$BUSseq$delta.est
         message("The output format is an N-dimensional vector.\n")
         return(.delta)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the intrinsic genes
intrinsic_genes_BUSseq <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         intrinsic_genes <- int_elementMetadata(sce_BUSseqfit)$BUSseq
         return(intrinsic_genes)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the BIC score
BIC_BUSseq <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         .BIC <- int_metadata(sce_BUSseqfit)$BUSseq$BIC
         message("BIC is ", .BIC, "\n")
         message("The output is a scalar.\n")
         return(.BIC)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the raw count data
raw_read_counts <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
         Read_raw <- assay(sce_BUSseqfit, "counts")
         message("The output format is a matrix, in which each row",
                 " represents a gene and each column does a cell.\n")
         return(Read_raw)
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the underlying true count data
imputed_read_counts <- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         Read_imputed <- assay(sce_BUSseqfit, "imputed_data")
         message("The output format is a matrix, in which each row",
                 " represents a gene and each column does a cell.\n")
         return(Read_imputed)
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

#obtain the corrected count data
corrected_read_counts<- function(sce_BUSseqfit){
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         Read_imputed <- assay(sce_BUSseqfit, "imputed_data")
         .B <- int_metadata(sce_BUSseqfit)$BUSseq$n.batch
         .nb <- int_metadata(sce_BUSseqfit)$BUSseq$n.perbatch
         .N <- int_metadata(sce_BUSseqfit)$BUSseq$n.cell
         .G <- int_metadata(sce_BUSseqfit)$BUSseq$n.gene
         .K <- int_metadata(sce_BUSseqfit)$BUSseq$n.celltype
         .gamma <- int_metadata(sce_BUSseqfit)$BUSseq$gamma.est
         .logmu <- int_metadata(sce_BUSseqfit)$BUSseq$alpha.est + 
            int_metadata(sce_BUSseqfit)$BUSseq$beta.est
         .nu <- int_metadata(sce_BUSseqfit)$BUSseq$nu.est
         .delta <- int_metadata(sce_BUSseqfit)$BUSseq$delta.est
         .phi <- int_metadata(sce_BUSseqfit)$BUSseq$phi.est
         .w <- int_metadata(sce_BUSseqfit)$BUSseq$w.est
         
         message("   correcting read counts...\n")
         Read_corrected <- Read_imputed
         cell_index <- 1
         for(b in seq_len(.B)){
            for(i in seq_len(.nb[b])){
               unc_mu <- exp(.logmu[,.w[cell_index]] + 
                                .nu[,b] + .delta[cell_index])
               p_x<-pnbinom(Read_imputed[,cell_index],size=.phi[,b], mu=unc_mu)
               p_xminus1<-pnbinom(Read_imputed[,cell_index]-1,size=.phi[,b], 
                                  mu=unc_mu)
               u <- runif(.G,min=p_xminus1,max=p_x)
               u <- apply(cbind(u,0.9999),1,min)
               cor_mu <- exp(.logmu[,.w[cell_index]])
               Read_corrected[,cell_index] <- 
                  qnbinom(u,size =.phi[,1], mu=cor_mu)
               cell_index <- cell_index + 1
            }
         }
         
         assay(sce_BUSseqfit, withDimnames = FALSE, 
               "corrected_data") <- Read_corrected
         message("The corrected read count matrix is added into the output",
                 " \"SingleCellExperiment\" object.\n")
         return(sce_BUSseqfit)
         
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}

########################################################################
# Visualization
########################################################################
#visualize the read counts data by stacking all gene expression matrices 
heatmap_data_BUSseq <- function(sce_BUSseqfit, 
                        data_type = c("Raw","Imputed","Corrected"), 
                        gene_set=NULL, 
                        project_name= paste0("BUSseq_heatmap_",data_type), 
                        image_dir=NULL, color_key_seq=NULL, 
                        image_width=1440, image_height=1080){
   
   if(is(sce_BUSseqfit,"SingleCellExperiment")){
      if(!is.null(int_metadata(sce_BUSseqfit)$BUSseq)){
         if(length(data_type) > 1){
            data_type <- data_type[1]
         }
         if(data_type == "Raw"){
            count_data <- assay(sce_BUSseqfit, "counts")
         }else if(data_type == "Imputed"){
            count_data <- assay(sce_BUSseqfit, "imputed_data")
         }else if(data_type == "Corrected"){
            if(sum(assayNames(sce_BUSseqfit) == "corrected_data") == 0){
               message("The corrected read count matrix",
                       " has not been generated.\n")
               stop("Please run the \"corrected_read_counts\" function!\n")
            }else{
               count_data <- assay(sce_BUSseqfit, "corrected_data")
            }

         }else{
            stop("Please select the data type from ",
                 "\"Raw\", \"Imputed\" or \"Corrected\" to draw the heatmap!")
         }
        
         
         
         .B <- int_metadata(sce_BUSseqfit)$BUSseq$n.batch
         .nb <- int_metadata(sce_BUSseqfit)$BUSseq$n.perbatch
         .G <- int_metadata(sce_BUSseqfit)$BUSseq$n.gene
         
         if(is.null(gene_set)) gene_set <- seq_len(.G)
         #heatmap cell colors
         colfunc <- colorRampPalette(c("grey", "black"))
         #batch colors
         color_batch_func <- colorRampPalette(
            c("#EB4334","#FBBD06","#35AA53","#4586F3"))
         color_batch <- color_batch_func(.B)
         color_batch2 <- rep(color_batch, .nb)
         
         log1p_mat <- log1p(count_data)
         log1p_mat_interest <- log1p_mat[gene_set, ]
         if(is.null(color_key_seq)){
            range_data <- range(log1p_mat_interest)
            color_key_seq <- seq(from=floor(range_data[1]) - 0.5, 
                                 to=ceiling(range_data[2]) + 0.5, length.out=11)
         }
         if(is.null(image_dir)) image_dir <- "./image"
         #create the folder
         dir.create(image_dir,showWarnings=FALSE)
         png(paste(image_dir,"/",project_name,"_log1p_data.png",sep=""),
             width=image_width, height=image_height)
         heatmap.2(log1p_mat_interest,
                   dendrogram="none",#with cluster tree
                   Rowv=FALSE, Colv=FALSE,
                   labRow=FALSE, labCol=FALSE,
                   ColSideColors=color_batch2,
                   col=colfunc(length(color_key_seq)-1),breaks=color_key_seq,
                   density.info="histogram",
                   hclustfun=function(c)hclust(c,method="average"),
                   keysize=0.8, cexRow=0.5,trace="none")#font size
         dev.off()
         
      }else{
         stop("The output of the \"BUSseq_MCMC\" function does not exist!\n")
      }
   }else{
      stop("The input must be a \"SingleCellExperiment\" object!\n")
   }
}
