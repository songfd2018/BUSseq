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
                        hyper_pi = 2, hyper_gamma0 = 3, hyper_gamma1 = c(0.001,0.01),
                        hyper_alpha = 5, tau1sq = 50,
                        hyper_p = c(1,3), hyper_tau0sq = c(2,0.01),
                        hyper_nu = 5, hyper_delta = 5,
                        hyper_phi = c(1,0.1)){
  
  Read <- NULL #input data Y
  nb <- NULL
  
  if(is(ObservedData, "list")){       #The input data format is list
    #Each row represents a gene and each column is a cell
    B <- length(ObservedData)# batch number   
    for(b in seq_len(B)){#input data
      if(is(ObservedData[[b]],"data.frame") | is(ObservedData[[b]],"matrix")){
        ObservedData[[b]] <- as.matrix(ObservedData[[b]])
      }else{
        stop(paste0("Each element of ObservedData must be a \"matrix\" or \"data.frame\" object!\n"))
      }
      Read <- cbind(Read,ObservedData[[b]])
      nb <- c(nb,ncol(ObservedData[[b]]))
    }
    ###################################
    ###Test the consistency of genes###
  }else{
    stop(paste0("ObservedData must be a \"list\" object!\n"))
  }
  
  if(B < 2){
    stop("The batch number must be greater than one.\n")
  }
  
  
  K <- n.celltypes
  G <- nrow(Read)
  N <- sum(nb)
  
  if(sum(K > nb) > 0){
    stop(paste0("The sample size in any batch must be greater",
                " than the assumed cell type number.\n"))
  }
  
  ####Record the posterior samplingo on the hard disk
  if(!dir.exists(working_dir)){
    dir.create(working_dir)
  }
  #setwd(working_dir)
  sampling_dir <- paste0(working_dir,"/MCMC_sampling_K",K)
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
  message(paste0("   The MCMC sampling takes: ", 
                 round(difftime(t.end, t.start,units="mins"), 3), " mins", "\n"))
  
  #######################
  # Posterior inference #
  #######################
  message("   conducting the posterior inferences...\n")
  
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
                       alpha_est = as.double(rep(0,G)), alpha_sd = as.double(rep(0,G)),
                       beta_est = as.double(matrix(0,K,G)), beta_sd = as.double(matrix(0,K,G)),
                       nu_est = as.double(matrix(0,B,G)), nu_sd = as.double(matrix(0,B,G)),
                       delta_est = as.double(rep(0,N)), delta_sd = as.double(rep(0,N)),
                       gamma_est = as.double(matrix(0, 2, B)), gamma_sd = as.double(matrix(0,2,B)),
                       phi_est = as.double(matrix(0,B,G)), phi_sd = as.double(matrix(0,B,G)),
                       pi_est = as.double(matrix(0, K, B)), pi_sd = as.double(matrix(0, K, B)),
                       tau0_est = as.double(0), tau0_sd = as.double(0),
                       p_est = as.double(0), p_sd = as.double(0),
                       w_est = as.integer(rep(0, N)), PPI_est = as.double(matrix(0, K, G)),
                       D_est = as.integer(rep(0,G)), BIC = as.double(0))
  
  # alpha
  alpha.est <- post_inference$alpha_est
  alpha.sd <- post_inference$alpha_sd
  
  # beta
  beta.est <- matrix(post_inference$beta_est, G, K, byrow = TRUE)
  beta.sd <- matrix(post_inference$beta_sd, G, K, byrow = TRUE)
  
  # nu
  nu.est <- matrix(post_inference$nu_est, G, B, byrow = TRUE)
  nu.sd <- matrix(post_inference$nu_sd, G, B, byrow = TRUE)
  
  # delta
  delta.est <- post_inference$delta_est
  delta.sd <- post_inference$delta_sd
  
  # gamma
  gamma.est <- matrix(post_inference$gamma_est, B, 2, byrow = TRUE)
  gamma.sd <- matrix(post_inference$gamma_sd, B, 2, byrow = TRUE)
  
  # phi
  phi.est <- matrix(post_inference$phi_est, G, B, byrow = TRUE)
  phi.sd <- matrix(post_inference$phi_sd, G, B, byrow = TRUE)
  
  # pi
  pi.est <- matrix(post_inference$pi_est, B, K, byrow = TRUE)
  pi.sd <- matrix(post_inference$pi_sd, B, K, byrow = TRUE)
  
  # tau0
  tau0.est <- post_inference$tau0_est
  tau0.sd <- post_inference$tau0_sd
  
  # p
  p.est <- post_inference$p_est
  p.sd <- post_inference$p_sd
  
  # w
  w.est <- post_inference$w_est + 1
  
  # PPI
  PPI.est <- matrix(post_inference$PPI_est, G, K, byrow = TRUE)
  
  # D
  D.est <- post_inference$D_est
  
  # BIC
  BIC <- post_inference$BIC
  
  t.end<-Sys.time()
  message(paste0("   calculating posterior means and posterior takes: ", 
                 round(difftime(t.end, t.start,units="mins"), 3),
                 " mins", "\n"))
  ###Generate the output "BUSseqfits" object
  #transfer Read, x.post, delta.est, w.est, delta.sd as a list
  Read_list <- list()
  Read_sim_list <- list()
  delta.est_list <- list()
  delta.sd_list <- list()
  w_list <- list()
  cell_index <- 0
  for(b in seq_len(B)){
    Read_list[[b]] <- Read[,cell_index + seq_len(nb[b])]
    Read_sim_list[[b]] <- x_imputed[,cell_index + seq_len(nb[b])]
    delta.est_list[[b]] <- delta.est[cell_index + seq_len(nb[b])]
    delta.sd_list[[b]] <- delta.sd[cell_index + seq_len(nb[b])]
    w_list[[b]] <- w.est[cell_index + seq_len(nb[b])]
    cell_index <- cell_index + nb[b]
  }
  
  output <- list(CountData_raw=Read_list, CountData_imputed=Read_sim_list,
                 #dimensions
                 n.cell=N, n.gene=G, n.batch=B, 
                 n.perbatch=nb, n.celltype=K, 
                 n.iter=n.iterations, n.burnin = n.burnin,
                 seed=seed,
                 #posterior mean or mode of parameters
                 gamma.est=gamma.est, alpha.est=alpha.est, 
                 beta.est=beta.est, nu.est=nu.est,
                 delta.est=delta.est_list, phi.est=phi.est, pi.est=pi.est,
                 w.est=w_list, p.est=p.est, tau0.est=tau0.est,
                 PPI.est=PPI.est, D.est = D.est,
                 #posterior sd of pararmeters
                 gamma.sd = gamma.sd, alpha.sd=alpha.sd, 
                 beta.sd=beta.sd, nu.sd=nu.sd,
                 delta.sd=delta.sd_list, phi.sd=phi.sd, pi.sd=pi.sd,
                 p.sd=p.sd, tau0.sd=tau0.sd,
                 BIC = BIC)
  
  class(output) <- "BUSseqfits"
  
  return(output)
  
}

##################################
# Useful Outputs from BUSseqfits #
##################################
#obtain the cell type indicators for samples
celltypes <- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .B <- BUSseqfits_obj$n.batch
   .w <- BUSseqfits_obj$w.est
   
   for(b in seq_len(.B)){
      
      message(paste0("Batch ", b, " cells' cell type indicators: ",
            .w[[b]][1],",",.w[[b]][2],",",.w[[b]][3], 
            "... ...\n"))
   }
   message(paste0("The output format is a list with length",
         " equal to the batch number.\n"))
   message(paste0("Each element of the list is a cell type indicator",
                " vector in that batch.\n"))
   return(.w)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the dropout intercept and odds ratio
dropout_coefficient_values <- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .gamma<-BUSseqfits_obj$gamma.est
   message("The output format is a matrix.\n")
   message(paste0("Each row represents a batch, the first column corresponds",
                " to intercept and the second column is the odd ratio.\n"))
   return(.gamma)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the log-scale baseline expression values
baseline_expression_values <- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .alpha<-BUSseqfits_obj$alpha.est
   message("The output format is a vector.\n")
   return(.alpha)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}
   
#obtain the cell type effects
celltype_effects <- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .beta<-BUSseqfits_obj$beta.est
   message("The output format is a matrix.\n")
   message(paste0("Each row represents a gene, and each column corresponds",
                " to a cell type.\n"))
   return(.beta)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the cell-tpye-specific mean expression levels
celltype_mean_expression <- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .alpha<-BUSseqfits_obj$alpha.est
   .beta<-BUSseqfits_obj$beta.est
   
   mu <- exp(.alpha+.beta)
   
   message("The output format is a matrix.\n")
   message(paste0("Each row represents a gene, and each column corresponds",
                " to a cell type.\n"))
   return(mu)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the location batch effects
location_batch_effects <- function(BUSseqfits_obj){
   
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .nu<-BUSseqfits_obj$nu.est
   message("The output format is a matrix.\n")
   message(paste0("Each row represents a gene, and each column",
         " corresponds to a batch.\n"))
   return(.nu)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the scale batch effects
overdispersions <- function(BUSseqfits_obj){
   
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .phi<-BUSseqfits_obj$phi.est
   message("The output format is a matrix.\n")
   message(paste0("Each row represents a gene, and each column",
         " corresponds to a batch.\n"))
   return(.phi)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the cell-specific size effects
cell_effect_values <- function(BUSseqfits_obj){
   
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .delta<-BUSseqfits_obj$delta.est
   message(paste0("The output format is a list with length equal to",
         " the batch number.\n"))
   message(paste0("Each element of the list is a cell-specific",
                " size factor vector of that batch.\n"))
   return(.delta)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}


#obtain the instrinsic genes
intrinsic_genes_BUSseq <- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   
   intrinsic_genes <- which(BUSseqfits_obj$D.est==1)
   return(intrinsic_genes)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the BIC score
BIC_BUSseq <- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   .BIC<-BUSseqfits_obj$BIC
   message("BIC is ", .BIC, "\n")
   message("The output is a scalar.\n")
   return(.BIC)
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the raw count data
raw_read_counts<- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   
   Read_raw <- BUSseqfits_obj$CountData_raw
   
   class(Read_raw) <- "CountData"
   
   message(paste0("The output format is a \"CountData\" object with length",
         " equal to the batch number.\n"))
   message("Each element of the object is the raw read count matrix.\n")
   message("In each matrix, each row represents a gene and each column",
         " correspods to a cell.\n")
   
   return(Read_raw)
   
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obtain the underlying true count data
imputed_read_counts<- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   
   Read_imputed <- BUSseqfits_obj$CountData_imputed
   
   class(Read_imputed) <- "CountData"
   
   message(paste0("The output format is a \"CountData\" object with length",
         " equal to the batch number.\n"))
   message("Each element of the object is the imputed read count matrix.\n")
   message(paste0("In each matrix, each row represents a gene and each column",
         " correspods to a cell.\n"))
   
   return(Read_imputed)
   
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

#obatin the corrected count data
corrected_read_counts<- function(BUSseqfits_obj){
   if(is(BUSseqfits_obj,"BUSseqfits")){
   
   Truereads_list <- BUSseqfits_obj$CountData_imputed
   .B <- BUSseqfits_obj$n.batch
   .nb <- BUSseqfits_obj$n.perbatch
   .N <- BUSseqfits_obj$n.cell
   .G <- BUSseqfits_obj$n.gene
   .K <- BUSseqfits_obj$n.celltype
   .gamma <- BUSseqfits_obj$gamma.est
   .logmu <- BUSseqfits_obj$alpha.est + BUSseqfits_obj$beta.est
   .nu <- BUSseqfits_obj$nu.est
   .delta <- BUSseqfits_obj$delta.est
   .phi <- BUSseqfits_obj$phi.est
   .w <- BUSseqfits_obj$w.est
   
   Truereads <- NULL
   delta_vec <- NULL
   w_vec <- NULL
   for(b in seq_len(.B)){
      Truereads <- cbind(Truereads,Truereads_list[[b]])
      delta_vec <- c(delta_vec, .delta[[b]])
      w_vec <- c(w_vec, .w[[b]])
   }
   
   t.start<-Sys.time()
   message("   correcting read counts...\n")
   
   read_corrected <- Truereads
   cell_index <- 1
   for(b in seq_len(.B)){
      for(i in seq_len(.nb[b])){
      unc_mu <- exp(.logmu[,w_vec[cell_index]] + 
                      .nu[,b] + delta_vec[cell_index])
      p_x<-pnbinom(Truereads[,cell_index],size=.phi[,b], 
                   mu=unc_mu)
      p_xminus1<-pnbinom(Truereads[,cell_index]-1,size=.phi[,b], 
                      mu=unc_mu)
      u <- runif(.G,min=p_xminus1,max=p_x)
      u <- apply(cbind(u,0.9999),1,min)
      cor_mu <- exp(.logmu[,w_vec[cell_index]])
      read_corrected[,cell_index] <- 
        qnbinom(u,size =.phi[,1], mu=cor_mu)
      
      cell_index <- cell_index + 1
      }
   }
   
   Read_corrected <- list()
   cell_index <- 0
   for(b in seq_len(.B)){
      Read_corrected[[b]] <- read_corrected[,cell_index + seq_len(.nb[b])]
      cell_index <- cell_index + .nb[b]
   }
   
   class(Read_corrected) <- "CountData"
   
   t.end<-Sys.time()
   message(paste0("   Correcting read counts takes: ", 
                round(difftime(t.end, t.start, units="mins"), 3), 
                " mins", "\n"))
   
   message(paste0("The output format is a \"CountData\" object with length",
         " equal to the batch number.\n"))
   message("Each element of the object is the corrected read count matrix.\n")
   message(paste0("In each matrix, each row represents a gene and",
         " each column correspods to a cell.\n"))
   
   return(Read_corrected)
   
   }else{
   stop("BUSseqfits_obj must be   a \"BUSseqfits\" object!\n")
   }
}

##############################################################################
# print and summary
##############################################################################
#print BUSseqfits
print.BUSseqfits <- function(x, ...){
   BUSseqfits <- x
   .G <- BUSseqfits$n.gene
   .B <- BUSseqfits$n.batch
   
   
   message("Cell type indicators:\n")
   .w <- BUSseqfits$w.est
   .nb <- BUSseqfits$n.perbatch
   for(b in seq_len(.B)){
   
   message(paste0("Batch ", b, " cells' cell type indicators: ",
                .w[[b]][1],",",.w[[b]][2],",",.w[[b]][3], 
                "... ...\n"))
   }
   message("\n")
   
   message("The estimated location batch effects:\n")
   .nu <- BUSseqfits$nu.est
   for(b in seq_len(.B)){
   message(paste0("    Batch ", b, " location batch effects are: ",
         .nu[1,b],",",.nu[2,b],",",.nu[3,b], "... ...\n"))
   }
   message("\n")
   
   message("The estimated overdispersions:\n")
   .phi <- BUSseqfits$phi.est
   for(b in seq_len(.B)){
   message(paste0("    Batch ", b, " scale batch effects are: ",
         .phi[1,b],",",.phi[2,b],",",.phi[3,b], "... ...\n"))
   }
   message("\n")
}

#summarize BUSseqfits
summary.BUSseqfits <- function(object, ...){
   BUSseqfits <- object
   .G <- BUSseqfits$n.gene
   .B <- BUSseqfits$n.batch
   .K <- BUSseqfits$n.celltype
   .N <- BUSseqfits$n.cell
   
   num_iters <- BUSseqfits$n.iter
   num_burnin <- BUSseqfits$n.burnin
   message(c("B=", .B, " batches\n"))
   message(c("G=", .G, " genes\n"))
   message(c("K=", .K, " cell types\n"))
   message(c("N=", .N, " cells in total\n"))
   message(c("Run ", num_iters," iterations with the first ",num_burnin,
             " iterations as burnin in the MCMC algorithm.\n\n"))
   message(paste0("BUSseqfits is an R list that contains",
             " the following main elements:\n\n"))
   message(paste0("    BUSseqfits$w.est : the estimated cell type indicators,",
          " a list with length equal to B.\n"))
   message(paste0("    BUSseqfits$pi.est : the estimated cell type proportions",
             " across batches, a K by B matrix.\n"))
   message(paste0("    BUSseqfits$gamma.est : the estimated the coefficients",
             " of the logistic regression for the dropout events,",
             " a B by 2 matrix\n"))
   message(paste0("    BUSseqfits$alpha.est : the estimated log-scale baseline",
             " expression levels, a vector with length G.\n"))
   message(paste0("    BUSseqfits$beta.est : the estimated cell type effects,",
          " a G by K matrix.\n"))
   message(paste0("    BUSseqfits$delta.est : the estimated cell-specific",
             " effects, a list with length equal to B.\n"))
   message(paste0("    BUSseqfits$nu.est : the estimated location batch",
                  " effects, a G by B matrix.\n"))
   message(paste0("    BUSseqfits$phi.est : the estimated overdispersion",
             " parameters, a G by B matrix.\n"))
   message("    BUSseqfits$BIC : the BIC, a scalar.\n")
   message(paste0("    BUSseqfits$D.est : the intrinsic gene indicators,",
             " a vector with length N.\n"))
   message("    For more output values, please use \"?BUSseq_MCMC\"\n")
   message("\n")
}

#print CountData
print.CountData <- function(x, ...){
   CountData <- x

   .B <- length(CountData)
   .G <- nrow(CountData[[1]])
   
   message(paste0("There are ", .B, " batches and ", .G, " genes.\n"))
   for(b in seq_len(.B)){
   .nperbatch <- ncol(CountData[[b]])
   message(paste0("Batch ", b, " contains ", .nperbatch, 
                " cells, and their read counts in all genes are: \n"))
   message(paste0("Gene 1: ", CountData[[b]][1,1],", ", CountData[[b]][1,2],
                ", ", CountData[[b]][1,3], ", ... ...\n"))
   message(paste0("Gene 2: ", CountData[[b]][2,1],", ", CountData[[b]][2,2],
                ", ", CountData[[b]][2,3], ", ... ...\n"))
   message(paste0("Gene 3: ", CountData[[b]][3,1],", ", CountData[[b]][3,2],
                ", ", CountData[[b]][3,3], ", ... ...\n"))
   message("    ... ...\n\n")
   }
   message("\n")
   
}

#summarize CountData
summary.CountData <- function(object, ...){

   CountData <- object
   
   .B <- length(CountData)
   .G <- nrow(CountData[[1]])
   
   message(paste0("There are ", .B, " batches and ", .G, " genes.\n"))
   for(b in seq_len(.B)){
   .nperbatch <- ncol(CountData[[b]])
   message(paste0("Batch ", b, " contains ", .nperbatch, " cells."))
   }
   message("\n")
}

########################################################################
# Visualization
########################################################################
#visualize the read counts data by stacking all gene expression matrices 
heatmap_data_BUSseq <- function(CountData_obj, gene_set=NULL, 
                        project_name="BUSseq_heatmap", 
                        image_dir=NULL, color_key_seq=NULL, 
                        image_width=1440, image_height=1080){
   
   if(is(CountData_obj,"CountData")){
   .B <- length(CountData_obj)
   .G <- nrow(CountData_obj[[1]])
   .nb <- rep(NA,.B)
   for(b in seq_len(.B)){
   .nb[b] <- ncol(CountData_obj[[b]])
   }

   if(is.null(gene_set)){
   gene_set <- seq_len(.G)
   }

   #heatmap cell colors
   colfunc <- colorRampPalette(c("grey", "black"))
   #batch colors
   color_batch_func <- colorRampPalette(
     c("#EB4334","#FBBD06","#35AA53","#4586F3"))
   
   color_batch <- color_batch_func(.B)
   
   color_batch2 <- NULL
   
   for(b in seq_len(.B)) {
   color_batch2 <- c(color_batch2, rep(color_batch[b], .nb[b]))
   }
   log1p_mat <- NULL
   for(b in seq_len(.B)){
   log1p_mat <- cbind(log1p_mat, log1p(CountData_obj[[b]]))
   }
   log1p_mat_interest <- log1p_mat[gene_set, ]
   
   if(is.null(color_key_seq)){
   range_data <- range(log1p_mat_interest)
   color_key_seq <- seq(from=floor(range_data[1]) - 0.5, 
                   to=ceiling(range_data[2]) + 0.5, length.out=11)
   }
   
   if(is.null(image_dir)){
   image_dir <- "./image"
   }
   #create the folder
   dir.create(image_dir,showWarnings=FALSE)
   
   png(paste(image_dir,"/",project_name,"_log1p_data.png",sep=""),
      width=image_width, height=image_height)
   heatmap.2(log1p_mat_interest,
         dendrogram="none",#with cluster tree
         Rowv=FALSE, Colv=FALSE,
         labRow=FALSE, labCol=FALSE,
         ColSideColors=color_batch2,
         #RowSideColors=genetype_color,
         col=colfunc(length(color_key_seq)-1),breaks=color_key_seq,
         density.info="histogram",
         hclustfun=function(c)hclust(c,method="average"),
         keysize=0.8, cexRow=0.5,trace="none")#font size
   dev.off()
   }else{
   stop("CountData_obj must be   a \"CountData\" object!\n")
   }
}
