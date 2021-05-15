\name{BUSseqfits_example}
\docType{data}
\alias{BUSseqfits_example}
\title{An external example of the BUSseqfits object}
\description{
  This data set is a BUSseqfits object obtained from 
  the simulatied data in the "Example" of \code{BUSseq-package}.
}
\usage{BUSseqfits_example}
\details{
  The simuated count data consist of two 150-cell batches. 
  In total, 300 genes are measured. 
  Moreover, all cells come from 4 cell types.
}
\format{A BUSseqfits object.}
\examples{
\donttest{
################################
# Set the Synthetic Parameters #
################################
rm(list=ls())
set.seed(1234)
#The number of batches
B<-2

#The number of cells per batch
nb<-c(150,150)

#The total number of cells
N<-sum(nb)

#The number of genes
G<-300

#The number of cell types
K<-4

# the project name
proj <- "demo"

#The first column of gamma.syn denotes the intercept of 
#the logistic regression for dropout events
#The second column of gamma.syn denotes the odds ratios 
#of the logistic regression for dropout events
gamma.syn<-matrix(0,B,2)
gamma.syn[1,]<-c(-1,-0.5)
gamma.syn[2,]<-c(-1,-0.5)

#the log-scale baseline expression levels
alpha.syn<-rep(0,G)
alpha.syn[1:(G*0.2)]<-2

#the cell-type effects 
beta.syn<-matrix(0,G,K)

#the first cell type is regarded as the reference cell type
beta.syn[,1] <- 0

#the effects of the second cell type
beta.syn[1:(G * 0.2),2] <- -2
beta.syn[(G * 0.20 + 1):(G * 0.4),2] <- 2
beta.syn[(G * 0.4 + 1):G,2] <- 0

#the effects of the third cell type
beta.syn[1:(G * 0.2),3]<- -2
beta.syn[(G * 0.2 + 1):(G * 0.4),3] <- 0
beta.syn[(G * 0.4 + 1):(G * 0.6),3] <- 2
beta.syn[(G * 0.6 + 1):G,3] <- 0

#the effects of the forth cell type
beta.syn[1:(G * 0.2),4]<- -2
beta.syn[(G * 0.2 + 1):(G * 0.6),4] <- 0
beta.syn[(G * 0.6 + 1):(G * 0.8),4] <- 2
beta.syn[(G * 0.8 + 1):G,4] <- 0

#the batch effects
nu.syn<-matrix(NA,G,B)

#the first batch is regarded as the reference batch
nu.syn[,1] <- 0

#the effect of the second batch
nu.syn[1:(G * 0.4),2]<-2
nu.syn[(G * 0.4 + 1):(G * 0.8),2]<-1
nu.syn[(G * 0.8 + 1):G,2]<-2

#the cell-specific size factors
delta.syn <- list()
for(b in 1:B){
  delta.syn[[b]] <- rep(NA, nb[b])
}

#The frist cell in each bathc is regarded as the reference cell
delta.syn[[1]][1:(nb[1] * 0.5)]<-0
delta.syn[[1]][(nb[1] * 0.5 + 1):(nb[1] * 0.9)]<-1
delta.syn[[1]][(nb[1] * 0.9 + 1):nb[1]]<-2

#the second batch
delta.syn[[2]][1:(nb[2] * 0.5)]<-0
delta.syn[[2]][(nb[2] * 0.5 + 1):(nb[2] * 0.7)]<-2
delta.syn[[2]][(nb[2] * 0.7 + 1):nb[2]]<--1

#the batch-specific and gene-specific overdispersion parameters
phi.syn<-matrix(10,B,G)

#the cell-type proportions in each batch
pi.syn <- matrix(NA,K,B)

#the frist batch
pi.syn[,1]<-c(0.3,0.3,0.2,0.2)

#the second batch
pi.syn[,2]<-c(0.2,0.24,0.2,0.36)


##############################################
# Simulate Latent Varibles and Observed data #
##############################################
#the cell-type indicators of each cell
w<-list()

#the first batch
nw<-nb[1] * pi.syn[,1]
w[[1]]<- rep(1:4,nw)

#the second batch
nw<-nb[2] * pi.syn[,2]
w[[2]]<- rep(1:4,nw)

#the indicators for dropout events
z<-list()

#the underlying true expression levels
x<-list()

#the observed expression levels
y<-list()

#the logarithm of mean expreesion level of each gene in each cell
log.mu<-list()

for(b in 1:B){
  z[[b]] <- matrix(NA, G, nb[b])
  x[[b]] <- matrix(NA, G, nb[b])
  y[[b]] <- matrix(NA, G, nb[b])
  log.mu[[b]] <- matrix(NA, G, nb[b])
}

#generate the latent variable and observed data
for(b in 1:B){
  for(i in 1:nb[b]){
    log.mu[[b]][,i] <- alpha.syn + beta.syn[,w[[b]][i]] 
    log.mu[[b]][,i] <- log.mu[[b]][,i] + nu.syn[,b]
    log.mu[[b]][,i] <- log.mu[[b]][,i] + delta.syn[[b]][i]
    
    for(j in 1:G){
      x[[b]][j,i]<-rnbinom(1,phi.syn[b,j],
                           mu=exp(log.mu[[b]][j,i]))
      logit_pi <- gamma.syn[b,1] + gamma.syn[b,2] * x[[b]][j,i]
      
      z[[b]][j,i]<-rbinom(1,1,prob = 1/(1+exp(-logit_pi)))
      if(z[[b]][j,i]==0){
        y[[b]][j,i]<- x[[b]][j,i]
      }else{
        y[[b]][j,i]<- 0
      }
    }
  }
}
BUSseqfits_example <- BUSseq_MCMC(ObservedData = y, 
                                  seed = 1234, n.cores = 2,
                                  n.celltypes = 4, n.iterations = 500)

}
}
\source{The "Example" of \code{BUSseq-package}}
\references{
Song, Fangda, Ga Ming Angus Chan, and Yingying Wei. 
Flexible experimental designs for valid single-cell RNA-sequencing 
experiments allowing batch effects correction. 
Nature communications 11, no. 1 (2020): 1-15.
}
\keyword{datasets}