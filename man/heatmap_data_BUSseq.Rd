\name{heatmap_data_BUSseq}
\alias{heatmap_data_BUSseq}
\title{
  Draw the Heatmap of the Log-scale Read Count Data for 
  the Output of the \code{BUSseq_MCMC} Function
}
\description{
  Plot the heatmap of the log-scale read count data across multiple
  batches, and then save the resulting images in the user's directory 
  as "png" format. 
}
\usage{
heatmap_data_BUSseq(sce_BUSseqfit, 
                    data_type = c("Raw","Imputed","Corrected"), 
                    gene_set = NULL,
                    project_name= paste0("BUSseq_heatmap_",data_type), 
                    image_dir = NULL, color_key_seq = NULL, 
                    image_width = 1440, image_height = 1080)
}
\arguments{
  \item{sce_BUSseqfit}{An output \code{SingleCellExperiment} object obtained 
  from the function \code{BUSseq_MCMC}.}
  \item{data_type}{
    A string to determine which count data matrix is used to draw the heatmap,
    "Raw" for the raw count data, "Imputed" for the imputed data, and 
    "Corrected" for the corrected data.
  }
  \item{gene_set}{
    A vector of gene indices indicating the gene set of interest to
    display in the heatmap. The default is all genes. 
    We also recommend displaying the intrinsic genes obtained from 
    \code{intrisic_genes_BUSseq(BUSseqfits_obj)}.
  }
  \item{project_name}{
    A string to name the "png" image. By default, the figure is named 
    as "BUSseq_heatmap_Raw_log1p_data.png."
  }
  \item{image_dir}{
    A directory to store the gnereated heatmap. The default is to 
    create a folder called "image" in the current directory and 
    save there.
  }
  \item{color_key_seq}{
    A numeric vector indicating the splitting points for binning 
    log-scale read counts into colors. The default is to space the 
    color key points equally between the minimum and maximum of 
    the log-scale read count data.
  }
  \item{image_width}{
    The width in pixels of the graphical device to plot. 
    The default is 1440 px.
  }
  \item{image_height}{
    The height in pixels of the graphical device to plot. 
    The default is 1080 px.
  }
}
\details{
  To cope with the zeros in the count data, we take the transformation
  log(1+x) on all count data, which corresponds to the R function
  \code{log1p()} instead of \code{log()}.
}
\value{
  Visualize the gene expression data matrix, where each row represents 
  a gene, and each column represents a sample.
}
\references{
Song, Fangda, Ga Ming Angus Chan, and Yingying Wei. 
Flexible experimental designs for valid single-cell RNA-sequencing 
experiments allowing batch effects correction. 
Nature communications 11, no. 1 (2020): 1-15.
}
\author{
  Fangda Song
}
\examples{
  library(SingleCellExperiment)
  # Plot the imputed read count data of the first 100 genes
  heatmap_data_BUSseq(BUSseqfits_example, data_type = "Imputed",
                      gene_set = 1:100)
}
