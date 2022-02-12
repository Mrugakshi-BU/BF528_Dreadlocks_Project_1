if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("sva")
BiocManager::install("AnnotationDbi")
BiocManager::install("hgu133plus2.db")
install.packages("ggfortify")
install.packages("tidyverse")
install.packages("plotly")
#--------------------------------------------
library("affy") ##load the affy package
library("affyPLM")
library("sva")
library("AnnotationDbi")
library("hgu133plus2.db")
library("ggfortify")
library("tidyverse")
library("plotly")

#1--------------------------------------------

CEL_data <- ReadAffy(celfile.path = "/projectnb/bf528/users/dreadlocks/project_1/samples")
norm_CEL_data <- rma(CEL_data)

#2--------------------------------------------

fit_data <- fitPLM(CEL_data, normalize=TRUE, background=TRUE)

RLE(fit_data)
NUSE(fit_data)

RLE_stats <- RLE(fit_data, type="stats")
RLE_median <- RLE_stats["median",]
hist(RLE_median)

NUSE_stats <- NUSE(fit_data, type="stats")
NUSE_median <- NUSE_stats["median",]
hist(NUSE_median)

#3--------------------------------------------

proj_metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
mod <- model.matrix(~as.factor(proj_metadata$normalizationcombatmod))

batch_corrected_data <- ComBat(exprs(norm_CEL_data), batch=proj_metadata$normalizationcombatbatch, mod=mod[,2])
write.csv(batch_corrected_data,'normalized_expression_matrix.csv', sep = ' ')

#my output has extra columns 100 columns in the start that the sample output doesn't. 
#The sample output only has the last 34 cols of my output.

#4--------------------------------------------

scaled_data <- t(scale(t(batch_corrected_data)))
pca_data <- prcomp(scaled_data, center = FALSE, scale = FALSE)

#5--------------------------------------------

percent_variance <- round(pca_data$sdev^2/sum(pca_data$sdev^2)*100, 2)
plot <- ggplot(pca_data$x[,1:2], aes(x=PC1,y=PC2)) + 
  geom_point(aes(x=pca_data$x[,1], y=pca_data$x[,2]), colour = "maroon", size = 0.1) +
  labs(x=paste0("PC1 (",percent_variance[1],"%)"),
       y=paste0("PC2 (",percent_variance[2],"%)")) +
  ggtitle('PCA plot')
ggplotly(plot)

