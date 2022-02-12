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

#3--------------------------------------------

CEL_data <- ReadAffy(celfile.path = "/projectnb/bf528/users/dreadlocks/project_1/samples")
norm_CEL_data <- rma(CEL_data)

#4--------------------------------------------

fit_data <- fitPLM(CEL_data, normalize=TRUE, background=TRUE)

RLE_median <- RLE(fit_data, type="stats")["median",]
jpeg(file='outputs/RLE_median.jpeg')
hist(RLE_median, xlab = 'Median RLE', main = 'Median of RLE Values')
dev.off()

NUSE_median <- NUSE(fit_data, type="stats")["median",]
jpeg(file='outputs/NUSE_median.jpeg')
hist(NUSE_median, xlab = 'Median Nuse', main = 'Median of Nuse Values')
dev.off()

#5--------------------------------------------

proj_metadata <- read.csv('/project/bf528/project_1/doc/proj_metadata.csv')
mod <- model.matrix(~as.factor(proj_metadata$normalizationcombatmod))
batch_corrected_data <- ComBat(exprs(norm_CEL_data), batch=proj_metadata$normalizationcombatbatch, mod=mod[,2])
write.csv(batch_corrected_data,'/projectnb/bf528/users/dreadlocks/project_1/normalized_adjusted_expression_values.csv')

#6--------------------------------------------

scaled_data <- t(scale(t(batch_corrected_data)))
pca_data <- prcomp(scaled_data, center = FALSE, scale = FALSE)

#7--------------------------------------------

percent_variance <- round(pca_data$sdev^2/sum(pca_data$sdev^2)*100, 2)
plot <- ggplot(pca_data$x[,1:2], aes(x=PC1,y=PC2)) +
  geom_point(aes(x=pca_data$x[,1], y=pca_data$x[,2]), colour = "maroon", size = 0.1) +
  labs(x=paste0("PC1 (",percent_variance[1],"%)"),
       y=paste0("PC2 (",percent_variance[2],"%)")) +
  ggtitle('PCA plot')
jpeg(file='outputs/PCA_plot.jpeg')
ggplotly(plot)
dev.off()


