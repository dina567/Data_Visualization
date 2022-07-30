# ----Plotting PCA (Principal Component Analysis) with probability ellipse or k-mean cluster ------ 

# Method 1 - DESeq2:plotPCA() for DESeq Dataset -----

# More info: https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/plotPCA

# Example:
# DEseq2 Data preparation
BiocManager::install("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)
setwd("C:/Rutgers/Github/Data_Visualization")
data <- read.csv("featurecounts.example.csv", sep = ",", header = TRUE, row.names=1)
data <- data[,c(6:35)] 
colnames(data) <- c(1:30)
cts <-as.matrix(data)

coldata <- read.csv("design table_RNA-seq-07-05-22.csv", sep = ",", header = TRUE, row.names = 1)
coldata <- coldata[,c("trt","wk")] # only keep the columns of conditions by which you want to compare
coldata$trt <- factor(coldata$trt)
coldata$wk <- factor(coldata$wk)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ wk + trt + wk:trt) #wk:trt interaction term
rld <- rlog(dds) #rlog transformation

# (1) plot PCA directly with plotPCA() 
pcaData <- plotPCA(rld, intgroup=c("trt", "wk"), returnData=TRUE)


# (2) plot PCA in conjunction with ggplot2
# add probability ellipse using stat_ellipse() (can specify grouping for ellipse)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$trt <- factor(pcaData$trt, levels = c("Control","BAP+TPA", "BAP+TPA+UA"))
ggplot(pcaData, aes(PC1, PC2, color=trt, shape=wk)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+
  theme(aspect.ratio=1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size =1.5),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size =10, face = "bold")) +
  scale_shape_manual(values = c(16, 17, 15, 10)) +
  stat_ellipse(aes(PC1, PC2, group = trt),type = "norm")


# Method 2 - ggfortify:autoplot()

library(ggplot2)
library(ggfortify)
library(cluster)

# add probability ellipse by pam(data, # of clustering)

autoplot(pam(pcaData, 3),frame = TRUE, , frame.type = 'norm')

# add k-means clustering by fanny(data, # of clustering)
# *background info from wiki: k-means clustering is a method of vector quantization:
# aims to partition n observations into k clusters in which each observation belongs to the cluster with the nearest mean (cluster centers or cluster centroid), serving as a prototype of the cluster. This results in a partitioning of the data space into Voronoi cells. 
# k-means clustering minimizes within-cluster variances (squared Euclidean distances), but not regular Euclidean distances

autoplot(fanny(pcaData, 3),frame = TRUE)

