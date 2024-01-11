# Packages
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(glmpca)
library(ggbeeswarm)
library(apeglm)
library(IHW)
library(genefilter)


# Load data into DEseq2
countdata <- read.table("Gacu_gut_counts.tsv", sep = '\t', header = TRUE, row.names = 1)
coldata <- read.table("Gacu_gut_metadata.tsv", sep = '\t', header = TRUE, stringsAsFactors = TRUE)

ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = coldata,
                                 design = ~ Population + Sex + Treatment)

dds <- DESeq(ddsMat)


# Pre-filter the dataset
nrow(dds) # 22456

keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds) # 20797


# rlog Transformation
rld <- rlog(dds, blind = TRUE)
head(assay(rld), 3)
colData(rld)

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(assay(rld)[,1:3]) %>%
    mutate(transformation = "rlog")
)

colnames(df)[1:3] <- c("x", "y", "z")

lvls <- c("rlog")
df$transformation <- factor(df$transformation, levels = lvls)

ggplot(df, aes(x=x, y=y)) + 
  geom_hex(bins=80) +
  coord_fixed() +
  facet_grid(. ~ transformation)


# Generate sample distance heatmaps with dendrograms

# Euclidean distance
sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Population, rld$Sex, rld$Treatment, sep = '-')
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         cllustering_distance_rows = sampleDists,
         col = colors)

# Poisson Distance
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(dds$Population, dds$Sex, dds$Treatment, sep = '-')
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


# PCA plot
plotPCA(rld, intgroup = c("Population", "Sex", "Treatment"))

# PCA plot using generalized PCA
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$pop <- dds$Population
gpca.dat$treatment <- dds$Treatment
gpca.dat$sex <- dds$Sex
gpca.dat$all <- paste(dds$Sex, dds$Treatment, dds$Population)

ggplot(gpca.dat, aes(x=dim1, y=dim2, color=all)) +
  geom_point(size=3) +
  coord_fixed() +
  ggtitle("glmcpa - Generalized PCA")


# MDS plot
mds <- as.data.frame(colData(rld)) %>%
  cbind(cmdscale(sampleDistMatrix))

ggplot(mds, aes(x= `1`, y= `2`, color = Population:Sex:Treatment)) +
  geom_point(size=3) + coord_fixed() + ggtitle("MDS with rld data")


# Results table
res <- results(dds)
res

mcols(res, use.names = TRUE)
summary(res)

res.01 <- results(dds, alpha = 0.01)
table(res.01$padj < 0.05)

res.s <- results(dds, alpha = 0.05)
table(res.s$padj < 0.1)

sum(res$pvalue < 0.1, na.rm = TRUE)


# Count plot
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup = c("Treatment"))

geneCounts <-plotCounts(dds, gene = topGene, intgroup = c("Sex", "Treatment"),
                        returnData = TRUE)
ggplot(geneCounts, aes(x=Sex, y=count, color=Treatment)) +
  scale_y_log10() +
  geom_beeswarm(cex=3)

ggplot(geneCounts, aes(x=Sex, y=count, color=Treatment, group=Treatment)) +
  scale_y_log10() + 
  geom_point(size=3) +
  geom_line()


# MA-plot
resultsNames(dds)

res <- lfcShrink(dds, coef = "Treatment_GF_vs_CV", type = "apeglm")
plotMA(res, ylim=c(-5,5))

res.noshr <- results(dds, name="Treatment_GF_vs_CV")
plotMA(res.noshr, ylim=c(-5,5))

plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


# Gene clustering
topvarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)

mat <- assay(rld)[topvarGenes,]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("Treatment", "Sex")])
pheatmap(mat, annotation_col = anno)