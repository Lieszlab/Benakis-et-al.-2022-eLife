

# 10x scRNA-sequencing MG-Tc project #

## Title: T cells modulate the microglial response to brain ischemia
## Authors: Corinne Benakis, Alba Simats, Sophie Tritschler, Steffanie Heindl, Simon Besson-Girard,Gemma Llovera, Kelsey Pinkham, Anna Kolz, Alessio Ricci, Fabian Theis, Stefan Bittner, Özgün Gökce, Anneli Peters, Arthur Liesz



library(SeuratDisk)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(airway)
library(magrittr)
library(EnhancedVolcano)

# Load data and create seurat objects
Raw_data117 <- Read10X(data.dir = "/Users/Asimatso/Documents/raw_filtered_count_matrices/117-S-WT")
sample117 <- CreateSeuratObject(counts = Raw_data117, project = "117", min.cells = 3, min.features = 200)

Raw_data118 <- Read10X(data.dir = "/Users/Asimatso/Documents/raw_filtered_count_matrices/118-S-RAG")
sample118 <- CreateSeuratObject(counts = Raw_data118, project = "118", min.cells = 3, min.features = 200)

Raw_data119 <- Read10X(data.dir = "/Users/Asimatso/Documents/raw_filtered_count_matrices/119-N-WT")
sample119 <- CreateSeuratObject(counts = Raw_data119, project = "119", min.cells = 3, min.features = 200)

Raw_data120 <- Read10X(data.dir = "/Users/Asimatso/Documents/raw_filtered_count_matrices/120-N-RAG")
sample120 <- CreateSeuratObject(counts = Raw_data120, project = "120", min.cells = 3, min.features = 200)


# Merge files
MG <- merge(sample117, y = c(sample118,sample119,sample120), add.cell.ids = c("117", "118", "119", "120"), project = "MGTc")

head(MG)
head(MG@meta.data)
table(MG$orig.ident)


# Normalizing the data
MG <- NormalizeData(object = MG, normalization.method = "LogNormalize", scale.factor = 1e4)

# Detection of variable genes across the single cells
MG <- FindVariableFeatures(object = MG, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = MG))

# Scaling the data and removing unwanted sources of variation
MG <- ScaleData(object = MG, features = rownames(x = MG), vars.to.regress = c("nCount_RNA", "percent.mito"))
saveRDS(MG, file = "MG.rds")


# Perform linear dimensional reduction
MG <- RunPCA(object = MG, features = VariableFeatures(object = MG), verbose = FALSE)

# Examine and visualize PCA results a few different ways
print(x = MG[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = MG, dims = 1:2)
DimPlot(object = MG)

# ProjectDim scores each feature in the dataset
MG <- ProjectDim(MG)
DimHeatmap(object = MG, dims = 1:15, cells = 500, balanced = TRUE)

# Find neighbors and clusters
MG <- FindNeighbors(object = MG, dims = 1:5)
MG <- FindClusters(object = MG, resolution = 1)

# Run Non-linear dimensional reduction (tSNE)
MG <- RunTSNE(object = MG, dims = 1:5)
DimPlot(MG, reduction = 'tsne', label = TRUE)
DimPlot(MG, reduction = 'tsne', label = TRUE, group.by = "organ")

# Run UMAP
MG <- RunUMAP(object = MG, dims = 1:5)
DimPlot(MG, reduction = 'umap', label = TRUE, group.by = "seurat_clusters")
DimPlot(MG, reduction = 'umap', label = TRUE, split.by = "group")

# Identify microglia clusters
MGmarkers <- FindAllMarkers(object = MG, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MGmarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table((MGmarkers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)), 'markers_by_cluster.MG.tsv', sep='\t')

FeaturePlot(MG, reduction = 'umap', features = c("Tmem119", "Cx3cr1", "Selplg", "Csf1r","Trem2", "Sall1", "P2ry12", "Apoe", "Lyz2"), cols = c("darkmagenta","darkcyan", "yellow"))

## Add info about the group
### Get batches based on cell names
sample <- sapply(colnames(GetAssayData(object = MG, slot = "counts")),
                 FUN=function(x){substr(x,1,3)})

head(sample)
tail(sample)

### Turn to numbers and add cell names to them
sample <- as.numeric(as.character(sample))
names(sample) <- colnames(GetAssayData(object = MG, slot = "counts"))

class(sample)
head(sample)

MG <- AddMetaData(MG, sample, "group")
head(MG)
tail(MG)


## Select only microglia clusters: clusters 0,1,2,3,8
Idents(MG) <- "seurat_clusters"
MGisolated <- WhichCells(object = MG, idents = c("0","1","2","3","8"))
MGisolated <- subset(x = MG, cells = MGisolated)

table(MGisolated$group)
DimPlot(MGisolated, reduction = 'umap', label = TRUE, group.by = "seurat_clusters")
DimPlot(MGisolated, reduction = 'umap', label = TRUE, group.by = "group", cols = c("blue", "darkblue", "red", "darkred"))
FeaturePlot(MGisolated, reduction = 'umap', features = c("Tmem119", "Cx3cr1", "Selplg", "Csf1r","Trem2", "Sall1", "P2ry12", "Apoe", "Lyz2"), cols = c("darkmagenta","darkcyan", "yellow"))

# Change names of group
MGisolated$group[MGisolated$group == "117"] <- "Stroke WT"
MGisolated$group[MGisolated$group == "118"] <- "Stroke RAG1KO"
MGisolated$group[MGisolated$group == "119"] <- "Naive WT"
MGisolated$group[MGisolated$group == "120"] <- "Naive RAG1KO"

DimPlot(MGisolated, reduction = 'umap', label = TRUE, group.by = "group")
DimPlot(MGisolated, reduction = 'umap', label = FALSE, group.by = "group", cols = c("blue", "darkblue", "red", "darkred"))
DimPlot(MGisolated, reduction = 'umap', label = TRUE, split.by = "group", cols = c("blue", "darkblue", "red", "darkred"))



# Normalization, scaling and recluster again only with MGisolated
MGisolated <- NormalizeData(object = MGisolated, normalization.method = "LogNormalize", scale.factor = 1e4)

## Detection of variable genes across the single cells
MGisolated <- FindVariableFeatures(object = MGisolated, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = MGisolated))

## Scaling the data and removing unwanted sources of variation
MGisolated <- ScaleData(object = MGisolated, features = rownames(x = MGisolated), vars.to.regress = c("nCount_RNA", "percent.mito"))

## Perform linear dimensional reduction
MGisolated <- RunPCA(object = MGisolated, features = VariableFeatures(object = MGisolated), verbose = FALSE)

## Examine and visualize PCA results a few different ways
print(x = MGisolated[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)
VizDimLoadings(object = MGisolated, dims = 1:2)
DimPlot(object = MGisolated)
DimHeatmap(object = MGisolated, dims = 1:15, cells = 500, balanced = TRUE)

## Find neighbors and clusters
MGisolated <- FindNeighbors(object = MGisolated, dims = 1:5)
MGisolated <- FindClusters(object = MGisolated, resolution = 0.8)

## Run Non-linear dimensional reduction (tSNE)
MGisolated <- RunTSNE(object = MGisolated, dims = 1:5)
DimPlot(MGisolated, reduction = 'tsne', label = TRUE)
DimPlot(MGisolated, reduction = 'tsne', label = TRUE, group.by = "organ")

## Run UMAP
MGisolated <- RunUMAP(object = MGisolated, dims = 1:5)
DimPlot(MGisolated, reduction = 'umap', label = TRUE, group.by = "seurat_clusters")
DimPlot(MGisolated, reduction = 'umap', label = TRUE, split.by = "group")
DimPlot(MGisolated, reduction = 'umap', label = FALSE, group.by = "group", cols = c("blue", "darkblue", "red", "darkred"))

# Compute DEG between conditions
Idents(object = MGisolated) <- "group"
MGisolated_RAGstroke_WTstroke <- FindMarkers(object = MGisolated, ident.1 = "Stroke RAG1KO", ident.2 = "Stroke WT", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(MGisolated_RAGstroke_WTstroke, 'MGisolated_RAGstroke_WTstroke.tsv', sep='\t')

MGisolated_RAGnaive_WTnaive <- FindMarkers(object = MGisolated, ident.1 = "Naive RAG1KO", ident.2 = "Naive WT", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(MGisolated_RAGnaive_WTnaive, 'MGisolated_RAGnaive_WTnaive.tsv', sep='\t')

MGisolated_RAGStroke_RAGnaive <- FindMarkers(object = MGisolated, ident.1 = "Stroke RAG1KO", ident.2 = "Naive RAG1KO", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(MGisolated_RAGStroke_RAGnaive, 'MGisolated_RAGStroke_RAGnaive.tsv', sep='\t')

MGisolated_WTStroke_WTnaive <- FindMarkers(object = MGisolated, ident.1 = "Stroke WT", ident.2 = "Naive WT", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(MGisolated_WTStroke_WTnaive, 'MGisolated_WTStroke_WTnaive.tsv', sep='\t')

# Re-check quality of the data #
## Mitocondrial genes
mito.features <- grep(pattern = "^mt-", x = rownames(x = MGisolated), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = MGisolated, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = MGisolated, slot = 'counts'))
mito.features

MGisolated[['percent.mito']] <- percent.mito
VlnPlot(object = MGisolated, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, group.by = "orig.ident")

## FeatureScatter 
plot1 <- FeatureScatter(object = MGisolated, feature1 = "nFeature_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(object = MGisolated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# save PROCESSED DATA #
saveRDS(MG, file = "MG.rds")
saveRDS(MGisolated, file = "MGisolated.rds")

MGisolated <- readRDS("MGisolated.rds")
MG <- readRDS("MG.rds")

# Volcano plot libraries


## MGisolated_WTStroke_WTnaive
EnhancedVolcano(MGisolated_WTStroke_WTnaive,
                lab = rownames(MGisolated_WTStroke_WTnaive),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'MGisolated_WTStroke_WTnaive',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

### Confirm genes down in Stroke vs Naive WT
RidgePlot(MGisolated, features = c("Spp1", "Malat1"), group.by = "group")+ NoLegend()
DotPlot(MGisolated, features = c("Spp1", "Apoe", "S100a9", "Malat1"), group.by = "group")+ NoLegend()

## MGisolated_RAGStroke_RAGnaive
EnhancedVolcano(MGisolated_RAGStroke_RAGnaive,
                lab = rownames(MGisolated_RAGStroke_RAGnaive),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'MGisolated_RAGStroke_RAGnaive',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

## MGisolated_RAGnaive_WTnaive
EnhancedVolcano(MGisolated_RAGnaive_WTnaive,
                lab = rownames(MGisolated_RAGnaive_WTnaive),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'MGisolated_RAGnaive_WTnaive',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)

RidgePlot(MGisolated, features = c("Spp1", "Lgals3", "S100a9", "Apoe"), group.by = "group")+ NoLegend()

## MGisolated_RAGstroke_WTstroke
EnhancedVolcano(MGisolated_RAGstroke_WTstroke,
                lab = rownames(MGisolated_RAGstroke_WTstroke),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'MGisolated_RAGstroke_WTstroke',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 1.0,
                labSize = 6.0,
                col=c('black', 'black', 'black', 'red3'),
                colAlpha = 1)


# Volcano plots for publication
EnhancedVolcano(MGisolated_RAGstroke_WTstroke_2,
                lab = rownames(MGisolated_RAGstroke_WTstroke_2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                title = "MGisolated_RAGstroke_WTstroke_2",
                labSize = 6.0,
                xlim = c(-2.5,2.5),
                ylim = c(0, 30),
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                arrowheads = FALSE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 40)

EnhancedVolcano(MGisolated_RAGnaive_WTnaive_2,
                lab = rownames(MGisolated_RAGnaive_WTnaive_2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0,
                title = "MGisolated_RAGnaive_WTnaive_2",
                xlim = c(-2,2),
                ylim = c(0, 30),
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                arrowheads = FALSE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 50)


# Heatmap per condition
Idents(object = MGisolated) <-  "group"
MGisolated.markers.group <- FindAllMarkers(object = MGisolated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MGisolated.markers.group %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

MGisolated.markers.group %>% group_by(cluster) %>% top_n(n=10, avg_log2FC) -> top10
heatmap <- DoHeatmap(object = MGisolated, features = top10$gene, group.colors = brewer.pal(n=9, "Paired")) + scale_fill_gradient(low="yellow", high="darkblue")
heatmap

















