#
# Making a seurat object of Tau THY-21 cortex samples - scRNA-seq 
#


library(Seurat) # Seurat_4.3.0
library(ggplot2) # ggplot2_3.4.2
library(cluster) # cluster_2.1.4
set.seed(1)


files = dir()[grep("DGE",dir())]
rmsum = grep("Summary",files)
files = files[-rmsum]

files
#01WTM_scRNA_MM48_NPG_S1_DGE.txt  02TGM_scRNA_MM48_NPG_S2_DGE.txt
#03WTF_scRNA_MM48_NPG_S3_DGE.txt  04TGF_scRNA_MM48_NPG_S4_DGE.txt
#05WTM_scRNA_MM48_NPG_S9_DGE.txt  06WTM_scRNA_MM48_NPG_S10_DGE.txt
#07WTM_scRNA_MM48_NPG_S11_DGE.txt 08TGM_scRNA_MM48_NPG_S12_DGE.txt
#09TGM_scRNA_MM48_NPG_S1_DGE.txt  10TGM_scRNA_MM48_NPG_S2_DGE.txt
#11WTF_scRNA_MM48_NPG_S3_DGE.txt  12WTF_scRNA_MM48_NPG_S4_DGE.txt
#13WTF_scRNA_MM48_NPG_S5_DGE.txt  14TGM_scRNA_MM48_NPG_S6_DGE.txt
#15TGM_scRNA_MM48_NPG_S7_DGE.txt  16TGM_scRNA_MM48_NPG_S8_DGE.txt
#17WTM_scRNA_MM48_NPG_S5_DGE.txt  18TGM_scRNA_MM48_NPG_S6_DGE.txt
#19WTF_scRNA_MM48_NPG_S7_DGE.txt  20TGF_scRNA_MM48_NPG_S8_DGE.txt 


if(exists("merged"))
	rm(merged)
for(file in files){
	if(!exists("merged")){
		counts <- read.table(paste0(count_path,"/",file), sep="\t", header=T, row.names=1)
		
		print(dim(counts))
		
		# add filename to column names
		colnames(counts) = paste(file,colnames(counts))
		
		
		# create a Seurat object containing the RNA adata
		merged <- CreateSeuratObject(counts = counts, assay = "RNA", project=file, min.cells = 3, min.features = 200)

	}else{
		counts <- read.table(paste0(count_path,"/",file), sep="\t", header=T, row.names=1)
				
		colnames(counts) = paste(file,colnames(counts))
		
		sobj = CreateSeuratObject(counts = counts, assay = "RNA", project=file, min.cells = 3, min.features = 200)
		merged <- merge(merged,sobj)
	}
}

#
# conditions from file "Tau22scRNAseq_sampels.xlsx" (23rd Nov. 2021) from Kamil
#
#Sample used	Digital unique names of scRNAseq libraries
#WT M	1WTM
#TG M	2TGM
#WT F	3WTF
#TG F	4TGF
#WT M	5WTM
#WT M	6WTM
#WT M	7WTM
#TG M	8TGM
#TG M	9TGM
#TG M	10TGM
#WT F	11WTF
#WT F	12WTF
#WT F	13WTF
#TG F	14TGM
#TG F	15TGM
#TG F	16TGM
#WT M	17WTM
#TG M	18TGM
#WT F	19WTF
#TG F	20TGF
#

condition_codes = c("WT M","TG M","WT F","TG F","WT M","WT M","WT M","TG M","TG M","TG M","WT F","WT F","WT F","TG F","TG F","TG F","WT M","TG M","WT F","TG F")
samples_codes = c("1WTM","2TGM","3WTF","4TGF","5WTM","6WTM","7WTM","8TGM","9TGM","10TGM","11WTF","12WTF","13WTF","14TGM","15TGM","16TGM","17WTM","18TGM","19WTF","20TGF")

pmatch(samples_codes, gsub("^0","",files))
# [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
# same order as in read files



#
# extract conditions from file names (not adequate here: three files have "TGM" abbreviation, although there are for transgenic females)
#

# extract conditions from condition codes
wt_mice = grep("WT",condition_codes)
tg_mice = grep("TG",condition_codes)
male_mice = grep(" M",condition_codes)
female_mice = grep(" F",condition_codes)



conditions_samples = rep("wt_male",length(files))
conditions_samples[intersect(wt_mice,female_mice)] = rep("wt_female",length(intersect(wt_mice,female_mice)))
conditions_samples[intersect(tg_mice,female_mice)] = rep("tg_female",length(intersect(tg_mice,female_mice)))
conditions_samples[intersect(tg_mice,male_mice)] = rep("tg_male",length(intersect(tg_mice,male_mice)))

table(conditions_samples)
#conditions_samples
#tg_female   tg_male wt_female   wt_male 
#        5         5         5         5 


# use merged colnames for assignmend


colcodes = sapply(gsub("^0","",colnames(merged)), function(x) strsplit(x, "_")[[1]][1])
mapids = match(colcodes, samples_codes)

# check for missing values
any(is.na(mapids))
#[1] FALSE

# assign conditions to individual cells
conditions_cells = condition_codes[mapids]
table(conditions_cells)
# TG F  TG M  WT F  WT M 
#11494 12047 12238  9131

# Add a label to the data containing the conditions
merged@meta.data$conditions = conditions_cells


save(merged, file="~/data/thy21_tau_mouse_cortex_scrnaseq.Rdata")



#
# Quality control and pre-processing
#



# Percentage mitochondrial genes


merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")

summary(merged@meta.data$percent.mt)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      0       0       0       0       0       0

# no mitochondrial contamination



pdf("~/figures/violinplot.pdf")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# We filter cells that have unique feature counts over 7000 or less than 200
# We filter cells that have >5% mitochondrial counts
dat_filt <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)


# use new pre-processing approach
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
sctrans2 = SCTransform(dat_filt, vst.flavor = "v2")



#
# Clustering
#

sctrans2 <- RunPCA(sctrans2, npcs = 30, verbose = FALSE) # only 30 PCs used in next analysis

sctrans2 <- RunUMAP(sctrans2, reduction = "pca", dims = 1:30, verbose = FALSE) # reduction = "pca" is default parameter

sctrans2 <- FindNeighbors(sctrans2, reduction = "pca", dims = 1:30, verbose = FALSE)
sctrans2 <- FindClusters(sctrans2, resolution = 0.8, verbose = FALSE) # resolution, Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities. default = 0.8
# consider adjusting the resolution parameter later for the clustering to obtain a smaller number of communities

pdf("~/figures/dimplot_umap_v2.pdf")
DimPlot(sctrans2, label = TRUE) + NoLegend()
dev.off()


#
# Evaluate clustering - choose optimal number of clusters automatically
#


# Elbow plot: a ranking of principle components based on the percentage of variance explained by each one
pdf("~/figures/ElbowPlot_PC30_mm_Cortex.pdf")
ElbowPlot(sctrans2)
dev.off()

#
# Do not choose the number of cluster manually from the Elbow plot, but use the following automated approach:
#
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html#:~:text=View%20on%20GitHub-,Elbow%20plot%3A%20quantitative%20approach,the%20majority%20of%20the%20variation.
#
# "A more quantitative approach may be a bit more reliable. We can calculate where the principal components start to elbow by taking the larger value of:
# The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
# The point where the percent change in variation between the consecutive PCs is less than 0.1%.
# We will start by calculating the first metric:"
#

# Determine percent of variation associated with each PC
pct <- sctrans2[["pca"]]@stdev / sum(sctrans2[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
# 25

#The first metric returns PC25 as the PC matching these requirements. Lets check the second metric, which identifies the PC where the percent change in variation between consecutive PCs is less than 0.1%:

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2
#15

#This second metric returns PC15. Usually, we would choose the minimum of these two metrics as the PCs covering the majority of the variation in the data.

# Minimum of the two calculations
pcs <- min(co1, co2)

pcs
# 15

pdf("~/figures/colored_elbow_plot.pdf")
# Create a dataframe with values
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))

# Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()


#
# Alternative: Silhouette width approach
#

# Make silhoutte plots at different clustering resolutions (of seurat) to identify best cluster number/size

reduction <- "pca"
dims <- 1:pcs
dist.matrix <- dist(x = Embeddings(object = sctrans2[[reduction]])[, dims])

# resolution <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
resolution <- c(0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5)

silhoutte_width <- data.frame()
# par(mfrow=c(4,5))
# png("silhoutte_width_plots.png", width=8, height=6, units="in", res=300)
for (i in 1:length(resolution)) {
	print(resolution[i])
	sctrans2 <- FindClusters(sctrans2, resolution = resolution[i])
	clusters <- eval(parse(text = paste0("sctrans2$SCT_snn_res.",resolution[i])))
	sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
	#plot(sil, main=resolution[i]) # silhoutte_5_cluster_res0.01 # Average silhoutte width: 0.63
	
	if(length(table(clusters)) < 2)
		next
	
	print(mean(sil[, "sil_width"])) # 0.5256185
	result <- cbind(resolution[i], length(unique(clusters)), mean(sil[, "sil_width"]))
	colnames(result) <- c("resolution", "clusters", "sil_width")
	silhoutte_width <- rbind(silhoutte_width, result)
}
# dev.off

silhoutte_width
#   resolution clusters sil_width
#1       0.005        3 0.3596401
#2       0.010        5 0.4795727
#3       0.020        6 0.4857564
#4       0.030        7 0.4977688
#5       0.040        7 0.4989802
#6       0.050        9 0.5112953
#7       0.060        9 0.5089220
#8       0.070       10 0.5157153 # max value --> 10 clusters
#9       0.080       10 0.5157000
#10      0.100       10 0.5122275 
#11      0.200       14 0.3443768
#12      0.300       16 0.2890529
#13      0.400       17 0.2906151
#14      0.500       19 0.2699902

pdf("~/figures/silhoutte_width_mm_Cortex.pdf")
ggplot(silhoutte_width, aes(x = clusters, y = sil_width)) +
  geom_line() + geom_point() + scale_x_continuous(breaks = 2:31) + theme_bw()
dev.off()

# NOTE: Best result are at "0.07" resolution: 10 clusters with avg. sil_width = 0.52 (highest)
bestres = 0.07

sctrans2 <- FindClusters(sctrans2, resolution = bestres)

#Maximum modularity in 10 random starts: 0.9813
#Number of communities: 10
#Elapsed time: 14 seconds



# Look at the number of cells in each cluster
table(sctrans2$seurat_clusters) # table(Idents(scNorm))
#     0     1     2     3     4     5     6     7     8     9
# 16769  9204  8605  6153   961   816   748   505   474   347


# Look at cluster IDs of the first 5 cells
head(Idents(sctrans2), 5)

head(sctrans2@meta.data$conditions)
sctrans2$stim <- sctrans2@meta.data$conditions

pdf("~/figures/DimPlot_umap_mm_Cortex.pdf")
DimPlot(sctrans2, reduction = "umap", split.by = "stim")
dev.off()

Idents(sctrans2) = sctrans2$seurat_clusters

png("~/figures/DimPlot_umap_mm_Cortex.png", width=10, height=6, units="in", res=300)
DimPlot(sctrans2, reduction = "umap", split.by = "stim")
dev.off()

sctrans2 <- RunTSNE(sctrans2, reduction = "pca", dims = 1:pcs)
pdf("~/figures/DimPlot_tSNE_mm_Cortex.pdf")
DimPlot(sctrans2, reduction = "tsne", split.by = "stim")
dev.off()

pdf("~/figures/DimPlot_umap_mm_Cortex_nosplit.pdf")
DimPlot(sctrans2, reduction = "umap", label = T, repel = T) + ggtitle("Unsupervised clustering")
dev.off()


save(sctrans2, file="~/data/sctrans2.RData")