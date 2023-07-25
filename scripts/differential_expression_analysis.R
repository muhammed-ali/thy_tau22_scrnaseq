#
# Check Marker Genes 
#


library(Seurat) # Seurat_4.3.0
library(ggplot2) # ggplot2_3.4.2
library(cluster) # cluster_2.1.4
library(tidyverse) # tidyverse_2.0.0
library(dplyr) # dplyr_1.1.2
library(HGNChelper) # HGNChelper_0.8.1
library(openxlsx) # openxlsx_4.2.5.2

set.seed(1)

load("~/data/sctrans2.RData")


DefaultAssay(sctrans2) <- "SCT"

cluster_markers <- FindAllMarkers(sctrans2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.05)

table(cluster_markers$cluster)
#   0   1   2   3   4   5   6   7   8   9
# 180 157 136 159  82 116 196 165 198 140



cluster_markers %>%
    group_by(cluster) %>%
    slice_max(n = 3, order_by = avg_log2FC) %>% print(n = 30)
# Groups:   cluster [10]
#        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene
#        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>
#  1 0              1.24  0.956 0.104 0         0       Ctss
#  2 0              1.14  0.934 0.098 0         0       Hexb
#  3 0              1.06  0.896 0.089 0         0       Cx3cr1
#  4 0              1.17  0.875 0.043 0         1       Flt1
#  5 0              1.14  0.888 0.147 0         1       Bsg
#  6 0              0.874 0.811 0.134 0         1       Sptbn1
#  7 0              1.29  0.977 0.201 0         2       Atp1a2
#  8 0              1.28  0.962 0.129 0         2       Slc1a2
#  9 0              1.17  0.919 0.075 0         2       Plpp3
# 10 0              2.10  1     0.191 0         3       Plp1
# 11 0              1.15  0.926 0.054 0         3       Mbp
# 12 0              1.11  0.908 0.071 0         3       Mal
# 13 0              1.77  0.879 0.131 0         4       Ttr
# 14 0              0.955 0.565 0.126 0         4       Enpp2
# 15 0              0.445 0.327 0.009 0         4       1500015O10Rik
# 16 0              0.998 0.65  0.02  0         5       Meg3
# 17 4.33e-282      0.535 0.564 0.142 8.09e-278 5       Ttc3
# 18 0              0.491 0.346 0.005 0         5       Syt1
# 19 0              0.966 0.699 0.05  0         6       Ccnd2
# 20 0              0.753 0.634 0.008 0         6       Sox11
# 21 0              0.697 0.642 0.092 0         6       Sox4
# 22 1.99e-203      1.00  0.861 0.416 3.72e-199 7       Apoe
# 23 0              0.987 0.77  0.055 0         7       Lyz2
# 24 0              0.836 0.699 0.005 0         7       Mrc1
# 25 0              1.20  0.962 0.23  0         8       Dbi
# 26 0              1.13  0.909 0.04  0         8       Nnat
# 27 0              1.01  0.943 0.299 0         8       Hsp90aa1
# 28 0              1.15  0.72  0.033 0         9       Mgp
# 29 3.09e-139      0.821 0.568 0.127 5.78e-135 9       Ptgds
# 30 0              0.784 0.55  0.022 0         9       Igfbp5


cluster_markers_top10 <- cluster_markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) %>% print(n = 100)
write.csv(as.data.frame(cluster_markers_top10), file="cluster_markers_top10.txt", quote=F)


cluster_markers %>%
    group_by(cluster) %>%
    top_n(n = 4, wt = avg_log2FC) -> top4

png("~/figures/Featureplot_top_marker_genes.png", width=12, height=14, units="in", res=300)
FeaturePlot(sctrans2, features = top4$gene, 
  label = TRUE,
  label.size = 4,
  label.color = "red",)
dev.off()


# remove genes that do not exist in scale.data as it will give error gene not found
gene_panel_intersect <- intersect(top4$gene, rownames(sctrans2@assays$SCT@scale.data))
length(gene_panel_intersect)

# Down-sample the cells because with all cells, heatmap does not work
# https://github.com/satijalab/seurat/issues/2724
png("~/figures/Heatmap_top_marker_genes.png", width=7, height=7, units="in", res=300)
DoHeatmap(subset(sctrans2, downsample = 100), 
          features = gene_panel_intersect, 
                  assay = "SCT")
dev.off()


# Verify Cluster annotation based on cluster marker genes
# Aim: Manually go through marker genes and find in which cell they are uniquely expressed and then label the cell types accordingly.
# http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.jsp
# https://panglaodb.se/markers.html
# http://mousebrain.org/adolescent/celltypes.html

# cluster0 = Ctss, Hexb, Cx3cr1, C1qb, C1qa, P2ry12 = Microglial cells
# cluster1 = Flt1, Ly6c1, Cldn5, Abcb1a = Endothelial cells
# cluster2 = Clu, Aldoc, Plpp3 = Astrocytes
# cluster3 = Cldn11, Plp1, Mbp, Ermn = Oligodendrocytes
# cluster4 = Ttr, Enpp2, 1500015O10Rik, Chchd10, Cox6c, Igf2, Cox8a, Stk39, Ndufa4, Arl6ip1 = 
# Enpp2 and Stk39 are Oligo, Igf2 is Mural, Chchd10 is ganglion neuron markers - Confusing but more confidence in Oligo.
# cluster5 = Meg3, Ttc3, Syt1, Ptprz1, Rtn1, Map1b, Snhg11 = Neurons 
# cluster6 = Ccnd2, Dlx1, Dlx6os1, Sox11, Sox4 = Neuroblasts
# cluster7 = Apoe, Lyz2, Mrc1, Pf4, Dab2, Cd74, Maf, Ccl7 = Macrophage
# cluster8 = Ccdc153, Gm973, Rarres2, Enkur, Tmem212 = Ependymal
# cluster9 = Mgp and Cald1 (mural), Ptgds (oligo), Apod (oligo), Ptn (astro) = 2 strong Mural cell marker - overall a mixed cluster.


# Automated annotation by SCType
# as.data.frame(sctype_scores_new)
#    cluster                           type    scores ncells
# 1        8                 Ependymal cell  9206.912    474
# 2        1               Endothelial cell 77221.396   9204
# 3        5                         Neuron  7764.588    816
# 4        3                Oligodendrocyte 95210.435   6153
# 5        2                      Astrocyte 62548.640   8605
# 6        9                     Mural cell  1445.235    347
# 7        6                     Neuroblast  6373.149    748
# 8        0                Microglial cell 90251.488  16769
# 9        7                     Macrophage  3508.650    505
# 10       4 Oligodendrocyte precursor cell  1185.298    961



#
# Automated cell-type mapping and annotation
#


# download SCType source code
#system('wget https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R')
#system('wget https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R')

# load downloaded source
source("~/data/gene_sets_prepare.R");
source("~/data/sctype_score_.R")

# DB file
# system('wget https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx')

db_ = "~/data/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain


# Load prepared mouse annotations from Cell Marker database
load(file="~/data/cellmarker_scrnaseq.Rdata")


#
# Adjusted version of the sctype_score function (avoiding human gene format in upper case letter format for compatibility with mouse MGI symbols)
#

sctype_score_noupper = function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){

  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }

  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T);
 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                      gene_ = names(marker_stat), stringsAsFactors = !1)
	# higher number of occurrence = lower sensitivity


  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]

  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }

  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]

  # combine scores (across marker genes per cell type)
  es = do.call("rbind", lapply(names(gs), function(gss_){
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

  es.max
}


gspositive_lst = sapply(cellmarker_scrnaseq$Cell.Marker, function(x) strsplit(x,", ")[[1]])
names(gspositive_lst) = cellmarker_scrnaseq$Cell.Type

# r merge named vectors with duplicate names
# https://stackoverflow.com/questions/57020599/r-list-combine-elements-with-same-name
gspositive_combined = tapply(unlist(gspositive_lst, use.names = FALSE), rep(names(gspositive_lst), lengths(gspositive_lst)), FUN = c)

# make unique
gspositive_combined  = sapply(gspositive_combined,function(x) unique(x))

es.max_new = sctype_score_noupper(scRNAseqData = sctrans2[["SCT"]]@scale.data, scaled = TRUE, gs = gspositive_combined, gs2 = NULL)
dim(es.max_new)
# 36 44582

cL_resutls_new = do.call("rbind", lapply(unique(sctrans2@meta.data$seurat_clusters), function(cl){
		#print(cl)
    es.max.cl = sort(rowSums(es.max_new[ ,rownames(sctrans2@meta.data[sctrans2@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sctrans2@meta.data$seurat_clusters==cl)), 10)
}))

dim(cL_resutls_new)
# 100 4

# show best annotations and corresponding scores for cluster 0
head(cL_resutls_new[which(cL_resutls_new$cluster == 0),], 10)


# filter clusters with too low scores as "unknown"
sctype_scores_new = cL_resutls_new %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

sctype_scores_new$type[as.numeric(as.character(sctype_scores_new$scores)) < sctype_scores_new$ncells/4] = "Unknown"


# show overview of cluster assignments
print(sctype_scores_new[,1:3])
#   cluster type                           scores
#   <fct>   <chr>                           <dbl>
# 1 8       Ependymal cell                  9207.
# 2 1       Endothelial cell               77191.
# 3 5       Neuron                          7776.
# 4 3       Oligodendrocyte                95241.
# 5 2       Astrocyte                      62535.
# 6 9       Mural cell                      1432.
# 7 6       Neuroblast                      6361.
# 8 0       Microglial cell                90291.
# 9 7       Macrophage                      3504.
#10 4       Oligodendrocyte precursor cell  1187.


as.data.frame(sctype_scores_new)
#    cluster                           type    scores ncells
# 1        8                 Ependymal cell  9206.912    474
# 2        1               Endothelial cell 77221.396   9204
# 3        5                         Neuron  7764.588    816
# 4        3                Oligodendrocyte 95210.435   6153
# 5        2                      Astrocyte 62548.640   8605
# 6        9                     Mural cell  1445.235    347
# 7        6                     Neuroblast  6373.149    748
# 8        0                Microglial cell 90251.488  16769
# 9        7                     Macrophage  3508.650    505
# 10       4 Oligodendrocyte precursor cell  1185.298    961


# Look at the number of cells in each cluster
table(sctrans2$seurat_clusters)
#     0     1     2     3     4     5     6     7     8     9
# 16769  9204  8605  6153   961   816   748   505   474   347



sctrans2@meta.data$classint = ""
for(j in unique(sctype_scores_new$cluster)){
  cl_type = sctype_scores_new[sctype_scores_new$cluster==j,]; 
  sctrans2@meta.data$classint[sctrans2@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}


# Dimplot with cluster annotation
pdf("~/figures/umap_plot_celltype_annot_2022-11-14.pdf")
DimPlot(sctrans2, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'classint')  
dev.off()


# Dimplot with cluster annotation
png("~/figures/Dimplot_SCtype_cluster_labelling.png", width=7, height=7, units="in", res=300)
DimPlot(sctrans2, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'classint')
dev.off()



# Automatically detect the tissue type of the dataset

# load auto-detection function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

# guess a tissue type
tissue_guess = auto_detect_tissue_type(path_to_db_file = db_, seuratObject = sctrans2, scaled = TRUE, assay = "SCT")  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used
dim(tissue_guess) # 15  2
tissue_guess
#           tissue     score
# 14       Stomach 17103.009
# 6          Brain 14331.818
# 10     Intestine 14274.502
# 4            Eye  9972.007
# 7           Lung  7580.909
# 3          Liver  6696.392
# 1  Immune system  6129.240
# 2       Pancreas  5978.701
# 8        Adrenal  5957.614
# 5         Kidney  5057.339
# 9          Heart  4981.028
# 11        Muscle  4901.251
# 12      Placenta  4481.925
# 15        Thymus  4257.373
# 13        Spleen  2663.855

# image saved as "Tissue_labelling.PNG"



#
# Differential Expression Analyses
#


# Set sample identities to celltype/condition combinations
orig.ident = sctrans2$orig.ident

sctrans2$orig.ident = paste(sctrans2$classint, sctrans2$conditions)

Idents(sctrans2) <- paste(sctrans2$classint, sctrans2$conditions)
table(Idents(sctrans2))
#
#Oligodendrocyte precursor cell WT M                     Mural cell WT F 
#                                 56                                  76 
#                    Mural cell TG F                     Mural cell TG M 
#                                 79                                  86 
#Oligodendrocyte precursor cell WT F                 Ependymal cell TG F 
#                                 96                                  97 
#                    Macrophage WT M                     Mural cell WT M 
#                                100                                 105 
#                Ependymal cell WT M                 Ependymal cell WT F 
#                                106                                 114 
#                    Neuroblast WT M                     Macrophage WT F 
#                                129                                 133 
#                    Macrophage TG M                     Macrophage TG F 
#                                135                                 135 
#                        Neuron TG F                 Ependymal cell TG M 
#                                155                                 157 
#                    Neuroblast TG F                     Neuroblast WT F 
#                                183                                 184 
#Oligodendrocyte precursor cell TG M                         Neuron TG M 
#                                194                                 202 
#                        Neuron WT M                     Neuroblast TG M 
#                                209                                 250 
#                        Neuron WT F Oligodendrocyte precursor cell TG F 
#                                252                                 616 
#               Oligodendrocyte TG F                Oligodendrocyte WT M 
#                               1171                                1319 
#               Oligodendrocyte WT F                Oligodendrocyte TG M 
#                               1758                                1911 
#                     Astrocyte WT M               Endothelial cell WT M 
#                               1982                                2025 
#                     Astrocyte WT F               Endothelial cell TG M 
#                               2088                                2108 
#                     Astrocyte TG F                      Astrocyte TG M 
#                               2210                                2310 
#              Endothelial cell TG F               Endothelial cell WT F 
#                               2372                                2694 
#               Microglial cell WT M                Microglial cell TG F 
#                               2991                                4406 
#               Microglial cell TG M                Microglial cell WT F 
#                               4621                                4767 



#
# Global differential analysis
#

table(sctrans2@meta.data$conditions)
#  TG F  TG M  WT F  WT M
# 11424 11974 12162  9022
# sctrans2$cond <- conditions_cells
sctrans2$cond <- sctrans2@meta.data$conditions

Idents(sctrans2) <- "cond"

# don't use the following parameters: min.pct = -Inf, logfc.threshold = -Inf, min.cells.feature = 0, min.cells.group = 0, 
DEG_Global_Poisson_F <- FindMarkers(sctrans2, ident.1 = "TG F", ident.2 = "WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos = FALSE, logfc.threshold = -Inf)

dim(DEG_Global_Poisson_F)
# 896   5
write.csv(DEG_Global_Poisson_F, file="~/data/DEG_Global_Poisson_Female.csv", quote=F)


DEG_Global_Poisson_M <- FindMarkers(sctrans2, ident.1 = "TG M", ident.2 = "WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos = FALSE, logfc.threshold = -Inf)

dim(DEG_Global_Poisson_M)
# 851   5
write.csv(DEG_Global_Poisson_M, file="~/data/DEG_Global_Poisson_Male.csv", quote=F)




# Determining gender-specific and gender-dimorphic DEGs using adjusted significance and a nominal p-value specificity filter
# a minimum absolute logFC threshold for gender-dimorphic genes + min. abs. logFC for the target gender for gender-specific genes
# (still consider as gender-shared if significant in both genders with shared logFC, and abs logFC only above 0.25 in one gender)
gender_spec_genes = function(oligodendrocyte_M, oligodendrocyte_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25){

	# filter by target gene FDR
	male_canddeg = rownames(oligodendrocyte_M)[which(oligodendrocyte_M$p_val_adj < target_fdr)]
		
	# filter other gender non-significance (and not close to significance) nominal p-value threshold
	male_spec_genes = male_canddeg[which(oligodendrocyte_F[match(male_canddeg, rownames(oligodendrocyte_F)),]$p_val > exclude_pval)]
	
	# filter by abs. logFc threshold in target gender
	male_spec_genes = male_spec_genes[which(abs(oligodendrocyte_M[match(male_spec_genes, rownames(oligodendrocyte_M)),]$avg_log2FC) >minabslog)]
	
	
	print("Male-specific genes:")
	print(male_spec_genes)
	# print(oligodendrocyte_M[match(male_spec_genes, rownames(oligodendrocyte_M)),])
	#	               p_val  avg_log2FC pct.1 pct.2    p_val_adj
	#Tsc22d3 1.394385e-11 -0.07631600 0.102 0.146 2.604990e-07
	#Glul    2.038288e-11 -0.11120746 0.413 0.476 3.807930e-07
	#Klf6    1.102510e-09 -0.10227687 0.420 0.480 2.059709e-05
	#Cyth4   9.061919e-08  0.07703651 0.269 0.223 1.692948e-03
	# low logFCs!
	
	female_canddeg = rownames(oligodendrocyte_F)[which(oligodendrocyte_F$p_val_adj < target_fdr)]
		
	female_spec_genes = female_canddeg[which(oligodendrocyte_M[match(female_canddeg, rownames(oligodendrocyte_M)),]$p_val > exclude_pval)]
	
	# filter by abs. logFc threshold in target gender
	female_spec_genes = female_spec_genes[which(abs(oligodendrocyte_F[match(female_spec_genes, rownames(oligodendrocyte_F)),]$avg_log2FC) >minabslog)]	
	
	print("Female-specific genes:")
	print(female_spec_genes)
	
	# shared DEGs == gender-dimorphic
	
	dimorphic_genes = NULL
	intdegs = intersect(male_canddeg, female_canddeg)
	
	if(length(intdegs)){
		# different sign of the fold-change
		logfcs_male = oligodendrocyte_M[match(intdegs, rownames(oligodendrocyte_M)),]$avg_log2FC
		logfcs_female = oligodendrocyte_F[match(intdegs, rownames(oligodendrocyte_F)),]$avg_log2FC	
		
		diff_fc = intersect(which(sign(logfcs_male) != sign(logfcs_female)), intersect(which(abs(logfcs_male)>minabslog), which(abs(logfcs_female)>minabslog)))
		dimorphic_genes = intdegs[diff_fc]	
		
		shared_fc = which(sign(logfcs_male) == sign(logfcs_female))
		# optionally add: minabslog fulfilled in at least one of the genders
		shared_genes = intdegs[shared_fc]		
	}
	
	print("Gender-dimorphic genes:")
	print(dimorphic_genes)
	
	print("Gender-shared genes:")
	print(shared_genes)
			
	dfres = data.frame("DEG type"=c(rep("male-specific",length(male_spec_genes)), rep("female-specific",length(female_spec_genes)), rep("gender-dimorphic",length(dimorphic_genes)), rep("gender-shared",length(shared_genes))), "Gene symbols"=c(male_spec_genes, female_spec_genes, dimorphic_genes, shared_genes), "Male avg. logFC"=c(oligodendrocyte_M[match(male_spec_genes, rownames(oligodendrocyte_M)),]$avg_log2FC, oligodendrocyte_M[match(female_spec_genes, rownames(oligodendrocyte_M)),]$avg_log2FC, oligodendrocyte_M[match(dimorphic_genes, rownames(oligodendrocyte_M)),]$avg_log2FC, oligodendrocyte_M[match(shared_genes, rownames(oligodendrocyte_M)),]$avg_log2FC), "Female avg. logFC"=c(oligodendrocyte_F[match(male_spec_genes, rownames(oligodendrocyte_F)),]$avg_log2FC, oligodendrocyte_F[match(female_spec_genes, rownames(oligodendrocyte_F)),]$avg_log2FC, oligodendrocyte_F[match(dimorphic_genes, rownames(oligodendrocyte_F)),]$avg_log2FC, oligodendrocyte_F[match(shared_genes, rownames(oligodendrocyte_F)),]$avg_log2FC), "Male FDR"=c(oligodendrocyte_M[match(male_spec_genes, rownames(oligodendrocyte_M)),]$p_val_adj, oligodendrocyte_M[match(female_spec_genes, rownames(oligodendrocyte_M)),]$p_val_adj, oligodendrocyte_M[match(dimorphic_genes, rownames(oligodendrocyte_M)),]$p_val_adj, oligodendrocyte_M[match(shared_genes, rownames(oligodendrocyte_M)),]$p_val_adj), "Female FDR"=c(oligodendrocyte_F[match(male_spec_genes, rownames(oligodendrocyte_F)),]$p_val_adj, oligodendrocyte_F[match(female_spec_genes, rownames(oligodendrocyte_F)),]$p_val_adj, oligodendrocyte_F[match(dimorphic_genes, rownames(oligodendrocyte_F)),]$p_val_adj, oligodendrocyte_F[match(shared_genes, rownames(oligodendrocyte_F)),]$p_val_adj))
	
	return(dfres)
	
}



# Variant with min abs. logFC threshold
#global_gender_spec_genes = gender_spec_genes(DEG_Global_Poisson_M, DEG_Global_Poisson_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#character(0)
#[1] "Gender-dimorphic genes:"
#character(0)
#[1] "Gender-shared genes:"
# [1] "Ttr"      "mt-Rnr2"  "mt-Rnr1"  "Gm42418"  "Ccl4"     "Flt1"    
# [7] "Enpp2"    "Hexb"     "Gm37376"  "Lgmn"     "Abcb1a"   "Ccl3"    
#[13] "Actb"     "Sptbn1"   "Fosb"     "Tmsb4x"   "Sdpr"     "Ptprb"   
#[19] "Gm10800"  "Csf1r"    "C1qb"     "Igf1r"    "C1qa"     "C1qc"    
#[25] "Siglech"  "Ctsb"     "Cd9"      "Utrn"     "Fth1"     "Ivns1abp"
#[31] "Txnip"    "Tyrobp"   "Golim4"   "Dclk1"    "Ctsl"     "Olfml3"  
#[37] "Tmem119"  "mt-Co1"   "Ctsd"     "Neat1"    "Calr"     "Ly86"    
#[43] "Ctsz"     "Npc2"

# only identifies gender-shared genes (min abs. logFC criterion is too restrictive)



# Variant without min abs. logFC threshold
global_gender_spec_genes = gender_spec_genes(DEG_Global_Poisson_M, DEG_Global_Poisson_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0)
#[1] "Male-specific genes:"
# [1] "Ctss"   "Klf4"   "Car2"   "Cyth4"  "Egr1"   "Qdpr"   "Ccr5"   "Laptm5"
# [9] "Epas1"  "Gm4617" "Pnisr"  "Vim"    "Ptp4a2" "Itgb5" 
#[1] "Female-specific genes:"
# [1] "H3f3b"      "Ubb"        "mt-Nd4"     "Mfge8"      "Gm14586"   
# [6] "Gpr37l1"    "Gm5963"     "Rps15a"     "Rps10"      "Rps19"     
#[11] "Aldoa"      "Rps24-ps3"  "Gm5805"     "Rplp0"      "Rps10-ps1" 
#[16] "Myl6"       "Chchd2"     "Rps18"      "Fcho2"      "Rps3"      
#[21] "Hspa5"      "Slc25a4"    "Myl12b"     "Ppp1r12a"   "Cox6a1"    
#[26] "Rdx"        "Kmt2e"      "Slc1a3"     "Atp5h"      "Prdx1"     
#[31] "Rpl3-ps1"   "Rps5"       "Rps9"       "Smc6"       "Cnbp"      
#[36] "Rpl13"      "Nfe2l2"     "Kmt2a"      "Oaz1"       "Mid1ip1"   
#[41] "Atp5a1"     "Rpl36a-ps2" "Rabac1"     "Serpine2"   "Nfkbiz"    
#[46] "Strn3"      "Itgb1"      "Btg1"       "Cfdp1"      "Zbtb20"    
#[51] "Tmem30a"    "Ctnna1"    
#[1] "Gender-dimorphic genes:"						# exclude the gender-dimorphic genes: no major logFC difference, insufficient evidence
# [1] "Cx3cr1"   "Aldoc"    "Cst3"     "Ly6c1"    "Nfkbia"   "Plp1"    
# [7] "Sparcl1"  "Mt1"      "Hsp90b1"  "Klf2"     "Maf"      "Ly6a"    
#[13] "Fos"      "Ermn"     "Hsp90aa1" "Epb41l2"  "mt-Nd2"   "Mt2"     
#[19] "Htra1"    "Mef2c"    "Fyb"      "Stmn4"    "Gm12346"  "Bsg"     
#[25] "Itm2a"    "Rrbp1"    "Zeb2"     "Bin1"     "Apoe"     "Id3"     
#[31] "Tgfbr1"   "Rbm25"    "Ttyh1"    "Atp1a2"   "Mag"      "Rtn4"    
#[37] "Prdx6"    "Zfhx3"    "Ccdc88a"  "Plxdc2"   "Ier2"     "Elmo1"   
#[43] "mt-Nd1"   "Qk"       "Tsc22d1"  "Tanc2"    "Prrc2c"   "Sgk1"    
#[49] "Phf14"    "Gstm5"    "Cpe"      "Gja1"     "Cldn5"    "Mbnl1"   
#[55] "Atrx"     "Gstm1"    "Son"      "Dusp1"    "Itgam"    "Id2"     
#[61] "Zfp638"   "Pmepa1"   "Ckb"      "Akap9"    "Chd9"    
#[1] "Gender-shared genes:"
# [1] "Ttr"      "mt-Rnr2"  "mt-Rnr1"  "Gm42418"  "Ccl4"     "Flt1"    
# [7] "Enpp2"    "Hexb"     "Gm37376"  "Lgmn"     "Abcb1a"   "Ccl3"    
#[13] "Actb"     "Sptbn1"   "Fosb"     "Tmsb4x"   "Sdpr"     "Ptprb"   
#[19] "Gm10800"  "Csf1r"    "C1qb"     "Igf1r"    "C1qa"     "C1qc"    
#[25] "Siglech"  "Ctsb"     "Cd9"      "Utrn"     "Fth1"     "Ivns1abp"
#[31] "Txnip"    "Tyrobp"   "Golim4"   "Dclk1"    "Ctsl"     "Olfml3"  
#[37] "Tmem119"  "mt-Co1"   "Ctsd"     "Neat1"    "Calr"     "Ly86"    
#[43] "Ctsz"     "Npc2" 


save(DEG_Global_Poisson_M, DEG_Global_Poisson_F, global_gender_spec_genes, file="~/data/Global_DEGs.RData")
DEG_Global_Poisson <- rbind(DEG_Global_Poisson_M, DEG_Global_Poisson_F)
DEG_Global_Poisson$SYMBOL <- rownames(DEG_Global_Poisson)

table(row.names(DEG_Global_Poisson_F) %in% row.names(DEG_Global_Poisson_M))
# FALSE  TRUE
#   156   740

table(global_gender_spec_genes$DEG.type)
#  female-specific gender-dimorphic    gender-shared    male-specific
#               52               65               44               14



#
# Check overlap of DEGS between scRNAseq and Bulk RNAseq data (thy-tau22 AD models)
#

cd /home/m.ali/Projects/UL/mice_snRNAseq_cortex/enrico_script/


# load single-cell data (fdr significant)
load("~/data/Global_DEGs.RData")

dim(global_gender_spec_genes) # 175   6
table(global_gender_spec_genes$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#              52               65               44               14
sc_male_DEGs <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "male-specific",]
dim(sc_male_DEGs) # 14  6
sc_female_DEGs <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "female-specific",]
dim(sc_female_DEGs) # 52  6
sc_dimorphic_DEGs <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "gender-dimorphic",]
dim(sc_dimorphic_DEGs) # 65  6
sc_neutral_DEGs <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "gender-shared",]
dim(sc_neutral_DEGs) # 44  6


# load bulk data (fdr significant DEGs)
load("~/data/bulk_thy_tau22_DEGs_fdr.RData")

names(res_gender)
# [1] "male_spec"      "female_spec"    "gender_dimorph" "gender_shared"
sapply(res_gender, length)
#     male_spec    female_spec gender_dimorph  gender_shared
#           186              0              0              0

table(sc_male_DEGs$Gene.symbols %in% res_gender$male_spec)
# FALSE  TRUE
#    13     1

sc_male_DEGs[sc_male_DEGs$Gene.symbols %in% res_gender$male_spec,]
#         DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR Female.FDR
# 11 male-specific        Pnisr      0.05559323       0.004420577 0.0001280967          1


# NOTE: In case of FDR, only 1 male-specific gene is common between single-cell and bulk DEGs


# load bulk data (nominally significant)
load("~/data/bulk_thy_tau22_DEGs_nominal.RData")

names(res_gender)
# [1] "male_spec"      "female_spec"    "gender_dimorph" "gender_shared"
sapply(res_gender, length)
#     male_spec    female_spec gender_dimorph  gender_shared
#          1678            493             20            120


table(sc_male_DEGs$Gene.symbols %in% res_gender$male_spec)
# FALSE  TRUE
#    11     3
sc_male_DEGs[sc_male_DEGs$Gene.symbols %in% res_gender$male_spec,] # single-cell stats
#         DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR Female.FDR
# 5  male-specific         Egr1     -0.08308093     -0.0027920161 1.490953e-12          1
# 11 male-specific        Pnisr      0.05559323      0.0044205771 1.280967e-04          1
# 13 male-specific       Ptp4a2      0.03323807     -0.0003169573 3.191378e-02          1

tau_vs_wt_male <- as.data.frame(tau_vs_wt_male)
tau_vs_wt_male[sc_male_DEGs[sc_male_DEGs$Gene.symbols %in% res_gender$male_spec,]$Gene.symbols,] # bulk stats
#         baseMean log2FoldChange      lfcSE      stat       pvalue        padj
# Egr1    5308.739    -0.29427156 0.14401035 -2.043406 4.101230e-02 0.299640598
# Pnisr  13954.236     0.35735916 0.07997208  4.468549 7.875189e-06 0.003200543
# Ptp4a2  7170.063    -0.06875479 0.03084658 -2.228927 2.581873e-02 0.245548669



table(sc_female_DEGs$Gene.symbols %in% res_gender$female_spec)
# FALSE  TRUE
#    51     1
sc_female_DEGs[sc_female_DEGs$Gene.symbols %in% res_gender$female_spec,] # single-cell stats
#           DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC Male.FDR Female.FDR
# 61 female-specific        Itgb1    -0.001047833       -0.03621102        1  0.0120905

tau_vs_wt_female <- as.data.frame(tau_vs_wt_female)
tau_vs_wt_female[sc_female_DEGs[sc_female_DEGs$Gene.symbols %in% res_gender$female_spec,]$Gene.symbols,] # bulk stats
#        baseMean log2FoldChange      lfcSE     stat   pvalue      padj
# Itgb1  3306.321      0.2004609 0.07089295 2.827657 0.004689 0.9997837

table(sc_dimorphic_DEGs$Gene.symbols %in% res_gender$gender_dimorph)
# FALSE
#    65

table(sc_neutral_DEGs$Gene.symbols %in% res_gender$gender_shared)
# FALSE  TRUE
#    41     3
sc_neutral_DEGs[sc_neutral_DEGs$Gene.symbols %in% res_gender$gender_shared,] # single-cell stats
#          DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR    Female.FDR
# 154 gender-shared         C1qa      0.08158415        0.20923356 3.387348e-12 2.821833e-105
# 155 gender-shared         C1qc      0.07709631        0.18593094 1.675023e-10  1.359628e-82
# 173 gender-shared         Ly86      0.03789761        0.05930501 7.732349e-03  9.297272e-12

tau_vs_wt_male[sc_neutral_DEGs[sc_neutral_DEGs$Gene.symbols %in% res_gender$gender_shared,]$Gene.symbols,] # bulk stats
#       baseMean log2FoldChange      lfcSE      stat       pvalue         padj
# C1qa 1961.5884     -0.4591156 0.09667953 -4.748839 2.045874e-06 0.0012756662
# C1qc 1806.6322     -0.4019390 0.08153790 -4.929475 8.245102e-07 0.0007152805
# Ly86  508.2798     -0.5144373 0.10610652 -4.848310 1.245175e-06 0.0008873203
tau_vs_wt_female[sc_neutral_DEGs[sc_neutral_DEGs$Gene.symbols %in% res_gender$gender_shared,]$Gene.symbols,] # bulk stats
#       baseMean log2FoldChange      lfcSE      stat      pvalue      padj
# C1qa 2041.4803     -0.2062161 0.07680513 -2.684926 0.007254592 0.9997837
# C1qc 1867.0993     -0.2036050 0.07774670 -2.618824 0.008823337 0.9997837
# Ly86  499.8666     -0.2186813 0.10164669 -2.151387 0.031445690 0.9997837



library(venn) # venn_1.11 {Venn was too complex to be understood}
library(dplyr) # dplyr_1.1.2
overlap = list(M_SC = unique(sc_male_DEGs$Gene.symbols), F_SC = unique(sc_female_DEGs$Gene.symbols), GN_SC = unique(sc_neutral_DEGs$Gene.symbols),
	M_Bulk = unique(res_gender$male_spec), F_Bulk = unique(res_gender$female_spec), GN_Bulk = unique(res_gender$gender_shared))
flag_palette <- c("M_SC" = "black",
                  "F_SC" = "#FF883E",
                  "GN_SC" = "#ce2b37",
                  "M_Bulk" = "#009246",
                  "F_Bulk" = "#FCD116",
                  "GN_Bulk" = "blue")
# png("Venn_Overlapping_BioProc_thy21.png", width=8, height=6, units="in", res=300)
overlap %>% magrittr::extract(!(names(.) %in% c("bars", "stripes", "no sections"))) %>%
  # the palette has to be in the same order as the list of sets, so we reorder it
 venn::venn(zcolor = flag_palette[names(.)], cexil = 5, cexsn = 5)
# dev.off()


library(UpSetR) # UpSetR_1.4.0
png("~/figures/Venn_Overlapping_DEGs_SC_Bulk.png", width=12, height=10, units="in", res=300)
upset(fromList(overlap), sets=c("M_SC", "F_SC", "GN_SC", "M_Bulk", "F_Bulk", "GN_Bulk"), order.by = "freq", point.size = 5, line.size = 2, text.scale = 3)
dev.off()




#
# Cell-type specific DEG analysis
#

save(sctrans2, file="sctrans2.RData")


Idents(sctrans2) <- paste(sctrans2$classint, sctrans2$conditions)

#microglial_cell_F_logfcfilt <- FindMarkers(sctrans2, ident.1 = "Microglial cell TG F", ident.2 = "Microglial cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE)
#dim(microglial_cell_F_logfcfilt) # 20  5

microglial_cell_F <- FindMarkers(sctrans2, ident.1 = "Microglial cell TG F", ident.2 = "Microglial cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_F) # 931   5
nrow(microglial_cell_F[microglial_cell_F$p_val_adj < 0.05,]) # 393
nrow(microglial_cell_F[microglial_cell_F$p_val < 0.05,]) # 709


microglial_cell_F
#                p_val avg_log2FC pct.1 pct.2     p_val_adj
#Cst3     0.000000e+00  0.2927972 0.999 0.998  0.000000e+00
#Gm42418  0.000000e+00 -0.7996983 0.396 0.684  0.000000e+00
#Malat1   0.000000e+00 -0.4822471 0.989 0.998  0.000000e+00
#mt-Rnr2  0.000000e+00 -0.8787733 0.992 0.999  0.000000e+00
#Ccl4    6.589522e-288  0.5437800 0.464 0.312 1.231054e-283
#mt-Rnr1 4.357102e-194 -0.4619767 0.378 0.578 8.139937e-190
#Ccl3    2.410441e-159  0.4117526 0.388 0.257 4.503186e-155
#H3f3b   8.234762e-124  0.3599896 0.499 0.365 1.538418e-119
#C1qb    3.780894e-118  0.3127882 0.892 0.825 7.063465e-114
#B2m     2.707839e-117  0.3331968 0.799 0.688 5.058784e-113
#C1qa    1.585712e-106  0.3107610 0.852 0.764 2.962428e-102
#Xist    1.891499e-100 -0.2869008 0.143 0.317  3.533699e-96
#Fth1    8.264168e-100  0.3184041 0.677 0.537  1.543912e-95
#Zfp36    6.437948e-90  0.3015275 0.595 0.482  1.202737e-85
#C1qc     1.372144e-84  0.2733614 0.869 0.771  2.563440e-80
#mt-Cytb  1.581851e-79  0.2863151 0.611 0.474  2.955214e-75
#mt-Nd1   3.386725e-74  0.2761482 0.585 0.437  6.327080e-70
#Itm2b    3.568941e-74  0.2677213 0.784 0.676  6.667496e-70
#Rplp1    4.196868e-70  0.2675211 0.666 0.538  7.840588e-66
#Klf2     2.121588e-65  0.2514903 0.320 0.236  3.963550e-61


microglial_cell_M <- FindMarkers(sctrans2, ident.1 = "Microglial cell TG M", ident.2 = "Microglial cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(microglial_cell_M) # 911   5
nrow(microglial_cell_M[microglial_cell_M$p_val_adj < 0.05,]) # 66
nrow(microglial_cell_M[microglial_cell_M$p_val < 0.05,]) # 340

head(microglial_cell_M, 20)
#                p_val  avg_log2FC pct.1 pct.2     p_val_adj
#mt-Rnr2  0.000000e+00 -0.26883433 0.999 0.999  0.000000e+00
#Cst3    1.505074e-218 -0.23635177 0.996 0.997 2.811780e-214
#Cx3cr1   6.577550e-69  0.24925408 0.929 0.910  1.228818e-64
#Nfkbia   1.287527e-66 -0.28370371 0.447 0.560  2.405359e-62
#Apoe     3.539242e-65  0.27041880 0.163 0.133  6.612012e-61
#Ttr      9.369726e-62  0.20619663 0.137 0.029  1.750452e-57
#Egr1     2.860531e-61 -0.26349723 0.463 0.572  5.344043e-57
#Fos      7.222697e-60 -0.24427690 0.538 0.679  1.349344e-55
#Fosb     8.308985e-42 -0.22046489 0.268 0.392  1.552285e-37
#Ccl4     7.306415e-39  0.22197590 0.433 0.378  1.364984e-34
#Gm10800  4.221309e-36  0.19224966 0.228 0.124  7.886249e-32
#Plp1     8.159537e-31  0.16694754 0.210 0.120  1.524365e-26
#Klf4     1.158018e-30 -0.16799181 0.164 0.246  2.163409e-26
#Jun      3.634048e-28 -0.14164907 0.781 0.856  6.789128e-24
#mt-Rnr1  4.622726e-25 -0.15945812 0.692 0.659  8.636176e-21
#Junb     1.450068e-19 -0.15107532 0.513 0.591  2.709017e-15
#Hsp90b1  7.598988e-19  0.14653822 0.756 0.702  1.419643e-14
#Rbm25    4.934186e-18  0.14814166 0.501 0.440  9.218046e-14
#Malat1   7.211847e-17  0.05881859 0.998 0.999  1.347317e-12
#Hspa5    7.658203e-17 -0.13800230 0.614 0.649  1.430706e-12


microglia_gender_spec_genes = gender_spec_genes(microglial_cell_M, microglial_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#[1] "C1qb"  "C1qa"  "C1qc"  "Itm2b"
#[1] "Gender-dimorphic genes:"
#character(0)
#[1] "Gender-shared genes:"
# [1] "mt-Rnr2" "Apoe"    "Ttr"     "Ccl4"    "Gm10800" "mt-Rnr1" "Hspa5"  
# [8] "Ccl3"    "Lgmn"    "Actb"

dim(microglia_gender_spec_genes) # 82  6
microglia_gender_spec_genes
#          DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR
#1  female-specific         C1qb   -0.0080042637         0.3127882 1.000000e+00
#2  female-specific         C1qa   -0.0010253728         0.3107610 1.000000e+00
#3  female-specific         C1qc   -0.0010513935         0.2733614 1.000000e+00
#4  female-specific        Itm2b   -0.0006332102         0.2677213 1.000000e+00
#5    gender-shared      mt-Rnr2   -0.2688343285        -0.8787733 0.000000e+00
#6    gender-shared         Apoe    0.2704187958         0.2083971 6.612012e-61
#7    gender-shared          Ttr    0.2061966289         0.2398774 1.750452e-57
#8    gender-shared         Ccl4    0.2219759018         0.5437800 1.364984e-34
#9    gender-shared      Gm10800    0.1922496593         0.1146926 7.886249e-32
#10   gender-shared      mt-Rnr1   -0.1594581241        -0.4619767 8.636176e-21
#11   gender-shared        Hspa5   -0.1380022975        -0.1749584 1.430706e-12
#12   gender-shared         Ccl3    0.0971932114         0.4117526 2.018675e-04
#13   gender-shared         Lgmn    0.0900992440         0.1957592 1.819856e-03
#14   gender-shared         Actb    0.0697981244         0.2445548 4.691103e-03
#      Female.FDR
#1  7.063465e-114
#2  2.962428e-102
#3   2.563440e-80
#4   6.667496e-70
#5   0.000000e+00
#6   7.123830e-48
#7   7.007098e-76
#8  1.231054e-283
#9   7.524420e-23
#10 8.139937e-190
#11  1.572642e-26
#12 4.503186e-155
#13  3.040685e-35
#14  3.947490e-95




oligodendrocyte_cell_F <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte TG F", ident.2 = "Oligodendrocyte WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_cell_F) # 1087    5
nrow(oligodendrocyte_cell_F[oligodendrocyte_cell_F$p_val_adj < 0.05,]) # 120
nrow(oligodendrocyte_cell_F[oligodendrocyte_cell_F$p_val < 0.05,]) # 577

oligodendrocyte_cell_M <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte TG M", ident.2 = "Oligodendrocyte WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_cell_M) # 960   5
nrow(oligodendrocyte_cell_M[oligodendrocyte_cell_M$p_val_adj < 0.05,]) # 34
nrow(oligodendrocyte_cell_M[oligodendrocyte_cell_M$p_val < 0.05,]) # 367

oligodendrocyte_gender_spec_genes = gender_spec_genes(oligodendrocyte_cell_M, oligodendrocyte_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#[1] "Mbp"    "mt-Nd1"
#[1] "Gender-dimorphic genes:"
#character(0)
#[1] "Gender-shared genes:"
# [1] "Gm42418" "mt-Rnr2" "mt-Rnr1" "Malat1"  "Ttr"     "Car2"    "Neat1"  
# [8] "Gatm"    "mt-Co1"  "Qdpr"    "Dbi"     "Plp1"    "Sept4"   "mt-Cytb"

dim(oligodendrocyte_gender_spec_genes) # 27  6
oligodendrocyte_gender_spec_genes
#          DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC      Male.FDR
#1  female-specific          Mbp     0.007883511        0.26006125  1.000000e+00
#2  female-specific       mt-Nd1     0.009306023        0.30024805  1.000000e+00
#3    gender-shared      Gm42418    -0.464208688       -0.74770674 3.579321e-138
#4    gender-shared      mt-Rnr2    -0.189659564       -0.95690487 3.974834e-131
#5    gender-shared      mt-Rnr1    -0.281693693       -0.53965603  9.678615e-37
#6    gender-shared       Malat1    -0.121691119       -0.45843193  8.643463e-34
#7    gender-shared          Ttr     0.162410671        0.19532810  2.696372e-17
#8    gender-shared         Car2     0.200562954        0.22234183  1.770905e-12
#9    gender-shared        Neat1     0.184096635        0.19276578  1.601312e-08
#10   gender-shared         Gatm     0.181622484        0.15561603  3.820175e-08
#11   gender-shared       mt-Co1     0.162211793        0.44130114  8.054102e-07
#12   gender-shared         Qdpr     0.156055178        0.23928198  3.892832e-05
#13   gender-shared          Dbi     0.139072583        0.22981955  1.162777e-03
#14   gender-shared         Plp1    -0.033707189       -0.08947719  2.130333e-03
#15   gender-shared        Sept4     0.131924166        0.23019614  7.464089e-03
#16   gender-shared      mt-Cytb     0.131895522        0.45219439  7.955292e-03
#      Female.FDR
#1   9.473757e-27
#2   3.277292e-25
#3  1.188307e-142
#4   0.000000e+00
#5   2.548937e-75
#6  2.692605e-252
#7   4.614527e-14
#8   5.208515e-17
#9   5.732043e-09
#10  5.042899e-05
#11  6.052377e-56
#12  2.725318e-15
#13  2.749146e-13
#14  3.660846e-33
#15  2.094739e-14
#16  3.998055e-65



neuron_cell_F <- FindMarkers(sctrans2, ident.1 = "Neuron TG F", ident.2 = "Neuron WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_cell_F) # 1082    5
nrow(neuron_cell_F[neuron_cell_F$p_val_adj < 0.05,]) # 37
nrow(neuron_cell_F[neuron_cell_F$p_val < 0.05,]) # 433

neuron_cell_M <- FindMarkers(sctrans2, ident.1 = "Neuron TG M", ident.2 = "Neuron WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuron_cell_M) # 1109    5
nrow(neuron_cell_M[neuron_cell_M$p_val_adj < 0.05,]) # 4
nrow(neuron_cell_M[neuron_cell_M$p_val < 0.05,]) # 94

neuron_gender_spec_genes = gender_spec_genes(neuron_cell_M, neuron_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
[1] "Male-specific genes:"
character(0)
[1] "Female-specific genes:"
[1] "Ppp3ca" "Atp1b1" "Mef2c"  "Map1a"  "Jund"  
[1] "Gender-dimorphic genes:"
character(0)
[1] "Gender-shared genes:"
[1] "Gm42418" "Ttr" 

dim(neuron_gender_spec_genes) # 9 6
neuron_gender_spec_genes
#        DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR
#1 female-specific       Ppp3ca      0.01395822        -0.5246821 1.000000e+00
#2 female-specific       Atp1b1     -0.01411724        -0.5168074 1.000000e+00
#3 female-specific        Mef2c     -0.02223698        -0.4485916 1.000000e+00
#4 female-specific        Map1a      0.02343048        -0.3772508 1.000000e+00
#5 female-specific         Jund     -0.02177414         0.3155776 1.000000e+00
#6   gender-shared      Gm42418     -0.39299508        -0.9139937 3.481636e-06
#7   gender-shared          Ttr      0.28294483         0.4287832 7.607243e-03
#    Female.FDR
#1 8.883821e-07
#2 1.656255e-06
#3 8.290837e-05
#4 4.742604e-03
#5 1.669542e-02
#6 9.242819e-24
#7 2.117414e-07



astrocyte_cell_F <- FindMarkers(sctrans2, ident.1 = "Astrocyte TG F", ident.2 = "Astrocyte WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_cell_F) # 808   5
nrow(astrocyte_cell_F[astrocyte_cell_F$p_val_adj < 0.05,]) # 103
nrow(astrocyte_cell_F[astrocyte_cell_F$p_val < 0.05,]) # 434

astrocyte_cell_M <- FindMarkers(sctrans2, ident.1 = "Astrocyte TG M", ident.2 = "Astrocyte WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(astrocyte_cell_M) # 742   5
nrow(astrocyte_cell_M[astrocyte_cell_M$p_val_adj < 0.05,]) # 49
nrow(astrocyte_cell_M[astrocyte_cell_M$p_val < 0.05,]) # 255

astrocyte_gender_spec_genes = gender_spec_genes(astrocyte_cell_M, astrocyte_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#character(0)
#[1] "Gender-dimorphic genes:"
#[1] "Aldoc"
#[1] "Gender-shared genes:"
#[1] "mt-Rnr2" "Ttr"     "Clu"     "Mfge8"   "Plp1"    "Gm10800" "mt-Rnr1"
#[8] "Calr"    "Mt3" 

dim(astrocyte_gender_spec_genes) # 27  6
astrocyte_gender_spec_genes
#           DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR
#1  gender-dimorphic        Aldoc     -0.28450294        0.26726914 7.940903e-39
#2     gender-shared      mt-Rnr2     -0.08495282       -1.03666376 2.151239e-59
#3     gender-shared          Ttr      0.18550132        0.23780094 5.059471e-31
#4     gender-shared          Clu      0.19674391        0.12536956 1.547653e-18
#5     gender-shared        Mfge8      0.13848626        0.18348407 4.298140e-06
#6     gender-shared         Plp1      0.10575646        0.12815582 3.367077e-04
#7     gender-shared      Gm10800      0.10975635        0.17309358 4.509817e-04
#8     gender-shared      mt-Rnr1     -0.08151970       -0.77269101 2.069834e-03
#9     gender-shared         Calr      0.09702138        0.09541161 3.679050e-03
#10    gender-shared          Mt3      0.11259625        0.25487398 6.939229e-03
#      Female.FDR
#1   6.435591e-40
#2   0.000000e+00
#3   3.079699e-31
#4   9.530001e-08
#5   5.120273e-13
#6   6.057639e-07
#7   4.605669e-23
#8  2.188730e-302
#9   2.502061e-02
#10  1.286003e-26


microglia_gender_spec_genes$cell_type <- "microglialCells"
oligodendrocyte_gender_spec_genes$cell_type <- "oligodendrocytes"
neuron_gender_spec_genes$cell_type <- "neurons"
astrocyte_gender_spec_genes$cell_type <- "astrocytes"

astrocytes_res_tauCortex <- astrocyte_gender_spec_genes
microglial_res_tauCortex <- microglia_gender_spec_genes
neurons_res_tauCortex <- neuron_gender_spec_genes
oligodendrocyte_res_tauCortex <- oligodendrocyte_gender_spec_genes
save(oligodendrocyte_res_tauCortex, microglial_res_tauCortex, astrocytes_res_tauCortex, neurons_res_tauCortex, 
  file="~/data/TauCortex_DEGs_4CT_astro_oligo_microg_neuron.RData")



endothelial_cell_F <- FindMarkers(sctrans2, ident.1 = "Endothelial cell TG F", ident.2 = "Endothelial cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_F) # 1045
nrow(endothelial_cell_F[endothelial_cell_F$p_val_adj < 0.05,]) # 186 
nrow(endothelial_cell_F[endothelial_cell_F$p_val < 0.05,]) # 594

endothelial_cell_M <- FindMarkers(sctrans2, ident.1 = "Endothelial cell TG M", ident.2 = "Endothelial cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(endothelial_cell_M) # 942   5
nrow(endothelial_cell_M[endothelial_cell_M$p_val_adj < 0.05,]) # 34
nrow(endothelial_cell_M[endothelial_cell_M$p_val < 0.05,]) # 267

endothelial_gender_spec_genes = gender_spec_genes(endothelial_cell_M, endothelial_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#character(0)
#[1] "Gender-dimorphic genes:"
#character(0)
#[1] "Gender-shared genes:"
#[1] "mt-Rnr2" "Gm42418" "mt-Rnr1" "Ttr"     "Tmsb10"

dim(endothelial_gender_spec_genes) # 47  6
endothelial_gender_spec_genes
#       DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC      Male.FDR
#1 gender-shared      mt-Rnr2     -0.16724157        -0.8462723 3.101494e-169
#2 gender-shared      Gm42418     -0.35671754        -0.6518393 9.460860e-106
#3 gender-shared      mt-Rnr1     -0.36140497        -0.5284986  8.300946e-87
#4 gender-shared          Ttr      0.18057248         0.2826311  5.371449e-30
#5 gender-shared       Tmsb10      0.09906632         0.1375786  2.203363e-04
#     Female.FDR
#1  0.000000e+00
#2 1.871766e-238
#3 3.544488e-150
#4  2.282290e-62
#5  5.768844e-10




neuroblast_cell_F <- FindMarkers(sctrans2, ident.1 = "Neuroblast TG F", ident.2 = "Neuroblast WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuroblast_cell_F) # 1211    5
nrow(neuroblast_cell_F[neuroblast_cell_F$p_val_adj < 0.05,]) # 8
nrow(neuroblast_cell_F[neuroblast_cell_F$p_val < 0.05,]) # 298

neuroblast_cell_M <- FindMarkers(sctrans2, ident.1 = "Neuroblast TG M", ident.2 = "Neuroblast WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(neuroblast_cell_M) # 1130
nrow(neuroblast_cell_M[neuroblast_cell_M$p_val_adj < 0.05,]) # 2
nrow(neuroblast_cell_M[neuroblast_cell_M$p_val < 0.05,]) # 118

neuroblast_gender_spec_genes = gender_spec_genes(neuroblast_cell_M, neuroblast_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#character(0)
#[1] "Gender-dimorphic genes:"
#character(0)
#[1] "Gender-shared genes:"
#[1] "mt-Rnr2" "Gm42418"

dim(neuroblast_gender_spec_genes) # 2 6
neuroblast_gender_spec_genes
#       DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR
#1 gender-shared      mt-Rnr2      -0.2825041        -0.7434547 1.160800e-30
#2 gender-shared      Gm42418      -0.4198996        -0.6216388 9.372428e-09
#     Female.FDR
#1 3.255496e-130
#2  1.689604e-12



oligodendrocyte_precursor_cell_F <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte precursor cell TG F", ident.2 = "Oligodendrocyte precursor cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_F) # 959
nrow(oligodendrocyte_precursor_cell_F[oligodendrocyte_precursor_cell_F$p_val_adj < 0.05,]) # 100 
nrow(oligodendrocyte_precursor_cell_F[oligodendrocyte_precursor_cell_F$p_val < 0.05,]) # 566

oligodendrocyte_precursor_cell_M <- FindMarkers(sctrans2, ident.1 = "Oligodendrocyte precursor cell TG M", ident.2 = "Oligodendrocyte precursor cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(oligodendrocyte_precursor_cell_M) # 1141
nrow(oligodendrocyte_precursor_cell_M[oligodendrocyte_precursor_cell_M$p_val_adj < 0.05,]) # 5
nrow(oligodendrocyte_precursor_cell_M[oligodendrocyte_precursor_cell_M$p_val < 0.05,]) # 122

oligodendrocyte_precursor_gender_spec_genes = gender_spec_genes(oligodendrocyte_precursor_cell_M, oligodendrocyte_precursor_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
# [1] "Igfbp2" "Folr1"  "Pcp4l1" "Arl3"   "Gm6265" "Cd63"   "Cdkn1c" "Ctsl"  
# [9] "Kl"     "Romo1"  "Cd59a"  "Dstn"  
#[1] "Gender-dimorphic genes:"
#[1] "mt-Rnr2" "Ttr"     "mt-Rnr1" "Malat1"  "Slc1a2" 
#[1] "Gender-shared genes:"
#character(0)

dim(oligodendrocyte_precursor_gender_spec_genes) # 25  6
oligodendrocyte_precursor_gender_spec_genes
#           DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC
#1   female-specific       Igfbp2   -0.0004697547        -0.6004985
#2   female-specific        Folr1    0.0147970019        -0.3794068
#3   female-specific       Pcp4l1   -0.0193287820        -0.5050209
#4   female-specific         Arl3    0.0006110957        -0.3642414
#5   female-specific       Gm6265    0.0370706309        -0.4565560
#6   female-specific         Cd63    0.0098815115        -0.3714069
#7   female-specific       Cdkn1c    0.0623489314        -0.3519581
#8   female-specific         Ctsl    0.0334126801        -0.3247821
#9   female-specific           Kl    0.0393193211        -0.2821271
#10  female-specific        Romo1    0.0176175210        -0.2874263
#11  female-specific        Cd59a    0.0328685749        -0.2730911
#12  female-specific         Dstn    0.0545675306        -0.2928735
#13 gender-dimorphic      mt-Rnr2   -0.7599459674         1.2185605
#14 gender-dimorphic          Ttr    0.7303901493        -1.7811761
#15 gender-dimorphic      mt-Rnr1   -0.8890659010         0.9804805
#16 gender-dimorphic       Malat1   -0.4423816439         1.4192250
#17 gender-dimorphic       Slc1a2   -0.5079429837         1.0167335
#        Male.FDR    Female.FDR
#1   1.000000e+00  1.317980e-15
#2   1.000000e+00  3.211095e-10
#3   1.000000e+00  7.665716e-10
#4   1.000000e+00  3.200278e-09
#5   1.000000e+00  1.560243e-07
#6   1.000000e+00  2.145034e-05
#7   1.000000e+00  4.100514e-05
#8   1.000000e+00  1.098515e-04
#9   1.000000e+00  1.006407e-03
#10  1.000000e+00  1.328087e-03
#11  1.000000e+00  2.188528e-03
#12  1.000000e+00  2.620715e-02
#13 6.797236e-158  0.000000e+00
#14 5.210689e-145  0.000000e+00
#15  1.282432e-21  4.017749e-24
#16  1.362932e-07 1.017452e-134
#17  2.917570e-04  8.217947e-14



mural_cell_F <- FindMarkers(sctrans2, ident.1 = "Mural cell TG F", ident.2 = "Mural cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_F) # 1422
nrow(mural_cell_F[mural_cell_F$p_val_adj < 0.05,]) # 15
nrow(mural_cell_F[mural_cell_F$p_val < 0.05,]) # 250

mural_cell_M <- FindMarkers(sctrans2, ident.1 = "Mural cell TG M", ident.2 = "Mural cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(mural_cell_M) # 1308
nrow(mural_cell_M[mural_cell_M$p_val_adj < 0.05,]) # 8
nrow(mural_cell_M[mural_cell_M$p_val < 0.05,]) # 137

mural_gender_spec_genes = gender_spec_genes(mural_cell_M, mural_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#[1] "mt-Cytb" "Col1a2"  "Sparcl1"
#[1] "Gender-dimorphic genes:"
#[1] "Atp1a2" "Mgp"   
#[1] "Gender-shared genes:"
#[1] "mt-Rnr2" "Gm42418" "mt-Rnr1"

dim(mural_gender_spec_genes) # 9 6
mural_gender_spec_genes
#          DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR
#1  female-specific      mt-Cytb     0.000178452         0.6124707 1.000000e+00
#2  female-specific       Col1a2     0.000178452         0.6434873 1.000000e+00
#3  female-specific      Sparcl1     0.014529741        -0.5901897 1.000000e+00
#4 gender-dimorphic       Atp1a2     0.704031035        -0.7132017 6.260753e-11
#5 gender-dimorphic          Mgp    -0.530240893         0.9928410 1.268099e-05
#6    gender-shared      mt-Rnr2    -0.173850332        -1.0493517 6.038696e-09
#7    gender-shared      Gm42418    -0.452090742        -0.9739420 1.276914e-04
#8    gender-shared      mt-Rnr1    -0.347862907        -1.0147480 2.114328e-02
#     Female.FDR
#1  3.181512e-06
#2  3.553994e-03
#3  1.271891e-02
#4  3.933838e-05
#5  1.724677e-24
#6 2.738248e-214
#7  7.276531e-12
#8  1.377016e-14




macrophage_cell_F <- FindMarkers(sctrans2, ident.1 = "Macrophage TG F", ident.2 = "Macrophage WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_cell_F) # 1080    5
nrow(macrophage_cell_F[macrophage_cell_F$p_val_adj < 0.05,]) # 9
nrow(macrophage_cell_F[macrophage_cell_F$p_val < 0.05,]) # 208

macrophage_cell_M <- FindMarkers(sctrans2, ident.1 = "Macrophage TG M", ident.2 = "Macrophage WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(macrophage_cell_M) # 1041
nrow(macrophage_cell_M[macrophage_cell_M$p_val_adj < 0.05,]) # 4
nrow(macrophage_cell_M[macrophage_cell_M$p_val < 0.05,]) # 104

macrophage_gender_spec_genes = gender_spec_genes(macrophage_cell_M, macrophage_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#[1] "Ttr"     "Lyz2"    "mt-Rnr1" "Prrc2c" 
#[1] "Gender-dimorphic genes:"
#character(0)
#[1] "Gender-shared genes:"
#[1] "mt-Rnr2"

dim(macrophage_gender_spec_genes) # 4 6
macrophage_gender_spec_genes
#         DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR
#1 female-specific          Ttr     0.208586622         0.7038445 1.000000e+00
#2 female-specific         Lyz2    -0.006299988         0.4942493 1.000000e+00
#3 female-specific      mt-Rnr1    -0.031378188        -0.5069600 1.000000e+00
#4 female-specific       Prrc2c     0.052467420        -0.3881303 1.000000e+00
#5   gender-shared      mt-Rnr2    -0.171283846        -0.6279812 4.084329e-06
#    Female.FDR
#1 9.431672e-09
#2 1.284876e-06
#3 1.601417e-04
#4 4.312942e-02
#5 3.927798e-82



ependymal_cell_F <- FindMarkers(sctrans2, ident.1 = "Ependymal cell TG F", ident.2 = "Ependymal cell WT F", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_F) # 1360    5
nrow(ependymal_cell_F[ependymal_cell_F$p_val_adj < 0.05,]) # 8
nrow(ependymal_cell_F[ependymal_cell_F$p_val < 0.05,]) # 167

ependymal_cell_M <- FindMarkers(sctrans2, ident.1 = "Ependymal cell TG M", ident.2 = "Ependymal cell WT M", verbose = FALSE, test.use = "poisson", assay="SCT", only.pos=FALSE, logfc.threshold = -Inf)
dim(ependymal_cell_M) # 1055
nrow(ependymal_cell_M[ependymal_cell_M$p_val_adj < 0.05,]) # 2
nrow(ependymal_cell_M[ependymal_cell_M$p_val < 0.05,]) # 106

ependymal_gender_spec_genes = gender_spec_genes(ependymal_cell_M, ependymal_cell_F, target_fdr = 0.05, exclude_pval = 0.5, minabslog=0.25)
#[1] "Male-specific genes:"
#character(0)
#[1] "Female-specific genes:"
#[1] "Mt3"
#[1] "Gender-dimorphic genes:"
#[1] "Nnat"
#[1] "Gender-shared genes:"
#[1] "Gm42418"
# 

dim(ependymal_gender_spec_genes) # 3 6
ependymal_gender_spec_genes
#          DEG.type Gene.symbols Male.avg..logFC Female.avg..logFC     Male.FDR
#1  female-specific          Mt3    -0.008821673         0.6165037 1.000000e+00
#2 gender-dimorphic         Nnat    -0.472548857         0.4805448 1.440311e-05
#3    gender-shared      Gm42418    -0.332317672        -0.5477098 1.340760e-03
#    Female.FDR
#1 6.744238e-06
#2 1.971671e-04
#3 3.938983e-04


microglia_gender_spec_genes$cell_type <- "microglialCells"
oligodendrocyte_gender_spec_genes$cell_type <- "oligodendrocytes"
neuron_gender_spec_genes$cell_type <- "neuron"
astrocyte_gender_spec_genes$cell_type <- "astrocytes"
endothelial_gender_spec_genes$cell_type <- "endothelial"
neuroblast_gender_spec_genes$cell_type <- "neuroblast"
oligodendrocyte_precursor_gender_spec_genes$cell_type <- "OPCs"
mural_gender_spec_genes$cell_type <- "mural"
macrophage_gender_spec_genes$cell_type <- "macrophages"
ependymal_gender_spec_genes$cell_type <- "ependymals"

astrocytes_res_tauCortex <- astrocyte_gender_spec_genes
microglial_res_tauCortex <- microglia_gender_spec_genes
neurons_res_tauCortex <- neuron_gender_spec_genes
oligodendrocyte_res_tauCortex <- oligodendrocyte_gender_spec_genes
endothelial_res_tauCortex <- endothelial_gender_spec_genes
neuroblast_res_tauCortex <- neuroblast_gender_spec_genes
OPC_res_tauCortex <- oligodendrocyte_precursor_gender_spec_genes
mural_res_tauCortex <- mural_gender_spec_genes
macrophage_res_tauCortex <- macrophage_gender_spec_genes
ependymal_res_tauCortex <- ependymal_gender_spec_genes

save(oligodendrocyte_res_tauCortex, microglial_res_tauCortex, astrocytes_res_tauCortex, neurons_res_tauCortex, endothelial_res_tauCortex,
	neuroblast_res_tauCortex, OPC_res_tauCortex, mural_res_tauCortex, macrophage_res_tauCortex, ependymal_res_tauCortex,
  file="~/data/TauCortex_DEGs_All10CTs.RData")

CT_specific_DEGs <- rbind(oligodendrocyte_res_tauCortex, microglial_res_tauCortex, astrocytes_res_tauCortex, neurons_res_tauCortex, 
	endothelial_res_tauCortex, neuroblast_res_tauCortex, OPC_res_tauCortex, mural_res_tauCortex, macrophage_res_tauCortex, ependymal_res_tauCortex)
dim(CT_specific_DEGs) # 235   7
write.csv(CT_specific_DEGs, file="~/data/TauCortex_DEGs_All10CTs.csv", row.names=F, quote=F)


table(astrocytes_res_tauCortex$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#               16                1                9                1
table(microglial_res_tauCortex$DEG.type)
# female-specific   gender-shared   male-specific
#              68              10               4
table(neurons_res_tauCortex$DEG.type)
# female-specific   gender-shared   male-specific
#               6               2               1
table(oligodendrocyte_res_tauCortex$DEG.type)
# female-specific   gender-shared   male-specific
#              11              14               2
table(endothelial_res_tauCortex$DEG.type)
# female-specific   gender-shared   male-specific
#              41               5               1
table(neuroblast_res_tauCortex$DEG.type)
# gender-shared
#             2
table(OPC_res_tauCortex$DEG.type)
# female-specific gender-dimorphic
#              20                5
table(mural_res_tauCortex$DEG.type)
# female-specific gender-dimorphic    gender-shared
#               3                3                3
table(macrophage_res_tauCortex$DEG.type)
# female-specific   gender-shared
#               3               1
table(ependymal_res_tauCortex$DEG.type)
# female-specific gender-dimorphic    gender-shared
#               1                1                1

		
write.xlsx(list("Microglial cell"=microglia_gender_spec_genes, "Neuron"=neuron_gender_spec_genes, "Astrocyte"=astrocyte_gender_spec_genes, 
	"Oligodendrocyte"=oligodendrocyte_gender_spec_genes, "Oligodendrocyte precursor"=oligodendrocyte_precursor_gender_spec_genes, 
	"Endothelial"=endothelial_gender_spec_genes, "Neuroblast"=neuroblast_gender_spec_genes, "Mural cell"=mural_gender_spec_genes, 
	"Macrophage"=macrophage_gender_spec_genes, "Ependymal cell"=ependymal_gender_spec_genes),'~/data/thy_tau22_single_cell_cortex_differential_analysis_gender.xlsx')

# count number of gender-specific genes per cell-type (male and female)
gen_degs = list("Microglial cell"=microglia_gender_spec_genes, "Neuron"=neuron_gender_spec_genes, "Astrocyte"=astrocyte_gender_spec_genes, 
	"Oligodendrocyte"=oligodendrocyte_gender_spec_genes, "Oligodendrocyte precursor"=oligodendrocyte_precursor_gender_spec_genes, 
	"Endothelial"=endothelial_gender_spec_genes, "Neuroblast"=neuroblast_gender_spec_genes, "Mural cell"=mural_gender_spec_genes, 
	"Macrophage"=macrophage_gender_spec_genes, "Ependymal cell"=ependymal_gender_spec_genes)


sapply(gen_degs, function(x) nrow(x))
#          Microglial cell                    Neuron                 Astrocyte 
#                       14                         7                        10 
#          Oligodendrocyte Oligodendrocyte precursor               Endothelial 
#                       16                        17                         5 
#               Neuroblast                Mural cell                Macrophage 
#                        2                         8                         5 
#           Ependymal cell 
#                        3

# majority of changes always in oligodendrocytes or oligodendrocyte precurosor cells in single-cell analyses (bias due to best marker genes for oligodendrocytes?)
# then microglia next most common


sapply(gen_degs, function(x) length(which(x$DEG.type=="female-specific")))
#          Microglial cell                    Neuron                 Astrocyte 
#                        4                         5                         0 
#          Oligodendrocyte Oligodendrocyte precursor               Endothelial 
#                        2                        12                         0 
#               Neuroblast                Mural cell                Macrophage 
#                        0                         3                         4 
#           Ependymal cell 
#                        1 

# most female-changes also in oligodendrocyte precursor cells

                        
sapply(gen_degs, function(x) length(which(x$DEG.type=="male-specific")))
# nothing

sapply(gen_degs, function(x) length(which(x$DEG.type=="gender-dimorphic")))
          Microglial cell                    Neuron                 Astrocyte 
                        0                         0                         1 
          Oligodendrocyte Oligodendrocyte precursor               Endothelial 
                        0                         5                         0 
               Neuroblast                Mural cell                Macrophage 
                        0                         2                         0 
           Ependymal cell 
                        1

# most gender-dimorphic changes also in oligodendrocyte precursor cells