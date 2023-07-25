#
# Gene set enrichment analysis
#

library(clusterProfiler) # clusterProfiler_4.8.1
library(org.Mm.eg.db) # org.Mm.eg.db_3.17.0
library(enrichplot) # enrichplot_1.20.0 {http://yulab-smu.top/clusterProfiler-book/chapter12.html}
library(tidyverse) # tidyverse_2.0.0
library(cowplot) # cowplot_1.1.1
library(venn) # venn_1.11
set.seed(1)


load("~/data/Global_DEGs.RData")

DEG_Global_Poisson <- rbind(DEG_Global_Poisson_M, DEG_Global_Poisson_F)
DEG_Global_Poisson$SYMBOL <- rownames(DEG_Global_Poisson)

table(row.names(DEG_Global_Poisson_F) %in% row.names(DEG_Global_Poisson_M))
# FALSE  TRUE
#   156   740

table(global_gender_spec_genes$DEG.type)
#  female-specific gender-dimorphic    gender-shared    male-specific
#               52               65               44               14


#gene.df <- bitr(rownames(DEG_Global_Poisson), fromType = "SYMBOL",
gene.df <- bitr(global_gender_spec_genes$Gene.symbols[which(global_gender_spec_genes$DEG.type=="gender-shared")], fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# Warning message: 11.36% of input gene IDs fail to map
dim(gene.df) # 39 2


# GOI <- left_join(gene.df, rownames(DEG_Global_Poisson), by=c("SYMBOL"="rownames(DEG_Global_Poisson)"))
GOI <- left_join(gene.df, DEG_Global_Poisson, by="SYMBOL")
dim(GOI) # 39  7


BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 132 9
#plot(BioProc)

png("~/figures/gender_shared_DEGs_enrichment_BP.png", width=12, height=8, units="in", res=300)
barplot(BioProc, showCategory=20)
dev.off()
png("~/figures/gender_shared_DEGs_enrichment_BP_v2.png", width=12, height=10, units="in", res=300)
barplot(BioProc, showCategory=17, font.size=20)
dev.off()

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 22 9
#yy <- enrichGO(GOI$ENTREZID, 'org.Mm.eg.db', ont="MF", pvalueCutoff=0.05, pAdjustMethod = "none")
#head(summary(yy))
#plot(yy)
png("~/figures/gender_shared_DEGs_enrichment_MF.png", width=12, height=8, units="in", res=300)
barplot(MolFun, showCategory=20)
dev.off()


kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'mmu',
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.05)
dim(kk) # 12 9
png("~/figures/gender_shared_DEGs_enrichment_KEGG.png", width=12, height=8, units="in", res=300)
barplot(kk, showCategory=10)
dev.off()
write.table(kk, "~/data/gender_shared_DEGs_enrichment_KEGG_FDRadj.txt", sep="\t", row.names=F, quote=F)

require(cowplot) # cowplot_1.1.1
p1 <- dotplot(BioProc, showCategory=18) + ggtitle("dotplot for BP enrichment")
# dotplot(BioProc, showCategory=15, split=".sign") + facet_grid(.~.sign)
#p2 <- dotplot(kk, showCategory=20) + ggtitle("dotplot for KEGG enrichment")
p2 <- dotplot(MolFun, showCategory=18) + ggtitle("dotplot for MF enrichment")
png("~/figures/gender_shared_DEGs_enrichment_BP_MF.png", width=16, height=10, units="in", res=300)
plot_grid(p1, p2, ncol=2)
dev.off()

png("~/figures/gender_shared_DEGs_enrichment_BP_MF_v2.png", width=16, height=8, units="in", res=300)
plot_grid(p1, p2, ncol=2)
dev.off()
save(BioProc, MolFun, file="~/data/gender_shared_DEGs_enrichment_BP_MF_v2.RData")


BP_GenderShared <- BioProc
MF_GenderShared <- MolFun
KEGG_GenderShared <- kk


# gseGO Fig. 5 if you want in "Activate/Suppressed" fashion
# we want the log2 fold change
DEA_padj <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "gender-shared",]
original_gene_list <- DEA_padj$Female.avg..logFC
# name the vector
names(original_gene_list) <- DEA_padj$Gene.symbols
# omit any NA values 
gene_list <- na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
length(gene_list) # 44

gse_BP <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             eps = 0, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none") # with "none" = 50 pathways, with "fdr" or "BH", we get 0 pathways)
dim(gse_BP) # 37 11

gse_MF <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "SYMBOL", 
             eps = 0, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none") # with "none" = 50 pathways, with "fdr" or "BH", we get 0 pathways)
dim(gse_MF) # 9 11

gse_BP@result[gse_BP@result$Description == "adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",]$Description <- "GO:0002460"
gse_BP@result[gse_BP@result$Description == "humoral immune response mediated by circulating immunoglobulin",]$Description <- "humoral immune response"

p1 <- dotplot(gse_BP, showCategory=8, split=".sign") + facet_grid(.~.sign) # Fig. 5 if you want in "Activate/Suppressed" fashion
p2 <- dotplot(gse_MF, showCategory=8, split=".sign") + facet_grid(.~.sign) # Fig. 5 if you want in "Activate/Suppressed" fashion
# write.table(gse_BP, file="gender_shared_DEGs_enrichment_gseBP_nominal.txt", sep="\t", row.names=F, quote=F)

png("~/figures/gender_shared_DEGs_enrichment_BP_MF_Splitted.png", width=16, height=8, units="in", res=300)
plot_grid(p1, p2, ncol=2)
dev.off()





#
# female-specific
#

gene.df <- bitr(global_gender_spec_genes$Gene.symbols[which(global_gender_spec_genes$DEG.type=="female-specific")], fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)
# 1.92% of input gene IDs fail to map


BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 153 9

png("~/figures/DEGs_enrichment_BP_female.png", width=12, height=8, units="in", res=300)
barplot(BioProc, showCategory=20)
dev.off()

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 6 9
#yy <- enrichGO(GOI$ENTREZID, 'org.Mm.eg.db', ont="MF", pvalueCutoff=0.05, pAdjustMethod = "none")
#head(summary(yy))
#plot(yy)
png("~/figures/DEGs_enrichment_MF_female.png", width=12, height=8, units="in", res=300)
barplot(MolFun, showCategory=20)
dev.off()


kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'mmu',
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.05)
dim(kk) # 9 9
png("DEGs_enrichment_KEGG_female.png", width=12, height=8, units="in", res=300)
barplot(kk, showCategory=10)
dev.off()
write.table(kk, "~/data/DEGs_enrichment_KEGG_FDRadj_female.txt", sep="\t", row.names=F, quote=F)


p1 <- dotplot(BioProc, showCategory=20, font.size=9) + ggtitle("dotplot for BP enrichment")
#p2 <- dotplot(kk, showCategory=20) + ggtitle("dotplot for KEGG enrichment")
p2 <- dotplot(MolFun, showCategory=20) + ggtitle("dotplot for MF enrichment")
png("~/figures/DEGs_enrichment_BP_MF_female.png", width=16, height=8, units="in", res=300)
plot_grid(p1, p2, ncol=2)
dev.off()

BP_Female <- BioProc
MF_Female <- MolFun
KEGG_Female <- kk


# gseGO

# we want the log2 fold change
DEA_padj <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "female-specific",]
original_gene_list <- DEA_padj$Female.avg..logFC
# name the vector
names(original_gene_list) <- DEA_padj$Gene.symbols
# omit any NA values 
gene_list <- na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
length(gene_list) # 52

gse_BP_female <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             eps = 0, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")
dim(gse_BP_female) # 57 11

gse_MF_female <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "SYMBOL", 
             eps = 0, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none") # with "none" = 50 pathways, with "fdr" or "BH", we get 0 pathways)
dim(gse_MF_female) # 5 11

p1 <- dotplot(gse_BP_female, showCategory=7, split=".sign") + facet_grid(.~.sign) # Fig. 5 if you want in "Activate/Suppressed" fashion
p2 <- dotplot(gse_MF_female, showCategory=7, split=".sign") + facet_grid(.~.sign) # Fig. 5 if you want in "Activate/Suppressed" fashion

png("~/figures/DEGs_enrichment_BP_MF_female_Splitted.png", width=16, height=8, units="in", res=300)
plot_grid(p1, p2, ncol=2)
dev.off()



#
# male-specific
#

gene.df <- bitr(global_gender_spec_genes$Gene.symbols[which(global_gender_spec_genes$DEG.type=="male-specific")], fromType = "SYMBOL",
        toType = c("ENTREZID", "SYMBOL"),
        OrgDb = org.Mm.eg.db)


BioProc <- enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "fdr", # replace "none" with "fdr" and "BH"
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(BioProc) # 33 9

png("~/figures/DEGs_enrichment_BP_male.png", width=12, height=8, units="in", res=300)
barplot(BioProc, showCategory=20)
dev.off()

MolFun <- enrichGO(gene        = gene.df$ENTREZID,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "MF",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
dim(MolFun) # 20 9
#yy <- enrichGO(GOI$ENTREZID, 'org.Mm.eg.db', ont="MF", pvalueCutoff=0.05, pAdjustMethod = "none")
#head(summary(yy))
#plot(yy)
png("~/figures/DEGs_enrichment_MF_male.png", width=12, height=8, units="in", res=300)
barplot(MolFun, showCategory=20)
dev.off()


kk <- enrichKEGG(gene         = gene.df$ENTREZID,
                 organism     = 'mmu',
                 pAdjustMethod = "fdr",
                 pvalueCutoff = 0.05)
dim(kk) # 0 9
png("~/figures/DEGs_enrichment_KEGG_male.png", width=12, height=8, units="in", res=300)
barplot(kk, showCategory=10)
dev.off()


p1 <- dotplot(BioProc, showCategory=20, font.size=8) + ggtitle("dotplot for BP enrichment")
#p2 <- dotplot(kk, showCategory=20) + ggtitle("dotplot for KEGG enrichment")
p2 <- dotplot(MolFun, showCategory=20, font.size=7.5) + ggtitle("dotplot for MF enrichment")
png("~/figures/DEGs_enrichment_BP_MF_male.png", width=16, height=8, units="in", res=300)
plot_grid(p1, p2, ncol=2)
dev.off()

BP_Male <- BioProc
MF_Male <- MolFun
KEGG_Male <- kk


# gseGO

# we want the log2 fold change
DEA_padj <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "male-specific",]
original_gene_list <- DEA_padj$Male.avg..logFC
# name the vector
names(original_gene_list) <- DEA_padj$Gene.symbols
# omit any NA values 
gene_list <- na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
length(gene_list) # 14

gse_BP_male <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "SYMBOL", 
             eps = 0, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")
dim(gse_BP_male) # 75 11
# write.table(gse, file="DEGs_enrichment_gseGO_nominal_male.txt", sep="\t", row.names=F, quote=F)


gse_MF_male <- gseGO(geneList=gene_list, 
             ont ="MF", 
             keyType = "SYMBOL", 
             eps = 0, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none") # with "none" = 50 pathways, with "fdr" or "BH", we get 0 pathways)
dim(gse_MF_male) # 18 11

p1 <- dotplot(gse_BP_male, showCategory=7, split=".sign") + facet_grid(.~.sign) # Fig. 5 if you want in "Activate/Suppressed" fashion
p2 <- dotplot(gse_MF_male, showCategory=7, split=".sign") + facet_grid(.~.sign) # Fig. 5 if you want in "Activate/Suppressed" fashion

png("~/figures/DEGs_enrichment_BP_MF_male_Splitted.png", width=16, height=8, units="in", res=300)
plot_grid(p1, p2, ncol=2)
dev.off()


save(BP_Male, MF_Male, KEGG_Male, KEGG_Female, BP_Female, MF_Female, BP_GenderShared, MF_GenderShared, KEGG_GenderShared, file="~/data/Pathway_Enrichment_Analysis.RData")
save(gse_BP, gse_MF, gse_BP_female, gse_MF_female, gse_BP_male, gse_MF_male, file="~/data/Pathway_Enrichment_Analysis_gseGO_Splitted.RData")



#####
## Get the number of genes in a BP/MF pathway that belong to up/down regulated DEGs
#####

library(ggplot2) # ggplot2_3.4.2
library(cowplot) # cowplot_1.1.1
library(clusterProfiler) # clusterProfiler_4.8.1
rm(list=ls())
set.seed(1)


load("Pathway_Enrichment_Analysis.RData")
load("Global_DEGs.RData")
DEG_Global_Poisson <- rbind(DEG_Global_Poisson_M, DEG_Global_Poisson_F)
DEG_Global_Poisson$SYMBOL <- rownames(DEG_Global_Poisson)


# head(MF_GenderShared, 2)

DEG_GenderShared <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "gender-shared",]
dim(DEG_GenderShared) # 44  6
DEG_Male <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "male-specific",]
dim(DEG_Male) # 14  6
DEG_Female <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "female-specific",]
dim(DEG_Female) # 52  6


# Gender-Shared
MF_GenderShared@result$up <- NA
MF_GenderShared@result$down <- NA
BP_GenderShared@result$up <- NA
BP_GenderShared@result$down <- NA
get_up_down_genes_GS <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC > 0 & deg_object$Female.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC < 0 & deg_object$Female.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_GenderShared <- get_up_down_genes_GS(MF_GenderShared, DEG_GenderShared)
BP_GenderShared <- get_up_down_genes_GS(BP_GenderShared, DEG_GenderShared)


p1 <- dotplot(BP_GenderShared, showCategory=18) + ggtitle("dotplot for BP enrichment")
p2 <- dotplot(MF_GenderShared, showCategory=18) + ggtitle("dotplot for MF enrichment")
plot_grid(p1, p2, ncol=2)



# Male
MF_Male@result$up <- NA
MF_Male@result$down <- NA
BP_Male@result$up <- NA
BP_Male@result$down <- NA
get_up_down_genes_Male <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Male.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_Male <- get_up_down_genes_Male(MF_Male, DEG_Male)
BP_Male <- get_up_down_genes_Male(BP_Male, DEG_Male)

p1 <- dotplot(BP_Male, showCategory=20, font.size=9) + ggtitle("dotplot for BP enrichment")
p2 <- dotplot(MF_Male, showCategory=20) + ggtitle("dotplot for MF enrichment")
plot_grid(p1, p2, ncol=2)


# Female
MF_Female@result$up <- NA
MF_Female@result$down <- NA
BP_Female@result$up <- NA
BP_Female@result$down <- NA
get_up_down_genes_female <- function(pathway_object, deg_object){
        for (i in 1:nrow(pathway_object)){
                pathway_genes <- pathway_object@result$geneID[i]
                pathway_genes <- as.character(strsplit(pathway_genes, "/")[[1]])
                up <- length(which(pathway_genes %in% deg_object[deg_object$Female.avg..logFC > 0,]$Gene.symbols))
                down <- length(which(pathway_genes %in% deg_object[deg_object$Female.avg..logFC < 0,]$Gene.symbols))
                pathway_object@result$up[i] <- up
                pathway_object@result$down[i] <- down
        }
        pathway_object@result$upDown <- ifelse(pathway_object@result$up >= pathway_object@result$down, "up", "down")
        pathway_object@result$upDown <- ifelse(pathway_object@result$up == pathway_object@result$down, "equal", pathway_object@result$upDown)
        return(pathway_object)
}

MF_Female <- get_up_down_genes_female(MF_Female, DEG_Female)
BP_Female <- get_up_down_genes_female(BP_Female, DEG_Female)


p1 <- dotplot(BP_Female, showCategory=20, font.size=9) + ggtitle("dotplot for BP enrichment")
p2 <- dotplot(MF_Female, showCategory=20) + ggtitle("dotplot for MF enrichment")
plot_grid(p1, p2, ncol=2)


png("~/figures/DEGs_enrichment_BP_MF_female_Splitted_upDown.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_Female, showCategory=20, font.size=9) + ggtitle("dotplot for BP enrichment (Female)") + facet_grid(.~upDown)
p2 <- dotplot(MF_Female, showCategory=20) + ggtitle("dotplot for MF enrichment (Female)") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2)
dev.off()

png("~/figures/DEGs_enrichment_BP_MF_male_Splitted_upDown.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_Male, showCategory=12, font.size=9) + ggtitle("dotplot for BP enrichment (Male)") + facet_grid(.~upDown)
p2 <- dotplot(MF_Male, showCategory=12) + ggtitle("dotplot for MF enrichment (Male)") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2)
dev.off()

png("~/figures/gender_shared_DEGs_enrichment_BP_MF_Splitted_upDown.png", width=16, height=8, units="in", res=300)
p1 <- dotplot(BP_GenderShared, showCategory=20, font.size=9) + ggtitle("dotplot for BP enrichment (Sex-Neutral)") + facet_grid(.~upDown)
p2 <- dotplot(MF_GenderShared, showCategory=16) + ggtitle("dotplot for MF enrichment (Sex-Neutral)") + facet_grid(.~upDown)
plot_grid(p1, p2, ncol=2)
dev.off()


save(BP_Male, MF_Male, BP_Female, MF_Female, BP_GenderShared, MF_GenderShared, 
        file="~/data/Pathway_Enrichment_Analysis_upDownDEGs.RData")

top_n <- 5
table_pathways <- rbind(head(BP_GenderShared, top_n), head(MF_GenderShared, top_n), head(BP_Female, top_n), 
        head(MF_Female, top_n), head(BP_Male, top_n), head(MF_Male, top_n))
table_pathways <- table_pathways[,c("ID", "Description", "p.adjust", "upDown")]
dim(table_pathways) # 30 4
write.table(table_pathways, file="~/data/top5_pathways_table.txt", sep="\t", row.names=F, quote=F)



#
# Check overlap of enriched pathways between male, female, and gender-shared DEGs
#

table(BP_Male$ID %in% BP_Female$ID)
# FALSE  TRUE
#    32     1
table(BP_Male$ID %in% BP_GenderShared$ID)
# FALSE  TRUE
#    30     3
table(BP_Female$ID %in% BP_GenderShared$ID)
# FALSE  TRUE
#   149     4

table(MF_Male$ID %in% MF_Female$ID)
# FALSE
#    20
table(MF_Male$ID %in% MF_GenderShared$ID)
# FALSE
#    20
table(MF_Female$ID %in% MF_GenderShared$ID)
# FALSE  TRUE
#     5     1

table(KEGG_Male$ID %in% KEGG_Female$ID) # Male have no KEGG pathway

table(KEGG_Male$ID %in% KEGG_GenderShared$ID) # Male have no KEGG pathway

table(KEGG_Female$ID %in% KEGG_GenderShared$ID)
# FALSE
#     9

na.omit(BP_Male[BP_Male$ID %in% BP_Female$ID,])
#                    ID                      Description GeneRatio   BgRatio      pvalue   p.adjust     qvalue           geneID Count
# GO:0022407 GO:0022407 regulation of cell-cell adhesion      3/13 468/28814 0.001078571 0.03379832 0.01634495 Klf4/Ccr5/Laptm5     3
na.omit(BP_Male[BP_Male$ID %in% BP_GenderShared$ID,])
#                    ID                                          Description GeneRatio  BgRatio       pvalue   p.adjust     qvalue    geneID Count
# GO:0032731 GO:0032731 positive regulation of interleukin-1 beta production      2/13 61/28814 0.0003387366 0.02736239 0.01323252 Egr1/Ccr5     2
# GO:0032732 GO:0032732      positive regulation of interleukin-1 production      2/13 73/28814 0.0004849649 0.02938079 0.01420862 Egr1/Ccr5     2
# GO:0045453 GO:0045453                                      bone resorption      2/13 78/28814 0.0005534620 0.03095130 0.01496812 Ctss/Car2     2
na.omit(BP_Female[BP_Female$ID %in% BP_GenderShared$ID,])
#                    ID                      Description GeneRatio   BgRatio       pvalue   p.adjust      qvalue                 geneID Count
# GO:0008360 GO:0008360         regulation of cell shape      4/42 157/28814 8.080397e-05 0.01093582 0.006795155 Aldoa/Myl12b/Rdx/Cfdp1     4
# GO:0042255 GO:0042255                ribosome assembly      3/42  64/28814 1.127748e-04 0.01160838 0.007213065       Rps19/Rplp0/Rps5     3
# GO:0022604 GO:0022604 regulation of cell morphogenesis      4/42 338/28814 1.465011e-03 0.02926847 0.018186459 Aldoa/Myl12b/Rdx/Cfdp1     4
# GO:0045216 GO:0045216  cell-cell junction organization      3/42 216/28814 3.845346e-03 0.04089495 0.025410763       Rdx/Itgb1/Ctnna1     3

na.omit(MF_Female[MF_Female$ID %in% MF_GenderShared$ID,])
#                    ID                    Description GeneRatio   BgRatio      pvalue   p.adjust     qvalue                 geneID Count
# GO:0050839 GO:0050839 cell adhesion molecule binding      4/43 307/28275 0.001204926 0.03735271 0.02600104 Mfge8/Rdx/Itgb1/Ctnna1     4


overlap = list(BP_Male = unique(BP_Male$ID), BP_Female = unique(BP_Female$ID), BP_GenderShared = unique(BP_GenderShared$ID))
flag_palette <- c("BP_Male" = "black",
                  "BP_Female" = "#FF883E",
                  "BP_GenderShared" = "#ce2b37")
png("~/figures/Venn_Overlapping_BioProc_thy21.png", width=8, height=6, units="in", res=300)
overlap %>% magrittr::extract(!(names(.) %in% c("bars", "stripes", "no sections"))) %>%
  # the palette has to be in the same order as the list of sets, so we reorder it
 venn::venn(zcolor = flag_palette[names(.)], cexil = 5, cexsn = 5)
dev.off()



#
# Check overlap of male and female enriched pathways between scRNAseq and Bulk RNAseq data
#

female_bulk_GO <- read.csv("~/data/Bulk_Female_TG_vs_WT_enrichGO_ALL.csv", stringsAsFactors=F)
dim(female_bulk_GO) # 52 10

male_bulk_GO <- read.csv("~/data/Bulk_Male_TG_vs_WT_enrichGO_ALL.csv", stringsAsFactors=F)
dim(male_bulk_GO) # 55 10



female_scRNA_GO <- rbind(as.data.frame(BP_Female), as.data.frame(MF_Female))
dim(female_scRNA_GO) # 159   9
male_scRNA_GO <- rbind(as.data.frame(BP_Male), as.data.frame(MF_Male))
dim(male_scRNA_GO) # 53  9

table(female_bulk_GO$ID %in% female_scRNA_GO$ID)
# FALSE  TRUE
#    49     3

table(male_bulk_GO$ID %in% male_scRNA_GO$ID)
# FALSE
#    55


na.omit(female_scRNA_GO[female_scRNA_GO$ID %in% female_bulk_GO$ID,])
#                    ID                        Description GeneRatio   BgRatio       pvalue   p.adjust     qvalue                geneID Count
# GO:0061041 GO:0061041        regulation of wound healing      3/42 133/28814 0.0009675587 0.02272621 0.01412131 Nfe2l2/Serpine2/Itgb1     3
# GO:1903034 GO:1903034 regulation of response to wounding      3/42 174/28814 0.0020899264 0.03123578 0.01940888 Nfe2l2/Serpine2/Itgb1     3
# GO:0030193 GO:0030193    regulation of blood coagulation      2/42  67/28814 0.0043188318 0.04384269 0.02724239       Nfe2l2/Serpine2     2

