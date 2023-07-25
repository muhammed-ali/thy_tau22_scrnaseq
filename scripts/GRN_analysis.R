#
# Add GRN network analysis for either all cell types combined (preferred version) or for one of the dominating cell types (less preferred)
#


# Female GRN (Global nominally significant DEGs; N=758)

DEG_Global_Poisson_F <- read.csv("~/data/DEG_Global_Poisson_Female.csv", header=T, row.names=1)
dim(DEG_Global_Poisson_F) # 896   5
head(DEG_Global_Poisson_F, 2)
DEA_pval <- DEG_Global_Poisson_F[DEG_Global_Poisson_F$p_val < 0.05,]
dim(DEA_pval) # 758   5
range(DEA_pval$avg_log2FC) # -0.7268129  0.6534306

nrow(DEG_Global_Poisson_F[DEG_Global_Poisson_F$p_val_adj < 0.05,]) # 541

nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.25,]) # 13
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.20,]) # 21
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.15,]) # 39
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.10,]) # 88
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.05,]) # 323

hist(abs(DEA_pval$avg_log2FC))

DEA_pval$Expression <- ifelse(DEA_pval$avg_log2FC > 0, 1, 0)
DEA_pval$Gene <- rownames(DEA_pval)
DEA_pval$Gene <- toupper(DEA_pval$Gene)
head(DEA_pval, 2)
#      p_val avg_log2FC pct.1 pct.2 p_val_adj Expression Gene
# Apoe     0  0.3516178 0.473 0.401         0          1 APOE
# Cst3     0  0.2823946 0.827 0.760         0          1 CST3

Expression <- DEA_pval[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 758   2

# write GRN input files
write.table(Expression, file="~/GRN/Female/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="~/GRN/Female/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="~/GRN/Female/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="~/GRN/Female/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)

head(Expression, 2)
Expression$Expression <- ifelse(Expression$Expression == 0, 1, 0)
write.table(Expression, file="~/GRN/Female/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Activation + Inhibitions only (no unspecified interactions).

cd ~/GRN/Female/

java -jar ~/GRN/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/GRN/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 498
wc -l NetworkPhenotype2.txt # 414

java -jar ~/GRN/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 61

java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 92


# NOTE: I tried network perturbation both with N1 (AD) and N2 (Control) phenotypes, reporting here the commands for N2 only.


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt

sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt


head -n 2 BruteForcePerturbationsUpdatedN1_*.txt

==> BruteForcePerturbationsUpdatedN1_1.txt <==
101     [JUN]
99      [RBX1]
==> BruteForcePerturbationsUpdatedN1_2.txt <==
107     [JUN, RBX1]
106     [JUN, UBB]
==> BruteForcePerturbationsUpdatedN1_3.txt <==
121     [JUN, MALAT1, FOS]
111     [JUN, FUS, RBX1]
==> BruteForcePerturbationsUpdatedN1_4.txt <==
128     [JUN, FUS, MALAT1, FOS]
124     [JUN, SMARCA5, MALAT1, FOS]


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 1 500000 BruteForcePerturbationsUpdatedN2_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 2 500000 BruteForcePerturbationsUpdatedN2_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 3 500000 BruteForcePerturbationsUpdatedN2_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 4 500000 BruteForcePerturbationsUpdatedN2_4.txt

sort -rn BruteForcePerturbationsUpdatedN2_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_1.txt
sort -rn BruteForcePerturbationsUpdatedN2_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_2.txt
sort -rn BruteForcePerturbationsUpdatedN2_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_3.txt
sort -rn BruteForcePerturbationsUpdatedN2_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_4.txt

head -n 2 BruteForcePerturbationsUpdatedN2_*.txt

==> BruteForcePerturbationsUpdatedN2_1.txt <==
83      [NR3C1]
68      [FOS]
==> BruteForcePerturbationsUpdatedN2_2.txt <==
89      [YBX1, NR3C1]
89      [NR3C1, FERMT2]
==> BruteForcePerturbationsUpdatedN2_3.txt <==
95      [HSP90AA1, YBX1, NR3C1]
95      [HSP90AA1, NR3C1, FERMT2]
==> BruteForcePerturbationsUpdatedN2_4.txt <==
102     [PKM, MAF, MALAT1, NR3C1]
97      [SFPQ, HSP90AA1, NR3C1, FERMT2]


# move everything in a new subdirectory
cd 
mkdir nominal_758DEGs
mv * ./nominal_758DEGs/



#
# Male GRN (Global nominally significant DEGs; N=530)
#

DEG_Global_Poisson_M <- read.csv("~/data/DEG_Global_Poisson_Male.csv", header=T, row.names=1)
dim(DEG_Global_Poisson_M) # 851   5
head(DEG_Global_Poisson_M, 2)


DEA_pval <- DEG_Global_Poisson_M[DEG_Global_Poisson_M$p_val < 0.05,]
dim(DEA_pval) # 530   5
range(DEA_pval$avg_log2FC) # -0.2546953  0.9496518

nrow(DEG_Global_Poisson_M[DEG_Global_Poisson_M$p_val_adj < 0.05,]) # 195

nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.25,]) # 3
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.20,]) # 5
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.15,]) # 13
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.10,]) # 28
nrow(DEA_pval[abs(DEA_pval$avg_log2FC) > 0.05,]) # 120

hist(abs(DEA_pval$avg_log2FC))



DEA_pval$Expression <- ifelse(DEA_pval$avg_log2FC > 0, 1, 0)
DEA_pval$Gene <- rownames(DEA_pval)
DEA_pval$Gene <- toupper(DEA_pval$Gene)
head(DEA_pval, 2)
#         p_val avg_log2FC pct.1 pct.2 p_val_adj Expression    Gene
# Ttr         0  0.9496518 0.135 0.025         0          1     TTR
# mt-Rnr2     0 -0.2227041 0.999 1.000         0          0 MT-RNR2

Expression <- DEA_pval[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 530   2

# write GRN input files
write.table(Expression, file="~/GRN/Male/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="~/GRN/Male/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="~/GRN/Male/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="~/GRN/Male/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)

head(Expression, 2)
Expression$Expression <- ifelse(Expression$Expression == 0, 1, 0)
write.table(Expression, file="~/GRN/Male/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Activation + Inhibitions only (no unspecified interactions).


cd ~/GRN/Male/

java -jar ~/GRN/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/GRN/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 264
wc -l NetworkPhenotype2.txt # 421

java -jar ~/GRN/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 74

java -jar PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 87


# NOTE: I tried network perturbation both with N1 (AD) and N2 (Control) phenotypes, reporting here the commands for N2 only.


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt


sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt


head -n 2 BruteForcePerturbationsUpdatedN1_*.txt

==> BruteForcePerturbationsUpdatedN1_1.txt <==
47      [CTNNB1]
41      [NR3C1]
==> BruteForcePerturbationsUpdatedN1_2.txt <==
52      [CTNNB1, IRF8]
49      [SPARC, CTNNB1]
==> BruteForcePerturbationsUpdatedN1_3.txt <==
54      [SPARC, IRF8, CTNNB1]
54      [IRF8, CTNNB1, ZBTB7A]
==> BruteForcePerturbationsUpdatedN1_4.txt <==
56      [SPARC, IRF8, CTNNB1, ZBTB7A]
56      [CD63, SPARC, IRF8, CTNNB1]


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 1 500000 BruteForcePerturbationsUpdatedN2_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 2 500000 BruteForcePerturbationsUpdatedN2_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 3 500000 BruteForcePerturbationsUpdatedN2_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression_N2.txt NetworkPhenotype2.txt 1 PerturbagenListGeneratorN2.txt 4 500000 BruteForcePerturbationsUpdatedN2_4.txt

sort -rn BruteForcePerturbationsUpdatedN2_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_1.txt
sort -rn BruteForcePerturbationsUpdatedN2_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_2.txt
sort -rn BruteForcePerturbationsUpdatedN2_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_3.txt
sort -rn BruteForcePerturbationsUpdatedN2_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN2_4.txt

head -n 2 BruteForcePerturbationsUpdatedN2_*.txt

==> BruteForcePerturbationsUpdatedN2_1.txt <==
119     [NR3C1]
97      [ABCB1A]
==> BruteForcePerturbationsUpdatedN2_2.txt <==
125     [JUN, NR3C1]
122     [NR4A1, CTNNB1]
==> BruteForcePerturbationsUpdatedN2_3.txt <==
129     [JUN, KLF4, NR3C1]
127     [ZEB2, JUN, NR3C1]
==> BruteForcePerturbationsUpdatedN2_4.txt <==
131     [JUN, SPARC, KLF4, NR3C1]
130     [ZEB2, JUN, NR3C1, KLF4]


# move everything in a new subdirectory
cd 
mkdir nominal_530DEGs
mv * ./nominal_530DEGs/



#
# Add GRN network analysis for one of the dominating cell types (less preferred)
#


load("~/data/TauCortex_DEGs_All10CTs.RData")

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


#
# Female GRN (Microglial)
#

microglial_female <- microglial_res_tauCortex[microglial_res_tauCortex$DEG.type == "female-specific",]
dim(microglial_female) # 68  7
max(microglial_female$Female.FDR) # 0.04734332

microglial_female$Expression <- ifelse(microglial_female$Female.avg..logFC > 0, 1, 0)
table(microglial_female$Expression)
#  0  1
# 32 36
microglial_female$Gene <- toupper(microglial_female$Gene.symbols)
Expression <- microglial_female[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 68  2
head(Expression, 2)
#   Gene Expression
# 5 C1QB          1
# 6 C1QA          1


# write GRN input files
write.table(Expression, file="~/GRN/Female_Microglial/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="~/GRN/Female_Microglial/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="~/GRN/Female_Microglial/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="~/GRN/Female_Microglial/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)

head(Expression_N2, 2)
#   Gene Expression
# 5 C1QB          0
# 6 C1QA          0

# Expression$Expression <- ifelse(Expression$Expression == 0, 1, 0)
write.table(Expression_N2, file="~/GRN/Female_Microglial/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional Interactions + Binding Interactions (Use for network building)


cd ~/GRN/Female_Microglial/

java -jar ~/GRN/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/GRN/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 37 interactions, 32 nodes
wc -l NetworkPhenotype2.txt # 40 interactions, 32 nodes

java -jar ~/GRN/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 5

java -jar PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 7


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt


sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt


head -n 2 BruteForcePerturbationsUpdatedN1_*.txt

==> BruteForcePerturbationsUpdatedN1_1.txt <==
15	[RPL3]
15	[DDX21]

==> BruteForcePerturbationsUpdatedN1_2.txt <==
20	[RPL3, TOP1]
20	[RPL3, HNRNPU]

==> BruteForcePerturbationsUpdatedN1_3.txt <==
20	[RPL3, TOP1]
20	[RPL3, HNRNPU]

==> BruteForcePerturbationsUpdatedN1_4.txt <==
20	[RPL3, TOP1]
20	[RPL3, HNRNPU]




#
# Female GRN (Endothelial)
#

endothelial_female <- endothelial_res_tauCortex[endothelial_res_tauCortex$DEG.type == "female-specific",]
dim(endothelial_female) # 41  7
max(endothelial_female$Female.FDR) # 0.04927151

endothelial_female$Expression <- ifelse(endothelial_female$Female.avg..logFC > 0, 1, 0)
table(endothelial_female$Expression)
#  0  1
# 15 26
endothelial_female$Gene <- toupper(endothelial_female$Gene.symbols)
Expression <- endothelial_female[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 41  2
head(Expression, 2)
#    Gene  Expression
# 5 PTPRB           0
# 6   UBB           1


# write GRN input files
write.table(Expression, file="~/GRN/Female_Endothelial/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="~/GRN/Female_Endothelial/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="~/GRN/Female_Endothelial/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="~/GRN/Female_Endothelial/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)

head(Expression_N2, 2)
#    Gene  Expression
# 5 PTPRB           1
# 6   UBB           0

# Expression$Expression <- ifelse(Expression$Expression == 0, 1, 0)
write.table(Expression_N2, file="~/GRN/Female_Endothelial/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional Interactions + Binding Interactions (Use for network building).


cd ~/GRN/Female_Endothelial/

java -jar ~/GRN/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/GRN/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 16 interactions, 16 nodes
wc -l NetworkPhenotype2.txt # 15 interactions, 16 nodes

java -jar ~/GRN/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 2

java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 0


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt


sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt


head -n 2 BruteForcePerturbationsUpdatedN1_*.txt
==> BruteForcePerturbationsUpdatedN1_1.txt <==
10      [CHCHD2]
9       [SFPQ]

==> BruteForcePerturbationsUpdatedN1_2.txt <==
12      [CHCHD2, SFPQ]
10      [CHCHD2]

==> BruteForcePerturbationsUpdatedN1_3.txt <==
12      [CHCHD2, SFPQ]
11      [SFPQ, CHCHD2, ITGA1]

==> BruteForcePerturbationsUpdatedN1_4.txt <==
12      [CHCHD2, SFPQ]
11      [SFPQ, CHCHD2, ITGA1]



#
# Add GRN network analysis for global male and female DEGs obtained using "gender_spec_genes' function
#

#
# Female (Gender-specific DEGs) GRN
#

load("~/data/Global_DEGs.RData")

dim(global_gender_spec_genes) # 175   6
table(global_gender_spec_genes$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#              52               65               44               14


global_female <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "female-specific",]
dim(global_female) # 41  7
max(global_female$Female.FDR) # 0.04141131

global_female$Expression <- ifelse(global_female$Female.avg..logFC > 0, 1, 0)
table(global_female$Expression)
#  0  1
# 14 38
global_female$Gene <- toupper(global_female$Gene.symbols)
Expression <- global_female[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 52  2
head(Expression, 2)
#     Gene  Expression
# 15 H3F3B           1
# 16   UBB           1


# write GRN input files
write.table(Expression, file="~/GRN/Female/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="~/GRN/Female/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="~/GRN/Female/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="~/GRN/Female/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)

head(Expression_N2, 2)
#     Gene  Expression
# 15 PTPRB           0
# 16   UBB           0

# Expression$Expression <- ifelse(Expression$Expression == 0, 1, 0)
write.table(Expression_N2, file="~/GRN/Female/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional Interactions + Binding Interactions (Use for network building).


cd ~/GRN/Female/

java -jar ~/GRN/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/GRN/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 83 interactions, 26 nodes
wc -l NetworkPhenotype2.txt # 94 interactions, 31 nodes

java -jar ~/GRN/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 19

java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 13


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt


sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt


head -n 2 BruteForcePerturbationsUpdatedN1_*.txt

==> BruteForcePerturbationsUpdatedN1_1.txt <==
14      [UBB]
8       [RDX]

==> BruteForcePerturbationsUpdatedN1_2.txt <==
19      [UBB, RPS19]
17      [UBB, ALDOA]

==> BruteForcePerturbationsUpdatedN1_3.txt <==
20      [UBB, RPS19, ALDOA]
19      [UBB, RPS19]

==> BruteForcePerturbationsUpdatedN1_4.txt <==
20      [UBB, RPS19, CNBP, ALDOA]
20      [UBB, RPS19, ALDOA]



#
# Male (Gender-specific DEGs) GRN
#


load("~/data/Global_DEGs.RData")

dim(global_gender_spec_genes) # 175   6
table(global_gender_spec_genes$DEG.type)
# female-specific gender-dimorphic    gender-shared    male-specific
#              52               65               44               14


global_male <- global_gender_spec_genes[global_gender_spec_genes$DEG.type == "male-specific",]
dim(global_male) # 14  6
max(global_male$Male.FDR) # 0.03976319

global_male$Expression <- ifelse(global_male$Male.avg..logFC > 0, 1, 0)
table(global_male$Expression)
# 0  1
# 4 10
global_male$Gene <- toupper(global_male$Gene.symbols)
Expression <- global_male[,c("Gene", "Expression")]
Expression <- unique(Expression)
Expression <- na.omit(Expression)
dim(Expression) # 14  2
head(Expression, 2)
#     Gene  Expression
# 15  CTSS           1
# 16  KLF4           0


# write GRN input files
write.table(Expression, file="~/GRN/Male/expression.txt", sep="\t", row.names=F, quote=F, col.names=F)
write.table(Expression$Gene, file="~/GRN/Male/geneList.txt", sep="\t", row.names=F, quote=F, col.names=F)

write.table(Expression, file="~/GRN/Male/N1_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)
Expression_N2 <- Expression
Expression_N2$Expression <- ifelse(Expression_N2$Expression < 1, 1, 0)
write.table(Expression_N2, file="~/GRN/Male/N2_nodeColorMapping.txt", sep="\t", row.names=F, quote=F, col.names=F)

head(Expression_N2, 2)
#     Gene  Expression
# 15  CTSS           0
# 16  KLF4           1

# Expression$Expression <- ifelse(Expression$Expression == 0, 1, 0)
write.table(Expression_N2, file="~/GRN/Male/expression_N2.txt", sep="\t", row.names=F, quote=F, col.names=F)


# Metacore Filter applied are:
# Species = Mus Musculus
# Interaction type = Binding, CrT, Regulation, influence on expression, transcriptional regulation
# Additional Filter: Functional Interactions + Binding Interactions (Use for network building).


cd ~/GRN/Male

java -jar ~/GRN/JARs/Preprocessor.jar . nodemap.txt interactions.txt geneList.txt
java -jar ~/GRN/JARs/DifferentialNetworkAnalysis.jar expression.txt adjacency.txt GAResult.txt 0 true 1000 100 .

wc -l NetworkPhenotype1.txt # 5 interactions, 4 nodes
wc -l NetworkPhenotype2.txt # 9 interactions, 9 nodes

java -jar ~/GRN/JARs/CommonNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt CommonNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/DifferentialNetworkGenerator.jar NetworkPhenotype1.txt NetworkPhenotype2.txt DifferentialNetworkGenerator_Output.txt
java -jar ~/GRN/JARs/ComputeCycles.jar CommonNetworkGenerator_Output.txt expression.txt pos.txt neg.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype1.txt 1 SteadyStateCalculatorN1.txt
java -jar ~/GRN/JARs/SteadyStateCalculator.jar expression.txt NetworkPhenotype2.txt 2 SteadyStateCalculatorN2.txt
java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN1.txt PerturbagenListGeneratorN1.txt

wc -l PerturbagenListGeneratorN1.txt # 8

java -jar ~/GRN/JARs/PerturbagenListGenerator.jar pos.txt neg.txt DifferentialNetworkGenerator_Output.txt SteadyStateCalculatorN2.txt PerturbagenListGeneratorN2.txt

wc -l PerturbagenListGeneratorN2.txt # 2


java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 1 500000 BruteForcePerturbationsUpdatedN1_1.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 2 500000 BruteForcePerturbationsUpdatedN1_2.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 3 500000 BruteForcePerturbationsUpdatedN1_3.txt
java -jar ~/GRN/JARs/BruteForcePerturbationsUpdated.jar expression.txt NetworkPhenotype1.txt 1 PerturbagenListGeneratorN1.txt 4 500000 BruteForcePerturbationsUpdatedN1_4.txt


sort -rn BruteForcePerturbationsUpdatedN1_1.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_1.txt
sort -rn BruteForcePerturbationsUpdatedN1_2.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_2.txt
sort -rn BruteForcePerturbationsUpdatedN1_3.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_3.txt
sort -rn BruteForcePerturbationsUpdatedN1_4.txt > Z.txt
mv Z.txt BruteForcePerturbationsUpdatedN1_4.txt


head -n 2 BruteForcePerturbationsUpdatedN1_*.txt

==> BruteForcePerturbationsUpdatedN1_1.txt <==
2       [KLF4]
2       [EGR1]

==> BruteForcePerturbationsUpdatedN1_2.txt <==
2       [VIM, KLF4]
2       [KLF4]

==> BruteForcePerturbationsUpdatedN1_3.txt <==
2       [VIM, KLF4]
2       [KLF4]

==> BruteForcePerturbationsUpdatedN1_4.txt <==
2       [VIM, KLF4]
2       [KLF4]