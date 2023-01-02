library("tximport")
library("readr")
pwd()
library("tximportData")
library("DESeq2")
library(GenomicFeatures)

names(files) <- samples$Name
# STEP 1: get transcript file to map tx to gene (R studio runs in home directory ~)
# tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
txdb <- makeTxDbFromGFF(file="gencode.v30.annotation.gff3.gz")
txdb
saveDb(x=txdb, file = "gencode.v30.annotation.TxDb")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
write.table(tx2gene, "tx2gene.gencode.v30.csv", sep = "\t", row.names = FALSE)
head(k)
head(tx2gene)
dim(tx2gene)
length(k)
tx2gene <- read_csv("tx2gene.gencode.v30.csv")
#tx2gene <- read.table("tx2gene.gencode.v30.csv", Header=T)

# STEP 2: import transcripts from quants directory (output of run_STAR.sh)
dir <- "/home/chienlab/Alan/quants" #this directory is for the location of samples.txt
samples <- read.table(file.path(dir,"samples.csv"), header=T)
# rownames(samples) <- samples$Run
# samples[,c("Run","cell_ID","condition","Group")]
files <- file.path(dir, samples$Run, "quant.sf")
txi <- tximport (files, type = "salmon", tx2gene = tx2gene)

rownames(samples) <- colnames(txi$counts)

# dds$culture.nested = factor(rep(rep(1:5,each=2),2))
#ddsTxi <- DESeqDataSetFromTximport(txi,colData=samples,design = ~condition )
#dds <- DESeq(ddsTxi)
#results <- results(dds)
#resultsNames(dds)
# subgroup <- factor(rep(1:10),each=3)

#create a factor "group" that contains the merged cell_ID and condition fields
# and use this as the design factor for DESeq2 analysis
dds$group <- factor(paste0(dds$cell_ID, dds$condition))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

H460 <- results(dds, contrast=c("group","H4605","H4600")) 
write.csv(H460,"~/Alan/VG_2019_RNASeq/pg545_H460-2.csv")

OVCAR5A <- results(dds, contrast=c("group","OVCAR55","OVCAR50")) 
write.csv(OVCAR5A,'~/Alan/VG_2019_RNASeq/pg545_OVCAR5-2.csv')

OVCAR5B <- results(dds, contrast=c("group","OVCAR5NTC","OVCAR5shRNA")) 
write.csv(OVCAR5B,'~/Alan/RNASeq-analysis/OVCAR5_shRNA.csv')

OVCAR7 <- results(dds, contrast=c("group","OVCAR7bl","OVCAR7ev")) 
write.csv(OVCAR7,'~/Alan/RNASeq-analysis/OVCAR7_bl-ev.csv')

OVCAR8 <- results(dds, contrast=c("group","OVCAR8ctr","OVCAR8ko")) 
write.csv(OVCAR8,'~/Alan/RNASeq-analysis/OVCAR8_ko.csv')

#BiocManager::install("org.Hs.eg.db")
library("pheatmap")
source("~/Alan/VG_2019_RNASeq/Setup.R")
source("~/Workspace/scripts/support_func.R")
source("~/Workspace/scripts/mk_pheatmap.R")
source("~/Workspace/scripts/mk_volcano.R")  
source("~/Workspace/scripts/mk_ma.R")  
library("AnnotationDbi")
#library("ensembldb")
library("AnnotationFilter")
#library("EnsDb.Hsapiens.v86")
library("BiocManager")

# BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")

matrix_H460 <- mk_matrix(H460,rlog_data,0.05)
dataFrame <- colnames(matrix_H460)
mk_pheatmap(matrix_H460,dataFrame,rownames(matrix_LH460),0.05,"TEST")

#relevant column (experiment) conditions for these experiments
pg545_dataFrames <- colData(rlog_data)[1:12,c("Name")]
LRRC15_dataFrames <- colData(rlog_data)[13:30,c("Name")]

# get the top 1000 (smallest pAdj) common to all three LRRC15 experiments
# Unfortunately, this heatmap highlights overlaps with high p-value that do not always agree 
# in terms of the direction of up/down regulation, ie a gene may be highly upregulated by 
# LRRC15 KO and highly downregulated by LRRC15 shRNA and it would make this list. 
OVCAR5B_top1k <- head(order(OVCAR5B$padj),1000)
OVCAR7_top1k <- head(order(OVCAR7$padj),1000)
OVCAR8_top1k <- head(order(OVCAR8$padj),1000)
LRRC_top_pvalue <- intersect(OVCAR5B_top1k,OVCAR7_top1k)
LRRC_top_pvalue <- intersect(LRRC_top_pvalue,OVCAR8_top1k)
LLRRC_top_pvalue # has 43 that are in all three! (?)
LRRC_top <- head(LRRC_top_pvalue,100) # reduce to top 100

#pheatmap(OVCAR_matrix, annotation_color=LRRC_top_pvalue, 
#         labels_col=LRRC15_dataFrames,cellheight = 4,cellwidth = 8,    fontsize = 4)
# mk_matrix(OVCAR_matrix,dds,0.05,LRRC_top_pvalue)

LRRC_matrix <- assay(rlog_data)[LRRC_top,] # select rows that match top p-values
LRRC_matrix <- LRRC_matrix - rowMeans(LRRC_matrix)
rownames(LRRC_matrix) <- map_symbols_adv(rownames(LRRC_matrix)) #maps symbols

LRRC15_dataFrames <- colData(rlog_data)[13:30,c("Name")] # LRRC15 columns
pheat_LRRC <- LRRC_matrix[,LRRC15_dataFrames]        # limit columns to LRRC15 experiments
pdf("LRRC15_intersection_top_pvalues.pdf")
pheatmap(pheat_LRRC)
dev.off()

# with Jeremy

cell_ID <- rep("H460", 6)
treatment <- c(0,0,0,5,5,5)
df <- cbind(cell_ID, treatment)
# sample_names <- colnames(dds)[1:6]
#  sample_names <- c("H460_0_1", "H460_0_2", "H460_0_3", "H460_5_1", "H460_5_2", "H460_5_3")
sample_names = paste(samples$cell_ID[1:12], samples$condition[1:12], "uM")
H460_rld <- rlog_data[,1:6] # columns 1-6 are pg545
mat <- mk_matrix(H460,H460_rld,p,)
pdf()
mk_pheatmap(mat, df, sample_names, "PG545 Treatment")
dev.off()

cell_ID <- rep("OVCAR5A", 6)
treatment <- c(0,0,0,5,5,5)
df <- cbind(cell_ID, treatment)
sample_names = paste(samples$cell_ID[7:12], samples$condition[7:12], "uM")
OVCAR5A_rld <- rlog_data[,7:12] # columns 1-6 are pg545
mat <- mk_matrix(OVCAR5A,OVCAR5A_rld,p,)
pdf("pg545_OVCAR5.pdf")
mk_pheatmap(mat, df, sample_names, "PG545 Treatment")
dev.off()


# do both pg545 experiments together in one heatmap
# both pg545 experiments use 0 and 5 in the condition column, so this can be used in contrast()
# requires new dds object created from design using only condition field
ddsCondTxi <- DESeqDataSetFromTximport(txi,colData=samples,design = ~condition )
dds_cond <- DESeq(ddsCondTxi)
results_cond <- results(dds_cond)
resultsNames(dds_cond)
pg545_results <- results(dds_cond,name="condition_5_vs_0")
write.csv(pg545_results,'~/Alan/VG_2019_RNASeq/pg545_bothExperiments.csv')

cell_ID <- c( rep("H460", 6),rep("OVCAR5",6)) 
treatment <- rep(c(0,0,0,5,5,5),2)
df <- cbind(cell_ID, treatment)
sample_names = paste(samples$cell_ID[1:12], samples$condition[1:12], "uM")
pg545_rld <- rlog_data[,1:12] # columns 1-6 are pg545
pg545_matrix <- mk_matrix(pg545_results,pg545_rld,0.05,)
pdf()
mk_pheatmap(pg545_matrix,df, sample_names, "PG545 H460+OVCAR5")
dev.off()

# do combined LRRC15 experiments using LRRC_top_pvalue intersection above
LRRC_top <- head(LRRC_top_pvalue,100) # top 100 p-value
pheat_LRRC <- OVCAR_matrix[,LRRC15_dataFrames] 
sample_names = paste(samples$cell_ID[13:30], samples$condition[13:30])

#top_rows <- assay(rlog_data)[LRRC_top, ]                  # creates gene matrix based on position
#top_rows <- top_rows - rowMeans(top_rows)                 # better highlight expression of genes
#rownames(top_rows) <- map_symbols_adv(rownames(top_rows)) # maps symbols

gene_names <- map_symbols_adv(rownames(rlog_data))

LRRC_topGenes <- gene_names[LRRC_top]

LRRC_matrix <- mk_matrix(dds,rlog_data[,13:30],,LRRC_topGenes)

LRRC15_df <- as.character(samples$cell_ID[13:30])
treatment <- rep( c( rep("reduced expression", 3),rep("increased expression",3)), 3)
df <- cbind(LRRC15_df, treatment)

pdf("LRRC_top_combined.pdf")
mk_pheatmap(rlog_data[,[13:30]],df, sample_names, "")
dev.off()

LRRC15_df <- as.character(samples$cell_ID[13:30])
treatment <- rep( c( rep("reduced expression", 3),rep("increased expression",3)), 3)
df <- cbind(LRRC15_df, treatment)

# Do OVCAR5B shRNA knockdown
cell_ID <- rep("OVCAR5", 6)
treatment <- c( rep("control",3),rep("shRNA",3))
ov5_df <- cbind(cell_ID, treatment)
sample_names = paste(samples$cell_ID[13:18], samples$condition[13:18])
OVCAR5_rld <- rlog_data[,13:18] # columns for OVCAR5shRNA
mat <- mk_matrix(OVCAR5B,OVCAR5_rld,p,) # created from results() above
pdf("LRRC15_shRNA_OVCAR5.pdf")
mk_pheatmap(mat, ov5_df, sample_names, "LRRC15 shRNA OVCAR5")
dev.off()

OVCAR5B <- results(dds, contrast=c("group","OVCAR5NTC","OVCAR5shRNA")) 
pdf("volcano_LRRC15_shRNA_OVCAR5.pdf")
mk_volcano(OVCAR5B,59000,0.01,,2,,,,"LRRC15 shRNA Gene Expression Changes in OVCAR5 Cells","Treatment/Control log2 fold-change")
dev.off()

# Do OVCAR7 overexpression experiment only 

cell_ID <- rep("OVCAR7", 6)
treatment <- c( rep("overExp",3),rep("emptyVec",3))
ov7_df <- cbind(cell_ID, treatment)
sample_names = paste(samples$cell_ID[19:24], samples$condition[19:24])
OVCAR7_rld <- rlog_data[,19:24] # columns for OVCAR7
mat7 <- mk_matrix(OVCAR7,OVCAR7_rld,p,) # created from results() above
pdf()
mk_pheatmap(mat7, ov7_df, sample_names, "OVCAR7")
dev.off()

# show all LRRC15 data using OVCAR7 pvalue subset to choose genes to display 
LRRC_rld <-rlog_data[,13:30]
matLRRC <- mk_matrix(OVCAR7,LRRC_rld,0.05,)
sample_names = paste(samples$cell_ID[13:30],samples$condition[13:30])
pdf("LRRC15_combined_OE-pvalues.pdf")
mk_pheatmap(matLRRC,ov7_df,sample_names, "LRRC15 all experiments (OverExpn pvalue)")
dev.off()

# Do OVCAR8
cell_ID <- rep("OVCAR8", 6)
treatment <- c( rep("control",3),rep("KO",3))
ov8_df <- cbind(cell_ID, treatment)
sample_names = paste(samples$cell_ID[25:30], samples$condition[25:30])
OVCAR8_rld <- rlog_data[,25:30] # columns for OVCAR8
mat <- mk_matrix(OVCAR8,OVCAR8_rld,p,) # created from results() above
pdf("heat_LRRC15-KO_OVCAR8.pdf")
mk_pheatmap(mat, ov8_df, sample_names, "OVCAR8")
dev.off()

# show all LRRC15 data using OVCAR8 pvalue subset to choose genes to display 
LRRC_rld <-rlog_data[,13:30]
matLRRC <- mk_matrix(OVCAR8,LRRC_rld,0.05,)
sample_names = paste(samples$cell_ID[13:30], samples$condition[13:30])
pdf("LRRC15_combined_OVCAR8_KO-pvalues.pdf")
mk_pheatmap(matLRRC,ov8_df,sample_names, "LRRC15 all experiments (OVCAR8 KO pvalues)")
dev.off()

# show all LRRC15 data using OVCAR5B pvalue subset to choose genes to display 
LRRC_rld <-rlog_data[,13:30]
matLRRC <- mk_matrix(OVCAR5B,LRRC_rld,0.05,)
sample_names = paste(samples$cell_ID[13:30], samples$condition[13:30])
pdf("LRRC15_combined_shRNA-pvalues.pdf")
mk_pheatmap(matLRRC,ov5_df,sample_names, "LRRC15 all experiments (OVCAR5 shRNA pvalues)")
dev.off()

H460_rev <- results(dds, contrast=c("group","H4605","H4600")) 
pdf("volcano_pg545_H460.pdf")
mk_volcano(H460_rev,30,0.05,,2,,,,"pg545-induced Gene Expression Changes in H460 Cells","pg545-treated/Control log2 fold-change")
dev.off()

OVCAR5A_rev <- results(dds, contrast=c("group","OVCAR55","OVCAR50")) 
pdf("volcano_pg545_OVCAR5.pdf")
mk_volcano(OVCAR5A_rev,30,0.01,,2,,,,"pg545-induced Gene Expression Changes in OVCAR5 Cells","pg545-treated/Control log2 fold-change")
dev.off()

# combined pg545 from above
pdf("volcano_pg545_combined.pdf")
mk_volcano(pg545_results,30,0.01,,2,,,,"pg545 Gene Expression Changes in H460 and OVCAR5 Cells","Control/pg545-treated (inverted) log2 fold-change")
dev.off()

OVCAR7 <- results(dds, contrast=c("group","OVCAR7bl","OVCAR7ev")) 
pdf("volcano_LRRC15_OverExp_OVCAR7.pdf")
mk_volcano(OVCAR7,50,0.01,,2,,,,"LRRC15 OverExpn: Gene Expression Changes in OVCAR7 Cells","Overexpression/Control log2 fold-change")
dev.off()

OVCAR8 <- results(dds, contrast=c("group","OVCAR8ctr","OVCAR8ko")) 
pdf("volcano_LRRC15_KO_OVCAR8.pdf")
mk_volcano(OVCAR7,40,0.01,,2,,,,"LRRC15 KO vs control: Gene Expression Changes in OVCAR8 Cells","LRRC15 KO/Control log2 fold-change")
dev.off()

# group shRNA,KO and OE-control conditions together as one group
# This correctly groups low expression versus high expression conditions 
# across shRNA,KO, and OverExpn experiments (ie, OverExpn control is grouped with shRNA + KO treatment)
low_high_factor <- factor(c( rep('na',12), rep( c( rep('high',3),rep('low',3)), 3) ))
samples2 <- cbind(samples,low_high)
ddsCondTxi <- DESeqDataSetFromTximport(txi,colData=samples2,design = ~low_high )
dds_lowHigh <- DESeq(ddsCondTxi)
resultsNames(dds)
LRRC15_combined_results <- results(dds_lowHigh,name="low_high_low_vs_high")
write.csv(LRRC15_combined_results,'~/Alan/VG_2019_RNASeq/LRRC15_combined_low-versus-high-expn.csv')

pdf("volcano_LRRC15_combined.pdf")
mk_volcano(LRRC15_combined_results,40,0.05,,2,,,,"LRRC15 Combined: Gene Expression Changes in OVCAR5/7/8 Cells","Treatment/Control log2 fold-change")
dev.off()

LRRC_rld <-rlog_data[,13:30]
matLRRC <- mk_matrix(LRRC15_combined_results,LRRC_rld,0.05,)
sample_names = paste(samples$cell_ID[13:30], samples$condition[13:30])
pdf("LRRC15_combined_low-high.pdf")
mk_pheatmap(matLRRC,ov7_df,sample_names, "LRRC15 all experiments high-vs-low")
dev.off()

library(dplyr)
library(data.table)
LRRC15_combined_results[ue$]
LRRC15_combined_results[ rownames(LRRC15_combined_results) %like% "ENSG00000164458",]
LRRC15_combined_results[ rownames(LRRC15_combined_results) == "ENSG00000164458.9",]
matLRRC[rownames(matLRRC) == "TBXT",]

#  mk_MAplot <- function(res, n, p, selectLabs, unit, title, subtitle, ylim, by){

H460_rev <- results(dds, contrast=c("group","H4605","H4600")) 
pdf("MAplot-pg545-H460.pdf")
mk_MAplot(H460_rev,100,0.05,,"pg545 / control","pg545 H460 cells")
dev.off()

OVCAR5A_rev <- results(dds, contrast=c("group","OVCAR55","OVCAR50")) 
pdf("MAplot-pg545-OVCAR5.pdf")
mk_MAplot(OVCAR5A_rev,100,0.05,,"pg545 / control","pg545 OVCAR5 cells")
dev.off()

# both pg545 experiments use 0 and 5 in the condition column, so this can be used in contrast()
# requires new dds object created from design using only condition field
ddsCondTxi <- DESeqDataSetFromTximport(txi,colData=samples,design = ~condition )
dds_cond <- DESeq(ddsCondTxi)
results_cond <- results(dds_cond)
resultsNames(dds_cond)
pg545_results <- results(dds_cond,name="condition_5_vs_0")
pdf("MAplot-pg545-combined-top100.pdf")
mk_MAplot(pg545_results,100,0.05,,"pg545 / control","pg545 combined results: H460+OVCAR5 cells")
dev.off()

OVCAR5B <- results(dds, contrast=c("group","OVCAR5shRNA","OVCAR5NTC")) 
pdf("MAplot_LRRC15_shRNA_OVCAR5.pdf")
mk_MAplot(OVCAR5B,50,0.05,,"LRRC15 shRNA / control","LRRC15 shRNA in OVCAR5 cells")
dev.off()

OVCAR7 <- results(dds, contrast=c("group","OVCAR7bl","OVCAR7ev")) 
pdf("MAplot_LRRC15_OverExp_OVCAR7.pdf")
mk_MAplot(OVCAR7,50,0.05,,"Overexpn / EmptyVec","LRRC15 OverExpn in OVCAR7 cells")
dev.off()

OVCAR8 <- results(dds, contrast=c("group","OVCAR8ctr","OVCAR8ko")) 
pdf("MAplot_LRRC15_KO_OVCAR8.pdf")
mk_MAplot(OVCAR8,50,0.05,,"Knockout / control","LRRC15 KO vs control in OVCAR8 Cells")
dev.off()

low_high_factor <- factor(c( rep('na',12), rep( c( rep('high',3),rep('low',3)), 3) ))
samples2 <- cbind(samples,low_high)
ddsCondTxi <- DESeqDataSetFromTximport(txi,colData=samples2,design = ~low_high )
dds_lowHigh <- DESeq(ddsCondTxi)
resultsNames(dds)
LRRC15_combined_results <- results(dds_lowHigh,name="low_high_low_vs_high")
pdf("MAplot_LRRC15_combined-results.pdf")
mk_MAplot(LRRC15_combined_results,50,0.05,,"Low expression / High expression","LRRC15 combined results OVCAR5/7/8 cells")
dev.off()

#
#  Look up values for individual genes (investigate/confirm top hits)
LRRC15_combined_results[ rownames(LRRC15_combined_results) %like% "ENSG00000164458",]
LRRC15_combined_results[ rownames(LRRC15_combined_results) == "ENSG00000164458.9",]
matLRRC[rownames(matLRRC) == "ZG16",]
matLRRC[rownames(matLRRC) == "OLR1",]
# all results
results[rownames(results) %like% "ENSG00000125968",]
matLRRC[rownames(matLRRC) == "PTPN7",]

filtered_8 <- head(res[order(OVCAR8$padj),],100)
filtered_8 <- dplyr::filter(as.data.frame(labels), log2FoldChange >= 20)
head_ovcar8[ head_ovcar8$log2FoldChange < -20, ]
# get gene name of ensembl gene ID
map_symbols_adv( rownames(head_ovcar8[ head_ovcar8$log2FoldChange < -20, ]))
head_ovcar8[ head_ovcar8$log2FoldChange > 2, ]

comb <- LRRC15_combined_results
upreg_comb <- comb[ !is.na(comb$log2FoldChange) & comb$log2FoldChange > 3 & !is.na(comb$padj) & comb$padj < 0.05, ]
map_symbols_adv( rownames(upreg_hits))
downreg_comb <- comb[ !is.na(comb$log2FoldChange) & comb$log2FoldChange < -3 & !is.na(comb$padj) & comb$padj < 0.05, ]
map_symbols_adv( rownames(downreg_comb))

# Viji request email for heatmap Jan 8th 2020
# saved .xls as .csv
# run "perl -pe 's/,.*\Z//' OVCAR5_Viji_genelist.csv > Viji_genelist.txt"
# run perl script convert_to_ENSG.pl to get ENSG prefix from gene name
fileList <- file("~/Alan/VG_2019_RNASeq/Viji_ENSG_list.txt")
geneList <- strsplit(readLines(fileList),"\n")
close(fileList)

#OVCAR5A <- results(dds, contrast=c("group","OVCAR55","OVCAR50")) 
# write.csv(OVCAR5A,'~/Alan/VG_2019_RNASeq/pg545_OVCAR5-2.csv')
# add gene names to file with perl script
# OVCAR5A_named <- read.csv('~/Alan/VG_2019_RNASeq/pg545_OVCAR5-2_named.csv', header=TRUE)

# remove ENSG suffix on OVCAR5A_rld to allow matching with ENSG prefix names
OVCAR5A_rld <- rlog_data[,7:12] # columns 1-6 are pg545
all_ENSG <- rownames(OVCAR5A_rld)
OVCAR5A_rlog_ENSGprefix <- sub( "(\\d+)\\.(\\d+)", "\\1",all_ENSG, perl=TRUE)
rownames(OVCAR5A_rld) <- OVCAR5A_rlog_ENSGprefix

cell_ID <- rep("OVCAR5A", 6)
treatment <- c(0,0,0,5,5,5)
df <- cbind(cell_ID, treatment)
sample_names = paste(samples$cell_ID[7:12], samples$condition[7:12], "uM")

# mk_matrix modified to not filter when p=0 (parameter 3)
mat <- mk_matrix(OVCAR5A_named,OVCAR5A_rld,0,geneList)

pdf("pg545_OVCAR5_Viji-geneList-Jan08-2020-2.pdf",5,28)
mk_pheatmap(mat, df, sample_names, "PG545 Treatment")
dev.off()

write.csv(mat,"~/Alan/VG_2019_RNASeq/pg545_geneList2.csv")


