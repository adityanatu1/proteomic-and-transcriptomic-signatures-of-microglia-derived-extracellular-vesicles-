
##############CHECK##############
#Comparison:Control VS LPS : Update your script name and comparison here + Double check the folder you are in to aviod mixing diff ex results!
rm(list = ls())
dev.off()
cat("\014")


##############CHECK##############
setwd("/Users/adityanatu/Library/Mobile Documents/com~apple~CloudDocs/Project_2023/Juliet_Data_Analysis/EV_Cell_RNASeq/5.BulkRNAseq_Responder Cells_RAW counts/LPS/")
getwd()
options(stringsAsFactors = FALSE)

library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)

#########################################################CHANGE##############################################################################
#Subsetting the cleanDat and traits file  was done in excel!
##############CHECK##############
cleanDat <- read.csv("LPS_BulkRNAseq_Responder Cells_RAW_counts.csv", header = TRUE, row.names = 1)
numericMeta <- read.csv("LPS_BulkRNAseq_Responder Cells_RAW_traits.csv", header = TRUE, row.names = 1)

colnames(cleanDat)
rownames(numericMeta)

##############CHECK##############
colnames(cleanDat) <- c("control_Evs_Tx_1","control_Evs_Tx_2","control_Evs_Tx_3","LPS_Evs_Tx_1","LPS_Evs_Tx_2","LPS_Evs_Tx_3")
#########################################################CHANGE##############################################################################

#sanity check 1 
all(colnames(cleanDat) %in% rownames(numericMeta))

#sanity check 2
all(colnames(cleanDat) == rownames(numericMeta))



#Deseq2 data matrix is needed for downstream analysis
dds <- DESeqDataSetFromMatrix(countData = cleanDat,
                              colData = numericMeta,
                              design = ~ Treatment)   



# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 15 (15 X 6 samples = 90) reads total across ALL the samples 
#check the threshold and how did you make this! 
keep <- rowSums(counts(dds)) >= 30 #check how! 
dds <- dds[keep,]
dds



########Symbol conversion for mouse IDs#######################################################################################
library(biomaRt)

#mart=useMart('ensembl')

#set up connection to ensembl mart for mouse
mouse<- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
tableSymbols<-getBM(attributes=c("mgi_symbol","ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(dds),mart=mouse) 

#sanity checks
dim(dds) #ControlVsDay28

dim(tableSymbols)
#put symbols in countData order

symbols<-tableSymbols$mgi_symbol[match(rownames(dds),tableSymbols$ensembl_gene_id)] #ControlVsDay28
length(symbols)
tail(symbols,100)

#replace NA with empty string where no lookup was possible.
symbols[is.na(symbols)]<-""
tail(symbols,100)

uniqueIDs<-paste0(symbols,"|",rownames(dds)) #ControlVsDay28
tail(uniqueIDs)
rownames(dds)<-uniqueIDs  #ControlVsDay28
####################################################################################################################
dds
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.csv(normalized_counts, "norm_counts.csv")
####################################################################################################################

#contrast <- c("condition", "level_to_compare", "base_level")

# set the factor level
#this allows us to set the default comparisons for control vs treatment condition! 
dds$Treatment <- relevel(dds$Treatment, ref = "Control")

#Run DESeq ----------------------
dds <- DESeq(dds)

################################################################################################################
##############CHECK############## : Here you need to snapshot the diff ex results that are printed to the screen - Make sure you grab it before buffer runs out!
#diff p adj value 
res05 <- results(dds) 
res05 #table w statistics

#write.csv(res05, "res05_withoutsettingALPHA.csv")

#diff ex analysis
summary(res05) #Summary deseq2 results 
dim(res05)

# contrasts
resultsNames(dds)


#up and down list created to pull out diff ex genes
up.1 <- res05[which(res05$log2FoldChange > 0 & res05$padj < .05),]
down.1 <- res05[which(res05$log2FoldChange < 0 & res05$padj < .05),]

#write.csv(up.1, file = "Up_dif_ex.csv")
#write.csv(down.1, file = "Down_dif_ex.csv")

#########################################################################CHANGE#################################################################################################

#diff p adj value 
#res10 <- results(dds, alpha = 0.10)
res10 <- results(dds, alpha = 0.10)
res10 #table w statistics

#diff ex analysis
summary(res10) #Summary deseq2 results 
dim(res10)
#write.csv(res10, "res10.csv")

# contrasts
resultsNames(dds)

#up and down list created to pull out diff ex genes
up.10 <- res10[which(res10$log2FoldChange > 0 & res10$padj < .10),]
down.10 <- res10[which(res10$log2FoldChange < 0 & res10$padj < .10),]

---------
  
up.10 <- res10[which(res10$log2FoldChange > 0 & res10$padj < .10),]
down.10 <- res10[which(res10$log2FoldChange < 0 & res10$padj < .10),]

#write.csv(up.10, file = "GO_Up_dif_ex.csv")
#write.csv(down.10, file = "GO_Down_dif_ex.csv")

####################################################VOLCANO############################################################
##############CHECK############## : Here you need to make sure that the comparison that you are running is accurately described in the volcano! 
#Verify what side is what using the resultsNames(dds) parameter for VOLC and Heatmap!
#FDR cutoff of alpha 
library(EnhancedVolcano)
#res05 = results parameter from deseq2 that needs to be plugged in here!

EnhancedVolcano(res05,
                lab = rownames(res05),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Treatment_LPSvsCtrl',
                pCutoff = 0.05,
                FCcutoff = 1,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 4.0,
                labSize = 6.0,
                colAlpha = 1,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0)

##############################################################################PCA##############################################################################
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup = c("Treatment"))

############################################################################################################################################################
AdjP <- as.data.frame(res05)
AdjP$padj[is.na(AdjP$padj)] <- 1
var <- order(AdjP$padj, decreasing=F)
matrix <- assay(rld)[var[1:50],]

#Set the thresh val here which you get from identifying the min and max range from your input matrix!
matrix <- matrix - rowMeans(matrix)
matrixx <- matrix
thresh_val <- min(abs(max(matrix)),abs(min(matrix)))
matrixx[matrixx < -thresh_val] <- -thresh_val
matrixx[matrixx > thresh_val] <- thresh_val

# select the 'contrast' you want
annotation_data <- as.data.frame(colData(rld)[c("Treatment")])
#########################################################CHANGE##############################################################################
##############CHECK#############: Change the main tag in Heatmap
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,cluster_cols= FALSE, main = "LPS Vs Control") #true / False to toggle supervised and unsupervised clustering
#########################################################CHANGE##############################################################################


#########################################################################IGNORE FOR NOW!###################################################################

# #Pull the norm counts! 
# VolcCounts <- counts(dds)
# VolcCounts.log2 <- log2(VolcCounts + 1)
# 
# write.csv(VolcCounts.log2,"VolcCounts_log2.csv")
# 
# #Parameters needed for ParANOVA! 
# cleanDat <- VolcCounts.log2
# Grouping <- numericMeta$Treatment
# parallelThreads=12
# source("./parANOVA.dex.R")
# ANOVAout <- parANOVA.dex()
# labelTop <- 25
# FCmin=1
# signifP=0.05 
# plotVolc()
# 
# head(ANOVAout)
######################################################################Volcano#######################################################################
#Sanity Check
# Find the rows that have a column with a value of 0
# zero_rows.1 <- VolcCounts.log2[rowSums(VolcCounts.log2 == 0) > 0,]
# # Print the zero rows
# print(sum(zero_rows.1))
######################################################################Volcano Code edit###################################################################

