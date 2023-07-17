
##############CHECK##############
#Comparison:Control VS TGFB : Update your script name and comparison here + Double check the folder you are in to avoid mixing diff ex results!
rm(list = ls())
dev.off()
cat("\014")


##############CHECK##############
setwd("/Users/adityanatu/Library/Mobile Documents/com~apple~CloudDocs/Project_2023/Juliet_Data_Analysis/EV_RNASeq/2.mRNAseq_EVs_RAW counts/TGFB/")
getwd()
options(stringsAsFactors = FALSE)

library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)

#########################################################CHANGE##############################################################################
#Subsetting the cleanDat and traits file  was done in excel!
##############CHECK##############
cleanDat <- read.csv("TGFB_mRNAseq_EVs_RAW counts.csv", header = TRUE)


#Need Unique rownames! Setting them now!

dat <-cleanDat
dim(dat)
df <-dat

#Making rownams unique
rownames(df) <- make.names(dat$Geneid, unique = TRUE)
write.csv(df, file = "Unique_TGFB_mRNAseq_EVs_RAW counts.csv")


df <- read.csv("Unique_TGFB_mRNAseq_EVs_RAW counts.csv", header = TRUE, row.names = 1)[,-c(1)]
cleanDat <- df
dim(cleanDat)

###Remove var that arent needed!
rm(dat)
rm(df)


numericMeta <- read.csv("TGFB_mRNAseq_EVs_RAW_Traits.csv", header = TRUE, row.names = 1)

colnames(cleanDat)
rownames(numericMeta)

##############CHECK##############
colnames(cleanDat) <- c("control_Evs_1","control_Evs_2","control_Evs_3","TGFB_Evs_1","TGFB_Evs_2","TGFB_Evs_3")
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


############Not needed as gene ids were provided for the raw counts##################
# ########Symbol conversion for mouse IDs#######################################################################################
# library(biomaRt)
# 
# #mart=useMart('ensembl')
# 
# #set up connection to ensembl mart for mouse
# mouse<- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org/")
# tableSymbols<-getBM(attributes=c("mgi_symbol","ensembl_gene_id"), filters="ensembl_gene_id", values=rownames(dds),mart=mouse) 
# 
# #sanity checks
# dim(dds) #ControlVsDay28
# 
# dim(tableSymbols)
# #put symbols in countData order
# 
# symbols<-tableSymbols$mgi_symbol[match(rownames(dds),tableSymbols$ensembl_gene_id)] #ControlVsDay28
# length(symbols)
# tail(symbols,100)
# 
# #replace NA with empty string where no lookup was possible.
# symbols[is.na(symbols)]<-""
# tail(symbols,100)
# 
# uniqueIDs<-paste0(symbols,"|",rownames(dds)) #ControlVsDay28
# tail(uniqueIDs)
#rownames(dds)<-uniqueIDs  #ControlVsDay28
####################################################################################################################
dds
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
res05 <- results(dds, alpha = 0.05) 
res05 #table w statistics

write.csv(res05, "res05.csv")

#diff ex analysis
summary(res05) #Summary deseq2 results 
dim(res05)

# contrasts
resultsNames(dds)


#up and down list created to pull out diff ex genes
up.1 <- res05[which(res05$log2FoldChange > 0 & res05$padj < .05),]
down.1 <- res05[which(res05$log2FoldChange < 0 & res05$padj < .05),]

write.csv(up.1, file = "Up_dif_ex.csv")
write.csv(down.1, file = "Down_dif_ex.csv")

#########################################################################CHANGE#################################################################################################

#diff p adj value 
res10 <- results(dds, alpha = 0.10) 
res10 #table w statistics

#diff ex analysis
summary(res10) #Summary deseq2 results 
dim(res10)

# contrasts
resultsNames(dds)

#up and down list created to pull out diff ex genes
up.10 <- res10[which(res10$log2FoldChange > 0 & res10$padj < .10),]
down.10 <- res10[which(res10$log2FoldChange < 0 & res10$padj < .10),]

write.csv(up.10, file = "GO_Up_dif_ex.csv")
write.csv(down.10, file = "GO_Down_dif_ex.csv")

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
                title = 'Treatment_TGFB_vs_Control',
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
pheatmap(matrixx, annotation_col=annotation_data, show_rownames = T,cluster_cols= FALSE, main = "Treatment_TGFB_vs_Control") #true / False to toggle supervised and unsupervised clustering
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



#Heatmap Explanation:
# This code performs several steps to generate a heatmap visualization of a gene expression matrix:
# The first line creates a data frame called "AdjP" by converting the results of a statistical analysis ("res05") into a data frame.
# The second line replaces any missing values in the "padj" column of "AdjP" with 1.
# The third line sorts the rows of "AdjP" by the "padj" column in ascending order and stores the resulting row order in a variable called "var".
# The fourth line extracts the top 50 rows of a gene expression matrix ("rld") based on the row order in "var" and stores the result in a variable called "matrix".
# The fifth line subtracts the row means of "matrix" from each row.
# The sixth and seventh lines set a threshold value for the absolute maximum and minimum values of "matrix" and clip any values outside this range to the threshold value.
# The eighth line creates a data frame called "annotation_data" by extracting a single column ("Treatment") from the column data of the "rld" object.
# The final line generates a heatmap visualization of "matrixx" using the "pheatmap" function. The "annotation_col" argument specifies the column data to use for color-coding the heatmap, "show_rownames" specifies whether to display row names, "cluster_cols" specifies whether to cluster columns (genes), and "main" specifies the main title of the plot.
# Overall, this code performs data preprocessing, selects a subset of genes based on statistical significance, and generates a heatmap visualization to help visualize gene expression patterns across samples.




####Rld vst and blind = false argument
##when performing clustering we must remove the dependence of the variance on the mean. 
#In other words we must remove the experiment wide trend in the data before clustering. 
#There are two functions within DEseq2 to transform the data in such a manner, the first is to use a regularized logarithm rlog() and the second is the variance stablizing transform vst(). 
#There are pros and cons to each method, we will use vst() here simply because it is much faster. 
#By default both rlog() and vst() are blind to the sample design formula given to DEseq2 in DESeqDataSetFromMatrix(). 
#However this is not appropriate if one expects large differences in counts, which can be explained by the differences in the experimental design. 
#In such cases the blind parameter should be set to FALSE.






###################################
#Heatmap theory: 
# 1. https://www.biostars.org/p/273199/ 
#2. https://genviz.org/module-04-expression/0004/02/01/DifferentialExpression/ 
#3. https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html#:~:text=NOTE%3A%20DESeq2%20doesn't%20actually,Generalized%20Linear%20Model%20(GLM).
