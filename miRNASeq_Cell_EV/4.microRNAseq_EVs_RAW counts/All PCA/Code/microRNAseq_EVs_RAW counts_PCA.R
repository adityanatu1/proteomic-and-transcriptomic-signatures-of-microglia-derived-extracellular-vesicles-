
##############CHECK##############
#ALL PCA!
rm(list = ls())
dev.off()
cat("\014")


##############CHECK##############
setwd("/Users/adityanatu/Library/Mobile Documents/com~apple~CloudDocs/Project_2023/Juliet_Data_Analysis/EV_RNASeq/4.microRNAseq_EVs_RAW counts/All PCA/")
getwd()
options(stringsAsFactors = FALSE)

library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)

#########################################################CHANGE##############################################################################
#Subsetting the cleanDat and traits file  was done in excel!
##############CHECK##############
cleanDat <- read.csv("microRNAseq_EVs_RAW counts.csv", header = TRUE, row.names = 1)

#not needed for these ids
# #Need Unique rownames! Setting them now!
# 
# dat <-cleanDat
# dim(dat)
# df <-dat
# 
# #Making rownams unique
# rownames(df) <- make.names(dat$Geneid, unique = TRUE)
# write.csv(df, file = "Unique_mRNAseq_EVs_RAW_counts.csv")
# 
# 
# df <- read.csv("Unique_mRNAseq_EVs_RAW_counts.csv", header = TRUE, row.names = 1)[,-c(1)]
# cleanDat <- df
# dim(cleanDat)
# 
# ###Remove var that arent needed!
# rm(dat)
# rm(df)


numericMeta <- read.csv("microRNAseq_EVs_RAW_Traits.csv", header = TRUE, row.names = 1)

colnames(cleanDat)
rownames(numericMeta)

##############CHECK##############
colnames(cleanDat) <- c("control_Evs_1","control_Evs_2","control_Evs_3","LPS_Evs_1" ,"LPS_Evs_2","LPS_Evs_3","IL10_Evs_1",   
                        "IL10_Evs_2","IL10_Evs_3","TGFB_Evs_1","TGFB_Evs_2","TGFB_Evs_3")
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

#12 X 5 = 60 - This case
keep <- rowSums(counts(dds)) >= 60  #check how! 
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


##############################################################################PCA##############################################################################
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup = c("Treatment"))
##############################################################################PCA##############################################################################

