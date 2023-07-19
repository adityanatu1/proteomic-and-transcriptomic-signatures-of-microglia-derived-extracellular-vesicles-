#####cell lysates from my BV2 polarization experiment Data Analysis - Aditya Natu - Rangaraju Lab
 

##################Setting the work space##################
cat("\014")
rm(list = ls())
dev.off()

library(M3C)
setwd("/Set/Your/Working/Dir")
getwd()
options(stringsAsFactors = FALSE)

##################Loading | Inspecting | cleaning the Data##################
#Loading and Cleaning the Input Data
dat.original <- read.csv("BV2 exosome perseus output normalized_9.15.22.csv", header = TRUE)
df <-dat.original

Clean.ids <- df[!grepl("CON_", df$Protein.IDs),]
Contaminated.Genename.CleanDat <- df[grepl("CON_", df$Protein.IDs),] #Duplicates are stored here


#write.csv(Duplicated.Genename.CleanDat, file = "Exosome_Duplicates_Dropped_BV2.csv")

Cleaned.Genename.CleanDat <- Clean.ids #CleanDat without duplicates are stoed here and used moving forward for analysis! 
#write.csv(Cleaned.Genename.CleanDat, file = "Exosome_Cleaned_Gene_Names_BV2_cells.csv") #file needed to be manipulated in excel before it could be used - someprotids were copied into genename to prevent blank gene names and other non imp columns were dropped 

###################START FROM BELOW COMMENT TO USE CLEANED DATA###################
###############################CLEANED DATA################################


#dat.cleaned <- read.csv("Exosome_Cleaned_Gene_Names_BV2_cells.csv", header = TRUE, row.names = 17)
dat.cleaned.reg <- read.csv("regressed_newNorm_Juliet.csv", header = TRUE, row.names = 1)

#[1] 535  16
#dim(dat.cleaned)
dim(dat.cleaned.reg)

#Traits
traits <- read.csv("Sample manifests.csv", header = TRUE, stringsAsFactors = TRUE, row.names = 2)
dim(traits)


rownames(traits) <- c("Norm_BV2_control_1_cell", 
                      "Norm_BV2_TGFB_1_cell",  
                      "Norm_BV2_IL10_1_cell", 
                      "Norm_BV2_LPS_1_cell", 
                      "Norm_BV2_control_2_cell",
                      "Norm_BV2_TGFB_2_cell",  
                      "Norm_BV2_IL10_2_cell",  
                      "Norm_BV2_LPS_2_cell", 
                      "Norm_BV2_control_3_cell",
                      "Norm_BV2_TGFB_3_cell",
                      "Norm_BV2_IL10_3_cell",  
                      "Norm_BV2_LPS_3_cell",
                      "Norm_BV2_control_4_cell",
                      "Norm_BV2_TGFB_4_cell",  
                      "Norm_BV2_IL10_4_cell",
                      "Norm_BV2_LPS_4_cell")


#Cleaning up columns
colnames(dat.cleaned.reg) <- c("Norm_BV2_control_1_cell", 
                           "Norm_BV2_TGFB_1_cell",  
                           "Norm_BV2_IL10_1_cell", 
                           "Norm_BV2_LPS_1_cell", 
                           "Norm_BV2_control_2_cell",
                           "Norm_BV2_TGFB_2_cell",  
                           "Norm_BV2_IL10_2_cell",  
                           "Norm_BV2_LPS_2_cell", 
                           "Norm_BV2_control_3_cell",
                           "Norm_BV2_TGFB_3_cell",
                           "Norm_BV2_IL10_3_cell",  
                           "Norm_BV2_LPS_3_cell",
                           "Norm_BV2_control_4_cell",
                           "Norm_BV2_TGFB_4_cell",  
                           "Norm_BV2_IL10_4_cell",
                           "Norm_BV2_LPS_4_cell")

#Rewriting to conventional variables used in analysis 
cleanDat <- dat.cleaned.reg
numericMeta <- traits 


table(rownames(numericMeta) == colnames(cleanDat))

#TRUE 
#16

#gives you the overall anova P value that you can plug into your heatmap and anovaP PCA! 
Grouping <- numericMeta$Treatment
#Choose columns to get unadjp or adjp volc
source("./parANOVA.dex.R")
ANOVAout <- parANOVA.dex()


#Nterm vs Cterm (n on the left and c on the right)
labelTop <- 25
outFileSuffix="All Comps"
signifP=0.05 
plotVolc()
##################Principle Component Analysis##################

cleanDat.original_test <- as.data.frame(cleanDat)
cleanDat.original_test <- as.matrix(cleanDat)

cleanDat.original_test <- apply(cleanDat.original_test, 2, as.numeric)
#PCA 
pca(cleanDat.original_test, colvec = c('skyblue', 'red', 'black', 'green'), labels = numericMeta$Treatment, controlscale=TRUE,scale=3,printres=TRUE,printwidth=25, legendtitle = "Reg_Original_Treatment")
##################Principle Component Analysis##################

#Control Vs LPS 

#options(ggrepel.max.overlaps = Inf)
cleanDat <- cleanDat[,c(1,5,9,13,4,8,12,16)]
Grouping.numericMeta <- numericMeta[c(1,5,9,13,4,8,12,16),]
Grouping <- Grouping.numericMeta$Treatment
table(rownames(Grouping.numericMeta) == colnames(cleanDat))
#Choose columns to get unadjp or adjp volc
source("./parANOVA.dex.R")
ANOVAout <- parANOVA.dex()
#Nterm vs Cterm (n on the left and c on the right)
labelTop <- 25
outFileSuffix="Control Vs LPS"
signifP=0.05 
plotVolc()


#Control vs TGF 
cleanDat <- cleanDat[,c(1,5,9,13,2,6,10,14)]
Grouping.numericMeta <- numericMeta[c(1,5,9,13,2,6,10,14),]
Grouping <- Grouping.numericMeta$Treatment
table(rownames(Grouping.numericMeta) == colnames(cleanDat))
#Choose columns to get unadjp or adjp volc
source("./parANOVA.dex.R")
ANOVAout <- parANOVA.dex()
#Nterm vs Cterm (n on the left and c on the right)
labelTop <- 25
outFileSuffix="Control vs TGF"
signifP=0.05 
plotVolc()



#Control vs IL10 
cleanDat <- cleanDat[,c(1,5,9,13,3,7,11,15)]
Grouping.numericMeta <- numericMeta[c(1,5,9,13,3,7,11,15),]
Grouping <- Grouping.numericMeta$Treatment
table(rownames(Grouping.numericMeta) == colnames(cleanDat))
#Choose columns to get unadjp or adjp volc
source("./parANOVA.dex.R")
ANOVAout <- parANOVA.dex()
#Nterm vs Cterm (n on the left and c on the right)
labelTop <- 25
outFileSuffix="Control vs IL10"
signifP=0.05 
plotVolc()


##################ANOVA Heatmap adjp < 0.05 ##################
cleanDat.double <- cleanDat
df_Frd_vec_double <- read.csv("Reg_ANOVA_diffEx-ALL-unspecified_study.csv", header = TRUE, row.names = 1)
df_Frd_vec_double <- df_Frd_vec_double[,c(1:2)]
colnames(df_Frd_vec_double) <- c("F-Test", "P_Overall")


##################Heatmap code##################
max_frd = 0.05 

#Genes that are less than ANOVAp overall < 0.05 are used in heatmap and that are significant accross any comparison!
df_Frd_vec_double$SupervisingKeep <- FALSE

table(df_Frd_vec_double$SupervisingKeep)
#535

df_Frd_vec_double$SupervisingKeep[which(df_Frd_vec_double$P_Overall < max_frd) ] <- TRUE
table(df_Frd_vec_double$SupervisingKeep)

toplot <- cleanDat.double[df_Frd_vec_double$SupervisingKeep,]
#toplot.z <- t(apply(toplot, 1,function(x){(x - mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE)}))

reorder <- c(1,5,9,13, 2,6,10,14, 3,7,11,15, 4,8,12,16)
dim(toplot)
metdat <- data.frame(Status=numericMeta$Treatment)


heatmapLegendColors=list('Status'=c("blue","purple","red", "orange"))

library(NMF)
pdf(file = "1test.pdf",width = 10, height = 20)
par(mfrow=c(1,1))
par(mar = c(2, 2, 2, 1));
aheatmap(x=na.omit(toplot[,reorder]), ## Numeric Matrix
         main="gene abound of overallP BV2 Exosome Polarization",
         annCol=metdat[reorder,],
         #annRow=data.frame(Modules=colnames(t(na.omit(t(MEs))))),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=WGCNA::blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA)



## Do not cluster columns - keep given order

aheatmap(x=na.omit(toplot[,reorder]), ## Numeric Matrix
         main="gene abund of overallP BV2 Exosome Polarization Clust Col TRUE",
         annCol=metdat[reorder,],
         #annRow=data.frame(Modules=colnames(t(na.omit(t(MEs))))),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=WGCNA::blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=TRUE) ## Do not cluster columns - keep given order

dev.off()



 

#### PCA With DeX proteins: coming from toplot!
cleanDat.PCA_Anova <- toplot
#PCA of Diff Ex prots
cleanDat.original_test <- as.data.frame(cleanDat.PCA_Anova)
cleanDat.original_test <- apply(cleanDat.original_test, 2, as.numeric)
#Sample PCA 
pca(cleanDat.original_test, colvec = c('skyblue', 'red', 'black', 'green'), labels = numericMeta$Treatment, controlscale=TRUE,scale=3,printres=TRUE,printwidth=25, legendtitle = "Dex_Treatment")















