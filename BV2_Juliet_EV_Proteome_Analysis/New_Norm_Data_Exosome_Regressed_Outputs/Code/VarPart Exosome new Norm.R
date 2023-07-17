
#####cell lysates from my BV2 polarization experiment Data Analysis - Aditya Natu - Rangaraju Lab

###Data: mass spec data on the cell lysates from my BV2 polarization experiment (with gene names on the far right column). From this I would like to see if the treatments I added (LPS, TGFB, IL-10) did in fact polarize my BV2 cells.

###Required: 

cat("\014")
rm(list = ls())
dev.off()

##################Setting the work space##################
library(M3C)
setwd("/Users/adityanatu/Downloads/PROJECTS/BV2_Juliet_Exosome_Analysis/New_Norm_Data_Exosome/")
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

#Dat dim 533 to 531: Duplicates and missing values are stored in Duplicated.Genename.CleanDat
dat.cleaned <- read.csv("Exosome_Cleaned_Gene_Names_BV2_cells.csv", header = TRUE, row.names = 17)
#[1] 535  16
dim(dat.cleaned)

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
colnames(dat.cleaned) <- c("Norm_BV2_control_1_cell", 
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
cleanDat <- dat.cleaned
numericMeta <- traits 


table(rownames(numericMeta) == colnames(cleanDat))





############################################
#Run VP to check how much batch in data 
##########################################
library(variancePartition)
regvars.vp<-data.frame(numericMeta)

regvars.vp$Treatment<-factor(regvars.vp$Treatment)
regvars.vp$Batch<-factor(regvars.vp$Batch)


#regvars.vp$Age<-as.numeric(regvars.vp$Age)
#regvars.vp$Group<-factor(regvars.vp$Group)
#regvars.vp$NIHSS<-as.numeric(regvars.vp$NIHSS)
#regvars.vp$GCS<-as.numeric(regvars.vp$GCS)

form <- ~ (1|Treatment) + (1|Batch) #Age + NIHSS + GCS + (1|Sex)

#Variance Partition to check for Batch 
#Status: Post-correction (or regression)

cleanDat.noNA<-impute::impute.knn(as.matrix(cleanDat))$data    

#regvars.vp.nc <- regvars.vp[match(colnames(cleanDat.nc),rownames(regvars.vp)),]
varPart2 <- fitExtractVarPartModel(cleanDat.noNA, form, regvars.vp) 

vp2 <- sortCols(varPart2,FUN=max,last= c("Treatment","Batch","Residuals"))

pdf(file="atop50.pdf")
par(mfrow=c(1,1))

plotVarPart( vp2, main="New Norm Exosome data VarPart" )

BatchSortOrder<-order(vp2[["Treatment"]],decreasing=TRUE)  
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][BatchSortOrder]; }
rownames(vp2)<-rownames(vp2)[BatchSortOrder]
plotPercentBars( vp2[1:50,]) + ggtitle( "Top Treatment-covariate Gene Products" ) 

BatchSortOrder<-order(vp2[["Batch"]],decreasing=TRUE)  
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][BatchSortOrder]; }
rownames(vp2)<-rownames(vp2)[BatchSortOrder]
plotPercentBars( vp2[1:50,]) + ggtitle( "Top Batch-covariate Gene Products" ) 

dev.off()





####################################################
######################################
#Setting up var!
cleanDat.original <- cleanDat
######################################
#Parallel Bootstrap Regression

cleanDat.unreg<-cleanDat.original
cleanDat <- cleanDat.unreg
dim(numericMeta)


library("doParallel")
parallelThreads=10 
#clusterLocal <- makeCluster(c(rep("haplotein.biochem.emory.edu",parallelThreads)), type = "SOCK", port=10191, user="edammer", rscript="/usr/bin/Rscript",rscript_args="OUT=/dev/null SNOWLIB=/usr/lib64/R/library",manual=FALSE)
##OR to run parallel processing threads for regression: :
#parallelThreads=8 #max is number of processes that can run on your computer at one time
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")

registerDoParallel(clusterLocal)

 ##OR for no parallel processing skip all doParallel functions (much slower):
#parallelThreads=1


library(boot)
boot <- TRUE #keep true 
numboot <- 1000
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows bootstrap function to select samples
  fit <- lm(formula, data=d)
  return(coef(fit))
}  

#condition <- as.factor(numericMeta$Diagnosis.Simple)  #$Group

#Define your variables: Variable formats needs to be specified 
#Contineous variable is defined as numeric 
#categorical variable is defined as a factor 


treatment=as.factor(numericMeta$Treatment)
batch<-as.factor(numericMeta$Batch)  

#batch<-as.factor(numericMeta$Batch)
batch<-relevel(batch,"1")  #Change reference level to a particular batch


#b15<-b14<-b13<-b12<-b11<-b10<-b09<-b08<-b07<-
#b06<-b05<-

b04<-b03<-b02<-b01<-rep(0,nrow(numericMeta))
b01[numericMeta$Batch=="1"]<- 1
b02[numericMeta$Batch=="2"]<- 1
b03[numericMeta$Batch=="3"]<- 1
b04[numericMeta$Batch=="4"]<- 1

#b05[numericMeta$Batch=="b05"]<- 1
#b06[numericMeta$Batch=="b06"]<- 1
#b07[numericMeta$Batch=="b07"]<- 1
#b08[numericMeta$Batch=="b08"]<- 1
#b09[numericMeta$Batch=="b09"]<- 1
#b10[numericMeta$Batch=="b10"]<- 1
#b11[numericMeta$Batch=="b11"]<- 1
#b12[numericMeta$Batch=="b12"]<- 1
#b13[numericMeta$Batch=="b13"]<- 1
#b14[numericMeta$Batch=="b14"]<- 1
#b15[numericMeta$Batch=="b15"]<- 1

#check here and update var!
regvars <- as.data.frame(cbind(treatment)) #condition,... #group in somalogic instance but global var doesnt change!
#regvars <- data.frame(batch=batch)


## Run the regression
normExpr.reg <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(cleanDat))
rownames(normExpr.reg) <- rownames(cleanDat)
colnames(normExpr.reg) <- colnames(cleanDat)

# + 2: 
coefmat <- matrix(NA,nrow=nrow(cleanDat),ncol=ncol(regvars) + 2) ##length(levels(batch))-1 (removed - includes the variable space for different batches)

## change this to ncol(regvars)+2 when condition has 2 levels if BOOT=TRUE, +1 if BOOT=FALSE
#Coefmat column 2, first level of batch, will be all NA because it is the reference level of this factored variable. Column 1 is from thisexp (cleanDat).
#another RNG seed set for reproducibility

set.seed(8675309);

if (parallelThreads > 1) {
  
  if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
    set.seed(8675309)
    cat('[bootstrap-PARALLEL] Working on ORDINARY NONPARAMETRIC BOOTSTRAP regression with ', parallelThreads, ' threads over ', nrow(cleanDat), ' iterations.\n Estimated time to complete:', round(120/parallelThreads*nrow(cleanDat)/2736,1), ' minutes.\n') #intermediate progress printouts would not be visible in parallel mode
    coefmat <- foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% {
      set.seed(8675309)
      options(stringsAsFactors=FALSE)
      library(boot)
      thisexp <- as.numeric(cleanDat[i,])
      #variables here are all that you need to include!
      df.reg<-model.matrix(~thisexp +treatment+batch, data=data.frame(thisexp,regvars)) 
      
      bs.results <- boot(data=as.data.frame(df.reg),statistic=bs,
                         R=numboot, formula=thisexp ~ .)  ## run 1000 resamplings
      
      apply(bs.results$t,2,function(x) median(x,na.rm=TRUE))  #this is bs.stats (output to coefmat 1 row at a time)
    }
    #batchwise-missing coefficients should be 0, not NA so that NAs do not propagate on the below line
    coefmat.noNA<-coefmat
    coefmat.noNA[is.na(coefmat)]<-0
    normExpr.reg <-foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% { (cleanDat[i,]  - coefmat.noNA[i,4]*b02 - coefmat.noNA[i,5]*b03 -coefmat.noNA[i,6]*b04) } # - coefmat.noNA[i,6]*b02 - coefmat.noNA[i,7]*b03 - coefmat.noNA[i,8]*b04 - coefmat.noNA[i,9]*b05 - coefmat.noNA[i,10]*b06 (removed as no batch)
    normExpr.regBatchOnly <- foreach (i=1:nrow(cleanDat), .combine=rbind) %dopar% { (cleanDat[i,] - coefmat.noNA[i,4]*b02 - coefmat.noNA[i,5]*b03 -coefmat.noNA[i,6]*b04) } #- coefmat.noNA[i,6]*b02 - coefmat.noNA[i,7]*b03 - coefmat.noNA[i,8]*b04 - coefmat.noNA[i,9]*b05 - coefmat.noNA[i,10]*b06
    #make sure your age and sex are numeric - else the multiplication will break!
    #Else <-  Boot = False
  } else { #linear model regression; faster but incomplete regression of Age, Sex, PMI effects, SO NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
    #NOTE: batch/factored not implemented for lm option BOOT=FALSE.
    coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(cleanDat[i,])~batch,data=regvars) #+PMI
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      normExpr.reg[i,] <- coef[1] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      #cat('Done for Protein ',i,'\n')
    }
  } #end parallel option -- Average run time estimate printed in console in minutes is calculated based on benchmark of a 2+ GHz intel Xeon with 30 threads 
} else {
  if (boot==TRUE) { #ORDINARY NONPARAMETRIC BOOTSTRAP
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      thisexp <- as.numeric(cleanDat[i,])
      bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                         R=numboot, formula=thisexp~ batch)  #+PMI) ## run 1000 resamplings
      #                       R=numboot, formula=thisexp~ condition +age+sex+PMI) ## run 1000 resamplings
      ## get the median - we can sometimes get NA values here... so let's exclude these - old code #bs.stats <- apply(bs.results$t,2,median) 
      bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
      for (n in 1:ncol(bs.results$t)) {
        bs.stats[n] <- median(na.omit(bs.results$t[,n]))
      }
      coefmat[i,] <- bs.stats
      normExpr.reg[i,] <- thisexp - bs.stats[2]*regvars[,"batch"]  #- bs.stats[4]*regvars[,"sex"]  #- bs.stats[5]*regvars[,"PMI"]
      cat('[bootstrap] Done for Protein ',i,'\n')
    }
  } else { #linear model regression; faster but NOT USED WITH boot=TRUE (requires changing coefmat matrix ncol to 1 less above)
    coefmat<-coefmat[,-ncol(coefmat)] #handles different column requirement for lm regression method
    for (i in 1:nrow(cleanDat)) {
      if (i%%1000 == 0) {print(i)}
      lmmod1 <- lm(as.numeric(cleanDat[i,])~batch,data=regvars)  #+PMI
      ##datpred <- predict(object=lmmod1,newdata=regvars)
      coef <- coef(lmmod1)
      coefmat[i,] <- coef
      normExpr.reg[i,] <- coef[1] + lmmod1$residuals ## The full data - the undesired covariates
      ## Also equivalent to <- thisexp - coef*var expression above
      cat('Done for Protein ',i,'\n')
    }
  } #end nonparallel option
}


## Sanity Check -- Did regression do something unexpected to our abundance data?
quantile(cleanDat.unreg[,ncol(cleanDat)],c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)
##paste below as comment the abundance data spread (quantiles) before regression for comparison of data spread to after-regression abundance data
#         0%        2.5%         25%         50%         75%       97.5%        100% 
#-3.85209815 -0.56198482 -0.09152705  0.00000000  0.08232791  0.35085274  2.98840540

quantile(normExpr.reg[,ncol(cleanDat)],c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)
##paste below as comment the abundance data spread (quantiles) after regression -- they often spread a bit but should not go out of a reasonable range
#          0%         2.5%          25%          50%          75%        97.5%         100% 
#-7.637627909 -0.801092300 -0.142569142 -0.008009563  0.127239972  0.545465455  3.903886569

quantile(normExpr.regBatchOnly[,ncol(cleanDat)],c(0,0.025,0.25,0.5,0.75,0.975,1),na.rm=TRUE)
##paste below as comment the abundance data spread (quantiles) after regression -- they often spread a bit but should not go out of a reasonable range
#         0%        2.5%         25%         50%         75%       97.5%        100% 
#-7.94097066 -0.86626968 -0.16654333 -0.01295554  0.13272330  0.62843667  3.87143607


##Overwrite cleanDat with regressed data
##(DO NOT RERUN OUT OF CONTEXT)
############################################
cleanDat<-normExpr.reg

cleanDat.regBatchOnly<-normExpr.regBatchOnly

rownames(cleanDat)<-rownames(cleanDat.unreg)

rownames(cleanDat.regBatchOnly)<-rownames(cleanDat.unreg)
############################################


############################################
#Run VP to check how much batch is regressed in data 
##########################################

library(variancePartition)
regvars.vp<-data.frame(numericMeta)
regvars.vp$Treatment<-factor(regvars.vp$Treatment)
regvars.vp$Batch<-factor(regvars.vp$Batch)


#regvars.vp$Age<-as.numeric(regvars.vp$Age)
#regvars.vp$Group<-factor(regvars.vp$Group)
#regvars.vp$NIHSS<-as.numeric(regvars.vp$NIHSS)
#regvars.vp$GCS<-as.numeric(regvars.vp$GCS)

form <- ~ (1|Treatment) + (1|Batch) #Age + NIHSS + GCS + (1|Sex)

#Variance Partition to check for Batch 
#Status: Post-correction (or regression)

cleanDat.noNA<-impute::impute.knn(as.matrix(cleanDat))$data    #removes NA's

#regvars.vp.nc <- regvars.vp[match(colnames(cleanDat.nc),rownames(regvars.vp)),]
varPart2 <- fitExtractVarPartModel(cleanDat.noNA, form, regvars.vp) 

vp2 <- sortCols(varPart2,FUN=max,last= c("Treatment","Batch","Residuals"))

pdf(file="atop50 Regressed.pdf")
par(mfrow=c(1,1))

plotVarPart( vp2, main="Regressed New Norm Exosome data VarPart" )

BatchSortOrder<-order(vp2[["Treatment"]],decreasing=TRUE)  
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][BatchSortOrder]; }
rownames(vp2)<-rownames(vp2)[BatchSortOrder]
plotPercentBars( vp2[1:50,]) + ggtitle( "Reg Top Treatment-covariate Gene Products" ) 

BatchSortOrder<-order(vp2[["Batch"]],decreasing=TRUE)  
for (i in ls(vp2)) { vp2[[i]]<-vp2[[i]][BatchSortOrder]; }
rownames(vp2)<-rownames(vp2)[BatchSortOrder]
plotPercentBars( vp2[1:50,]) + ggtitle( "Reg Top Batch-covariate Gene Products" ) 

dev.off()


write.csv(cleanDat, file = "regressed_newNorm_Juliet.csv")



