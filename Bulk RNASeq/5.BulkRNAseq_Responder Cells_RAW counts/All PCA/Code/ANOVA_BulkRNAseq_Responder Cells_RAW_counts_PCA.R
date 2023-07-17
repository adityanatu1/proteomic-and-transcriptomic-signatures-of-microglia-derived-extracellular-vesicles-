
##############CHECK##############
#ALL PCA!
rm(list = ls())
dev.off()
cat("\014")


##############CHECK##############
setwd("/Users/adityanatu/Library/Mobile Documents/com~apple~CloudDocs/Project_2023/Juliet_Data_Analysis/EV_RNASeq/5.BulkRNAseq_Responder Cells_RAW counts/All PCA/")
getwd()
options(stringsAsFactors = FALSE)

library(DESeq2)
library(tidyverse)
library(biomaRt)
library(pheatmap)

#########################################################CHANGE##############################################################################
#Subsetting the cleanDat and traits file  was done in excel!
##############CHECK##############
cleanDat <- read.csv("BulkRNAseq_Responder Cells_RAW_counts.csv", header = TRUE, row.names = 1)

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


numericMeta <- read.csv("BulkRNAseq_Responder Cells_RAW_traits.csv", header = TRUE, row.names = 1)

colnames(cleanDat)
rownames(numericMeta)

##############CHECK##############
colnames(cleanDat) <- c("control_Evs_Tx_1","control_Evs_Tx_2","control_Evs_Tx_3" ,"LPS_Evs_Tx_1", "LPS_Evs_Tx_2", "LPS_Evs_Tx_3", "TGFB_Evs_Tx_1"   
                        ,"TGFB_Evs_Tx_2"  ,  "TGFB_Evs_Tx_3"  ,  "IL10_Evs_Tx_1"  ,  "IL10_Evs_Tx_2"   , "IL10_Evs_Tx_3")
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


###########Not needed as gene ids were provided for the raw counts##################
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
####################################################################################################################
dds.LTR <- dds


#Overall Model Values: 
#Run DESeq ----------------------
#ANOVA
#This p-value indicates whether the inclusion of the treatment variable significantly improves the fit of the model 
#compared to an intercept-only model. Therefore, it is not a p-value for all the comparisons, 
#but rather for the overall significance of the treatment variable across all the comparisons.


#Here, we are using the "LRT" (Likelihood Ratio Test) option in the DESeq function to perform ANOVA analysis. The "reduced" formula specifies a null model with no experimental variables. The full model includes all the variables in the design matrix.
dds.LTR <- DESeq(dds.LTR, test="LRT", reduced=~1)

#dds <- DESeq(dds)

#diff p adj value 
#Here, we are extracting the ANOVA results using the "results" function and ordering them by the adjusted p-value. The "padj" column contains the adjusted p-value, which accounts for multiple testing.

res05 <- results(dds.LTR, alpha = 0.05) 
res05 #table w statistics

write.csv(res05, "Overall_res05.csv")

#diff ex analysis
#These are the basic steps for performing ANOVA analysis using DESeq2 in R. Note that there are other options for performing ANOVA analysis, such as using the "glm" or "nbinomLRT" test options in the DESeq function. Additionally, it is important to carefully consider the experimental design and other factors that may influence the results of the analysis.
summary(res05) #Summary deseq2 results 
dim(res05)
#############Individual Model

dds <- DESeq(dds)

results_TGFBvsControl <- results(dds, contrast=c("Treatment", "TGFB", "Control"), alpha = 0.05)
results_TGFBvsControl

write.csv(results_TGFBvsControl, "TFGBres05.csv")

#diff ex analysis
summary(results_TGFBvsControl) #Summary deseq2 results 
dim(results_TGFBvsControl)

##############################################################################PCA##############################################################################
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup = c("Treatment"))
##############################################################################PCA##############################################################################
#FAQ Section for Code & Analysis: 
#LRT p-value: '~ Treatment' vs '~ 1'  : What does this mean?
#In statistics, the likelihood ratio test (LRT) is a statistical test used to compare two nested models, where one model is a simpler version of the other. The LRT compares the likelihood of the data under the simpler model (the null hypothesis) to the likelihood of the data under the more complex model (the alternative hypothesis).
  # 
  # In the context of '~ Treatment' vs '~ 1', this means we are comparing two linear regression models:
  # The null model, which has an intercept term only (~ 1).
  # The alternative model, which includes a treatment variable in addition to the intercept term (~ Treatment).
  # The LRT tests whether the treatment variable significantly improves the fit of the model compared to the null model. In other words, it tests whether the treatment variable explains a significant amount of the variation in the response variable.
  # The p-value of the LRT tells us the probability of obtaining a test statistic as extreme as, or more extreme than, the observed test statistic if the null hypothesis were true (i.e., if the treatment variable did not significantly improve the model fit). A small p-value (typically less than 0.05) indicates strong evidence against the null hypothesis, and suggests that the treatment variable does significantly improve the model fit.
  # So, in summary, the LRT p-value for '~ Treatment' vs '~ 1' tests whether including a treatment variable significantly improves the fit of the linear regression model compared to an intercept-only model. A small p-value indicates strong evidence that the treatment variable is a significant predictor of the response variable.

##whether including a treatment variable significantly improves the fit of the linear regression model compared to an intercept-only model - explain this
  # In linear regression, the goal is to model the relationship between a response variable (also called the dependent variable) and one or more predictor variables (also called independent variables). The linear regression model assumes that the relationship between the response variable and the predictor variables can be described by a straight line.
  # The simplest linear regression model is an intercept-only model, which assumes that the response variable is only related to a constant intercept term. In other words, the intercept-only model assumes that the response variable is constant for all levels of the predictor variable(s), and does not depend on any specific values of the predictor variable(s).
  # However, in many cases, the response variable is influenced by one or more predictor variables. For example, in a medical study, the response variable might be a patient's blood pressure, and the predictor variable might be a treatment that the patient receives (e.g., a drug or a placebo). In this case, the treatment variable is likely to have an effect on the patient's blood pressure, and including the treatment variable in the regression model can improve the model's ability to explain the variation in the response variable.
  # The likelihood ratio test (LRT) is a statistical test that compares the fit of two nested models: a simpler model (the null model) and a more complex model (the alternative model). In the context of '~ Treatment' vs '~ 1', the null model is an intercept-only model, and the alternative model includes a treatment variable in addition to the intercept term.
  # The LRT tests whether including the treatment variable in the model significantly improves the model fit compared to the intercept-only model. If the p-value of the LRT is less than a pre-specified significance level (e.g., 0.05), then we reject the null hypothesis (that the intercept-only model is a better fit for the data) and conclude that the treatment variable significantly improves the model fit.
  # Therefore, when we say that "including a treatment variable significantly improves the fit of the linear regression model compared to an intercept-only model", we mean that the treatment variable is a significant predictor of the response variable, and that including it in the model improves the model's ability to explain the variation in the response variable.


##resultsNames(dds) explain this
  #In R, the resultsNames() function is used in conjunction with the DESeq2 package to obtain the names of the differentially expressed gene (DEG) statistics generated by the differential gene expression analysis.
  #After running the differential expression analysis using DESeq2, the resulting dds object contains a variety of statistical measures for each gene, including fold change, p-values, and adjusted p-values. resultsNames() is a function that returns a character vector of the names of these statistics.
  #The resultsNames() function is useful because it allows users to easily extract specific results from the dds object. For example, if a user wants to extract the p-values for all DEGs in their analysis, they can use resultsNames() to obtain the name of the p-value statistic and then use that name to extract the p-values for all genes in the dds object.
  #Overall, resultsNames() is a convenient function for working with the output of the differential gene expression analysis and for quickly extracting specific results of interest.


#[1] "Intercept"                 "Treatment_IL10_vs_Control" "Treatment_LPS_vs_Control"  "Treatment_TGFB_vs_Control" Explain this output from results names after LTR test
  #The output you provided from resultsNames() after performing an LRT test in DESeq2 indicates the names of the statistical tests that were performed in the differential gene expression analysis.
  #"Intercept" is the name of the statistic for the base model or the intercept-only model. This model assumes that there are no differences in gene expression between the groups being compared, and it is used as a reference point to compare with more complex models that include additional variables.
  #"Treatment_IL10_vs_Control", "Treatment_LPS_vs_Control", and "Treatment_TGFB_vs_Control" are the names of the statistics for the treatment variables being compared with the control group. These variables represent the experimental groups being compared in the analysis, where IL10, LPS, and TGFB are treatments being compared with a control group.
  #Each of these statistics contains information about differential expression between the experimental groups and the control group for each gene, such as fold change, p-values, and adjusted p-values. These statistics can be used to identify differentially expressed genes and to perform downstream analyses.

# 
# Do the P value that I get after LTR test - is that the p value for all the comparisons or just 1 specific comparison above
  # The p-value that you obtain after performing a LRT (Likelihood Ratio Test) in DESeq2 is a single p-value that represents the significance of the comparison between two models: the full model (including the treatment variable) and the reduced model (without the treatment variable).
  # This p-value indicates whether the inclusion of the treatment variable significantly improves the fit of the model compared to an intercept-only model. Therefore, it is not a p-value for all the comparisons, but rather for the overall significance of the treatment variable across all the comparisons.
  # After obtaining a significant p-value for the LRT, you can then use other functions in DESeq2, such as results(), to obtain the results of individual comparisons (e.g., Treatment_IL10_vs_Control, Treatment_LPS_vs_Control, and Treatment_TGFB_vs_Control) and to obtain p-values and adjusted p-values for each of these comparisons.
  # log2 fold change (MLE): Treatment TGFB vs Control - explain this

# The log2 fold change (MLE) for Treatment TGFB vs Control is a statistic that is generated during differential expression analysis using DESeq2.
 #The log2 fold change (LFC) represents the difference in gene expression levels between two groups, in this case, Treatment TGFB and Control. Specifically, it represents the log2-transformed ratio of the mean expression level in Treatment TGFB versus Control.
  #The MLE (maximum likelihood estimate) is the estimate of the LFC that maximizes the likelihood function of the model, given the observed data. This estimate is based on the counts of RNA-seq reads that are mapped to each gene, as well as on the normalization and dispersion estimates that are calculated by DESeq2.
  # A positive LFC indicates that the gene is upregulated in the Treatment TGFB group compared to the Control group, while a negative LFC indicates downregulation. The absolute value of the LFC represents the magnitude of the difference in expression level between the two groups, with larger absolute values indicating greater differences.
  # The log2 transformation is used to make the LFCs more interpretable and easier to compare between genes. The log2 scale is commonly used in genomics because it provides a linear relationship between fold changes and their corresponding log2-transformed values, and because it can be interpreted as a fold change of a binary variable that is 2-fold different.

# #how to use results function in deseq2 to obtain the results of individual comparisons
  #You can use the results() function in DESeq2 to obtain the results of individual comparisons between treatment groups. The results() function takes two arguments: the DESeqDataSet object (dds) and a character vector indicating the name of the column containing the comparison of interest.
  # Here's an example of how to use the results() function to obtain the results of the comparison between Treatment TGFB and Control:
  # In this example, contrast=c("Treatment", "TGFB", "Control") specifies that we want to compare the expression levels of the "TGFB" treatment group to the "Control" group, with "Treatment" being the column containing the treatment information in the DESeqDataSet.
  # The results() function returns a DataFrame object with the results of the comparison, including the log2 fold change (LFC), p-value, and adjusted p-value (i.e. false discovery rate, FDR) for each gene. By default, the results are sorted by adjusted p-value, with the most significant genes listed first.
  # You can repeat this process with different comparison contrasts to obtain the results for different pairwise comparisons between treatment groups.

# Explain how to do an LTR Test in DEseq2
  # The likelihood ratio test (LRT) is a statistical test that can be used in DESeq2 to compare the fit of two nested models. The purpose of the test is to determine whether the inclusion of a particular variable (e.g., a treatment variable) significantly improves the fit of the model compared to a simpler model without that variable.
  # Here are the general steps to perform an LRT in DESeq2:
  # First, you need to have already loaded and preprocessed your count data using DESeq2. This involves steps such as normalization, filtering, and estimating size factors and dispersion values.
  # Once you have your DESeqDataSet object, you can use the DESeq() function to estimate the size factors and dispersions and perform the differential expression analysis.
  # After running DESeq(), you can use the results() function to obtain the log2 fold change and p-values for each gene or transcript.
  # To perform the LRT, you need to create two models: one with the variable you are interested in (e.g., treatment) and one without it. You can do this using the formula interface, where the '~' separates the formula into the dependent variable (on the left) and the independent variable(s) (on the right).
  # For example, to create a model with treatment, you would use:blah blah blah as shown above!
# 
# Explain the following terms: 
  # 1. compare the fit of two nested models
  # 2. simpler model without that variable

  # Compare the fit of two nested models: This refers to comparing the performance of two statistical models, one of which is a more complex version of the other. 
  #The more complex model contains additional variables or parameters compared to the simpler model. Comparing the fit of these models is done to determine if the additional variables or parameters in the more complex model provide a statistically significant improvement in the model's ability to explain the data. The process of comparing the models involves performing a likelihood ratio test (LRT), which is a statistical test that compares the likelihood of the data given the two models and calculates a p-value indicating whether the difference in fit between the models is statistically significant.
  
  # Simpler model without that variable: This refers to a statistical model that is less complex than another model it is being compared to, by excluding one or more variables or parameters. 
  #The simpler model is referred to as a "reduced model" or "null model" and is often used as a baseline to compare to a more complex model that includes additional variables or parameters. 
  #The simpler model without a particular variable is obtained by fitting the model with all variables except the one being excluded from the analysis. This is done to determine whether the excluded variable is statistically significant in explaining the variability in the response variable, and whether its inclusion in the more complex model is warranted.

# dds <- DESeq(dds, test="LRT", reduced=~1) Explain the reduced parameter
  # In DESeq2, the reduced parameter is used in the likelihood ratio test (LRT) to specify the reduced model. The reduced model is a simpler model with fewer explanatory variables than the full model. The LRT compares the fit of the full model with the fit of the reduced model to determine if including the additional explanatory variable(s) in the full model significantly improves the fit.
  # In the example given, the reduced parameter is set to ~1, which means that the reduced model includes only an intercept term. This is equivalent to fitting a model with no explanatory variables, and the LRT will test whether any of the explanatory variables in the full model are significant in predicting the response variable.
  # The null hypothesis of the LRT is that the reduced model is a better fit to the data than the full model, while the alternative hypothesis is that the full model provides a significantly better fit. By default, DESeq2 uses the chi-squared distribution to calculate the p-values for the LRT.