##################Principle Component Analysis##################
setwd("/Users/adityanatu/Library/Mobile Documents/com~apple~CloudDocs/Project_2023/Juliet_Data_Analysis/BV2_Juliet_Cell_Proteome_Analysis/R Code File/")
cleanDat.original_test <- as.data.frame(cleanDat)
cleanDat.original_test <- as.matrix(cleanDat)
typeof(cleanDat.original_test)
cleanDat.original_test <- apply(cleanDat.original_test, 2, as.numeric)
#Sample PCA SomaData
pca(cleanDat.original_test, colvec = c('skyblue', 'red', 'yellow', 'green'), labels = numericMeta$Treatment, controlscale=TRUE,scale=3,printres=TRUE,printwidth=25, legendtitle = "Treatment")



#### PCA With 333 DeX proteins: 



#Sanity check and just doing a rowmatch leads to 128 unique across all comparisons!
###Creating new cleanDat 
#Kmeans_Expr_Matrix <- backg_nano_kmeans[backg_nano_kmeans$X %in% match_diff_ex_kmeans$V1,]

write.csv(cleanDat, file = "cleanDat.333.csv")

write.csv(toplot, file = "toplot.csv")

cleanDat.333 <- read.csv("cleanDat.333.csv", header = TRUE)
toplot.PCA <- read.csv("toplot.csv", header = TRUE)

cleanDat.PCA.333 <- cleanDat.333[cleanDat.333$X %in% toplot.PCA$X,]


write.csv(cleanDat.PCA.333, file = "PCA_updated.csv")

cleanDat.PCA.333 <- read.csv("PCA_updated.csv", header = TRUE, row.names = 1)


cleanDat.original_test <- as.data.frame(cleanDat.PCA.333)
cleanDat.original_test <- as.matrix(cleanDat.PCA.333)
typeof(cleanDat.original_test)
cleanDat.original_test <- apply(cleanDat.original_test, 2, as.numeric)
#Sample PCA SomaData
pca(cleanDat.original_test, colvec = c('skyblue', 'red', 'black', 'green'), labels = numericMeta$Treatment, controlscale=TRUE,scale=3,printres=TRUE,printwidth=25, legendtitle = "333_Treatment")

