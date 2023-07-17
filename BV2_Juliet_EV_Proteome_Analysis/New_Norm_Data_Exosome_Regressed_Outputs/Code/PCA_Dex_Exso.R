

#### PCA With DeX proteins: coming from toplot!


#Sanity check and just doing a rowmatch leads to 128 unique across all comparisons!
###Creating new cleanDat 
#Kmeans_Expr_Matrix <- backg_nano_kmeans[backg_nano_kmeans$X %in% match_diff_ex_kmeans$V1,]

write.csv(cleanDat, file = "ALL_cleanDat_used4PCA.csv")

write.csv(toplot, file = "toplot_thePCA_prots.csv")

cleanDat.ALL <- read.csv("ALL_cleanDat_used4PCA.csv", header = TRUE)
toplot.PCA <- read.csv("toplot_thePCA_prots.csv", header = TRUE)

cleanDat.PCA_Anova <- cleanDat.ALL[cleanDat.ALL$X %in% toplot.PCA$X,]


write.csv(cleanDat.PCA_Anova, file = "PCA_prots_final.csv")

cleanDat.PCA_Anova <- read.csv("PCA_prots_final.csv", header = TRUE, row.names = 1)

#PCA of Diff Ex prots
cleanDat.original_test <- as.data.frame(cleanDat.PCA_Anova)
cleanDat.original_test <- as.matrix(cleanDat.PCA_Anova)
typeof(cleanDat.original_test)
cleanDat.original_test <- apply(cleanDat.original_test, 2, as.numeric)
#Sample PCA SomaData
pca(cleanDat.original_test, colvec = c('skyblue', 'red', 'yellow', 'green'), labels = numericMeta$Treatment, controlscale=TRUE,scale=3,printres=TRUE,printwidth=25, legendtitle = "Dex_Treatment")

