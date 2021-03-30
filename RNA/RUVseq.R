# RUV pipeline used in the EE study
# Example script. Can be applied to RNAseq or any other dataset

library(ggbiplot)
library(RUVSeq)
library(devtools)
library(factoextra)
library(dplyr)
library(EDASeq)
library(digest)
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
library("wesanderson")
library(DESeq2)

data_PCA<-as.data.frame(EE1[,-1])

rownames(data_PCA) <- make.names(EE1[,1], unique = TRUE)
length(data_PCA)

log.ir <- log(data_PCA[, 1:6])
filter<- apply(log.ir, 1, function(x) length(x[x>0])>=6)
filtered<-log.ir[filter,]
PCA_data_t<-as.data.frame(t(filtered))
which(apply(PCA_data_t, 2, var)==0)
test<-PCA_data_t[ , apply(PCA_data_t, 2, var) != 0]  # var estimation, if interger (0) no concerns


res.pca <- prcomp(test, scale = TRUE, center=TRUE)
fviz_eig(res.pca)

# fviz_pca_ind(res.pca,
#             col.ind = "cos2", # Colorer par le cos2
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE
#)



res.pca$group[(1:3)] <- "EE"
res.pca$group[(4:6)] <- "CTL"

fviz_pca_ind(res.pca,
            col.ind = res.pca$group, # colorer par groupes
            addEllipses = TRUE, # Ellipse de concentration
            ellipse.type = "confidence",
            legend.title = "Groups",
            repel = TRUE,
            geom="point"
)



# correlation

library(corrplot)
library("PerformanceAnalytics")

flattenCorrMatrix <- function(cormat, pmat) {
 ut <- upper.tri(cormat)
 data.frame(
   row = rownames(cormat)[row(cormat)[ut]],
   column = rownames(cormat)[col(cormat)[ut]],
   cor  =(cormat)[ut],
   p = pmat[ut]
   )
}

library(Hmisc)
res2<-rcorr(as.matrix(data_PCA[,1:9]))
flattenCorrMatrix(res2$r, res2$P)

res <- cor(data_PCA)
round(res, 2)
corrplot(res, type = "upper", order = "hclust",
        tl.col = "black", tl.srt = 45)
# pal <- wes_palette("Zissou1", length(mat_breaks) - 1, type = "continuous")
corrplot(res, method = "square", type = "upper", tl.col = "black", order = "hclust", addrect = 1 , rect.col= "white", col = wes_palette("Darjeeling1"))
# chart.Correlation(data_PCA, histogram=TRUE, pch=19)

# as there is a sustantial difference in day 0, print the number of transcripts that they survive count > 0

day0<-as.data.frame(EE1[c("input_WEN1","input_WEN2","input_WEN3")])
filter<- apply(day0, 1, function(x) length(x[x>0])>=3)
filtered.day0<-day0[filter,]
nr.day0<-nrow(filtered.day0)
# 22038

day1<-as.data.frame(EE1[c("input_WNN1","input_WNN2","input_WNN3")])
filter<- apply(day1, 1, function(x) length(x[x>0])>=3)
filtered.day1<-day1[filter,]
nr.day1<-nrow(filtered.day1)

b<-as.data.frame(c(nr.day0,nr.day1))
b$day<-c("EE","CTL")
colnames(b)<-c("N_transcripts","day")
pal <- wes_palette("Zissou1", length(b$day) , type = "continuous")
p<-ggplot(data=b,aes(x=day,y=N_transcripts)) +
 geom_bar(stat="identity", fill=pal, colour="black")


# differential analysis day 1 to day 0
data_PCA<-as.data.frame(EE1[c("MSC510d" , "MSC5_2_0d" , "MSC610d", "MSC6_2_0d", "MSC511d" , "MSC5_2_1d", "MSC611d", "MSC6_2_1d", "MSC513d" , "MSC6_2_3d", "MSC613d",  "MSC5_2_3d", "MSC5_1_7d" , "MSC5_2_7d", "MSC6_1_7d", "MSC6_2_7d", "MSC5_1_14d",  "MSC5_2_14d" , "MSC6_1_14d", "MSC6_2_14d" )])
rownames(data_PCA) <- make.names(EE1[,1], unique = TRUE)

EEa<-data_PCA
dir.create("enhancers")
OUTF1<-"./enhancers/"
filter<- apply(EEa, 1, function(x) length(x[x>4])>=5)
filtered<-EEa[filter,]
genes<-rownames(filtered)
x<-as.factor(rep(c("EE","CTL"),  c(3,3)))
set<-newSeqExpressionSet(as.matrix(filtered),
                        phenoData=data.frame(x,row.names=colnames(filtered)))
set



#########################
### UNNORMALIZE COUNTS###
#########################

# exploring the data without normalization
library(RColorBrewer)
colors<-brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF1, "boxplot_unnormalizecounts.pdf", sep=""))
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()
#saving
pdf(paste(OUTF1, "PCA_unnormalize_all.pdf", sep=""))
plotPCA(set, col=colors[x], cex=1.2)
dev.off()

##upper-quartile normalization
setn<-betweenLaneNormalization(set, which="upper")
plotRLE(setn, outline=FALSE, ylim=c(-2,2), col=colors[x])
plotPCA(setn, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF1, "boxplot_normalizecounts.pdf", sep=""))
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()
#saving
pdf(paste(OUTF1, "PCA_normalize.pdf", sep=""))
plotPCA(set, col=colors[x], cex=1.2)
dev.off()

################################################
######## RUVg USING CONTROL GENES###############
################################################


## no ERCC-spikes so then, using empirical controls :
design<-model.matrix(~x, data=pData(set))   # create a design for the study based in x matrix defined previously
y_empirical<-DGEList(counts=counts(set), group=x)
y_empirical<-calcNormFactors(y_empirical, method="upperquartile")
y_empirical<- estimateGLMCommonDisp(y_empirical, design)
y_empirical<- estimateGLMTagwiseDisp(y_empirical, design)

fit_empirical<-glmFit(y_empirical, design)

lrt_empirical<- glmLRT(fit_empirical, coef=2)
# saving unnormalized DE genes
tab_ruvg_10k <- topTags(lrt_empirical, n=length(rownames(lrt_empirical)))
write.csv(tab_ruvg_10k, paste(OUTF1, "lrt_empirical2.csv", sep=""))



# reporting the number of genes significantly up-regulated or downregulated at 5% FDR :

summary(dt_empirical<-decideTestsDGE(lrt_empirical))
write.table(summary(dt_empirical<-decideTestsDGE(lrt_empirical)), paste(OUTF1, "summary_uncorrected_FDR005.txt", sep=""))
isDE_empirical<-as.logical(dt_empirical)
DEnames_empirical<-rownames(y_empirical)[isDE_empirical]
plotSmear(lrt_empirical, de.tags=DEnames_empirical)
abline(h= c(-1,1), col="blue")

# plotting unnormalized DE genes
pdf(paste(OUTF1, "unbatched_DE_genes.pdf", sep=""))
plotSmear(lrt_empirical, de.tags=DEnames_empirical)
abline(h= c(-1,1), col="blue")
dev.off()

top<- topTags(lrt_empirical, n=nrow(set))$table
empirical<- rownames(set)[which(!(rownames(set) %in% rownames(top) [1:15000]))]

# we will use the most significant DE genes found by lrt_empirical that are less influenced by the variation
set2<- RUVg(set, empirical, k=1 )
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-2,2), col=colors[x])

#saving
pdf(paste(OUTF1, "boxplot_RUVg.pdf", sep=""))
plotRLE(set2, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()

plotPCA(set2, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF1, "PCA_RUVg.pdf", sep=""))
plotPCA(set2, col=colors[x], cex=1.2)
dev.off()

### reporting the DE genes in set2
design<-model.matrix(~x + W_1, data=pData(set2))
y_empirical_ruvg<-DGEList(counts=counts(set), group=x)
y_empirical_ruvg<-calcNormFactors(y_empirical_ruvg, method="upperquartile")
y_empirical_ruvg<-estimateGLMCommonDisp(y_empirical_ruvg, design)
y_empirical_ruvg<-estimateGLMTagwiseDisp(y_empirical_ruvg, design)

fit_empirical_ruvg<- glmFit(y_empirical_ruvg, design)
lrt_empirical_ruvg<-glmLRT(fit_empirical_ruvg, coef=2)

topTags(lrt_empirical_ruvg)

summary(dt_empirical_ruvg<-decideTestsDGE(lrt_empirical_ruvg))
write.table(summary(dt_empirical_ruvg<-decideTestsDGE(lrt_empirical_ruvg)), paste(OUTF1, "summary_RUVg_FDR005.txt", sep=""))
isDE_empirical_ruvg<-as.logical(dt_empirical_ruvg)
DEnames_empirical_ruvg<-rownames(y_empirical)[isDE_empirical_ruvg]
plotSmear(lrt_empirical_ruvg, de.tags=DEnames_empirical_ruvg)
abline(h= c(-1,1), col="blue")

#saving
pdf(paste(OUTF1, "DE_genes_RUVg.pdf", sep=""))
plotSmear(lrt_empirical_ruvg, de.tags=DEnames_empirical_ruvg)
abline(h= c(-1,1), col="blue")
dev.off()

#save into a table

tab_ruvg_200 <- topTags(lrt_empirical_ruvg, n=200)
tab_ruvg_500 <- topTags(lrt_empirical_ruvg, n=500)
tab_ruvg_10k <- topTags(lrt_empirical_ruvg, n=10000)
tab_ruvg_length <- topTags(lrt_empirical_ruvg, n=rownames(length(lrt_empirical_ruvg)))

write.table(tab_ruvg_length, paste(OUTF1, "tab_ruvg_test.txt", sep=""))
write.table(tab_ruvg_length, paste(OUTF1, "tab_ruvg_test.csv", sep=""))

dds <- DESeqDataSetFromMatrix(countData = counts(set2),colData = pData(set2),design = ~ W_1 + x)
dds <- DESeq(dds)
res  <- results(dds)
res
write.table(res, paste(OUTF1, "dds.DESEQ2.Wald.csv", sep=""))

dds2 <- DESeq(dds, test="LRT", reduced=as.formula("~ W_1"))
res <- results(dds2)
res

write.table(res, paste(OUTF1, "dds.DESEQ2.likehood.csv", sep=""))

################################################
######## RUVs USING REPLICATE SAMPLES ##########
################################################


#creating a matrix
differences<- matrix(data=c(1:2, 3:4), byrow=TRUE, nrow=2)
set3<-RUVs(set, genes, k=1, differences)
pData(set3)

plotRLE(set3, outline=FALSE, ylim=c(-2,2), col=colors[x])

#saving
pdf(paste(OUTF1, "boxplot_RUVs.pdf", sep=""))
plotRLE(set3, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()

plotPCA(set3, col=colors[x], cex=1.2)

pdf(paste(OUTF1, "PCA_RUVs.pdf", sep=""))
plotPCA(set3, col=colors[x], cex=1.2)
dev.off()

### reporting the DE genes in set3
design<-model.matrix(~x + W_1, data=pData(set3))
y_replicate_ruvs<-DGEList(counts=counts(set), group=x)
y_replicate_ruvs<-calcNormFactors(y_replicate_ruvs, method="upperquartile")
y_replicate_ruvs<-estimateGLMCommonDisp(y_replicate_ruvs, design)
y_replicate_ruvs<-estimateGLMTagwiseDisp(y_replicate_ruvs, design)

fit_replicate_ruvs<- glmFit(y_replicate_ruvs, design)
lrt_replicate_ruvs<-glmLRT(fit_replicate_ruvs, coef=2)

topTags(lrt_replicate_ruvs)

summary(dt_replicate_ruvs<-decideTestsDGE(lrt_replicate_ruvs))
write.table(summary(dt_replicate_ruvs<-decideTestsDGE(lrt_replicate_ruvs)), paste(OUTF1, "summary_RUVs_FDR005.txt", sep=""))
isDE_replicate_ruvs<-as.logical(dt_replicate_ruvs)
DEnames_replicate_ruvs<-rownames(y_empirical)[isDE_replicate_ruvs]
plotSmear(lrt_replicate_ruvs, de.tags=DEnames_replicate_ruvs)
abline(h= c(-1,1), col="blue")

#saving
pdf(paste(OUTF1, "DE_genes_RUVs.pdf", sep=""))
plotSmear(lrt_replicate_ruvs, de.tags=DEnames_replicate_ruvs)
abline(h= c(-1,1), col="blue")
dev.off()


#save into a table

tab_ruvs_200 <- topTags(lrt_replicate_ruvs, n=200)
tab_ruvs_500 <- topTags(lrt_replicate_ruvs, n=500)
tab_ruvs_10k <- topTags(lrt_replicate_ruvs, n=10000)
tab_ruvs_length <- topTags(lrt_replicate_ruvs, n=rownames(length(lrt_replicate_ruvs)))


# total togtags

write.table(tab_ruvs_length, paste(OUTF1, "tab_ruvs.txt", sep=""))
write.table(tab_ruvs_length, paste(OUTF1, "tab_ruvs.csv", sep=""))



################################################
######## RUVr USING RESIDUALS ##########
################################################



design<-model.matrix(~x, data=pData(set))
y<-DGEList(counts=counts(set), group=x)
y<-calcNormFactors(y, method="upperquartile")
y<-estimateGLMCommonDisp(y, design)
y<-estimateGLMTagwiseDisp(y, design)

fit<-glmFit(y, design)
res<- residuals(fit, type="deviance")

set4<-RUVr(set, genes, k=1, res)
pData(set4)


plotRLE(set4, outline=FALSE, ylim=c(-2,2), col=colors[x])

plotPCA(set4, col=colors[x], cex=1.2)

pdf(paste(OUTF1, "PCA_RUVr.pdf", sep=""))
plotPCA(set4, col=colors[x], cex=1.2)
dev.off()

### reporting the DE genes in set4
design<-model.matrix(~x + W_1, data=pData(set4))
y_residuals_ruvr<-DGEList(counts=counts(set), group=x)
y_residuals_ruvr<-calcNormFactors(y_residuals_ruvr, method="upperquartile")
y_residuals_ruvr<-estimateGLMCommonDisp(y_residuals_ruvr, design)
y_residuals_ruvr<-estimateGLMTagwiseDisp(y_residuals_ruvr, design)

fit_residuals_ruvr<- glmFit(y_residuals_ruvr, design)
lrt_residuals_ruvr<-glmLRT(fit_residuals_ruvr, coef=2)

topTags(lrt_residuals_ruvr)

summary(dt_residuals_ruvr<-decideTestsDGE(lrt_residuals_ruvr))
write.table(summary(dt_residuals_ruvr<-decideTestsDGE(lrt_residuals_ruvr)), paste(OUTF1, "summary_RUVr_FDR005.txt", sep=""))
isDE_residuals_ruvr<-as.logical(dt_residuals_ruvr)
DEnames_residuals_ruvr<-rownames(y_residuals_ruvr)[isDE_residuals_ruvr]
plotSmear(lrt_residuals_ruvr, de.tags=DEnames_residuals_ruvr)
abline(h= c(-1,1), col="blue")

# saving

pdf(paste(OUTF1, "DE_genes_RUVr.pdf", sep=""))
plotSmear(lrt_residuals_ruvr, de.tags=DEnames_residuals_ruvr)
abline(h= c(-1,1), col="blue")
dev.off()

#save into a table

tab_ruvr_200 <- topTags(lrt_residuals_ruvr, n=200)
tab_ruvr_500 <- topTags(lrt_residuals_ruvr, n=500)
tab_ruvr_10k<- topTags(lrt_residuals_ruvr, n=10000)
tab_ruvr_length <- topTags(lrt_residuals_ruvr, n=rownames(length(lrt_residuals_ruvr)))

write.table(tab_ruvr_length, paste(OUTF1, "tab_ruvr.txt", sep=""))
write.table(tab_ruvr_length, paste(OUTF1, "tab_ruvr.csv", sep=""))



#####

promoters

#########



setwd("")

# reading the table counts (done by featureCounts)

EE1 <- read.table("table.txt",row.names=NULL, header=TRUE)


library(ggbiplot)
library(RUVSeq)
library(devtools)
library(factoextra)
library(dplyr)
library(EDASeq)
library(digest)
library(Hmisc)
library(corrplot)
library("PerformanceAnalytics")
library("wesanderson")
library(DESeq2)

data_PCA<-as.data.frame(EE1[,-1])

rownames(data_PCA) <- make.names(EE1[,1], unique = TRUE)
length(data_PCA)

log.ir <- log(data_PCA[, 1:6])
filter<- apply(log.ir, 1, function(x) length(x[x>0])>=6)
filtered<-log.ir[filter,]
PCA_data_t<-as.data.frame(t(filtered))
which(apply(PCA_data_t, 2, var)==0)
test<-PCA_data_t[ , apply(PCA_data_t, 2, var) != 0]  # var estimation, if interger (0) no concerns


res.pca <- prcomp(test, scale = TRUE, center=TRUE)
fviz_eig(res.pca)

# fviz_pca_ind(res.pca,
#             col.ind = "cos2", # Colorer par le cos2
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE
#)



res.pca$group[(1:3)] <- "EE"
res.pca$group[(4:6)] <- "CTL"

fviz_pca_ind(res.pca,
            col.ind = res.pca$group, # colorer par groupes
            addEllipses = TRUE, # Ellipse de concentration
            ellipse.type = "confidence",
            legend.title = "Groups",
            repel = TRUE,
            geom="point"
)



# correlation

library(corrplot)
library("PerformanceAnalytics")

flattenCorrMatrix <- function(cormat, pmat) {
 ut <- upper.tri(cormat)
 data.frame(
   row = rownames(cormat)[row(cormat)[ut]],
   column = rownames(cormat)[col(cormat)[ut]],
   cor  =(cormat)[ut],
   p = pmat[ut]
   )
}

library(Hmisc)
res2<-rcorr(as.matrix(data_PCA[,1:9]))
flattenCorrMatrix(res2$r, res2$P)

res <- cor(data_PCA)
round(res, 2)
corrplot(res, type = "upper", order = "hclust",
        tl.col = "black", tl.srt = 45)
# pal <- wes_palette("Zissou1", length(mat_breaks) - 1, type = "continuous")
corrplot(res, method = "square", type = "upper", tl.col = "black", order = "hclust", addrect = 1 , rect.col= "white", col = wes_palette("Darjeeling1"))
# chart.Correlation(data_PCA, histogram=TRUE, pch=19)

# as there is a sustantial difference in day 0, print the number of transcripts that they survive count > 0

day0<-as.data.frame(EE1[c("input_WEN1","input_WEN2","input_WEN3")])
filter<- apply(day0, 1, function(x) length(x[x>0])>=3)
filtered.day0<-day0[filter,]
nr.day0<-nrow(filtered.day0)
# 22038

day1<-as.data.frame(EE1[c("input_WNN1","input_WNN2","input_WNN3")])
filter<- apply(day1, 1, function(x) length(x[x>0])>=3)
filtered.day1<-day1[filter,]
nr.day1<-nrow(filtered.day1)

b<-as.data.frame(c(nr.day0,nr.day1))
b$day<-c("EE","CTL")
colnames(b)<-c("N_transcripts","day")
pal <- wes_palette("Zissou1", length(b$day) , type = "continuous")
p<-ggplot(data=b,aes(x=day,y=N_transcripts)) +
 geom_bar(stat="identity", fill=pal, colour="black")



EEa<-data_PCA
dir.create("promoters")
OUTF1<-"./promoters/"
filter<- apply(EEa, 1, function(x) length(x[x>4])>=5)
filtered<-EEa[filter,]
genes<-rownames(filtered)
x<-as.factor(rep(c("EE","CTL"),  c(3,3)))
set<-newSeqExpressionSet(as.matrix(filtered),
                        phenoData=data.frame(x,row.names=colnames(filtered)))
set



#########################
### UNNORMALIZE COUNTS###
#########################

# exploring the data without normalization
library(RColorBrewer)
colors<-brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF1, "boxplot_unnormalizecounts.pdf", sep=""))
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()
#saving
pdf(paste(OUTF1, "PCA_unnormalize_all.pdf", sep=""))
plotPCA(set, col=colors[x], cex=1.2)
dev.off()

##upper-quartile normalization
setn<-betweenLaneNormalization(set, which="upper")
plotRLE(setn, outline=FALSE, ylim=c(-2,2), col=colors[x])
plotPCA(setn, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF1, "boxplot_normalizecounts.pdf", sep=""))
plotRLE(set, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()
#saving
pdf(paste(OUTF1, "PCA_normalize.pdf", sep=""))
plotPCA(set, col=colors[x], cex=1.2)
dev.off()

################################################
######## RUVg USING CONTROL GENES###############
################################################


## no ERCC-spikes so then, using empirical controls :
design<-model.matrix(~x, data=pData(set))   # create a design for the study based in x matrix defined previously
y_empirical<-DGEList(counts=counts(set), group=x)
y_empirical<-calcNormFactors(y_empirical, method="upperquartile")
y_empirical<- estimateGLMCommonDisp(y_empirical, design)
y_empirical<- estimateGLMTagwiseDisp(y_empirical, design)

fit_empirical<-glmFit(y_empirical, design)

lrt_empirical<- glmLRT(fit_empirical, coef=2)
# saving unnormalized DE genes
tab_ruvg_10k <- topTags(lrt_empirical, n=length(rownames(lrt_empirical)))
write.csv(tab_ruvg_10k, paste(OUTF1, "lrt_empirical2.csv", sep=""))



# reporting the number of genes significantly up-regulated or downregulated at 5% FDR :

summary(dt_empirical<-decideTestsDGE(lrt_empirical))
write.table(summary(dt_empirical<-decideTestsDGE(lrt_empirical)), paste(OUTF1, "summary_uncorrected_FDR005.txt", sep=""))
isDE_empirical<-as.logical(dt_empirical)
DEnames_empirical<-rownames(y_empirical)[isDE_empirical]
plotSmear(lrt_empirical, de.tags=DEnames_empirical)
abline(h= c(-1,1), col="blue")

# plotting unnormalized DE genes
pdf(paste(OUTF1, "unbatched_DE_genes.pdf", sep=""))
plotSmear(lrt_empirical, de.tags=DEnames_empirical)
abline(h= c(-1,1), col="blue")
dev.off()

top<- topTags(lrt_empirical, n=nrow(set))$table
empirical<- rownames(set)[which(!(rownames(set) %in% rownames(top) [1:15000]))]

# we will use the most significant DE genes found by lrt_empirical that are less influenced by the variation
set2<- RUVg(set, empirical, k=1 )
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-2,2), col=colors[x])

#saving
pdf(paste(OUTF1, "boxplot_RUVg.pdf", sep=""))
plotRLE(set2, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()

plotPCA(set2, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF1, "PCA_RUVg.pdf", sep=""))
plotPCA(set2, col=colors[x], cex=1.2)
dev.off()

### reporting the DE genes in set2
design<-model.matrix(~x + W_1, data=pData(set2))
y_empirical_ruvg<-DGEList(counts=counts(set), group=x)
y_empirical_ruvg<-calcNormFactors(y_empirical_ruvg, method="upperquartile")
y_empirical_ruvg<-estimateGLMCommonDisp(y_empirical_ruvg, design)
y_empirical_ruvg<-estimateGLMTagwiseDisp(y_empirical_ruvg, design)

fit_empirical_ruvg<- glmFit(y_empirical_ruvg, design)
lrt_empirical_ruvg<-glmLRT(fit_empirical_ruvg, coef=2)

topTags(lrt_empirical_ruvg)

summary(dt_empirical_ruvg<-decideTestsDGE(lrt_empirical_ruvg))
write.table(summary(dt_empirical_ruvg<-decideTestsDGE(lrt_empirical_ruvg)), paste(OUTF1, "summary_RUVg_FDR005.txt", sep=""))
isDE_empirical_ruvg<-as.logical(dt_empirical_ruvg)
DEnames_empirical_ruvg<-rownames(y_empirical)[isDE_empirical_ruvg]
plotSmear(lrt_empirical_ruvg, de.tags=DEnames_empirical_ruvg)
abline(h= c(-1,1), col="blue")

#saving
pdf(paste(OUTF1, "DE_genes_RUVg.pdf", sep=""))
plotSmear(lrt_empirical_ruvg, de.tags=DEnames_empirical_ruvg)
abline(h= c(-1,1), col="blue")
dev.off()

#save into a table

tab_ruvg_200 <- topTags(lrt_empirical_ruvg, n=200)
tab_ruvg_500 <- topTags(lrt_empirical_ruvg, n=500)
tab_ruvg_10k <- topTags(lrt_empirical_ruvg, n=10000)
tab_ruvg_length <- topTags(lrt_empirical_ruvg, n=rownames(length(lrt_empirical_ruvg)))

write.table(tab_ruvg_length, paste(OUTF1, "tab_ruvg_test.txt", sep=""))
write.table(tab_ruvg_length, paste(OUTF1, "tab_ruvg_test.csv", sep=""))

dds <- DESeqDataSetFromMatrix(countData = counts(set2),colData = pData(set2),design = ~ W_1 + x)
dds <- DESeq(dds)
res  <- results(dds)
res
write.table(res, paste(OUTF1, "dds.DESEQ2.Wald.csv", sep=""))

dds2 <- DESeq(dds, test="LRT", reduced=as.formula("~ W_1"))
res <- results(dds2)
res

write.table(res, paste(OUTF1, "dds.DESEQ2.likehood.csv", sep=""))

################################################
######## RUVs USING REPLICATE SAMPLES ##########
################################################


#creating a matrix
differences<- matrix(data=c(1:2, 3:4), byrow=TRUE, nrow=2)
set3<-RUVs(set, genes, k=1, differences)
pData(set3)

plotRLE(set3, outline=FALSE, ylim=c(-2,2), col=colors[x])

#saving
pdf(paste(OUTF1, "boxplot_RUVs.pdf", sep=""))
plotRLE(set3, outline=FALSE, ylim=c(-2,2), col=colors[x])
dev.off()

plotPCA(set3, col=colors[x], cex=1.2)

pdf(paste(OUTF1, "PCA_RUVs.pdf", sep=""))
plotPCA(set3, col=colors[x], cex=1.2)
dev.off()

### reporting the DE genes in set3
design<-model.matrix(~x + W_1, data=pData(set3))
y_replicate_ruvs<-DGEList(counts=counts(set), group=x)
y_replicate_ruvs<-calcNormFactors(y_replicate_ruvs, method="upperquartile")
y_replicate_ruvs<-estimateGLMCommonDisp(y_replicate_ruvs, design)
y_replicate_ruvs<-estimateGLMTagwiseDisp(y_replicate_ruvs, design)

fit_replicate_ruvs<- glmFit(y_replicate_ruvs, design)
lrt_replicate_ruvs<-glmLRT(fit_replicate_ruvs, coef=2)

topTags(lrt_replicate_ruvs)

summary(dt_replicate_ruvs<-decideTestsDGE(lrt_replicate_ruvs))
write.table(summary(dt_replicate_ruvs<-decideTestsDGE(lrt_replicate_ruvs)), paste(OUTF1, "summary_RUVs_FDR005.txt", sep=""))
isDE_replicate_ruvs<-as.logical(dt_replicate_ruvs)
DEnames_replicate_ruvs<-rownames(y_empirical)[isDE_replicate_ruvs]
plotSmear(lrt_replicate_ruvs, de.tags=DEnames_replicate_ruvs)
abline(h= c(-1,1), col="blue")

#saving
pdf(paste(OUTF1, "DE_genes_RUVs.pdf", sep=""))
plotSmear(lrt_replicate_ruvs, de.tags=DEnames_replicate_ruvs)
abline(h= c(-1,1), col="blue")
dev.off()


#save into a table

tab_ruvs_200 <- topTags(lrt_replicate_ruvs, n=200)
tab_ruvs_500 <- topTags(lrt_replicate_ruvs, n=500)
tab_ruvs_10k <- topTags(lrt_replicate_ruvs, n=10000)
tab_ruvs_length <- topTags(lrt_replicate_ruvs, n=rownames(length(lrt_replicate_ruvs)))


# total togtags

write.table(tab_ruvs_length, paste(OUTF1, "tab_ruvs.txt", sep=""))
write.table(tab_ruvs_length, paste(OUTF1, "tab_ruvs.csv", sep=""))



################################################
######## RUVr USING RESIDUALS ##########
################################################



design<-model.matrix(~x, data=pData(set))
y<-DGEList(counts=counts(set), group=x)
y<-calcNormFactors(y, method="upperquartile")
y<-estimateGLMCommonDisp(y, design)
y<-estimateGLMTagwiseDisp(y, design)

fit<-glmFit(y, design)
res<- residuals(fit, type="deviance")

set4<-RUVr(set, genes, k=1, res)
pData(set4)


plotRLE(set4, outline=FALSE, ylim=c(-2,2), col=colors[x])

plotPCA(set4, col=colors[x], cex=1.2)

pdf(paste(OUTF1, "PCA_RUVr.pdf", sep=""))
plotPCA(set4, col=colors[x], cex=1.2)
dev.off()

### reporting the DE genes in set4
design<-model.matrix(~x + W_1, data=pData(set4))
y_residuals_ruvr<-DGEList(counts=counts(set), group=x)
y_residuals_ruvr<-calcNormFactors(y_residuals_ruvr, method="upperquartile")
y_residuals_ruvr<-estimateGLMCommonDisp(y_residuals_ruvr, design)
y_residuals_ruvr<-estimateGLMTagwiseDisp(y_residuals_ruvr, design)

fit_residuals_ruvr<- glmFit(y_residuals_ruvr, design)
lrt_residuals_ruvr<-glmLRT(fit_residuals_ruvr, coef=2)

topTags(lrt_residuals_ruvr)

summary(dt_residuals_ruvr<-decideTestsDGE(lrt_residuals_ruvr))
write.table(summary(dt_residuals_ruvr<-decideTestsDGE(lrt_residuals_ruvr)), paste(OUTF1, "summary_RUVr_FDR005.txt", sep=""))
isDE_residuals_ruvr<-as.logical(dt_residuals_ruvr)
DEnames_residuals_ruvr<-rownames(y_residuals_ruvr)[isDE_residuals_ruvr]
plotSmear(lrt_residuals_ruvr, de.tags=DEnames_residuals_ruvr)
abline(h= c(-1,1), col="blue")

# saving

pdf(paste(OUTF1, "DE_genes_RUVr.pdf", sep=""))
plotSmear(lrt_residuals_ruvr, de.tags=DEnames_residuals_ruvr)
abline(h= c(-1,1), col="blue")
dev.off()

#save into a table

tab_ruvr_200 <- topTags(lrt_residuals_ruvr, n=200)
tab_ruvr_500 <- topTags(lrt_residuals_ruvr, n=500)
tab_ruvr_10k<- topTags(lrt_residuals_ruvr, n=10000)
tab_ruvr_length <- topTags(lrt_residuals_ruvr, n=rownames(length(lrt_residuals_ruvr)))

write.table(tab_ruvr_length, paste(OUTF1, "tab_ruvr.txt", sep=""))
write.table(tab_ruvr_length, paste(OUTF1, "tab_ruvr.csv", sep=""))
