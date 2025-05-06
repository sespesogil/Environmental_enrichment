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


######################################################
######################################################
###### RUV-SEQ: removing unwanted variation. #########
######################################################
######################################################

## R version 3.2.0
## website: http://www.bioconductor.org/packages/release/bioc/html/RUVSeq.html

## current script: sergio.espeso@crg.eu


system("mkdir /RUV/WNN_vs_WEN/")
OUTF<-"/analysis/mRNA/RUV/WNN_vs_WEN/"

library(RUVSeq)

# reading the table counts (done by featureCounts)
EE<-read.table("/transcriptome/analysis/mRNA/featureCounts/WNN_vs_WEN/counts.txt", row.names=1 ,  header=T )

# changing the format of the table
colnames(EE)<-c( "Chr", "Start", "End", "Strand", "Length", "WNN1", "WNN2", "WEN1", "WEN2")
EEb <- EE[c('WNN1','WNN2','WEN1','WEN2')]

## filtering and exploratory analysis
filter<- apply(EEb, 1, function(x) length(x[x>5])>=2)
filtered<-EEb[filter,]
genes<-rownames(filtered)

## if ERCC spikes are present, then:
# genes<-rownames(filtered)[grep("^ENS", rownames(filtered))]
# spikes<-rownames(filtered)[grep("^ERCC", rownames(filtered))]

## creating a matrix of conditions and replicates
x<-as.factor(rep(c("WNN","WEN"), each=2))

## storing the data into a S4 object from the EDAseq package. Usfull for plotting
set<-newSeqExpressionSet(as.matrix(filtered),
phenoData=data.frame(x,row.names=colnames(filtered)))
set

#########################
### UNNORMALIZE COUNTS###
#########################

# exploring the data without normalization
library(RColorBrewer)
colors<-brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#saving 
pdf(paste(OUTF, "boxplot_unnormalizecounts.pdf", sep=""))
plotRLE(set, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
dev.off()
#saving
pdf(paste(OUTF, "PCA_unnormalize.pdf", sep=""))
plotPCA(set, col=colors[x], cex=1.2)
dev.off()

##upper-quartile normalization
setn<-betweenLaneNormalization(set, which="upper")
plotRLE(setn, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
plotPCA(setn, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF, "boxplot_normalizecounts.pdf", sep=""))
plotRLE(set, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
dev.off()
#saving
pdf(paste(OUTF, "PCA_normalize.pdf", sep=""))
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

# reporting the number of genes significantly up-regulated or downregulated at 5% FDR :

summary(dt_empirical<-decideTestsDGE(lrt_empirical))
isDE_empirical<-as.logical(dt_empirical)
DEnames_empirical<-rownames(y_empirical)[isDE_empirical]
plotSmear(lrt_empirical, de.tags=DEnames_empirical)
abline(h= c(-1,1), col="blue")

# plotting unnormalized DE genes
pdf(paste(OUTF, "unbatched_DE_genes.pdf", sep=""))
plotSmear(lrt_empirical, de.tags=DEnames_empirical)
abline(h= c(-1,1), col="blue")
dev.off()

top<- topTags(lrt_empirical, n=nrow(set))$table
empirical<- rownames(set)[which(!(rownames(set) %in% rownames(top) [1:15000]))]

# we will use the most significant DE genes found by lrt_empirical that are less influenced by the variation 
set2<- RUVg(set, empirical, k=1 )
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])

#saving
pdf(paste(OUTF, "boxplot_RUVg.pdf", sep=""))
plotRLE(set2, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
dev.off()

plotPCA(set2, col=colors[x], cex=1.2)

#saving
pdf(paste(OUTF, "PCA_RUVg.pdf", sep=""))
plotRLE(set2, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
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
isDE_empirical_ruvg<-as.logical(dt_empirical_ruvg)
DEnames_empirical_ruvg<-rownames(y_empirical)[isDE_empirical_ruvg]
plotSmear(lrt_empirical_ruvg, de.tags=DEnames_empirical_ruvg)
abline(h= c(-1,1), col="blue")

#saving
pdf(paste(OUTF, "DE_genes_RUVg.pdf", sep=""))
plotSmear(lrt_empirical_ruvg, de.tags=DEnames_empirical_ruvg)
abline(h= c(-1,1), col="blue")
dev.off()

#annotating
library(org.Mm.eg.db)
idfound<-rownames(lrt_empirical_ruvg) %in% mappedRkeys(org.Mm.egENSEMBL)
list_ruvg<-lrt_empirical_ruvg[idfound,]
dim(list_ruvg)

egENSEMBL<-toTable(org.Mm.egENSEMBL)
head(egENSEMBL)
m_ruvg<- match(rownames(lrt_empirical_ruvg), egENSEMBL$ensembl_id)
lrt_empirical_ruvg$genes$EntrezGene <- egENSEMBL$gene_id[m_ruvg]

egSYMBOL<- toTable(org.Mm.egSYMBOL)
head(egSYMBOL)
m_ruvg<-match(lrt_empirical_ruvg$genes$EntrezGene, egSYMBOL$gene_id)
lrt_empirical_ruvg$genes$Symbol<- egSYMBOL$symbol[m_ruvg]
head(lrt_empirical_ruvg$genes)

#save into a table

tab_ruvg_200 <- topTags(lrt_empirical_ruvg, n=200)
tab_ruvg_500 <- topTags(lrt_empirical_ruvg, n=500)

# first 200 more significant

write.table(tab_ruvg_200, file="/RUV/WNN_vs_WEN/FIRST_200_RUVg_mygenelist.txt")

write.csv(tab_ruvg_200, file="/RUV/WNN_vs_WEN/FIRST_200_RUVg_mygenelist.csv")

# first 500 more significant

write.table(tab_ruvg_500, file="/WNN_vs_WEN/FIRST_500_RUVg_mygenelist.txt")

write.csv(tab_ruvg_500, file="/WNN_vs_WEN/FIRST_500_RUVg_mygenelist.csv")

#full list 
write.table(lrt_empirical_ruvg, file="/WNN_vs_WEN/RUVg_mygenelist.txt")

write.csv(lrt_empirical_ruvg, file="/WNN_vs_WEN/RUVg_mygenelist.csv")



################################################
######## RUVs USING REPLICATE SAMPLES ##########
################################################


#creating a matrix
differences<- matrix(data=c(1:2, 3:4), byrow=TRUE, nrow=2)
set3<-RUVs(set, genes, k=1, differences)
pData(set3)

plotRLE(set3, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])

#saving
pdf(paste(OUTF, "boxplot_RUVs.pdf", sep=""))
plotRLE(set3, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
dev.off()

plotPCA(set3, col=colors[x], cex=1.2)

pdf(paste(OUTF, "PCA_RUVs.pdf", sep=""))
plotRLE(set3, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
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
isDE_replicate_ruvs<-as.logical(dt_replicate_ruvs)
DEnames_replicate_ruvs<-rownames(y_empirical)[isDE_replicate_ruvs]
plotSmear(lrt_replicate_ruvs, de.tags=DEnames_replicate_ruvs)
abline(h= c(-1,1), col="blue")

#saving
pdf(paste(OUTF, "DE_genes_RUVg.pdf", sep=""))
plotSmear(lrt_replicate_ruvs, de.tags=DEnames_replicate_ruvs)
abline(h= c(-1,1), col="blue")
dev.off()

#annotating
library(org.Mm.eg.db)
idfound<-rownames(lrt_replicate_ruvs) %in% mappedRkeys(org.Mm.egENSEMBL)
list_ruvg<-lrt_replicate_ruvs[idfound,]
dim(list_ruvg)

egENSEMBL<-toTable(org.Mm.egENSEMBL)
head(egENSEMBL)
m_ruvg<- match(rownames(lrt_empirical_ruvg), egENSEMBL$ensembl_id)
lrt_empirical_ruvg$genes$EntrezGene <- egENSEMBL$gene_id[m_ruvg]

egSYMBOL<- toTable(org.Mm.egSYMBOL)
head(egSYMBOL)
m_ruvg<-match(lrt_empirical_ruvg$genes$EntrezGene, egSYMBOL$gene_id)
lrt_empirical_ruvg$genes$Symbol<- egSYMBOL$symbol[m_ruvg]
head(lrt_empirical_ruvg$genes)

#save into a table

tab_ruvg_200 <- topTags(lrt_empirical_ruvg, n=200)
tab_ruvg_500 <- topTags(lrt_empirical_ruvg, n=500)

# first 200 more significant

write.table(tab_ruvg_200, file="/WNN_vs_WEN/FIRST_200_RUVg_mygenelist.txt")

write.csv(tab_ruvg_200, file="/WNN_vs_WEN/FIRST_200_RUVg_mygenelist.csv")

# first 500 more significant

write.table(tab_ruvg_500, file="/WNN_vs_WEN/FIRST_500_RUVg_mygenelist.txt")

write.csv(tab_ruvg_500, file="/WNN_vs_WEN/FIRST_500_RUVg_mygenelist.csv")

#full list 
write.table(lrt_empirical_ruvg, file="/WNN_vs_WEN/RUVg_mygenelist.txt")

write.csv(lrt_empirical_ruvg, file="/WNN_vs_WEN/RUVg_mygenelist.csv")

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


plotRLE(set4, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])

#saving
OUTF<-"/analysis/mRNA/RUV/"
pdf(paste(OUTF, "boxplot_RUVr.pdf", sep=""))
plotRLE(set4, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
dev.off()

plotPCA(set4, col=colors[x], cex=1.2)

OUTF<-"/analysis/mRNA/RUV/"
pdf(paste(OUTF, "PCA_RUVr.pdf", sep=""))
plotRLE(set4, outline=FALSE, ylim=c(-0.4,0.4), col=colors[x])
dev.off()



OUTF<-"/transcriptome/analysis/mRNA/RUV/"
pdf(paste(OUTF, "DE_plotsmear.pdf", sep=""))
plotSmear(lrt, de.tags=DEnames)
abline(h=c(-1,1), col="blue")
dev.off()


