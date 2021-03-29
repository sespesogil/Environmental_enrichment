#' Computing Pearson correlation to asses reproducibilty 
#' Directly from FeatureCounts read counts table

library(ggplot2)

dir<-""
files<-list.files(dir, "total_counts.txt", recursive = F)
files

for (i in 1:length(files)){
  file<-files[i]
  filename<-basename(file)
  print(filename)
  if ( i == 1) {
  tmp<-read.table(file, header=T, sep = "\t")
  tmp.samples<-tmp[, -c(1:6)]
  head(tmp.samples)
  cor.matrix<-cor(tmp.samples)
  cor.df<-as.data.frame(as.table(cor.matrix))
  cor.df.pw<-subset(cor.df, Freq < 1)
  cor.df.pw$group<-gsub(".total_counts.txt", "", filename)
  } else {
    tmp<-read.table(file, header=T, sep = "\t")
    tmp.samples<-tmp[, -c(1:6)]
    head(tmp.samples)
    cor.matrix<-cor(tmp.samples)
    cor.df<-as.data.frame(as.table(cor.matrix))
    tmp.cor.df.pw<-subset(cor.df, Freq < 1)
    tmp.cor.df.pw$group<-gsub(".total_counts.txt", "", filename)
    cor.df.pw<-rbind(cor.df.pw,tmp.cor.df.pw)
  }
}

library(ggplot2)
p <- ggplot(cor.df.pw, aes(x=group, y=Freq, color=group)) + 
  geom_boxplot() + ylim(0,1)
p
