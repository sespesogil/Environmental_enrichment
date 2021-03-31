#Example of MarkerGeneProfile estimation 

devtools::install_github('oganm/markerGeneProfile')

library('ggplot2')
library('gplots')
library('viridis')
library('dplyr')
library('knitr')


data(mouseMarkerGenes)
names(mouseMarkerGenes)
library(markerGeneProfile)

lapply(mouseMarkerGenes$Cortex[1:3],head, 14)

library(dplyr)

mgp_LesnickParkinsonsExp <- read.table("H3K79me2_neuNpos_neg_rpkms.txt", header=T)
mgp_LesnickParkinsonsMeta<- read.table("cell_types.txt", header=T)

mgp_LesnickParkinsonsMeta




    unfilteredParkinsonsExp = mgp_LesnickParkinsonsExp # keep this for later
medExp = mgp_LesnickParkinsonsExp %>% 
    sepExpr() %>% {.[[2]]} %>%
    unlist %>% median

# mostVariable function is part of this package that does probe selection and filtering for you
mgp_LesnickParkinsonsExp = mostVariable(mgp_LesnickParkinsonsExp, 
                                        threshold = medExp, 
                                        threshFun= median)



estimations =  mgpEstimate(exprData=mgp_LesnickParkinsonsExp,
                           genes=mouseMarkerGenes$Cortex,
                           geneColName='Gene.Symbol',
                           outlierSampleRemove=F, # should outlier samples removed. This is done using boxplot stats.
                           geneTransform =NULL, # this is the default option for geneTransform
                           groups=mgp_LesnickParkinsonsMeta$disease, #if there are experimental groups provide them here. if not desired set to NULL
                           seekConsensus = FALSE, # ensures gene rotations are positive in both of the groups
                           removeMinority = TRUE) 

ls(estimations$estimates)

OUTF<-"/CHIP-seq/"

OligoFrame =
    data.frame(`OLigo MGP` = estimations$estimates$Oligo, 
               state = estimations$groups$Oligo, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
               check.names=FALSE)


ggplot2::ggplot(OligoFrame, 
                aes(x = state, y = `OLigo MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t

AstrocyteFrame =
data.frame(`Astrocyte MGP` = estimations$estimates$Astrocyte, 
state = estimations$groups$Astrocyte, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(AstrocyteFrame, 
                aes(x = state, y = `Astrocyte MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t



GabaPVFrame =
data.frame(`GabaPV MGP` = estimations$estimates$GabaPV, 
state = estimations$groups$GabaPV, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(GabaPVFrame, 
                aes(x = state, y = `GabaPV MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t



GabaVIPRelnFrame =
data.frame(`GabaVIPReln MGP` = estimations$estimates$GabaVIPReln, 
state = estimations$groups$GabaVIPReln, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)



ggplot2::ggplot(GabaVIPRelnFrame, 
                aes(x = state, y = `GabaVIPReln MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t
dev.off()


PyramidalFrame =
data.frame(`Pyramidal MGP` = estimations$estimates$Pyramidal, 
state = estimations$groups$Pyramidal, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(PyramidalFrame, 
                aes(x = state, y = `Pyramidal MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t




GabaRelnCalbFrame =
data.frame(`PGabaRelnCalb MGP` = estimations$estimates$GabaRelnCalb, 
state = estimations$groups$GabaRelnCalb, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(GabaRelnCalbFrame, 
                aes(x = state, y = `GabaRelnCalb MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t



Endothelial =
data.frame(`Endothelial MGP` = estimations$estimates$Endothelial, 
state = estimations$groups$Endothelial, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(Endothelial, 
                aes(x = state, y = `Endothelial MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t



Microglia_activation =
data.frame(`Microglia_activation MGP` = estimations$estimates$Microglia_activation, 
state = estimations$groups$Microglia_activation, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(Microglia_activation, 
                aes(x = state, y = `Microglia_activation MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t



Microglia_deactivation =
data.frame(`Microglia_deactivation MGP` = estimations$estimates$Microglia_deactivation, 
state = estimations$groups$Microglia_deactivation, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(Microglia_deactivation, 
                aes(x = state, y = `Microglia_deactivation MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t



Microglia =
data.frame(`Microglia MGP` = estimations$estimates$Microglia, 
state = estimations$groups$Microglia, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(Microglia, 
                aes(x = state, y = `Microglia MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t


OligoPrecursors =
data.frame(`OligoPrecursors MGP` = estimations$estimates$OligoPrecursors, 
state = estimations$groups$OligoPrecursors, # note that unless outlierSampleRemove is TRUE this will be always the same as the groups input
check.names=FALSE)


ggplot2::ggplot(OligoPrecursors, 
                aes(x = state, y = `OligoPrecursors MGP`)) + 
    geom_boxvio() + geom_jitter(width = .05) # t


