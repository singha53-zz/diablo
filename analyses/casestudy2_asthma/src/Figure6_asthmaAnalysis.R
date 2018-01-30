#################################################
#
# Figure6_asthmaAnalysis.R
# Author: Amrit Singh
# Date: April 01, 2016
#
#################################################
WhereAmI <- "~/Dropbox/Manuscript/mixOmics.org:DIABLO/"

## load libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(mixOmicsv6)
library(corrplot)
library(annotate)
library(pathview)
library(pROC)
library("hugene10sttranscriptcluster.db")
source("~/Dropbox/Manuscript/mixOmics.org:DIABLO/functions/visualizationFunctions.R")
source("~/Dropbox/Manuscript/mixOmics.org:DIABLO/functions/functions.R")
source(paste0(WhereAmI, "functions/diablo.perf_diablo.R"))

## load data
load("~/Dropbox/Manuscript/mixOmics.org:DIABLO/1-Data & Preprocessing/asthma/data/asthmaDatasets.RDATA")

## clean datasets
cells <- cells[, c("Relative.Neutrophils",
                   "Relative.Lymphocytes","Relative.Monocytes",
                   "Relative.Eosinophils","Relative.Basophils",
                   "Treg","Tcells","Bcells","Th17")]
metExp <- metExp[apply(metExp, 1, mad) >0, ]
rownames(metExp) <- metabolites[rownames(metExp), "BIOCHEMICAL"]

##-----------------------------------------------------------------------------
#
# 1) Figure 2A: FEV1 profiles
#
##-----------------------------------------------------------------------------
fev1 <- read.csv("/Users/asingh/Dropbox/manuscript/mixOmics.org:DIABLO/1-Data & Preprocessing/asthma/fev1Data.csv", row.names = 1)[, c("BLFEV","F10L","F20L","F30L","F45L","F60L","F90L","F120L","F180L","F240L","F300L","F360L","F420L")]

p <- scale(t(fev1), center = fev1$BLFEV, scale = fev1$BLFEV) %>% as.data.frame() %>%
  tbl_df() %>%
  mutate(Time = c(0, 10, 20, 30, 45, 60, 90, 120, 180 ,240, 300, 360, 420)) %>%
  tidyr::gather(Subject, FEV1, -Time) %>%
  filter(Time %in% c(0, 10, 20, 30, 45, 60, 90, 120)) %>%
  mutate(fev1 = 100*FEV1) %>%
  ggplot(aes(x = Time, y = fev1, fill = Subject, color = Subject)) +
  geom_point() + geom_line(color="black")
  
pdf(paste0(WhereAmI, "7-Figure6_asthmaAnalysis/Figure6A.pdf"))
p + scale_y_continuous(expression('Percent drop in '~ FEV[1])) + theme_bw() +
  scale_x_continuous(expression('Time (minutes)')) +
  theme(axis.text.y = element_text(size = 15, hjust = 1)) + theme(axis.text.x = element_text(size = 15, hjust = 0.5))+
  theme(axis.title.x=element_text(size = 15)) + theme(axis.title.y=element_text(size = 15,angle = 90))+ 
  theme(legend.key = element_rect(colour = "black", fill="white"))  +
  theme(plot.background = element_rect()) +  
  theme(strip.background = element_rect(colour = "black", fill = "white",
                                        size = 1), strip.text.x = element_text(size=20)) +
  geom_hline(yintercept = 0, colour="yellow3", linetype = "longdash") +
  theme(legend.position="none")
dev.off()

###########################################################################
#
# run diablo and generate visuals
#    pre vs. post
#
###########################################################################
##-------------------------------------
##------------- unpaired
##-------------------------------------
Y <- relevel(demo$Time, ref = "pre")
X <- list(Cells = cells, geneMods = genMEs0, metMods = metMEs0)
dim(X[[1]]); dim(X[[2]]); dim(X[[3]]); 

## compute design matrix
## print design matrix to pdf
x.xList <- list()
for(i in 1:length(X)){
  corDat <- rep(0, length(X))
  names(corDat) <- paste("cor", names(X)[i], names(X), sep = "_")
  for(j in 1:length(X)){
    result <- pls(X = X[[i]], Y = X[[j]], ncomp = 1)
    corDat[j] <- as.numeric(cor(result$variates$X, result$variates$Y))
  }
  x.xList[[i]] <- corDat
}
corMat <- do.call(rbind, x.xList)
rownames(corMat) <- colnames(corMat) <- names(X)
pdf(paste0(WhereAmI, "7-Figure6_asthmaAnalysis/unpairedDesignMatrixFigureS3.pdf"), width = 5, height = 5)
corrplot(corMat)
corrplot(corMat,add=TRUE, type="upper", method="ell",order="original",
  diag=FALSE,tl.pos="d", cl.pos="n")
corrplot(corMat,add=TRUE, type="lower", method="number",order="original",
  diag=TRUE,tl.pos="d", cl.pos="n")
dev.off()


## run DIABLO
ncomp <- rep(2, length(X))
design <- matrix(c(0, 0, 1,
 				   0, 0, 1,
 				   1, 1, 0), nrow = 3, ncol = 3)
keepX = list(rep(1, ncomp[1]), rep(5, ncomp[1]), rep(5, ncomp[1]))
result = diablo(X = X, Y = Y, ncomp = ncomp, 
                            keepX = keepX, design = design,
                            mode = "regression", bias = TRUE)
feat1 <- lapply(result$loadings, function(x)
  apply(x, 2, function(i) names(i)[which(i != 0)]))
plotDiablo(result)

plotIndiv_diablo(object=result, ncomp = 1, groupOrder = c("pre", "post"))

## average error rate
perf=perf_diablo(result, validation="loo")
perf$AverageScore.error.rate

## plot area under the curve
time <- relevel(demo$Time, ref = "pre")
names(time) <- rownames(demo)
predictScores <- Reduce("+", lapply(perf$predict, function(i) i[[1]]))/length(X)
all(rownames(predictScores) == names(time))

roc.score = roc(response = as.character(time), predictor = predictScores[, "post"], plot = TRUE, percent = TRUE, na.rm =TRUE)
roc.res = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.res$Specificity = 100 - as.numeric(roc.res$Specificity)
roc.res$Sensitivity = as.numeric(roc.res$Sensitivity)
roc.score$auc

##-------------------------------------
##------------- paired
##-------------------------------------
Cov <- data.frame(sample = rep(1:14, 2), stimul = Y)
A = lapply(X, function(i) suppressMessages(withinVariation(X = i, design = Cov)))

##-----------------------------------------------------------------------------
#
# plot design matrix
#
##-----------------------------------------------------------------------------
## determine design
## compute design matrix
## print design matrix to pdf
x.xList <- list()
for(i in 1:length(A)){
  corDat <- rep(0, length(A))
  names(corDat) <- paste("cor", names(A)[i], names(A), sep = "_")
  for(j in 1:length(A)){
    result <- pls(X = A[[i]], Y = A[[j]], ncomp = 1)
    corDat[j] <- as.numeric(cor(result$variates$X, result$variates$Y))
  }
  x.xList[[i]] <- corDat
}
corMat <- do.call(rbind, x.xList)
rownames(corMat) <- colnames(corMat) <- names(A)
pdf(paste0(WhereAmI, "7-Figure6_asthmaAnalysis/pairedDesignMatrixFigureS3.pdf"), width = 5, height = 5)
corrplot(corMat)
corrplot(corMat,add=TRUE, type="upper", method="ell",order="original",
  diag=FALSE,tl.pos="d", cl.pos="n")
corrplot(corMat,add=TRUE, type="lower", method="number",order="original",
  diag=TRUE,tl.pos="d", cl.pos="n")
dev.off()

##-----------------------------------------------------------------------------
#
# Figure 6B: path diagram (made in ppt)
#
##-----------------------------------------------------------------------------

## run Diablo
ncomp <- c(2, 2, 2)
keepX = list(rep(1, ncomp[1]), rep(5, ncomp[2]), rep(5, ncomp[3]))
result2 = diablo(X = A, Y = Y, ncomp = ncomp, 
                            keepX = keepX, design = design,
                            mode = "regression", bias = TRUE)
feat2 <- lapply(result2$loadings, function(x)
  apply(x, 2, function(i) names(i)[which(i != 0)]))

##-----------------------------------------------------------------------------
#
# Figure 6C: component plots
#
##-----------------------------------------------------------------------------
pdf(paste0(WhereAmI, "7-Figure6_asthmaAnalysis/Figure6C.pdf"), width = 7, height = 5)
plotIndiv_diablo(object=result2, ncomp = 1, groupOrder = c("pre", "post"))
dev.off()


##-----------------------------------------------------------------------------
#
# circos plot
#
##-----------------------------------------------------------------------------
plotDiablo3(result2, ncomp = 1)


circosPlot_diablo(result2, corThreshold=0.75, cex.label = 0.5)


##-----------------------------------------------------------------------------
#
# Figure 6D: heatmap (correlation between features)
#
##-----------------------------------------------------------------------------
object = result2
margins = c(2, 15)
pos.legend = "topright"
cex.legend = 1.5
pdf(paste0(WhereAmI, "7-Figure6_asthmaAnalysis/Figure6D.pdf"), width = 10)
#cimDiablo2(result2, cex.legend = 0.9, margins = c(2, 17))
X <- object$X
Y <- object$Y
keepA <- lapply(object$loadings, function(x)
  apply(x, 2, function(i) names(i)[which(i != 0)]))
XDatList <- mapply(function(x, y) {
  x[, y]
}, x = X, y = keepA[-length(keepA)], SIMPLIFY = FALSE)

XDat <- do.call(cbind, XDatList)
XDat[which(XDat > 2)] <- 2
XDat[which(XDat < -2)] <- -2
dark <- brewer.pal(n = 12, name = "Paired")[seq(2, 12, by = 2)]
col = rep(1:length(X), unlist(lapply(keepA[-length(keepA)], length)))+3

#VarLabels <- factor(rep(names(X), lapply(keepA[-length(keepA)],
#sum)), levels = names(X)[order(names(X))])
opar <- par()[!names(par()) %in% c("cin", "cra", "csi", "cxy",
  "din", "page")]
par(mfrow = c(1, 1))
cim(cor(XDat), col.names = rep("", ncol(XDat)), row.sideColors = color.mixo(col),
  col.sideColors = color.mixo(col), margins = margins)
legend("topright", names(X), col = color.mixo(4:6), pch = 19, bty = "n")
par(opar)
dev.off()


##-----------------------------------------------------------------------------
## Compute error rate
cvPaired <- perf_diablo(result2, validation = "loo")
cvPaired$AveragePredict.error.rate

## plot area under the curve
time <- relevel(demo$Time, ref = "pre")
names(time) <- rownames(demo)
predictScores <- Reduce("+", lapply(cvPaired$predict, function(i) i[[2]]))/length(X)
all(rownames(predictScores) == names(time))

roc.score = roc(response = as.character(time), predictor = predictScores[, "post"], plot = TRUE, percent = TRUE, na.rm =TRUE)
roc.resPaired = data.frame(Specificity = rev(roc.score$specificities), Sensitivity = rev(roc.score$sensitivities))
roc.resPaired$Specificity = 100 - as.numeric(roc.resPaired$Specificity)
roc.resPaired$Sensitivity = as.numeric(roc.resPaired$Sensitivity)
roc.resPaired$auc


##-----------------------------------------------------------------------------
#
# Figure 6E: plot AUCs
#
##-----------------------------------------------------------------------------
df.auc <- data.frame(Specificity = c(roc.res$Specificity, roc.resPaired$Specificity),
  Sensitivity = c(roc.res$Sensitivity, roc.resPaired$Sensitivity))
df.auc$AnalysisType <- rep(c("Unpaired", "Paired"), each = length(roc.res$Sensitivity))

pdf(paste0(WhereAmI, "7-Figure6_asthmaAnalysis/Figure6E.pdf"), width = 7, height = 5)
ggplot(df.auc, aes(y = Sensitivity, x = Specificity, color = AnalysisType)) + geom_point(size = 3) + geom_line(size = 2) +
  customTheme(sizeStripFont = 15, xAngle = 0, hjust = 0.5, vjust = 0.5, xSize = 15, ySize = 15, xAxisSize = 15, yAxisSize = 15) +
  xlab("100 - Specificity") + 
  annotate("text", x = 30, y = 95, label = "AUC = 99.5%", size = 5) +
  annotate("text", x = 30, y = 74, label = "AUC = 85.2%", size = 5)
dev.off()




##-----------------------------------------------------------------------------
#
# Supplementary Figure 6: Asthma KEGG diagram
#
##-----------------------------------------------------------------------------
## pre vs post
fc <- rowMeans(genExp[, demo$Time == "post"])-rowMeans(genExp[, demo$Time == "pre"])
geneSymbols <- getSYMBOL(names(fc), "hugene10sttranscriptcluster")
eID <- id2eg(geneSymbols)
names(fc) <- eID[,2]
pv.out <- pathview(gene.data = fc, pathway.id = "05310",
                   species = "hsa",  kegg.dir= "/Users/asingh/Dropbox/manuscript/mixOmics.org:DIABLO/5-SupplementaryMaterial/",
                   #out.suffix = "SupplementaryFig.6A_PrevsPost", 
                   kegg.native = T,
                   limit = list(c(-0.5, 0.5)))

## perform differential expression
library(org.Hs.eg.db)
library(KEGG.db)

## obtain KEGG id for the Asthma pathway
xx <- as.list(KEGGPATHID2NAME)
hsaID <- unlist(xx)[unlist(xx)=="Asthma"]

## extract the entrez gene ids for members of the Asthma pathway
kegg <- org.Hs.egPATH2EG
entrez <- unlist(as.list(kegg[names(hsaID)]))

## map entrez ids to affymetrix probe set ids
eg <- getEG(rownames(genExp), "hugene10sttranscriptcluster")
geneIDs <- eg[eg %in% entrez]

## extract expression data for genes in the asthma pathway
genesDat <- genExp[names(geneIDs), ]
dim(genesDat)

pval <- apply(genesDat, 1, function(i){
  t.test(i[demo$Time == "post"], i[demo$Time == "pre"], paired = TRUE)$p.value
})
padj <- p.adjust(pval, "BH")
fc <- apply(genesDat, 1, function(i){
 mean(i[demo$Time == "post"]) - mean(i[demo$Time == "pre"])
})

#pdf(paste0(WhereAmI, "5-SupplementaryMaterial/SupplementaryFigure6B.pdf"))
plot(-log10(pval) ~ fc, pch = 19, col = 2, xlim = c(-0.17, 0.155),
     ylab=expression("-log"[10]~~~"p-value"), xlab = expression("log"[2]~~~"fold-change"),
     main = "Asthma KEGG pathway")
points(-log10(pval) ~ fc, pch = 21)
abline(h = -log10(0.05), lty = 2)
text(x = fc, y = -log10(pval), labels = getSYMBOL(names(pval), "hugene10sttranscriptcluster"),
     cex = 0.6, pos = 1:4)
#dev.off()

## using limma
library(limma)
time <- relevel(demo$Time, ref = "pre")
subj <- factor(rep(1:14, 2))
design <- model.matrix(~time+subj)
fit <- eBayes(lmFit(genesDat, design))
top <- topTable(fit, coef = 2, adjust.method = "BH", n = nrow(fit))

pval <- top$P.Value
names(pval) <- rownames(top)
fc <- top$logFC
pdf(paste0(WhereAmI, "5-SupplementaryMaterial/12-SupplementaryFigure6B.pdf"), height = 5.5, width = 5.5)
plot(-log10(pval) ~ fc, pch = 19, col = 2, xlim = c(-0.17, 0.155),
     ylab=expression("-log"[10]~~~"p-value"), xlab = expression("log"[2]~~~"fold-change"),
     main = "Asthma KEGG pathway")
points(-log10(pval) ~ fc, pch = 21)
abline(h = -log10(0.05), lty = 2)
text(x = 0.1, y = -log10(0.045), labels = "p=0.05")
text(x = fc, y = -log10(pval), labels = getSYMBOL(names(pval), "hugene10sttranscriptcluster"),
     cex = 0.6, pos = 1:4)
text(x = -0.1, y = -log10(0.01616927), labels = "BH-FDR=0.46", cex = 0.7, col=2)
dev.off()

##-----------------------------------------------------------------------------
#
# Supplementary Figure 7: Valine, leucine and isoleucine
#
##-----------------------------------------------------------------------------
## how does the valine pathway response in genes and metabolites
grep("Valine", colnames(genMEs0), value = TRUE)
grep("Valine", colnames(metMEs0), value = TRUE)

apply(genMEs0[, grep("Valine", colnames(genMEs0), value = TRUE)], 2, function(i){
  t.test(i[demo$Time == "post"], i[demo$Time == "pre"], paired = TRUE)
})
apply(metMEs0[, grep("Valine", colnames(metMEs0), value = TRUE), drop=FALSE], 2, function(i){
  t.test(i[demo$Time == "post"], i[demo$Time == "pre"], paired = TRUE)
})

## plot pathways
segmentDat <- data.frame(x = rep(1, 3), x1 = rep(2, 3),
                         y = c(2, 3.9, 1.5), y1 = c(2, 3.9, 1.5),
                         Pathway = c("Valine, leucine and isoleucine biosynthesis",
                                     "Valine, leucine and isoleucine degradation",
                                     "Valine, leucine and isoleucine metabolism"))
pvalDat <- data.frame(x = rep(1.5, 3), y = c(1.7, 3.4, 1.2),
                      labels = c("0.003", "0.15", "0.0009"),
                      Pathway = c("Valine, leucine and isoleucine biosynthesis",
                                  "Valine, leucine and isoleucine degradation",
                                  "Valine, leucine and isoleucine metabolism"))

pdf(paste0(WhereAmI, "5-SupplementaryMaterial/13-SupplementaryFigure7a.pdf"), width = 4)
cbind(genMEs0[, grep("Valine", colnames(genMEs0), value = TRUE)], metMEs0[, grep("Valine", colnames(metMEs0), value = TRUE), drop=FALSE]) %>% 
  as.data.frame %>% 
  mutate(Time = relevel(demo$Time, "pre"), Subject = rep(1:14, 2)) %>% 
  gather(Pathway, Exp, -c(Time:Subject)) %>% 
  ggplot(aes(x = Time, y = Exp, fill = Pathway)) + geom_point() +
  geom_line(aes(group = Subject)) + facet_wrap(~Pathway, ncol = 1, scales = "free") +
  customTheme(sizeStripFont=10, xAngle=0, hjust=0.5, vjust=0.5, xSize=10, ySize=10, xAxisSize=10, yAxisSize=10) +
  ylab("Pathway activity") + xlab("Time") + theme(legend.position = "none") +
  geom_segment(data = segmentDat, aes(x=x, y=y, xend = x1, yend = y1)) +
  geom_text(aes(x=x, y=y, label = paste0("p ==", labels), size = 1), 
            data = pvalDat, parse = TRUE)
dev.off()


## correlate gene and metabolite modules
geneMod <- genMEs0[demo$Time == "post", "Valine, leucine and isoleucine biosynthesis"]-genMEs0[demo$Time == "pre", "Valine, leucine and isoleucine biosynthesis"]
metMod <- metMEs0[demo$Time == "post", "Valine, leucine and isoleucine metabolism"]-metMEs0[demo$Time == "pre", "Valine, leucine and isoleucine metabolism"]

plot(geneMod ~ metMod, pch = 19)
abline(lm(geneMod ~ metMod))

## extract the entrez gene ids for members of the Asthma pathway
kegg <- org.Hs.egPATH2EG
entrez <- unlist(as.list(kegg["00290"]))

## map entrez ids to affymetrix probe set ids
eg <- getEG(rownames(genExp), "hugene10sttranscriptcluster")
geneIDs <- eg[eg %in% entrez]

Cov <- data.frame(sample = rep(1:14, 2), stimul = demo$Time)
## extract expression data for genes in the valine, leucine and isoleucine KEGG pathway
bcaaGenes <- t(genExp[names(geneIDs), ])
dim(bcaaGenes)
bcaaGenes.Xw <- withinVariation(X = bcaaGenes, design = Cov)
bcaaGenesLabels <- getSYMBOL(colnames(bcaaGenes), "hugene10sttranscriptcluster")

## metabolomics datasets
bcaaMet <- t(metExp[as.character(metabolites$BIOCHEMICAL[metabolites$SUB_PATHWAY == "Valine, leucine and isoleucine metabolism"]), ])
dim(bcaaMet)
bcaaMet.Xw <- withinVariation(X = bcaaMet, design = Cov)


result <- pls(Y = bcaaGenes.Xw, X = bcaaMet.Xw)
pdf(paste0(WhereAmI, "5-SupplementaryMaterial/14-SupplementaryFigure7b.pdf"), width = 6)
cim(result, col.names = bcaaGenesLabels,
    row.names = colnames(bcaaMet), margins = c(5, 13),
    xlab = "Valine, leucine and isoleucine biosynthesis",
    ylab = "Valine, leucine and isoleucine metabolism")
dev.off()




















###########################################################################
#
# Extras
#
###########################################################################
###################
## pre vs post
rownames(genExp) <- getEG(rownames(genExp), "hugene10sttranscriptcluster")
fc <- rowMeans(genExp[, demo$Time == "post"])-rowMeans(genExp[, demo$Time == "pre"])
geneSymbols <- getSYMBOL(names(fc), "hugene10sttranscriptcluster")
eID <- id2eg(geneSymbols)
names(fc) <- eID[,2]
pv.out <- pathview(gene.data = fc, pathway.id = "05310",
                   species = "hsa", 
                   out.suffix = "PrevsPost", 
                   kegg.native = T,
                   limit = list(c(-0.5, 0.5)))

dat <- as.data.frame(t(genExp[c("7923907", "7944152", "8068254"), ]))
colnames(dat) <- c("IL10", "IL10RA", "IL10RB")
dat$Time <- as.character(demo$Time) 
dat$Time[dat$Time == "pre"] <- "Pre"                     
dat$Time[dat$Time == "post"] <- "Post"    
dat$Time <- factor(dat$Time, levels = c("Pre", "Post")) 
dat$Response <- demo$Response
dat$Group <- rep(1:14, 2)

pdf(paste0(WhereAmI, "Figure 3/IL10.pdf"))
ggplot(dat, aes(x = Time, y = IL10, fill = Response)) + geom_boxplot() +
  geom_line(aes(group = Group, color = Response)) +
  ylab(expression("log"[2]~"expression")) +
  ggtitle("IL-10 (Microarray)") +
  customTheme(sizeStripFont=15, xAngle=0, hjust=0.5, vjust=0.5, 
              xSize=15, ySize=15, xAxisSize=15, yAxisSize=15)
dev.off()

## changes x-axis to Response
ggplot(dat, aes(x = Response, y = IL10, fill = Time)) + geom_boxplot() +
  geom_point(position=position_dodge(width=0.5) ,size = 2.5) +
  geom_line(aes(group = Group, color = Time)) +
  ylab(expression("log"[2]~"expression")) +
  ggtitle("IL-10 (Microarray)") +
  customTheme(sizeStripFont=15, xAngle=0, hjust=0.5, vjust=0.5, 
              xSize=15, ySize=15, xAxisSize=15, yAxisSize=15)


int.dis <- subset(dat, Time == "post")[,1:3]-subset(dat, Time == "pre")[,1:3]
apply(int.dis, 2, function(i){
  t.test(i[dat$Response[dat$Time == "pre"] == "DR"],
         i[dat$Response[dat$Time == "pre"] == "ER"])
})

## test this gene in the RNA-Seq data
## import raw data
load("~/Documents/Asthma/biomarkerPanels/data/discovery/rnaseq/allRnaseqDatasets_rawData.RDATA")
load("~/Documents/Asthma/biomarkerPanels/data/discovery/rnaseq/allRnaseqDatasets_normalized.RDATA")
## normalize by library size
ensemblDat <- normalizelibSum(starEnsemblExp[-(60156:nrow(starEnsemblExp)), colnames(genDats$starEnsemblExp)])
all(rownames(demo) == colnames(ensemblDat))

library(biomaRt)
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", 
               dataset="hsapiens_gene_ensembl")
filterList <- 'hgnc_symbol'
attr = c('ensembl_gene_id', 'hgnc_symbol')

hk.known <- getBM(attributes = attr,
                  filters = filterList,
                  values = c("IL10", "IL10RA", "IL10RB"),
                  mart = mart)


rowNames <- unlist(lapply(strsplit(rownames(ensemblDat), "\\."), function(i) i[1]))
dat2 <- as.data.frame(t(ensemblDat[rowNames %in% hk.known$ensembl_gene_id, ]))
colnames(dat2) <- c("IL10", "IL10RA", "IL10RB")
dat2$Time <- relevel(demo$Time, ref = "Pre") 
dat2$Response <- demo$calculated_Response
dat2$Group <- as.numeric(droplevels(demo$UniqueID))

pdf(paste0(WhereAmI, "Figure 3/IL10_rnaseq.pdf"))
ggplot(dat2, aes(x = Time, y = IL10, fill = Response)) + geom_boxplot() +
  geom_line(aes(group = Group, color = Response)) +
  ylab(expression("log"[2]~"cpm")) +
  ggtitle("IL-10 (RNA-Seq)") +
  customTheme(sizeStripFont=15, xAngle=0, hjust=0.5, vjust=0.5, 
              xSize=15, ySize=15, xAxisSize=15, yAxisSize=15)
dev.off()


# ER
genExpER <- genExp[, demo$Response == "ER"]
time <- demo$Time[demo$Response == "ER"]
fc <- rowMeans(genExpER[, time == "post"])-rowMeans(genExpER[, time == "pre"])
geneSymbols <- getSYMBOL(names(fc), "hugene10sttranscriptcluster")
eID <- id2eg(geneSymbols)
names(fc) <- eID[,2]
pv.out <- pathview(gene.data = fc, pathway.id = "05310",
                   species = "hsa", 
                   out.suffix = "ERprevspost", 
                   kegg.native = T,
                   limit = list(c(-0.5, 0.5)))
genExpER.asthma <- genExpER[names(fc) %in% unlist(strsplit(pv.out$plot.data.gene$all.mapped, ",")), ]
erGenSym <- apply(genExpER.asthma, 1, function(i){
  t.test(i[time == "post"], i[time == "pre"], paired = TRUE)$p.value
})
names(erGenSym) <- as.character(getSYMBOL(names(erGenSym), "hugene10sttranscriptcluster"))

# DR
genExpDR <- genExp[, demo$Response == "DR"]
time <- demo$Time[demo$Response == "DR"]
fc <- rowMeans(genExpDR[, time == "post"])-rowMeans(genExpDR[, time == "pre"])
geneSymbols <- getSYMBOL(names(fc), "hugene10sttranscriptcluster")
eID <- id2eg(geneSymbols)
names(fc) <- eID[,2]
pv.out <- pathview(gene.data = fc, pathway.id = "05310",
                   species = "hsa", 
                   out.suffix = "DRprevspost", 
                   kegg.native = T,
                   limit = list(c(-0.5, 0.5)))
genExpDR.asthma <- genExpDR[names(fc) %in% unlist(strsplit(pv.out$plot.data.gene$all.mapped, ",")), ]
drGenSym <- apply(genExpDR.asthma, 1, function(i){
  t.test(i[time == "post"], i[time == "pre"], paired = TRUE)$p.value
})
names(drGenSym) <- as.character(getSYMBOL(names(drGenSym), "hugene10sttranscriptcluster"))


## interaction
diff <- genExp[, demo$Time == "post"]-genExp[, demo$Time == "pre"]
res <- demo$Response[demo$Time == "post"]
genExp.asthma <- diff[names(fc) %in% unlist(strsplit(pv.out$plot.data.gene$all.mapped, ",")), ]
genSym <- apply(genExp.asthma, 1, function(i){
  t.test(i[res == "DR"], i[res == "ER"])$p.value
})
names(genSym) <- as.character(getSYMBOL(names(genSym), "hugene10sttranscriptcluster"))
p.adjust(genSym, "BH")

## pre vs post
prepost <- genExp[names(fc) %in% unlist(strsplit(pv.out$plot.data.gene$all.mapped, ",")), ]
genSym <- apply(prepost, 1, function(i){
  t.test(i[demo$Time == "post"], i[demo$Time == "pre"], paired = TRUE)$p.value
})
names(genSym) <- as.character(getSYMBOL(names(genSym), "hugene10sttranscriptcluster"))
p.adjust(genSym, "BH")

