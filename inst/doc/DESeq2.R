## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,dev="png",fig.show="hide",
               fig.width=4,fig.height=4.5,
               message=FALSE)

## ----style, eval=TRUE, echo=FALSE, results="asis"---------------------------------------
BiocStyle::latex()

## ----loadDESeq2, echo=FALSE-------------------------------------------------------------
library("DESeq2")

## ----options, results="hide", echo=FALSE--------------------------------------
options(digits=3, width=80, prompt=" ", continue=" ")

## ----quick, eval=FALSE--------------------------------------------------------
#  dds <- DESeqDataSet(se = se, design = ~ condition)
#  dds <- DESeq(dds)
#  res <- results(dds)

## ----loadSumExp,cache=TRUE----------------------------------------------------
library("airway")
data("airway")
se <- airway

## ----sumExpInput, cache=TRUE--------------------------------------------------
library("DESeq2")
ddsSE <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE

## ----loadPasilla--------------------------------------------------------------
library("pasilla")
library("Biobase")
data("pasillaGenes")
countData <- counts(pasillaGenes)
colData <- pData(pasillaGenes)[,c("condition","type")]

## ----matrixInput--------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
dds

## ----addFeatureData-----------------------------------------------------------
featureData <- data.frame(gene=rownames(pasillaGenes))
(mcols(dds) <- DataFrame(mcols(dds), featureData))

## ----htseqDirI, eval=FALSE----------------------------------------------------
#  directory <- "/path/to/your/files/"

## ----htseqDirII---------------------------------------------------------------
directory <- system.file("extdata", package="pasilla", mustWork=TRUE)

## ----htseqInput, cache=TRUE---------------------------------------------------
sampleFiles <- grep("treated",list.files(directory),value=TRUE)
sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

## ----relevel------------------------------------------------------------------
dds$condition <- relevel(dds$condition, "untreated")

## ----droplevels---------------------------------------------------------------
dds$condition <- droplevels(dds$condition)

## ----deseq, cache=TRUE--------------------------------------------------------
dds <- DESeq(dds)
res <- results(dds)

## ----parallel, eval=FALSE-----------------------------------------------------
#  library("BiocParallel")
#  register(MulticoreParam(4))

## ----resOrder-----------------------------------------------------------------
resOrdered <- res[order(res$padj),]
head(resOrdered)

## ----sumRes-------------------------------------------------------------------
summary(res)

## ----MA, fig.width=4.5, fig.height=4.5----------------------------------------
plotMA(res, main="DESeq2", ylim=c(-2,2))

## ----resMLE-------------------------------------------------------------------
resMLE <- results(dds, addMLE=TRUE)
head(resMLE, 4)

## ----MANoPrior, echo=FALSE, fig.width=4.5, fig.height=4.5---------------------
df <- data.frame(resMLE$baseMean, resMLE$lfcMLE, ifelse(is.na(res$padj), FALSE, res$padj < .1))
plotMA(df, main=expression(unshrunken~log[2]~fold~changes), ylim=c(-2,2))

## ----MAidentify, eval=FALSE---------------------------------------------------
#  identify(res$baseMean, res$log2FoldChange)

## ----plotCounts, dev="pdf", fig.width=4.5, fig.height=5-----------------------
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

## ----plotCountsAdv, dev="pdf", fig.width=3.5, fig.height=3.5------------------
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

## ----metadata-----------------------------------------------------------------
mcols(res)$description

## ----addMLE-------------------------------------------------------------------
head(results(dds, addMLE=TRUE),4)

## ----export, eval=FALSE-------------------------------------------------------
#  write.csv(as.data.frame(resOrdered),
#            file="condition_treated_results.csv")

## ----subset-------------------------------------------------------------------
resSig <- subset(resOrdered, padj < 0.1)
resSig

## ----multifactor--------------------------------------------------------------
colData(dds)

## ----copyMultifactor----------------------------------------------------------
ddsMF <- dds

## ----replaceDesign,cache=TRUE-------------------------------------------------
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)

## ----multiResults-------------------------------------------------------------
resMF <- results(ddsMF)
head(resMF)

## ----multiTypeResults---------------------------------------------------------
resMFType <- results(ddsMF, contrast=c("type","single-read","paired-end"))
head(resMFType)

## ----rlogAndVST---------------------------------------------------------------
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
rlogMat <- assay(rld)
vstMat <- assay(vsd)

## ----vsd1, echo=FALSE, fig.width=4.5, fig.height=4.5--------------------------
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord],
        cbind(assay(vsd)[, 1], log2(px))[ord, ],
        type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright",
       legend = c(
        expression("variance stabilizing transformation"),
        expression(log[2](n/s[1]))),
       fill=vstcol)

## ----vsd2, fig.width=8, fig.height=3------------------------------------------
library("vsn")
par(mfrow=c(1,3))
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))

## ----heatmap------------------------------------------------------------------
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

## ----figHeatmap2a, dev="pdf", fig.width=7, fig.height=10----------------------
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))

## ----figHeatmap2b, dev="pdf", fig.width=7, fig.height=10----------------------
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

## ----figHeatmap2c, dev="pdf", fig.width=7, fig.height=10----------------------
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))

## ----sampleClust--------------------------------------------------------------
distsRL <- dist(t(assay(rld)))

## ----figHeatmapSamples, dev="pdf", fig.width=7, fig.height=7------------------
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition, type, sep=" : "))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none", 
          col = rev(hmcol), margin=c(13, 13))

## ----figPCA, dev="pdf", fig.width=5, fig.height=3-----------------------------
plotPCA(rld, intgroup=c("condition", "type"))

## ----figPCA2, dev="pdf", fig.width=5, fig.height=3----------------------------
data <- plotPCA(rld, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

## ----WaldTest, eval=FALSE-----------------------------------------------------
#  dds <- estimateSizeFactors(dds)
#  dds <- estimateDispersions(dds)
#  dds <- nbinomWaldTest(dds)

## ----dispFit------------------------------------------------------------------
plotDispEsts(dds)

## ----dispFitCustom------------------------------------------------------------
ddsCustom <- dds
useForMedian <- mcols(ddsCustom)$dispGeneEst > 1e-7
medianDisp <- median(mcols(ddsCustom)$dispGeneEst[useForMedian],na.rm=TRUE)
dispersionFunction(ddsCustom) <- function(mu) medianDisp
ddsCustom <- estimateDispersionsMAP(ddsCustom)

## ----filtByMean, dev="pdf"----------------------------------------------------
attr(res,"filterThreshold")
plot(attr(res,"filterNumRej"),type="b",
     ylab="number of rejections")

## ----noFilt-------------------------------------------------------------------
resNoFilt <- results(dds, independentFiltering=FALSE)
addmargins(table(filtering=(res$padj < .1), noFiltering=(resNoFilt$padj < .1)))

## ----ddsNoPrior---------------------------------------------------------------
ddsNoPrior <- DESeq(dds, betaPrior=FALSE)

## ----lfcThresh----------------------------------------------------------------
par(mfrow=c(2,2),mar=c(2,2,1,1))
yl <- c(-2.5,2.5)

resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(ddsNoPrior, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")

plotMA(resGA, ylim=yl)
abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resLA, ylim=yl)
abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resG, ylim=yl)
abline(h=.5,col="dodgerblue",lwd=2)
plotMA(resL, ylim=yl)
abline(h=-.5,col="dodgerblue",lwd=2)

## ----mcols--------------------------------------------------------------------
mcols(dds,use.names=TRUE)[1:4,1:4]
# here using substr() only for display purposes
substr(names(mcols(dds)),1,10) 
mcols(mcols(dds), use.names=TRUE)[1:4,]

## ----muAndCooks---------------------------------------------------------------
head(assays(dds)[["mu"]])
head(assays(dds)[["cooks"]])

## ----dispersions--------------------------------------------------------------
head(dispersions(dds))
# which is the same as 
head(mcols(dds)$dispersion)

## ----sizefactors--------------------------------------------------------------
sizeFactors(dds)

## ----coef---------------------------------------------------------------------
head(coef(dds))

## ----betaPriorVar-------------------------------------------------------------
attr(dds, "betaPriorVar")

## ----dispPriorVar-------------------------------------------------------------
dispersionFunction(dds)
attr(dispersionFunction(dds), "dispPriorVar")

## ----normFactors, eval=FALSE--------------------------------------------------
#  normFactors <- normFactors / exp(rowMeans(log(normFactors)))
#  normalizationFactors(dds) <- normFactors

## ----offsetTransform, eval=FALSE----------------------------------------------
#  cqnOffset <- cqnObject$glm.offset
#  cqnNormFactors <- exp(cqnOffset)
#  EDASeqNormFactors <- exp(-1 * EDASeqOffset)

## ----cooksPlot----------------------------------------------------------------
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic", 
     ylab="maximum Cook's distance per gene",
     ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

## ----indFilt------------------------------------------------------------------
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

## ----histindepfilt, dev="pdf", fig.width=7, fig.height=5----------------------
use <- res$baseMean > attr(res,"filterThreshold")
table(use)
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

## ----fighistindepfilt---------------------------------------------------------
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

## ----sortP, cache=TRUE--------------------------------------------------------
resFilt <- res[use & !is.na(res$pvalue),]
orderInPlot <- order(resFilt$pvalue)
showInPlot <- (resFilt$pvalue[orderInPlot] <= 0.08)
alpha <- 0.1

## ----sortedP, fig.width=4.5, fig.height=4.5-----------------------------------
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)

## ----doBH, echo=FALSE, results="hide"-----------------------------------------
whichBH <- which(resFilt$pvalue[orderInPlot] <= alpha * seq(along=resFilt$pvalue)/length(resFilt$pvalue))
whichBH <- seq_len(max(whichBH))
## Test some assertions:
## - whichBH is a contiguous set of integers from 1 to length(whichBH)
## - the genes selected by this graphical method coincide with those
##   from p.adjust (i.e. padjFilt)
stopifnot(length(whichBH)>0,
          identical(whichBH, seq(along=whichBH)),
          resFilt$padj[orderInPlot][ whichBH] <= alpha,
          resFilt$padj[orderInPlot][-whichBH]  > alpha)

## ----SchwSpjot, echo=FALSE, results="hide"------------------------------------
j  <- round(length(resFilt$pvalue)*c(1, .66))
px <- (1-resFilt$pvalue[orderInPlot[j]])
py <- ((length(resFilt$pvalue)-1):0)[j]
slope <- diff(py)/diff(px)

## ----SchwederSpjotvoll, fig.width=4.5, fig.height=4.5-------------------------
plot(1-resFilt$pvalue[orderInPlot],
     (length(resFilt$pvalue)-1):0, pch=".",
     xlab=expression(1-p[i]), ylab=expression(N(p[i])))
abline(a=0, slope, col="red3", lwd=2)

## ----vanillaDESeq, eval=FALSE-------------------------------------------------
#  dds <- DESeq(dds, minReplicatesForReplace=Inf)
#  res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)

## ----sessInfo, results="asis", echo=FALSE-------------------------------------
toLatex(sessionInfo())

## ----resetOptions, results="hide", echo=FALSE---------------------------------
options(prompt="> ", continue="+ ")

