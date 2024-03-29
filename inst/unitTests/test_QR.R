test_QR <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=100,betaSD=1)
  dds <- DESeq(dds)
  ddsNoQR <- nbinomWaldTest(dds,useQR=FALSE)
  res <- results(dds)
  resNoQR <- results(ddsNoQR)
  checkEqualsNumeric(res$log2FoldChange, resNoQR$log2FoldChange, tolerance=1e-6)
}
