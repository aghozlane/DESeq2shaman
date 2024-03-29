test_interactions <- function() {
  dds <- makeExampleDESeqDataSet(n=100,m=8)
  colData(dds)$group <- factor(rep(c("x","y"),times=ncol(dds)/2))
  design(dds) <- ~ condition + condition:group
  dds <- DESeq(dds)
  design(dds) <- ~ condition + group + condition:group
  dds <- DESeq(dds)
}
