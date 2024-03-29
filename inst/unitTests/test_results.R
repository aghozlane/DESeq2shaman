test_contrasts <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200,m=12)
  dds$condition <- factor(rep(1:3,each=4))
  dds$group <- factor(rep(1:2,length=ncol(dds)))
  counts(dds)[1,] <- rep(c(100L,200L,800L),each=4)
  
  design(dds) <- ~ group + condition
  dds <- DESeq(dds)
  
  # check to see if the contrasts with expanded model matrix
  # are close to expected (although shrunk due to the beta prior)
  lfc31 <- results(dds,contrast=c("condition","3","1"))[1,2]
  lfc21 <- results(dds,contrast=c("condition","2","1"))[1,2]
  lfc32 <- results(dds,contrast=c("condition","3","2"))[1,2]
  checkEqualsNumeric(lfc31, 3, tolerance=.1)
  checkEqualsNumeric(lfc21, 1, tolerance=.1)
  checkEqualsNumeric(lfc32, 2, tolerance=.1)
  checkEqualsNumeric(results(dds,contrast=c("condition","1","3"))[1,2], -3, tolerance=.1)
  checkEqualsNumeric(results(dds,contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  checkEqualsNumeric(results(dds,contrast=c("condition","2","3"))[1,2], -2, tolerance=.1)

  # check that results are not changed by releveling
  dds2 <- dds
  colData(dds2)$condition <- relevel(colData(dds2)$condition, "2")
  dds2 <- DESeq(dds2)
  checkEqualsNumeric(results(dds2,contrast=c("condition","3","1"))[1,2], lfc31, tolerance=1e-6)
  checkEqualsNumeric(results(dds2,contrast=c("condition","2","1"))[1,2], lfc21, tolerance=1e-6)
  checkEqualsNumeric(results(dds2,contrast=c("condition","3","2"))[1,2], lfc32, tolerance=1e-6)

  # check the default prior  variance on the intercept and group LFC's
  dds3 <- dds
  design(dds3) <- ~ group + condition + condition:group
  dds3 <- nbinomWaldTest(dds3)
  checkEquals(attr(dds3,"betaPriorVar")[1:6],
              c(Intercept=1e6, group1=1e3, group2=1e3,
                condition1=1e3, condition2=1e3, condition3=1e3))

  # test a number of contrast as list options
  checkEqualsNumeric(results(dds, contrast=list("condition3","condition1"))[1,2], lfc31, tolerance=1e-6)
  results(dds, contrast=list("condition3","condition1"), listValues=c(.5,-.5))
  results(dds, contrast=list("condition3",character()))
  results(dds, contrast=list("condition3",character()), listValues=c(.5,-.5))
  results(dds, contrast=list(character(),"condition1"))
  results(dds, contrast=list(character(),"condition1"), listValues=c(.5,-.5))

  # test no prior on intercept
  checkEqualsNumeric(attr(dds,"betaPriorVar")[1], 1e6)  
}


test_designs_zero_intercept <- function() {
  # test some special cases for results()
  # using designs with a + 0 
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200,m=12)
  dds$condition <- factor(rep(1:3,each=4))
  dds$group <- factor(rep(1:2,length=ncol(dds)))

  counts(dds)[1,] <- rep(c(100L,200L,400L),each=4)
  
  design(dds) <- ~ condition + 0
  dds <- DESeq(dds, betaPrior=FALSE)

  checkEqualsNumeric(results(dds)[1,2], 2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","1"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","3","2"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","3"))[1,2], -2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","3"))[1,2], -1, tolerance=.1)
  
  design(dds) <- ~ group + condition + 0
  dds <- DESeq(dds, betaPrior=FALSE)

  checkEqualsNumeric(results(dds)[1,2], 2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","1"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","3","2"))[1,2], 1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","3"))[1,2], -2, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","1","2"))[1,2], -1, tolerance=.1)
  checkEqualsNumeric(results(dds, contrast=c("condition","2","3"))[1,2], -1, tolerance=.1)
}

test_LRT_then_Wald <- function() {
  set.seed(1)
  dds <- makeExampleDESeqDataSet(n=200)
  dds$group <- factor(rep(1:2,6))
  design(dds) <- ~ group + condition
  dds <- DESeq(dds, test="LRT", reduced=~group)
  checkTrue(!all(results(dds,name="condition_B_vs_A")$stat ==
                 results(dds,name="condition_B_vs_A",test="Wald")$stat))
}
