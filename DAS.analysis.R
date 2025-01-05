# Analysis for differential alternative splicing

sig.cut <- function(test.mat){
  test.mat <- test.mat[as.double(test.mat[,"FDR"]) < cut.fdr & abs(as.double(test.mat[,"IncLevelDifference"])) > cut.fc,]
  sample1.r.IC <- do.call(rbind,strsplit(test.mat[,is.element(colnames(test.mat),c("IC_SAMPLE_1","IJC_SAMPLE_1"))],","))
  sample1.r.SC <- do.call(rbind,strsplit(test.mat[,is.element(colnames(test.mat),c("SC_SAMPLE_1","SJC_SAMPLE_1"))],","))
  sample2.r.IC <- do.call(rbind,strsplit(test.mat[,is.element(colnames(test.mat),c("IC_SAMPLE_2","IJC_SAMPLE_2"))],","))
  sample2.r.SC <- do.call(rbind,strsplit(test.mat[,is.element(colnames(test.mat),c("SC_SAMPLE_2","SJC_SAMPLE_2"))],","))
  group1.exp <- do.call(rbind,strsplit(test.mat[,"IncLevel1"],","))
  group2.exp <- do.call(rbind,strsplit(test.mat[,"IncLevel2"],","))
  group1.sd <- lapply(1:length(group1.exp[,1]),function(ea.nm){
    sd(as.double(group1.exp[ea.nm,as.double(sample1.r.IC[ea.nm,]) + as.double(sample1.r.SC[ea.nm,]) >= 4]))
  })
  group2.sd <- lapply(1:length(group2.exp[,1]),function(ea.nm){
    sd(as.double(group2.exp[ea.nm,as.double(sample2.r.IC[ea.nm,]) + as.double(sample2.r.SC[ea.nm,]) >= 4]))
  })
  group1.sd <- unlist(group1.sd)
  group2.sd <- unlist(group2.sd)
  #test.mat <- test.mat[group1.sd < 0.1 & group2.sd < 0.1,]
  return (test.mat)
}
