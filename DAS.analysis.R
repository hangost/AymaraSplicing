# Analysis for differential alternative splicing


sig.cut <- function(test.mat){
  test.mat <- test.mat[as.double(test.mat[,"FDR"]) < cut.fdr & abs(as.double(test.mat[,"IncLevelDifference"])) > cut.fc,]
  return (test.mat)
}


tarMet <- "JunctionCountOnly"
GNC.A3SS.f <- paste("./5.rMATs.result/GNC/MATS_output/A3SS.MATS.",tarMet,".txt",sep="")
GNC.A5SS.f <- paste("./5.rMATs.result/GNC/MATS_output/A5SS.MATS.",tarMet,".txt",sep="")
GNC.ES.f <- paste("./5.rMATs.result/GNC/MATS_output/SE.MATS.",tarMet,".txt",sep="")
GNC.IR.f <- paste("./5.rMATs.result/GNC/MATS_output/RI.MATS.",tarMet,".txt",sep="")
GNC.MXE.f <- paste("./5.rMATs.result/GNC/MATS_output/MXE.MATS.",tarMet,".txt",sep="")

GNC.A3SS.re <- gsub(" ","",as.matrix(read.table(GNC.A3SS.f,sep='\t',header=T)))
GNC.A5SS.re <- gsub(" ","",as.matrix(read.table(GNC.A5SS.f,sep='\t',header=T)))
GNC.ES.re <- gsub(" ","",as.matrix(read.table(GNC.ES.f,sep='\t',header=T)))
GNC.IR.re <- gsub(" ","",as.matrix(read.table(GNC.IR.f,sep='\t',header=T)))
GNC.MXE.re <- gsub(" ","",as.matrix(read.table(GNC.MXE.f,sep='\t',header=T)))

cut.fdr <- 0.05
cut.fc <- 0.1

sig.GNC.A3SS.re <- sig.cut(GNC.A3SS.re)
sig.GNC.A5SS.re <- sig.cut(GNC.A5SS.re)
sig.GNC.ES.re <- sig.cut(GNC.ES.re)
sig.GNC.IR.re <- sig.cut(GNC.IR.re)
sig.GNC.MXE.re <- sig.cut(GNC.MXE.re)

a3ss <- sig.GNC.A3SS.re[,c("geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE",
	"flankingES","flankingEE","PValue","FDR","IncLevelDifference")]
a5ss <- sig.GNC.A5SS.re[,c("geneSymbol","chr","strand","longExonStart_0base","longExonEnd","shortES","shortEE","flankingES",
	"flankingEE","PValue","FDR","IncLevelDifference")]
es <- sig.GNC.ES.re[,c("geneSymbol","chr","strand","exonStart_0base","exonEnd","upstreamES","upstreamEE","downstreamES",
	"downstreamEE","PValue","FDR","IncLevelDifference")]
ir <- sig.GNC.IR.re[,c("geneSymbol","chr","strand","riExonStart_0base","riExonEnd","upstreamES","upstreamEE","downstreamES",
	"downstreamEE","PValue","FDR","IncLevelDifference")]
mxe <- sig.GNC.MXE.re[,c("geneSymbol","chr","strand","X1stExonStart_0base","X1stExonEnd","X2ndExonStart_0base",
	"X2ndExonEnd","upstreamES","upstreamEE","downstreamES","downstreamEE","PValue","FDR","IncLevelDifference")]

write.table(a3ss,"./results/A3SS.re",sep='\t',quote=F,row.names=F)
write.table(a5ss,"./results/A5SS.re",sep='\t',quote=F,row.names=F)
write.table(es,"./results/ES.re",sep='\t',quote=F,row.names=F)
write.table(ir,"./results/IR.re",sep='\t',quote=F,row.names=F)
write.table(mxe,"./results/MXE.re",sep='\t',quote=F,row.names=F)

