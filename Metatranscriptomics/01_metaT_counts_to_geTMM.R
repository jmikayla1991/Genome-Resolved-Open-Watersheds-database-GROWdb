#input metaT mapping featureCounts output
#be sure to remove the first line on server before transferring using sed -i '1d' <name of file>
x=read.delim("INPUT.out")
#set row names
row.names(x)=x$Geneid
x=x[,-1]
#remove unnecessary columns
x=x[,-c(1:4)]
#convert to rpk (length normalization= divide each count of gene by gene length in column 1)
rpk <- (x[,2:ncol(x)]/x[,1])
#TMM normalization
library(edgeR)
rpk.norm <- DGEList(counts=rpk)
rpk.norm <- calcNormFactors(rpk.norm)
norm.counts.rpk_edger <- cpm(rpk.norm)
write.csv(norm.counts.rpk_edger,file="OUTPUT.csv")

