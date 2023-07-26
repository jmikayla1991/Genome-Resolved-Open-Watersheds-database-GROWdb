#input metaT mapping featureCounts output
#be sure to remove the first line on server before transferring using sed -i '1d' <name of file>
x=read.csv("metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd.csv")
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
write.csv(norm.counts.rpk_edger,file="geTMM_norm.counts.rpk_edger_genes.csv")

#join in annotations
anno=read.delim("annotations.tsv")
row.names(anno)=anno[,1]
anno=anno[,-1]
norm.counts.rpk_edger_withanno=merge(norm.counts.rpk_edger,anno, by = 'row.names')
write.csv(norm.counts.rpk_edger_withanno,file="geTMM_norm.counts.rpk_edger_genes_withanno.csv")

#now need to collapse by genome
library(dplyr)
norm.counts.rpk_edger.zerosremoved_for_bins <- as.data.frame(norm.counts.rpk_edger)
norm.counts.rpk_edger.zerosremoved_for_bins <- tibble::rownames_to_column(norm.counts.rpk_edger.zerosremoved_for_bins, "rn")
norm.counts.rpk_edger.zerosremoved_for_bins$rn<-gsub('_Ga(.*)','',norm.counts.rpk_edger.zerosremoved_for_bins$rn)
norm.counts.rpk_edger.zerosremoved_for_bins$rn<-gsub('_k121(.*)','',norm.counts.rpk_edger.zerosremoved_for_bins$rn)
norm.counts.rpk_edger.zerosremoved_for_bins$rn<-gsub('_scaffold(.*)','',norm.counts.rpk_edger.zerosremoved_for_bins$rn)
#collapse rows by sum and median
norm.counts.rpk_edger.bins_mean<-norm.counts.rpk_edger.zerosremoved_for_bins%>%group_by(rn)%>%summarise_all(funs(mean))
norm.counts.rpk_edger.bins_sum<-norm.counts.rpk_edger.zerosremoved_for_bins%>%group_by(rn)%>%summarise_all(funs(sum))
norm.counts.rpk_edger.bins_median<-norm.counts.rpk_edger.zerosremoved_for_bins%>%group_by(rn)%>%summarise_all(funs(median))
write.csv(norm.counts.rpk_edger.bins_mean,"norm.counts.rpk_edger.bins_mean.csv")
#add the adjectives 
adjectives=read.delim("adjectives.tsv")
merged_data_adjectives_median_bins<-inner_join(adjectives,norm.counts.rpk_edger.bins_median,by='rn')
write.csv(merged_data_adjectives_median_bins,"merged_data_adjectives_median_bins.csv")




