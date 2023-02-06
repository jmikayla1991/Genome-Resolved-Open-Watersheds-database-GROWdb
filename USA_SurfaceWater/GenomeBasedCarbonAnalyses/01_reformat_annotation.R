library(dplyr)
library(tidyr)

annotations=read.delim("41827_metaT_genes_annotations.tsv",sep="\t",header=TRUE)
geTMM=read.delim("metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd_geTMM.csv",sep=",",header=TRUE)

gene_bam=geTMM %>% gather(-X,key="fasta",value="geTMM")%>%filter(geTMM>0)%>% select(X,fasta)
  
annotations_2=left_join(gene_bam,annotations,by="X")%>%select(-fasta.y)%>%rename(fasta=fasta.x)
write.table(annotations_2,file="41827genes_58bams_annotations.tsv",sep="\t",row.names=FALSE)
