# Mapping for GROWdb surface water USA

### Samples- there are 158 samples with metagenomes 

	List_of_samples.txt

### Mapping 

Build bowtie database

	bowtie2-build 2093_dRep99_MAG_scaffolds.fa 2093_dRep99_MAG_scaffolds_DB --threads 7 
	
Output database

	/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaG_mapping/mapping_bowtie_012722/2093_dRep99_MAG_scaffolds_DB

Mapping all 158 metaG (sickle trimmed read sets) with bowtie to the scaffolds derived from DRAM annotations of 2093 MAGs (532750 scaffolds) 

	bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x 2093_dRep99_MAG_scaffolds_DB -S <sample name>_bowtie.sam -1 <sample name>_R1_trimmed.fastq -2 <sample name>_2_trimmed.fastq

