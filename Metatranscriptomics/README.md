# GROWdb Metatranscriptomics Workflow

## Samples

### Surface water

**1. WHONDRS** - there are 61 metaT total from WHONDRS surface water S19S, of those 57 are USA based and will be in the USA focused surface water story
	- see ListOfWHONDRS_sw_metaT.xlsx for list and locations  


**2. WROL** - coming soon

### Sediment

**1. WHONDRS** - coming soon

## USA based metatranscriptomic mapping and analysis

### Trimming and filtering 

Trim adapters and quality trim

	bbduk.sh in=<sampleid>.fastq out=<sampleid>_trimmed.fastq ktrim=r k=23 Mink=11 hdist=1 qtrim=rl trimq=20 minlen=75 Maq=10 t=7 ref=/opt/bbtools/bbmap/resources/adapters.fa 

Zip reads for rqcfilter

	gzip -c <sampleid>_trimmed.fastq > <sampleid>_trimmed.fastq.gz 
	
Filter reads using rqcfilter 

	rqcfilter2.sh jni=t in=<sampleid>_trimmed.fastq.gz path=/home/projects-wrighton-2/GROWdb/metatranscriptomes/trimmed_reads/<sampleid>_filter  rqcfilterdata=/home/opt/RQCFilterData rna=t trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t removemicrobes=t sketch kapa=t clumpify=t barcodefilter=f trimpolyg=5 usejni=f

### Mapping 

Unzip reads into mapping dir

	gunzip -c /home/projects-wrighton-2/GROWdb/metatranscriptomes/trimmed_reads/<sampleid>_filter/<sampleid>_trimmed.anqrpht.fastq.gz > <sampleid>_trimmed.anqrpht.fastq 

Deinterleave trimmed and filtered fastqs
	
	/ORG-Data/scripts/deinterleave_fastq.sh < <sampleid>_trimmed.anqrpht.fastq <sampleid>_trimmed.anqrpht_R1.fastq <sampleid>_trimmed.anqrpht_R2.fastq

Mapping with bowtie
	
	bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x /home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaG_mapping/mapping_bowtie_012722/2093_dRep99_MAG_scaffolds_DB -S <sampleid>_bowtie.sam -1 <sampleid>_trimmed.anqrpht_R1.fastq -2 <sampleid>_trimmed.anqrpht_R2.fastq 

Sam to bam file
	
	samtools view -@ 30 -bS <sampleid>_bowtie.sam > <sampleid>_bowtie.bam 

Reformat bam file to retain mapped reads at 97% id
	
	reformat.sh in=<sampleid>_bowtie.sorted.bam out=<sampleid>_bowtie.reformat97.bam minidfilter=0.97 primaryonly=t pairedonly=f 

Name sort bam files
	
	samtools sort -n -o <sampleid>_bowtie.reformat97_sorted.bam <sampleid>_bowtie.reformat97.bam -@ 30 

### Counting and filtering
