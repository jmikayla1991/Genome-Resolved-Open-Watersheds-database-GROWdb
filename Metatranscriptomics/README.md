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

### Counting and filtering
