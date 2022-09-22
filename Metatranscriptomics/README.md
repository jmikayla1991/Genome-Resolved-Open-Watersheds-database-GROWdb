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

Mapping all 61 metaT with bowtie to the gff derived from DRAM annotations of 
	
	bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 30 -x /home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaG_mapping/mapping_bowtie_012722/2093_dRep99_MAG_scaffolds_DB -S <sampleid>_bowtie.sam -1 <sampleid>_trimmed.anqrpht_R1.fastq -2 <sampleid>_trimmed.anqrpht_R2.fastq 

Sam to bam file
	
	samtools view -@ 30 -bS <sampleid>_bowtie.sam > <sampleid>_bowtie.bam 

Reformat bam file to retain mapped reads at 97% id
	
	reformat.sh in=<sampleid>_bowtie.sorted.bam out=<sampleid>_bowtie.reformat97.bam minidfilter=0.97 primaryonly=t pairedonly=f 

Name sort bam files
	
	samtools sort -n -o <sampleid>_bowtie.reformat97_sorted.bam <sampleid>_bowtie.reformat97.bam -@ 30 

### Counting

Feature counts 

	featureCounts -T 20 -t CDS -g ID -s 2 -p -a /home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/all_bins/dRep_v2.6.2/dereplicated_genomes/merged_annotations_drep_bins/genes_no_rna_fix.gff -o metaT_output_feature_counts_revstranded_070722_97perc.out *_MT_bowtie.reformat97_sorted.bam

Server location of feature counts output metaT_output_feature_counts_revstranded_070722_97perc.out (3991706 lines) 

	/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/USA/metaT_output_feature_counts_revstranded_07192022_97perc.out

	
### Filtered and transformed data tables 

**1. USA only metaT gene table transformed to geTMM (not filtered to genes in 3 or more samples)**  

Remove first line of file before putting into R 

	sed -i '1d' metaT_output_feature_counts_revstranded_07192022_97perc.out 

Now metaT_output_feature_counts_revstranded_07192022_97perc.out is 3991705 lines

Run R script 01_metaT_counts_to_geTMM.R to obtain geTMMs

Output

	/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/USA/geTMM_tables/geTMM_norm.counts.rpk_edger_genes_nofilter_092122.csv

**2. USA only metaT gene table transformed to geTMM filtered to genes in 3 or more samples** 

Filter feature count output to genes that are in 3 or more samples 

	with open("./metaT_output_feature_counts_revstranded_07192022_97perc.out", 'r') as counts:
    with open("./metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd.csv", 'w') as out_data:
        for i, l in enumerate(counts):
            if  i < 2:
                _ = out_data.write(l)
                continue
            if len([j for j in l.split('\t')[6:] if float(j)> 0]) > 5:
                ln = out_data.write(l) 


Server location of filtered feature counts output metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd.csv (41829 lines)

	/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/USA/metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd.csv

Remove first line of file before putting into R 

	sed -i '1d' metaT_output_feature_counts_revstranded_07192022_97perc_sample_filterd.csv 

Run R script 01_metaT_counts_to_geTMM.R to obtain geTMMs



