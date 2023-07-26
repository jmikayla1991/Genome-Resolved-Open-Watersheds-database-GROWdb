# GROWdb Metagenomic Pipeline for generating metagenome assembled genomes (MAGs) from raw reads

![workflow](https://github.com/jmikayla1991/Genome-Resolved-Open-Watersheds-database-GROWdb/blob/main/USA_SurfaceWater/Metagenomic_Pipeline_GROWdb/Analysis_flowchart_v3.jpg)

For each set of metagenomic reads in the GROW database, 4 types of assemblies will be performed:

1. Reads trimmed with sickle (v1.33), assembled with megahit (v1.2.9), and binned with metabat2 (2.12.1)
2. Reads trimmed with sickle (v1.33), filtered randomly to 25%, assembled with idba-ud (1.1.0), and binned with metabat2 (2.12.1)
3. Reads trimmed and filtered with rqc2filter, filtered using bbcms, assembled with metaspades (v3.13.0), and binned with metabat2 (2.12.1). ** Note this is the exact JGI pipeline. For samples sequenced at JGI, MAGs were downloaded directly from from JGI. rqc2filter and bbcms are under bbtools version 38.89.

### Commands for each pathway can be found below. 

Pathway 1: 

Read trimming with sickle:

sickle pe -f <forward untrimmed reads>  -r <reverse untrimmed reads>   -t <sequencing platform> -o <sample_name>_R1_trimmed.fastq -p <sample_name>_R2_trimmed.fastq -s discared_R1R2.fastq

Assembly with megahit:

megahit -1 <sample_name>_R1_trimmed.fastq -2 <sample_name>_R2_trimmed.fastq --k-min 31 --k-max 121 --k-step 10 -m <max memory in byte, for W2 cluster, generally do 0.4, 40% of machines total mem>  -t <number of threads>

Binning with metabat2:

pullseq.py -i <assembly fasta> -m 2500 -o <sample_name>_assembly_2500.fa

bbmap.sh -Xmx48G threads=<number of threads> minid=85 overwrite=t ref=<sample_name>_assembly_2500.fa in1=<sample_name>_R1_trimmed.fastq in2=<sample_name>_R2_trimmed.fastq out= <sample_name>_mapped.sam

samtools view -@ <number of threads> -bS <sample_name>_mapped.sam > <sample_name>_mapped.bam

samtools sort -T <sample_name> 2500.sorted -o =<sample_name>_2500.sorted.bam <sample_name>_mapped.bam -@ <number of threads>

runMetaBat.sh <sample_name>_assembly_2500.fa <sample_name> 2500.sorted.bam

Pathway 2: 

Read trimming with sickle:

sickle pe -f <forward untrimmed reads>  -r <reverse untrimmed reads>   -t <sequencing platform> -o <sample_name>_R1_trimmed.fastq -p <sample_name>_R2_trimmed.fastq -s discared_R1R2.fastq

Assembly with idba-ud:
reformat.sh in1=<forward trimmed reads> in2=<reverse trimmed reads> out1=<forward trimmed reads>_25pcnt.fastq out2=<reverse trimmed reads>_25pcnt.fastq samplerate=0.25 -sampleseed=1234


fq2fa --merge --filter <forward trimmed reads>_25pcnt.fastq <reverse trimmed reads>_25pcnt.fastq R1R2_ALL_trimmed_25pcnt.fa

idba_ud -r R1R2_ALL_trimmed_25pcnt.fa -o idba_assembled_output_25pcnt --num_threads 20

Binning with metabat2:

pullseq.py -i <assembly fasta> -m 2500 -o <sample_name>_assembly_2500.fa

bbmap.sh -Xmx48G threads=<number of threads> minid=85 overwrite=t ref=<sample_name>_assembly_2500.fa in1=<sample_name>_R1_trimmed.fastq in2=<sample_name>_R2_trimmed.fastq out= <sample_name>_mapped.sam

samtools view -@ <number of threads> -bS <sample_name>_mapped.sam > <sample_name>_mapped.bam

samtools sort -T <sample_name> 2500.sorted -o =<sample_name>_2500.sorted.bam <sample_name>_mapped.bam -@ <number of threads>

runMetaBat.sh <sample_name>_assembly_2500.fa <sample_name> 2500.sorted.bam

Pathway 3 (Note this is the JGI pipeline, if your sequencing is from JGI, you do not need to run this pathway): 

Read trimming and filtering with rqc2filter:

rqcfilter2.sh jni=t in=<interleaved untrimmed reads>.fastq.gz path=filter  rqcfilterdata=/home/opt/RQCFilterData rna=f trimfragadapter=t qtrim=r trimq=0 maxns=3 maq=3 minlen=51 mlf=0.33 phix=t removehuman=t removedog=t removecat=t removemouse=t khist=t removemicrobes=t sketch kapa=t clumpify=t barcodefilter=f trimpolyg=5 usejni=f

Assembly with metaspades:
bbcms.sh mincount=2 highcountfraction=0.6 in=<interleaved trimmed and filtered reads>.fastq.gz out=<interleaved trimmed and filtered reads>_bbcms.fastq.gz

spades.py -o <sample name>_metaspades_assembly --only-assembler -k 33,55,77,99,127 --meta -t <number of threads> --12 <interleaved trimmed and filtered reads>_bbcms.fastq.gz

Binning with metabat2:

pullseq.py -i <assembly fasta> -m 2500 -o <sample_name>_assembly_2500.fa

bbmap.sh -Xmx48G threads=<number of threads> minid=85 overwrite=t ref=<sample_name>_assembly_2500.fa in1=<sample_name>_R1_trimmed.fastq in2=<sample_name>_R2_trimmed.fastq out= <sample_name>_mapped.sam

samtools view -@ <number of threads> -bS <sample_name>_mapped.sam > <sample_name>_mapped.bam

samtools sort -T <sample_name> 2500.sorted -o =<sample_name>_2500.sorted.bam <sample_name>_mapped.bam -@ <number of threads>

runMetaBat.sh <sample_name>_assembly_2500.fa <sample_name> 2500.sorted.bam

Pathway 4 (Note this should only be done once you have assembled, binned, and quality checked bins from pathways 1-3 above- see notes for quality checking bins below):

Concatenate all the bins from a single river

cat *fa > all_bins_<river_name>.fasta

Concatenate forward and reverse reads for a particular river and sample type (e.g. surfacewater versus porewater)

cat <forward SW reads 1>.fastq <forward SW reads 2>.fastq    > < all forward SW reads>.fastq
cat <reverse SW reads 1>.fastq < reverse SW reads 2>.fastq    > < all reverse SW reads>.fastq
cat <forward PW reads 1>.fastq <forward PW reads 2>.fastq    > < all forward PW reads>.fastq
cat <reverse PW reads 1>.fastq < reverse PW reads 2>.fastq    > < all reverse PW reads>.fastq

Run the following steps with PW and SW separately (below with PW as example)

Map concatenated reads to medium and high quality bins

bbmap.sh -Xmx48G threads=<number of threads> minid=99 overwrite=t ref=all_bins_<river_name>.fasta in1=< all forward SW reads>.fastq in2=< all reverse SW reads>.fastq outu1=MapToBins_R1_<number of bins>_bins.fastq outu2=MapToBins_R2_<number of bins>_bins.fastq ambiguous=all

Record number of reads that are unmapped

fq2fa --merge --filter MapToBins_R1_<number of bins>_bins.fastq MapToBins_R2_<number of bins>_bins.fastq R1R2_MaptoBins_<number of bins>.fa

Assembly with idba-ud

idba_ud -r R1R2_MaptoBins_<number of bins>.fa -o idba_output_ReadsThatDidNotMap_<number of bins> --num_threads <number of threads>

Binning with metabat2:

pullseq.py -i <assembly fasta> -m 2500 -o <sample_name>_assembly_2500.fa

bbmap.sh -Xmx48G threads=<number of threads> minid=85 overwrite=t ref=<sample_name>_assembly_2500.fa in1=MapToBins_R1_<number of bins>_bins.fastq in2=MapToBins_R2_<number of bins>_bins.fastq out= <sample_name>_mapped.sam

samtools view -@ <number of threads> -bS <sample_name>_mapped.sam > <sample_name>_mapped.bam

samtools sort -T <sample_name>_2500.sorted -o =<sample_name>_2500.sorted.bam <sample_name>_mapped.bam -@ <number of threads>

runMetaBat.sh <sample_name>_assembly_2500.fa <sample_name> 2500.sorted.bam

Quality check bins with checkM and then repeat pathway 4 with new bin set that includes new medium and high quality bins recovered from coassembly and pathways 1-3.

Bin analyses

Quality checking (checkM)

checkm lineage_wf -t <number of threads>  -x fa . ./checkM 
checkm qa checkM/lineage.ms checkM -o 1 -f analyze_bins.txt --tab_table

Taxonomy assignment

gtdbtk classify_wf --extension fa --genome_dir <genome dir> --out_dir <genome dir>/gtdb --min_perc_aa 0 –cpus <number of threads>

Annotation

DRAM.py annotate -i <bins> -o <output dir name>  --min_contig_size 2500 --threads <number of threads>
