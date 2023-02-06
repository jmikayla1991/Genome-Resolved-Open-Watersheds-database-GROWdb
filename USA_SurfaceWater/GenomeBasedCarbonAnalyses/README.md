## Add CAMPER annotations to existing MAG set
Used MAG annotations located here 
```
/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/all_bins/dRep_v2.6.2/dereplicated_genomes/merged_annotations_drep_bins
```

```
conda activate CAMPER
camper_annotate -i genes.faa -a annotations.tsv -o add_camper
camper_distill  -a add_camper/annotations.tsv -o add_camper/camper_distillate.tsv
```
---
## functional analysis of metaT genes

Working here:
```
/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/USA/carbon_analyses
```

pull gene seqs of "on" genes
```
pullseq_header_name.py -i /home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/all_bins/dRep_v2.6.2/dereplicated_genomes/merged_annotations_drep_bins/genes.faa -o 41827_metaT_genes.faa -n on_genes.txt -e F
```

annotate genes using camper, fe, s, dram 1.4.4
```
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
DRAM.py annotate_genes -i 41827_metaT_genes.faa -o DRAM_1.4.4_02062023 --use_camper --use_fegenie --use_sulfur --threads 10
```
log here: 
```
/home/projects-wrighton-2/GROWdb/USAfocus_FinalBins110121/dereplicated_bin_analyses/metaT_mapping/bowtie_mapping_final/USA/carbon_analyses/slurm-57468.out
```

ran R script XX to reformat annotations file with fasta = bam
> note, remove "X" from first column, remove "NA" in CAMPER_bitScore column first

distilled new annotations file to get distillate by site
```
DRAM.py distill -i 41827genes_58bams_annotations.tsv -o distilled_by_BAM
```
