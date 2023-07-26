
The analyses shown here were used to create Figure 5.

---
1. Pull DRAM annotations of genes "on" in metaT

```
pullseq_header_name.py -i genes.faa -o 41827_metaT_genes.faa -n on_genes.txt -e F
```

2. Annotate with DRAM v1.4.4
```
source /opt/Miniconda2/miniconda2/bin/activate DRAM1.4.4
DRAM.py annotate_genes -i 41827_metaT_genes.faa -o DRAM_1.4.4_02062023 --use_camper --use_fegenie --use_sulfur --threads 10
```

3. run R script **01_reformat_annotation.R** to reformat annotations file with fasta = bam
> note, remove "X" from first column, remove "NA" in CAMPER_bitScore column first

4. distill new annotations file to get distillate by site
```
DRAM.py distill -i 41827genes_58bams_annotations.tsv -o distilled_by_BAM
```

5. run R script **02_XYZ.R** to curate annotations into microbial lifestyles, and generate Figure 5.

