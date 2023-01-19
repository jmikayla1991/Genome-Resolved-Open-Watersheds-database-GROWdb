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
