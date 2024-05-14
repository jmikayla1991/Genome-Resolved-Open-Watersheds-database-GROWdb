## This respository contains the code and commands for "A functional microbiome catalog crowdsourced from North American rivers" currently on [BioArchive](https://www.biorxiv.org/content/10.1101/2023.07.22.550117v1). 

### Information in this repository
#### Directories
1. Figures- R scripts for data analysis and figure generation
2. GenomeBasedCarbonAnalyses- R scripts and commands for assigning microbial genomes to carbon usage patterns, as shown in Fig. 5 of the manuscript
3. MetaG_mapping- commands used for mapping metagenomes to dereplicated MAG database
4. MetaT_mapping- commands used for mapping metatranscriptomes to dereplicated MAG database and R scripts for data normalization
5. Metagenomic_Pipeline_GROWdb- commands used for read trimming, assembly, and binning of GROWdb metagenome-assembled genomes
6. Geospatial_Data- R code and associated geospatial data generated for each sampling location in this manuscript

Other relevant GitHub respositories: 
- https://github.com/rossyndicate/GROWdb for geospatial data mining for GROWdb

Command line tools and associated versions used here are included in requirements.txt

### Data accessibility 
The data underlying GROWdb are accessible across multiple platforms to ensure many levels of data use and structure are widely available. First, all reads and MAGs are publicly hosted on National Center for Biotechnology (NCBI) under Bioproject PRJNA946291. Second, all data presented in this manuscript including MAG annotations, phylogenetic tree files, antibiotic resistance gene database files, and expression data tables are available in Zenodo (https://doi.org/10.5281/zenodo.8173287). 

Beyond the content listed above, our aim for GROWdb was to maximize data use by making the data available in searchable and interactive platforms including the National Microbiome Data Collaborative (NMDC) data portal, the Department of Energy’s Systems Biology Knowledgebase (KBase), and a GROW specific user interface released here, GROWdb Explorer. Each platform provides different ways to interact with data in the GROWdb: 
- NMDC GROWdb was a flagship project for the newly formed NMDC. Specifically, individual GROWdb datasets (metagenomes, metatranscriptomes, etc) are easily accessible and searchable through the NMDC data portal (https://data.microbiomedata.org/), where they are systematically connected to each other and to a rich suite of sample information, other data collected on the same samples, and standard analysis results, following Findable, Accessible, Interoperable, and Reusable (FAIR) data practices.
- KBase GROWdb is a publicly available collection within KBase3, with samples, MAGs, and corresponding genome scale metabolic models found in the KBase narrative structure (https://doi.org/10.25982/109073.30/1895615). Access within KBase allows for immediate access and reuse of data, including comparison to private data analyses using KBase’s 500+ analysis tools, in a point and click format.
- GROWdb Explorer is a graphical user interface built through the Colorado State University Geospatial Centroid (https://geocentroid.shinyapps.io/GROWdatabase/), allowing users to search and graph microbial and spatial data simultaneously. Here the microbial data was distilled into functional gene information, so that biogeochemical contributions and the microorganisms catalyzing them can be assessed and visualized rapidly across the dataset. 

In summary, GROWdb represents the first publicly available genome collection from rivers and offers data that can be leveraged across microbiome studies. GROWdb is an expanding repository to incorporate and unify global river multi-omic data for the future.


