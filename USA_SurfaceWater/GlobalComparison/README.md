## Comparison of US Surface Water GROW MAGs to global MAGs
### Sources of global MAGs
1. GROW global sites to date - This includes 1,286 additional MAGs derived from 23 metagenomes released in this study, including 11 sites beyond the United States (UK, Canada, Italy, Germany, Israel, Republic of the Congo).

>> [Link to list of 1,286 MAGs used here.](https://github.com/jmikayla1991/Genome-Resolved-Open-Watersheds-database-GROWdb/blob/main/USA_SurfaceWater/GlobalComparison/GROW_global_1286.csv) - GROW_global_1286.csv
   
2. [GROW Erpe river 48-hr study, Rodriguez-Ramos, et al.](https://www.frontiersin.org/articles/10.3389/frmbi.2023.1199766/full) - As a part of the GROW effort, the Erpe River in Germany was sequenced every 3 hours over a 48 hour period in 2018 resulting in 32 metagenomes (surface water and sediment). In total, 1,033 MAGs were reconstructed and dereplicated into 125 unique MAGs. Here we retain 75 MAGs from surface waters for our global comparison.

>> [Link to list of 75 MAGs used here.]() - Rodriguez-Ramos_GROWcompare_75.csv

3. [Nayfach, et al.](https://www.nature.com/articles/s41587-020-0718-6) - In an effort to catalog Earth's microbiomes, this study performed metagenomics on samples collected from diverse habitats covering all of Earth's continent and oceans. This comprehensive catalog includes 52,215 MAGs, of which we retain 5,336 MAGs from **aquatic freshwater** environments for our global comparison. Specifically these include **Aquatic, Environmental, Freshwater lake, Freshwater, Freshwater lake sediment, watersheds, Freshwater lentic, freshwater microbial mat, freshwater lake hypolimnion, Anoxic lake water, Surface water, Lake, Lentic,** and **Wetland** habitats. We have excluded soil, sediment, and wastewater related habitats within the aquatic freshwater ecosystem to more directly compare to surface water GROW.

>> [Link to list of 5,336 MAGs used here.](https://github.com/jmikayla1991/Genome-Resolved-Open-Watersheds-database-GROWdb/blob/main/USA_SurfaceWater/GlobalComparison/Nayfach_GROWcompare_5336.csv) - Nayfach_GROWcompare_5336.csv
   
4. [Garner, et al.](https://www.nature.com/articles/s41564-023-01435-6) - This large scale freshwater metagenome study of 308 Canadian lakes resulted in a catalog of 1,008 MAGs. In an interest to compare to other inland waterbodies, we have retained 1,008 MAGs for our global comparison.

>> [Link to list of 1,008 MAGs used here.](https://github.com/jmikayla1991/Genome-Resolved-Open-Watersheds-database-GROWdb/blob/main/USA_SurfaceWater/GlobalComparison/Garner_GROWcompare_1008.csv) - Garner_GROWcompare_1008.csv

### Compilation of MAGs from other freshwater global studies resulted in 7,705 MAGs to compare to 2,093 United States surface water GROW MAGs. 

#### How similar are the set of 9,798 MAGs recovered from various studies at 99% ID? 
We will use drep_3.0.0 to which genomes are similar at 99% ID. 

**99% ID**
```
dRep dereplicate dRep_v3.0.0_9798MAGs -p 20 -comp 50 -con 10 -g ./*fa 
```


