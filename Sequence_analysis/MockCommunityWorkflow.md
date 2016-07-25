#Part 1:  Mock community sequence processing
### PATH to Mock community data from Centralia 2014 16S amplicons
```
research/ShadeLab/Shade/20141230_16Stag_Centralia/20141230_B_16S_PE/Mock_Cmty_TCCTCTGTCGAC_L001_R*
```

### File names for forward (R1) and reverse (R2) read files
* Mock_Cmty_TCCTCTGTCGAC_L001_R1_001.fastq.gz
* Mock_Cmty_TCCTCTGTCGAC_L001_R2_001.fastq.gz

### Merge paired ends and quality filtering with usearch 8.1
```
usearch -fastq_mergepairs Mock_Cmty_TCCTCTGTCGAC_L001_R1_001.fastq -fastqout Mock.fastq -relabel @ -fastq_merge_maxee 1.0 -fastq_minmergelen 250 -fastq_maxmergelen 274 -fastq_nostagger

usearch v8.1.1803_i86linux32, 4.0Gb RAM (264Gb total), 20 cores
(C) Copyright 2013-15 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

#output files: Mock.fastq
```

### Deplication with UNOISE
```
usearch -derep_fulllength Mock.fastq -fastqout uniques_Mock.fastq -sizeout

#output files: uniques_Mock.fastq

#step 1.5 - remove singletons
usearch -sortbysize uniques_Mock.fastq -fastqout uniques_Mock_nosigs.fastq -minsize 2

#output files:   uniques_Mock_nosigs.fastq
```

### Denoise (pre-cluster) with parent sequences with USEARCH
```
usearch -cluster_fast uniques_Mock_nosigs.fastq -centroids_fastq denoised.fq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

#output files:  denoised.fq
```

### Pick OTUs with UPARSE, includes chimera detection
```
usearch -cluster_otus denoised.fq -otus mock_denoised_otus.fa -relabel OTU_ -sizeout -uparseout results.txt   

#output files:  mock_denoised_otus.fa, results.txt   
```

### Additional reference-based chimera detection (only for mock)
```
usearch -uchime_ref mock_denoised_otus.fa -db gold.fa -strand plus -nonchimeras mock_denoised_NoChimeraRef_otus.fa

output files: mock_denoised_NoChimeraRef_otus.fa
```

### Mapping de-noised reads to defined OTUs.  Anything that does not hit at 97% identity or greater (e.g., chimeras) to the OTUs will be discarded.   
```
usearch -usearch_global denoised.fq -db mock_denoised_NoChimeraRef_otus.fa -strand plus -id 0.97 -uc map_denoised.uc -otutabout Mock_OTU_table.txt

#output files:  map.uc, Mock_OTU_table.txt
```

### Assigning taxonomy with the RDP Classifier
```
module load RDPClassifier/2.9

java -jar $RDP_JAR_PATH/classifier.jar classify -c 0.5 -o mock_denoised_classified.txt -h mock_hier.txt mock_denoised_NoChimeraRef_otus.fa

#output files: mock_hier.txt, mock_denoised_classified.txt
```

### Make a database of contaminant OTUs detected in mock community by removing real rep. OTU sequences (mock type stains) from `mock_denoised_NoChimeraRef_otus.fa`.  New db file is: `mock_craptaminant_OTU_db.fa`


#Part 2.  Community sequence processing using USEARCH tools
### PATH to raw fastq (zipped) on ShadeLab HPCC research space

```
research/ShadeLab/Shade/20141230_16Stag_Centralia/20141230_A_16S_PE/
research/ShadeLab/Shade/20141230_16Stag_Centralia/20141230_B_16S_PE/
```

### i.  Merge fastq reads - executed as a qsub 'merge.qsub'

```
mkdir mergedfastq

for file in $(<merge_fq_list.txt)
do

    usearch -fastq_mergepairs ${file} -fastqout mergedfastq/${file}_merged.fastq -relabel @ -fastq_merge_maxee 1.0 -fastq_minmergelen 250 -fastq_maxmergelen 274 -fastq_nostagger

done
```

### ii.  Pool all merged-paired ends across samples

```
find ./mergedfastq -type f -name '*fastq' -exec cat '{}' > combined_merged.fastq ';'

## output file:  combined_merged.fastq
```

### iii.  De-replicate

```
usearch -derep_fulllength combined_merged.fastq -fastqout uniques_combined_merged.fastq -sizeout

## output file: uniques_combined_merged.fastq
```

### iv.  Remove single sequences

```
usearch -sortbysize uniques_combined_merged.fastq -fastqout nosigs_uniques_combined_merged.fastq -minsize 2
## output file: nosigs_uniques_combined_merged.fastq
```

### v. Precluster (denoise) - for our dataset on the 64-bit usearch

```
usearch -cluster_fast nosigs_uniques_combined_merged.fastq -centroids_fastq denoised_nosigs_uniques_combined_merged.fastq -id 0.9 -maxdiffs 5 -abskew 10 -sizein -sizeout -sort size

## output file: denoised_nosigs_uniques_combined_merged.fastq
```

### vi.  Remove sequences that match 100% to our craptaminant database (using closed-reference 100% OTU picking)

```
usearch -search_exact denoised_nosigs_uniques_combined_merged.fastq -db mock_craptaminant_OTU_db.fa -otus craptaminantOTUs_denoised_nosigs_uniques_combined_merged.fa -notmatchedfq nocrap_denoised_nosigs_uniques_combined_merged.fastq -matchedfq craptaminantSeqs_denoised_nosigs_uniques_combined_merged.fastq -strand plus


## output files:
### craptaminantOTUs_denoised_nosigs_uniques_combined_merged.fa;
### nocrap_denoised_nosigs_uniques_combined_merged.fastq
```

### vi.  Reference-based OTU picking using usearch_global: Cluster sequences at 97% identity to the greengenes database, version 13.8

```
usearch -usearch_global nocrap_denoised_nosigs_uniques_combined_merged.fastq -id 0.97 -db gg_13_8_otus/rep_set/97_otus.fasta -notmatchedfq RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fastq -strand plus -uc RefMatchOTUMap_nocrap_denoised_nosigs_uniques_combined_merged.uc -dbmatched gg_97_rep_set_matched.fa

## output files:
### RefMatchOTUMap_nocrap_denoised_nosigs_uniques_combined_merged.uc (usearch standard results table);
### RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fastq (sequences that did not hit the gg db and need to be clustered de novo)
### gg_97_rep_set_matched.fa (gg 97_rep_set matches to our dataset - add to MASTER OTU rep. sequences)
```

### vii.  De novo OTU picking using usearch cluster_otus:  cluster sequences at 97% identity (includes chimera checking with uparse) This step removes singletons.

```
usearch -cluster_otus RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fastq -minsize 2 -otus DeNovoUclustOTUs_RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fa -relabel OTU_dn_ -uparseout DeNovoUclustResults_RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.up

## output files:
### DeNovoUclustOTUs_RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fa (representative sequences for de novo OTUs)
### DeNovoUclustResults_RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.up (uparse standard results table)
```

### viii.  Combine ref-based and de novo representative OTU sequences into one master OTU "db" file.

```
cat gg_97_rep_set_matched.fa DeNovoUclustOTUs_RefNoMatch_nocrap_denoised_nosigs_uniques_combined_merged.fa > RepSeqs.fa

## output files
### RepSeqs.fa
```

### ix.  Map all sequences (pre-dereplication) back to OTU definitions using usearch_global.  Any sequences that do not hit the new OTU database are discarded.

```
usearch -usearch_global combined_merged.fastq -db RepSeqs.fa  -strand plus -id 0.97 -uc OTU_map.uc -otutabout OTU_table.txt -biomout OTU_jsn.biom

## output files
### OTU_jsn.biom (Biom OTU table)
### OTU_table.txt (classic OTU table)
```


#Part 3:  Continued sequence processing in QIIME 1.9.1
### i.  Convert jsn biom format to HD5 and add QIIME metadata
```
biom convert -i OTU_jsn.biom -o OTU_hdf5.biom --table-type="OTU table" --to-hdf5

## output files
### OTU_hdf5.biom
```

### ii.  Alignment against silva v123 and taxonomic assignment with RDP Classifier 2.2.  We chose silva because the alignment was higher quality; however, greengenes and silva alignments both resulted in similar proportions of failures
```
align_seqs.py -i RepSeqs.fa -t core_alignment_SILVA123.fasta -o qiime191_pynast_silva123/

## output
### /qiime191_pynast_silva123/RepSeqs_failures.fasta
### /qiime191_pynast_silva123/RepSeqs_aligned.fasta
### /qiime191_pynast_silva123/RepSeqs_log.txt
```

### iii. Filter failed alignments from OTU table
```
#from OTU table
filter_otus_from_otu_table.py -i OTU_hdf5.biom -o OTU_hdf5_filteredfailedalignments.biom -e qiime191_pynast_silva123/RepSeqs_failures.fasta

#from RepSeqs file
filter_fasta.py -f RepSeqs.fa -o RepSeqs_filteredfailedalignments.fa -a qiime191_pynast_silva123/RepSeqs_aligned.fasta

## output files
### OTU_hdf5_filteredfailedalignments.biom
### RepSeqs_filteredfailedalignments.fa
```

### iv. Assign taxonomy using RDP Classifier 2.2 against gg 97 database v. 13.8

```
assign_taxonomy.py -i RepSeqs_filteredfailedalignments.fa -m rdp -c 0.8 -t gg_13_8_otus/taxonomy/97_otu_taxonomy.txt -r gg_13_8_otus/rep_set/97_otus.fasta -o rdp_assigned_taxonomy22/

#Add taxonomy as metadata to .biom table

echo "#OTUID"$'\t'"taxonomy"$'\t'"confidence" > templine.txt

cat  templine.txt rdp_assigned_taxonomy22/RepSeqs_filteredfailedalignments_tax_assignments.txt >> rdp_assigned_taxonomy22/RepSeqs_filteredfailedalignments_tax_assignments_header.txt

#source /mnt/research/ShadeLab/software/loadanaconda2.sh

biom add-metadata -i OTU_hdf5_filteredfailedalignments.biom -o OTU_hdf5_filteredfailedalignments_rdp.biom --observation-metadata-fp rdp_assigned_taxonomy22/RepSeqs_filteredfailedalignments_tax_assignments_header.txt --sc-separated taxonomy --observation-header OTUID,taxonomy

rm templine.txt

## output files
### rdp_assigned_taxonomy/MASTER_RepSeqs_filteredfailedalignments_tax_assignments.log
### rdp_assigned_taxonomy/MASTER_RepSeqs_filteredfailedalignments_tax_assignments.txt
### rdp_assigned_taxonomy/MASTER_RepSeqs_filteredfailedalignments_tax_assignment_header.txt
### OTU_hdf5_filteredfailedalignments_rdp.biom
```

### vi.  Filter chloroplast and mitochondria from OTU table and RepSeqs file
```
#remove chloroplasts and mitochondria (keep cyanobacteria that are not chloroplasts)
filter_taxa_from_otu_table.py -i OTU_hdf5_filteredfailedalignments_rdp.biom -o OTU_hdf5_filteredfailedalignments_rdp_rmCM.biom -n o__Streptophyta,o__Chlorophyta,f__mitochondria

#remove same Mito and Chloro sequences from RepSeqs file
filter_fasta.py -f RepSeqs_filteredfailedalignments.fa -o MASTER_RepSeqs_filteredfailedalignments_rmCM.fa -b OTU_hdf5_filteredfailedalignments_rdp_rmCM.biom

#optional sanity check: summarize new biom table and check that it has fewer sequences than original
#original
biom summarize-table -i OTU_hdf5_filteredfailedalignments_rdp.biom
#filtered
biom summarize-table -i OTU_hdf5_filteredfailedalignments_rdp_rmCM.biom

#optional sanity check:  count seqs in new fasta, and check that it has fewer than original
#orginal
count_seqs.py -i RepSeqs_filteredfailedalignments.fa
#filtered
count_seqs.py -i MASTER_RepSeqs_filteredfailedalignments_rmCM.fa


## output files
### MASTER_RepSeqs_filteredfailedalignments_rmCM.fa
### OTU_hdf5_filteredfailedalignments_rdp_rmCM.biom
```

### vii.  Make phylogeny using FastTree (version)
```
#First, clean alignment by omitting highly variable regions before tree building - will make tree building more efficient
filter_alignment.py -i qiime191_pynast_silva123/RepSeqs_aligned.fasta -o clean_alignment/

#make phylogeny and root tree
make_phylogeny.py -i clean_alignment/RepSeqs_aligned_pfiltered.fasta -o MASTER_RepSeqs_aligned_clean.tre -r tree_method_default

## output files
### clean_alignment/RepSeqs_aligned_pfiltered.fasta
### MASTER_RepSeqs_aligned_clean.tre
```

### viii. Subsampling, collapsing the OTU table, exploratory rarefaction
```
#subsampling for analysis of technical variability in DNA extraction replicates
single_rarefaction.py -i OTU_hdf5_filteredfailedalignments_rdp_rmCM.biom -o OTU_hdf5_filteredfailedalignments_rdp_rmCM_even53000.biom -d 53000

#collapse table across technical (DNA extraction) replicates
collapse_samples.py -b OTU_hdf5_filteredfailedalignments_rdp_rmCM.biom -m Centralia_Full_Map.txt --output_biom_fp OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse.biom --output_mapping_fp Centralia_Collapsed_Map.txt --collapse_mode sum --collapse_fields Sample

#optional summary
biom summarize-table -i OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse.biom

#subsample to include mock community
single_rarefaction.py -i OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse.biom -o OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even227100.biom -d 227100

#subsample to omit mock community
single_rarefaction.py -i OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse.biom -o MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.biom -d 321000

#Exploratory rarefaction to determine sampling efforts
alpha_rarefaction.py -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.biom -m Centralia_Collapsed_Map.txt -a -n 15 -o alpha_rarefaction_collapsed/ -t MASTER_RepSeqs_aligned_clean.tre

### output files
## OTU_hdf5_filteredfailedalignments_rdp_rmCM_even53000.biom
## OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse.biom
## OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even227100.biom
## MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.biom
## alpha_rarefaction_collapsed/
```

### xiv.  Diversity calculations
```
#alpha diversity analyses on even, collapsed dataset
alpha_diversity.py -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.biom -m PD_whole_tree,observed_otus -t MASTER_RepSeqs_aligned_clean.tre -o MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_alphadiv.txt

#beta diversity analysis on even, collapsed dataset
beta_diversity.py -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.biom -t MASTER_RepSeqs_aligned_clean.tre -o betadiv_even321000/ -m unweighted_unifrac,weighted_unifrac

#summarize at phylum level
summarize_taxa.py -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.biom -L 2 -o phylum_summary/
```

### xv.  Convert biom for export
```
biom convert -i MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.biom -o MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage"
```

### Move files to local computer for analysis in R (remember to format to remove hash (#) from first row header
* OTU table: MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt
* phylum summary: phylum_summary/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_L2.txt
* alpha diversity: MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000_alphadiv.txt
* betadiv_even321000/weighted_unifrac_MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt
