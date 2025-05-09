# Optimizing Methods for Insect Metabarcoding

![](./docs/Figure_41.png)

## Introduction
This project goal was identifying the number of subsamples required to detect total estimated insect OTU diversity in a light trap.

## Codes
* Prepare databases.
  
  * [Download](https://mhenso.github.io/insect_metabarcoding/1_bold/all_insecta.nb.html) from BOLD repository
  
  * [Subsetting](https://mhenso.github.io/insect_metabarcoding/1_bold/db5.nb.html) database entries
    

* Reshape and compare the databases.
  
  * [Shell_scripting](https://github.com/mhenso/insect_metabarcoding/tree/main/2_databases/2_databases_revised.sh) Using seqkit, mafft, usearch, blast+ and R
    
  * [Remove NUMTS and perform blastn/tblastx](https://mhenso.github.io/insect_metabarcoding/2_databases/db5_nonumts.nb.html) 
 
  * [Blast results](https://mhenso.github.io/insect_metabarcoding/2_databases/blastn300.nb.html) 
 
  * [Figure 2](https://mhenso.github.io/insect_metabarcoding/2_databases/bitscore.nb.html) and [Figure 3](https://mhenso.github.io/insect_metabarcoding/2_databases/blastn300_venn.nb.html)
 
  * [Saturation]()
  

* OTUs pipeline

  * [Shell_scripting](https://github.com/mhenso/insect_metabarcoding/blob/main/3_cutadapt/3_cutadapt_revised.sh) Using usearch and R
 
  * [Remove NUMTS](https://mhenso.github.io/insect_metabarcoding/3_cutadapt/otus_numts.nb.html)
 
  * [Sample Coverage](https://mhenso.github.io/insect_metabarcoding/3_cutadapt/otus300.nb.html)


* OTUs pipeline for external dataset.

  * [Shell_scripting](https://github.com/mhenso/insect_metabarcoding/blob/main/4_cutadapt_external_dataset/notes_zizka_cutadapt_subsample.sh) Using usearch and R
 
  * [Remove NUMTS](https://mhenso.github.io/insect_metabarcoding/4_cutadapt_external_dataset/otus_numts.nb.html)
 
  * [Sample Coverage](https://mhenso.github.io/insect_metabarcoding/4_cutadapt_external_dataset/otus300.nb.html)


## Reference
Sihaloho HF, Azhar I, Gani M, Senawi J, Yong LS, Ayub Q, De R, Kingston T, Phillips CD. in press. Optimizing Methods for Insect Metabarcoding. 

## SRA
The raw reads generated in this study are available at GenBank Sequence Read Archive, BioProject accession number [PRJNA1260500](http://www.ncbi.nlm.nih.gov/bioproject/1260500)

## Folder structure
```bash
 ./Sihaloho_etal
 ├── 1_bold
 │   ├── rds
 │   ├── results
 ├── 2_databases
 │   ├── blast_db
 │   ├── blast_db_results
 │   ├── figures
 │   └── results
 │   │   ├── global
 │   │   ├── local
 │   │   └── regional
 ├── 3_cutadapt
 │   ├── filtered
 │   ├── lt_db
 │   ├── output
 │   └── stitched
 ├── 3_cutadapt_zizka_subsample
 │   ├── cut
 │   ├── fastq
 │   ├── fastq_subsample
 │   ├── filtered
 │   ├── output
 │   ├── plotting
 │   ├── rds
 │   └── stitched
 ├── COI_252_samples
 └── docs
```


