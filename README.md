# Optimizing Methods for Insect Metabarcoding

![](./docs/Figure_41.png)

## Introduction
This project goal was identifying the number of subsamples required to detect total estimated insect OTU diversity in a light trap.

## Methods
* Prepare databases.
  
  * [Download](https://mhenso.github.io/insect_metabarcoding/1_bold/all_insecta.nb.html) from BOLD repository
  
  * [Cleaning](https://github.com/mhenso/insect_metabarcoding/1_bold/db5.nb.html) database entries
    

* [Step 2](https://mhenso.github.io/public/docs/db5.nb.html) - Reshape and compare the databases.
  
  * [Shell_scripting](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) Using seqkit, mafft, usearch, blast+ and R
    
  * [Remove NUMTS and perform blastn/tblastx](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) db5_nonumts.nb.html
 
  * [Blast reseults](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) blastn300.nb.html
 
  * [Figure 2](https://github.com/mhenso/insect_metabarcoding/b) bitscore.nb.html  [Figure 3](https://github.com/mhenso/insect_metabarcoding/b) blastn300_venn.nb.html
 
  * [Saturation](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh)
  

* [OTUs pipeline](https://mhenso.github.io/public/docs/db5.nb.html)

  * [Shell_scripting](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) Using usearch and R
 
  * [Remove NUMTS](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh)
 
  * [Rarefaction](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) otus300.nb.html


* [OTUs pipeline](https://mhenso.github.io/public/docs/db5.nb.html) - OTUs pipeline for external dataset.

  * [Shell_scripting](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) Using usearch and R # notes_zizka_cutadapt_subsample.sh
 
  * [Remove NUMTS](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) # otus_numts.nb.html
 
  * [Rarefaction](https://github.com/mhenso/insect_metabarcoding/blob/main/2_databases/2_databases_revised.sh) # otus300.nb.html


## Reference
Sihaloho HF, Azhar I, Gani M, Senawi J, Yong LS, Ayub Q, De R, Kingston T, Phillips CD. in press. Optimizing Methods for Insect Metabarcoding. 

## SRA
The raw reads generated in this study are available at GenBank Sequence Read Archive, BioProject accession number [PRJNA1260500](http://www.ncbi.nlm.nih.gov/bioproject/1260500)

```bash
./Sihaloho_etal
├──1_bold
│   ├──rds
│   └──results
├──2_databases
│   ├──data
│   ├──results
│   │   ├──global
│   │   │   ├──blast300
│   │   │   ├──confirm_coordinate
│   │   │   └──usearch
│   │   ├──local
│   │   │   ├──blast300
│   │   │   └──usearch
│   │   └──regional
│   │       ├──blast300
│   │       └──usearch
│   ├──results_blastn300
│   ├──results_figures
│   └──results_tblastx300
├──3_cutadapt
│   ├──filtered
│   ├──lt_db
│   ├──output
│   │   └──otus_finalize
│   └──stitched
├──4_saturation
│   └──data
│       ├──dambe
│       └──global
│           └──dambe
└──COI_252_samples
```




