# sratoolkit
## Download sra info for creating sample list
esearch -db sra -query PRJNA883590 | esummary > temp_summary

# BASH
## create sample list
cat temp_summary | xtract -pattern DocumentSummary -element Biosample, LIBRARY_NAME, Submitter@acc, Experiment@acc, Sample@acc, Run@acc, Platform@instrument_model > PRJNA883590-info4.tsv
awk '{print $5}' PRJNA883590-info3.tsv > srr_list


# sratoolkit
## Download the SRA
prefetch --option-file srr_list

# Extract / split files
for s in $(cat srr_list); 
do 
fasterq-dump --threads 64 ./${s} \
-O ./fastq;
done


# SEQTK
## subsample
for x in $(cat srr_list); 
do
seqtk sample -s 8 ../zizka2/fastq/${x}_1.fastq 23119 > ./fastq_subsample/${x}_1.fastq
seqtk sample -s 8 ../zizka2/fastq/${x}_2.fastq 23119 > ./fastq_subsample/${x}_2.fastq;
done


# BASH
## rawreads
for x in $(cat srr_list);
do 
cat ./fastq_subsample/${x}_1.fastq | grep "^@SRR" | wc -l  >> R1_reads.txt; 
cat ./fastq_subsample/${x}_2.fastq | grep "^@SRR" | wc -l  >> R2_reads.txt; 
done
## comment: exact number of reads were sampled

mkdir cut
nano cutadapt.sh

# cutadapt.sh ---
OMP_NUM_THREADS=128
conda activate env_cutadapt

# forward
for x in $(cat srr_list);
do
# forward GGDACWGGWTGAACWGTWTAYCCHCC
cutadapt -j 128 -g GGDACWGGWTGAACWGTWTAYCCHCC ./fastq_subsample/${x}_1.fastq.gz -o ./cut/${x}_1.fastq
# reverse TANACYTCNGGRTGNCCRAARAAYCA
cutadapt -j 128 -g TANACYTCNGGRTGNCCRAARAAYCA ./fastq_subsample/${x}_2.fastq.gz -o ./cut/${x}_2.fastq;
done
# cutadapt.sh ===

### Notes: you need to run "conda activate env_cutadapt" before nohup it
conda activate env_cutadapt
nohup bash cutadapt_lt.sh &> cutadapt_lt.sh.out &


# PEAR
## stitch consensus read pairs for each sample
mkdir stitched
nano pear.sh
# pear.sh ---
for x in $(cat srr_list);
do
pear --threads 128 --memory 100G \
-f ./cut/${x}_1.fastq \
-r ./cut/${x}_2.fastq \
-o ./stitched/${x}.fastq;
done
# pear.sh ===
nohup bash pear.sh &> pear.sh.out &


# USEARCH
## filter consensus reads to only include consensus reads with estimated 1 error max
mkdir filtered && nano usearch_fastq_filter.sh

# usearch_fastq_filter.sh ---
export OMP_NUM_THREADS=128

for s in $(cat srr_list)
do
$usearch -fastq_filter stitched/${s}.fastq.assembled.fastq \
-fastaout filtered/${s}_filtered.fa \
-fastq_maxee 1.0 -threads 128;
done
# usearch_fastq_filter.sh ===
nohup bash usearch_fastq_filter.sh &> usearch_fastq_filter.sh.out &


# BASH
## relabel the reads so we can keep track after pooling of which sample each came from
nano relabel.sh
# relabel.sh ---
for s in $(cat srr_list)
do
sed "-es/^>\(.*\)/>\1;barcodelabel=${s};/" < filtered/${s}_filtered.fa > filtered/${s}.fa;
done
# relabel.sh ====
nohup bash relabel.sh &> relabel.sh.out &

## pool the relabeled concensus sequences from each sample
mkdir output 
cat ./filtered/*{0..9}.fa > ./output/combined.fa

cat output/combined.fa | grep -c "^>"
# 4017978

# USEARCH
## make a fasta with only unique sequences but keeping track of abundance

nano usearch_fastx_uniques.sh

# usearch_fastx_uniques.sh ---
export OMP_NUM_THREADS=128

$usearch -fastx_uniques ./output/combined.fa \
-fastaout ./output/uniques.fa \
-sizeout -relabel Uniq \
-tabbedout ./output/report_fastx_unique.txt \
-threads 128
# usearch_fastx_uniques.sh ===
nohup bash usearch_fastx_uniques.sh &> usearch_fastx_uniques.sh.out &


cat ./output/uniques.fa | grep -c "^>"
# 1,427,773

## identify otus
nano usearch_cluster_otus.sh
# usearch_cluster_otus.sh ---
export OMP_NUM_THREADS=128

$usearch -cluster_otus ./output/uniques.fa \
-otus ./output/otus97.fa \
-relabel Otu -threads 128
# usearch_cluster_otus.sh ===
nohup bash usearch_cluster_otus.sh &> usearch_cluster_otus.sh.out &

cat ./output/otus97.fa | grep -c "^>"
# 4,820


## build the otu table from the identified zotu and their counts across samples
nano usearch_otutab.sh
# usearch_otutab.sh ---
export OMP_NUM_THREADS=128

$usearch -threads 128 \
-otutab ./output/combined.fa \
-otus ./output/otus97.fa \
-otutabout ./output/otutab.txt \
-mapout ./output/otumap.txt 
# usearch_otutab.sh ===
nohup bash usearch_otutab.sh &> usearch_otutab.sh.out &


# OTUS finalize ------------------------------------------------------------ DONE
pwd
# /lustre/work/hsihaloh/zizka_cutadapt_subsample/output/otus_finalise

seqkit seq -g -w 0 ../otus97.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > otus97.fa.tab

## subsample 200 id
grep "^>" ../otus97.fa | head -n 200 | sed 's/^>//' > otus200.list
seqkit grep -n -f ./otus200.list ../otus97.fa > ./otus97.200.fa
## align the subsample otus
mafft --auto --thread 8 ./otus97.200.fa > ./otus97.200.afa


## align rest otus
## filter OUT otus that have been aligned (200otus)
seqkit grep -v -n -f ./otus200.list ../otus97.fa > ./otus97rest.fa
grep -c "^>" otus97rest.fa
# 4,620  ### Notes: the number is correct

mafft --auto --addfull ./otus97rest.fa --keeplength --thread 32 otus97.200.afa > otus97.200.rest.afa

grep -c "^>" otus97.200.rest.afa
# 4,820

## see the frequencies of length
seqkit seq -g -w 0 otus97.200.rest.afa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > otus97.200.rest.afa.tab

## obtain min 300 length --------------------
## create a list have seq length min 300
seqkit seq -g --min-len 300 -w 0 ./otus97.200.rest.afa  | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 4,505 list.min300

## obtain min300 length
seqkit grep -n -f list.min300 otus97.200.rest.afa | seqkit seq --upper-case -w 0 > otus97.200.rest.min300.afa

grep -c "^>" otus97.200.rest.min300.afa
# 4,505 ==> ### Notes: 4,505/4,820 * 100% = 93.46% ===> lost 315 OTUs

## see the frequencies of length again
seqkit seq -g -w 0 otus97.200.rest.min300.afa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > otus97.200.rest.min300.afa.tab
cat otus97.200.rest.min300.afa.tab
   # 3950 313
    # 147 312
     # 50 311
     # 60 310
     # 11 309
     # 17 308
    # 183 307
     # 11 306
     # 16 305
     # 25 304
      # 7 303
      # 6 302
     # 11 301
     # 11 300

## Ungap/remove the gaps from the local database
seqkit seq -g ./otus97.200.rest.min300.afa > otus97.final.300.fa
grep -c "^>" otus97.final.300.fa
# 4,505
## obtain min 300 length ===================

# OTU TABLE
## subset otutable to mins 300bp otus refseq 

# create id from final refseq
cat otus97.final.300.fa | grep "^>" | sed 's/>//' > otus97.final.300.id
# add extra id for first line for sample name
sed -i '1 i\#OTU' otus97.final.300.id
## grep the otus from otutab.txt and save it as otutab.300.txt
grep -w -f otus97.final.300.id ../otutab.txt > ./otutab.300.txt

wc -l ./otutab.300.txt
# 4506 ./otutab.300.txt; 4,505+1 colnames
wc -l ../otutab.txt
# 4821 ../otutab.txt 4,820+1 colnames
# OTUS finalize ============================================================


# Find NUMTS
# R: temp3_cutadapt.Rproj --> otus_numts.Rmd ==> produced "otus97.final.300.nonumts.fa"




