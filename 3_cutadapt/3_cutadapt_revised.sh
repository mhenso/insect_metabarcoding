# Sihaloho etal Optimizing Insect Metabarcoding
# OTU pipeline 
# Several softwares that we used are cutadapt


# create sample list file
ls -1 COI_252_samples/ | sed -E 's/_L001_R1_001.fastq|_L001_R2_001.fastq//' | uniq > samples

# CUTADAPT
pwd
# /COI_LT

for x in $(cat samples);
do
# forward
cutadapt -j 64 -g GGWACWGGWTGAACWGTWTAYCCYCC ./COI_252_samples/${x}_L001_R1_001.fastq -o ./cutadapt/lt_db/${x}_L001_R1_001.fastq
# reverse
cutadapt -j 64 -g TAAACTTCAGGGTGACCAAAAAATCA ./COI_252_samples/${x}_L001_R2_001.fastq -o ./cutadapt/lt_db/${x}_L001_R2_001.fastq;
done

# RAW READS
for x in $(cat samples);
do
cat ./lt_db/${x}_L001_R1_001.fastq | wc -l |  awk '{s+=$1} END {print s/4}' >> list_reads.txt;
done &
paste samples list_reads.txt | column -s $'\t' -t > list_reads_merged.txt


pwd 
# /COI_LT/cutadapt
cp ../samples ./

# PEAR
## stitch consensus read pairs for each sample
for x in $(cat samples);
do
pear --threads 64 --memory 100G \
-f lt_db/${x}_L001_R1_001.fastq \
-r lt_db/${x}_L001_R2_001.fastq \
-o stitched/${x}.fastq;
done


# USEARCH
## -fastq_filter filter consensus reads to only include consensus reads with estimated 1 error max
usearch='/home/caphilli/usearch11_64'
export OMP_NUM_THREADS=64

for x in $(cat samples);
do
$usearch -fastq_filter stitched/${x}.fastq.assembled.fastq \
-fastaout filtered/${x}_filtered.fa \
-fastq_maxee 1.0 -threads 64;
done

# BASH
## relabel the reads so we can keep track after pooling of which sample each came from
for x in $(cat samples);
do
sed "-es/^>\(.*\)/>\1;barcodelabel=${x};/" < filtered/${x}_filtered.fa > filtered/${x}.fa;
done

## pool the relabeled concensus sequences from each sample
mkdir output
cat filtered/*{0..9}.fa > ./output/combined.fa
cat output/combined.fa | grep -c "^>"
# 1,613,781

# USEARCH
## fastx_uniques make a fasta with only unique sequences but keeping track of abundances
$usearch -fastx_uniques ./output/combined.fa -fastaout ./output/uniques.fa -sizeout -relabel Uniq -threads 64
cat ./output/uniques.fa | grep -c "^>"
# 423,545

## -cluster_otus to identify otus
$usearch -cluster_otus ./output/uniques.fa -otus ./output/otus97.fa -relabel Otu -threads 64
cat ./output/otus97.fa | grep -c "^>"
# 3310

## -otutab build the otu table from the identified zotu and their counts across samples
usearch='/home/caphilli/usearch11_64'
$usearch -threads 128 \
-otutab ./output/combined.fa \
-otus ./output/otus97.fa \
-otutabout ./output/otutab.txt \
-mapout ./output/otumap.txt 

# we got our otu_table, we need to finalize the otu table
pwd 
# /COI_LT/cutadapt/output/otus_finalize
# OTUS finalize 
seqkit seq -g -w 0 ../otus97.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > otus97.fa.tab

## subsample 200 id
grep "^>" ../otus97.fa | head -n 200 | sed 's/^>//' > otus200.list
seqkit grep -n -f ./otus200.list ../otus97.fa > ./otus97.200.fa
## align the subsample otus
mafft --auto --thread 8 ./otus97.200.fa > ./otus97.200.afa


## filter OUT otus that have been aligned (200otus)
seqkit grep -v -n -f ./otus200.list ../otus97.fa > ./otus97rest.fa
grep -c "^>" otus97rest.fa
# 3,110  ### Notes: the number is correct

## align rest otus
mafft --auto --addfull ./otus97rest.fa --keeplength --thread 16 otus97.200.afa > otus97.200.rest.afa

grep -c "^>" otus97.200.rest.afa
# 3310

## see the frequencies of length
seqkit seq -g -w 0 otus97.200.rest.afa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > otus97.200.rest.afa.tab


## obtain min 300 length
seqkit seq -g --min-len 300 -w 0 ./otus97.200.rest.afa  | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 3247 list.min300

## obtain min300 length
seqkit grep -n -f list.min300 otus97.200.rest.afa | seqkit seq --upper-case -w 0 > otus97.200.rest.min300.afa
# [INFO] 3247 patterns loaded from file
grep -c "^>" otus97.200.rest.min300.afa
# 3247 ==> ### Notes: 3247/3310 * 100% = 98.09%

## see the frequencies of length again
seqkit seq -g -w 0 otus97.200.rest.min300.afa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > otus97.200.rest.min300.afa.tab

## Ungap/remove the gaps from the local database
seqkit seq -g ./otus97.200.rest.min300.afa > otus97.final.300.fa
grep -c "^>" otus97.final.300.fa
# 3247

# retrive otu_table with the final id
grep "^>" otus97.final.300.fa | sed 's/^>//' > otus97.final.300.id
sed -i "1i #OTU" otus97.final.300.id
grep -w -f otus97.final.300.id ../otutab.txt > ./otutab.300.txt

# So, our final files
## otutab file is otutab.300.txt, and 
## refseq are otus97.final.300.fa (ungapped) and otus97.200.rest.min300.afa (gapped)


# Find NUMTS
# R: temp3_cutadapt.Rproj --> otus_numts.Rmd ==> produced "otus97.final.300.nonumts.fa"

~/seqkit seq -w 0 ./otus97.final.300.nonumts.fa > temp && mv temp ./otus97.final.300.nonumts.fa


# R: temp3_cutadapt.Rproj --> otus300.Rmd phyloseq, srs, inext plotting ==> produced Fig.4 & Fig. 5

























