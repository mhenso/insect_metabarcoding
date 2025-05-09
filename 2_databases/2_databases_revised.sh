# This db5_results3 was prepared for independent pipeline for each region. It removed the NUMTS and applied MAPUC and it extract from 393 to 711.

pwd
# /Sihaloho_etal/1_bold/results
ls -1
# db5.global.fa
# db5.local.fa
# db5.regional.fa

grep -c "^>" db5.local.fa
# 60342
grep -c "^>" db5.regional.fa
# 113895
grep -c "^>" db5.global.fa
# 4241641

# Removing leading/trailing Ns
### Notes: This will not delete the entries
perl -pe '/^>/ ? print "\n" : chomp' db5.local.fa | tail -n +2 | sed -r '/^>/! s/N+$|^N+//g' > temp && mv temp db5.local.fa
perl -pe '/^>/ ? print "\n" : chomp' db5.regional.fa | tail -n +2 | sed -r '/^>/! s/N+$|^N+//g' > temp && mv temp db5.regional.fa
perl -pe '/^>/ ? print "\n" : chomp' db5.global.fa | tail -n +2 | sed -r '/^>/! s/N+$|^N+//g' > temp && mv temp db5.global.fa

pwd
# /Sihaloho_etal/2_databases/results
cp ../../1_bold/results/*.fa ./

# Assessing the homopolimer ACTG >=12 & the ambiguous Ns >=6, and save
~/seqkit locate -P -f ./hopol.fa ./db5.local.fa --threads 8 | sort -nk 1,1 > db5.local.hopol.fa
~/seqkit locate -P -f ./hopol.fa ./db5.regional.fa --threads 8 | sort -nk 1,1 > db5.regional.hopol.fa
~/seqkit locate -P -f ./hopol.fa ./db5.global.fa --threads 8 | sort -nk 1,1 > db5.global.hopol.fa

# create list of sequence of hopol to be removed
cat db5.local.hopol.fa | awk '!seen[$1]++' | cut -f 1 > db5.local.hopol.fa.list
cat db5.regional.hopol.fa | awk '!seen[$1]++' | cut -f 1 > db5.regional.hopol.fa.list
cat db5.global.hopol.fa | awk '!seen[$1]++' | cut -f 1 > db5.global.hopol.fa.list

# remove the sequences
mkdir ./local ./regional ./global 

~/seqkit grep -v --pattern-file ./db5.local.hopol.fa.list ./db5.local.fa > ./local/db5.local.nohopol.fa
# [INFO] 116 patterns loaded from file
grep -c "^>" ./local/db5.local.nohopol.fa
# 60227
~/seqkit grep -v --pattern-file ./db5.regional.hopol.fa.list ./db5.regional.fa > ./regional/db5.regional.nohopol.fa
# [INFO] 468 patterns loaded from file
grep -c "^>" ./regional/db5.regional.nohopol.fa
# 113428
~/seqkit grep -v --pattern-file ./db5.global.hopol.fa.list ./db5.global.fa > ./global/db5.global.nohopol.fa
# [INFO] 12039 patterns loaded from file
grep -c "^>" ./global/db5.global.nohopol.fa
# 4229603

rm db5.global.fa db5.local.fa db5.regional.fa

# LOCAL -------------------------------------------------------------------------------------------------- #
pwd
# /Sihaloho_etal/2_databases/results/local

# use previously work
## align the rest of db5.local.nohopol.fa
## filter OUT the sequences to include only those NOT in the droplist file
wc -l ../droplist
# 2570 ../droplist

~/seqkit -w 0 grep -v --pattern-file ../droplist ./db5.local.nohopol.fa > ./db5.local.nohopol.rest.fa
# [INFO] 2570 patterns loaded from file
grep -c "^>" ./db5.local.nohopol.rest.fa
# 57657 --> 57657+2570 = 60,227 

## Align the rest db with previous aligned filtered db sequences
mafft --auto --addfull ./db5.local.nohopol.rest.fa --keeplength --thread 8 ../db.local2.nohopol.filt.afa.ler_fol.ref > db5.local.ler_fol.rest.afa &

grep -c "^>" db5.local.ler_fol.rest.afa 
# 60229

## add primer again
mafft --multipair --addfragments ../primers_ler_fol.fa --keeplength --thread 8 --mapout ./db5.local.ler_fol.rest.afa > ./db5.local.ler_fol.rest.ler_fol.afa
### Notes: 367-392, 712-737 ==> 711-393+1 = 319 it is correct as previous

## extract aligned region
python3 ./extract_alignment_region.py \
-i ./db5.local.ler_fol.rest.ler_fol.afa \
-o ./db5.local.ler_fol.rest.ler_fol.trimmed.afa \
-s 393 -e 711


## create a list for primer
grep "leray_f\|folmer_rc" --color -A 1 db5.local.ler_fol.rest.ler_fol.trimmed.afa
grep "leray_f\|folmer_rc" --color -A 1 db5.local.ler_fol.rest.ler_fol.trimmed.afa | grep "^>" | sed 's/>//' > drop.primers

## remove the primers
~/seqkit grep -v -n -f drop.primers db5.local.ler_fol.rest.ler_fol.trimmed.afa | ~/seqkit seq --upper-case -w 0 > db5.local.afa
grep -c "^>" db5.local.afa
# 60227

## see the frequencies of length
~/seqkit seq -g -w 0 ./db5.local.afa | grep -v "^>" | awk '{print length ($1)}' | sort | uniq -c | sort -nk 2,2 -r > ./db5.local.afa.tab


#--- local min 300bp --------------------------------
## create a list have seq length min 300
~/seqkit seq -g --min-len 300 -w 0 ./db5.local.afa | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 18809 list.min300

## obtain min300 length
~/seqkit grep -n -f list.min300 db5.local.afa | ~/seqkit seq --upper-case -w 0 > db5.local.min300.afa
# [INFO] 18809 patterns loaded from file

grep -c "^>"  db5.local.min300.afa
# 18809 ==> ### Notes: 18809/60227 * 100% = 31.23%

## see the frequencies of length again
~/seqkit seq -g -w 0 ./db5.local.min300.afa | grep -v "^>" | awk '{print length ($1)}' | sort | uniq -c | sort -nk 2,2 -r > ./db5.local.min300.afa.tab

## Ungap/remove the gaps from the local database
~/seqkit seq -g -w 0 ./db5.local.min300.afa > db5.local.min300.degap.fa
grep -c "^>" db5.local.min300.degap.fa
# 18809

# R check NUMTS --> db5.local.nonumts.fa
grep -c "^>" db5.local.nonumts.fa
# 17,888

## dereplicate ### Notes: using usearch
~/usearch11_64 --fastx_uniques db5.local.nonumts.fa -sizeout -fastaout ./db5.local.nonumts.uniq.fa -uc ./db5.local.nonumts.uniq.uc

grep -c "^>" ./db5.local.nonumts.uniq.fa
# 11,421

# R check MAPUC
# Produce "/results/local/db5.local.nonumts.uniq.mapuc.fa"

# BASH
pwd
# /Sihaloho_etal/2_databases/results/local
# convert to unix format and linearize
~/seqkit seq -w 0 ./db5.local.nonumts.uniq.mapuc.fa > db5.local.final.fa
grep -c "^>" db5.local.final.fa
# 11,372
#=== local min 300bp ================================

# Bash: blast
## Create blast database
pwd
# /results3/local
makeblastdb -in ./db5.local.final.fa -parse_seqids -dbtype nucl -out ../../blast_db/db5.local
# LOCAL ================================================================================================== # done


# REGIONAL --------------------------------------------------------------------------------------------------
pwd
# /Sihaloho_etal/2_databases/results/regional
grep -c "^>" db5.regional.nohopol.fa
# 113,428

### Notes: Use local results
## align the rest of db5.regional.nohopol.fa
## filter OUT the sequences to include only those NOT in the droplist file
seqkit -w 0 grep -v --pattern-file ../droplist ./db5.regional.nohopol.fa > ./db5.regional.nohopol.rest.fa
# [INFO] 2570 patterns loaded from file

grep -c "^>" db5.regional.nohopol.rest.fa
# 110,858 --> 110,858 + 2570 = 113,428 (same as on line 187)

## Align the rest dbregional with previous aligned filtered db sequences
# HPCC
mafft --auto --addfull ./db5.regional.nohopol.rest.fa --keeplength --thread 16 ../db.local2.nohopol.filt.afa.ler_fol.ref > db5.regional.nohopol.ler_fol.rest.afa &

grep -c "^>" db5.regional.nohopol.ler_fol.rest.afa
# 113430

## add primer again by subsetting
grep "^>" ./db5.regional.nohopol.ler_fol.rest.afa | head -n 2000 | sed 's/^>//' > subset_giant_regional.list
seqkit grep -n -f subset_giant_regional.list ./db5.regional.nohopol.ler_fol.rest.afa | seqkit seq --upper-case -w 0 > ./db5.regional.nohopol.ler_fol.rest.subset.afa
### align again
mafft --multipair --addfragments primers_ler_fol.fa --keeplength --thread 32 --mapout ./db5.regional.nohopol.ler_fol.rest.subset.afa > ./db5.regional.nohopol.ler_fol.rest.subset.ler_fol.afa
### Notes: 367-392, 712-737 ==> 711-393+1 = 319 it is correct as previous

# local_pc
## extract aligned region
python3 ./extract_alignment_region.py \
-i ./db5.regional.nohopol.ler_fol.rest.afa \
-o ./db5.regional.nohopol.ler_fol.rest.trimmed.afa \
-s 393 -e 711


## create a list for primer
grep "leray_f\|folmer_rc" --color -A 1 ./db5.regional.nohopol.ler_fol.rest.trimmed.afa
grep "leray_f\|folmer_rc" ./db5.regional.nohopol.ler_fol.rest.trimmed.afa | grep "^>" | sed 's/>//' > drop.primers

## remove the primers
~/seqkit grep -v -n -f drop.primers ./db5.regional.nohopol.ler_fol.rest.trimmed.afa | ~/seqkit seq --upper-case -w 0 > ./db5.regional.afa
# [INFO] 2 patterns loaded from file
grep -c "^>" ./db5.regional.afa
# 113,428
### Notes: the primers are removed, the same number as line 44

## see the frequencies of length
~/seqkit seq -g -w 0 ./db5.regional.afa | grep -v "^>" | awk '{print length ($1)}' | sort | uniq -c | sort -nk 2,2 -r > ./db5.regional.afa.tab

#--- regional min 300bp -----------------------------
## create a list have seq length min 300
~/seqkit seq -g --min-len 300 -w 0 ./db5.regional.afa | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 39461 list.min300

## obtain min300 length
~/seqkit grep -n -f list.min300 db5.regional.afa | seqkit seq --upper-case -w 0 > db5.regional.min300.afa
# [INFO] 39461 patterns loaded from file

grep -c "^>" db5.regional.min300.afa
# 39461 ==> ### Notes: 39,461/113,428 * 100% = 34.34.79%

## see the frequencies of length again
~/seqkit seq -g -w 0 ./db5.regional.min300.afa | grep -v "^>" | awk '{print length ($1)}' | sort | uniq -c | sort -nk 2,2 -r > ./db5.regional.min300.afa.tab

## Ungap/remove the gaps from the local database
~/seqkit seq -g -w 0 ./db5.regional.min300.afa > db5.regional.min300.degap.fa
grep -c "^>" db5.regional.min300.degap.fa
# 39,461

# R check NUMTS --> db5.regional.nonumts.fa
grep -c "^>" db5.regional.nonumts.fa
# 38,256

# usearch
## dereplicate ### Notes: using usearch
~/usearch11_64 --fastx_uniques db5.regional.nonumts.fa -sizeout -fastaout ./db5.regional.nonumts.uniq.fa -uc ./db5.regional.nonumts.uniq.uc

grep -c "^>" ./db5.regional.nonumts.uniq.fa
# 23061

# R check MAPUC
# Produce "/results/regional/db5.regional.nonumts.uniq.mapuc.fa"

# BASH
pwd
# /Sihaloho_etal/2_databases/results/regional
# convert to unix format and linearize
~/seqkit seq -w 0 ./db5.regional.nonumts.uniq.mapuc.fa > db5.regional.final.fa
grep -c "^>"  db5.regional.final.fa
# 22,738
#=== regional min 300bp =============================

# Bash: blast
## Create regional blast database
pwd
# /results3/regional/
makeblastdb -in ./db5.regional.final.fa -parse_seqids -dbtype nucl -out ../../blast_db/db5.regional
# REGIONAL ==================================================================================================


# GLOBAL ----------------------------------------------------------------------------------------------------
pwd
# /Sihaloho_etal/2_databases/results/global
grep -c "^>" db5.global.nohopol.fa
# 4229603

### Notes: Use local results
## align the rest of db5.regional.nohopol.fa
## filter OUT the sequences to include only those NOT in the droplist file
seqkit -w 0 grep -v --pattern-file ../droplist ./db5.global.nohopol.fa > ./db5.global.nohopol.rest.fa
# [INFO] 2570 patterns loaded from file

grep -c "^>" ./db5.global.nohopol.rest.fa
# 4227033 --> 4,227,033 + 2570 = 4229603 (same as on line 294)

## Align the rest dbregional with previous aligned filtered db sequences
mafft --auto --addfull ./db5.global.nohopol.rest.fa --keeplength --thread 16 ../db.local2.nohopol.filt.afa.ler_fol.ref > db5.global.nohopol.ler_fol.rest.afa 

## add primer again by subsetting
grep "^>" ./db5.global.nohopol.ler_fol.rest.afa | head -n 2000 | sed 's/^>//' > subset_giant_global.list
seqkit grep -n -f subset_giant_global.list ./db5.global.nohopol.ler_fol.rest.afa | seqkit seq --upper-case -w 0 > ./db5.global.nohopol.ler_fol.rest.subset.afa
### align again
mafft --multipair --addfragments primers_ler_fol.fa --keeplength --thread 32 --mapout ./db5.global.nohopol.ler_fol.rest.subset.afa > ./db5.global.nohopol.ler_fol.rest.subset.ler_fol.afa
### Notes: 367-392, 712-737 ==> 711-393+1 = 319 it is correct as previous

# local_pc
python3 ./extract_alignment_region.py \
-i ./db5.global.nohopol.ler_fol.rest.afa \
-o ./db5.global.nohopol.ler_fol.rest.trimmed.afa \
-s 393 -e 711

## create a list for primer
grep "leray_f\|folmer_rc" --color -A 1 ./db5.global.nohopol.ler_fol.rest.trimmed.afa
grep "leray_f\|folmer_rc" ./db5.global.nohopol.ler_fol.rest.trimmed.afa | grep "^>" | sed 's/>//' > drop.primers

## remove the primers
~/seqkit grep -v -n -f drop.primers ./db5.global.nohopol.ler_fol.rest.trimmed.afa | ~/seqkit seq --upper-case -w 0 > ./db5.global.afa
# [INFO] 2 patterns loaded from file
grep -c "^>" ./db5.global.afa
# 4,229,603
### Notes: the primers are removed, the same number as line 48

## see the frequencies of length
~/seqkit seq -g -w 0 ./db5.global.afa | grep -v "^>" | awk '{print length ($1)}' | sort | uniq -c | sort -nk 2,2 -r > ./db5.global.afa.tab


#--- regional min 300bp -----------------------------
## create a list have seq length min 300
~/seqkit seq -g --min-len 300 -w 0 ./db5.global.afa | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 2484568 list.min300

## obtain min300 length
~/seqkit grep -n -f list.min300 db5.global.afa | seqkit seq --upper-case -w 0 > db5.global.min300.afa
# [INFO] 2484568 patterns loaded from file

grep -c "^>" db5.global.min300.afa
# 2484568 ==> ### Notes: 2,484,568/4,229,603 * 100% = 58.74%

## see the frequencies of length again
~/seqkit seq -g -w 0 ./db5.global.min300.afa | grep -v "^>" | awk '{print length ($1)}' | sort | uniq -c | sort -nk 2,2 -r > ./db5.global.min300.afa.tab

## Ungap/remove the gaps from the local database
~/seqkit seq -g -w 0 ./db5.global.min300.afa > db5.global.min300.degap.fa
grep -c "^>" db5.global.min300.degap.fa
# 2,484,568

# R check NUMTS --> db5.global.nonumts.fa
grep -c "^>" db5.global.nonumts.fa
# 2,453,660

# usearch
## dereplicate ### Notes: using usearch
~/usearch11_64 --fastx_uniques db5.global.nonumts.fa -sizeout -fastaout ./db5.global.nonumts.uniq.fa -uc ./db5.global.nonumts.uniq.uc

grep -c "^>" ./db5.global.nonumts.uniq.fa
# 745,041

# R check MAPUC
# Produce "/results/global/db5.global.nonumts.uniq.mapuc.fa"

# BASH
pwd
# /Sihaloho_etal/2_databases/results/global
# convert to unix format and linearize
~/seqkit seq -w 0 ./db5.global.nonumts.uniq.mapuc.fa > db5.global.final.fa
grep -c "^>"  db5.global.final.fa
# 726,899
#=== regional min 300bp =============================

# Bash: blast
## Create regional blast database
pwd
# /results3/global
makeblastdb -in ./db5.global.final.fa -parse_seqids -dbtype nucl -out ../../blast_db/db5.global
# GLOBAL ====================================================================================================

