pwd
# /Sihaloho_etal/2_databases/results
ls -1 ../../1_bold/results/
# db.global2.fa
# db.local2.fa
# db.regional2.fa

cat ../../1_bold/results/db.local2.fa | grep -c "^>"
# 7,262
cat ../../1_bold/results/db.regional2.fa | grep -c "^>"
# 29,418
cat ../../1_bold/results/db.global2.fa | grep -c "^>"
# 1,741,530



# Removing leading/trailing Ns
### Notes: This will not delete the entries
perl -pe '/^>/ ? print "\n" : chomp' ../data/db.local2.fa | tail -n +2 | sed -r '/^>/! s/N+$|^N+//g' > db.local2.fa
perl -pe '/^>/ ? print "\n" : chomp' ../data/db.regional2.fa | tail -n +2 | sed -r '/^>/! s/N+$|^N+//g' > db.regional2.fa
perl -pe '/^>/ ? print "\n" : chomp' ../data/db.global2.fa | tail -n +2 | sed -r '/^>/! s/N+$|^N+//g' > db.global2.fa


# Assessing the homopolimer ACTG >=12, just check
seqkit locate -P -f ./hopolACTG.fa ./db.local2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 3
seqkit locate -P -f ./hopolACTG.fa ./db.regional2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 21
seqkit locate -P -f ./hopolACTG.fa ./db.global2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 432

# Assessing the ambiguous Ns >=6, just check
seqkit locate -P -f ./hopolNs.fa ./db.local2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 107
seqkit locate -P -f ./hopolNs.fa ./db.regional2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 432
seqkit locate -P -f ./hopolNs.fa ./db.global2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 10,443

# Assessing the homopolimer ACTG >=12 & the ambiguous Ns >=6, just check
seqkit locate -P -f ./hopol.fa ./db.local2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 109
seqkit locate -P -f ./hopol.fa ./db.regional2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 452
seqkit locate -P -f ./hopol.fa ./db.global2.fa --threads 8 | sort -nk 1,1 | awk '!seen[$1]++' | wc -l
# 10869

# Assessing the homopolimer ACTG >=12 & the ambiguous Ns >=6, and save
seqkit locate -P -f ./hopol.fa ./db.local2.fa --threads 8 | sort -nk 1,1 > db.local2.hopol.fa
seqkit locate -P -f ./hopol.fa ./db.regional2.fa --threads 8 | sort -nk 1,1 > db.regional2.hopol.fa
seqkit locate -P -f ./hopol.fa ./db.global2.fa --threads 8 | sort -nk 1,1 > db.global2.hopol.fa


# create list of sequence of hopol to be removed
wc -l db.local2.hopol.fa
# 3935 db.local2.hopol.fa
cat db.local2.hopol.fa | awk '!seen[$1]++' | wc -l
# 109

cat db.local2.hopol.fa | awk '!seen[$1]++' | cut -f 1 > db.local2.hopol.fa.list
cat db.regional2.hopol.fa | awk '!seen[$1]++' | cut -f 1 > db.regional2.hopol.fa.list
cat db.global2.hopol.fa | awk '!seen[$1]++' | cut -f 1 > db.global2.hopol.fa.list

# remove the sequences
mkdir ./local
mkdir ./regional
mkdir ./global

seqkit grep -v --pattern-file ./db.local2.hopol.fa.list ./db.local2.fa > ./local/db.local2.nohopol.fa
# [INFO] 109 patterns loaded from file
seqkit grep -v --pattern-file ./db.regional2.hopol.fa.list ./db.regional2.fa > ./regional/db.regional2.nohopol.fa
# [INFO] 452 patterns loaded from file
seqkit grep -v --pattern-file ./db.global2.hopol.fa.list ./db.global2.fa > ./global/db.global2.nohopol.fa
# [INFO] 10869 patterns loaded from file

# LOCAL -------------------------------------------------------------------------------------------------- #
pwd
# /Sihaloho_etal/2_databases/results/local
grep -c "^>" db.local2.nohopol.fa
# 7154
seqkit seq -g -w 0 db.local2.nohopol.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.local2.nohopol.fa.tab

## subset 655-1500bp
seqkit seq --min-len 655 --max-len 1500 -w 0 ./db.local2.nohopol.fa > ./db.local2.nohopol.filt.fa
grep -c ^">" db.local2.nohopol.filt.fa
# 2570

## align
mafft --auto --thread 8 ./db.local2.nohopol.filt.fa > ./db.local2.nohopol.filt.afa

## align the primers to db.local2.nohopol.filt.afa
mafft --multipair --addfragments primers_ler_fol.fa --keeplength --thread 8 --mapout --reorder db.local2.nohopol.filt.afa > db.local2.nohopol.filt.afa.ler_fol.ref
### Notes: 367-392, 712-737 ==> 711-393+1 = 319

## align the rest of db.local2.nohopol.fa
## create a droplist from 
grep "^>" ./db.local2.nohopol.filt.fa | sed 's/^>//' > droplist
## filter OUT the sequences to include only those NOT in the droplist file
seqkit grep -v --pattern-file ./droplist ./db.local2.nohopol.fa > ./db.local2.nohopol.rest.fa
grep -c "^>" db.local2.nohopol.rest.fa
# 4584
### Notes: Therea are 4,584 (7154-2570) entries that have < 655 nucleotides

## Align the rest db with previous aligned filtered db sequences
mafft --auto --addfull ./db.local2.nohopol.rest.fa --keeplength --thread 8 db.local2.nohopol.filt.afa.ler_fol.ref > db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa

grep -c "^>" ./db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa
# 7156
### Notes : It has 2 more entries (primers) than db.local2.nohopol.fa, we're on track

## aligned with primers_ler_fol again, be carefull here, cos you may replace previous map file
mv primers_ler_fol.fa.map primers_ler_fol.fa.map.1st

mafft --multipair --addfragments primers_ler_fol.fa --keeplength --thread 8 --mapout db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa > db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref
### Notes: 367-392, 712-737 ==> 711-393 = 319, it is the same as before

## extract aligned region
python3 ./extract_alignment_region.py \
-i ./db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref \
-o ./db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa \
-s 393 -e 711

grep -c "^>" db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa
# 7158 ### Notes: added 2 primers

## create a list for primer
grep "leray_f\|folmer_rc" --color -A 1 db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa
grep "leray_f\|folmer_rc" db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa | grep "^>" | sed 's/>//' > drop.primers

## remove the primers
seqkit grep -v -n -f drop.primers db.local2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa | seqkit seq --upper-case -w 0 > db.local2.final.fa
grep -c "^>" db.local2.final.fa
# 7154
### Notes: the primers are removed


## see the frequencies of length
seqkit seq -g -w 0 db.local2.final.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.local2.final.fa.tab

#--- local min 300bp --------------------------------
## create a list have seq length min 300
seqkit seq -g --min-len 300 -w 0 ./db.local2.final.fa | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 3733 list.min300

## obtain min300 length
seqkit grep -n -f list.min300 db.local2.final.fa | seqkit seq --upper-case -w 0 > db.local2.final.min300.afa
# [INFO] 3733 patterns loaded from file

grep -c "^>" db.local2.final.min300.afa
# 6888 ==> ### Notes: 3733/7154 * 100% = 52.18%

## see the frequencies of length again
seqkit seq -g -w 0 db.local2.final.min300.afa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.local2.final.min300.afa.tab

## Ungap/remove the gaps from the local database
seqkit seq -g ./db.local2.final.min300.afa > db.local2.final.min300.ungap.afa
grep -c "^>" db.local2.final.min300.ungap.afa
# 3733

## dereplicate ### Notes: using usearch
usearch11 --fastx_uniques db.local2.final.min300.ungap.afa -sizeout \
-fastaout ./usearch/db.local2.final.min300.ungap.unique.fa -uc ./usearch/local_map300.uc 

grep -c "^>" ./usearch/db.local2.final.min300.ungap.unique.fa
# 2468 ==> 2468/7154=34.49%

# R: "/Sihaloho_etal/2_databases/mapuc_300.rmd" ==> # MAPUC & # Dropping_sequences
# Produce "/Sihaloho_etal/2_databases/results/local/usearch/db.local2.final.min300.ungap.unique.clean.fa"

# BASH
pwd
# /Sihaloho_etal/2_databases/results/local/usearch
# convert to unix format and linearize
seqkit seq -w 0 ./db.local2.final.min300.ungap.unique.clean.fa > mydb300.fa
grep -c "^>" mydb300.fa
# 2439
#=== local min 300bp ================================


# Bash: blast
## Create blast database
pwd
# /Sihaloho_etal/2_databases/results/local/usearch
makeblastdb -in ./mydb300.fa -parse_seqids -dbtype nucl -out ../blast300/mydb300
# LOCAL ================================================================================================== #


# REGIONAL --------------------------------------------------------------------------------------------------
pwd
# /Sihaloho_etal/2_databases/results/regional
grep -c "^>" db.regional2.nohopol.fa
# 28967

seqkit seq -g -w 0 db.local2.nohopol.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.local2.nohopol.fa.tab


### Notes: Use local results

## align the rest of db.regional.nohopol.fa
## copy a droplist from local
cp ../local/droplist ./
wc -l droplist
# 2570 droplist ### Notes: This filtered db.local2.

## copy filtered aligned reference from local
cp ../local/db.local2.nohopol.filt.afa.ler_fol.ref ./

## filter OUT the sequences to include only those NOT in the droplist file
seqkit grep -v --pattern-file ./droplist ./db.regional2.nohopol.fa > ./db.regional2.nohopol.rest.fa
# [INFO] 2570 patterns loaded from file
grep -c "^>" db.regional2.nohopol.rest.fa
# 26397
### Notes: Therea are 26,397 (28,967-2,570) entries 

## Align the rest db with previous aligned filtered db sequences
mafft --auto --addfull ./db.regional2.nohopol.rest.fa --keeplength --thread 8 db.local2.nohopol.filt.afa.ler_fol.ref > db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa

grep -c "^>" ./db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa
# 28,969
### Notes : It has 2 more entries (primers) than db.regional2.nohopol.fa, we're on track

## aligned with primers_ler_fol again, be carefull here, cos you may replace previous map file, not needed 
mafft --multipair --addfragments primers_ler_fol.fa --keeplength --thread 8 --mapout db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa > db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref
### Notes: Check primers_ler_fol.fa.map
### Notes:367-392, 712-737 ==> 711-393 = 319, it is the same as before on local, we're on good track

grep -c "^>" db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref
# 28,971

## extract aligned region
cp ../local/extract_alignment_region.py ./

python3 ./extract_alignment_region.py \
-i ./db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref \
-o ./db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa \
-s 393 -e 711

grep -c "^>" db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa
# 28971 ### Notes: added 2 primers

## create a list for primer
grep "leray_f\|folmer_rc" --color -A 1 db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa
grep "leray_f\|folmer_rc" db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa | grep "^>" | sed 's/>//' > drop.primers

## remove the primers
seqkit grep -v -n -f drop.primers db.regional2.nohopol.filt.afa.ler_fol.ref.rest.afa.ref.afa | seqkit seq --upper-case -w 0 > db.regional2.final.fa

grep -c "^>" db.regional2.final.fa
# 28967
### Notes: the primers are removed, the same number as db.regional2.nohopol.fa


## see the frequencies of length
seqkit seq -g -w 0 db.regional2.final.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.regional2.final.fa.tab

#--- regional min 300bp -----------------------------
pwd
# /mnt/c/Docs/R/databases/results/regional/usearch
## create a list have seq length min 300
seqkit seq -g --min-len 300 -w 0 ./db.regional2.final.fa | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 18951 list.min300

## obtain min300 length
seqkit grep -n -f list.min300 db.regional2.final.fa | seqkit seq --upper-case -w 0 > db.regional2.final.min300.afa
# [INFO] 18951 patterns loaded from file

grep -c "^>" db.regional2.final.min300.afa
# 18951 ==> ### Notes: 18951/28967 * 100% = 65.42

## see the frequencies of length again
seqkit seq -g -w 0 db.regional2.final.min300.afa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.regional2.final.min300.afa.tab

## Ungap/remove the gaps from the local database
seqkit seq -g ./db.regional2.final.min300.afa > db.regional2.final.min300.ungap.afa

grep -c "^>" db.regional2.final.min300.ungap.afa
# 18951

## dereplicate ### Notes: using usearch
usearch11 --fastx_uniques db.regional2.final.min300.ungap.afa -sizeout \
-fastaout ./usearch/db.regional2.final.min300.ungap.unique.fa -uc ./usearch/regional_map300.uc 

grep -c "^>" ./usearch/db.regional2.final.min300.ungap.unique.fa
# 11,043


# R: "/Sihaloho_etal/2_databases/mapuc_300.rmd" ==> # MAPUC & # Dropping_sequences
# Produce: "/Sihaloho_etal/2_databases/results/regional/usearch/db.regional2.final.min300.ungap.unique.clean.fa"

# BASH
pwd
# /Sihaloho_etal/2_databases/results/regional/usearch
# convert to unix format and linearize
seqkit seq -w 0 ./db.regional2.final.min300.ungap.unique.clean.fa > imdb300.fa
grep -c "^>" imdb300.fa
# 10,830
#=== regional min 300bp =============================

# Bash: blast
## Create regional blast database
pwd
# /Sihaloho_etal/2_databases/results/regional/usearch
makeblastdb -in ./imdb300.fa -parse_seqids -dbtype nucl -out ../blast300/imdb300
# REGIONAL ==================================================================================================


# GLOBAL ----------------------------------------------------------------------------------------------------
pwd
# /Sihaloho_etal/2_databases/results/global
grep -c "^>" db.global2.nohopol.fa
# 1730662

seqkit seq -g -w 0 db.global2.nohopol.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.global2.nohopol.fa.tab

### Notes: Use local results
## align the rest of db.global.nohopol.fa
## copy a droplist from local
cp ../local/droplist ./
wc -l droplist
# 2570 droplist ### Notes: This filtered db.local2.

## copy filtered aligned reference from local
cp ../local/db.local2.nohopol.filt.afa.ler_fol.ref ./

## filter OUT the sequences to include only those NOT in the droplist file
seqkit grep -v --pattern-file ./droplist ./db.global2.nohopol.fa > ./db.global2.nohopol.rest.fa
# [INFO] 2570 patterns loaded from file
grep -c "^>" db.global2.nohopol.rest.fa
# 1728092
### Notes: Therea are 1,728,092 (1730662-2,570) entries 

## Align the rest db with previous aligned filtered db sequences
mafft --auto --addfull ./db.global2.nohopol.rest.fa \
--keeplength --thread 16 \
db.local2.nohopol.filt.afa.ler_fol.ref > db.global2.nohopol.filt.afa.ler_fol.ref.rest.afa

grep -c "^>" ./db.global2.nohopol.filt.afa.ler_fol.ref.rest.afa
# 1,730,664
### Notes : It has 2 more entries (primers) than db.regional2.nohopol.fa, we're on track


## aligned with primers_ler_fol again, be carefull here, cos you may replace previous map file, not needed 
### Notes: We will run out of memory by adding primers to giant alignment, so we need to subset the giant alignment then add primers

# ---- confirmation of the coordinates with the subset of giant alignment ----
cd confirm_coordinate/
cp ../primers_ler_fol.fa ./

grep "^>" ../db.global2.nohopol.filt.afa.ler_fol.ref.rest.afa | head -n 2000 | sed 's/^>//' > subset_giant.list

seqkit grep -n -f subset_giant.list ../db.global2.nohopol.filt.afa.ler_fol.ref.rest.afa | seqkit seq --upper-case -w 0 > ./db.global2.nohopol.filt.afa.ler_fol.ref.rest.2000.afa

mafft --multipair --addfragments primers_ler_fol.fa --keeplength --thread 128 --mapout \
./db.global2.nohopol.filt.afa.ler_fol.ref.rest.2000.afa > ./db.global2.nohopol.filt.afa.ler_fol.ref.rest.2000.afa.ref

### Notes: Check primers_ler_fol.fa.map
### Notes:367-392, 712-737 ==> 711-393 = 319, it is the same as before on local, we're on good track
# ==== confirmation of the coordinates with the subset of giant alignment ====

### Notes: Since the position are the same as we did with db.local / db.regional, we can continue extract the aligned region from the giant alignment.

pwd
# /Sihaloho_etal/2_databases/results/global
## extract aligned region
cp ../local/extract_alignment_region.py ./

python3 ./extract_alignment_region.py \
-i ./db.global2.nohopol.filt.afa.ler_fol.ref.rest.afa \
-o ./db.global2.nohopol.filt.afa.ler_fol.ref.rest.extract.afa \
-s 393 -e 711


grep -c "^>" ./db.global2.nohopol.filt.afa.ler_fol.ref.rest.extract.afa
# 1,730,664 ### Notes: added 2 primers

## create a list for primer
grep "leray_f\|folmer_rc" --color -A 1 db.global2.nohopol.filt.afa.ler_fol.ref.rest.extract.afa
grep "leray_f\|folmer_rc" db.global2.nohopol.filt.afa.ler_fol.ref.rest.extract.afa | grep "^>" | sed 's/>//' > drop.primers

## remove the primers
seqkit grep -v -n -f drop.primers db.global2.nohopol.filt.afa.ler_fol.ref.rest.extract.afa | seqkit seq --upper-case -w 0 > db.global2.final.fa

grep -c "^>" db.global2.final.fa
# 1,730,662
### Notes: the primers are removed, the same number as db.global2.nohopol.fa


## see the frequencies of length
seqkit seq -g -w 0 db.global2.final.fa | paste - - | awk '{print length ($2)}' | sort | uniq -c | sort -nk 2,2 -r > db.global2.final.fa.tab

# --- global min 300bp ------------------------------
## create a list have seq length min 300
seqkit seq -g --min-len 300 -w 0 ./db.global2.final.fa | grep "^>" | sed 's/>//' > list.min300
wc -l list.min300
# 1,127,866 list.min300

## obtain min300 length
seqkit grep -n -f list.min300 db.global2.final.fa | seqkit seq --upper-case -w 0 > db.global2.final.min300.afa
# [INFO] 1127866 patterns loaded from file

grep -c "^>" db.global2.final.min300.afa
# 1127866 ==> ### Notes: 1127866/1730662 * 100% = 65.17

## see the frequencies of length again
grep -v "^>" db.global2.final.min300.afa | sed 's/-//g' | awk '{print length}' | sort | uniq -c | sort -nk 2,2 -r > db.global2.final.min300.afa.tab


## Ungap/remove the gaps from the local database
seqkit seq -g -w 0 ./db.global2.final.min300.afa  > db.global2.final.min300.ungap.afa
grep -c "^>" db.global2.final.min300.ungap.afa
# 1,127,866

## see the frequencies of length of ungap
grep -v "^>" db.global2.final.min300.ungap.afa | sed 's/-//g' | awk '{print length}' | sort | uniq -c | sort -nk 2,2 -r > db.global2.final.min300.ungap.afa.tab


## dereplicate ### Notes: using usearch
mkdir usearch
usearch11 --fastx_uniques db.global2.final.min300.ungap.afa -sizeout \
-fastaout ./usearch/db.global2.final.min300.ungap.unique.fa -uc ./usearch/global_map300.uc 

grep -c "^>" ./usearch/db.global2.final.min300.ungap.unique.fa
# 371,183


# R: "/Sihaloho_etal/2_databases/mapuc_300.rmd" ==> # MAPUC & # Dropping_sequences
# Produce: "/Sihaloho_etal/2_databases/results/regional/usearch/db.regional2.final.min300.ungap.unique.clean.fa"

# BASH
pwd
# /Sihaloho_etal/2_databases/results/regional/usearch
seqkit seq -w 0 ./db.global2.final.min300.ungap.unique.clean.fa > global300.fa
grep -c "^>" global300.fa
# 360,843
# === global min 300bp ==============================


# Bash: blast
## Create blast database
pwd
# /Sihaloho_etal/2_databases/results/regional/usearch
makeblastdb -in ./global300.fa -parse_seqids -dbtype nucl -out ../blast300/global300
# GLOBAL ====================================================================================================

