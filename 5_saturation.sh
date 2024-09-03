pwd
# /mnt/c/Docs/R/satur_global/data

# Create id list of de-replicated databases
cat ../../databases/results/local/usearch/mydb300.fa | grep "^>" | sed 's/\ /;/' | sed 's/\ /,/g' | sed 's/>//' | sed 's/,$/; /' > mydb_id
wc -l mydb_id
# 2439 mydb_id

cat ../../databases/results/regional/usearch/imdb300.fa | grep "^>" | sed 's/\ /;/' | sed 's/\ /,/g' | sed 's/>//' | sed 's/,$/; /' > imdb_id
wc -l imdb_id
# 10830 imdb_id


# obtain alignment that related to mydb_id
seqkit grep -n -f mydb_id ../../databases/results/local/db.local2.final.fa | seqkit seq --upper-case -w 0 > db.local2.satur.afa
# [INFO] 2439 patterns loaded from file
grep -c "^>" db.local2.satur.afa
# 2439
seqkit grep -n -f imdb_id ../../databases/results/regional/db.regional2.final.fa | seqkit seq --upper-case -w 0 > db.regional2.satur.afa
# [INFO] 10830 patterns loaded from file
grep -c "^>" db.regional2.satur.afa
# 10830

# Use DAMBE to generate the 1st and 3rd codon for local and regional and save the files on data/dambe folder
# R: satur_codon.rmd ==> beta coefficient for local and regional



# GLOBAL
pwd
# /mnt/c/Docs/R/satur_global/data/global

## Create id list of de-replicated databases
cat ../../../databases/results/global/usearch/global300.fa | grep "^>" | sed 's/\ /;/' | sed 's/\ /,/g' | sed 's/>//' | sed 's/,$/; /' > global_id
wc -l global_id
# 360843 global_id

## obtain alignment that related to global_id
### before obtaining, need to shorten global id so DAMBE will not crash
seqkit replace -p ";tax=.*" -r ";tax " -w 0 ../../../git_databases/results/global/db.global2.final.fa > db.global2.final_singleid.fa
sed 's/;tax.*/;tax /' global_id > global_id2

## obtaining using global_id2
seqkit grep -n -f global_id2 db.global2.final_singleid.fa | seqkit seq --upper-case -w 0 > db.global2.satur.afa
# [INFO] 360843 patterns loaded from file

grep -c "^>" db.global2.satur.afa
# 360843

## subsample 1st 10,000
seqkit sample --rand-seed 101 --number 10000 --two-pass db.global2.satur.afa -w 0 > global_101.afa

# reduced the db.global2.satur.afa > db.global2.satur_after101.afa before subsample 102
# --
grep "^>" global_101.afa | sed 's/>//' > global_101.id
## filter OUT the sequences to include only those NOT in the global_101.id
seqkit grep -v -n -f ./global_101.id ./db.global2.satur.afa | seqkit seq --upper-case -w 0 > ./db.global2.satur_after101.afa
# [INFO] 10000 patterns loaded from file
grep -c "^>" ./db.global2.satur_after101.afa
# 350843, it is 10,000 less than db.global2.satur.afa
# ==


## subsample 2nd 10,000
seqkit sample --rand-seed 102 --number 10000 --two-pass ./db.global2.satur_after101.afa -w 0 > global_102.afa

# reduced again
grep "^>" global_102.afa | sed 's/>//' > global_102.id
seqkit grep -v -n -f ./global_102.id ./db.global2.satur_after101.afa | seqkit seq --upper-case -w 0 > ./db.global2.satur_after102.afa
grep -c "^>" ./db.global2.satur_after102.afa
# 340843


## subsample 3rd 10,000
seqkit sample --rand-seed 103 --number 10000 --two-pass ./db.global2.satur_after102.afa -w 0 > global_103.afa

# Use DAMBE to generate the 1st and 3rd codon for global and save the files on data/global/dambe folder
# R: satur_codon.rmd ==> beta coefficient for local and regional




































