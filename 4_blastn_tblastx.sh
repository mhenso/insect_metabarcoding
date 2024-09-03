# BASH: BLASTN
# "otus97.final.300.fa" blast to "min300" databases

pwd
# /mnt/c/Docs/R/databases/results_blastn300
grep -c "^>" ../../cutadapt/output/otus_finalize/otus97.final.300.fa
# 3247


## local
blastn -num_threads 10 \
-query ../../cutadapt/output/otus_finalize/otus97.final.300.fa \
-db ../results/local/blast300/mydb300 -evalue 1e-4 \
-out ./local300.out \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" &

wc -l local300.out
# 965919 local300.out
cat ./local300.out | sort -nk 12 -r -nk 3 | awk '!seen[$1]++' > local.out.top
wc -l local.out.top
# 3197 ### Notes: 50 otus were not assigned on local.db


## regional
blastn -num_threads 10 \
-query ../../cutadapt/output/otus_finalize/otus97.final.300.fa \
-db ../results/regional/blast300/imdb300 -evalue 1e-4 \
-out ./regional300.out \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" &

wc -l regional300.out
# 1296483 regional300.out
cat ./regional300.out | sort -nk 12 -r -nk 3 | awk '!seen[$1]++' > regional.out.top
wc -l regional.out.top
# 3230 ### Notes: 17 Otus were not assigned on regional.db


## global
blastn -num_threads 10 \
-query ../../cutadapt/output/otus_finalize/otus97.final.300.fa \
-db ../results/global/blast300/global300 -evalue 1e-4 \
-out ./global300.out \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" &

wc -l global.out
# 1601378 global300.out
cat global300.out | sort -nk 12 -r -nk 3 | awk '!seen[$1]++' > global.out.top
wc -l global.out.top
# 3243 global.out.top


# BASH: TBLASTX
# "otus97.final.300.fa" tblastx to "min300" databases

pwd
# /mnt/c/Docs/R/databases/results_blastn300
grep -c "^>" ../../cutadapt/output/otus_finalize/otus97.final.300.fa
# 3247

## local
tblastx -num_threads 10 \
-query ../../cutadapt/output/otus_finalize/otus97.final.300.fa \
-db ../results/local/blast300/mydb300 -evalue 1e-4 \
-out ./local300.out \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" &

wc -l local300.out
# 10,785,979 local300.out
cat ./local300.out | sort -nk 12 -r -nk 3 | awk '!seen[$1]++' > local.out.top
wc -l local.out.top
# 3247 local.out.top


## regional
tblastx -num_threads 10 \
-query ../../cutadapt/output/otus_finalize/otus97.final.300.fa \
-db ../results/regional/blast300/imdb300 -evalue 1e-4 \
-out ./regional300.out \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" &

cat regional300.out | wc -l
# 10828462
cat ./regional300.out | sort -nk 12 -r -nk 3 | awk '!seen[$1]++' > regional.out.top
wc -l regional.out.top
# 3247 regional.out.top


## global
tblastx -num_threads 10 \
-query ../../cutadapt/output/otus_finalize/otus97.final.300.fa \
-db ../results/global/blast300/global300 -evalue 1e-4 \
-out ./global300.out \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore salltitles" &

cat global300.out | wc -l 
# 11,467,601
cat global300.out | sort -nk 12 -r -nk 3 | awk '!seen[$1]++' > global.out.top
wc -l global.out.top
# 3247 global.out.top