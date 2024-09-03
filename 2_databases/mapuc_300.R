## .scroll-200 {
##   max-height: 200px;
##   overflow-y: auto;
##   background-color: inherit;
## }

## .tocify .tocify-header {
##   #position: fixed;
##   #top: 50px;
##   #left: 50px;
##   width: 350px;
##   #height: 400px;
##   }

suppressPackageStartupMessages({
  library(tidyverse)
  library(reshape2)
  library(Biostrings)
})

mism <- function(x,sh) { # x is object taxa level, sh is df from map uc
  uc_mism =c()
  a = rev(x)[1]
  # Species
  temp = sh
  temp$V9_2 = str_replace(temp$V9_2, paste0(".*", a), "") %>% str_replace(";:.*", "")
  temp$V10_2 = str_replace(temp$V10_2, paste0(".*", a), "") %>% str_replace(";:.*", "")
  b = length(which(temp$V10_2 != temp$V9_2)) # mismatch dataframe for each taxa level
  c = sub(",", "", a)
  c = sub(":", "", c)
  d = data.frame("tax_lev"=c, "nos.mismatch" = b)
  uc_mism = rbind(uc_mism, d)
  uc_mism_s = temp[which(temp$V10_2 != temp$V9_2),]
  
  # Genus
  a = rev(x)[2]
  temp = sh
  temp$V9_2 = str_replace(temp$V9_2, paste0(".*", a), "") %>% str_replace(",s:.*", "")
  temp$V10_2 = str_replace(temp$V10_2, paste0(".*", a), "") %>% str_replace(",s:.*", "")
  b = length(which(temp$V10_2 != temp$V9_2))
  c = sub(",", "", a)
  c = sub(":", "", c)
  d = data.frame("tax_lev"=c, "nos.mismatch" = b)
  uc_mism = rbind(uc_mism, d)
  uc_mism_g = temp[which(temp$V10_2 != temp$V9_2), ]
  
  # Family
  a = rev(x)[3]
  temp = sh
  temp$V9_2 = str_replace(temp$V9_2, paste0(".*", a), "") %>% str_replace(",g:.*", "")
  temp$V10_2 = str_replace(temp$V10_2, paste0(".*", a), "") %>% str_replace(",g:.*", "")
  b = length(which(temp$V10_2 != temp$V9_2))
  c = sub(",", "", a)
  c = sub(":", "", c)
  d = data.frame("tax_lev"=c, "nos.mismatch" = b)
  uc_mism = rbind(uc_mism, d)
  uc_mism_f = temp[which(temp$V10_2 != temp$V9_2), ]
  
  # Order
  a = rev(x)[4]
  temp = sh
  temp$V9_2 = str_replace(temp$V9_2, paste0(".*", a), "") %>% str_replace(",f:.*", "")
  temp$V10_2 = str_replace(temp$V10_2, paste0(".*", a), "") %>% str_replace(",f:.*", "")
  b = length(which(temp$V10_2 != temp$V9_2))
  c = sub(",", "", a)
  c = sub(":", "", c)
  d = data.frame("tax_lev"=c, "nos.mismatch" = b)
  uc_mism = rbind(uc_mism, d)
  uc_mism_o = temp[which(temp$V10_2 != temp$V9_2), ]
  
  uc_mism_list = list("mism_df"      = uc_mism,
                      "mism_species" = uc_mism_s,
                      "mism_genus"   = uc_mism_g,
                      "mism_family"  = uc_mism_f,
                      "mism_order"   = uc_mism_o)

  print(uc_mism)
  return (uc_mism_list)
  # rm(uc_mism, temp, a,b,c,d, uc_mism_s, uc_mism_g, uc_mism_f, uc_mism_o, uc_mism_list)
}

uc_mydb = read.delim("./results/local/usearch/local_map300.uc", header = FALSE)
uc_mydb_seed = uc_mydb %>% select(V1, V9) %>% filter(V1 == "S")
unique(uc_mydb_seed$V1) # check what they are, there are only S not H and C
uc_mydb_hit = uc_mydb %>% select(V1, V9, V10) %>% filter(V1 == "H")
unique(uc_mydb_hit$V1) # check what they are, there are only H

# V1 = S/H/C
# v9 = label_taxa
# v10 = Label of target sequence (H records only)
# V9_2 = V9
# V10_2 = V10
uc_mydb_hit$V9_2 = sub(".*tax=", "", uc_mydb_hit$V9)
uc_mydb_hit$V10_2 = sub(".*tax=", "", uc_mydb_hit$V10)
taxa = c(",o:", ",f:", ",g:", ",s:")
uc_mydb_mism = mism(taxa, uc_mydb_hit)

uc_mydb_mism$mism_species
uc_mydb_mism$mism_genus
uc_mydb_mism$mism_family
uc_mydb_mism$mism_order

uc_mydb_mism$mism_df
temp = unique(uc_mydb_mism$mism_species$V10)
str(temp)
head(temp)

mydb_unique = readDNAStringSet("./results/local/usearch/db.local2.final.min300.ungap.unique.fa") # usearch fastx output
mydb_unique_df = data.frame(mydb_unique)
mydb_unique_df$taxa = mydb_unique@ranges@NAMES

head(mydb_unique_df$taxa)
mydb_unique_df$taxa = str_replace(mydb_unique_df$taxa, "size=.*", "")
mydb_unique_df$taxa = str_replace(mydb_unique_df$taxa, " ;", " ")
head(mydb_unique_df$taxa)
head(temp)

length(mydb_unique_df$taxa) ; length(temp)
mydb_unique_df = subset(mydb_unique_df, !(taxa %in% temp)) # dropping bad sequences
length(mydb_unique_df$taxa)
head(mydb_unique_df$taxa)
head(mydb_unique_df$mydb_unique)
length(mydb_unique) ; dim(mydb_unique_df)[1]

# reformat for blast ---
head(mydb_unique_df$taxa)
mydb_unique_df$taxa = str_replace_all(mydb_unique_df$taxa, "; ", " ") # replace semicolumn+space with space
mydb_unique_df$taxa = str_replace_all(mydb_unique_df$taxa, ";", " ") # replace semicolumn with space
mydb_unique_df$taxa = str_replace_all(mydb_unique_df$taxa, ",", " ") # replace comma with space
head(mydb_unique_df$taxa)
# reformat for blast ===

temp2 = DNAStringSet(mydb_unique_df$mydb_unique)
names(temp2) = mydb_unique_df$taxa
writeXStringSet(temp2, "./results/local/usearch/db.local2.final.min300.ungap.unique.clean.fa")
rm(temp, temp2)

uc_imdb = read.delim("./results/regional/usearch/regional_map300.uc", header = FALSE)
uc_imdb_seed = uc_imdb %>% select(V1, V9) %>% filter(V1 == "S")
unique(uc_imdb_seed$V1) # check what they are, there are only S not H and C
uc_imdb_hit = uc_imdb %>% select(V1, V9, V10) %>% filter(V1 == "H")
unique(uc_imdb_hit$V1) # check what they are, there are only H

uc_imdb_hit$V9_2 = sub(".*tax=", "", uc_imdb_hit$V9)
uc_imdb_hit$V10_2 = sub(".*tax=", "", uc_imdb_hit$V10)
taxa = c(",o:", ",f:", ",g:", ",s:")
uc_imdb_mism = mism(taxa, uc_imdb_hit)

uc_imdb_mism$mism_species
uc_imdb_mism$mism_genus
uc_imdb_mism$mism_family
uc_imdb_mism$mism_order

uc_imdb_mism$mism_df
temp = unique(uc_imdb_mism$mism_species$V10)
str(temp)
head(temp)

imdb_unique = readDNAStringSet("./results/regional/usearch/db.regional2.final.min300.ungap.unique.fa") # usearch fastx output
imdb_unique_df = data.frame(imdb_unique)
imdb_unique_df$taxa = imdb_unique@ranges@NAMES

head(imdb_unique_df$taxa)
imdb_unique_df$taxa = str_replace(imdb_unique_df$taxa, "size=.*", "")
imdb_unique_df$taxa = str_replace(imdb_unique_df$taxa, " ;", " ")
head(imdb_unique_df$taxa)
head(temp)

length(imdb_unique_df$taxa) ; length(temp)
imdb_unique_df = subset(imdb_unique_df, !(taxa %in% temp)) # dropping bad sequences
length(imdb_unique_df$taxa)
head(imdb_unique_df$taxa)
head(imdb_unique_df$imdb_unique)
length(imdb_unique) ; dim(imdb_unique_df)[1]

# reformat for blast
# ---
head(imdb_unique_df$taxa)
imdb_unique_df$taxa = str_replace_all(imdb_unique_df$taxa, "; ", " ") # replace semicolumn+space with space
imdb_unique_df$taxa = str_replace_all(imdb_unique_df$taxa, ";", " ") # replace semicolumn with space
imdb_unique_df$taxa = str_replace_all(imdb_unique_df$taxa, ",", " ") # replace comma with space
head(imdb_unique_df$taxa)
# ===

temp2 = DNAStringSet(imdb_unique_df$imdb_unique)
names(temp2) = imdb_unique_df$taxa
writeXStringSet(temp2, "./results/regional/usearch/db.regional2.final.min300.ungap.unique.clean.fa")
rm(temp, temp2)

uc_global = read.delim("./results/global/usearch/global_map300.uc", header = FALSE)
uc_global_seed = uc_imdb %>% select(V1, V9) %>% filter(V1 == "S")
unique(uc_global_seed$V1) # check what they are, there are only S not H and C
uc_global_hit = uc_global %>% select(V1, V9, V10) %>% filter(V1 == "H")
unique(uc_global_hit$V1) # check what they are, there are only H

uc_global_hit$V9_2 = sub(".*tax=", "", uc_global_hit$V9)
uc_global_hit$V10_2= sub(".*tax=", "", uc_global_hit$V10)
taxa = c(",o:", ",f:", ",g:", ",s:")
uc_global_mism = mism(taxa, uc_global_hit)

uc_global_mism$mism_species
uc_global_mism$mism_genus
uc_global_mism$mism_family
uc_global_mism$mism_order

uc_global_mism$mism_df
temp = unique(uc_global_mism$mism_species$V10)
str(temp)
head(temp)

global_unique = readDNAStringSet("./results/global/usearch/db.global2.final.min300.ungap.unique.fa") # usearch fastx output
global_unique_df = data.frame(global_unique)
global_unique_df$taxa = global_unique@ranges@NAMES

head(global_unique_df$taxa)
global_unique_df$taxa = str_replace(global_unique_df$taxa, "size=.*", "")
global_unique_df$taxa = str_replace(global_unique_df$taxa, " ;", " ")
head(global_unique_df$taxa)
head(temp)

length(global_unique_df$taxa) ; length(temp)
global_unique_df = subset(global_unique_df, !(taxa %in% temp)) # dropping bad sequences
length(global_unique_df$taxa)
head(global_unique_df$taxa)
head(global_unique_df$global_unique)

length(global_unique) ; dim(global_unique_df)[1]

# reformat for blast ---
head(global_unique_df$taxa)
global_unique_df$taxa = str_replace_all(global_unique_df$taxa, "; ", " ") # replace semicolumn+space with space
global_unique_df$taxa = str_replace_all(global_unique_df$taxa, ";", " ") # replace semicolumn with space
global_unique_df$taxa = str_replace_all(global_unique_df$taxa, ",", " ") # replace comma with space
head(global_unique_df$taxa)
# reformat for blast ===

temp2 = DNAStringSet(global_unique_df$global_unique)
names(temp2) = global_unique_df$taxa
writeXStringSet(temp2, "./results/global/usearch/db.global2.final.min300.ungap.unique.clean.fa")
rm(temp, temp2)

gc()
save.image("./mapuc_300.RData")

sessioninfo::session_info(pkgs="attached")
