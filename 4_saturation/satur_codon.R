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
  library(ape)
})

local <- read.dna(file="./data/db.local2.satur.afa", format="fasta")
local.d_TN93 = data.frame(d_TN93 = c(dist.dna(local, model="TN93", as.matrix = TRUE)
                                     [upper.tri(dist.dna(local, model="TN93", as.matrix = TRUE))]))
length(which(local.d_TN93$d_TN93 > 0))

local_1 <- read.dna(file="./data/dambe/db.local2.satur_1st.FAS", format="fasta")
dim(local_1)
class(local_1)

local_1.TS = data.frame(d1.TS = c(dist.dna(local_1, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(local_1, model="TS", as.matrix = TRUE))]))

local_1.TV = data.frame(d1.TV = c(dist.dna(local_1, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(local_1, model="TV", as.matrix = TRUE))]))

# 3rd codon
local_3 <- read.dna(file="./data/dambe/db.local2.satur_3rd.FAS", format="fasta")
local_3.TS = data.frame(d3.TS = c(dist.dna(local_3, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(local_3, model="TS", as.matrix = TRUE))]))
local_3.TV = data.frame(d3.TV = c(dist.dna(local_3, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(local_3, model="TV", as.matrix = TRUE))]))

local.d_merg = data.frame(d_TN93=local.d_TN93$d_TN93, 
                          d1_TS=local_1.TS$d1.TS/106, 
                          d1_TV=local_1.TV$d1.TV/106,
                          d3_TS=local_3.TS$d3.TS/107, 
                          d3_TV=local_3.TV$d3.TV/107)
dim(local.d_merg)
local.d_merg %>% drop_na(d_TN93) %>% dim()
local.d_merg %>% drop_na(d1_TS) %>% dim()
local.d_merg %>% drop_na(d1_TV) %>% dim()
local.d_merg %>% drop_na(d3_TS) %>% dim()
local.d_merg %>% drop_na(d3_TV) %>% dim()

regional <- read.dna(file="./data/db.regional2.satur.afa", format="fasta")
regional.d_TN93 = data.frame(d_TN93 = c(dist.dna(regional, model="TN93", as.matrix = TRUE)
                                        [upper.tri(dist.dna(regional, model="TN93", as.matrix = TRUE))]))

library(ape)
regional_1 <- read.dna(file="./data/dambe/db.regional2.satur_1st.FAS", format="fasta")

regional_1.TS = data.frame(d1.TS = c(dist.dna(regional_1, model="TS", as.matrix = TRUE)
                                   [upper.tri(dist.dna(regional_1, model="TS", as.matrix = TRUE))]))

regional_1.TV = data.frame(d1.TV = c(dist.dna(regional_1, model="TV", as.matrix = TRUE)
                                   [upper.tri(dist.dna(regional_1, model="TV", as.matrix = TRUE))]))

regional_3 <- read.dna(file="./data/dambe/db.regional2.satur_3rd.FAS", format="fasta")
regional_3.TS = data.frame(d3.TS = c(dist.dna(regional_3, model="TS", as.matrix = TRUE)
                                   [upper.tri(dist.dna(regional_3, model="TS", as.matrix = TRUE))]))
regional_3.TV = data.frame(d3.TV = c(dist.dna(regional_3, model="TV", as.matrix = TRUE)
                                   [upper.tri(dist.dna(regional_3, model="TV", as.matrix = TRUE))]))

regional.d_merg = data.frame(d_TN93=regional.d_TN93$d_TN93, 
                          d1_TS=regional_1.TS$d1.TS/106, 
                          d1_TV=regional_1.TV$d1.TV/106,
                          d3_TS=regional_3.TS$d3.TS/107, 
                          d3_TV=regional_3.TV$d3.TV/107)
dim(regional.d_merg)
regional.d_merg %>% drop_na(d_TN93) %>% dim()
regional.d_merg %>% drop_na(d1_TS) %>% dim() 
regional.d_merg %>% drop_na(d1_TV) %>% dim()
regional.d_merg %>% drop_na(d3_TS) %>% dim()
regional.d_merg %>% drop_na(d3_TV) %>% dim()

regional.d_merg = regional.d_merg %>% drop_na(d_TN93)
dim(regional.d_merg)

rm(list=grep("local.d_merg|regional.d_merg", ls(), value = TRUE, invert = TRUE))
gc()

summary(lm(d1_TS ~ I(d_TN93^1), data=local.d_merg))

summary(lm(d1_TV ~ I(d_TN93^1), data=local.d_merg))

summary(lm(d3_TS ~ I(d_TN93^1), data=local.d_merg))

summary(lm(d3_TV ~ I(d_TN93^1), data=local.d_merg))

summary(lm(d1_TS ~ I(d_TN93^1), data=regional.d_merg))

summary(lm(d1_TV ~ I(d_TN93^1), data=regional.d_merg))

summary(lm(d3_TS ~ I(d_TN93^1), data=regional.d_merg))

summary(lm(d3_TV ~ I(d_TN93^1), data=regional.d_merg))

rm(list=ls())
gc()

global_101 <- read.dna(file="./data/global/global_101.afa", format="fasta")
global_101.d_TN93 = data.frame(d101_TN93 = c(dist.dna(global_101, model="TN93", as.matrix = TRUE)
                                             [upper.tri(dist.dna(global_101, model="TN93", as.matrix = TRUE))]))

global_102 <- read.dna(file="./data/global/global_102.afa", format="fasta")
global_102.d_TN93 = data.frame(d102_TN93 = c(dist.dna(global_102, model="TN93", as.matrix = TRUE)
                                          [upper.tri(dist.dna(global_102, model="TN93", as.matrix = TRUE))]))

global_103 <- read.dna(file="./data/global/global_103.afa", format="fasta")
global_103.d_TN93 = data.frame(d103_TN93 = c(dist.dna(global_103, model="TN93", as.matrix = TRUE)
                                          [upper.tri(dist.dna(global_103, model="TN93", as.matrix = TRUE))]))

# All sequences of same length: 106 
global_101_1 <- read.dna(file="./data/global/dambe/global_101_1st.afa", format="fasta")

global_101_1.TS = data.frame(d101.1.TS = c(dist.dna(global_101_1, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_101_1, model="TS", as.matrix = TRUE))]))

global_101_1.TV = data.frame(d101.1.TV = c(dist.dna(global_101_1, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_101_1, model="TV", as.matrix = TRUE))]))

global_102_1 <- read.dna(file="./data/global/dambe/global_102_1st.afa", format="fasta")

global_102_1.TS = data.frame(d102.1.TS = c(dist.dna(global_102_1, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_102_1, model="TS", as.matrix = TRUE))]))

global_102_1.TV = data.frame(d102.1.TV = c(dist.dna(global_102_1, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_102_1, model="TV", as.matrix = TRUE))]))

global_103_1 <- read.dna(file="./data/global/dambe/global_103_1st.afa", format="fasta")

global_103_1.TS = data.frame(d103.1.TS = c(dist.dna(global_103_1, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_103_1, model="TS", as.matrix = TRUE))]))

global_103_1.TV = data.frame(d103.1.TV = c(dist.dna(global_103_1, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_103_1, model="TV", as.matrix = TRUE))]))

global_101_3 <- read.dna(file="./data/global/dambe/global_101_3rd.afa", format="fasta")

global_101_3.TS = data.frame(d101.3.TS = c(dist.dna(global_101_3, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_101_3, model="TS", as.matrix = TRUE))]))

global_101_3.TV = data.frame(d101.3.TV = c(dist.dna(global_101_3, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_101_3, model="TV", as.matrix = TRUE))]))

global_102_3 <- read.dna(file="./data/global/dambe/global_102_3rd.afa", format="fasta")

global_102_3.TS = data.frame(d102.3.TS = c(dist.dna(global_102_3, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_102_3, model="TS", as.matrix = TRUE))]))

global_102_3.TV = data.frame(d102.3.TV = c(dist.dna(global_102_3, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_102_3, model="TV", as.matrix = TRUE))]))

global_103_3 <- read.dna(file="./data/global/dambe/global_103_3rd.afa", format="fasta")

global_103_3.TS = data.frame(d103.3.TS = c(dist.dna(global_103_3, model="TS", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_103_3, model="TS", as.matrix = TRUE))]))

global_103_3.TV = data.frame(d103.3.TV = c(dist.dna(global_103_3, model="TV", as.matrix = TRUE)
                                [upper.tri(dist.dna(global_103_3, model="TV", as.matrix = TRUE))]))

global_101.d_merg = data.frame(d101_TN93=global_101.d_TN93$d101_TN93, # all seq
                               d101.1_TS=global_101_1.TS$d101.1.TS/106, # 1st codon
                               d101.1_TV=global_101_1.TV$d101.1.TV/106, 
                               d101.3_TS=global_101_3.TS$d101.3.TS/106, # 3rd codon
                               d101.3_TV=global_101_3.TV$d101.3.TV/106)
dim(global_101.d_merg)[1]

global_101.d_merg %>% drop_na(d101_TN93) %>% dim() 
global_101.d_merg %>% drop_na(d101.1_TS) %>% dim() 
global_101.d_merg %>% drop_na(d101.1_TV) %>% dim()
global_101.d_merg %>% drop_na(d101.3_TS) %>% dim() 
global_101.d_merg %>% drop_na(d101.3_TV) %>% dim()

global_101.d_merg = global_101.d_merg %>% drop_na(d101_TN93)
dim(global_101.d_merg)[1]
rm(global_101.d_TN93,
   global_101_1.TS, 
   global_101_1.TV,
   global_101_3.TS,
   global_101_3.TV)
gc()

global_102.d_merg = data.frame(d102_TN93=global_102.d_TN93$d102_TN93, # all seq
                               d102.1_TS=global_102_1.TS$d102.1.TS/106, # 1st codon
                               d102.1_TV=global_102_1.TV$d102.1.TV/106, 
                               d102.3_TS=global_102_3.TS$d102.3.TS/106, # 3rd codon
                               d102.3_TV=global_102_3.TV$d102.3.TV/106)
dim(global_102.d_merg)[1]
head(global_102.d_merg)

global_102.d_merg %>% drop_na(d102_TN93) %>% dim() 
global_102.d_merg %>% drop_na(d102.1_TS) %>% dim() 
global_102.d_merg %>% drop_na(d102.1_TV) %>% dim()
global_102.d_merg %>% drop_na(d102.3_TS) %>% dim() 
global_102.d_merg %>% drop_na(d102.3_TV) %>% dim()

global_102.d_merg = global_102.d_merg %>% drop_na(d102_TN93)

dim(global_102.d_merg)[1]

rm(global_102.d_TN93,
   global_102_1.TS, 
   global_102_1.TV,
   global_102_3.TS,
   global_102_3.TV)
gc()

global_103.d_merg = data.frame(d103_TN93=global_103.d_TN93$d103_TN93, # all seq
                               d103.1_TS=global_103_1.TS$d103.1.TS/106, # 1st codon
                               d103.1_TV=global_103_1.TV$d103.1.TV/106, 
                               d103.3_TS=global_103_3.TS$d103.3.TS/106, # 3rd codon
                               d103.3_TV=global_103_3.TV$d103.3.TV/106)
dim(global_103.d_merg)[1]

global_103.d_merg %>% drop_na(d103_TN93) %>% dim() 
global_103.d_merg %>% drop_na(d103.1_TS) %>% dim() 
global_103.d_merg %>% drop_na(d103.1_TV) %>% dim()
global_103.d_merg %>% drop_na(d103.3_TS) %>% dim() 
global_103.d_merg %>% drop_na(d103.3_TV) %>% dim()

global_103.d_merg = global_103.d_merg %>% drop_na(d103_TN93)

dim(global_103.d_merg)[1]

rm(global_103.d_TN93,
   global_103_1.TS, 
   global_103_1.TV,
   global_103_3.TS,
   global_103_3.TV)
gc()

names(global_101.d_merg)
names(global_102.d_merg)
names(global_103.d_merg)
# renames column so we can rbind
names(global_102.d_merg) = names(global_101.d_merg)
names(global_103.d_merg) = names(global_101.d_merg)

names(global_101.d_merg)
names(global_102.d_merg)
names(global_103.d_merg)

global.d = rbind(global_101.d_merg, global_102.d_merg, global_103.d_merg)
names(global.d)
names(global.d) = c("TN93", "TS1", "TV1", "TS3", "TV3")
head(global.d)

rm(list=grep("global.d", ls(), value = TRUE, invert = TRUE))
gc()

summary(lm(TS1 ~ I(TN93^1), data=global.d))
gc()

summary(lm(TV1 ~ I(TN93^1), data=global.d))
gc()

summary(lm(TS3 ~ I(TN93^1), data=global.d))
gc()

summary(lm(TV3 ~ I(TN93^1), data=global.d))
gc()

sessioninfo::session_info(pkgs="attached")
