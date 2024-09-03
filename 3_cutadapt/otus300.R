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
  library(Biostrings)
  library(reshape2)
  library(phyloseq)
  library(iNEXT)
  library(ggpubr)
  library(DT)
})

otu = read.table("./output/otus_finalize/otutab.300.txt", row.names = 1)
dim(otu)

## grep "^#" ./output/otus_finalize/otutab.300.txt > temp_colnames
## sed -E 's/#OTU ID\t//' temp_colnames > temp_colnames2
## ### Using Shell Terminal ===

temp_colnames = read.table("temp_colnames2")
temp_colnames = gsub("COI-", "", temp_colnames)
temp_colnames = gsub("_S.*", "", temp_colnames)
temp_colnames = data.frame(t(temp_colnames))
colnames(otu)= temp_colnames[1,]

temp = c("TK05", "TK06", "TK07", "TK08", "TK09", "TK10")
meta = data.frame("TK_number"=colnames(otu))
meta$lt_code = meta$TK_number
meta$lt_code = str_sub(meta$lt_code, 1, 4)

# check ligth trap names
unique(meta$lt_code)

rownames(meta) = meta$TK_number
meta %>% datatable()

ref_seq = Biostrings::readDNAStringSet("./output/otus_finalize/otus97.final.300.fa")
ref_seq
range(ref_seq@ranges@width)
ref_seq2 = data.frame(ref_seq)
ref_seq2$length = ref_seq@ranges@width

otu2 = otu_table(as.matrix(otu), taxa_are_rows = T)

dim(meta)
dim(otu2)

length(ref_seq@ranges@group)

otu_phy = merge_phyloseq(otu2, ref_seq, sample_data(meta))
otu_phy

# check if otu is correct
length(intersect(rownames(otu_table(otu_phy)), rownames(otu)))

otu_info = DataInfo(otu)
otu_info = otu_info %>% arrange(site)

head(otu_info)
otu_info$site2 = str_sub(otu_info$site, 1,4)

library(doParallel)
ncores <- 6
cl <- makeCluster(ncores)
registerDoParallel(cl)

est <- function(TK_sample){
  require(iNEXT)
  require(tidyverse)
  require(phyloseq)
  a = data.frame(otu_table(otu_phy))
  b = data.frame(sample_data(otu_phy)) %>% subset(lt_code == TK_sample) %>% select(TK_number)
  b = sort(b$TK_number)
  a = a %>% select(all_of(b))
  
  c = estimateD(a, datatype="abundance", base="coverage", level=NULL, conf=0.95) 
  return(c)
}

temp = unique(otu_info$site2)

system.time(est_otu_all <- foreach(i=temp) %dopar% est(i)) # another way
stopCluster(cl)

est_otu_all2 = as.data.frame(do.call(rbind, est_otu_all))
est_otu_all2 = est_otu_all2 %>% filter(order == 0)

temp = est_otu_all2$m # get the number of abundance of each subsample that has same coverage across each sample (rdept)
names(temp) = est_otu_all2$site

tmp_otu = otu[,-c(1:72)] # create empty otu_table as the original otu_table

for (i in 1:length(temp)){
  a = names(temp[i])
  rdepth = temp[[i]]
  tmp = otu %>% select(all_of(a))
  tmp = SRS::SRS(tmp, Cmin = rdepth, set_seed = T, seed = 42) # SRS::SRS
  tmp_otu = cbind(tmp_otu, tmp)
}
head(tmp_otu)

# comparing the randomized subsample otu table (tmp_otu) with the intended sample number(temp)
# colSums(otu) # original otu_table
# colSums(tmp_otu) # subsample otu_table
identical(colSums(tmp_otu), temp)

# make temporary phyloseq otu_phy_srs
otu_phy_srs = otu_phy

# replacing the original otu_table with the srs one
tmp_otu2 = otu_table(as.matrix(tmp_otu), taxa_are_rows = T)

otu_table(otu_phy_srs) = tmp_otu2
otu_phy_srs

# colSums(otu_table(otu_phy_srs))
# colSums(tmp_otu)
# check if it is identical
identical(colSums(otu_table(otu_phy_srs)), colSums(tmp_otu))

temp = unique(sample_data(otu_phy_srs)$lt_code)

otu_srs = data.frame(otu_table(otu_phy_srs))
otu_srs_list = list()

for (i in temp){
  a = subset(sample_data(otu_phy_srs), subset = lt_code %in% i)
  b = a$TK_number
  c = otu_srs[,b]
  c[c>0]=1
  otu_srs_list[[i]] = c
}
rm(a, b, c, i, temp)

library(iNEXT)
otu_srs_list_out = iNEXT(otu_srs_list, q=c(0,1,2), datatype="incidence_raw", nboot = 100)

po1 = ggiNEXT(otu_srs_list_out, facet.var = "order", type=1) + theme_bw(base_size=25) +
  ylab("OTU") + xlab("Number of subsamples") +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) +
  scale_y_continuous(breaks = seq(300, 1300, by = 100)) +
  theme(legend.position="none", axis.text.x = element_text(size=rel(0.6), angle=35, hjust=0.5, vjust=0.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA))


po2 = ggiNEXT(otu_srs_list_out, type=2) + theme_bw(base_size=25) + xlab("Number of subsamples") +
  scale_x_continuous(breaks = seq(0, 24, by = 1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.05)) +
  theme(legend.position = c(0.8,0.42),
        legend.text = element_text(size =25),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=rel(0.6), angle=35, hjust=0.5, vjust=0.8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank()) +
  geom_vline(xintercept = 4, color = "blue", linetype= "dashed", size=1.2)
ggarrange(po1, po2, nrow = 1, widths = c(0.65, 0.35))

est_otu_srs_size12 = estimateD(otu_srs_list, datatype="incidence_raw", base="size", level=12, conf=0.95)
est_otu_srs_size12 %>% filter(order == 0)
est_otu_srs_size12 %>% filter(order == 0) %>% select(SC) %>% summarise(mean= mean(SC), sd = sd(SC))

est_otu_srs_cover90 = estimateD(otu_srs_list, datatype="incidence_raw", base="coverage", level=0.90, conf=0.95)
est_otu_srs_cover90 %>% filter(order == 0)
est_otu_srs_cover90 %>% filter(order == 0) %>% select(SC) %>% summarise(mean= round(mean(SC),2), sd = round(sd(SC),2))

est_otu_srs_cover95 = estimateD(otu_srs_list, datatype="incidence_raw", base="coverage", level=0.95, conf=0.95)
est_otu_srs_cover95 %>% filter(order == 0)
est_otu_srs_cover95 %>% filter(order == 0) %>% select(SC) %>% summarise(mean= round(mean(SC),3), sd = round(sd(SC),2))

est_otu_srs_size4 = estimateD(otu_srs_list, datatype="incidence_raw", base="size", level=4, conf=0.95)
est_otu_srs_size4 %>% filter(order == 0) 
est_otu_srs_size4 %>% filter(order == 0) %>% select(SC) %>% summarise(mean= mean(SC), sd = sd(SC))
est_otu_srs_size4 %>% filter(order == 0) %>% reframe(SC= round(SC,2))
est_otu_srs_size4 %>% filter(order == 0) %>% select(SC) %>% summarise(mean= round(mean(SC),2), sd = round(sd(SC),2))

otu_phy_srs
temp = unique(sample_data(otu_phy_srs)$lt_code)

# set.seed(12) #original seed number
sample4 = function(temp, numberofsubsample){
  set.seed(12)
  hasil = c()
  
  for (x in temp){
    hasil_TK = c()
    a = sample(LETTERS[1:numberofsubsample], 4)
    
    for (y in a){
      b = paste0(x, y)
      c = data.frame("TK" = x, "tk_sub"= b)
      hasil_TK = rbind(hasil_TK, c)
    }
    
    d = hasil_TK
    hasil = rbind(hasil, d)
  }
  
  return(hasil)
}
chosen4 = sample4(temp, 12)

chosen4
temp = unique(sample_data(otu_phy_srs)$lt_code)
e = otu_srs[,-c(1:72)] # otu level

for (i in temp){
  a = subset(chosen4, subset = TK %in% i)
  b = a$tk_sub

  c = otu_srs[,b] # night_otu
  c$tk = rowSums(c)
  d = data.frame(c[,5])
  rownames(d) = rownames(c)
  names(d) = i
  d[d>0]=1
  e = cbind(d,e)
  # prun_srs_night_otu_list = data.frame(d)
  # assign(paste0("night", i) , d)
}

otu_srs_night_list = list(e)
names(otu_srs_night_list) = "TKs"
rm(a, b, c, d, e, i)

otu_srs_night_list_out = iNEXT(otu_srs_night_list, q=c(0,1,2), datatype="incidence_raw", nboot = 100)

est_otu_srs_night_size = estimateD(otu_srs_night_list, datatype="incidence_raw", base="size", conf=0.95) # , level=6
est_otu_srs_night_size

est_otu_srs_night_cover95 = estimateD(otu_srs_night_list, datatype="incidence_raw", base="coverage", level=0.95, conf=0.95)
est_otu_srs_night_cover95

est_otu_srs_night_cover90 = estimateD(otu_srs_night_list, datatype="incidence_raw", base="coverage", level=0.90, conf=0.95)
est_otu_srs_night_cover90

po1 = ggiNEXT(otu_srs_night_list_out, facet.var = "order", type=1) + 
  theme_bw(base_size=25) + 
  ylab("OTU") + xlab("Number of Nights") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(breaks = seq(1000, 5000, by = 1000)) +
  theme(legend.position="none", 
        # axis.text.x = element_text(size=rel(0.6), angle=35, hjust=0.5, vjust=0.8),
        axis.text.x = element_text(size=rel(1.2), angle=35, hjust=0.5, vjust=0.8),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour = "black", fill=NA))

po2 = ggiNEXT(otu_srs_night_list_out, type=2) + 
  theme_bw(base_size=25) +  xlab("Number of Nights") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  geom_vline(xintercept = 6, color = "blue", linetype= "dashed", size=1.2) +
  theme(legend.position = c(0.8,0.42),
        # axis.text.x = element_text(size=rel(0.6), angle=35, hjust=0.5, vjust=0.8),
        axis.text.x = element_text(size=rel(1.2), angle=35, hjust=0.5, vjust=0.8),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())

ggarrange(po1, po2, widths = c(0.65, 0.35))

otu_srs

kombina2 = function(samp_id, num_subsam, combi, tab_otu){
  a=c()
  for (z in LETTERS[1:num_subsam]){
    b = paste0(samp_id, z)
    a = c(a, b)
  }
  komb = data.frame(combn(a, combi))
  temp = c()
  
  for (y in 1:495){
    c = tab_otu %>% select(komb[,y])
    c[c>0] = 1
    d = data.frame("jum"=colSums(c))
    d$sub = rownames(d)
    e = data.frame("sam" = samp_id, "sam_comb"= paste0(samp_id, "_", y), "samsub"=d$sub, "jum"=d$jum)
    temp = rbind(temp,e)
    }
  return(temp)
}

jum.tk05 = kombina2("TK05", 12, 4, otu_srs)
jum.tk06 = kombina2("TK06", 12, 4, otu_srs)
jum.tk07 = kombina2("TK07", 12, 4, otu_srs)
jum.tk08 = kombina2("TK08", 12, 4, otu_srs)
jum.tk09 = kombina2("TK09", 12, 4, otu_srs)
jum.tk10 = kombina2("TK10", 12, 4, otu_srs)

jum.all = rbind(jum.tk05, jum.tk06, jum.tk07, jum.tk08, jum.tk09, jum.tk10)
rm(list=ls(pattern = "jum.tk"))

system.time(mod <- aov(jum ~ factor(sam)/factor(sam_comb), data = jum.all))

summary(mod)

system.time(tukey.mod <- TukeyHSD(mod))
# tukey.mod

sessioninfo::session_info(pkgs="attached")
