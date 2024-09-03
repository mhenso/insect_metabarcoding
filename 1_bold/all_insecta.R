## .tocify .tocify-header {
##   #position: fixed;
##   #top: 50px;
##   #left: 50px;
##   width: 350px;
##   #height: 400px;
##   }

suppressPackageStartupMessages({
  library(bold)
  library(taxize)
  library(dplyr)
})

allArthropod_names <- downstream("Arthropoda", db = "bold", downto = "class")
otherArth_names <- allArthropod_names$Arthropoda %>% filter(name != "Insecta") %>% select(name)
otherArth_list <- lapply(otherArth_names, bold_seqspec)

gc()
excludeNames <- c("Coleoptera","Diptera","Hymenoptera","Lepidoptera")
allInsect_names <- downstream("Insecta", db = "bold", downto = "order")
otherInsects_names <- allInsect_names$Insecta %>% filter(!name %in% excludeNames) %>% select(name)
otherInsects_list <- lapply(otherInsects_names, bold_seqspec)

# split into: Carabidae,Chrysomelidae,Curculionidae,Staphylinidae,(remaining others)
gc()
save.image()
Col_list <- downstream("Coleoptera", db = "bold", downto = "family")

# Carabidae
Col_Carabidae_names <- Col_list$Coleoptera %>% filter(name=="Carabidae") %>% select(name)
Col_Carabidae_list <- lapply(Col_Carabidae_names, bold_seqspec)

# Chrysomelidae
gc()
save.image()
Col_Chrysomelidae_names <- Col_list$Coleoptera %>% filter(name=="Chrysomelidae") %>% select(name)
Col_Chrysomelidae_list <- lapply(Col_Chrysomelidae_names, bold_seqspec)
# Col_Chrysomelidae_df <- gatherBOLDdat_function(Col_Chrysomelidae_list)

# Curculionidae
gc()
save.image()
Col_Curculionidae_names <- Col_list$Coleoptera %>% filter(name=="Curculionidae") %>% select(name)
Col_Curculionidae_list <- lapply(Col_Curculionidae_names, bold_seqspec)
# Col_Curculionidae_df <- gatherBOLDdat_function(Col_Curculionidae_list)

# Staphylinidae
gc()
save.image()
Col_Staphylinidae_names <- Col_list$Coleoptera %>% filter(name=="Staphylinidae") %>% select(name)
Col_Staphylinidae_list <- lapply(Col_Staphylinidae_names, bold_seqspec)
# Col_Staphylinidae_df <- gatherBOLDdat_function(Col_Staphylinidae_list)

# All others Coleoptera
gc()
save.image()
excludeColNames <- c("Carabidae","Chrysomelidae","Curculionidae","Staphylinidae")
Col_allother_names <- Col_list$Coleoptera %>% filter(!name %in% excludeColNames) %>% select(name)
Col_allothers_list <- lapply(Col_allother_names, bold_seqspec)

#### split into: Sciaridae,Cecidomyiidae,Chironomidae,(remaining others)
gc()
save.image()
diptera_list <- downstream("Diptera", db = "bold", downto = "family")

# Sciaridae
Dip_Sciaridae_names <- diptera_list$Diptera %>% filter(name=="Sciaridae") %>% select(name)
Dip_Sciaridae_list <- lapply(Dip_Sciaridae_names, bold_seqspec)

# Cecidomyiidae
Dip_Cecidomyiidae_names <- diptera_list$Diptera %>% filter(name=="Cecidomyiidae") %>% select(name)
Dip_Cecidomyiidae_list <- lapply(Dip_Cecidomyiidae_names, bold_seqspec)

# Chironomidae
Dip_Chironomidae_names <- diptera_list$Diptera %>% filter(name=="Chironomidae") %>% select(name)
Dip_Chironomidae_list <- lapply(Dip_Chironomidae_names, bold_seqspec)

# Phoridae
Dip_Phoridae_names <- diptera_list$Diptera %>% filter(name=="Phoridae") %>% select(name)
Dip_Phoridae_list <- lapply(Dip_Phoridae_names, bold_seqspec)

# Ceratopogonidae
Dip_Ceratopogonidae_names <- diptera_list$Diptera %>% filter(name=="Ceratopogonidae") %>% select(name)
Dip_Ceratopogonidae_list <- lapply(Dip_Ceratopogonidae_names, bold_seqspec)

# All others Diptera
excludeDipNames <- c("Sciaridae", "Cecidomyiidae", "Chironomidae", "Phoridae", "Ceratopogonidae")
Dip_allother_names <- diptera_list$Diptera %>% filter(!name %in% excludeDipNames) %>% select(name)
Dip_allothers_list <- lapply(Dip_allother_names, bold_seqspec)

gc()
save.image()

#### split into: Braconidae, Formicidae, Ichneumonidae, Platygastridae,(remaining others)
gc()

Hym_list <- downstream("Hymenoptera", db = "bold", downto = "family")
# Braconidae
Hym_Braconidae_names <- Hym_list$Hymenoptera %>% filter(name=="Braconidae") %>% select(name)
Hym_Braconidae_list <- lapply(Hym_Braconidae_names, bold_seqspec)

# Formicidae
gc()
save.image()
Hym_Formicidae_names <- Hym_list$Hymenoptera %>% filter(name=="Formicidae") %>% select(name)
Hym_Formicidae_list <- lapply(Hym_Formicidae_names, bold_seqspec)

# Ichneumonidae
gc()
save.image()
Hym_Ichneumonidae_names <- Hym_list$Hymenoptera %>% filter(name=="Ichneumonidae") %>% select(name)
Hym_Ichneumonidae_list <- lapply(Hym_Ichneumonidae_names, bold_seqspec)

# Platygastridae
gc()
save.image()
Hym_Platygastridae_names <- Hym_list$Hymenoptera %>% filter(name=="Platygastridae") %>% select(name)
Hym_Platygastridae_list <- lapply(Hym_Platygastridae_names, bold_seqspec)

# All others Hymenoptera
excludeHymNames <- c("Braconidae","Formicidae","Ichneumonidae","Platygastridae")
Hym_allother_names <- Hym_list$Hymenoptera %>% filter(!name %in% excludeHymNames) %>% select(name)
Hym_allothers_list <- lapply(Hym_allother_names, bold_seqspec)

#### split into: Noctuidae,Erebidae,Sphingidae,Geometridae,(remaining others)
gc()
save.image()
Lep_list <- downstream("Lepidoptera", db = "bold", downto = "family")
# Noctuidae
Lep_Noctuidae_names <- Lep_list$Lepidoptera %>% filter(name=="Noctuidae") %>% select(name)
Lep_Noctuidae_list <- lapply(Lep_Noctuidae_names, bold_seqspec)

# Erebidae
gc()
save.image()
Lep_Erebidae_names <- Lep_list$Lepidoptera %>% filter(name=="Erebidae") %>% select(name)
Lep_Erebidae_list <- lapply(Lep_Erebidae_names, bold_seqspec)

# Sphingidae
gc()
save.image()
Lep_Sphingidae_names <- Lep_list$Lepidoptera %>% filter(name=="Sphingidae") %>% select(name)
Lep_Sphingidae_list <- lapply(Lep_Sphingidae_names, bold_seqspec)

# Geometridae
gc()
save.image()
Lep_Geometridae_names <- Lep_list$Lepidoptera %>% filter(name=="Geometridae") %>% select(name)
Lep_Geometridae_list <- lapply(Lep_Geometridae_names, bold_seqspec)

# Crambidae
gc()
Lep_Crambidae_names <- Lep_list$Lepidoptera %>% filter(name=="Crambidae") %>% select(name)
Lep_Crambidae_list <- lapply(Lep_Crambidae_names, bold_seqspec)

# Gelechiidae
gc()
Lep_Gelechiidae_names <- Lep_list$Lepidoptera %>% filter(name=="Gelechiidae") %>% select(name)
Lep_Gelechiidae_list <- lapply(Lep_Gelechiidae_names, bold_seqspec)

# Saturniidae
gc()
Lep_Saturniidae_names <- Lep_list$Lepidoptera %>% filter(name=="Saturniidae") %>% select(name)
Lep_Saturniidae_list <- lapply(Lep_Saturniidae_names, bold_seqspec)

# All other Lepidoptera
gc()
save.image()
excludeLepNames <- c("Noctuidae","Erebidae","Sphingidae","Geometridae","Crambidae", "Gelechiidae", "Saturniidae")
Lep_allother_names <- Lep_list$Lepidoptera %>% filter(!name %in% excludeLepNames) %>% select(name)
Lep_allothers_list <- lapply(Lep_allother_names, bold_seqspec)
