### lama refseq shotgun plotting
library(dplyr)
library(tidyr)
library(ggplot2)
library(analogue)

setwd("C:/Users/ugcabuk/Desktop/test_barbara")

### kraken outputs
apmg <- read.delim("Lama_APMG_nt2022_0.8.txt", header = FALSE)
names(apmg) <- c("samples", "percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")

bhv <- read.delim("Lama_BHV_nt2022_0.8.txt", header = FALSE)
names(bhv) <- c("samples", "percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")

reseq <- read.delim("Lama_resequenced_nt2022_0.8.txt", header = FALSE)
names(reseq) <- c("samples", "percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")

apmg_lin <- read.delim2("lineage_lama_apmg3738_nt_0.8.csv", sep = ",", header = TRUE)
bhv_lin <- read.delim2("lineage_lama_bhv1_nt_0.8.csv", sep = ",", header = TRUE)
reseq_lin <- read.delim2("lineage_lama_bhv2_nt_0.8.csv", sep = ",", header = TRUE)


apmg.lineage <- plyr::join(apmg, apmg_lin, by = "taxID")
bhv.lineage <- plyr::join(bhv, bhv_lin, by = "taxID")
reseq.lineage <- plyr::join(reseq, reseq_lin, by = "taxID")

### load metadata containing sample name, age and depth
metadata_apmg <- openxlsx::read.xlsx(xlsxFile = "APMG3738_Lama_all_shotgun_sample_name_change.xlsx")
metadata_bhv <- openxlsx::read.xlsx(xlsxFile = "BHV1_shotgun_sample_name_change.xlsx")

metadata_apmg$samples %in% apmg$samples 
metadata_bhv$samples %in% bhv$samples

metadata_reseq <- openxlsx::read.xlsx(xlsxFile = "BHV2_shotgun_sample_name_change.xlsx")

### check if there are any missing samples
metadata_reseq$samples %in% reseq$samples 

shotgun.join_apmg <- plyr::join(apmg.lineage, metadata_apmg, by = "samples")
shotgun.join_apmg <- data.frame(shotgun.join_apmg)
shotgun.join_apmg$age <- as.numeric(shotgun.join_apmg$age)

shotgun.join_bhv <- plyr::join(bhv.lineage, metadata_bhv, by = "samples")
shotgun.join_bhv <- data.frame(shotgun.join_bhv)
shotgun.join_bhv$age <- as.numeric(shotgun.join_bhv$age)

shotgun.join_reseq <- plyr::join(reseq.lineage, metadata_reseq, by = "samples")
shotgun.join_reseq <- data.frame(shotgun.join_reseq)
shotgun.join_reseq$age <- as.numeric(shotgun.join_reseq$age)

shotgun_all <- rbind(shotgun.join_apmg, shotgun.join_bhv, shotgun.join_reseq)

### bacteria abundance as stratigraphic plot
### filter to bacteria superkingdom
bacteria <- filter(shotgun_all, type == "sample" & superkingdom == "Bacteria")
### extract at the genus and species rank
bact_g_s <- filter(bacteria, Rank == "G" | Rank == "S")

selected.bact_g_s <- dplyr::select(bact_g_s, c("age", "Name", "TaxCount"))

### remove rows with zero counts in TaxCount
selected.bact_g_s$Name <- stringr::str_replace_all(selected.bact_g_s$Name,"_", " ") 

bact_contam <- read.delim2("unique_bacteria_gen_spe_count20_assigned.csv", header = T, sep = ";")
bact_contam$Name <- stringr::str_replace_all(bact_contam$Name,"_", " ") 

bact_gen_spe_assign <- merge(selected.bact_g_s, bact_contam, by = "Name", all = TRUE)

### drop aquatic , NA and contaminant
bact_gen_spe_assign_clean <- bact_gen_spe_assign[!(bact_gen_spe_assign$element_cycle == "aquatic" | bact_gen_spe_assign$element_cycle == "contaminant" |
                                                     bact_gen_spe_assign$element_cycle == "NA" |bact_gen_spe_assign$element_cycle == ""),]

### make new Name_assigned with Name, assignment, real_5
bact_gen_spe_assign_clean$Name_assigned <- paste(bact_gen_spe_assign_clean$Name, bact_gen_spe_assign_clean$element_cycle, bact_gen_spe_assign_clean$real_5, sep ="_")
### select columns to convert
selected.columns_bac <- dplyr::select(bact_gen_spe_assign_clean, c("age", "Name_assigned", "TaxCount"))

### make the data from long to wide
wide_bact <- pivot_wider(selected.columns_bac, names_from = Name_assigned, names_sep = "_", values_from = TaxCount, values_fn = sum)
wide_bact[is.na(wide_bact)] <- 0
wide_bact <- wide_bact[-45,]
wide_bact = wide_bact[-which(wide_bact$age == 0),]

### table for resampling
resampling_bact <- as.data.frame(wide_bact[, -1])
resampling_bact_occ3 <- chooseTaxa(resampling_bact, n.occ = 3)
resampling_bact_occ3$age <- wide_bact$age
resampling_bact_occ3 <- resampling_bact_occ3 %>%
  dplyr::select(age, everything())

write.csv2(resampling_bact_occ3, "bacteria_assigned_resampling_occ3.csv")


### FUNGI
fungi <- filter(shotgun_all, type == "sample" & kingdom == "Fungi")

selected.fungi <- dplyr::select(fungi, c("age", "Name", "TaxCount"))

### remove rows with zero counts in TaxCount
selected.fungi$Name <- stringr::str_replace_all(selected.fungi$Name,"_", " ") 

### for contaminant selection
fungi_unique <- unique(fungi$Name)

fung_contam <- read.delim2("unique_fungi_all_assigned.csv", header = T, sep = ";")
fung_contam$Name <- stringr::str_replace_all(fung_contam$Name,"_", " ") 

fungi_assign <- merge(selected.fungi, fung_contam, by = "Name", all = TRUE)
fungi_assignments <- unique(fungi_assign$ecology)

### drop aquatic , NA and contaminant
fungi_assign_clean <- fungi_assign[!(fungi_assign$ecology == "aquatic" | fungi_assign$ecology == "contaminant" |
                                       fungi_assign$ecology == "NA" | fungi_assign$ecology == ""),]


### make new Name_assigned with Name, assignment, real_5
fungi_assign_clean$Name_assigned <- paste(fungi_assign_clean$Name, fungi_assign_clean$ecology, fungi_assign_clean$real_5, sep ="_")
selected.columns_fungi <- dplyr::select(fungi_assign_clean, c("age", "Name_assigned", "TaxCount")) # select columns to convert


### make the data from long to wide
wide_fungi <- pivot_wider(selected.columns_fungi, names_from = Name_assigned, names_sep = "_", values_from = TaxCount, values_fn = sum)
wide_fungi[is.na(wide_fungi)] <- 0
wide_fungi <- wide_fungi[-45,]
wide_fungi = wide_fungi[-which(wide_fungi$age == 0),]

### table for resampling
resampling_fung <- as.data.frame(wide_fungi[, -1])
resampling_fung3 <- chooseTaxa(resampling_fung, n.occ = 3)
resampling_fung3$age <- wide_fungi$age
resampling_fung3 <- resampling_fung3 %>%
  dplyr::select(age, everything())

write.csv2(resampling_fung3, "fungi_assigned_resampling_occ3.csv")

### PLANT
vegetation <- filter(shotgun_all, type == "sample" & kingdom == "Viridiplantae")

selected.vege <- dplyr::select(vegetation, c("age", "Name", "TaxCount"))

### remove rows with zero counts in TaxCount
selected.vege$Name <- stringr::str_replace_all(selected.vege$Name,"_", " ") 

### for contaminant selection
vege_unique <- unique(vegetation$Name)

vege_contam <- read.delim2("unique_plant_gen_spe_assigned.csv", header = T, sep = ";")
vege_contam$Name <- stringr::str_replace_all(vege_contam$Name,"_", " ") 

vege_assign <- merge(selected.vege, vege_contam, by = "Name", all = TRUE)
vege_assignments <- unique(vege_assign$assignment)

### drop aquatic , NA and contaminant
vege_assign_clean <- vege_assign[!(vege_assign$assignment == "aquatic" | vege_assign$assignment == "contaminant" |
                                     vege_assign$assignment == "NA" | vege_assign$assignment == ""),]

### make new Name_assigned with Name, assignment, real_5
vege_assign_clean$Name_assigned <- paste(vege_assign_clean$Name, vege_assign_clean$assignment, vege_assign_clean$real_5, sep ="_")
selected.columns_veg <- dplyr::select(vege_assign_clean, c("age", "Name_assigned", "TaxCount")) # select columns to convert

#### make the data from long to wide
wide_veg <- pivot_wider(selected.columns_veg, names_from = Name_assigned, names_sep = "_", values_from = TaxCount, values_fn = sum)
wide_veg[is.na(wide_veg)] <- 0
wide_veg <- wide_veg[-45,]
wide_veg = wide_veg[-which(wide_veg$age == 0),]

### table for resampling
resampling_veg <- as.data.frame(wide_veg[, -1])
resampling_veg3 <- chooseTaxa(resampling_veg, n.occ = 3)
resampling_veg3$age <- wide_veg$age
resampling_veg3 <- resampling_veg3 %>%
  dplyr::select(age, everything())

write.csv2(resampling_veg3, "plant_assigned_resampling_occ3.csv")
