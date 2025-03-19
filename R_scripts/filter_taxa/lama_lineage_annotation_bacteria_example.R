###lama refseq shotgun plotting
setwd("...")
library(dplyr)
library(tidyr)
library(ggplot2)

apmg <- read.delim("Lama_APMG_nt2022_0.8.txt", header = FALSE)
names(apmg) <- c("samples", "percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")

bhv <- read.delim("Lama_BHV_nt2022_0.8.txt", header = FALSE)
names(bhv) <- c("samples", "percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")

reseq <- read.delim("Lama_resequenced_nt2022_0.8.txt", header = FALSE)
names(reseq) <- c("samples", "percentage", "CladeCount", "TaxCount", "Rank", "taxID", "Name")

apmg_lin <- read.delim2("lineage_lama_apmg3738_nt_0.8.csv", sep = ",", header = TRUE)
bhv_lin <- read.delim2("lineage_lama_bhv1_nt_0.8.csv", sep = ",", header = TRUE)
reseq_lin <- read.delim2("lineage_lama_bhv2_nt_0.8.csv", sep = ",", header = TRUE)


###memory groesse reicht nicht, deshalb muss das limit angehoben werden
memory.limit(size=56000)

apmg.lineage <- plyr::join(apmg, apmg_lin, by = "taxID")
bhv.lineage <- plyr::join(bhv, bhv_lin, by = "taxID")
reseq.lineage <- plyr::join(reseq, reseq_lin, by = "taxID")


# load metadata containing sample name, age and depth
metadata_apmg <- openxlsx::read.xlsx(xlsxFile = "APMG3738_Lama_all_shotgun_sample_name_change.xlsx")
metadata_bhv <- openxlsx::read.xlsx(xlsxFile = "BHV1_shotgun_sample_name_change.xlsx")
metadata_apmg$samples %in% apmg$samples # check if there are any missing samples
metadata_bhv$samples %in% bhv$samples # check if there are any missing samples

metadata_reseq <- openxlsx::read.xlsx(xlsxFile = "BHV2_shotgun_sample_name_change.xlsx")
metadata_reseq$samples %in% reseq$samples # check if there are any missing samples

##
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
#save(shotgun_all, file = "lama_shotgun_join_0.8.rdata")


# bacteria abundance as stratigraphic plot ----------------------------------
# filter to bacteria superkingdom
bacteria <- filter(shotgun_all, type == "sample" & superkingdom == "Bacteria")
bact_g_s <- filter(bacteria, Rank == "G" | Rank == "S")  # extract at the genus and species rank

selected.bact_g_s <- dplyr::select(bact_g_s, c("age", "Name", "TaxCount"))

###remove rows with zero counts in TaxCount
#selected.bact_g_s <- selected.bact_g_s[!(selected.bact_g_s$TaxCount=="0"),]
selected.bact_g_s$Name <- stringr::str_replace_all(selected.bact_g_s$Name,"_", " ") 

# count records per species
#species_counts <- bact_g_s %>%
 # group_by(Name) %>%
  #tally

#head(species_counts)

# get names of the species with counts >= 3
#frequent_species <-  species_counts %>%
 # filter(n >= 20) 

# filter out the less-frequent species
#bact_g_s_occ20 <- bact_g_s %>%
 # filter(bact_g_s$Name %in% frequent_species$Name)


##for contaminant selection
gen_spe_unique <- unique(bact_g_s_occ20$Name)
#write.csv2(gen_spe_unique, "unique_bacteria_gen_spe.csv")

bact_contam <- read.delim2("unique_bacteria_gen_spe_count20_assigned.csv", header = T, sep = ";")
bact_contam$Name <- stringr::str_replace_all(bact_contam$Name,"_", " ") 

bact_gen_spe_assign <- merge(selected.bact_g_s, bact_contam, by = "Name", all = TRUE)

##drop aquatic , NA and contaminant
bact_gen_spe_assign_clean <- bact_gen_spe_assign[!(bact_gen_spe_assign$element_cycle == "aquatic" | bact_gen_spe_assign$element_cycle == "contaminant" |
                                                     bact_gen_spe_assign$element_cycle == "NA" |bact_gen_spe_assign$element_cycle == ""),]


###make new Name_assigned with Name, assignment, real_5
bact_gen_spe_assign_clean$Name_assigned <- paste(bact_gen_spe_assign_clean$Name, bact_gen_spe_assign_clean$element_cycle, bact_gen_spe_assign_clean$real_5, sep ="_")
selected.columns <- dplyr::select(bact_gen_spe_assign_clean, c("age", "Name_assigned", "TaxCount")) # select columns to convert
#write.csv2(selected.columns, "bacteria_age_count_clean.csv")
wide_bact <- pivot_wider(selected.columns, names_from = Name_assigned, names_sep = "_", values_from = TaxCount, values_fn = sum) # make the data from long to wide
wide_bact[is.na(wide_bact)] <- 0 # replace NA with zero
wide_bact <- wide_bact[-45,] ##irgendwie ist zus?tzliche reihe hier aufgetaucht

###table for resampling
resampling_bact <- as.data.frame(wide_bact[, -1])
resampling_bact_occ3 <- chooseTaxa(resampling_bact, n.occ = 3)
resampling_bact_occ3$age <- wide_bact$age ###bring das alter wieder dazu
resampling_bact_occ3 <- resampling_bact_occ3 %>%
 dplyr::select(age, everything())#alter soll erste column sein

write.csv2(resampling_bact_occ3, "bacteria_assigned_resampling_occ3_taxCount.csv")
