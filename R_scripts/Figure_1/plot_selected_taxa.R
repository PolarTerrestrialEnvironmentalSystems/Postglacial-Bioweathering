###lama shotgun plotting

setwd("Figure_I")
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tidypaleo)

################
#######
################
###trial mit geresampled
fungi_resampl <- read.delim("Resampling/Output/2023-03-07_fungi_occ3_median_resampled_resampled_specieslevel_Sampleeffort275_aggregated_pcainput.csv", sep = ";", dec = ",")

long.convert_fungi_res <- fungi_resampl %>% # convert into parameter-long form
  gather(-X, key = Name, value = percentage) %>%
  group_by(X) %>%
  ungroup()
long_fungi_res <- arrange(long.convert_fungi_res, desc(percentage)) # sort the taxa by relative abundance
long_fungi_res$taxa <- factor(long_fungi_res$Name, levels = unique(long_fungi_res$Name)) # force the taxa order

####split the assigned name into three columns
long_fungi_res[c('Name', 'assignment', "real_5")] <- stringr::str_split_fixed(long_fungi_res$Name, '_', 3)
long_fungi_res$real_percent <- (long_fungi_res$percentage)/275*100 ###real percentage is count/resampling count*100

######extract only the mentioned fungi
long_fungi_res_gg <- long_fungi_res %>%
  filter(stringr::str_detect(Name, "Ascomycota|Penicillium|Mortierella|Malassezia|Komagataella|Verrucariales|Peltigera|Suillineae|Glomeraceae|Rhizophagus|Laccaria|Hyaloscypha|Tuberaceae"))

long_fungi_res_gg$age <- as.numeric(long_fungi_res_gg$X) / 1000 # show age as ka
max.age <- round(max(long_fungi_res_gg$age)) # check the oldest age and adjust the y axis in the plot

###merge genera together
long_fungi_res_gg$Name_clean <- gsub("\\..*", "", long_fungi_res_gg$Name)        # Apply gsub with \\
long_fungi_res_gg$Name_clean[long_fungi_res_gg$Name_clean == 'unclassified'] <- 'Mortierella'
long_fungi_res_gg <- long_fungi_res_gg %>%
  group_by(Name_clean, age) %>%
  summarise(merged_percent = (sum(real_percent)))


################
#######
################
###trial mit geresampled
bact_resampl <- read.delim("Resampling/Output/2023-03-08_bacteria_clean_occ3_median_resampled_resampled_specieslevel_Sampleeffort24285_aggregated_pcainput.csv", sep = ";", dec = ",")

long.convert_bact_res <- bact_resampl %>% # convert into parameter-long form
  gather(-X, key = Name, value = percentage) %>%
  group_by(X) %>%
  ungroup()
long_bact_res <- arrange(long.convert_bact_res, desc(percentage)) # sort the taxa by relative abundance
long_bact_res$taxa <- factor(long_bact_res$Name, levels = unique(long_bact_res$Name)) # force the taxa order

####split the assigned name into three columns
long_bact_res[c('Name', 'assignment', "real_5")] <- stringr::str_split_fixed(long_bact_res$Name, '_', 3)
long_bact_res$real_percent <- (long_bact_res$percentage)/24285*100 ###real percentage is count/resampling count*100

long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"As..Sb", "As") ##As, Sb to As
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"S..As", "As") ##S, As to As
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Au", "metalls") ##Au to metalls
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C..Cl", "C..halogene")# C, Cl to C, halogene
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Fe..S", "Fe")# Fe, S to Fe
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Fe..Mn", "Fe")# Fe, Mn to Fe
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Mn", "C, Mn")# Mn to C, Mn
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C, C, Mn", "C, Mn")# Mn to C, Mn
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"symbiosis", "PGPB")# symbiosis to PGPB
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C, Mn", "C")# symbiosis to PGPB
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C, Mo", "C")# symbiosis to PGPB
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C..C", "C")# symbiosis to PGPB



##############
######extract only the mentioned bacteria
long_bact_res_gg <- long_bact_res %>%
  filter(stringr::str_detect(Name, "Brevundimonas|Hydrogenophaga|Herminiimonas|Bradyrhizobium|Ferrigenium|Sideroxydans|Pseudolabrys|Delftia"))

long_bact_res_gg$age <- as.numeric(long_bact_res_gg$X) / 1000 # show age as ka
max.age <- round(max(long_bact_res_gg$age)) # check the oldest age and adjust the y axis in the plot

###merge genera together
long_bact_res_gg$Name_clean <- gsub("\\..*", "", long_bact_res_gg$Name)        # Apply gsub with \\
long_bact_res_gg <- long_bact_res_gg %>%
  group_by(Name_clean, age) %>%
  summarise(merged_percent = (sum(real_percent)))

##########plants resampled
plant_resampl <- read.delim("Resampling/Output/2023-03-06_plants_median_resampled_resampled_specieslevel_Sampleeffort3533_aggregated_pcainput.csv", sep = ";", dec = ",")

long.convert_plant_res <- plant_resampl %>% # convert into parameter-long form
  gather(-X, key = Name, value = percentage) %>%
  group_by(X) %>%
  ungroup()
long_plant_res <- arrange(long.convert_plant_res, desc(percentage)) # sort the taxa by relative abundance
long_plant_res$taxa <- factor(long_plant_res$Name, levels = unique(long_plant_res$Name)) # force the taxa order

####split the assigned name into three columns
long_plant_res[c('Name', 'assignment', "real_5")] <- stringr::str_split_fixed(long_plant_res$Name, '_', 3)
long_plant_res$real_percent <- (long_plant_res$percentage)/3533*100 ###real percentage is count/resampling count*100

######extract only the mentioned plants
long_plant_res_gg <- long_plant_res %>%
  filter(stringr::str_detect(Name, "Dryas|Saxifragaceae|Salix|Betula|Alnus|Pinaceae|Larix|Picea|Ericaceae|
                             Asteraceae|Pyrola rotundifolia|Vaccinium|Pyrola|Saxifraga|Potentilla"))

long_plant_res_gg$age <- as.numeric(long_plant_res_gg$X) / 1000 # show age as ka
max.age <- round(max(long_plant_res_gg$age)) # check the oldest age and adjust the y axis in the plot

###merge genera together
long_plant_res_gg$Name_clean <- gsub("\\..*", "", long_plant_res_gg$Name)        # Apply gsub with \\
long_plant_res_gg <- long_plant_res_gg %>%
  group_by(Name_clean, age) %>%
  summarise(merged_percent = (sum(real_percent)))

#######gg all
fungi_all_selected_input <- long_fungi_res_gg
plant_all_selected_input <- long_plant_res_gg
bact_all_selected_input <- long_bact_res_gg

fungi_all_selected_input$ecology <- "fungi"
plant_all_selected_input$ecology <- "plants"
bact_all_selected_input$ecology <- "bacteria"

fungi_all_selected_input$order <- "2"
plant_all_selected_input$order <- "1"
bact_all_selected_input$order <- "3"

all_selected_gg <- rbind(fungi_all_selected_input, plant_all_selected_input, bact_all_selected_input)

####new_order
new_order <- read.delim("Figure_I/ggplot_order_new.csv", sep = ";")

####
all_selected_gg_order <- merge(all_selected_gg, new_order, by = "Name_clean")

####force order
my_taxa_order=c(unique(plant_all_selected_input$Name_clean), unique(fungi_all_selected_input$Name_clean), unique(bact_all_selected_input$Name_clean)) # it is a vector with unique taxa from plants to fungi to bacteria.
my_taxa_order2 = all_selected_gg_order$order_new

all_selected_gg$my_taxa_order = factor(all_selected_gg$Name_clean, levels=my_taxa_order)

#added new line (144) because this should be unique !
my_taxa_order2 = unique(all_selected_gg_order$order_new)
all_selected_gg_order$my_taxa_order_new = factor(all_selected_gg_order$Name_clean, levels=my_taxa_order2)

all_selected_gg_order <- arrange(all_selected_gg_order, order_new)
all_selected_gg_order$Name_clean <- factor(all_selected_gg_order$Name_clean, levels = unique(all_selected_gg_order$Name_clean))

all_gg <- ggplot(all_selected_gg_order, aes(x = merged_percent, y = age), fill = ecology) +
  geom_areah_exaggerate(data = all_selected_gg, exaggerate_x = 5, fill = "lightgrey") +
  geom_areah(aes(fill = ecology)) + 
  scale_fill_manual(values = c("#1C629D", "#784C23", "#748700")) +
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) + # exaggeration by 5
  facet_grid(~ Name_clean, scales = "free", space = "fixed") + # facet by taxon
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 5))) + # reverse the y axis for age
  xlab(paste0("Relative abundance (%) of most important taxa")) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 14, angle = 90),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        panel.spacing.x = unit(2, "mm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks.length = unit(2, "mm"))
svg("2023-10-23_stratplot_lama_nt22_0.8_selected_all_merged_resampled_all_ordered_final_exx_trial_new.svg", width = 30, height = 6, pointsize = 8)
plot(all_gg)
dev.off()

#### ordered
all_gg <- ggplot(all_selected_gg_order, aes(x = merged_percent, y = age), fill = ecology) +
  #coord_flip() + # to create horizontal plots
  geom_areah_exaggerate(data = all_selected_gg_order, exaggerate_x = 5, fill = "lightgrey") +
  geom_areah(aes(fill = ecology)) + 
  scale_fill_manual(values = c("#1C629D", "#784C23", "#748700")) +
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) + # exaggeration by 5
  facet_grid(~ Name_clean, scales = "free", space = "fixed") + # facet by taxon
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 5))) + # reverse the y axis for age
  xlab(paste0("Relative abundance (%) of most important taxa")) +
  theme_bw() +
  theme(strip.text.x = element_text(size = 14, angle = 90),
        strip.background = element_blank(),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        panel.spacing.x = unit(2, "mm"),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks.length = unit(2, "mm"))
svg("2023-10-26_stratplot_lama_nt22_0.8_selected_all_merged_resampled_all_ordered_final_exx_trial_coloured_new.svg", width = 30, height = 6, pointsize = 8)
plot(all_gg)
dev.off()
