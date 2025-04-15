### lama shotgun plotting
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(tidypaleo)
require(gridExtra)
library(grid)
library(gtable)

setwd("../Figure_6")

### load fungi data
fungi_resampl <- read.delim("../02_resampling/2023-03-08_fungi_clean_occ3_median_resampled_resampled_specieslevel_Sampleeffort274_aggregated_pcainput.csv", sep = ";", dec = ",")

### convert into parameter-long form
long.convert_fungi_res <- fungi_resampl %>% 
  gather(-X, key = Name, value = percentage) %>%
  group_by(X) %>%
  ungroup()

### sort the taxa by relative abundance
long_fungi_res <- arrange(long.convert_fungi_res, desc(percentage))
### force the taxa order
long_fungi_res$taxa <- factor(long_fungi_res$Name, levels = unique(long_fungi_res$Name))

### split the assigned name into three columns
long_fungi_res[c('Name', 'assignment', "real_5")] <- stringr::str_split_fixed(long_fungi_res$Name, '_', 3)
### real percentage is count/resampling count*100
long_fungi_res$real_percent <- (long_fungi_res$percentage)/275*100 


### extract only the selected fungi
long_fungi_weath_gg <- long_fungi_res %>%
  filter(stringr::str_detect(assignment, "lichen|mycorrhizae"))

### show age as ka
long_fungi_weath_gg$age <- as.numeric(long_fungi_weath_gg$X) / 1000
### check the oldest age and adjust the y axis in the plot
max.age <- round(max(long_fungi_weath_gg$age)) 

### merge genera together
long_fungi_weath_gg$Name_clean <- gsub("\\..*", "", long_fungi_weath_gg$Name)
long_fungi_weath_gg <- long_fungi_weath_gg %>%
  group_by(assignment, age) %>%
  summarise(merged_percent = (sum(real_percent)))

weath_gg <- ggplot(long_fungi_weath_gg, aes(x = merged_percent, y = age), fill = assignment) +
  coord_flip() +
  geom_areah() + 
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) +
  facet_grid(assignment ~ ., scales = "free", space = "fixed") + 
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 1))) +
  xlab(paste0("Relative abundance (%)")) +
  theme(panel.background = element_blank())+ theme(axis.line.x = element_line(color="black", size = 0.5),
                                                   axis.line.y = element_line(color="black", size = 0.5), legend.position = "none") 

### delete those with pH "out" or rename
fungi_pH_clean_res <- long_fungi_res[!(long_fungi_res$real_5 == "out" | long_fungi_res$real_5 == "pH tolerant" | long_fungi_res$real_5 == "unknown"),]

fungi_pH_gg_res <- fungi_pH_clean_res %>%
  group_by(X) %>%
  mutate(pH_percent= percentage/sum(percentage)*100)

fungi_pH_gg_res <- fungi_pH_gg_res %>%
  group_by(real_5, X) %>%
  summarise(ecol_percent = (sum(pH_percent)))
fungi_pH_clean_res$real_5 <- factor(fungi_pH_clean_res$real_5, levels = unique(fungi_pH_clean_res$real_5))

fungi_pH_gg_res$X <- as.numeric(fungi_pH_gg_res$X) / 1000 # show age as ka

fungi_pH_gg_res$levels <- ordered(fungi_pH_gg_res$real_5, levels=c(5, 4, 3, 2, 1))

### ph only alkaline and acidic
fungi_pH_gg_res$real_5 <- stringr::str_replace_all(fungi_pH_gg_res$real_5,"5", "4")
fungi_pH_gg_res$levels <- stringr::str_replace_all(fungi_pH_gg_res$levels,"5", "4")
fungi_pH_gg_res_acAl <- fungi_pH_gg_res[!(fungi_pH_gg_res$real_5 == "2" | fungi_pH_gg_res$real_5 == "3"),]
fungi_pH_gg_res_acAl <- fungi_pH_gg_res_acAl[,-4]

fungi_pH_gg_res_acAl <- fungi_pH_gg_res_acAl %>%
  group_by(real_5, X) %>%
  summarise(across(c(ecol_percent), sum))

fungi_pH_plot <- ggplot(fungi_pH_gg_res_acAl, aes(x = ecol_percent, y = X, fill = real_5))+
  coord_flip() +
  geom_areah() + 
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) + # exaggeration by 5
  facet_grid(real_5 ~ ., scales = "free", space = "fixed") + # facet by taxon
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 1))) + # reverse the y axis for age
  xlab(paste0("Relative abundance (%)")) +
  scale_fill_viridis_d() +
  #theme_bw() +
  theme(panel.background = element_blank())+ theme(axis.line.x = element_line(color="black", size = 0.5),
                                                   axis.line.y = element_line(color="black", size = 0.5), legend.position = "none")

### print the plot
print(fungi_pH_plot)

### Load bacterial data
bact_resampl <- read.delim("../02_resampling/2023-03-08_bacteria_clean_occ3_median_resampled_resampled_specieslevel_Sampleeffort24285_aggregated_pcainput.csv", sep = ";", dec = ",")

### convert into parameter-long form
long.convert_bact_res <- bact_resampl %>% 
  gather(-X, key = Name, value = percentage) %>%
  group_by(X) %>%
  ungroup()
### sort the taxa by relative abundance
long_bact_res <- arrange(long.convert_bact_res, desc(percentage))
### force the taxa order
long_bact_res$taxa <- factor(long_bact_res$Name, levels = unique(long_bact_res$Name))

### split the assigned name into three columns
long_bact_res[c('Name', 'assignment', "real_5")] <- stringr::str_split_fixed(long_bact_res$Name, '_', 3)
long_bact_res$real_percent <- (long_bact_res$percentage)/24285*100

### rename some of the initial assignments to narrow them down
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"As..Sb", "As")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"S..As", "As")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Au", "metalls")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C..Cl", "C..halogene")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Fe..S", "Fe")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Fe..Mn", "Fe")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"Mn", "C, Mn")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C, C, Mn", "C, Mn")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"symbiosis", "PGPB")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C, Mn", "C")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C, Mo", "C")
long_bact_res$assignment <- stringr::str_replace_all(long_bact_res$assignment,"C..C", "C")

## arsenic
as_gg_res <- long_bact_res %>%
  group_by(X) %>%
  mutate(element_percent= percentage/sum(percentage)*100)

as_gg_res <- as_gg_res %>%
  group_by(assignment, X) %>%
  summarise(tot_element_percent = (sum(element_percent)))
as_gg_res$assignment <- factor(as_gg_res$assignment, levels = unique(as_gg_res$assignment))

as_gg_res$X <- as.numeric(as_gg_res$X) / 1000 # show age as ka

as_gg_res_plot <- as_gg_res[(as_gg_res$assignment == "As" ),]

as_gg <- ggplot(as_gg_res_plot, aes(x = tot_element_percent, y = X, fill = assignment))+
  coord_flip() +
  geom_areah() + 
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) + # exaggeration by 5
  facet_grid(assignment ~ ., scales = "free", space = "fixed") + # facet by taxon
  #scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 1))) + # reverse the y axis for age
  xlab(paste0("Relative abundance (%)")) +
  #scale_fill_viridis_d() +
  #theme_bw() +
  theme(panel.background = element_blank())+ theme(axis.line.x = element_line(color="black", size = 0.5),
                                                   axis.line.y = element_line(color="black", size = 0.5), legend.position = "none") 


### delete those with pH "out" or rename
bact_pH_clean_res <- long_bact_res[!(long_bact_res$real_5 == "out" | long_bact_res$real_5 == "" | long_bact_res$real_5 == "unknown"),]

bact_pH_gg_res <- bact_pH_clean_res %>%
  group_by(X) %>%
  mutate(pH_percent= percentage/sum(percentage)*100)

bact_pH_gg_res <- bact_pH_gg_res %>%
  group_by(real_5, X) %>%
  summarise(ecol_percent = (sum(pH_percent)))
bact_pH_gg_res$real_5 <- factor(bact_pH_gg_res$real_5, levels = unique(bact_pH_gg_res$real_5))

bact_pH_gg_res$X <- as.numeric(bact_pH_gg_res$X) / 1000 # show age as ka

bact_pH_gg_res$levels <- ordered(bact_pH_gg_res$real_5, levels=c(5,4,3,2,1))

### pH only alkaline and acidic
bact_pH_gg_res$real_5 <- stringr::str_replace_all(bact_pH_gg_res$real_5,"5", "4")
bact_pH_gg_res$levels <- stringr::str_replace_all(bact_pH_gg_res$levels,"5", "4")
bact_pH_gg_res_acAl <- bact_pH_gg_res[!(bact_pH_gg_res$real_5 == "2" | bact_pH_gg_res$real_5 == "3"),]
bact_pH_gg_res_acAl <- bact_pH_gg_res_acAl[,-4]

bact_pH_gg_res_acAl <- bact_pH_gg_res_acAl %>%
  group_by(real_5, X) %>%
  summarise(across(c(ecol_percent), sum))

bact_pH_plot <- ggplot(bact_pH_gg_res_acAl, aes(x = ecol_percent, y = X, fill = real_5))+
  coord_flip() +
  geom_areah() + 
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) + # exaggeration by 5
  facet_grid(real_5 ~ ., scales = "free", space = "fixed") + # facet by taxon
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 1))) + # reverse the y axis for age
  xlab(paste0("Relative abundance (%)")) +
  scale_fill_viridis_d() +
  #theme_bw() +
  theme(panel.background = element_blank())+ theme(axis.line.x = element_line(color="black", size = 0.5),
                                                   axis.line.y = element_line(color="black", size = 0.5), legend.position = "none") 

### nutrient cyclers
bact_nutr <- long_bact_res[!(long_bact_res$assignment == "As" | long_bact_res$assignment == "Fe" | long_bact_res$assignment == "unknown" |
                               long_bact_res$assignment == "lichen" | long_bact_res$assignment == "metalls" | long_bact_res$assignment == "PGPB" 
),]

### merge further the assignments
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C..Fe", "C")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C..halogene", "C")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C..Mo", "C")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C..N..Cl|C..N..Mo", "C, N")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C..N", "C, N")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C, N..P", "C, N, P")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C..P", "C, P")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"C..S", "C, S")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"N..P", "N, P")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"N..S", "N, S")
bact_nutr$assignment <- stringr::str_replace_all(bact_nutr$assignment,"P..C..S..Cl..oxalates..silicates", "C, P, S")

nutr_gg_res <- bact_nutr %>%
  group_by(X) %>%
  mutate(nutr_percent= percentage/sum(percentage)*100)

nutr_gg_res <- nutr_gg_res %>%
  group_by(assignment, X) %>%
  summarise(ecol_percent = (sum(nutr_percent)))
nutr_gg_res$assignment <- factor(nutr_gg_res$assignment, levels = unique(nutr_gg_res$assignment))

nutr_gg_res$X <- as.numeric(nutr_gg_res$X) / 1000

nutr_gg_res_plot <- nutr_gg_res[(nutr_gg_res$assignment == "C" | nutr_gg_res$assignment == "N" | nutr_gg_res$assignment == "S"
),]

nutr_gg <- ggplot(nutr_gg_res_plot, aes(x = ecol_percent, y = X, fill = assignment))+
  coord_flip() +
  geom_areah() + 
  geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) + # exaggeration by 5
  facet_grid(assignment ~ ., scales = "free", space = "fixed") + # facet by taxon
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 1))) + # reverse the y axis for age
  xlab(paste0("Relative abundance (%)")) +
  scale_fill_viridis_d() +
  #theme_bw() +
  theme(panel.background = element_blank())+ theme(axis.line.x = element_line(color="black", size = 0.5),
                                                   axis.line.y = element_line(color="black", size = 0.5), legend.position = "none") 

### load plant data 
plant_resampl <- read.delim("../02_resampling/2023-03-06_plants_median_resampled_resampled_specieslevel_Sampleeffort3533_aggregated_pcainput.csv", sep = ";", dec = ",")

### convert into parameter-long form
long.convert_plant_res <- plant_resampl %>%
  gather(-X, key = Name, value = percentage) %>%
  group_by(X) %>%
  ungroup()

### sort the taxa by relative abundance
long_plant_res <- arrange(long.convert_plant_res, desc(percentage))
### force the taxa order
long_plant_res$taxa <- factor(long_plant_res$Name, levels = unique(long_plant_res$Name))

### split the assigned name into three columns
long_plant_res[c('Name', 'assignment', "real_5")] <- stringr::str_split_fixed(long_plant_res$Name, '_', 3)
long_plant_res$real_percent <- (long_plant_res$percentage)/3533*100 ###real percentage is count/resampling count*100

### delete those with pH "out" or rename
plant_pH_clean_res <- long_plant_res[!(long_plant_res$real_5 == "out" | long_plant_res$real_5 == "" | long_plant_res$real_5 == "unknown"
                                       | long_plant_res$real_5 == "pH.tolerant"),]

plant_pH_gg_res <- plant_pH_clean_res %>%
  group_by(X) %>%
  mutate(pH_percent= percentage/sum(percentage)*100)

plant_pH_gg_res <- plant_pH_gg_res %>%
  group_by(real_5, X) %>%
  summarise(ecol_percent = (sum(pH_percent)))
plant_pH_gg_res$real_5 <- factor(plant_pH_gg_res$real_5, levels = unique(plant_pH_gg_res$real_5))

plant_pH_gg_res$X <- as.numeric(plant_pH_gg_res$X) / 1000 # show age as ka

plant_pH_gg_res$levels <- ordered(plant_pH_gg_res$real_5, levels=c(5,4,3,2,1))

### ph only alkaline and acidic
plant_pH_gg_res$real_5 <- stringr::str_replace_all(plant_pH_gg_res$real_5,"5", "4")
plant_pH_gg_res$levels <- stringr::str_replace_all(plant_pH_gg_res$levels,"5", "4")

plant_pH_clean_res_acAl <- plant_pH_gg_res[!(plant_pH_gg_res$real_5 == "2" | plant_pH_gg_res$real_5 == "3"),]

plant_pH_clean_res_acAl <- plant_pH_clean_res_acAl[,-4]

plant_pH_clean_res_acAl <- plant_pH_clean_res_acAl %>%
  group_by(real_5, X) %>%
  summarise(across(c(ecol_percent), sum))

plant_pH_plot <- ggplot(plant_pH_clean_res_acAl, aes(x = ecol_percent, y = X, fill = real_5))+
  coord_flip() +
  geom_areah() + 
  #geom_lineh_exaggerate(exaggerate_x = 5, col = "grey70", lty = 2, linewidth = 0.6) + # exaggeration by 5
  facet_grid(real_5 ~ ., scales = "free", space = "fixed") + # facet by taxon
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 1))) + # reverse the y axis for age
  xlab(paste0("Relative abundance (%)")) +
  scale_fill_viridis_d() +
  #theme_bw() +
  theme(panel.background = element_blank())+ theme(axis.line.x = element_line(color="black", size = 0.5),
                                                   axis.line.y = element_line(color="black", size = 0.5), legend.position = "none") 

### merge all selected assignments and import xrf
xrf_new <- read.delim("2022-08-15_xrf_for_smoothing.csv", sep=";")
age <- read.delim("PG1341_woody_seg6_minus145cm_boundary_168_ages.txt", sep="\t",  header=TRUE, stringsAsFactors=FALSE, dec=",")

### combine age and xrf_clean data
xrf_age <- merge(xrf_new, age, by = "depth")
y <- xrf_age$mean ### smoothing of all values 

K_l <- xrf_age$K
Ti_l <- xrf_age$Ti

K_loess <- predict(loess(y ~ K_l, xrf_age, span = 0.75), se = T)
Ti_loess <- predict(loess(y ~ Ti_l, xrf_age, span = 0.75), se = T)

loess_xrf <- data.frame(y, K_loess$fit, Ti_loess$fit)
loess_xrf <- loess_xrf[!(loess_xrf$y=="-47"),] 
loess_xrf <- loess_xrf[!(loess_xrf$y=="-7"),] 
loess_xrf$facet <- "weathering"

K_Ti <- ggplot(data=loess_xrf, aes(y=y/1000, x = (K_loess.fit/Ti_loess.fit))) +
  coord_flip() +
  geom_lineh(group = 1, colour = "black") +
  geom_smooth(group = 1, colour = "darkred", orientation = "y") +
  facet_grid(facet ~. , scales = "free", space = "fixed") +
  scale_y_reverse(name = "Age (ka)", breaks = rev(seq(0, max.age, by = 1))) + # reverse the y axis for age
  xlab(paste0("K/Ti")) +
  #theme_bw() +
  theme(panel.background = element_blank())+ theme(axis.line.x = element_line(color="black", size = 0.5),
                                                   axis.line.y = element_line(color="black", size = 0.5), legend.position = "none") 


### Convert ggplot objects to gtables
K_Ti_grob <- ggplotGrob(K_Ti)
weath_gg_grob <- ggplotGrob(weath_gg)
as_gg_grob <- ggplotGrob(as_gg)
nutr_gg_grob <- ggplotGrob(nutr_gg)
plant_pH_plot_grob <- ggplotGrob(plant_pH_plot)
fungi_pH_plot_grob <- ggplotGrob(fungi_pH_plot)
bact_pH_plot_grob <- ggplotGrob(bact_pH_plot)

### Combine grobs vertically
all_combined <- gtable_rbind(K_Ti_grob, 
                             as_gg_grob,
                             weath_gg_grob,
                             nutr_gg_grob, 
                             plant_pH_plot_grob, 
                             fungi_pH_plot_grob, 
                             bact_pH_plot_grob)

grid.newpage()
grid.draw(all_combined)
svg(filename = "2023-04-06_all_weathering2.svg", width = 8, height = 20)
grid.draw(all_combined)
dev.off()
