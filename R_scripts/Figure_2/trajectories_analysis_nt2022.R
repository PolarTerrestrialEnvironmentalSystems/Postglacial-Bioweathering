#########trajectories analysis
setwd("R_script/Figure_2")

library(vegan)
library(BiodiversityR)
library(ppcor)
library(ggplot2)
library(graphics)
#install.packages("VennDiagram")   # Install & load VennDiagram package
library("VennDiagram")

vegetation <- read.delim("2023-03-06_plants_median_resampled_resampled_specieslevel_Sampleeffort3533_aggregated_pcainput.csv", sep=";", header=TRUE, stringsAsFactors=FALSE, dec=",")
fungi <- read.delim("2023-03-07_fungi_occ3_median_resampled_resampled_specieslevel_Sampleeffort275_aggregated_pcainput.csv", sep = ";", dec = ",")
bacteria <- read.delim("2023-03-08_bacteria_clean_occ3_median_resampled_resampled_specieslevel_Sampleeffort24285_aggregated_pcainput.csv", sep = ";", dec = ",")
temperature <- read.delim("tagg_temperatures_paleo.csv", sep = ";", dec = ",")
time <- read.delim("time_variable.txt", sep = "\t")

########
# the temperature needs the same format as the other variables
names(vegetation)[1] <- "sample" ##rename first column of vegetation
names(fungi)[1] <- "sample" ##rename first column of fungi
names(bacteria)[1] <- "sample" ##rename first column of bacteria

temp_age <- merge(temperature, vegetation, by = "sample", all = TRUE) 
y <- temp_age$sample
tr <- temp_age$tr

vegetation2 <- vegetation[,-1]
rownames(vegetation2) <- vegetation[,1]
rowsumsnotzero=which(apply(vegetation2,1,sum)>0)
colsumsnotzero=which(apply(vegetation2,2,sum)>0)
vegetation2=vegetation2[rowsumsnotzero,colsumsnotzero]

fungi2 <- fungi[,-1]
rownames(fungi2) <- fungi[,1]
fungi2=fungi2[rowsumsnotzero,colsumsnotzero]

bacteria2 <- bacteria[,-1]
rownames(bacteria2) <- bacteria[,1]
bacteria2=bacteria2[rowsumsnotzero,colsumsnotzero]

time2 <- time[,-1]

###interpolation of the missing tr values
library("imputeTS")
temp_int <- na_interpolation(temp_age)

temp_final <- merge(temp_int, bacteria, all.y = TRUE)
temp_final <- temp_final[,-c(3:1421)]

rownames(temp_final) <- temp_final[,1]
temp_input <- temp_final[,-1]


veg.pca <- rda(sqrt(sqrt(vegetation2)))

PCAsignificance(veg.pca, axes = 8) ##first two pc-axis are important

veg.pca.site1 = as.data.frame(veg.pca$CA$u)[,1]
veg.pca.site2 = as.data.frame(veg.pca$CA$u)[,2]

veg.pca.site = as.data.frame(veg.pca$CA$u)[,1:2]

pc1veg=sort(veg.pca$CA$v[,"PC1"])
par(mar=c(4,12,2,2),las=1);barplot(c(pc1veg[1:10],rev(pc1veg)[1:10]), horiz=TRUE)

pc2veg=sort(veg.pca$CA$v[,"PC2"])
par(mar=c(4,12,2,2),las=1);barplot(c(pc2veg[1:10],rev(pc2veg)[1:10]), horiz=TRUE)


########variation partitioning https://r.qcbs.ca/workshop10/book-en/variation-partitioning.html#variation-partitioning-in-r
##fungi with vegetation , temperature variation and time
fungi_vp <- varpart(sqrt(sqrt(fungi2)), veg.pca.site, temp_input, time2)
fungi_vp$part ###access results: Total variation (SS): 838.14, Variance: 19.492

###significance of the variables
anova.cca(rda(fungi2, veg.pca.site)) #Pr(>F) 0.001 ***
anova.cca(rda(fungi2, time2)) #Pr(>F) 0.001 ***
anova.cca(rda(fungi2, temp_input)) #Pr(>F) 0.001 ***
anova.cca(rda(fungi2, veg.pca.site, temp_input)) #Pr(>F) 0.001 *** 
anova.cca(rda(fungi2, veg.pca.site, time2)) #Pr(>F) 0.001 ***
anova.cca(rda(fungi2, time2, temp_input)) #Pr(>F) 0.049 *

# plot the variation partitioning Venn diagram
fungi_varpar <- plot(fungi_vp,
                     Xnames = c("Vegetation", "Temperature", "Time"), # name the partitions
                     bg = c("seagreen3", "yellow", "blue", "orange"), 
                     alpha = 50, # colour the circles
                     digits = 2, # only show 2 digits
                     cex = 1.5,
                     lty = 1,
                     bty = "n",
                     pty = "m")

ggsave("2023-03-08_fungi_var_par_withTime_vegeMerged.svg", width = 10, height = 10)

#####bacteria
###bacteria vege merged
bact_vp <- varpart(sqrt(sqrt(bacteria2)), veg.pca.site, temp_input, time2)
bact_vp$part ###access results: Total variation (SS): 1181.5, Variance: 27.476

###significance
anova.cca(rda(bacteria2, veg.pca.site)) #Pr(>F) 0.002 **
anova.cca(rda(bacteria2, time2)) #Pr(>F) 0.004 **
anova.cca(rda(bacteria2, temp_input)) #Pr(>F) 0.001 ***
anova.cca(rda(bacteria2, veg.pca.site, temp_input)) #Pr(>F) 0.005 **
anova.cca(rda(bacteria2, veg.pca.site, time2)) #Pr(>F) 0.022 * 
anova.cca(rda(bacteria2, time2, temp_input)) #Pr(>F) 0.46

bact_varpar <- plot(venneuler(bact_vp,
                    Xnames = c("Vegetation", "Temperature", "Time"), # name the partitions
                    bg = c("seagreen3", "mediumpurple", "red", "orange"), alpha = 80, # colour the circles
                    digits = 2, # only show 2 digits
                    cex = 1.5))

svg("2023-03-08_bacteria_var_par_withTime_vegeMerged.svg", bact_varpar, width = 15, height = 10, pointsize = 8)
plot(bact_varpar)
dev.off()
