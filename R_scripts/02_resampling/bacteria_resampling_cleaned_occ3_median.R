library(tidyverse)
options(stringsAsFactors=FALSE)
### select main folder containing folders named exactly "data"
setwd("~/Postglacial-Bioweathering-main/R_scripts/Resampling/02_resampling")

rm(list=ls())
### data table was exported from xlsx here csv2: sep=";", dec=","
specseq=read.delim("bacteria_assigned_resampling_occ3.csv", header = TRUE, sep = ";")
specseq <- specseq[,-1]

### check data is of required format, count data as integer/double, taxa names as strings
str(specseq)

### make sure additional rows or NA values are converted to zeros
specseq[is.na(specseq)]=0

### dataframe still needs to be transformed
t_specseq <- as.data.frame(t(specseq), header = F) 
names(t_specseq) <- t_specseq[1,] 
t_specseq <- t_specseq[-1,] 
t_specseq$name <- row.names(t_specseq)
t_specseq <- t_specseq %>%
  dplyr::select(name, everything())

### define columns that contain raw count data
names(t_specseq)
colstart=2
colend=45
names(t_specseq)[colstart:colend]# check
COLUMNNAMESAREYEARS=TRUE

### species must have unique sample names and need to be merged with the family (for technical reasons)
SPECIESNAMECOLUMN=1# position of the species assignments column in the data frame
names(t_specseq)[SPECIESNAMECOLUMN]="scientific_name"# change species name column to "scientific_name"

FAMILYNAMEPRESENT=FALSE
FAMILYNAMECOLUMN=1# position of the family assignments column in the data frame
if(FAMILYNAMEPRESENT)
{
  names(t_specseq)[FAMILYNAMECOLUMN]="element_cycle"# change family name column to "family_name"
} else
{
  t_specseq$family_name="NoFamilyName"# add empty family name column
}

### make sure to have individual names for each species/taxa entry
specseq_final_name=paste(gsub(" ", "_", t_specseq$family_name), make.unique(t_specseq$scientific_name))

### get mean and median counts
sample_counts <- colSums(t_specseq[colstart:colend])

### 90639.05
mean(sample_counts) 
### 24285
median(sample_counts) 


### resample loop for each sample/year present in the data table
### determine min. read counts for rarefaction, here automatic procedure to find the minimum within the data set
nsampleff=24285
### set here the number of resamplings, standard==100
resamplingnumber=100
genrare=list()

### in case of column/sample names are numbers corresponding to depths, an X was given in front of the number on data table import, in the loop this X is deleted and the year is assigned as the name of the entry in the filled loop
yr=NULL 

### this counter is needed to reorganize the list in case samples without any reads are present. These will be skipped (should be not the case otherwise the minimal counts would be zero what will produce in an error so please check before running the script if you have such samples in your data set end exclude them prior analyses).
missing=0

for(yrcoli in colstart:colend)
{
  print(yrcoli)
  
  allspec=specseq_final_name[which(t_specseq[,yrcoli]>0)]
  allspec_counts=t_specseq[which(t_specseq[,yrcoli]>0),yrcoli]
  
  if(length(allspec)>0)
  {
    if(COLUMNNAMESAREYEARS)
    {# convert to number
      yr=c(yr, as.numeric(gsub(",", ".", gsub("X","",names(t_specseq)[yrcoli]))))
    } else
    {# take the name as it is
      yr=c(yr, names(t_specseq)[yrcoli])
    }
    
    sampleeffort=list()
    for(nsampleeffi in nsampleff)
    {
      repeatsample=list()
      for(repi in 1:resamplingnumber)
      {
        repeatsample[[repi]]=sample(allspec,nsampleeffi,replace=TRUE, prob=allspec_counts/sum(allspec_counts))
      }
      sampleeffort[[which(nsampleff==nsampleeffi)]]=repeatsample
    }
    genrare[[yrcoli-missing-(colstart-1)]]=list(allspec,sampleeffort)
  } else
  {
    missing=missing+1
  }
}
names(genrare)=yr

### how to access the data in a resampled dataset (list called genrare)
### ... show how many taxa in the first newly generated data set were drawn
length(unique(genrare[[1]][[2]][[1]][[1]]))

### ... genrare[[1...number of columns]]
### ... genrare[[1...number of columns]][[1]] == potential unique taxa names in the corresponding sample
### ... genrare[[1...number of columns]][[2]] == resampled data
### ... genrare[[1...number of columns]][[2]][[1]] == resampled data of nsampleff, usually only 1 present
### ... genrare[[1...number of columns]][[2]][[1]][[1:...resamplingnumber]] == resampled data (individual taxa) as determined earlier

### processing of the resampled data set
### count total species/family number and reads of individual families
famorig=specseq$family_name
familylevels=names(rev(sort(table(famorig))))

totspec=NULL
totfam=NULL
for(li in 1:length(genrare))
{
  for(li2 in 1:length(genrare[[li]][[2]]))
  {
    print(paste0(li," - ",li2))
    spectot=NULL
    spectot4fam=NULL
    for(repi in 1:resamplingnumber)
    {
      pei=unique(genrare[[li]][[2]][[li2]][[repi]])
      spectot=c(spectot,length(pei))
      spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split=" ",pei),function(x)return(x[1]))),levels=familylevels)))
    }
    totspec=rbind(totspec, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),Nspecies=spectot))
    totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
  }
}

### simple plots of the processed data
### modify column names ifnot years that can be coerced to numbers
if(!COLUMNNAMESAREYEARS)
{
  totspec$T=factor(totspec$T, levels=names(genrare))
  totfam$T=factor(totfam$T, levels=names(genrare))
}
png(paste0("resampling/2023-03-08_bacteria_clean_occ3_median_resampled_totalspecies_Sampleeffort",nsampleff,"_plot.png"), width=480,height=480)
par(mar=c(8,4,3,1),las=2)
with(totspec,plot(Nspecies~T, main="number of species per sample"))
dev.off()

### save processed data
write.csv2(totspec, paste0("output/2023-03-08_bacteria_clean_occ3_median_resampled_resampled_totalspecies_Sampleeffort",nsampleff,".csv"), row.names=FALSE)	
write.csv2(totfam, paste0("output/2023-03-08_bacteria_clean_occ3_median_resampled_totalfamilies_Sampleeffort",nsampleff,".csv"), row.names=FALSE)	



### aggregate data on species level
famorig=t_specseq$scientific_name
familylevels=names(rev(sort(table(famorig))))

totfam=NULL
for(li in 1:length(genrare))
{
  for(li2 in 1:length(genrare[[li]][[2]]))
  {
    print(paste0(li," - ",li2))
    
    spectot4fam=NULL
    for(repi in 1:resamplingnumber)
    {
      pei=genrare[[li]][[2]][[li2]][[repi]]
      spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split=" ",pei),function(x)return(paste(x[-1],collapse = " ")))),levels=familylevels)))
    }
    totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
  }
}

### modify column names ifnot years that can be coerced to numbers
if(!COLUMNNAMESAREYEARS)
{
  totfam$T=factor(totfam$T, levels=names(genrare))
}

### calculate mean values for each species/taxa
speciesfamiliesdf_totfam=NULL
pdf(paste0("output/2023-03-08_bacteria_clean_occ3_median_resampled_resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated.pdf"))
par(mar=c(8,4,3,1),las=2)
for(fami in names(totfam)[3:dim(totfam)[2]])
{
  plot(totfam[,fami]~totfam$T,col=rainbow(length(names(totfam)[3:dim(totfam)[2]]),s=0.6)[which(names(totfam)[3:dim(totfam)[2]]==fami)],type="n",lwd=2,main=fami,ylab="Sequence counts",xlab="Depth (m)")
  
  mn=aggregate(totfam[,fami],list(as.numeric(totfam$T)),mean)
  ti=mn$Group.1
  mn=mn$x
  sd=aggregate(totfam[,fami],list(as.numeric(totfam$T)),sd)$x*1.96
  
  polygon(y=c(mn+sd,rev(mn-sd)), x=c(ti,rev(ti)), col="gray", border=NA)
  lines(mn~ti,col="steelblue2",lwd=3)
  
  if(!COLUMNNAMESAREYEARS)
  {
    ti=levels(totfam$T)
  }
  speciesfamiliesdf_totfam=rbind(speciesfamiliesdf_totfam, data.frame(Species=fami,TBP=ti,Mean=mn,CI95=sd))
}
dev.off()

str(speciesfamiliesdf_totfam)

### save processed data
write.csv2(speciesfamiliesdf_totfam, paste0("output/2023-03-08_bacteria_clean_occ3_median_resampled_resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated.csv"), row.names=FALSE)

### post processing	
### reformat for pca input ... rows are samples, cols are species/taxa
ordidf=NULL
for(fi in levels(factor(speciesfamiliesdf_totfam$Species)))
{
  numberi=speciesfamiliesdf_totfam[speciesfamiliesdf_totfam$Species==fi, ]$Mean
  ordidf=rbind(ordidf, data.frame(Spec=fi,t(numberi)))
}
names(ordidf)[2:dim(ordidf)[2]]=speciesfamiliesdf_totfam[speciesfamiliesdf_totfam$Species==fi, ]$TBP
str(ordidf)

### remove species/taxa column
row.names(ordidf)=ordidf$Spec
ordidf=ordidf[,-1]
str(ordidf)
head(ordidf)

### exclude samples or species that have no records
rowsumsnotzero=which(apply(ordidf,1,sum)>0)
colsumsnotzero=which(apply(ordidf,2,sum)>0)
ordidf=ordidf[rowsumsnotzero,colsumsnotzero]

### export data
write.csv2(t(ordidf), paste0("output/2023-03-08_bacteria_clean_occ3_median_resampled_resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated_pcainput.csv"))

### comparison of original and resampled data set
ordiorigdf=t_specseq[colstart:colend]
row.names(ordiorigdf)=make.unique(t_specseq$scientific_name)
rowsumsnotzero=which(apply(ordiorigdf,1,sum)>0)
colsumsnotzero=which(apply(ordiorigdf,2,sum)>0)
ordiorigdf=ordiorigdf[rowsumsnotzero,colsumsnotzero]

### species count (counts >= 1)
png(paste0("output/2023-03-08_bacteria_clean_occ3_median_resampled_resampled_speciesnumber_Sampleeffort",nsampleff,"_aggregated_comparisonplot.png"), width=480,height=480)
par(mar=c(8,4,3,1),las=2)
barplot(apply(ordiorigdf,2,function(x)length(which(x>=1))), col="tomato", border=FALSE)
barplot(apply(ordidf,2,function(x)length(which(x>=1))), add=TRUE, col="skyblue", border=FALSE)
legend("topright", c("original","resampled"), col=c("tomato","skyblue"), pch=15, pt.cex=1.5, title="Species number >= 1 count")
dev.off()

### ordination
pca_original=prcomp(sqrt(sqrt(t(ordiorigdf))))
pca_resampled=prcomp(sqrt(sqrt(t(ordidf))))

png(paste0("output/2023-03-08_bacteria_clean_occ3_median_resampled_specieslevel_Sampleeffort",nsampleff,"_aggregated_pca_comparisonplot.png"), width=960,height=480)
par(mfrow=c(1,2))
biplot(pca_original, main="original data")
biplot(pca_resampled, main="rarefied data")
dev.off()

