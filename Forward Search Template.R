#Script developed by Tyler, all functions are within linkagemapping
#load the linkagemapping package

devtools::install_github("AndersenLab/linkagemapping")
library("linkagemapping")

#insert cross data
data("N2xCB4856cross")
cross <- N2xCB4856cross


#insert pheno data
pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/ForTrip/RIAILs2_processed.rds")

#map from subset of phenotype or load map object if mapping has been saved
## only mapping bleomycin from the RIAILs2 data set (reprocessed)
#map <- fsearch(cross, iterations = 1000, phenotype="bleomycin")

#save(map, file="~/Dropbox/AndersenLab/LabFolders/Shannon/Scripts/Linkage Mapping/forwardsearchbleo.Rda")

#this mapping was already done, loading the mapping info now under name "map"
load("~/Dropbox/AndersenLab/LabFolders/Shannon/Scripts/Linkage Mapping/forwardsearchbleo.Rda")

#create completed cross object with pheno data set
cross <- mergepheno(cross, pheno, set = 2)

#annotate lods to find CI, VE, ES
annotatedlods <- annotate_lods(map, cross)

#subset your trait of interest
mean.TOF <- filter(annotatedlods, trait=="bleomycin.mean.TOF")

#plot the max lods
lodplot(mean.TOF)

#plot the effect size
effectplot(cross, mean.TOF)

#plot the phenotype by genotype split
pxgplot(cross, mean.TOF)
