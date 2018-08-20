#Initial scripts for taxonomic harmonization and phylogeny stuff

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FIA data prep

#For FIA, use "new scientific name" where not empty
FIA_taxa <- read.csv("C:/Users/Brian/Google Drive/DNH_scale/FIA_taxa.csv",na.strings = c("",NA))
FIA_taxa$ID<-1:nrow(FIA_taxa)

FIA_taxa$SPECIES<-as.character(FIA_taxa$SPECIES)
FIA_taxa$SCIENTIFIC_NAME<-as.character(FIA_taxa$SCIENTIFIC_NAME)
FIA_taxa$NEW_SCIENTIFIC_NAME<-as.character(FIA_taxa$NEW_SCIENTIFIC_NAME)
FIA_taxa$FAMILY<-as.character(FIA_taxa$FAMILY)
FIA_taxa$GENUS<-as.character(FIA_taxa$GENUS)

#populate missing species
FIA_taxa$SCIENTIFIC_NAME[which(is.na(FIA_taxa$SCIENTIFIC_NAME) & !is.na(FIA_taxa$GENUS) & !is.na(FIA_taxa$SPECIES))]<-
      paste(FIA_taxa$GENUS[which(is.na(FIA_taxa$SCIENTIFIC_NAME) & !is.na(FIA_taxa$GENUS) & !is.na(FIA_taxa$SPECIES))]," ",
            FIA_taxa$SPECIES[which(is.na(FIA_taxa$SCIENTIFIC_NAME) & !is.na(FIA_taxa$GENUS) & !is.na(FIA_taxa$SPECIES))],sep = "")

#populate missing genus and species
for(i in 1:nrow(FIA_taxa)){
  
if(is.na(FIA_taxa$GENUS[i]) & !is.na(FIA_taxa$SCIENTIFIC_NAME[i]) ){

FIA_taxa$GENUS[i]<-strsplit(FIA_taxa$SCIENTIFIC_NAME[i],split = " ")[[1]][1]

  
} #if scientific name populated, but not genus 
  
if(is.na(FIA_taxa$SPECIES[i]) & !is.na(FIA_taxa$SCIENTIFIC_NAME[i]) & length(unlist(strsplit(FIA_taxa$SCIENTIFIC_NAME[i],split = " ")))>=2
   ){
    
    FIA_taxa$SPECIES[i]<-strsplit(FIA_taxa$SCIENTIFIC_NAME[i],split = " ")[[1]][2]
    
    
  } #if scientific name populated, but not spp 
  
  
  
    
  
}
rm(i)


#populate missing families from FIA data

for(i in 1:nrow(FIA_taxa)){
if(is.na(FIA_taxa$FAMILY[i])){
new_fam<-unique(na.omit(FIA_taxa$FAMILY[which(FIA_taxa$GENUS==FIA_taxa$GENUS[i])]))


if(length(new_fam)!=1){

  
  bien_tax_i<-BIEN_taxonomy_genus(genus = as.character(FIA_taxa$GENUS[i]))  
  new_fam<-unique(bien_tax_i$scrubbed_family)
  new_fam<-new_fam[which(new_fam!="Unknown")]
  
}#try bien    


if(length(new_fam)==1){
FIA_taxa$FAMILY[i]<-new_fam  
}


if(length(new_fam)!=1){}


  
}#if family is missing
  
  
  
  
}
rm(i, new_fam,bien_tax_i)

#Manually correct spelling error in Apocynaceae
FIA_taxa$FAMILY[which(FIA_taxa$FAMILY=="Apocynacea")]<-"Apocynaceae"

#Prune to needed fields
FIA_taxa<-FIA_taxa[c('SCIENTIFIC_NAME',"NEW_SCIENTIFIC_NAME","FAMILY","GENUS","ID")]

#Replace any scientific name with the "new" scientific name
FIA_taxa$SCIENTIFIC_NAME[which(FIA_taxa$NEW_SCIENTIFIC_NAME!="")] <- 
    FIA_taxa$NEW_SCIENTIFIC_NAME[which(FIA_taxa$NEW_SCIENTIFIC_NAME!="")]

#Toss hybrids
FIA_taxa <- FIA_taxa[grep(pattern = " x ",x = FIA_taxa$SCIENTIFIC_NAME,invert = T),]

#Toss anything only ID'ed to the genus
FIA_taxa<-FIA_taxa[which(unlist(lapply(X = FIA_taxa$SCIENTIFIC_NAME,FUN = function(x){
  length(unlist(strsplit(x,split = " ")))>=2
}))),]


#Only keep genus and species
FIA_taxa$SCIENTIFIC_NAME<-unlist(lapply(X = FIA_taxa$SCIENTIFIC_NAME,FUN = function(x){
  paste(strsplit(x = x, split = " ")[[1]][1:2],collapse =" "  )
  
}))

FIA_taxa<-FIA_taxa[c("SCIENTIFIC_NAME","GENUS","FAMILY","ID")]
FIA_taxa<-unique(FIA_taxa)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#NEON data prep  

NEON_taxa <- read.csv("C:/Users/Brian/Google Drive/DNH_scale/NEON_taxa.csv")
NEON_taxa$scientificName <- as.character(NEON_taxa$scientificName)
#filter by taxon rank to only include species

NEON_taxa <- NEON_taxa[which(NEON_taxa$taxonRank == "species"),]

#Toss hybrids
NEON_taxa<-NEON_taxa[grep(pattern = "Ã",x = NEON_taxa$scientificName,invert = T),]

#Only keep genus and species
NEON_taxa$scientificName<-unlist(lapply(X = NEON_taxa$scientificName,FUN = function(x){
  paste(strsplit(x = x, split = " ")[[1]][1:2],collapse =" "  )
  
}))

NEON_taxa <- NEON_taxa[c("scientificName","family","taxonID")]
NEON_taxa <- unique(NEON_taxa)
NEON_taxa$genus<-unlist(lapply(X = NEON_taxa$scientificName,FUN = function(x){
  paste(strsplit(x = x, split = " ")[[1]][1],collapse =" "  )
  
}))

NEON_taxa<-NEON_taxa[c("scientificName","genus","family","taxonID")]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#USDA county lists
USDA_taxa<-read.csv(file = "C:/Users/Brian/Google Drive/DNH_scale/USDA/unique_USDA.csv")
USDA_taxa$species<-gsub(pattern = "_",replacement = " ",x = USDA_taxa$binomial)
USDA_taxa$ID<-paste("USDA",1:nrow(USDA_taxa),sep = "")
USDA_taxa<-USDA_taxa[c("species","Genus","Family","ID")]
#Toss the hybrids
USDA_taxa<-USDA_taxa[grep(pattern = "Ã",x = USDA_taxa$species,invert = T),]




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#combining sets of taxa

colnames(NEON_taxa)<-colnames(FIA_taxa)
colnames(USDA_taxa)<-colnames(FIA_taxa)
combined_taxa<- rbind(FIA_taxa,NEON_taxa,USDA_taxa)
rownames(combined_taxa)<-combined_taxa$ID
combined_taxa<-combined_taxa[-4]
combined_taxa<-unique(combined_taxa)
dup_taxa<-combined_taxa$SCIENTIFIC_NAME[which(duplicated(x = combined_taxa$SCIENTIFIC_NAME))]#Looks like a few duplicated sci names, many with wrong genus
dup_stuff<-combined_taxa[which(combined_taxa$SCIENTIFIC_NAME%in%dup_taxa),]

#Looks like species have been reassigned new genera, causing issues.

for(i in 1:nrow(combined_taxa)){

inferred_genus<-strsplit(combined_taxa$SCIENTIFIC_NAME[i],split = " ")[[1]][1]  
asserted_genus<-combined_taxa$GENUS[i]    
  
if(inferred_genus!=asserted_genus){
combined_taxa$GENUS[i]<-inferred_genus  
  
}# if the genus from the binomial and the genus field don't match, use the binomial  
  
}
combined_taxa<-unique(combined_taxa)


dup_taxa<-combined_taxa$SCIENTIFIC_NAME[which(duplicated(x = combined_taxa$SCIENTIFIC_NAME))]#Looks like a few duplicated sci names, many with wrong genus
dup_stuff<-combined_taxa[which(combined_taxa$SCIENTIFIC_NAME%in%dup_taxa),]

#Three species are assigned to different families.  Defer to BIEN

for(i in 1:length(unique(dup_stuff$SCIENTIFIC_NAME))){
  
sp_i<-unique(dup_stuff$SCIENTIFIC_NAME)[i]  
BIEN_tax<-BIEN_taxonomy_species(species = sp_i)  
fam_i<-unique(BIEN_tax$scrubbed_family)
if(length(fam_i)!=1){stop()}

combined_taxa$FAMILY[which(combined_taxa$SCIENTIFIC_NAME==sp_i) ]<-fam_i
  
}
combined_taxa<-unique(combined_taxa)


dup_taxa<-combined_taxa$SCIENTIFIC_NAME[which(duplicated(x = combined_taxa$SCIENTIFIC_NAME))]#Looks like a few duplicated sci names, many with wrong genus
dup_stuff<-combined_taxa[which(combined_taxa$SCIENTIFIC_NAME%in%dup_taxa),]

if(length(dup_taxa)!=0){stop("Duplicated shit, yo")}

rm(BIEN_tax,dup_stuff,dup_taxa,asserted_genus,fam_i,i,inferred_genus,sp_i)

combined_species<-combined_taxa$SCIENTIFIC_NAME


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Checking out phylogenies for taxon overlap
library(ape)


#BIEN tree first
BIEN_tree<-BIEN::BIEN_phylogeny_complete(n_phylogenies = 1)
paste(length(which(combined_species %in% gsub(pattern = "_",replacement = " ",x = BIEN_tree$tip.label)))/length(combined_species)*100, "percent coverage")
#about 63%

BIEN_cons_tree<-BIEN_phylogeny_conservative()
paste(length(which(combined_species %in% gsub(pattern = "_",replacement = " ",x = BIEN_cons_tree$tip.label)))/length(combined_species)*100, "percent coverage")
# 31%

#Smith 2017 all Open Tree backbone
allotb_tree <- read.tree("C:/Users/Brian/Desktop/phylogenies/Smith_big_seed_plant_trees_v0.1/ALLOTB.tre")
paste(length(which(combined_species %in% gsub(pattern = "_",replacement = " ",x = allotb_tree$tip.label)))/length(combined_species)*100, "percent coverage")
#about 76%

#Smith 2017 all Magallon backbone
allmb_tree <- read.tree("C:/Users/Brian/Desktop/phylogenies/Smith_big_seed_plant_trees_v0.1/ALLMB.tre")
paste(length(which(combined_species %in% gsub(pattern = "_",replacement = " ",x = allmb_tree$tip.label)))/length(combined_species)*100, "percent coverage")
#about 77%


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Looks like about 20% of species aren't matching any phylogeny
bad_names<-combined_species[which(!combined_species %in% gsub(pattern = "_",replacement = " ",x = allotb_tree$tip.label)  )]

library(BIEN)
source("r_scripts/sloppy_tnrs.R")
tnrsed_bad_names<-sloppy_tnrs(species_list = bad_names)

mismatched_names<-tnrsed_bad_names[which(tnrsed_bad_names$name_supplied != tnrsed_bad_names$scrubbed_species_binomial),]


#Zanne 2014

#qian hong has an updated version of Zanne in plotone      


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(ape)
library(phytools)
source("C:/Users/Brian/Desktop/current_projects/S.PhyloMaker/R_codes for S.PhyloMaker")
phytophylo<-read.tree(file = "C:/Users/Brian/Desktop/current_projects/S.PhyloMaker/PhytoPhylo")
nodes<-read.table("C:/Users/Brian/Desktop/current_projects/S.PhyloMaker/nodes",header=T,sep = "\t") # read in the nodes information of the megaphylogeny.
colnames(combined_taxa)<-c("species","genus","family")


example<-read.csv("C:/Users/Brian/Desktop/current_projects/S.PhyloMaker/example.splist",header=T) # read in the nodes information of the megaphylogeny.



#splist: a user specified data frame which contains a list of seed plant species for which S.PhyloMaker generates phylogenies. 
#This data frame should include three columns with field names of 
#species, genus, and family. 
#The three columns are, respectively, for species name (e.g. Acacia berlandieri), genus name (e.g. Acacia), and family name (e.g. Fabaceae).
combined_taxa$species<-gsub(pattern = " ",replacement = "_",x = combined_taxa$species)
results<-S.PhyloMaker(tree = phytophylo,spList = combined_taxa[1:3],nodes = nodes)

