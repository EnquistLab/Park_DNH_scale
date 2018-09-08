#Initial scripts for taxonomic harmonization and phylogeny stuff

library(BIEN)
library(ape)


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
NEON_taxa$taxonID<-paste(NEON_taxa$taxonID,1:nrow(NEON_taxa),sep = "")
NEON_taxa$scientificName <- as.character(NEON_taxa$scientificName)
#filter by taxon rank to only include species

NEON_taxa <- NEON_taxa[which(NEON_taxa$taxonRank == "species"),]

#Toss hybrids
NEON_taxa<-NEON_taxa[grep(pattern = "Ã",x = NEON_taxa$scientificName,invert = T),]

#Only keep genus and species
NEON_taxa$scientificName<-unlist(lapply(X = NEON_taxa$scientificName,FUN = function(x){
  paste(strsplit(x = x, split = " ")[[1]][1:2],collapse =" "  )
  
}))

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


#Remove taxa with X as species

combined_taxa <- combined_taxa[which(combined_taxa$GENUS != "x"),]
combined_taxa <- combined_taxa[which(combined_taxa$SCIENTIFIC_NAME != "(Valeton) Hosok."),]

rm(BIEN_tax,dup_stuff,dup_taxa,asserted_genus,fam_i,i,inferred_genus,sp_i)

rm(FIA_taxa,NEON_taxa,USDA_taxa)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Add anything Dan has that is currently missing
#annotate Dan's master list with species data

master_taxon<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/master_taxon_binomial.csv")
master_taxon<-master_taxon[c('binomial',"Genus","Family")]
colnames(master_taxon)<-colnames(combined_taxa)
master_taxon$SCIENTIFIC_NAME<-gsub(pattern = "_",replacement = " ",x = master_taxon$SCIENTIFIC_NAME)
combined_taxa<-rbind(combined_taxa,master_taxon)
combined_taxa<-unique(combined_taxa)
rm(master_taxon)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Checking out phylogenies for taxon overlap


#GB = genbank phylos

#Smith 2017 cons Open Tree backbone
gbotb_tree <- read.tree("C:/Users/Brian/Desktop/phylogenies/Smith_big_seed_plant_trees_v0.1/GBOTB.tre")
#Strip existing node labels
gbotb_tree$node.label[which(gbotb_tree$node.label!="")]<-""
summary(gbotb_tree)
paste(length(which(combined_taxa$SCIENTIFIC_NAME %in% gsub(pattern = "_",replacement = " ",x = gbotb_tree$tip.label)))/length(combined_taxa$SCIENTIFIC_NAME)*100, "percent coverage")
#about 45%
source("r_scripts/collapse_subsp.R")
gbotb_tree<-collapse_subsp(gbotb_tree)


#Double check that combined_taxa has no duplicated species

if(length(unique(combined_taxa$SCIENTIFIC_NAME))!=nrow(combined_taxa)){stop("Duplicated Taxa")}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Windows code for creating sunplin.spn object needed by sunplin code
  #cd "C:\Users\Brian\Desktop\current_projects\Park_DNH_scale"
  #R CMD SHLIB sunplin-r.cpp --output=sunplin.spn


source("C:/Users/Brian/Desktop/current_projects/misc_R_code/sunplin/sunplin/sunplin-functions.r")
dyn.load("C:/Users/Brian/Desktop/current_projects/misc_R_code/sunplin/sunplin/sunplin.spn")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prepping Smith tree:

#Need to add node labels where needed, and also identify which how species will be assigned to put (ie node based or tip based)

combined_taxa$put_level<-NA
combined_taxa$put<-NA

Smith_taxonomy<-BIEN_taxonomy_species(species = gsub(pattern = "_",replacement = " ",x = gbotb_tree$tip.label))

for( i in 1:nrow(combined_taxa)){
  
  #Check whether genus is present multiple times
  if(length(grep(pattern = paste(combined_taxa$GENUS[i],"_",sep = ""),x = gbotb_tree$tip.label))>1 ){
  combined_taxa$put_level[i]<-"mrca_genus"
  combined_taxa$put[i]<-combined_taxa$GENUS[i] }#if
  
  #Check whether genus is present at least once
  if(length(grep(pattern = paste(combined_taxa$GENUS[i],"_",sep = ""),x = gbotb_tree$tip.label))==1 ){
    combined_taxa$put_level[i]<-"congener"
    combined_taxa$put[i]<-gbotb_tree$tip.label[grep(pattern = paste(combined_taxa$GENUS[i],"_",sep = ""),x = gbotb_tree$tip.label)] }#if
  
  
  #Check whether family is present multiple times (and genus is not present)
  
    fam_i<-combined_taxa$FAMILY[i]
    
    if(length(fam_i)>=1 & is.na(combined_taxa$put_level[i])){
    
    if(length(fam_i)>1 | fam_i == "Unknown"){stop("too many families or only unknown")}  
    
    
    
    spp_in_fam_i<-gbotb_tree$tip.label[which(gbotb_tree$tip.label %in%   gsub(pattern = " ",replacement = "_",Smith_taxonomy$scrubbed_species_binomial[which(Smith_taxonomy$scrubbed_family==fam_i)]) )]
    
    
    #if genus attempts havent worked and there are multiple confamilials, use family mrca
    if( length(spp_in_fam_i)>1 &  is.na(combined_taxa$put_level[i])){
      combined_taxa$put_level[i]<-"mrca_family"
      combined_taxa$put[i]<-fam_i  }#if
    
      
    if( length(spp_in_fam_i)==1 &  is.na(combined_taxa$put_level[i])){
      combined_taxa$put_level[i]<-"confamilial"
      combined_taxa$put[i]<-spp_in_fam_i  }#if
    
    
    #If none of these have worked, the species will retain an NA, and be removed later
    
    rm(spp_in_fam_i,fam_i)
    }#only go through family stuff if we have a family in the taxonomy
      
}#for i loop
rm(fam_i,i)

#drop taxa that have NA for put level (mostly mosses and ferns, but a few weird seed plants (Rafflesia))
combined_taxa_og<-combined_taxa
combined_taxa <- combined_taxa[which(!is.na(combined_taxa$put_level)),]
write.csv(x = combined_taxa_og,file = "combined_taxa_BSM.csv",col.names = F)

#remove taxa that are already in the tree
combined_taxa<-combined_taxa[which(!gsub(pattern = " ",replacement = "_", x = combined_taxa$SCIENTIFIC_NAME) %in% gbotb_tree$tip.label),]

#2) Need to label all clades that will have species inserted into them (ie genera or families if that fails)


genera_to_add<-unique(combined_taxa$GENUS[which(combined_taxa$put_level=="mrca_genus")])
families_to_add<-unique(combined_taxa$FAMILY[which(combined_taxa$put_level=="mrca_family")])

#Label genera nodes
for(i in 1:length(genera_to_add)){
  
genus_i<-genera_to_add[i]  

mrca_i<-getMRCA(phy = gbotb_tree,
        tip = grep(pattern = paste(genus_i,"_",sep = ""),x = gbotb_tree$tip.label)   )

#if the node isnt labelled yet, label it.
if(gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]==""){

gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]<-genus_i

}else{

#if the node IS labelled, modify the put of the genus in question to match the node label
if(gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]!=""){
  
  label_i<-gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]
  combined_taxa$put[which(combined_taxa$put_level=="mrca_genus" & combined_taxa$put==genus_i)]<-label_i
  rm(label_i)
    
}else{stop("Something weird happened")}

}#if the node is labelled

rm(genus_i,mrca_i)

}


#Label family nodes

for(i in 1:length(families_to_add)){
  
  fam_i<-families_to_add[i]  
  spp_in_family <- Smith_taxonomy$scrubbed_species_binomial[which(Smith_taxonomy$scrubbed_family == fam_i)]
  
  
  mrca_i<-getMRCA(phy = gbotb_tree,
                  tip = which(gbotb_tree$tip.label %in% gsub(pattern = " ",replacement = "_", x = spp_in_family   ) )) 
  

  
  #if the node isnt labelled yet, label it.
  if(gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]==""){
    
    gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]<-fam_i
    
  }else{
    
    #if the node IS labelled, modify the put of the genus in question to match the node label
    if(gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]!=""){
      
      label_i<-gbotb_tree$node.label[mrca_i-length(gbotb_tree$tip.label)]
      combined_taxa$put[which(combined_taxa$put_level=="mrca_family" & combined_taxa$put==fam_i)]<-label_i
      rm(label_i)
      
    }else{stop("Something weird happened")}
    
  }#if the node is labelled
  
  
  rm(fam_i,mrca_i,spp_in_family)
  
}
rm(i,genera_to_add,families_to_add)



#3) Need to prepare a .puts (phylogenetically uncertain taxa) file containing all species to be added

file.create("gbotb_tree.puts")
writeLines(text = paste(gsub(pattern = " ",replacement = "_",x = combined_taxa$SCIENTIFIC_NAME)," ",combined_taxa$put,sep = ""),
           con = "gbotb_tree.puts")  



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#write current version of tree (with node labels)

write.tree(phy = gbotb_tree,file = "Smith_2017_gbotb.tre")
nreps=1000
for(i in 1:nreps){
    
  trees <- sunplin.expd("Smith_2017_gbotb.tre","gbotb_tree.puts",numTree =  1,method = 2)  
  t1<-strsplit(x = trees,split = ";")[[1]][2]
  t2<-read.tree(text = paste(t1,";",sep = " "))
  t2<-drop.tip(phy = t2,
           tip = t2$tip.label[which(!t2$tip.label %in% gsub(pattern = " ",replacement = "_",x = combined_taxa_og$SCIENTIFIC_NAME))] ,
           trim.internal = T, 
           collapse.singles = T
               
           )
  
  
  
  write.tree(phy = t2,file = paste("sunplin_trees/gbotb_tree_sunplin_method_2_rep",i,".tre",sep = ""))  
  print(paste(i/nreps*100, "percent done"))
  
  
}
rm(nreps)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

