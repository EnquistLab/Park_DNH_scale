#Making sunplin trees v2

library(ape)
source("r_scripts/get_genera.R")


#1 Check Smith and Brown's multiple trees for species and genus coverage and polytomies. 
#Pick one that has the best coverage and least polytomies. 
#Check how SUNPLIN deals with polytomies. 
#OTB series have better resolution for taxa without molecular data, 
#MB has less conflicts (more monophyly) but more uncertainty for grafted groups.

#For species names that are not matched, scrub the Smith tree tips with TNRS for synonyms and see if there are any more matches. 
#INCLUDE SUBSPECIES AND VARIETIES. 
#Spot check a few names to ensure the TNRS is not doing something weird again.

#Check for existing family nodes. Get MRCA of genera by grepping genus names (beware case and overlap with family names - e.g., Aster, Asteraceae, Cotoneaster).
#Generate 1) 1000 trees with taxa glued on to genus and family nodes; 2)1000 trees with taxa glued on to genus nodes only; and 3) the raw scrubbed, pruned base phylogeny.
#Use genus nodes whenever possible. Use pre-existing family nodes if they exist. If not, information on the species/genera in each family can be found at: http://www.theplantlist.org/1.1/browse/A/
  
  
#####################################################################################
#Step 1:

# Check Smith and Brown's multiple trees for species and genus coverage and polytomies. 
# Pick one that has the best coverage and least polytomies. Check how SUNPLIN deals with polytomies. 
# OTB series have better resolution for taxa without molecular data, MB has less conflicts (more monophyly) but more uncertainty for grafted groups.

#load phylos:
allmb <- read.tree("smith_and_brown_2018_trees/v0.1/v0.1/ALLMB.tre")
allotb <- read.tree("smith_and_brown_2018_trees/v0.1/v0.1/ALLOTB.tre")
gbotb <- read.tree("smith_and_brown_2018_trees/v0.1/v0.1/GBOTB.tre")
gbmb <- read.tree("smith_and_brown_2018_trees/v0.1/v0.1/GBMB.tre")

#load taxa
taxa<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/L48_taxa.csv",stringsAsFactors = F)

#allmb performance
  #species
  length(intersect(allmb$tip.label,taxa$binomial))/nrow(taxa) #13986, 80%
  (length(allmb$tip.label)-1)-allmb$Nnode #271,897 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(allmb)))/length(unique(taxa$Genus)) #95%
  
  
#allotb performance
  #species
  length(intersect(allotb$tip.label,taxa$binomial))/nrow(taxa) #13903, 79% 
  (length(allotb$tip.label)-1)-allotb$Nnode #267,505 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(allotb)))/length(unique(taxa$Genus)) #95%
  
#gbmb performance
  #species
  length(intersect(gbmb$tip.label,taxa$binomial))/nrow(taxa) #8398, 48%
  (length(gbmb$tip.label)-1)-gbmb$Nnode #60 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(gbmb)))/length(unique(taxa$Genus))#88%
  
#gbotb performance
  #species
  length(intersect(gbotb$tip.label,taxa$binomial))/nrow(taxa) #8398, 48%
  (length(gbotb$tip.label)-1)-gbotb$Nnode #62 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(gbotb)))/length(unique(taxa$Genus)) #88%
  

########################################################################################################
  
#Step 2
  
#convert all tips to species-level  

#allmb
allmb <- make_species_level(phylogeny = allmb)  

#allotb
allotb <- make_species_level(phylogeny = allotb)  

#gbmb  
gbmb <- make_species_level(phylogeny = gbmb)

#gbotb
gbotb <- make_species_level(phylogeny = gbotb)




#########################################################################################################


#Step 3

#fix synonyms

#allmb

allmb <- check_synonyms(phylogeny = allmb, species_list = taxa$binomial)

#allotb

allotb <- check_synonyms(phylogeny = allotb, species_list = taxa$binomial)

#gbmb  

gbmb <- check_synonyms(phylogeny = gbmb, species_list = taxa$binomial)

#gbotb

gbotb <- check_synonyms(phylogeny = gbotb, species_list = taxa$binomial)


#save modified phylogeny
write.tree(phy = allmb,"trees_w_updated_names/allmb.tre")
write.tree(phy = allotb,"trees_w_updated_names/allotb.tre")
write.tree(phy = gbmb,"trees_w_updated_names/gbmb.tre")
write.tree(phy = gbotb,"trees_w_updated_names/gbotb.tre")

#########################################################################################################

#Step 4

#sunplin trees w/ family and genera additions

#make sp_fam dataframe needed for sunplin fx
sp_fam <- taxa[c("binomial","Family")]


allmb$node.label[grep(pattern = "eae",x = allmb$node.label)]#19
allotb$node.label[grep(pattern = "eae",x = allotb$node.label)]#13

gbmb$node.label[grep(pattern = "eae",x = gbmb$node.label)]#3
gbotb$node.label[grep(pattern = "eae",x = gbotb$node.label)]#3


#########################################################################################################

#Step 5

#sunplin trees w/ genera additions (no family-level additions)

allmb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = allmb,genus_only_addition = T)
length(which(allmb_puts_genus_only$put_level!="remove"))/nrow(allmb_puts_genus_only) #99 percent

allotb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = allotb,genus_only_addition = T)
length(which(allotb_puts_genus_only$put_level!="remove"))/nrow(allotb_puts_genus_only) #99 percent

gbotb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = gbotb,genus_only_addition = T)
length(which(gbotb_puts_genus_only$put_level!="remove"))/nrow(gbotb_puts_genus_only) #94 percent

gbmb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = gbmb,genus_only_addition = T)
length(which(gbmb_puts_genus_only$put_level!="remove"))/nrow(gbmb_puts_genus_only) # percent



