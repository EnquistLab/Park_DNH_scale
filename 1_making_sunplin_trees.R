#Making sunplin trees v2
#install.packages("remotes")
#remotes::install_github("davidnipperess/PDcalc")

library(ape)
source("r_scripts/get_genera.R")
source("sunplin-functions.r")
source("r_scripts/sunplin_fxs_node_labels.R")
source("r_scripts/collapse_subsp.R")

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
write.csv(x = taxa,file = "L48_taxa.csv",row.names = F)

#allmb performance
  #species
  length(intersect(allmb$tip.label,taxa$binomial))/nrow(taxa) #13986, 80%
  (length(allmb$tip.label)-1)-allmb$Nnode #271,897 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(allmb)))/length(unique(taxa$Genus)) #95%
  #unique tip labels
  length(unique(allmb$tip.label))/length(allmb$tip.label)
  
  
#allotb performance
  #species
  length(intersect(allotb$tip.label,taxa$binomial))/nrow(taxa) #13903, 79% 
  (length(allotb$tip.label)-1)-allotb$Nnode #267,505 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(allotb)))/length(unique(taxa$Genus)) #95%
  #unique tip labels
  length(unique(allotb$tip.label))/length(allotb$tip.label)
  
    
#gbmb performance
  #species
  length(intersect(gbmb$tip.label,taxa$binomial))/nrow(taxa) #8398, 48%
  (length(gbmb$tip.label)-1)-gbmb$Nnode #60 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(gbmb)))/length(unique(taxa$Genus))#88%
  #unique tip labels
  length(unique(gbmb$tip.label))/length(gbmb$tip.label)
  
  
#gbotb performance
  #species
  length(intersect(gbotb$tip.label,taxa$binomial))/nrow(taxa) #8398, 48%
  (length(gbotb$tip.label)-1)-gbotb$Nnode #62 spp in polytomies
  #genera
  length(which(unique(taxa$Genus)  %in% get_genera(gbotb)))/length(unique(taxa$Genus)) #88%
  #unique tip labels
  length(unique(gbotb$tip.label))/length(gbotb$tip.label)
  

########################################################################################################
  
#Step 2
  
#convert all tips to species-level  

#allmb
allmb<- collapse_subsp(phylogeny = allmb)
  #unique tip labels
  length(unique(allmb$tip.label))/length(allmb$tip.label)
  
#allotb
allotb <- collapse_subsp(phylogeny = allotb)  
  #unique tip labels
  length(unique(allotb$tip.label))/length(allotb$tip.label)

#gbmb  
gbmb <-collapse_subsp(phylogeny = gbmb)
  #unique tip labels
  length(unique(gbmb$tip.label))/length(gbmb$tip.label)


#gbotb
gbotb <- collapse_subsp(phylogeny = gbotb)
  #unique tip labels
  length(unique(gbotb$tip.label))/length(gbotb$tip.label)




#########################################################################################################


#Step 3

#fix synonyms

#allmb

allmb <- check_synonyms(phylogeny = allmb, species_list = taxa$binomial)
  #unique tip labels
  length(unique(allmb$tip.label))/length(allmb$tip.label)

#allotb

allotb <- check_synonyms(phylogeny = allotb, species_list = taxa$binomial)
  #unique tip labels
  length(unique(allotb$tip.label))/length(allotb$tip.label)

#gbmb  

gbmb <- check_synonyms(phylogeny = gbmb, species_list = taxa$binomial)
  #unique tip labels
  length(unique(gbmb$tip.label))/length(gbmb$tip.label)

#gbotb

gbotb <- check_synonyms(phylogeny = gbotb, species_list = taxa$binomial)
  #unique tip labels
  length(unique(gbotb$tip.label))/length(gbotb$tip.label)


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

#allotb
#allmb

#gbotb
gbotb_genera <- get_genera(phylogeny = gbotb)
gbotb_genera_to_add <- setdiff(taxa$Genus,gbotb_genera)

#gbmb
gbmb_genera <- get_genera(phylogeny = gbmb)
gbmb_genera_to_add <- setdiff(taxa$Genus,gbmb_genera)

#which families do the missing genera fall in?  I'll pull full lists of genera from Kew for these fams
unique(taxa$Family[which(taxa$Genus %in% union(gbmb_genera_to_add,gbotb_genera_to_add))])
#write.csv(x = unique(taxa$Family[which(taxa$Genus %in% union(gbmb_genera_to_add,gbotb_genera_to_add))]),file = "fams_needed.csv")

#Dan has requested we follow POTW for placing missing genera.   I've combed through POTW online to get all the families that represent the missing genera
families_potw<-read.csv("families_to_add.csv",stringsAsFactors = F)
for(i in 1:length(families_potw$genus)){families_potw$genus[i]<-strsplit(x = families_potw$genus[i],split = " ")[[1]][1]}
families_potw <- families_potw[which(families_potw$genus!="x"),]

#missing anything?
missing_genera<-setdiff(x = union(gbmb_genera_to_add,gbotb_genera_to_add),y = families_potw$genus)
setdiff(x = taxa$Family[which(taxa$Genus %in% missing_genera)],y = families_potw$ï..Family)
#so missing a few genera due to taxonomic revisions (e.g. renamed species).  I'll graft these on with the original names.
colnames(families_potw)<-c("Family","Genus")
families_potw <- rbind(families_potw,taxa[c("Family","Genus")][which(taxa$Genus %in% missing_genera),])

#4.1 prepare puts info files
#allmb_puts <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = allmb,genus_only_addition = F)
#length(which(allmb_puts$put_level!="remove"))/nrow(allmb_puts_genus_only) #99 percent

#allotb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = allotb,genus_only_addition = T)
#length(which(allotb_puts_genus_only$put_level!="remove"))/nrow(allotb_puts_genus_only) #99 percent

gbotb_puts_genus_and_family <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],
                                                        phylogeny = gbotb, 
                                                        genus_only_addition = F,
                                                  genera_in_grafting_fams = families_potw)

length(which(gbotb_puts_genus_and_family$put_level!="remove"))/nrow(gbotb_puts_genus_and_family) #99.99 percent

gbmb_puts_genus_and_family <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],
                                                       phylogeny = gbmb,
                                                       genus_only_addition = F,
                                                       genera_in_grafting_fams = families_potw)

length(which(gbmb_puts_genus_and_family$put_level!="remove"))/nrow(gbmb_puts_genus_and_family) # 99.99 percent

#5.2 prepare puts phylo, .puts files

#make_puts_input_node_labels(puts_info = allmb_puts_genus_only,
#                            phylogeny = allmb,
#                            phylogeny_filename = "allmb_genus_only_puts_phylo.tre",
#                            puts_filename = "allmb_genus_only.puts")

#make_puts_input_node_labels(puts_info = allotb_puts_genus_only,
#                            phylogeny = allotb,
#                            phylogeny_filename = "allotb_genus_only_puts_phylo.tre",
#                            puts_filename = "allotb_genus_only.puts")

make_puts_input_node_labels(puts_info = gbmb_puts_genus_and_family,
                            phylogeny = gbmb,
                            phylogeny_filename = "gbmb_genus_and_family_phylo.tre",
                            puts_filename = "gbmb_genus_and_family.puts",
                            genera_in_grafting_fams = families_potw )

make_puts_input_node_labels(puts_info = gbotb_puts_genus_only,
                            phylogeny = gbotb,
                            phylogeny_filename = "gbotb_genus_and_family_phylo.tre",
                            puts_filename = "gbotb_genus_and_family.puts",
                            genera_in_grafting_fams = families_potw)

#4.3 make replicated phylogenies using puts info

sunplin_phylo_replicates(put_file = "gbmb_genus_and_family.puts",
                         phylogeny_file = "gbmb_genus_and_family_phylo.tre",
                         output_directory = "sunplin_trees/genus_and_family_additions/gbmb_genus_and_family_additions/",
                         output_base_filename = "gbmb_and_family_additions",
                         nrep = 1000)

sunplin_phylo_replicates(put_file = "gbotb_genus_and_family.puts",
                         phylogeny_file = "gbotb_genus_and_family_phylo.tre",
                         output_directory = "sunplin_trees/genus_and_family_additions/gbotb_genus_and_family_additions/",
                         output_base_filename = "gbotb_genus_and_family_additions",
                         nrep = 1000)
  
#########################################################################################################

#Step 5

#sunplin trees w/ genera additions (no family-level additions)

#5.1 prepare puts info files
allmb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = allmb,genus_only_addition = T)
length(which(allmb_puts_genus_only$put_level!="remove"))/nrow(allmb_puts_genus_only) #99 percent

allotb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = allotb,genus_only_addition = T)
length(which(allotb_puts_genus_only$put_level!="remove"))/nrow(allotb_puts_genus_only) #99 percent

gbotb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = gbotb,genus_only_addition = T)
length(which(gbotb_puts_genus_only$put_level!="remove"))/nrow(gbotb_puts_genus_only) #94 percent

gbmb_puts_genus_only <- get_put_info_node_labels(sp_fam = taxa[c("binomial","Family")],phylogeny = gbmb,genus_only_addition = T)
length(which(gbmb_puts_genus_only$put_level!="remove"))/nrow(gbmb_puts_genus_only) # 94 percent

#5.2 prepare puts phylo, .puts files

make_puts_input_node_labels(puts_info = allmb_puts_genus_only,
                            phylogeny = allmb,
                            phylogeny_filename = "allmb_genus_only_puts_phylo.tre",
                            puts_filename = "allmb_genus_only.puts")

make_puts_input_node_labels(puts_info = allotb_puts_genus_only,
                            phylogeny = allotb,
                            phylogeny_filename = "allotb_genus_only_puts_phylo.tre",
                            puts_filename = "allotb_genus_only.puts")

make_puts_input_node_labels(puts_info = gbmb_puts_genus_only,
                            phylogeny = gbmb,
                            phylogeny_filename = "gbmb_genus_only_puts_phylo.tre",
                            puts_filename = "gbmb_genus_only.puts")

make_puts_input_node_labels(puts_info = gbotb_puts_genus_only,
                            phylogeny = gbotb,
                            phylogeny_filename = "gbotb_genus_only_puts_phylo.tre",
                            puts_filename = "gbotb_genus_only.puts")


#5.2 make replicated phylogenies using puts info

#temp <- read.tree(file = "allmb_genus_only_puts_phylo.tre")
#temp <- multi2di(phy = temp)
#write.tree(phy = temp,file = "temp.tre")


#Crashes V
#sunplin_phylo_replicates(put_file = "allmb_genus_only.puts",
#                         phylogeny_file = "temp.tre",
#                         output_directory = "sunplin_trees/genus_addition_only/allmb_genus_additions_only/",
#                         output_base_filename = "allmb_genus_addition_only",
#                         nrep = 1000)


#temp <- drop.tip(phy = temp,tip = setdiff(x = temp$tip.label,y = taxa$binomial))

#sunplin_phylo_replicates(put_file = "allmb_genus_only.puts",
#                         phylogeny_file = "temp.tre",
#                         output_directory = "sunplin_trees/genus_addition_only/allmb_genus_additions_only/",
#                         output_base_filename = "allmb_genus_addition_only",
#                         nrep = 1000)



##############################

sunplin_phylo_replicates(put_file = "allmb_genus_only.puts",
                         phylogeny_file = "allmb_genus_only_puts_phylo.tre",
                         output_directory = "sunplin_trees/genus_addition_only/allmb_genus_additions_only/",
                         output_base_filename = "allmb_genus_addition_only",
                         nrep = 1000)

sunplin_phylo_replicates(put_file = "allotb_genus_only.puts",
                         phylogeny_file = "allotb_genus_only_puts_phylo.tre",
                         output_directory = "sunplin_trees/genus_addition_only/allotb_genus_additions_only/",
                         output_base_filename = "allotb_genus_addition_only",
                         nrep = 1000)

#sunplin_phylo_replicates(put_file = "gbmb_genus_only.puts",
#                         phylogeny_file = "gbmb_genus_only_puts_phylo.tre",
#                         output_directory = "sunplin_trees/genus_addition_only/gbmb_genus_additions_only/",
#                         output_base_filename = "gbmb_genus_addition_only",
#                         nrep = 1000)

#sunplin_phylo_replicates(put_file = "gbotb_genus_only.puts",
#                         phylogeny_file = "gbotb_genus_only_puts_phylo.tre",
#                         output_directory = "sunplin_trees/genus_addition_only/gbotb_genus_additions_only/",
#                         output_base_filename = "gbotb_genus_addition_only",
#                         nrep = 1000)





#####################################

#Copy files to google drive (since Github is running crazy slow..well, internet is crazy slow)

file.copy(from = list.files("sunplin_trees/genus_addition_only/gbmb_genus_additions_only/",full.names = T),
          to = "C:/Users/Brian/Google Drive/Park_DNH_trees/sunplin_trees/genus_addition_only/gbmb_genus_additions_only/")

file.copy(from = list.files("sunplin_trees/genus_addition_only/gbotb_genus_additions_only/",full.names = T),
          to = "C:/Users/Brian/Google Drive/Park_DNH_trees/sunplin_trees/genus_addition_only/gbotb_genus_additions_only/")

file.copy(from = list.files("sunplin_trees/genus_and_family_additions/gbmb_genus_and_family_additions/",full.names = T),
          to = "C:/Users/Brian/Google Drive/Park_DNH_trees/sunplin_trees/genus_and_family_additions/gbmb_genus_and_family_additions/")


file.copy(from = list.files("sunplin_trees/genus_and_family_additions/gbotb_genus_and_family_additions/",full.names = T),
          to = "C:/Users/Brian/Google Drive/Park_DNH_trees/sunplin_trees/genus_and_family_additions/gbotb_genus_and_family_additions/")
