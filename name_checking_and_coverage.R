library(ape)
library(phytools)
gbotb<-read.tree("C:/Users/Brian/Desktop/current_projects/Park_DNH_scale/Smith_2017_gbotb.tre")
taxa<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/L48_taxa.csv")

length(unique(taxa$binomial)) #initial dataset for L48 had 17,649 species
length(setdiff(x = taxa$binomial,y = gbotb$tip.label))#8958 species needed to be added, the remainder were present


#Of the 17,649 species in our dataset, approximately half (8691) were already present in the phylogeny
8691/17649

#how many genera?

tree_genera<-NULL
for(i in gbotb$tip.label){
  
  tree_genera<-unique(c(tree_genera,strsplit(x = i,split = "_")[[1]][1]))

}

taxa$Genus
length(unique(intersect(as.character(taxa$Genus),tree_genera)))/length(unique(taxa$Genus))

#~10 percent of genera missing

setdiff(as.character(taxa$Genus),tree_genera)
length(which(taxa$Genus %in% intersect(as.character(taxa$Genus),tree_genera)))/nrow(taxa) #6 pct of species

spp_to_add <- setdiff(x = gsub(pattern = " ",replacement = "_",x = taxa$binomial),y = gbotb$tip.label)

taxhead<-BIEN:::.BIEN_sql("SELECT * FROM bien_taxonomy LIMIT 1;")


sloppy_tnrs<-function(species, ...){
  
  query = paste("SELECT * FROM bien_taxonomy 
              WHERE scrubbed_species_binomial IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ")  ; ")  
  
  
  return(BIEN:::.BIEN_sql(query, ...))
  
  
  
}

tnrs_out <- sloppy_tnrs(species = gsub(pattern = "_",replacement = " ",x = spp_to_add))

not_scrubbed <- setdiff(x = gsub(pattern = "_",replacement = " ",x = spp_to_add),y = unique(tnrs_out$scrubbed_species_binomial))

# only 84 species were not matched as a "scrubbed name", predominantly accepted names but some "no opinions"

sloppy_tnrs<-function(species, ...){
  
  query = paste("SELECT * FROM bien_taxonomy 
              WHERE scrubbed_species_binomial IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ")  ; ")  
  
  
  return(BIEN:::.BIEN_sql(query, ...))
  
  
  
}


sloppy_tnrs_pt2 <- function(species, ...){
  
  query = paste("SELECT scrubbed_species_binomial, name_matched, verbatim_scientific_name 
                FROM view_full_occurrence_individual 
                WHERE name_matched IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ") 
                OR verbatim_scientific_name IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ") ; ")  
  
  
  return(BIEN:::.BIEN_sql(query, ...))
  
  
  
}



tnrs_out2 <- sloppy_tnrs_pt2(species = not_scrubbed)
tnrs_out2<-unique(tnrs_out2)


tnrs_out2$scrubbed_species_binomial

length(which(as.character(sapply(X = tnrs_out2$scrubbed_species_binomial,FUN = function(x){strsplit(x,split = " ")[[1]][1]}))==
               as.character(sapply(X = tnrs_out2$name_matched, FUN = function(x){strsplit(x,split = " ")[[1]][1]}))))

tnrs_3<-tnrs_out2[which(as.character(sapply(X = tnrs_out2$scrubbed_species_binomial,FUN = function(x){strsplit(x,split = " ")[[1]][1]}))!=
                          as.character(sapply(X = tnrs_out2$name_matched, FUN = function(x){strsplit(x,split = " ")[[1]][1]}))),]


unique(tnrs_out2$scrubbed_species_binomial[which(as.character(sapply(X = tnrs_out2$scrubbed_species_binomial,FUN = function(x){strsplit(x,split = " ")[[1]][1]}))==
                                                   as.character(sapply(X = tnrs_out2$name_matched, FUN = function(x){strsplit(x,split = " ")[[1]][1]})))])

#10 spelling variants


unique(tnrs_3$scrubbed_species_binomial)
#22 species that have been re-assigned


not_scrubbed[which(!not_scrubbed %in% c(tnrs_out2$name_matched,tnrs_out2$verbatim_scientific_name))]
#49 species not matched at all
length(unique(tnrs_out$scrubbed_species_binomial))/length(unique(spp_to_add))*100


#How many species overall match to bien taxonomy?
alltnrs<-sloppy_tnrs(species = gsub(pattern = "_",replacement = " ",x = taxa$binomial))

length(which(!gsub(pattern = "_",replacement = " ",x = taxa$binomial) %in%alltnrs$scrubbed_species_binomial))/length(unique(alltnrs$scrubbed_species_binomial))*100
#less than one half of one percent of names do not match to TNRS scrubbed names

#####################################################
#####################################################
#figure
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

set.seed(2005)
random_tree <- read.tree(list.files(path = "C:/Users/Brian/Desktop/current_projects/DNH_Park/phylos/sunplin_trees/",pattern = ".tre",full.names = T)[runif(n = 1,min = 1,max = 1000)])
plot(random_tree,show.tip.label = F)



library(ggtree)
library(ggplot2)
library(cowplot)

status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/status4.csv",stringsAsFactors = F)

pruned_tree<-drop.tip(phy = random_tree,tip = setdiff(x = random_tree$tip.label,y = status$binomial))

angios <- getMRCA(phy = pruned_tree, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = pruned_tree$tip.label) )
gymnos <- getMRCA(phy = pruned_tree, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = pruned_tree$tip.label) )

#phylo, invasions, angios, scale bar
ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = 5) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+theme(legend.position="right")+
  geom_cladelabel(node = angios, label="", offset = 20,barsize = 3,color="grey")



#We should have a scale bar for the branch lengths.
#A legend might be nice (red = non-native, blue = native).
#Maybe make the colored points a bit bigger (or depict them as bars).
#If it's not too much trouble, maybe label some lower level groups to show we care?
#Save as pdf (or other vector graphic) so that people can zoom in.



l1 <- get_legend(ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
plot(l1)

t1<-ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = angios, label="Angiosperms", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")


plot(t1)

combined<- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)


ggsave(filename = "DNH_scale_phylo",plot = combined,device = "svg",path = "C:/Users/Brian/Desktop/",width = 20,height = 20,units = "cm")

bien_taxonomy_status<-BIEN::BIEN_taxonomy_species(species = gsub(pattern = "_",replacement = " ",x = status$binomial))
unique(bien_taxonomy_status$order)

orders<-read.table("list_of_orders.txt")

status<-merge(x = status,y = orders,by.x = "Family",by.y = "V2")
colnames(x = status)[length(status)]<-"Order"

bien_taxonomy_status















#Modify above to add most common families
library(BIEN)
bien_taxonomy <- BIEN_taxonomy_species(species = gsub(pattern = "_",replacement = " ",x = pruned_tree$tip.label))
which(bien_taxonomy$scrubbed_species_binomial %in% gsub(pattern = "_",replacement = " ",x = pruned_tree$tip.label))
bien_taxonomy<-unique(bien_taxonomy[c("scrubbed_family","scrubbed_species_binomial")])
  #are any species listed in different families?
dups<- bien_taxonomy$scrubbed_species_binomial[which(duplicated(bien_taxonomy$scrubbed_species_binomial))]
#drop these conflicting species from the taxonomy

bien_dups <- bien_taxonomy[which(bien_taxonomy$scrubbed_species_binomial %in% dups),]
bien_dups<-bien_dups[c('scrubbed_species_binomial','scrubbed_family')]
bien_dups<-unique(bien_dups)
bien_taxonomy <- bien_taxonomy[which(!bien_taxonomy$scrubbed_species_binomial %in% dups),]

fam_table<- sort(table(bien_taxonomy$scrubbed_family),decreasing = T)
names(fam_table[1:3])

length(intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Asteraceae")]) 
          ,y = pruned_tree$tip.label))

length(pruned_tree$tip.label)


aster <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Asteraceae")]) 
                                                      ,y = pruned_tree$tip.label) )


fab <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Fabaceae")]) 
                                                    ,y = pruned_tree$tip.label) )


poa <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Poaceae")]) 
                                                  ,y = pruned_tree$tip.label) )



#Asteraceae, Fabaceae, Fagaceae, Rosaceae, Poaceae, Orchidaceae, Apiaceae, Lamiaceae, Brassicaeae, etc.



  ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = angios, label="Angiosperms", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = aster, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")


poa <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Poaceae")]) 
                                                    ,y = pruned_tree$tip.label) )
  


################################


gbotb_tnrs_out<-sloppy_tnrs(species= gsub(pattern = "_",replacement = " ",x = gbotb$tip.label))
unique(gbotb_tnrs_out$scrubbed_taxonomic_status)
length(unique(gbotb_tnrs_out$scrubbed_species_binomial))
length(gbotb$tip.label)
gbotb_unmatched<-setdiff(gsub(pattern = "_",replacement = " ",x = gbotb$tip.label),unique(gbotb_tnrs_out$scrubbed_species_binomial))


gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched)

length(which(sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2))/length(gbotb_unmatched)
#33% sub species, varieties, etc.

length(grep(pattern = "\\.",x = gbotb_unmatched))/length(gbotb_unmatched)
#35% .sp, etc.

#47% subspecies, other obviously bad names


gbotb_unmatched_species <- gbotb_unmatched[which(!sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2 & !grepl(pattern = "\\.",x = gbotb_unmatched))]

gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched_species)

#Do any non-matched taxa correspond to matched taxa?


missing_spp<- setdiff(x = taxa$binomial,y = gbotb$tip.label)

missing_spp_tnrs_out<- sloppy_tnrs(species = gsub(pattern = "_",replacement = " ",x = missing_spp))
#how many were matched?
length(unique(missing_spp_tnrs_out$scrubbed_species_binomial))#8874/8958 # all good names

phylo_spp<- setdiff(y = as.character(taxa$binomial),x = gbotb$tip.label)
length(grep(pattern = "\\.",x = gbotb_unmatched))/length(gbotb_unmatched)\

phylo_spp <-  phylo_spp[which(!sapply(X = phylo_spp,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2 & !grepl(pattern = "\\.",x = phylo_spp))]


#how many species from the phylo (that didn't match)





gbotb_tnrs_out<-sloppy_tnrs(species= gsub(pattern = "_",replacement = " ",x = gbotb$tip.label))
unique(gbotb_tnrs_out$scrubbed_taxonomic_status)
length(unique(gbotb_tnrs_out$scrubbed_species_binomial))#69960 species in the phylogeny match to TNRS scrubbed names
length(unique(gbotb_tnrs_out$scrubbed_species_binomial))/length(gbotb$tip.label) #89%

gbotb_unmatched<-setdiff(gsub(pattern = "_",replacement = " ",x = gbotb$tip.label),unique(gbotb_tnrs_out$scrubbed_species_binomial))
#of the remaining ~11%:

#gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched)

length(which(sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2))/length(gbotb_unmatched)
#33% sub species, varieties, etc.

length(grep(pattern = "\\.",x = gbotb_unmatched))/length(gbotb_unmatched)
#35% .sp, etc.

#37% subspecies, other obviously bad names

length(gbotb_unmatched_species)/length(gbotb_unmatched) #63 percent potential names 
gbotb_unmatched_species <- gbotb_unmatched[which(!sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2 & !grepl(pattern = "\\.",x = gbotb_unmatched))]


#gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched_species)

gbotb_unmatched_tnrs2
gbotb_unmatched_species
gbotb_unmatched_species_tnrs2 <- gbotb_unmatched_tnrs2[which(gbotb_unmatched_tnrs2$name_matched %in% gbotb_unmatched_species | gbotb_unmatched_tnrs2$verbatim_scientific_name %in% gbotb_unmatched_species),]

length(which(gbotb_unmatched_species %in% gbotb_unmatched_tnrs2$name_matched | gbotb_unmatched_species %in% gbotb_unmatched_tnrs2$verbatim_scientific_name))


length(gbotb_unmatched_species_tnrs2$scrubbed_species_binomial[which(unique(gbotb_unmatched_species_tnrs2$scrubbed_species_binomial) %in% gsub(pattern = "_",replacement = " ",x = spp_to_add))])/length(taxa$binomial)

#plot BIEN taxonomy





taxa_sloppy_tnrs <- sloppy_tnrs(species = taxa$Family)
#how many species added outside fams?






###########################################
fam_mrcas<-NULL
for(i in 1:length(unique(taxa$Family))){

  fam_i<-as.character(unique(taxa$Family)[i])
  taxa_mrca <- getMRCA(phy = pruned_tree,tip = as.character(taxa$binomial[which(taxa$Family ==fam_i)]))  
  bien_mrca <- getMRCA(phy = pruned_tree,tip = gsub(pattern = " ",replacement = "_",x = unique(as.character(bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family ==fam_i)]))))
  og_mrca <- getMRCA(phy = pruned_tree,tip = intersect(x = gbotb$tip.label,y = as.character(taxa$binomial[which(taxa$Family ==fam_i)])))
  
  if(is.null(bien_mrca)){bien_mrca<-NA}
  if(is.null(taxa_mrca)){taxa_mrca<-NA}
  if(is.null(og_mrca)){og_mrca<-NA}

  fam_mrcas<-rbind(fam_mrcas,cbind(fam_i,taxa_mrca,bien_mrca,og_mrca))    
  
}
rm(fam_i,taxa_mrca,bien_mrca,og_mrca)
fam_mrcas<-as.data.frame(fam_mrcas)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

placement_out<-NULL
for(i in 17219:length(pruned_tree$tip.label)){
print(i)
  tip_i<-pruned_tree$tip.label[i]

fam_taxa_i<-as.character(taxa$Family[which(gsub(pattern = " ",replacement = "_",x = taxa$binomial)==tip_i)])


if(is.na(as.numeric(as.character(fam_mrcas$taxa_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))){in_taxa_family<-NA}else{
  if(tip_i %in% extract.clade(phy = pruned_tree,node = as.numeric(as.character(fam_mrcas$taxa_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))$tip.label){
    in_taxa_family<-1
  }else(in_taxa_family<-0 )
  
  }




if(is.na(as.numeric(as.character(fam_mrcas$og_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))){in_og_family<-NA}else{
  
if(tip_i %in% extract.clade(phy = pruned_tree,node = as.numeric(as.character(fam_mrcas$og_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))$tip.label){
  in_og_family<-1
}else(in_og_family<-0 )

}


if(length(which(gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial)==tip_i))==0){in_bien_family <- NA}else{
fam_bien_i<-na.omit(unique(bien_taxonomy$scrubbed_family[which(gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial)==tip_i)]))
if(length(fam_bien_i)>1){in_bien_family<-0}else{ 
if(fam_bien_i%in%fam_mrcas$fam_i){
if(is.na(as.numeric(as.character(fam_mrcas$bien_mrca[which(fam_mrcas$fam_i==fam_bien_i)])))){in_bien_family<-NA}else{
  
if(tip_i %in% extract.clade(phy = pruned_tree,node = as.numeric(as.character(fam_mrcas$bien_mrca[which(fam_mrcas$fam_i==fam_bien_i)])))$tip.label){
  in_bien_family<-1
}else(in_bien_family<-0 )
}
}else{in_bien_family<-NA}
}}
placement_out<-rbind(placement_out,cbind(tip_i,in_taxa_family,in_og_family,in_bien_family))


}
saveRDS(object = placement_out,file = "placement_in_fams.rds")


placement_out <- as.data.frame(placement_out)
placement_out$in_taxa_family<-as.character(placement_out$in_taxa_family)
placement_out$in_og_family<-as.character(placement_out$in_og_family)
placement_out$in_bien_family<-as.character(placement_out$in_bien_family)




length(which(placement_out$in_og_family=="0"))






#joining old PUTS input to info on which species are correctly placed should help figure out issues



#plotting without introduced species outside their genera
library(ggtree)
library(ggplot2)
library(cowplot)
library(ape)

set.seed(2005)
random_tree <- read.tree(list.files(path = "C:/Users/Brian/Desktop/current_projects/DNH_Park/phylos/sunplin_trees/",pattern = ".tre",full.names = T)[runif(n = 1,min = 1,max = 1000)])
plot(random_tree,show.tip.label = F)

status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/status4.csv",stringsAsFactors = F)

pruned_tree<-drop.tip(phy = random_tree,tip = setdiff(x = random_tree$tip.label,y = status$binomial))




tree_genera<-NULL
for(i in gbotb$tip.label){
  
  tree_genera<-unique(c(tree_genera,strsplit(x = i,split = "_")[[1]][1]))
  
}


l1 <- get_legend(ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
plot(l1)

status$Genus %in% tree_genera

t1<-ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=(status %in% c('I') & Genus %in% tree_genera) ), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=(status %in% c('N')) & Genus %in% tree_genera) , color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )


plot(t1)
combined<- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)


ggsave(filename = "DNH_scale_phylo_2.pdf",plot = combined,device = "pdf",path = "C:/Users/Brian/Desktop/",width = 20,height = 20,units = "cm")
ggsave(filename = "DNH_scale_phylo_2.jpg",plot = combined,device = "jpg",path = "C:/Users/Brian/Desktop/",width = 20,height = 20,units = "cm")








##################################################################################

placement_out<-readRDS(file = "placement_in_fams.rds")
placement_out <- as.data.frame(placement_out)
placement_out$in_taxa_family<-as.character(placement_out$in_taxa_family)
placement_out$in_og_family<-as.character(placement_out$in_og_family)
placement_out$in_bien_family<-as.character(placement_out$in_bien_family)
placement_out$tip_i <- as.character(placement_out$tip_i)


puts<-read.table("gbotb_tree.puts")
puts$V1<-as.character(puts$V1)
puts$V2<-as.character(puts$V2)

placement_out$tip_i

merged_puts<-merge(x = placement_out,y = puts,by.x = "tip_i", by.y = "V1",all = T)


problem_placements<-merged_puts[which(merged_puts$in_og_family==0),]






###########################################
library(ape)
library(phytools)
gbotb<-read.tree("C:/Users/Brian/Desktop/current_projects/Park_DNH_scale/Smith_2017_gbotb.tre")
taxa<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/L48_taxa.csv")

length(unique(taxa$binomial)) #initial dataset for L48 had 17,649 species
length(setdiff(x = taxa$binomial,y = gbotb$tip.label))#8958 species needed to be added, the remainder were present


#Of the 17,649 species in our dataset, approximately half (8691) were already present in the phylogeny
8691/17649

#how many genera?

tree_genera<-NULL
for(i in gbotb$tip.label){
  
  tree_genera<-unique(c(tree_genera,strsplit(x = i,split = "_")[[1]][1]))
  
}

taxa$Genus
length(unique(intersect(as.character(taxa$Genus),tree_genera)))/length(unique(taxa$Genus))

#~10 percent of genera missing

setdiff(as.character(taxa$Genus),tree_genera)
length(which(taxa$Genus %in% intersect(as.character(taxa$Genus),tree_genera)))/nrow(taxa) #6 pct of species

spp_to_add <- setdiff(x = gsub(pattern = " ",replacement = "_",x = taxa$binomial),y = gbotb$tip.label)

taxhead<-BIEN:::.BIEN_sql("SELECT * FROM bien_taxonomy LIMIT 1;")


sloppy_tnrs<-function(species, ...){
  
  query = paste("SELECT * FROM bien_taxonomy 
              WHERE scrubbed_species_binomial IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ")  ; ")  
  
  
  return(BIEN:::.BIEN_sql(query, ...))
  
  
  
}

tnrs_out <- sloppy_tnrs(species = gsub(pattern = "_",replacement = " ",x = spp_to_add))

not_scrubbed <- setdiff(x = gsub(pattern = "_",replacement = " ",x = spp_to_add),y = unique(tnrs_out$scrubbed_species_binomial))

# only 84 species were not matched as a "scrubbed name", predominantly accepted names but some "no opinions"

sloppy_tnrs<-function(species, ...){
  
  query = paste("SELECT * FROM bien_taxonomy 
              WHERE scrubbed_species_binomial IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ")  ; ")  
  
  
  return(BIEN:::.BIEN_sql(query, ...))
  
  
  
}


sloppy_tnrs_pt2 <- function(species, ...){
  
  query = paste("SELECT scrubbed_species_binomial, name_matched, verbatim_scientific_name 
                FROM view_full_occurrence_individual 
                WHERE name_matched IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ") 
                OR verbatim_scientific_name IN (", paste(shQuote(species, type = "sh"),collapse = ', '), ") ; ")  
  
  
  return(BIEN:::.BIEN_sql(query, ...))
  
  
  
}



tnrs_out2 <- sloppy_tnrs_pt2(species = not_scrubbed)
tnrs_out2<-unique(tnrs_out2)


tnrs_out2$scrubbed_species_binomial

length(which(as.character(sapply(X = tnrs_out2$scrubbed_species_binomial,FUN = function(x){strsplit(x,split = " ")[[1]][1]}))==
               as.character(sapply(X = tnrs_out2$name_matched, FUN = function(x){strsplit(x,split = " ")[[1]][1]}))))

tnrs_3<-tnrs_out2[which(as.character(sapply(X = tnrs_out2$scrubbed_species_binomial,FUN = function(x){strsplit(x,split = " ")[[1]][1]}))!=
                          as.character(sapply(X = tnrs_out2$name_matched, FUN = function(x){strsplit(x,split = " ")[[1]][1]}))),]


unique(tnrs_out2$scrubbed_species_binomial[which(as.character(sapply(X = tnrs_out2$scrubbed_species_binomial,FUN = function(x){strsplit(x,split = " ")[[1]][1]}))==
                                                   as.character(sapply(X = tnrs_out2$name_matched, FUN = function(x){strsplit(x,split = " ")[[1]][1]})))])

#10 spelling variants


unique(tnrs_3$scrubbed_species_binomial)
#22 species that have been re-assigned


not_scrubbed[which(!not_scrubbed %in% c(tnrs_out2$name_matched,tnrs_out2$verbatim_scientific_name))]
#49 species not matched at all
length(unique(tnrs_out$scrubbed_species_binomial))/length(unique(spp_to_add))*100


#How many species overall match to bien taxonomy?
alltnrs<-sloppy_tnrs(species = gsub(pattern = "_",replacement = " ",x = taxa$binomial))

length(which(!gsub(pattern = "_",replacement = " ",x = taxa$binomial) %in%alltnrs$scrubbed_species_binomial))/length(unique(alltnrs$scrubbed_species_binomial))*100
#less than one half of one percent of names do not match to TNRS scrubbed names

#####################################################
#####################################################
#figure
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

set.seed(2005)
random_tree <- read.tree(list.files(path = "C:/Users/Brian/Desktop/current_projects/DNH_Park/phylos/sunplin_trees/",pattern = ".tre",full.names = T)[runif(n = 1,min = 1,max = 1000)])
plot(random_tree,show.tip.label = F)



library(ggtree)
library(ggplot2)
library(cowplot)

status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/status4.csv",stringsAsFactors = F)

pruned_tree<-drop.tip(phy = random_tree,tip = setdiff(x = random_tree$tip.label,y = status$binomial))

angios <- getMRCA(phy = pruned_tree, tip = intersect(x = status$binomial[which(status$Group=="Angiosperms")],y = pruned_tree$tip.label) )
gymnos <- getMRCA(phy = pruned_tree, tip = intersect(x = status$binomial[which(status$Group=="Gymnosperms")],y = pruned_tree$tip.label) )

#phylo, invasions, angios, scale bar
ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = 5) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+theme(legend.position="right")+
  geom_cladelabel(node = angios, label="", offset = 20,barsize = 3,color="grey")



#We should have a scale bar for the branch lengths.
#A legend might be nice (red = non-native, blue = native).
#Maybe make the colored points a bit bigger (or depict them as bars).
#If it's not too much trouble, maybe label some lower level groups to show we care?
#Save as pdf (or other vector graphic) so that people can zoom in.



l1 <- get_legend(ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
                   geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,position = ggplot2::position_nudge(x = 10))+
                   geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5  )+
                   geom_cladelabel(node = angios, label="Angiosperms", offset = 20,barsize = 3,color="grey")+
                   geom_tippoint(aes(color=status),show.legend = T)+theme(legend.position = "right")+labs(color = "Native Status\n")+
                   scale_color_manual(labels = c("Native", "Non-native"), values = c("blue", "red")))
plot(l1)

t1<-ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = angios, label="Angiosperms", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")


plot(t1)

combined<- ggdraw(t1) +
  draw_plot(plot = l1,x = .8,y = .8,height = .2,width = .2)

library(ggplot2)
ggsave(filename = "DNH_scale_phylo",plot = combined,device = "svg",path = "C:/Users/Brian/Desktop/",width = 40,height = 40,units = "cm")

bien_taxonomy_status<-BIEN::BIEN_taxonomy_species(species = gsub(pattern = "_",replacement = " ",x = status$binomial))
unique(bien_taxonomy_status$order)

orders<-read.table("list_of_orders.txt")

status<-merge(x = status,y = orders,by.x = "Family",by.y = "V2")
colnames(x = status)[length(status)]<-"Order"

bien_taxonomy_status















#Modify above to add most common families
library(BIEN)
bien_taxonomy <- BIEN_taxonomy_species(species = gsub(pattern = "_",replacement = " ",x = pruned_tree$tip.label))
which(bien_taxonomy$scrubbed_species_binomial %in% gsub(pattern = "_",replacement = " ",x = pruned_tree$tip.label))
bien_taxonomy<-unique(bien_taxonomy[c("scrubbed_family","scrubbed_species_binomial")])
#are any species listed in different families?
dups<- bien_taxonomy$scrubbed_species_binomial[which(duplicated(bien_taxonomy$scrubbed_species_binomial))]
#drop these conflicting species from the taxonomy

bien_dups <- bien_taxonomy[which(bien_taxonomy$scrubbed_species_binomial %in% dups),]
bien_dups<-bien_dups[c('scrubbed_species_binomial','scrubbed_family')]
bien_dups<-unique(bien_dups)
bien_taxonomy <- bien_taxonomy[which(!bien_taxonomy$scrubbed_species_binomial %in% dups),]

fam_table<- sort(table(bien_taxonomy$scrubbed_family),decreasing = T)
names(fam_table[1:3])

length(intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Asteraceae")]) 
                 ,y = pruned_tree$tip.label))

length(pruned_tree$tip.label)


aster <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Asteraceae")]) 
                                                    ,y = pruned_tree$tip.label) )


fab <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Fabaceae")]) 
                                                  ,y = pruned_tree$tip.label) )


poa <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Poaceae")]) 
                                                  ,y = pruned_tree$tip.label) )



#Asteraceae, Fabaceae, Fagaceae, Rosaceae, Poaceae, Orchidaceae, Apiaceae, Lamiaceae, Brassicaeae, etc.



ggtree(tr = pruned_tree,layout = "circular") %<+% status+geom_treescale(y = length(pruned_tree$tip.label)-1000,offset = 100,width = 100) + 
  geom_tippoint(aes(subset=status %in% c('I')), color='red', alpha=.8,size=2,position = ggplot2::position_nudge(x = 25))+
  geom_tippoint(aes(subset=status %in% c('N')), color='blue', alpha=.5,size=2, position = ggplot2::position_nudge(x = 10)  )+
  geom_cladelabel(node = angios, label="Angiosperms", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="grey")+
  geom_cladelabel(node = aster, label="Asteraceae", offset = 35,angle = 90,offset.text = 50,barsize = 3,color="brown")


poa <- getMRCA(phy = pruned_tree, tip = intersect(x = gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family == "Poaceae")]) 
                                                  ,y = pruned_tree$tip.label) )



################################


gbotb_tnrs_out<-sloppy_tnrs(species= gsub(pattern = "_",replacement = " ",x = gbotb$tip.label))
unique(gbotb_tnrs_out$scrubbed_taxonomic_status)
length(unique(gbotb_tnrs_out$scrubbed_species_binomial))
length(gbotb$tip.label)
gbotb_unmatched<-setdiff(gsub(pattern = "_",replacement = " ",x = gbotb$tip.label),unique(gbotb_tnrs_out$scrubbed_species_binomial))


gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched)

length(which(sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2))/length(gbotb_unmatched)
#33% sub species, varieties, etc.

length(grep(pattern = "\\.",x = gbotb_unmatched))/length(gbotb_unmatched)
#35% .sp, etc.

#47% subspecies, other obviously bad names


gbotb_unmatched_species <- gbotb_unmatched[which(!sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2 & !grepl(pattern = "\\.",x = gbotb_unmatched))]

gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched_species)

#Do any non-matched taxa correspond to matched taxa?


missing_spp<- setdiff(x = taxa$binomial,y = gbotb$tip.label)

missing_spp_tnrs_out<- sloppy_tnrs(species = gsub(pattern = "_",replacement = " ",x = missing_spp))
#how many were matched?
length(unique(missing_spp_tnrs_out$scrubbed_species_binomial))#8874/8958 # all good names

phylo_spp<- setdiff(y = as.character(taxa$binomial),x = gbotb$tip.label)
length(grep(pattern = "\\.",x = gbotb_unmatched))/length(gbotb_unmatched)\

phylo_spp <-  phylo_spp[which(!sapply(X = phylo_spp,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2 & !grepl(pattern = "\\.",x = phylo_spp))]


#how many species from the phylo (that didn't match)





gbotb_tnrs_out<-sloppy_tnrs(species= gsub(pattern = "_",replacement = " ",x = gbotb$tip.label))
unique(gbotb_tnrs_out$scrubbed_taxonomic_status)
length(unique(gbotb_tnrs_out$scrubbed_species_binomial))#69960 species in the phylogeny match to TNRS scrubbed names
length(unique(gbotb_tnrs_out$scrubbed_species_binomial))/length(gbotb$tip.label) #89%

gbotb_unmatched<-setdiff(gsub(pattern = "_",replacement = " ",x = gbotb$tip.label),unique(gbotb_tnrs_out$scrubbed_species_binomial))
#of the remaining ~11%:

#gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched)

length(which(sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2))/length(gbotb_unmatched)
#33% sub species, varieties, etc.

length(grep(pattern = "\\.",x = gbotb_unmatched))/length(gbotb_unmatched)
#35% .sp, etc.

#37% subspecies, other obviously bad names

length(gbotb_unmatched_species)/length(gbotb_unmatched) #63 percent potential names 
gbotb_unmatched_species <- gbotb_unmatched[which(!sapply(X = gbotb_unmatched,FUN = function(x){length(strsplit(x = x,split = " ")[[1]])})>2 & !grepl(pattern = "\\.",x = gbotb_unmatched))]


#gbotb_unmatched_tnrs2<- sloppy_tnrs_pt2(species = gbotb_unmatched_species)

gbotb_unmatched_tnrs2
gbotb_unmatched_species
gbotb_unmatched_species_tnrs2 <- gbotb_unmatched_tnrs2[which(gbotb_unmatched_tnrs2$name_matched %in% gbotb_unmatched_species | gbotb_unmatched_tnrs2$verbatim_scientific_name %in% gbotb_unmatched_species),]

length(which(gbotb_unmatched_species %in% gbotb_unmatched_tnrs2$name_matched | gbotb_unmatched_species %in% gbotb_unmatched_tnrs2$verbatim_scientific_name))


length(gbotb_unmatched_species_tnrs2$scrubbed_species_binomial[which(unique(gbotb_unmatched_species_tnrs2$scrubbed_species_binomial) %in% gsub(pattern = "_",replacement = " ",x = spp_to_add))])/length(taxa$binomial)

#plot BIEN taxonomy





taxa_sloppy_tnrs <- sloppy_tnrs(species = taxa$Family)
#how many species added outside fams?






###########################################
fam_mrcas<-NULL
for(i in 1:length(unique(taxa$Family))){
  
  fam_i<-as.character(unique(taxa$Family)[i])
  taxa_mrca <- getMRCA(phy = pruned_tree,tip = as.character(taxa$binomial[which(taxa$Family ==fam_i)]))  
  bien_mrca <- getMRCA(phy = pruned_tree,tip = gsub(pattern = " ",replacement = "_",x = unique(as.character(bien_taxonomy$scrubbed_species_binomial[which(bien_taxonomy$scrubbed_family ==fam_i)]))))
  og_mrca <- getMRCA(phy = pruned_tree,tip = intersect(x = gbotb$tip.label,y = as.character(taxa$binomial[which(taxa$Family ==fam_i)])))
  
  if(is.null(bien_mrca)){bien_mrca<-NA}
  if(is.null(taxa_mrca)){taxa_mrca<-NA}
  if(is.null(og_mrca)){og_mrca<-NA}
  
  fam_mrcas<-rbind(fam_mrcas,cbind(fam_i,taxa_mrca,bien_mrca,og_mrca))    
  
}
rm(fam_i,taxa_mrca,bien_mrca,og_mrca)
fam_mrcas<-as.data.frame(fam_mrcas)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

placement_out<-NULL
for(i in 17219:length(pruned_tree$tip.label)){
  print(i)
  tip_i<-pruned_tree$tip.label[i]
  
  fam_taxa_i<-as.character(taxa$Family[which(gsub(pattern = " ",replacement = "_",x = taxa$binomial)==tip_i)])
  
  
  if(is.na(as.numeric(as.character(fam_mrcas$taxa_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))){in_taxa_family<-NA}else{
    if(tip_i %in% extract.clade(phy = pruned_tree,node = as.numeric(as.character(fam_mrcas$taxa_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))$tip.label){
      in_taxa_family<-1
    }else(in_taxa_family<-0 )
    
  }
  
  
  
  
  if(is.na(as.numeric(as.character(fam_mrcas$og_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))){in_og_family<-NA}else{
    
    if(tip_i %in% extract.clade(phy = pruned_tree,node = as.numeric(as.character(fam_mrcas$og_mrca[which(fam_mrcas$fam_i==fam_taxa_i)])))$tip.label){
      in_og_family<-1
    }else(in_og_family<-0 )
    
  }
  
  
  if(length(which(gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial)==tip_i))==0){in_bien_family <- NA}else{
    fam_bien_i<-na.omit(unique(bien_taxonomy$scrubbed_family[which(gsub(pattern = " ",replacement = "_",x = bien_taxonomy$scrubbed_species_binomial)==tip_i)]))
    if(length(fam_bien_i)>1){in_bien_family<-0}else{ 
      if(fam_bien_i%in%fam_mrcas$fam_i){
        if(is.na(as.numeric(as.character(fam_mrcas$bien_mrca[which(fam_mrcas$fam_i==fam_bien_i)])))){in_bien_family<-NA}else{
          
          if(tip_i %in% extract.clade(phy = pruned_tree,node = as.numeric(as.character(fam_mrcas$bien_mrca[which(fam_mrcas$fam_i==fam_bien_i)])))$tip.label){
            in_bien_family<-1
          }else(in_bien_family<-0 )
        }
      }else{in_bien_family<-NA}
    }}
  placement_out<-rbind(placement_out,cbind(tip_i,in_taxa_family,in_og_family,in_bien_family))
  
  
}
saveRDS(object = placement_out,file = "placement_in_fams.rds")


placement_out <- as.data.frame(placement_out)
placement_out$in_taxa_family<-as.character(placement_out$in_taxa_family)
placement_out$in_og_family<-as.character(placement_out$in_og_family)
placement_out$in_bien_family<-as.character(placement_out$in_bien_family)




length(which(placement_out$in_og_family=="0"))













