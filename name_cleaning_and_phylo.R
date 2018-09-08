#Initial scripts for taxonomic harmonization and phylogeny stuff

library(BIEN)
library(ape)

source("r_scripts/collapse_subsp.R")
source("r_scripts/sloppy_tnrs.R")
source("r_scripts/sunplin_fxs.R")
source("sunplin-functions.r")
gbotb<-read.tree("Smith_2017_gbotb.tre")
taxa<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/L48_taxa.csv")
taxa$binomial<-gsub(pattern = "_",replacement = " ",x = taxa$binomial)
taxa$Family<-as.character(taxa$Family)
taxa$Family[which(taxa$binomial=="Hesperoyucca whipplei")]<-"Asparagaceae"



put_info<-get_put_info(sp_fam = taxa[c("binomial","Family")],phylogeny = gbotb)
problem_sp<-put_info$binomial[which(put_info$put=="remove")]#2 taxa with issues
problem_taxonomy<-BIEN_taxonomy_species(problem_sp)
taxa$Family[which(taxa$binomial%in%problem_sp)]
problem_taxonomy$scrubbed_family




make_puts_input(puts_info = put_info,phylogeny = gbotb,phylogeny_filename = "gbotb_tree.tre",puts_filename = "gbotb_tree.puts")


length(which(put_info$put_level!="present"))


testsputs<-read.csv("gbotb_tree.puts",header = F)
nrow(testsputs)#only one species missing, and its a weird family/species, so this is expected


sunplin_phylo_replicates(put_file = "gbotb_tree.puts",phylogeny_file = "gbotb_tree.tre",
                         output_directory = "sunplin_trees/",
                         output_base_filename = "gbotb_base_sunplin_method2_rep",
                         nrep = 1000, taxa_to_keep = put_info$binomial
                          )
