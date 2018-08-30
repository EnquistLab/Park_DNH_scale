#r script to collapse monophyletic subspecies, etc.
#this assumes that any species with more than 1 underscor is a subspecies, variety, cultivar, etc., which may not be correct

collapse_subsp<-function(phylogeny){
    
  
subsp<-phylogeny$tip.label[which(unlist( lapply(X = phylogeny$tip.label,function(x){
if(length(strsplit(x = x,split = "_")[[1]])>2){out<-TRUE}
if(length(strsplit(x = x,split = "_")[[1]])<=2){out<-FALSE}  
out})))]

subsp<-as.data.frame(subsp)
subsp$genus<-unlist(lapply(X = as.character(subsp$subsp),FUN = function(x){
    strsplit(x = x,split = "_")[[1]][1]  
  }))

subsp$species<-unlist(lapply(X = as.character(subsp$subsp),FUN = function(x){
  paste(unlist(strsplit(x = x,split = "_")[[1]][1:2]),collapse = "_"  )
  }))


for(i in 1:length(unique(subsp$species))){

species_i<-unique(subsp$species)[i]
taxa<-as.character(subsp$subsp[which(subsp$species==species_i)])    
taxa<-unique(c(taxa,species_i))
taxa<-taxa[which(taxa%in%phylogeny$tip.label)]#make sure the species and the subspecies are present

if(length(taxa)==1){
  
phylogeny$tip.label[which(phylogeny$tip.label==taxa) ] <-species_i
  
}

if(length(taxa)>1){
  

mrca<-getMRCA(phy = phylogeny,tip = taxa)
clade_i<-extract.clade(phy = phylogeny,node = mrca)

if(length(taxa)==length(clade_i$tip.label) ){

  
tips_to_drop<-taxa[2:length(taxa)]  
tips_to_rename<-taxa[1]      
phylogeny<-drop.tip(phy = phylogeny,tip = tips_to_drop)
phylogeny$tip.label[which(phylogeny$tip.label==tips_to_rename)]<-species_i

}#if can collapse
  
#if you can't collapse it, due to lack of monophyly, do nothing

}#if multiple species to consider
print(paste(round(i/length(unique(subsp$species)*100),digits = 2)," percent done."))

  
}#for i loop

return(phylogeny)


}#function 
  
  
