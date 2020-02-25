#sunplin functions using node labels


#helper functions for sunplin and phylogenies
#library(BIEN)

#sp_fam: dataframe where 1st column is species, second is family

#function to generate puts file
get_put_info_node_labels<-function(sp_fam,phylogeny,genus_only_addition=FALSE){
  
  sp_fam$put_level<-NA
  sp_fam$put<-NA
  sp_fam[,1] <- gsub(pattern = "_",replacement = " ",x = sp_fam[,1]) #in case underscores are used instead of spaces
  
  for( i in 1:nrow(sp_fam)){
    genus<-unlist(strsplit(x = as.character(sp_fam[i,1]),split = " "))[1]
    
    #check whether species is in phylogeny
    if(gsub(pattern = " ",replacement = "_",x = as.character(sp_fam[i,1])) %in% phylogeny$tip.label  ){
      sp_fam$put_level[i]<-"present"
      sp_fam$put[i]<-"present"
      next
    }
    
    
    #Check whether genus is present multiple times
    if(length(grep(pattern = paste(genus,"_",sep = ""),x = phylogeny$tip.label))>1 & is.na(sp_fam$put[i])){
      sp_fam$put_level[i]<-"mrca_genus"
      sp_fam$put[i]<-genus 
      next
      }#if
    
    #Check whether genus is present at least once
    if(length(grep(pattern = paste(genus,"_",sep = ""),x = phylogeny$tip.label))==1 & is.na(sp_fam$put[i])){
      sp_fam$put_level[i]<-"congener"
      sp_fam$put[i]<-phylogeny$tip.label[grep(pattern = paste(genus,"_",sep = ""),x = phylogeny$tip.label)] 
      next
      }#if
    
    
    
    if(genus_only_addition==T){
      if(is.na(sp_fam$put_level[i])){
        sp_fam$put_level[i]<-"remove"
        sp_fam$put[i]<-"remove"
      }
      next
      
    }
      
    
    
    
    
    
    
    
    
    #Check whether family is present multiple times (and genus is not present)
    
    fam_i <- as.character(sp_fam[i,2])
    
    #sanity check for # families != 1
    if(length(fam_i)>=1 & is.na(sp_fam$put_level[i])){  if(length(fam_i)>1 | fam_i == "Unknown"){stop("too many families or only unknown")}  
      
    #first, we'll try to assign to families based on nodes
    if(sp_fam[i,2] %in% phylogeny$node.label){stop("Brian, write family node code")}
      
    #failing that, we'll use 
      
        stop("The code below here is broken")
      
      
      spp_in_fam_i <- 
        phylogeny$tip.label[which(phylogeny$tip.label %in% 
                                    gsub(pattern = " ",replacement = "_",x = sp_fam$binomial[which(sp_fam$Family==fam_i)]) )]
      
      #if genus attempts havent worked and there are multiple confamilials, use family mrca
      if( length(spp_in_fam_i)>1 &  is.na(sp_fam$put_level[i])){
        sp_fam$put_level[i]<-"mrca_family"
        sp_fam$put[i]<-fam_i  }#if
      
      
      if( length(spp_in_fam_i)==1 &  is.na(sp_fam$put_level[i])){
        sp_fam$put_level[i]<-"confamilial"
        sp_fam$put[i]<-spp_in_fam_i  }#if
      
      
      #If none of these have worked, the species will retain an NA, and be removed later
      
      rm(spp_in_fam_i,fam_i)
    }#only go through family stuff if we have a family in the taxonomy
    
    if(is.na(sp_fam$put_level[i])){
      sp_fam$put_level[i]<-"remove"
      sp_fam$put[i]<-"remove"
    }
    
    
    
  }#for i loop
  #rm(fam_i,i)
  
  return(sp_fam)
  
}#end function

