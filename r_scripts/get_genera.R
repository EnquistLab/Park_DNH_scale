get_genera <- function(phylogeny){
  
genera <- NULL

    for(i in 1:length(phylogeny$tip.label)){
        
      genera <- unique(c(genera,strsplit(x = phylogeny$tip.label[i],split = "_")[[1]][1]))  
      
    }  
return(genera)    
  
}


###########################

make_species_level <- function(phylogeny){

  for(i in 1:length(phylogeny$tip.label)){
    
    split_tip <- strsplit(x = phylogeny$tip.label[i],split = "_")[[1]] 
        if(length(split_tip)>2){
          
      phylogeny$tip.label[i]  <-paste(split_tip[1:2],collapse = "_")
          
        }
          
      
  }  
  
  return(phylogeny)  
  
  
  
}



##########################


#check synonyms
  #input species_list is species you want on tree, but aren't present

check_synonyms <- function(phylogeny, species_list){
  
  good_names <- intersect(x = species_list,y = phylogeny$tip.label)  
  to_check <- setdiff(species_list,phylogeny$tip.label)
  to_check <- gsub(pattern = "_",replacement = " ",x = to_check)
  
  #this bit checks for the names you're looking for in BIEN
  query = paste("SELECT scrubbed_species_binomial, name_matched
                FROM view_full_occurrence_individual 
                WHERE name_matched IN (", paste(shQuote(to_check, type = "sh"),collapse = ', '), ")  ; ")  

  bien_out<- BIEN:::.BIEN_sql(query)

  # look for the "scrubbed names" in the tree
  possible_matches <- bien_out[which(gsub(pattern = " ",replacement = "_",x = bien_out$scrubbed_species_binomial) %in% phylogeny$tip.label),]
  possible_matches <- unique(possible_matches)
  possible_matches$name_matched <- gsub(pattern = " ",replacement = "_",x = possible_matches$name_matched)
  possible_matches$scrubbed_species_binomial <- gsub(pattern = " ",replacement = "_",x = possible_matches$scrubbed_species_binomial)
  
  #Toss any possible synonyms that match to species we already have (since they would be overwritten in tree)
  possible_matches <- possible_matches[which(!possible_matches$scrubbed_species_binomial %in% species_list),]
  
  for(i in 1:nrow(possible_matches)){
    
    possible_matches$name_matched[i] #species we want in tree
    possible_matches$scrubbed_species_binomial[i] #name that MIGHT be in tree
  
    if(possible_matches$scrubbed_species_binomial[i] %in% phylogeny$tip.label){
      if(!possible_matches$name_matched[i] %in% phylogeny$tip.label){ #this line is just to make sure we aren't adding a species that is already present (shouldn't be necessary, but better to be safe)
        phylogeny$tip.label[which(phylogeny$tip.label==possible_matches$scrubbed_species_binomial[i])] <- possible_matches$name_matched[i]
        
      }} #if statements
  
  } #i loop
  
  return(phylogeny)  


  
  
  
  
}





