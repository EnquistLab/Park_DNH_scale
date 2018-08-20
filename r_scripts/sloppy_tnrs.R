#Taxonomic scrubbing code for Kenji

#Brian Maitner, 2018

#First,load the BIEN library
library(BIEN)


#Next, load the following function:

#########################
sloppy_tnrs<-function(species_list){
  names_out<-NULL
  
  
  
  for(i in 1:length(species_list)){
    print(i/length(species_list))
    
    name_i<-NULL  
    
    query = paste("SELECT scrubbed_species_binomial, name_matched, verbatim_scientific_name 
              FROM view_full_occurrence_individual 
              WHERE name_matched IN (", paste(shQuote(species_list[i], type = "sh"),collapse = ', '), ") 
              OR verbatim_scientific_name IN (", paste(shQuote(species_list[i], type = "sh"),collapse = ', '), ") ; ")  
    
    
    while(is.null(name_i)){   try(  name_i<- BIEN:::.BIEN_sql(query = query,limit=1  ) ) }  #Using a while() since the sloppy_tnrs function sometimes times out 
    
    if(nrow(name_i)<1){
      name_i<-cbind(NA,NA,NA)
      name_i<-as.data.frame(name_i)
      colnames(name_i)<-colnames(names_out[2:4])
      
    }#if statement
    
    names_out<-rbind(names_out,cbind(species_list[i],name_i))    
  }#for loop
  
  names(names_out)[1]<-"name_supplied"
  
  return(names_out)
  
  
  
}
