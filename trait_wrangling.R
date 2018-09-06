#trait_wrangling
library(BIEN)
status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/status.csv")


status_species<-as.character(status$binomial[which(is.na(status$status))])

occurrence_spp<-BIEN_occurrence_species_2(species = gsub(pattern = "_",replacement = " ",x = status_species))
occurrence_spp<-occurrence_spp[c("scrubbed_species_binomial","native_status_country")]
occurrence_spp<-unique(occurrence_spp)
occurrence_spp$scrubbed_species_binomial
occurrence_spp$scrubbed_species_binomial<-gsub(pattern = " ",replacement = "_",x = occurrence_spp$scrubbed_species_binomial)
merged_status<-merge(x = status,y = occurrence_spp,by.x = "binomial",by.y = "scrubbed_species_binomial",all.x = T)
status_species[which(!status_species%in%occurrence_spp$scrubbed_species_binomial)]#species missing from BIEN native assessment

merged_status$native_status_country
colnames(merged_status)[which(colnames(merged_status)=="native_status_country")]<-"BIEN_native_status_USA"

#standardizing BIEN codes to native vs introduced

unique(merged_status$BIEN_native_status_USA)
merged_status$BIEN_native_status_USA[which(merged_status$BIEN_native_status_USA%in%c('N',"Ne"))]<-"N"
merged_status$BIEN_native_status_USA[which(merged_status$BIEN_native_status_USA%in%c('I',"Ie"))]<-"I"

write.csv(x = merged_status,file = "C:/Users/Brian/Google Drive/DNH_scale/status_w_BIEN.csv",row.names = F)


BIEN_occurrence_species_2<-function(species, native.status=TRUE,natives.only=FALSE,...){
  
  native_<-.native_check(native.status)
  natives_<-.natives_check(natives.only)
  
  # set the query
  query <- paste("SELECT DISTINCT scrubbed_species_binomial",native_$select,"native_status_country 
                 FROM view_full_occurrence_individual 
                 WHERE scrubbed_species_binomial in (", paste(shQuote(species, type = "sh"),collapse = ', '), ")",
                 natives_$query,  "
                 AND country in ('United States') AND higher_plant_group IS NOT NULL AND (is_geovalid = 1 OR is_geovalid IS NULL) 
                 ORDER BY scrubbed_species_binomial ;")
  
  return(.BIEN_sql(query, ...))
  
}



#######################################

#Growth form, etc.
status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/status_w_BIEN.csv")
FIA_status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/FIA/FIA_taxa_status_clean.csv")
USDA_status<-read.csv("C:/Users/Brian/Google Drive/DNH_scale/USDA/unique_USDA.csv")
colnames(FIA_status)
FIA_status$SCIENTIFIC_NAME<-gsub(pattern = " ",replacement = "_",x = FIA_status$SCIENTIFIC_NAME)
FIA_status<-FIA_status[c("SCIENTIFIC_NAME","Tree","Subshrub","Shrub","Vine","Herb","Graminoid","Perennial","Annual","Biennial")]

#Populate columns from FIA data:
status<-merge(x = status,y = FIA_status,by.x = "binomial",by.y = "SCIENTIFIC_NAME",all.x = T)
rm(FIA_status)

#Populate columns from USDA data
USDA_status$Scientific.Name<-gsub(pattern = " ",replacement = "_",x = USDA_status$Scientific.Name)
USDA_status<-USDA_status[which(USDA_status$Scientific.Name%in%status$binomial),]#toss unneeded records

#inelegant for loop

for(i in 1:nrow(USDA_status)){
species_i<-USDA_status$Scientific.Name[i]  
duration_i<-as.character(USDA_status$Duration[i]  )
growth_i<-as.character(USDA_status$Growth.Habit[i])

#fill in annual vs perrennial vs bienniel
try(if(grep(pattern = "AN",duration_i,ignore.case = T)>0){ status$Annual[which(status$binomial==species_i)]<-1  },silent = T)
try(if(grep(pattern = "per",duration_i,ignore.case = T)>0){ status$Perennial[which(status$binomial==species_i)]<-1  },silent = T)
try(if(grep(pattern = "bien",duration_i,ignore.case = T)>0){ status$Biennial[which(status$binomial==species_i)]<-1  },silent = T)

#Fill in growth form
try(if(grep(pattern = "Tree",growth_i,ignore.case = T)>0){ status$Tree[which(status$binomial==species_i)]<-1  },silent = T)
try(if(grep(pattern = "Subshrub",growth_i,ignore.case = T)>0){ status$Subshrub[which(status$binomial==species_i)]<-1  },silent = T)
try(if(grep(pattern = "Shrub",growth_i,ignore.case = F)>0){ status$Shrub[which(status$binomial==species_i)]<-1  },silent = T)
try(if(grep(pattern = "Vine",growth_i,ignore.case = T)>0){ status$Vine[which(status$binomial==species_i)]<-1  },silent = T)
try(if(grep(pattern = "herb",growth_i,ignore.case = T)>0){ status$Herb[which(status$binomial==species_i)]<-1  },silent = T)
try(if(grep(pattern = "graminoid",growth_i,ignore.case = T)>0){ status$Graminoid[which(status$binomial==species_i)]<-1  },silent = T)

print(paste(round(i/nrow(USDA_status)*100,digits = 2)," percent done",sep = ""))
rm(duration_i,growth_i,species_i)  
}#end for loop
rm(i,USDA_status)  
  

#Pull growth form data from BIEN


BIEN_trait_species_2<-function(species, trait, ...){
  
  
  # set the query
  query <- paste("SELECT DISTINCT scrubbed_species_binomial, trait_name, trait_value 
                  FROM agg_traits 
                 WHERE scrubbed_species_binomial in (", paste(shQuote(species, type = "sh"),collapse = ', '), ") and 
                 trait_name in (", paste(shQuote(trait, type = "sh"),collapse = ', '), ")
                 AND higher_plant_group IS NOT NULL
                 ORDER BY scrubbed_species_binomial ;")
  
  return(BIEN:::.BIEN_sql(query, ...))
  
}


bien_growth_forms<-BIEN_trait_species_2(species = gsub(pattern = "_",replacement = " ",x = status$binomial),
                                        trait = c("whole plant growth form diversity"))


#inelegant for loop

for(i in 1:nrow(bien_growth_forms)){
  species_i<-gsub(pattern = " ",replacement = "_",x = bien_growth_forms$scrubbed_species_binomial[i]  )
  growth_i<-as.character(bien_growth_forms$trait_value[i])
  
  ##fill in annual vs perrennial vs bienniel
  #try(if(grep(pattern = "AN",duration_i,ignore.case = T)>0){ status$Annual[which(status$binomial==species_i)]<-1  },silent = T)
  #try(if(grep(pattern = "per",duration_i,ignore.case = T)>0){ status$Perennial[which(status$binomial==species_i)]<-1  },silent = T)
  #try(if(grep(pattern = "bien",duration_i,ignore.case = T)>0){ status$Biennial[which(status$binomial==species_i)]<-1  },silent = T)
  
  #Fill in growth form
  try(if(grep(pattern = "Tree",growth_i,ignore.case = T)>0){ status$Tree[which(status$binomial==species_i)]<-1  },silent = T)
  try(if(grep(pattern = "Shrub",growth_i,ignore.case = T)>0){ status$Shrub[which(status$binomial==species_i)]<-1  },silent = T)
  try(if(grep(pattern = "Vine",growth_i,ignore.case = T)>0){ status$Vine[which(status$binomial==species_i)]<-1  },silent = T)
  try(if(grep(pattern = "herb",growth_i,ignore.case = T)>0){ status$Herb[which(status$binomial==species_i)]<-1  },silent = T)
  #try(if(grep(pattern = "graminoid",growth_i,ignore.case = T)>0){ status$Graminoid[which(status$binomial==species_i)]<-1  },silent = T)
  
  print(paste(round(i/nrow(bien_growth_forms)*100,digits = 2)," percent done",sep = ""))
  rm(growth_i,species_i)  
}#end for loop
rm(i,bien_growth_forms)  

write.csv(x = status,file = "C:/Users/Brian/Google Drive/DNH_scale/status_w_BIEN.csv",row.names = F)



  
  