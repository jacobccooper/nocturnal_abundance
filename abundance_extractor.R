abundance_extracter <- function(rds_file,abundance_rast,species_name){
  
  observations <- read_rds(rds_file)
  
  observations$time_observations_started <- 
    as.POSIXct(observations$time_observations_started, 
               format = "%H:%M:%S")
  
  abundance <- rast(abundance_rast)
  
  # convert to same coordinate system
  abundance <- abundance %>% project("epsg:4326")
  
  occ <- observations[,c("longitude","latitude")]
  
  vals <- extract(abundance,occ,cells = T)
  colnames(vals)[2] <- "abundance"
  vals <- vals %>%
    cbind(observations$Time_of_Day,
          observations$species_observed) %>%
    rename('time_of_day' = 'observations$Time_of_Day',
           'species_observed' = 'observations$species_observed')
  
  write_rds(vals,paste0(filepath,species_name,"_abundance_data.rds"))
}
