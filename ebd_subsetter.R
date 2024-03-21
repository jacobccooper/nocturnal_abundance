# filepath is where ebd files are contained, where new data will be written
# species should be scientific, genus - species

ebd_subsetter <- function(bird_file,sample_file,
                          species,filepath,filename){
  ebd <- auk_ebd(paste0(filepath,bird_file), 
                 file_sampling = paste0(filepath,sample_file))
  
  # filter data
  # declare filters here
  ebd_filters <- ebd %>%
    # define species here
    # English or Scientific Name
    auk_species(species) %>% 
    # restrict to the standard traveling and stationary count protocols
    auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
    # restrict to complete checklists
    auk_complete()
  
  # rename files for each species
  # here, "rtha" for Red-tailed Hawk
  f_ebd <- file.path(filepath,paste0("ebd_",filename))
  f_sampling <- file.path(filepath,paste0("smp_",filename))
  
  # only run if the files don't already exist
    auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling,
               overwrite = T)
  
  # observed and not observed
  # set to more recent years
  zerofill_file <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE) %>% 
    mutate(
      # convert X to NA
      observation_count = if_else(observation_count == "X", 
                                  NA_character_, observation_count),
      observation_count = as.integer(observation_count),
      # effort_distance_km to 0 for non-travelling counts
      effort_distance_km = if_else(protocol_type != "Traveling", 
                                   0, effort_distance_km),
      # convert time to decimal hours since midnight
      # kept here - needed for modeling
      # time_observations_started = time_to_decimal(time_observations_started),
      # split date into year and day of year
      # we keep date different - kept here for notes
      # year = year(observation_date),
      # day_of_year = yday(observation_date)
    ) %>%
    filter(year(observation_date) >= 2000,
           year(observation_date) <= 2023) %>%
    # remove na from duration; some exist!
    drop_na(duration_minutes) %>% 
    drop_na(time_observations_started) %>%
    drop_na(observation_count) %>% 
    ebd_time_data()
  
  write_rds(zerofill_file,file=paste0(filepath,filename,".rds"))
}

