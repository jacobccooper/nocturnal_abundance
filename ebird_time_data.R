# creating a function for time
# function gets sunrise, sunset
# categorizes data by time of day

ebd_time_data <- function(zerofill_file,maxduration=720){
  # get date-time
  # times saved as UTC into file if POSIXct
  
  ebd_zerofill <- zerofill_file %>%
    filter(duration_minutes <= maxduration)
  longs <- ebd_zerofill$longitude
  lats <- ebd_zerofill$latitude
  
  # get timezone - accurate
  tz_ebd <- tz_lookup_coords(lon = longs,
                             lat = lats,
                             method = "accurate")
  ebd_zerofill$Time_Zone <- tz_ebd
  
  # start date, start time
  x.date <- ebd_zerofill$observation_date
  start.time <- ebd_zerofill$time_observations_started
  
  sun_data <- as.data.frame(cbind(lats,longs))
  # add date to make sure format correct
  sun_data$date <- as.Date(x.date)
  
  sun_data <- sun_data%>%
    select(date,lats,longs) %>%
    rename(lat = lats, lon = longs)
  
  # get sun data
  tz_unique <- unique(tz_ebd)
  
  ebd_zerofill$start_condition <- "Diurnal"
  ebd_zerofill$end_condition <- "Diurnal"
  
  ebd_zerofill$Time_Date <- NA %>% as_datetime()
  
  ebd_zerofill$start_condition <- "Diurnal"
  ebd_zerofill$end_condition <- "Diurnal"
  ebd_zerofill$Time_of_Day <- NA
  
  for(i in 1:length(tz_unique)){
    index <- which(tz_ebd==tz_unique[i])
    sun_times <- getSunlightTimes(data = sun_data[index,],
                                  tz = tz_unique[i])
    all_start <- ymd_hms(paste0(x.date[index]," ",start.time[index]),
                         tz = tz_unique[i])
    ebd_zerofill$Time_Date[index] <- all_start
    
    # get end time
    duration <- ebd_zerofill$duration_minutes[index]
    duration[which(is.na(duration))] <- 0
    duration_hours <- trunc(duration/60,digits=0)
    duration_minutes <- (duration - (duration_hours*60))
    # next step can be problematic; should be ok
    duration_hms <- paste0(duration_hours,":",duration_minutes,":00") %>%
      hms::as_hms()
    # remember - saved as UTC
    end.time <- all_start + duration_hms
    
    # get nocturnal start and ends
    sub_index_start <- which(all_start < sun_times$sunrise|
                               all_start > sun_times$sunset)
    sub_index_end <- which(end.time < sun_times$sunrise|
                             end.time > sun_times$sunset)
    
    ebd_zerofill$start_condition[index[sub_index_start]] <- "Nocturnal"
    ebd_zerofill$end_condition[index[sub_index_end]] <- "Nocturnal"
  }
  
  night <- which(ebd_zerofill$start_condition==ebd_zerofill$end_condition&
                   ebd_zerofill$start_condition=="Nocturnal")
  dawn <- which(ebd_zerofill$start_condition!=ebd_zerofill$end_condition&
                  ebd_zerofill$start_condition=="Nocturnal")
  day <- which(ebd_zerofill$start_condition==ebd_zerofill$end_condition&
                 ebd_zerofill$start_condition=="Diurnal")
  dusk <- which(ebd_zerofill$start_condition!=ebd_zerofill$end_condition&
                  ebd_zerofill$start_condition=="Diurnal")
  
  ebd_zerofill$Time_of_Day[night] <- "Night"
  ebd_zerofill$Time_of_Day[dawn] <- "Dawn"
  ebd_zerofill$Time_of_Day[day] <- "Day"
  ebd_zerofill$Time_of_Day[dusk] <- "Dusk"
  
  ebd_zerofill$Time_Zone <- as.factor(ebd_zerofill$Time_Zone)
  ebd_zerofill$Time_of_Day <- as.factor(ebd_zerofill$Time_of_Day)
  ebd_zerofill$start_condition <- as.factor(ebd_zerofill$start_condition)
  ebd_zerofill$end_condition<- as.factor(ebd_zerofill$end_condition)
  
  return(ebd_zerofill)
}
