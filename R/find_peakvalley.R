require(tidyverse)
# Find peaks and valleys takes a series of values in "valuecol" ordered by 
# dates in "datecol" and identifies peaks and valleys in it looking forward
# and backwards a maximum number of timesteps defined by "max.offset"
# Note that it ONLY counts the timesteps included in the dataset. Therefore,
# if a month is missing (e.g. no Mar between Feb and April), then the 
# algorithm treats this as a time difference of 1, not 2. To have it consider
# all months in between, add in the missing months.
find_peakvalley <- function(x, valuecol, datecol, max.offset){
  xvalue <- setNames(x[,c(datecol, valuecol)], c("date", "value")) %>% ungroup()
  first.date <- min(xvalue$date)
  xvalue <- xvalue %>% 
    mutate(value_lag6 = lag(value, n=6, order_by = date)) %>%
    mutate(value_lag5 = lag(value, n=5, order_by = date)) %>%
    mutate(value_lag4 = lag(value, n=4, order_by = date)) %>%
    mutate(value_lag3 = lag(value, n=3, order_by = date)) %>%
    mutate(value_lag2 = lag(value, n=2, order_by = date)) %>%
    mutate(value_lag1 = lag(value, n=1, order_by = date)) %>%
    mutate(value_lead1 = lead(value, n=1, order_by = date)) %>%
    mutate(value_lead2 = lead(value, n=2, order_by = date)) %>%
    mutate(value_lead3 = lead(value, n=3, order_by = date)) %>%
    mutate(value_lead4 = lead(value, n=4, order_by = date)) %>%
    mutate(value_lead5 = lead(value, n=5, order_by = date)) %>%
    mutate(value_lead6 = lead(value, n=6, order_by = date)) %>%
    rownames_to_column(var = "ID")
  
  xvalue_long <- xvalue %>%
    pivot_longer(
      cols = value_lag6:value_lead6, 
      names_to = c("dir", "offset"),
      names_pattern = "value_(lag|lead)([1-6])",
      values_to = "offset_value"
      ) %>% rowwise() %>%
    mutate(incr = if_else(dir == "lag", value > offset_value, offset_value > value)) %>%
    mutate(decr = if_else(dir == "lag", value <= offset_value, offset_value <= value)) %>%
    ungroup()
    
  sumoffsets <- xvalue_long %>% filter(offset <= max.offset) %>%
    group_by(ID, dir) %>%
    summarize(incr = sum(incr, na.rm = TRUE), decr = sum(decr, na.rm = TRUE)) %>%
    ungroup() %>%
    pivot_wider(id_cols = ID, names_from = dir, values_from = c(incr, decr),
                names_glue = "{dir}_{.value}") %>%
    left_join(xvalue %>% dplyr::select(value, value_lag1, date, ID), by = "ID") %>%
    mutate(valley = ifelse(!is.na(value), lag_decr + lead_incr, NA)) %>%
    mutate(peak = ifelse(!is.na(value), lag_incr + lead_decr, NA)) %>%
    rowwise() %>%
    mutate(takeoff = peak == 0 & valley >= 6) %>%
    mutate(takeoff.num = ifelse(is.na(takeoff), 0, as.numeric(takeoff))) %>%
    mutate(end = ifelse(!is.na(value_lag1), value > value_lag1, TRUE) & 
             valley == 0 & peak >= 6 & date != last.date) %>%
    mutate(end.num = ifelse(is.na(end), 0, as.numeric(end))) %>%
    ungroup()
  
  return(sumoffsets %>% dplyr::select(date, takeoff, takeoff.num, end, end.num))
}

# Finds the "runs" of increasing or decreasing data by assuming a valley followed
# by a peak is a run. "pvdata" is the output of the "find_peakvalley" function but
# custom peak and value columns can be identified with the "peak" and "valley"
# arguments. THese columns must consist of 0 and 1s indicating presence of the
# event described by the column.
find_runs <- function(pvdata, start.run = "takeoff.num", end.run = "end.num",  run.prefix = "",
                      include.end = TRUE){
  rundata <- pvdata %>% 
    mutate(startrun = order_by(date, cumsum(!!sym(start.run)))) %>%
    mutate(endrun = order_by(date, cumsum(!!sym(end.run)))) %>%
    { if(include.end) mutate(., endrun = endrun - (!!sym(end.run))) else . } %>%  # subtract 1 from peaks to start the run after the peak
    mutate(runID = paste(run.prefix, str_pad(startrun + endrun, 2, pad = "0"), sep = "")) %>%
    group_by(runID) %>% 
    { if(include.end){
      mutate(., runtype = ifelse(sum(!!sym(start.run), !!sym(end.run))==2, "growth", 
                              ifelse(sum(!!sym(end.run), na.rm = TRUE)==0, "decrease", "other")))
    }else{
      mutate(., runtype = ifelse(sum(!!sym(start.run), na.rm = TRUE) == 1,
                              "growth",
                              ifelse(sum(!!sym(end.run), na.rm = TRUE) == 1,
                                     "decrease", "other"))
      )
    }} %>%
    ungroup()
  return(rundata)
}

# Find the "middle" date needed for defining the first and second halfs of the
# slowing-down regressions. half=1 gives the last date in the first half and 
# half=2 gives the first date of the second half. Halves are defined as
# including the "middle" month if the period has an odd number of months
findmiddle <- function(daterange, half){
  wholemonths <- (((min(daterange) + max(daterange))*12)%/%2)/12
  middlemonth <- (((min(daterange) + max(daterange))*12)%%2)/12
  if(middlemonth==0){
    return(wholemonths + middlemonth)
  }else if(middlemonth>0){
    if(half==1){
      return(wholemonths)
    }else if(half==2){
      return(wholemonths + (1/12))
    }
  }
}