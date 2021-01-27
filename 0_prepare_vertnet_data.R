library(tidyverse)
library(rvertnet)
library(spThin)
library(rnaturalearth)
library(rnaturalearthdata)

# Load, thin, map species occurrence data --------------------------------------
# For better control, define these here:
spp_abbrv <- c('bbwo','mart','tahu')

# Request and download VertNet occurrence data  
# *************************** Completed on 9/29 ********************************
# To retrieve more than 1000 records, must submit request and then download CSV file
# bigsearch(genus = "Picoides", specificepithet = "arcticus*", rfile = 'bbwo_vert_data.csv', email = 'hoecker@wisc.edu') # arcticus
# bigsearch(genus = "Martes", specificepithet = "americana*", rfile = 'mart_vert_data.csv', email = 'hoecker@wisc.edu')
# bigsearch(genus = "Tamiasciurus", specificepithet = "hudsonicus*", rfile = 'tahu_vert_data.csv', email = 'hoecker@wisc.edu') # hudsonicus

# Get citation information
ref <- map_df(paste0('vertnet_data/',spp_abbrv,'_vert_data.csv'), read_csv)

citations <- ref %>% 
  select(bibliographiccitation, references) %>% 
  separate(references, into = c('url',NA), sep = '\\?') %>% 
  separate(bibliographiccitation, into = c('museum','collection','x1','x2'), sep = '\\.') %>% 
  separate(x1, into = c(NA,'part2','part3','part4','part5', 'part6', 'part7','part8','part9','part10'), sep = ' ') %>% 
  mutate(across(part3:part10, ~ ifelse(part2 == 'Record', ' ', .x))) %>% 
  mutate(across(part3:part10, ~ ifelse(is.na(.x), ' ', .x))) %>% 
  # mutate(part3 = ifelse(part2 == 'Record', NA, part2),
  #        part4 = ifelse(part2 == 'Record', NA, part4)) %>% 
  distinct(museum, collection, part2, .keep_all = TRUE) %>% 
  mutate(museum = ifelse(part2 != 'Record', paste(museum,collection), museum),
         collection = ifelse(part2 != 'Record', paste(part2,part3,part4,part5,part6,part7,part8,part9,part10), collection)) %>% 
  distinct(museum, collection, url)

write_csv(citations, 'vertnet_citations.csv')

# North American map outline


# # Run the function over the set of prefixes
# occ_records <- map(spp_abbrv, thin_vert) %>% 
#   set_names(spp_abbrv)

# Function to read in VertNet data, thin it to 1 km, and return the optimal thinned dataset as a sf object
# 'north_america' is hard-coded into this function and performs intersection of points to remove bad points
occ_records <- 
  purrr::map(spp_abbrv, function(species){
    
    world <- ne_countries(scale = "medium", returnclass = "sf")
    north_america <- world %>% 
      filter(continent == 'North America')
    
    vert_thin_list <- read_csv(paste0('vertnet_data/',species,'_vert_data.csv')) %>% 
      # Remove other species of 
      # Remove records with positional accuracy < 5000 m 
      filter(coordinateuncertaintyinmeters < 5000) %>% 
      # Spatially thin records within 1 km (the resolution of the climate data)
      thin(., lat.col = 'decimallatitude', long.col = 'decimallongitude',
           spec.col = 'genus',
           # KEY ARG: THIN RECORDS WITHIN 250 m (0.25 km)
           thin.par = 1,
           reps = 3,
           write.files = F,
           locs.thinned.list.return = T)
    
    # Pick the longest rep (most records within min buffer)
    max_recs <- max(map_dbl(vert_thin_list, ~ length(.x[,1])))
    first_longest <-
      which(map_dbl(vert_thin_list, ~ length(.x[,1]) == max_recs) == 1) %>% 
      min()
    
    # Use the largest thinned dataset
    thin_df <- vert_thin_list[[first_longest]] %>% 
      filter(!is.na(Longitude)) %>% 
      # Convert to a simple features (spatial) object
      st_as_sf(., coords = c('Longitude', 'Latitude')) %>% 
      `st_crs<-`(4326)  %>% 
      # st_transform(crs = '+init=EPSG:26912')
      # Remove points that don't overlap with North America
      st_intersection(., north_america)
    
    return(thin_df)}) %>% 
  set_names(spp_abbrv)


saveRDS(occ_records, 'vertnet_data/vertnet_occ_data_thin_1km.rds')
