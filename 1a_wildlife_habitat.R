# THIS MAKES A DECADAL TIME SERIES OF HABITAT AND CLIMATE MODELS TOGETHER
library(raster)
library(tidyverse)
library(furrr)
library(tmap)
library(RSQLite)
library(dbplyr)
library(viridis)

setwd("/Volumes/thexternal/iLand_gye/wildlife_pc")

# Make sure select is from dplyr package
select <- dplyr::select

# Where are the untarred iLand outputs located?
output_dir <- file.path('untar_output')

# Things to specify
species <- c('bbwo','mart','tahu')
decades = seq(5,95,10)
# Re ran 2010 after including spin-up fires (1984-2005)
# decades = decades[1]
simulations <- seq(0,19)

# This is a dataframe listing all landscape-GCM-RCP-iteration combinations for function to iterate over
scenarios <- 
  expand.grid('landscape' = c('grte','nynp','synp','wynp','sgye'), # 'grte','nynp','nynp','synp','wynp','sgye'
              'gcm'= c('canesm2','hadgem2es'), # 'canesm2', 'hadgem2cc','hadgem2es',,'hadgem2cc',,'canesm2'
              'rcp' = c('45','85')) #'45',

# TESTING
# landscape = c('grte')
# gcm= c('canesm2')
# rcp = c('85')

# A set of nested functions produce a set of binary rasters of habitat suitability
# One raster per Landscape/GCM/RCP/Simulation/Decade/Species (!!!)
habitat_fn <- function(landscape,gcm,rcp,sim){
  
  # Status update... check if scenario is already complete, if so, skip
  scenario_name <- paste0(landscape,'_',gcm,'_',rcp)
  
  # if(length(
  #   list.files(
  #     path = 'habitat_results/decadal',
  #     pattern = c(paste0('habitat_',landscape,'_',gcm,'_',rcp), '.*'))
  #   ) ==  600){
  # 
  #   print(paste0(paste0(landscape,'_',gcm,'_',rcp), ' already complete'))
  #   return(NULL)
  # }
  
  print(paste('Running:', scenario_name))
  
  # Read in iLand environment raster
  iland_env <- raster(paste0('../stand_grids/',landscape,'_landscape_env_grid.asc'), crs = '+init=EPSG:26912')
  iland_poly <- rgdal::readOGR(paste0('../stand_grids/',landscape,'_polygon/',landscape,'_polygon.shp'))
  
  # Identify cells burned by 'planned fires' - save these IDs for use in BBWO rules
  # B/c how the modeling was set up, this will include only fires that happened 2006-2016. 
  imposed_1984_2005 <- raster('imposed_fires_1984_2005.txt', crs = '+init=EPSG:26912')
  imposed_2006_2016 <- raster('imposed_fires_2006_2016.txt', crs = '+init=EPSG:26912')
  
  planned <- mosaic(imposed_1984_2005,imposed_2006_2016, fun = max)
  # Concert raster of years to spatial dataframe
  fire_year_cells <- rasterToPoints(planned, function(x) x > 0, spatial = T)
  # Use spatial dataframe to idenitfy RIDs
  burned_rid <- raster::extract(iland_env, fire_year_cells)
  # Wrangling
  planned_fire_cells <- fire_year_cells %>% 
    as_tibble() %>% 
    mutate(rid = burned_rid) %>% 
    filter(!is.na(rid)) %>% 
    dplyr::select(rid, fire_year = layer) 
  # Remove things that are no longer needed for memory sake
  rm(planned, fire_year_cells, burned_rid)
  
  # The meat: loop through simulations and identify habitat based on rules
  sim_fn <- function(sim){
    
    # if(length(
    #   list.files(
    #     path = 'habitat_results/decadal',
    #     pattern = c(paste0('habitat_',landscape,'_',gcm,'_',rcp,'_',sim), '.*'))
    #   ) ==  30){
    #   
    #   print(paste0(paste0(landscape,'_',gcm,'_',rcp), ' already complete'))
    #   return(NULL) 
    # }
    # 
    # Define paths to outputs etc
    output_path <- file.path('/Volumes/thexternal/iLand_gye',output_dir,
                             paste(scenario_name,'output', sim, sep = '_'),
                             'NR','output')
    database <- file.path(output_path,'output.sqlite')
    
    # Check if the output file exists, if not, skip
    if( !file.exists(database) ){
      print(paste('Scenario', landscape, gcm, rcp, sim, 
                  'output sqlite table missing!'))
      return(NULL)
      # Otherwise, connect to the database for scenario_i
    } else {
      # Connect to the sqlite database
      output <- DBI::dbConnect(RSQLite::SQLite(), database)
      print(paste('Analyzing simulation', sim))
    }
    
    # For each scenario, iterate through each timestep (decades)
    decade_fn <- function(model_year){
      
      # if(length(
      #   list.files(
      #     path = 'habitat_results/decadal',
      #     pattern = c(paste0('habitat_',landscape,'_',gcm,'_',rcp,'_',sim,'_',model_year+2005),'.*'))
      #   ) ==  3){
      #   
      #   print(paste0(paste0(landscape,'_',gcm,'_',rcp), ' already complete'))
      #   return(NULL) 
      # }
      
      
      # For BBWO, iterate each year of the decade proceeding each timestep
      bbwo_fun <- function(year_x){
        
        print(paste0('Running year ',year_x,' (',year_x+2005,')'))
        
        # Pull out the planned fires for this year
        planned_cells <- planned_fire_cells %>% 
          filter(fire_year == year_x+2005)
        
        # List the iLand fire raster files for this year
        fire_files <- list.files(path = file.path(output_path,'fire_rasters'), 
                                 pattern = c(paste0('crownkill_',year_x,'_'), '[:digit:]+'), 
                                 full.names = T)
        
        # Get RIDs from simulated fires this year
        if (length(fire_files) > 0) {
          # Read in those rasters
          fire_rast <- map(fire_files, raster, crs = '+init=EPSG:26912')
          
          # Merge all fire rasters from same year into 1 raster 
          # Use mosaic to specify how to handle overlapping values (use max)
          if ( length(fire_rast) > 1 ){
            #print(paste(length(fire_rast), 'fires in year', year_x))
            fire_rast.mosaicargs <- fire_rast
            fire_rast.mosaicargs$fun <- max
            fire_rast <- do.call(mosaic, fire_rast.mosaicargs)
          } else {
            fire_rast <- fire_rast[[1]]
          }
          # Create spatial points dataframe listing the cells that burned at high severity
          burned_cells <- rasterToPoints(fire_rast, function(x) x > 0.5, spatial = T)
          # Remove fire_rast now, for memory
          rm(fire_rast)
          
          # Extract the resource unit ids of cells that burned at high severity (from the landscape environment grid)
          sim_fire_rids <- raster::extract(iland_env, burned_cells)
        } else {
          sim_fire_rids <- NA
        }
        
        # Get RIDs from planned fires this year
        if (nrow(planned_cells[,'rid']) > 0) {
          # Combine them with simulated fire IDs (they may be NULL)
          burned_rids <- c(sim_fire_rids, as.numeric(planned_cells[['rid']]))
        } else {
          # Store just the simulate fire IDs if not
          burned_rids <- sim_fire_rids
        }
        
        
        # If there were fires of any kind, use them to filter data, otherwise save no data
        if (length(burned_rids) > 0){
          output_tbl_raw <- tbl(output, 'output') %>%
            # BURNED; LIVE 1 YEAR BEFORE FIRE;
            filter(rid %in% burned_rids, year == year_x-1) %>%
            # Create some new columns summarizing conifers
            mutate(con_dens = tot_count - potr_dens) %>%
            select(rid, year, tot_ba, con_dens, contains('dbh'), -contains('dead')) %>%
            # CONIFER DENSITY >200/ha; AVG DBH > 20 cm
            filter(con_dens >= 200) %>%
            collect()
          
          # Do more calculatin and filterin
          output_tbl <- output_tbl_raw %>%
            gather(species, dbh, contains('dbh')) %>%
            # Convert zeros to NA (they are not zeros!!)
            mutate(dbh = ifelse(dbh == 0, NA, dbh)) %>%
            group_by(rid, year) %>%
            summarise(dbh =  mean(dbh, na.rm = T), .groups = 'drop_last') %>%
            ungroup() %>%
            # AVG DBH > 20 cm
            filter(dbh > 20) 
          
          bbwo_rids <- as.numeric(output_tbl[['rid']])
        } else {
          bbwo_rids <- NA
        }
        rm(output_tbl_raw)
        return(burned_rids)
      } 
      
      # Each year of the decade proceeding each timestep
      bbwo_years <- seq(model_year-10,model_year) # 11,95
      
      # Run the BBWO habitat function over the 10 years prior to the map year
      print(paste('Starting BBWO habitat module for decade ', model_year))
      bbwo_list <- map(bbwo_years, bbwo_fun)
      bbwo_rids <- as.integer(unlist(bbwo_list))
      
      # Other species 
      print(paste('Running marten and squirrel habitat module ', model_year))
      # Bring the output table into R. This is the slow part... not sure how else
      output_tbl_raw <- tbl(output, 'output') %>% 
        filter(year == model_year) %>% 
        # Create some new columns summarizing conifers 
        mutate(con_dens = tot_count - potr_dens) %>% 
        dplyr::select(rid, year, tot_ba, tot_count, con_dens, contains('dbh'), -contains('dead'), age_p90) %>% 
        collect() 
      
      # Do more calclulatin and filterin
      habitat_results <- output_tbl_raw %>%
        gather(species, dbh, contains('dbh')) %>% 
        # Filter only rows with DBH > 0 (zero for dbh means this species was not present)
        # Remove POTR 
        filter(dbh > 0, species != 'potr_dbh') %>% 
        # TAHU RULE #1: Cone production by species
        mutate(tahu = ifelse(species == 'abla_dbh' & dbh >= 30, 1,
                             ifelse(species == 'pien_dbh' & dbh >= 38, 1, 
                                    ifelse(species == 'psme_dbh' & dbh >= 20, 1, 
                                           ifelse(species == 'pics_dbh' & dbh >= 20, 1,
                                                  ifelse(species == 'pial_dbh', 1, 0)))))) %>%
        # MARTEN filtering
        # RULE 1:  >50 trees/ha; RULE 2: >18 m2/ha basal area; RULE 3: stand age > 75 years old
        mutate(mart = ifelse(tot_count >= 50 & tot_ba > 18 & age_p90 >= 75, 1, 0),
               # TAHU RULE 2: tree density > 50
               tahu = ifelse(tahu == 1 & tot_count >= 50, 1, 0),
               # BBWO rids from above will be filtered below... keep for now
               bbwo = ifelse(rid %in% bbwo_rids, 1, 0)) %>%
        group_by(rid) %>% 
        # Calculate the mean dbh among all species, and tally the number of species in the stand
        summarise(mean_dbh =  mean(dbh, na.rm = T),
                  con_dens = mean(con_dens, na.rm = T), # THIS MIGHT BE BBWO ISSUE
                  n_species = n(),
                  tahu = max(tahu),
                  mart = max(mart),
                  bbwo = max(bbwo),
                  .groups = 'drop_last') %>% 
        # MARTEN RULE #4: AVG DBH > 20 cm
        mutate(mart = ifelse(mart == 1 & mean_dbh >= 20, 1, 0),
               # TAHU RULE #3: >1 species of conifer       
               tahu = ifelse(tahu == 1 & n_species > 1, 1, 0))
      
      # For memory?
      rm(output_tbl_raw)
      
      # # Extract just the RIDs of habitat for each species - this will interface with climate output
      # habitat_rids <- map(species, function(x) {
      #   habitat_results[['rid']][habitat_results[,x] == 1]}) %>%
      #   set_names(species)
      
      # Create and save the rasters for each iteration
      combined_rast <- walk(species, function(x) {
        hab_rast <- reclassify(iland_env, cbind(habitat_results[['rid']], habitat_results[[x]]))
        hab_rast[hab_rast > 1] = 0
        writeRaster(hab_rast, 
                    paste0('habitat_results/decadal/habitat_',
                           scenario_name,'_',sim,'_',model_year+2005,'_',x,'.tif'), 
                    overwrite = T)}) 
    }
    
    # decades <-  c(5,15)  
    walk(decades, decade_fn) 
    # Disconnect sqlite table now
    dbDisconnect(output)
  }
  
  # Map the sim_fn across all simulations
  walk(simulations, sim_fn) 
}

pwalk(scenarios, habitat_fn)
