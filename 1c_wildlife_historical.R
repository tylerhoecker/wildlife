# SCRIPT TO CREATE COMBINED NICHE RASTERS FOR 'HISTORICAL' PERIOD,
# USING HABITAT FROM 2020 AND HISTORICAL CLIMATE
library(raster)
library(tidyverse)

# Make sure select is from dplyr package
select <- dplyr::select

# Things to specify
species <- c('bbwo','mart','tahu')
simulations <- seq(0,19)


# TESTING
# ls = c('grte')
# gcm= c('canesm2')
# rcp = c('45')

# Function to read in the results of the habitat (1a) and climate (1b) envelope modeling scripts
landscapes = c('grte','nynp','synp','wynp','sgye')

walk(landscapes, function(ls){
  
  # Read iLand environment
  iland_env <- raster(paste0('../stand_grids/', ls,'_landscape_env_grid.asc'),
                      crs = '+init=EPSG:26912')
  
  # Extract points as sp and just RIDs
  iland_pts <- rasterToPoints(iland_env, spatial = T)
  iland_rids <- rasterToPoints(iland_env)[,3]
  
  # Run over all species
  walk(species, function(sp_name){
    
    # Climate from historical projection
    biomod_hist <- raster(paste0('biomod2ing/',sp_name,
                                 '/proj_hist_gye/proj_hist_gye_',
                                 sp_name,'_ensemble.grd')) %>%
      projectRaster(., crs = '+init=EPSG:26912') %>%
      raster::extract(., iland_pts)
    
    # Reclassify iLand grid with CCASF climate data
    biomod_iland <- reclassify(iland_env, cbind(iland_rids, biomod_hist))
    # Write out this raster
    writeRaster(biomod_iland,
                filename = paste0('biomod2ing/historical_rasters/',ls,'_',sp_name,'.tif'),
                overwrite = T)
  
  # Habitat from earliest timestep
  # Just take all habitat rasters for this species in 2010 and average them
  habitat_files <- list.files(
    path = 'habitat_results/decadal',
    pattern = paste0('habitat_',ls,'.*',2010,'_',sp_name,'.*'), full.names = T)
  habitat_hist <- map(habitat_files, raster) %>% 
    stack() %>% 
    mean(., na.rm = T)
  
  # Quicker than reclassify - just round! Not very consequential here since all scenarios are very similar at this point
  habitat_bin <- round(habitat_hist)
  writeRaster(habitat_bin,
              filename = paste0('habitat_results/historical/',ls,'_',sp_name,'.tif'),
              overwrite = T)
  
  # Combine climate and habitat
  combined_hist <- biomod_iland * habitat_bin
  
  writeRaster(combined_hist,
              filename = paste0('combined_spatial/',ls,'_historical_',sp_name,'.tif'),
              overwrite = T)
  
  })
})

cutoffs <- readRDS('biomod2ing/cutoff_values.Rds')

