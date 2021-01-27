library(tidyverse)
library(raster)


# Crop and transform future climate data to iland landscape envrionments
# This only needs to be done once. 
# Completed 10/20/2020 
fut_clim_files <- list.files('biomod2ing/bioclim_data/ccafs', pattern = '.*.asc', recursive = T, full.names = T)
landscapes <- c('nynp','synp','wynp','sgye') # ,'nynp','synp','wynp','sgye'
bioclim_vars <- seq(1,19,1)


walk(landscapes, function(ls){
  
  # Read iLand environment
  iland_env <- raster(paste0('../stand_grids/', ls,'_landscape_env_grid.asc'),
                      crs = '+init=EPSG:26912')
  
  # Extract points as sp and just RIDs
  iland_pts <- rasterToPoints(iland_env, spatial = T)
  iland_rids <- rasterToPoints(iland_env)[,3]
  
  
  # Transform iLand env to WGS for cropping
  iland_wgs <- projectRaster(iland_env, crs = '+init=EPSG:4326')
  
  walk(fut_clim_files, function(a_file){
    
    # A gross nest of splits to get the pertinent info...
    info1 <- str_split(str_split(a_file, 'biomod2ing/bioclim_data/ccafs/')[[1]], '_bio/')[[2]][1]
    info2 <- str_split(str_split(str_split(a_file, 'biomod2ing/bioclim_data/ccafs/')[[1]], '_bio/')[[2]][2],'.asc')[[1]][1]
    
    if( !info2 %in% paste0('bio_',bioclim_vars)){
      return(NULL)
    } else {
      print(paste('Running', ls, info1, info2))
    }
    
    # Read in, crop to landscape, transform to UTM, extract data to iLand stands
    ccasf <- raster(a_file, crs = '+init=EPSG:4326') %>%
      crop(., iland_wgs, snap = 'near') %>%
      projectRaster(., crs = '+init=EPSG:26912') %>%
      raster::extract(., iland_pts)
    
    # Reclassify iLand grid with CCASF climate data
    ccasf_iland <- reclassify(iland_env, cbind(iland_rids, ccasf))
    # Write out this raster
    writeRaster(ccasf_iland,
                filename = paste0('biomod2ing/bioclim_data/ccafs_iland/',ls,'_',info1,'_',info2,'.tif'),
                overwrite = T)
    
  })
})
