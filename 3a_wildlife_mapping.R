library(raster)
library(tmap)
library(tidyverse)
library(sf)
library(stars)

species <- c('bbwo','mart','tahu')
snaps <- c('historical',2050,2080)
landscapes = c('nynp','wynp','synp','grte','sgye')

# Labels
gcm_labs = c("canesm2" = "CanESM2","hadgem2es" = "HadGEM2-ES")
spp_labs <- c('bbwo' = 'P. arcticus', 'mart' = 'M. americana', 'tahu' = 'T. hudsonicus')

# Data sources
# Model projections
combined_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/combined_spatial/'
habitat_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/habitat_results/snapshots/'
climate_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/biomod2ing/future_rasters/'

# Landscapes
ls_outlines <- map(landscapes, ~ rgdal::readOGR(paste0('../stand_grids/',.x,'_polygon/',.x,'_polygon.shp')))
names(ls_outlines) <- landscapes

ls_outside <- map(landscapes, function(x){
  read_sf(paste0('/Volumes/thexternal/iLand_gye/stand_grids/',x,'_polygon/',x,'_polygon.shp')) %>% 
    nngeo::st_remove_holes()}) 
names(ls_outside) <- landscapes

# Bounding box for landscape maps
ls_bboxes <- map(ls_outside, st_bbox) 


crop_box <- extent(491277.8, 569877.8, 4742803.9,5010703.9)
     

dem_background <- raster('/Users/sprague/Box/Work/PhD/GIS/GYE/DEM/DEM_hillshade.tif')

dem_crop <- crop(dem_background, crop_box) 

dem_ls <- map(ls_outside, function(x){
  dem_crop <- crop(dem_background, x)
  #dem_mask  <- mask(dem_crop, x)
})

# Maxent thresholds 
cutoffs <- readRDS('biomod2ing/cutoff_values.Rds')
names(cutoffs) <- species



# TESTING
# spp = 'mart'
# gcm = 'hadgem2es'
# rcp = '45'
# ls = landscapes[[1]]
# 
mapping_sp <- function(spp){
  
  # Historical habitat
  hist_hab_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/habitat_results/historical/'
  hab_hist <- map(landscapes, ~  raster(paste0(hist_hab_dir,.x,'_',spp,'.tif'))) 
  names(hab_hist) <- landscapes
  # Historical climate
  hist_clim_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/biomod2ing/historical_rasters/'
  clim_hist <- map(landscapes, ~  raster(paste0(hist_clim_dir,.x,'_',spp,'.tif')))
  names(clim_hist) <- landscapes
  # Historical combined
  comb_hist <- map(landscapes, ~ raster(paste0(combined_dir,.x,'_historical_',spp,'.tif')))
  names(comb_hist) <- landscapes
  
  # Threshold
  thresh <- cutoffs[[spp]]
  
  # Mapping function 
  map_fn <- function(ls,hab,clim,comb){
    tm_shape(clim,
             bbox = ls_bboxes[[ls]]) +
      tm_raster(breaks = c(0,thresh,1000),
                palette = c("#525252",'#35978f'),
                legend.show = FALSE) +
      tm_shape(hab,
               bbox = ls_bboxes[[ls]]) +
      tm_raster(breaks = c(0,0.5,1),
                palette = c("#00000000", '#bf812d'),
                legend.show = FALSE) +
      tm_shape(comb,
               bbox = ls_bboxes[[ls]]) +
      tm_raster(breaks = c(0,thresh,1000),
                palette = c("#00000000", '#762a83'),
                legend.show = FALSE) +
      
      tm_shape(dem_crop, raster.downsample = T) +
      tm_raster(palette = gray(0:10 / 10), style = "cont", 
                legend.show = FALSE,
                alpha = 0.25) +
      tm_shape(ls_outlines[[ls]]) +
      tm_borders(col = 'black', lwd = 0.5) +
      tm_shape(ls_outside[[ls]]) +
      tm_borders(col = 'black', lwd = 1.5) 
  }
  
  maps_hist <- pmap(list(ls = landscapes,
                         hab = hab_hist,
                         clim = clim_hist,
                         comb = comb_hist), 
                    map_fn)
  maps_hist_col <- tmap_arrange(maps_hist, ncol = 1, asp = 1)
  tmap_save(tm = maps_hist_col, 
            filename = paste0('maps/temp_pngs/map_col_hist_',spp,'.png'), 
            width = 2, height = 8, units = 'in',
            dpi = 500)
  
  mapping_scenario <- function(gcm, rcp){
    
    # 2050 maps
    # - habitat
    hab_2050 <- map(landscapes, ~ 
                      raster(paste0(
                        habitat_dir,.x,'_',gcm,'_',rcp,'_2050_',spp,'.tif'))) 
    names(hab_2050) <- landscapes
    # - climate
    clim_2050 <- map(landscapes, ~ 
                       stack(paste0(climate_dir,.x,'_',gcm,'_',rcp,'_2050.tif'), 
                             bands = which(species == spp)))
    names(clim_2050) <- landscapes
    # - combined
    comb_2050 <- map(landscapes, ~ 
                       raster(paste0(
                         combined_dir,.x,'_',gcm,'_',rcp,'_2050_',spp,'.tif'))) 
    names(comb_2050) <- landscapes
    # - maps
    maps_2050 <- pmap(list(ls = landscapes,
                           hab = hab_2050,
                           clim = clim_2050,
                           comb = comb_2050), 
                      map_fn)
    maps_2050_col <- tmap_arrange(maps_2050, ncol = 1, asp = 1)
    tmap_save(tm = maps_2050_col, 
              filename = paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_2050_',spp,'.png'), 
              width = 2, height = 8, units = 'in',
              dpi = 500)

    # - habitat 
    hab_2080 <- map(landscapes, ~ 
                      raster(paste0(
                        habitat_dir,.x,'_',gcm,'_',rcp,'_2080_',spp,'.tif'))) 
    names(hab_2080) <- landscapes
    # - climate
    clim_2080 <- map(landscapes, ~ 
                       stack(paste0(climate_dir,.x,'_',gcm,'_',rcp,'_2080.tif'), 
                             bands = which(species == spp)))
    names(clim_2050) <- landscapes
    # - comb
    comb_2080 <- map(landscapes, ~ 
                      raster(paste0(
                        combined_dir,.x,'_',gcm,'_',rcp,'_2080_',spp,'.tif'))) 
    names(comb_2080) <- landscapes
    # - maps
    maps_2080 <- pmap(list(ls = landscapes,
                           hab = hab_2080,
                           clim = clim_2080,
                           comb = comb_2080), 
                      map_fn)
    maps_2080_col <- tmap_arrange(maps_2080, ncol = 1, asp = 1)
    
    tmap_save(tm = maps_2080_col, 
              filename = paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_2080_',spp,'.png'),
              width = 2, height = 8, units = 'in',
              dpi = 500)
    
    library(grid)
    library(gridExtra)
    library(png)
    
    col_hist <- readPNG(paste0('maps/temp_pngs/map_col_hist_',spp,'.png'))
    col_2050 <- readPNG(paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_2050_',spp,'.png'))
    col_2080 <- readPNG(paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_2080_',spp,'.png'))
    
    lay <- rbind(c(1,3,5),
                 c(2,4,6))
    
    all_grobs <- list(textGrob('Historical'), rasterGrob(col_hist),
                   textGrob('2050'), rasterGrob(col_2050),
                   textGrob('2080'), rasterGrob(col_2080))
    
    png(filename = paste0('maps/',gcm,'_',rcp,'_',spp,'.png'), 
        width = 4, height = 6, units = 'in', res = 500)
    grid.arrange(grobs = all_grobs,
                 heights = c(0.05,1),
                 widths = c(1.2,1.2,1.2),
                 layout_matrix = lay,
                 top = paste0(gcm_labs[[gcm]],' RCP-',rcp,'  ',spp_labs[[spp]]))
    dev.off()
  }
  scenarios <- 
    expand.grid('gcm'= c('canesm2','hadgem2es'),
                'rcp' = c('45','85'))
  pwalk(scenarios, mapping_scenario)
}

walk(species, mapping_sp)



