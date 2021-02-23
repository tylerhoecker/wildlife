library(raster)
library(tmap)
library(tidyverse)
library(sf)
library(stars)

species <- c('bbwo','mart','tahu')
periods <- c('hist',2030,2050,2070)
landscapes = c('nynp','wynp','synp','grte','sgye')

# Labels
gcm_labs = c("canesm2" = "Wet","hadgem2es" = "Dry")
spp_labs <- c('bbwo' = 'Black-backed Woodpecker', 'mart' = 'North American marten', 'tahu' = 'Red squirrel')
rcp_labs <- c('45' = 'Warm', '85' = 'Hot')

# Data sources
# Model projections
combined_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/combined_spatial/'
habitat_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/habitat_results/decadal_means/'
climate_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/biomod2ing/future_rasters/'

# Landscapes polygons, used for cropping and in the maps
ls_polygons <- 
  map(landscapes, ~ 
        st_read(paste0('../stand_grids/',.x,'_polygon/',.x,'_polygon.shp'))) %>% 
  set_names(landscapes) 
# Just the outter perimeter of landscapes, no holes.
ls_perimeters <- map(ls_polygons, ~ nngeo::st_remove_holes(.x)) 

ls_bboxes <- map(ls_polygons, ~ 
                   st_buffer(.x, 1000) %>% 
                   st_bbox())

# An extent object around ALL landscapes that can be used for cropping... 
# use raster::extent instead of st_bbox for compatability with raster
crop_box <-
  ls_polygons %>% 
  bind_rows() %>% 
  st_buffer(1000) %>% 
  extent() 

# DEM layer for hillshade
dem_background <- raster('/Users/sprague/Box/Work/PhD/GIS/GYE/DEM/DEM_hillshade.tif')
# Make it smaller, size of the GYE
dem_crop <- crop(dem_background, crop_box) 
# Make it just the size of the landscapes BB
dem_ls <- map(ls_perimeters, function(x){
  dem_crop <- crop(dem_background, x)
  #dem_mask  <- mask(dem_crop, x)
})

# Maxent thresholds 
cutoffs <- readRDS('biomod2ing/cutoff_values.Rds')

# Now the functions
# TESTING
# spp = 'mart'
# gcm = 'hadgem2es'
# rcp = '45'

# Outter-most wrapper function, iterates over each species  
mapping_sp <- function(spp){
  # Threshold
  thresh <- cutoffs[[spp]]$mean
  
  # Mapping function 
  map_fn <- function(ls,hab,clim,comb){
    tm_shape(clim,
             bbox = ls_bboxes[[ls]]) +
      tm_raster(breaks = c(0,thresh,1000),
                palette = c("#9b9d99",'#35978f'),
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
      tm_shape(ls_polygons[[ls]]) +
      tm_borders(col = 'black', lwd = 0.5) +
      tm_shape(ls_perimeters[[ls]]) +
      tm_borders(col = 'black', lwd = 1.5) 
  }

  mapping_scenario <- function(gcm, rcp){
    
    # First make historical column, since paths are a bit different
    # Historical habitat
    hist_hab_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/habitat_results/decadal_means/'
    hab_hist <- map(landscapes, ~  raster(paste0(hist_hab_dir,.x,'_',gcm,'_',rcp,'_hist_',spp,'.tif'))) 
    names(hab_hist) <- landscapes
    # Historical climate
    hist_clim_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/biomod2ing/historical_rasters/'
    clim_hist <- map(landscapes, ~  raster(paste0(hist_clim_dir,.x,'_',spp,'.tif')))
    names(clim_hist) <- landscapes
    # Historical combined
    comb_hist <- map(landscapes, ~ raster(paste0(combined_dir,.x,'_',gcm,'_',rcp,'_hist_',spp,'.tif')))
    names(comb_hist) <- landscapes
    
    maps_hist <- pmap(list(ls = landscapes,
                             hab = hab_hist,
                             clim = clim_hist,
                             comb = comb_hist), 
                        map_fn)
    maps_hist_col <- tmap_arrange(maps_hist, ncol = 1, asp = 1, outer.margins = 0)
    tmap_save(tm = maps_hist_col, 
              filename = paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_hist_',spp,'.png'), 
              width = 2, height = 8, units = 'in',
              dpi = 600)
    
    # Now create maps for each future period
    walk(periods[-1], function(period){
      # - habitat
      habitat <- map(landscapes, ~ 
                        raster(paste0(
                          habitat_dir,.x,'_',gcm,'_',rcp,'_',period,'_',spp,'.tif'))) 
      names(habitat) <- landscapes
      # - climate
      climate <- map(landscapes, ~ 
                         stack(paste0(climate_dir,.x,'_',gcm,'_',rcp,'_',period,'.tif'), 
                               bands = which(species == spp)))
      names(climate) <- landscapes
      # - combined
      combined <- map(landscapes, ~ 
                         raster(paste0(
                           combined_dir,.x,'_',gcm,'_',rcp,'_',period,'_',spp,'.tif'))) 
      names(combined) <- landscapes
      # - maps
      maps_period <- pmap(list(ls = landscapes,
                               hab = habitat,
                               clim = climate,
                               comb = combined), 
                          map_fn)
      maps_period_col <- tmap_arrange(maps_period, ncol = 1, asp = 1, outer.margins = 0)
      tmap_save(tm = maps_period_col, 
                filename = paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_',period,'_',spp,'.png'), 
                width = 2, height = 8, units = 'in',
                dpi = 600)
    })
    
    library(grid)
    library(gridExtra)
    library(png)
    
    map_cols <- map(periods, ~ 
                      readPNG(paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_',.x,'_',spp,'.png')))
    map_grobs <- map(map_cols, ~
                       rasterGrob(.x,
                                  interpolate = T))
    
    # col_hist <- readPNG(paste0('maps/temp_pngs/map_col_hist_',spp,'.png'))
    # col_2030 <- readPNG(paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_2030_',spp,'.png'))
    # col_2080 <- readPNG(paste0('maps/temp_pngs/map_col_',gcm,'_',rcp,'_2080_',spp,'.png'))
    # 
    lay <- rbind(c(1,3,5,7),
                 c(2,4,6,8))
    
    all_grobs <- list(textGrob('Historical'), map_grobs[[1]],
                      textGrob('2030'), map_grobs[[2]],
                      textGrob('2050'), map_grobs[[3]],
                      textGrob('2070'), map_grobs[[4]])
    
    png(filename = paste0('maps/',gcm,'_',rcp,'_',spp,'.png'), 
        width = 8, height = 9, units = 'in', res = 600)
    grid.arrange(grobs = all_grobs,
                 ncol = 4, nrow = 2,
                 heights = unit(c(0.3,8),'in'),
                 #widths = unit(rep(1.0,4),'in'),
                 layout_matrix = lay,
                 top = textGrob(paste0(gcm_labs[[gcm]],'-',rcp_labs[[rcp]],'  ',spp_labs[[spp]]), gp=gpar(fontsize=10)),
                 padding = unit(0, 'mm'))
    dev.off()
  }
  scenarios <- 
    expand.grid('gcm'= c('canesm2','hadgem2es'),
                'rcp' = c('45','85'))
  pwalk(scenarios, mapping_scenario)
}

walk(species, mapping_sp)



