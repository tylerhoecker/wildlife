library(raster)
library(tidyverse)
#library(furrr)
library(tmap)

# Things to specify
species <- c('bbwo','mart','tahu')
snaps <- c(2050,2080)
#simulations <- seq(0,19)


# TESTING
# landscape = c('grte')
# gcm= c('canesm2')
# rcp = c('45')

# Function to read in the results of the habitat (1a) and climate (1b) envelope modeling scripts
wildlife_space <- function(landscape,gcm,rcp){
  
  print(paste(landscape,gcm,rcp, sep = ' '))
  ## CLIMATE
  # List the results files for this scenario
  climate_files <- list.files(
    path = 'biomod2ing/future_rasters',
    pattern = paste0(landscape,'_',gcm,'_',rcp,'.*'), full.names = T) 
  
  these_files <- map(snaps, function(x){
    climate_files[grepl(paste0('.*',x,'.*'), climate_files)]}) 
  
  if(length(these_files) > 1){
    climate_snaps <- map(these_files, ~ stack(.x))
  } else {
    climate_snaps <- stack(these_files[[1]])
  }
  
  ## HABITAT
  habitat_files <- list.files(
    path = 'habitat_results/decadal',
    pattern = c(paste0('habitat_',landscape,'_',gcm,'_',rcp), '.*'), full.names = T) 
  
  # Select on files for certain simulations (not neccesary after all simulations are run)
  # these_files <- paste0(landscape,'_',gcm,'_',rcp,'_[',simulations[1],'-',simulations[length(simulations)],']_.*') 
  # habitat_files <- habitat_files[grepl(these_files, habitat_files)]
  
  decade_files <- map(snaps, function(x){
    habitat_files[grepl(paste0('.*',x,'.*'), habitat_files)]}) 
  
  habitat_snaps <- map(decade_files, function(x){
    map(species, function (sp){
      x[grepl(paste0('.*_[0-9]{4}_',sp,'.*'), x)] %>% 
        stack() %>% 
        mean() %>% 
        # Quicker than reclassify - just round! Not very consequential here since all scenarios are very similar at this point
        round()
     })%>% 
      set_names(species) %>% 
      stack()}) %>% 
    set_names(snaps)
  
  # LEAVING OFF HERE NEED TO MAKE BOTH A STACK OR SOMETING
  combined_snaps <- 
    map2(habitat_snaps, climate_snaps, ~ .x*.y)
  
  walk2(combined_snaps, snaps, ~
         writeRaster(.x,
                     filename = paste0('combined_spatial/',landscape,'_',gcm,'_',rcp,'_',.y,'.tif'),
                     bylayer = T, suffix = species,
                     overwrite = T)
       )
  
  
  # habitat = habitat_snaps[[1]]
  # climate = climate_snaps[[1]]
  # snap_year = snaps[[1]]
  # Can produce PNG maps... ----------------------------------------------------
  # iland_poly <- rgdal::readOGR(paste0('../stand_grids/',landscape,'_polygon/',landscape,'_polygon.shp'))
  # map_niches <- function(snap_yr,habitat,climate){
  #   
  #   combined <- stack(habitat) * stack(climate)
  #   combined_l <- unstack(combined) %>% 
  #     set_names(species)
  #   
  #   # Plot the maps, save as pdf
  #   col_breaks = seq(0, 1, by = 0.1)
  #   
  #   #layer = habitat[[1]]
  #   #plot_title = 'BBWO'
  #   #col_pal <-'Oranges'
  #   
  #   map_fn <- function(layer, plot_title, col_pal){
  #     map <- 
  #       tm_shape(layer) +
  #       tm_raster(title = 'Prob. habitat',
  #                 breaks = col_breaks,
  #                 palette = col_pal,
  #                 midpoint = 0.50) +
  #       tm_compass(type = "arrow", position = c("left")) +
  #       tm_scale_bar(breaks = c(0, 5, 10), text.size = 0.75) +
  #       tm_layout(main.title =
  #                   toupper(paste(plot_title,landscape,gcm,rcp,snap_yr)),
  #                 main.title.size = 1,
  #                 frame = T, frame.lwd = 3,
  #                 legend.outside = T,
  #                 # legend.position = c('right','top'),
  #                 inner.margins = c(0.02,0.02,0.02,0.02),
  #                 bg.color = "grey85") +
  #       tm_legend(show=TRUE) +
  #       tm_shape(iland_poly) +
  #       tm_borders()
  #     
  #     tmap_save(map, file = paste0('habitat_maps/',plot_title,'_',landscape,'_',gcm,'_',rcp,'_',snap_yr,'.png'),
  #               width = 7, height = 7, units = 'in', dpi = 300)
  #   }
  #   
  #   habitat_maps <- list(layer = habitat,
  #                        plot_title = paste0('habitat-',species),
  #                        col_pal = c('Oranges','Greens','Purples'))
  #   pwalk(habitat_maps, map_fn)
  #   
  #   climate_maps <- list(layer = climate,
  #                        plot_title = paste0('climate-',species),
  #                        col_pal = c('Oranges','Greens','Purples'))
  #   pwalk(climate_maps, map_fn)
  #   
  #   combined_maps <- list(layer = combined_l,
  #                         plot_title = paste0('combined-',species),
  #                         col_pal = c('Oranges','Greens','Purples'))
  #   pwalk(combined_maps, map_fn)
  #   
  # }
  # 
  # to_map <- list(snap_yr = snaps,
  #                habitat = habitat_snaps,
  #                climate = climate_snaps)
  # 
  # pwalk(to_map, map_niches)
  # ----------------------------------------------------------------------------
}

# This is a dataframe listing all landscape-GCM-RCP-iteration combinations for function to iterate over
scenarios <- 
  expand.grid('landscape' = c('grte','nynp','synp','wynp','sgye'), # 'grte','nynp','nynp','synp','wynp','sgye'
              'gcm'= c('canesm2','hadgem2es'), # 'canesm2', 'hadgem2cc','hadgem2es',,'hadgem2cc',,'canesm2'
              'rcp' = c('45','85'))

pwalk(scenarios, wildlife_space)














# OLD?
# Produce maps by decade ---------------------------------------------------
# Plot the maps, save as pdf
col_breaks = seq(0, 1, by = 0.1)

#layer = habitat[[1]]
#plot_title = 'BBWO'
#col_pal <-'Oranges'

map_fn <- function(layer, plot_title, col_pal){
  map <- 
    tm_shape(layer) +
    tm_raster(title = 'Prob. habitat',
              breaks = col_breaks,
              palette = col_pal,
              midpoint = 0.50) +
    tm_compass(type = "arrow", position = c("left")) +
    tm_scale_bar(breaks = c(0, 5, 10), text.size = 0.75) +
    tm_layout(main.title =
                toupper(paste(plot_title,landscape,gcm,rcp,decade_name)),
              main.title.size = 1,
              frame = T, frame.lwd = 3,
              legend.outside = T,
              # legend.position = c('right','top'),
              inner.margins = c(0.02,0.02,0.02,0.02),
              bg.color = "grey85") +
    tm_legend(show=TRUE) +
    tm_shape(iland_poly) +
    tm_borders()
  
  tmap_save(map, file = paste0('habitat_maps/',plot_title,'_',landscape,'_',gcm,'_',rcp,'_',decade_name,'.png'),
            width = 7, height = 7, units = 'in', dpi = 300)
}

habitat_maps <- list(layer = habitat,
                     plot_title = paste0('habitat-',species),
                     col_pal = c('Oranges','Greens','Purples'))
pwalk(habitat_maps, map_fn)

climate_maps <- list(layer = climate,
                     plot_title = paste0('climate-',species),
                     col_pal = c('Oranges','Greens','Purples'))
pwalk(climate_maps, map_fn)

combined_maps <- list(layer = combined_l,
                      plot_title = paste0('combined-',species),
                      col_pal = c('Oranges','Greens','Purples'))
pwalk(combined_maps, map_fn)

}

to_map <- list(decade_name = decades,
               habitat = habitat_decades,
               climate = climate_decades)

pwalk(to_map, map_niches)
