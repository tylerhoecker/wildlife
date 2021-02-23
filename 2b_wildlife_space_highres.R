library(raster)
library(tidyverse)

# Things to specify ------------------------------------------------------------
species <- c('bbwo','mart','tahu')
# simulations <- c(0,1,5,10,15,19)
simulations <- seq(0,19)
# This is a bit funky, but these are repeated so that climate periods are associated with each decade of habitat
periods <- c(2030,2050,2070) # This will get revised later in the script to be c('historical','2030','2050','2070')
# Note here that 2010 is in the historical slot
dec_groups <- list(2010, c(2020,2030,2040),c(2040,2050,2060),c(2060,2070,2080)) 

# decades <- seq(2030,2100,10)
cutoffs <- readRDS('biomod2ing/cutoff_values.Rds')
# ------------------------------------------------------------------------------


landscape = c('grte')
gcm= c('canesm2')
rcp = c('45')

wildlife_space <- function(landscape,gcm,rcp){
  
  if(length(
    list.files(path = 'combined_spatial/',
               pattern = paste0(landscape,'_',gcm,'_',rcp,'.*'))
    ) == length(species)*length(dec_groups) ){
    print(paste(landscape,gcm,rcp,'already complete!', sep = ' '))
    return(NULL)
  } else {
    print(paste(landscape,gcm,rcp, sep = ' '))
  }
  if(landscape == 'sgye' & gcm == 'hadgem2es' & rcp == '85'){
    simulations <- c(0:13,15:19)
  }
  
  ## CLIMATE -------------------------------------------------------------------
  # List the results files for this scenario
  climate_files <- list.files(
    path = 'biomod2ing/future_rasters',
    pattern = paste0(landscape,'_',gcm,'_',rcp,'.*'), full.names = T) 
  
  # # Some weird stuff to sort this list of files by period, then species
  # these_files <- map(periods, function(x){
  #   climate_files[grepl(paste0('.*',x,'.*'), climate_files)]}) 
  
  # A list of future climate rasters now sorted by decade and then by species
  fut_clim <- map(climate_files, ~
                    stack(.x) %>%  
                    unstack(.) %>% 
                    set_names(., species)) %>% 
    set_names(periods)
  
  # Add historical climate rasters to this list
  hist_files <- list.files(
    path = 'biomod2ing/historical_rasters',
    pattern = paste0(landscape,'.*'), full.names = T)
  
  # Read in the raster for each species
  hist_clim <- map(hist_files, raster) %>% 
    set_names(species)
  # Make it a list named by decade, like the future list
  # hist_decade <- list('2010' = hist_clim, '2020' = hist_clim)
  hist_clim <- list('hist' = hist_clim)
  
  # Combine the historical and future into decade-element list
  period_climate <- c(hist_clim, fut_clim)
  
  # Now rename `periods` so that it includes these new historical decades
  periods <- names(period_climate)
  
  ## HABITAT  ------------------------------------------------------------------
  # List habitat files
  habitat_files <- list.files(
    path = 'habitat_results/decadal',
    pattern = c(paste0('habitat_',landscape,'_',gcm,'_',rcp), '.*'), full.names = T) 
  # Select only files for certain simulations (not necessary after all simulations are run)
  these_files <- paste0(landscape,'_',gcm,'_',rcp,'_[',paste(simulations, collapse = '|'),']+_.*')
  habitat_files <- habitat_files[grepl(these_files, habitat_files)]
  
  
  # Check that everything is there!
  if(length(habitat_files)/length(simulations) < 30){
    stop(paste0('Habitat results for ', paste0(landscape,'_',gcm,'_',rcp), ' missing: ', 
                length(habitat_files), ' of ', length(simulations) * 30,' present'))
  }
  
  # Figuring this out... must average all simulations first, then do mosaic within each climate period-group
  all_habitat <- habitat_files %>% 
    map(., raster) 
  
  # Organize the rasters by decades, since all sims for each decade will get averaged
  decade_hab <- map(seq(2010,2100,10), function(x){
    all_habitat[grepl(paste0('.*',x,'.*'), all_habitat)]}) %>% 
    set_names(seq(2010,2100,10))
  
  # For each species, stack all sims in each decade and calculate mean,
  # then reclassify the mean to binary
  mean_hab <- map(decade_hab, function(x){
    map(species, function (sp){
      x[grepl(paste0('.*_[0-9]{4}_',sp,'.*'), x)] %>% 
        stack() %>% 
        mean() %>% 
        reclassify(matrix(c(0,0.5,0, 0.5,1,1), ncol=3, byrow=TRUE))}) %>% 
      set_names(species)}) %>% 
    unlist()
  
  name_index <- names(mean_hab)
  
  period_hab <- 
    map(dec_groups, function(group){ 
       map(species, function(sp){
        period_hab <- 
          mean_hab[grepl(paste0(paste(group, collapse = paste0('.',sp,'|')),'.',sp),
                         name_index)]
        # For climate periods that include more than 1 decade (all but historical [==2010]),
        # Calculate the mean
        period_hab <- mean(stack(period_hab))
        # Then reclassify the mean (thus, 2/3 = 1 and 1/3 = 0)
        period_hab <- reclassify(period_hab, matrix(c(0,0.5,0, 0.5,1,1), ncol=3, byrow=TRUE))
        return(period_hab)
      }) %>% 
        set_names(species)
    }) %>% 
    set_names(periods)
  
  # habitat = period_hab[[2]]
  # climate = period_climate[[2]]
  # period = periods[[2]]
  pwalk(list(period_hab, period_climate, periods), 
        function(habitat, climate, period){
          combined <- unstack(stack(habitat) * stack(climate))
          to_write <- list(habitat, combined)
          walk2(to_write, c('habitat_results/decadal_means/','combined_spatial/'),
                function(x,write_path){
                  writeRaster(stack(x),
                              filename = paste0(write_path,landscape,'_',gcm,'_',rcp,'_',period,'.tif'),
                              bylayer = T, suffix = species,
                              overwrite = T)})
          })
}

# This is a dataframe listing all landscape-GCM-RCP-iteration combinations for function to iterate over
scenarios <- 
  expand.grid('landscape' = c('grte','nynp','synp','wynp','sgye'), # 'grte','nynp','nynp','synp','wynp','sgye'
              'gcm'= c('canesm2','hadgem2es'), # 'canesm2', 'hadgem2cc','hadgem2es',,'hadgem2cc',,'canesm2'
              'rcp' = c('45','85'))

pwalk(scenarios, wildlife_space)









# OLD ---------
# TESTING

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
