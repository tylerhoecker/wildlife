library(raster)
library(landscapemetrics)
library(landscapetools)
library(tidyverse)

raster_dir <- '/Volumes/thexternal/iLand_gye/wildlife_pc/combined_spatial/'
species <- c('bbwo','mart','tahu')
periods <- c('hist',2030,2050,2070)
simulations <- seq(0,19)

# Maxent thresholds 
cutoffs <- readRDS('biomod2ing/cutoff_values.Rds')

# List of patch metrics
met_list <- list('area','enn','shape')


# TESTING
# landscape = c('grte')
# gcm = c('hadgem2es')
# rcp = c('85')


wildlife_ls_metrics <- function(landscape, gcm, rcp){
  
  # Names of combined rasters
  spatial_files <- list.files(
    path = 'combined_spatial',
    pattern = paste0(landscape,'_',gcm,'_',rcp,'.*\\.tif$'), full.names = T)
 
  decade_files <- map(periods, function(x){
    spatial_files[grepl(paste0('.*',x,'.*'), spatial_files)]}) 
  
  # Read in the files, do some things
  spatial_periods <- map(decade_files, function(x){
    map(species, function (sp){
      r <- raster(x[grepl(paste0('.*',sp,'.*'), x)]) %>%
        reclassify(., matrix(c(0,cutoffs[[sp]]$mean, 0, cutoffs[[sp]]$mean,1000, 1), ncol=3, byrow=TRUE)) 
      names(r) <- sp
      return(r) }) %>% 
      set_names(species) }) %>% 
    set_names(periods)
  
  # ------------------------------------------------------------------------------
  # Calculate metrics... these are patch-level metrics that are then summarized to class-level (just habitat patches)
  # ------------------------------------------------------------------------------
  # Function to calculate an area-weighted patch metric
  area_weighted <- function(metric, surface) {
    
    surface <- stack(surface) 
    metric_fn  <- get(paste0('lsm_p_', metric))
    metric_patch <- metric_fn(surface)
    area_patch <- lsm_p_area(surface) 
    
    # Weight by area
    result <- 
      left_join(x = metric_patch, y = area_patch, 
                by = c("layer", "level", "class", "id")) %>%
      mutate(value.w = value.x * value.y) %>%
      group_by(layer, class) %>%
      summarise(value = sum(value.w, na.rm = T) / sum(value.y, na.rm = T), .groups = 'drop') %>% 
      filter(class == 1) %>% select(-class) %>% 
      mutate(metric = paste0(metric,'_am'),
             species = case_when(layer == 1 ~ species[1],
                                 layer == 2 ~ species[2],
                                 layer == 3 ~ species[3])) %>% 
      select(species, metric, value)
    
    return(result)
  }
  
  patch_met_df <- 
    map_df(spatial_periods, function(x){
      map_df(met_list, area_weighted, surface = x) %>%
      mutate(landscape = landscape,
             gcm = gcm,
             rcp = rcp)}, .id = 'year')
}

# Build dataframe of scenarios to analyze
scenarios <- 
  expand.grid('landscape' = c('grte','synp','nynp','wynp','sgye'), #,'grte','nynp','wynp','synp','sgye'
              'gcm'= c('canesm2','hadgem2es'), # ,'hadgem2es','canesm2'
              'rcp' = c('45','85')) # ,2050,2100 

wildlife_metrics_df <- 
  pmap_df(scenarios, wildlife_ls_metrics) %>% 
  mutate(year = factor(year))

# Fill in NAs
complete_df <- expand.grid('landscape' = c('grte','synp','nynp','wynp','sgye'), #,'grte','nynp','wynp','synp','sgye'
                           'gcm'= c('canesm2','hadgem2es'), # ,'hadgem2es','canesm2'
                           'rcp' = c('45','85'), 
                           'year'  = factor(periods, levels = periods),
                           'species' = species,
                           'metric' = paste0(unlist(met_list),'_am')) %>% 
  full_join(wildlife_metrics_df) %>% 
  group_by(gcm, rcp, year, species, metric) %>% 
  summarise(value = mean(value, na.rm = T))

# ------------------------------------------------------------------------------
# Plotting
plot_df <- complete_df %>% filter(gcm == 'hadgem2es' & rcp == '85' | gcm == 'canesm2' & rcp == '45') %>%
  mutate(value = ifelse(value > 10000, 10000, value))

fig_cols <- c('#f16913','#7f2704','#41ab5d','#00441b','#807dba','#3f007d')

gcm_labs <- c('canesm2' = 'Wet (CanESM2)', 'hadgem2es' = 'Dry (HadGEM2-ES)')
rcp_labs <- c('45' = 'Warm-Wet', '85' = 'Hot-Dry')
metric_labs <- c('area_am' = 'Patch size (ha)', 
                 'enn_am' = 'NN-Distance (m)', 
                 'shape_am' = 'Shape index')

ggplot(plot_df) +
  geom_col(aes(x = year, y = value, fill = species), 
           position = 'dodge2', color = 'black') +
  scale_fill_manual('Species', values = fig_cols[c(5,3,1)],
                    labels = c('Black-backed Woodpecker','North American marten', 'Red squirrel')) +
  facet_grid(metric~rcp, 
             labeller = labeller(rcp = rcp_labs, metric = metric_labs), 
             scales = 'free') + # landscape~gcm+
  scale_x_discrete(labels = c('Hist.','2030','2050','2070')) +
  labs(y = 'Landscape metric value') +
  theme_bw() %+replace%
  theme(
    # axes
    axis.text = element_text(color="black",size=6),
    axis.title = element_text(size=8),
    axis.title.x = element_blank(),
    axis.ticks = element_line(color="black",size=0.25),
    # legend
    legend.position ="bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.height = unit(0.25, 'cm'),
    legend.key.width = unit(0.25, 'cm'),
    # facets
    strip.background =  element_blank(),
    strip.text = element_text(size = 8))  

ggsave(filename = 'landscape_metrics.png',
       height = 120, width = 100, units = 'mm',
       dpi = 600)



















# Previous version that is based just on the habitat rasters
# ------------------------------------------------------------------------------
# Set the directory where the rasters are stored
raster_dir <- "/Volumes/thexternal/iLand_gye/wildlife_pc/habitat_results/decadal/"

# Things to specify
species <- c('bbwo','mart','tahu')
snaps <- c(2020,2050,2100)
simulations <- seq(0,5)

# Build dataframe of scenarios to analyze
scenarios <- 
  expand.grid('landscape' = c('grte','synp','nynp'), #,'grte','nynp','wynp','synp','sgye'
             'gcm'= c('canesm2','hadgem2es'), # ,'hadgem2es','canesm2'
             'rcp' = c('45','85')) # ,2050,2100 

# Set important parameters for landscape metric analysis
# Threshold of proportion of simulations it was habitat
hab_threshold <- 0.7
# List of patch metrics
met_list <- list('area','enn','shape')

# landscape = c('grte')
# gcm= c('hadgem2es')
# rcp = c('85')
# year = 202

wildlife_ls_metrics <- function(landscape, gcm, rcp){
  
  habitat_files <- list.files(
    path = 'habitat_results/decadal',
    pattern = c(paste0('habitat_',landscape,'_',gcm,'_',rcp), '.*'), full.names = T) 
  
  # Select on files for certain simulations (not neccesary after all simulations are run)
  these_files <- paste0(landscape,'_',gcm,'_',rcp,'_[',simulations[1],'-',simulations[length(simulations)],']_.*') 
  
  habitat_files <- habitat_files[grepl(these_files, habitat_files)]
  
  decade_files <- map(snaps, function(x){
    habitat_files[grepl(paste0('.*',x,'.*'), habitat_files)]}) 
  
  habitat_snaps <- map(decade_files, function(x){
    map(species, function (sp){
      r <- mean(stack(x[grepl(paste0('.*_[0-9]{4}_',sp,'.*'), x)])) %>% 
        reclassify(., matrix(c(0,hab_threshold, 0, hab_threshold,1, 1), ncol=3, byrow=TRUE)) 
      names(r) <- sp
      return(r) }) %>% 
      set_names(species) }) %>% 
    set_names(snaps)
  
  
  combined_files <- list.files(
    path = 'habitat_results/decadal',
    pattern = c(paste0('combined_',landscape,'_',gcm,'_',rcp), '.*'), full.names = T) 

  # ------------------------------------------------------------------------------
  # Calculate metrics... these are patch-level metrics that are then summarized to class-level (just habitat patches)
  # ------------------------------------------------------------------------------
  # Function to calculate an area-weighted patch metric
  # 'mart' and 'tahu' are currently hard-coded here...
  area_weighted <- function(metric, surface) {
    
    surface <- stack(surface) 
    
    metric_fn  <- get(paste0('lsm_p_', metric))
    metric_patch <- metric_fn(surface)
    area_patch <- lsm_p_area(surface) 
    
    # Weight by area
    result <- 
      left_join(x = metric_patch, y = area_patch, 
                by = c("layer", "level", "class", "id")) %>%
      mutate(value.w = value.x * value.y) %>%
      group_by(layer, class) %>%
      summarise(value = sum(value.w) / sum(value.y), .groups = 'drop') %>% 
      filter(class == 1) %>% select(-class) %>% 
      mutate(metric = paste0(metric,'_am'),
             species = case_when(layer == 1 ~ species[1],
                                 layer == 2 ~ species[2],
                                 layer == 3 ~ species[3])) %>% 
      select(species, metric, value)
    
    return(result)
  }
  
  patch_met_df <- map_df(habitat_snaps, function(x){
    map_df(met_list, area_weighted, surface = x) %>%
      mutate(landscape = landscape,
             gcm = gcm,
             rcp = rcp,
             year = year)}, .id = 'year')
  
  
}


wildlife_metrics_df <- 
  pmap_df(scenarios, wildlife_ls_metrics) %>% 
  mutate(year = factor(year))

# Fill in NAs
complete_df <- expand.grid('landscape' = c('grte','synp','nynp'), #,'grte','nynp','wynp','synp','sgye'
                        'gcm'= c('canesm2','hadgem2es'), # ,'hadgem2es','canesm2'
                        'rcp' = c('45','85'), 
                        'year'  = factor(snaps),
                        'species' = species,
                        'metric' = paste0(unlist(met_list),'_am')) %>% 
  full_join(wildlife_metrics_df) %>% 
  group_by(gcm, rcp, year, species, metric) %>% 
  summarise(value = mean(value, na.rm = T))
                      
plot_df <- complete_df %>% filter(rcp == '85') %>% 
  mutate(value = ifelse(value > 700, 700, value))


fig_cols <- c('#f16913','#7f2704','#41ab5d','#00441b','#a6a3ce','#7570b3')
gcm_labs <- c('Warm-wet (CanESM2)', 'Warm-dry (HadGEM2-ES)')
names(gcm_labs) <- unique(plot_df$gcm)

b_plot <- 
ggplot(plot_df) +
  geom_col(aes(x = year, y = value, fill = species), position = 'dodge2') +
  scale_fill_manual('Species', values = fig_cols[c(2,4,6)]) +
  facet_grid(metric~gcm, , labeller = labeller(gcm = gcm_labs), scales = 'free') + # landscape~gcm+
  labs(x = 'Year', y = 'Landscape metric value') +
  theme_bw(base_size = 16) +
  theme(strip.background =  element_rect('transparent'),
        strip.text = element_text(face = 'bold', size = 16))
b_plot

ggsave(b_plot, filename = 'habitat_landscape_metrics.png', dpi = 300, height = 8, width = 12)



# TO DELETE
# # Read in rasters
# all_files <- list.files(raster_dir, full.names = T)
# habitat_pattern <- paste0('habitat_results/decadal/habitat_',
#                           scenario_name,'_',sim,'_',model_year+2005,'_',x,'.tif')
# 
# paste0('.*_',landscape,'_',gcm,'_rcp85_yr',year-2005)
# files <- all_files[grepl(pattern, all_files)]
# habitats <- map(files, raster) 
# 
# # Create raster stack, removing BBWO for now... since threshold will be very different (and maybe irrelevant)
# cont_hab <- habitats[-1] %>% 
#   stack() %>% 
#   `names<-` (c('mart','tahu'))
# # Take a look
# show_landscape(cont_hab)

# # Convert stack to binary: 1 = habitat, 0 = not
# binary_hab <- cont_hab > hab_threshold
# Take a look
# show_landscape(habitat_snaps[[1]], discrete = T)



