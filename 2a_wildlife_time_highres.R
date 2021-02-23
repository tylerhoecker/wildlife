library(raster)
library(tidyverse)

setwd("/Volumes/thexternal/iLand_gye/wildlife_pc")

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

# TESTING
landscape = c('grte')
gcm= c('hadgem2es')
rcp = c('85')
 
# Function to read in the results of the habitat (1a) and climate (1b) envelope modeling scripts
wildlife_combined <- function(landscape,gcm,rcp){
  
  if(file.exists(paste0('combined_results_csv/decadal/',paste(landscape,gcm,rcp, sep = '_'), '.csv'))){
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
    pattern = paste0(landscape,'_',gcm,'_',rcp,'.*\\.tif$'), full.names = T) 
  
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
    pattern = paste0(landscape,'.*\\.tif$'), full.names = T)
  
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
  # ----------------------------------------------------------------------------
  
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
  
  # Weird stuff to sort the lists by simulations > decades > species
  all_habitat <- habitat_files %>% 
    map(., raster) 
  
  sim_habitat <- 
    map(simulations, function(x){
      all_habitat[grepl(paste0('.*_',x,'_[0-9]{4}.*'), habitat_files)]}) 
  
  sim_list_index <- 
    map(sim_habitat, function(x){ map_chr(x, names)})[[1]]
  
  sim_period_habitat <- 
    map(sim_habitat, function(sim){ 
      map(dec_groups, function(group){ 
        map(species, function(sp){
          period_hab <- 
            sim[grepl(paste0('[0-9]+_',paste(group, collapse = paste0('_',sp,'|')),'_',sp), 
                      sim_list_index)]
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
    }) 
  # ----------------------------------------------------------------------------
  
  # habitat = sim_period_habitat[[1]][[2]]
  # climate = period_climate[[2]]
  
  calc_niche <- function(habitat, climate, period, sim_num){
    
    # First, save the continuous version of the joint SDM for other uses
    combined_cont <- stack(habitat) * stack(climate)
    writeRaster(combined_cont,
                filename = paste0('combined_spatial/per_sim/', 
                                  paste(landscape,gcm,rcp,sim_num,period, sep = '_'),
                                  '.tif'),
                bylayer = T, suffix = species,
                overwrite = T)
    
    # Convert to binary
    climate_med <- 
      map2(climate, cutoffs, ~ 
             reclassify(.x, matrix(c(0,.y[['mean']],0, .y[['mean']],1000,1), ncol=3, byrow=TRUE)))
    climate_low <- 
      map2(climate, cutoffs, ~ 
             reclassify(.x, matrix(c(0,.y[['lower']],0, .y[['lower']],1000,1), ncol=3, byrow=TRUE)))
    climate_high <- 
      map2(climate, cutoffs, ~ 
             reclassify(.x, matrix(c(0,.y[['upper']],0, .y[['upper']],1000,1), ncol=3, byrow=TRUE)))
    
    
    # Multiply the rasters together (easily done as a stack, then unstack)
    climates <- list(climate_med, climate_low, climate_high)
    
    combined <- map(climates, ~ unstack(stack(habitat) * stack(.))) 
    
    # combined <- stack(habitat) * stack(climate)
    # combined_l <- unstack(combined)
    # names(combined_l) <- species

    # Use raster::freq to count suitable/unsuitable
    result <- data.frame(
      'combined_med' = map_dbl(combined[[1]], freq, value = 1),
      'combined_low' = map_dbl(combined[[2]], freq, value = 1),
      'combined_high' = map_dbl(combined[[3]], freq, value = 1),
      # 'comb_unsuitable' = map_dbl(combined_l, freq, value = 0),
      'habitat' = map_dbl(habitat, freq, value = 1),
      # 'hab_unsuitable' = map_dbl(habitat, freq, value = 0),
      'climate_med' = map_dbl(climate_med, freq, value = 1),
      'climate_low' = map_dbl(climate_low, freq, value = 1),
      'climate_high' = map_dbl(climate_high, freq, value = 1)) %>% 
      # 'clim_unsuitable'= map_dbl(climate, freq, value = 0)) %>% 
      mutate(species = rownames(.))
    
    return(result)
  }
  
  niche_df <- 
    map2_df(sim_period_habitat, simulations, function(period_habitat, sim_num){
      pmap_df(list(period_habitat, period_climate, periods, sim_num), calc_niche, .id = 'decade')}, 
      .id = 'simulation') %>% 
    mutate(landscape = landscape,
           gcm = gcm,
           rcp = rcp) %>% 
    write_csv(., paste0('combined_results_csv/decadal/',paste(landscape,gcm,rcp, sep = '_'), '.csv'))
  
}

scenarios <- 
  expand.grid('landscape' = c('grte','nynp','synp','wynp','sgye'), # 'grte','nynp','nynp','synp','wynp','sgye'
              'gcm'= c('canesm2','hadgem2es'), # 'canesm2', 'hadgem2cc','hadgem2es',,'hadgem2cc',,'canesm2'
              'rcp' = c('45','85'))

pwalk(scenarios, wildlife_combined)


# ------------------------------------------------------------------------------
# Plotting ---------------------------------------------------------------------

# Read in all datas
full_wildlife_df <- list.files('combined_results_csv/decadal', full.names = T) %>% 
  map_df(., read_csv) %>% 
  mutate(decade = factor(decade, 
                         levels = c('hist','2030','2050','2070'),
                         labels = c('Hist.','2030s','2050s','2070s')),
         decade_num = as.numeric(decade))

# Sum climate-only columns across landscapes, 
# pick first simulation since there is no variation among simulations 
climate_df <- full_wildlife_df %>% 
  group_by(simulation, gcm, rcp, species, decade_num) %>% 
  summarise(across(starts_with('climate'), sum)) %>% 
  group_by(gcm, rcp, species, decade_num) %>% 
  summarise(across(starts_with('climate'), first)) 

# climate_df %>%  filter(species == 'tahu', decade == 'Hist.') %>% 
#   group_by(species) %>% 
#   summarise(across(where(is.numeric), mean)) %>% 
#   view()
# 
# climate_df %>%  filter(species == 'tahu') %>% 
#   group_by(gcm, rcp, species, decade) %>% 
#   summarise(across(where(is.numeric), mean)) %>% 
#   view()

# Habitat and combined 
hab_comb_df <- full_wildlife_df %>% 
  dplyr::select(-starts_with('climate'), -decade) %>% 
  # Sum values across all landscapes within each simulation
  group_by(simulation, gcm, rcp, species, decade_num) %>% 
  summarise(across(everything() & !landscape, sum)) %>% 
  # Caclulate mean across simulations
  group_by(gcm, rcp, species, decade_num) %>% 
  summarise(across(everything() & !simulation, 
                   list(med = median, 
                        low = ~ quantile(.x, 0.025), 
                        high = ~ quantile(.x, 0.975)))) 

# hab_comb_df %>%  filter(species == 'tahu', decade == 'Hist.') %>% 
#   group_by(species) %>% 
#   summarise(across(where(is.numeric), mean)) %>% 
#   view()
# 
# hab_comb_df %>%  filter(species == 'tahu') %>% 
#   group_by(gcm, rcp, species, decade) %>% 
#   summarise(across(where(is.numeric), mean)) %>% 
#   view()

# Plot colors
fig_cols <- c('#807dba','#3f007d','#41ab5d','#00441b','#f16913','#7f2704')
gcm_labs <- c("canesm2" = 'Wet', "hadgem2es" = 'Dry')
sppName_labs <- c('bbwo' = 'Black-backed Woodpecker', 'mart' = 'North American marten', 'tahu' = 'Red squirrel')
sppLat_labs <- c('bbwo' = 'P. arcticus', 'mart' = 'M. americana', 'tahu' = 'T. hudsonicus')
rcp_labs <- c('45' = 'Warm', '85' = 'Hot')
# names(gcm_labs) <- unique(habitat_df$gcm)
# lands_labs <- c('Grand Teton NP', 'Custer-Gallatin NF')
# names(lands_labs) <- unique(plot_df$landscape)

# Theme elements that will be common to all plots
theme_figs = theme_bw() %+replace%
  theme(
    axis.text = element_text(color="black",size=6),
    axis.title = element_text(size=8),
    axis.ticks = element_line(color="black",size=0.25),
    # facet titles
    strip.background = element_blank(),
    strip.text = element_text(size=7),
    # legend
    #legend.position="right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    # plot title
    plot.title=element_text(hjust=0.5,vjust=0.5,size=8,face="plain"),
    plot.background = element_rect(color = 'transparent'))

# New habitat vs climate plots - bar plots
plot_df <- full_join(climate_df, hab_comb_df) %>% 
  mutate(combined_med = combined_med_med,
         combined_low = combined_high_low,
         combined_high = combined_low_high) %>% 
  pivot_longer(cols = -c(gcm,rcp,species,decade_num), 
               names_to = c("set", ".value"),
               names_pattern = "(.+)_(.+)") %>% 
  filter(set %in% c('climate','habitat','combined')) %>%
  mutate(set = factor(set, levels = c('climate','habitat','combined')))


hab_area_plots <- map(species, function(spp){
  df <- filter(plot_df, species == !!spp) 
  
  if(spp == 'bbwo'){
    theme_spp <-
      theme(axis.title.x = element_blank())}
  if(spp == 'mart'){
    theme_spp <-
      theme(axis.title.x = element_blank(),
            strip.text.x = element_blank())}
  if(spp == 'tahu'){
    theme_spp <-
      theme(strip.text.x = element_blank(),
            legend.position = 'bottom')} 

  p <- 
    ggplot(df, aes(x = decade_num, 
                   y = med/279488, 
                   ymin = low/279488,
                   ymax = high/279488
                   )) +
    # geom_area(aes(fill = set), position = 'identity', color = 'black', alpha = 0.7) +
    # geom_line(aes(color = set)) +
    # geom_pointrange(aes(fill = set), shape = 21, size = 0.8, fatten = 2) +
    geom_bar(aes(fill = set),
             color = 'black', stat = 'identity', position = "dodge", width = 0.9,
             size = 0.3) +
    geom_errorbar(aes(group = set),
                  position = position_dodge(0.9), width = 0.5) +
    facet_grid(gcm~rcp, 
               labeller = labeller(rcp = rcp_labs, gcm = gcm_labs),
               scales = 'free_x', space = 'free_x') +
    scale_fill_manual('Model', values = c('#35978f','#bf812d','#762a83'),
                      labels = c('Climate','Vegetation','Joint habitat')) +
    scale_color_manual('Model', values = c('#35978f','#bf812d','#762a83'),
                      labels = c('Climate','Vegetation','Joint habitat')) +
    scale_y_continuous(breaks = c(0,0.5, 1), limits = c(0,1.1)) +
    # scale_y_continuous(breaks = seq(0,250000, 100000), labels = seq(0,250,100),
    #                    sec.axis = sec_axis(~ ./279488, 
    #                                        name = 'Proportion of study area',
    #                                        breaks = c(0,0.5))) +
    scale_x_continuous(breaks = c(1,2,3,4),
                       labels = c('Hist.','2030','2050','2070')) +
    labs(y = 'Proportion of focal landscapes', x = 'Time period') + #
    theme_figs +
    theme(legend.position = 'none',
          plot.margin = margin(1,5.5,1,5.5)) +
    theme_spp 
     
  return(p)
})

cowplot::plot_grid(hab_area_plots[[1]],hab_area_plots[[2]],hab_area_plots[[3]], 
                   nrow = 3, rel_heights = c(1.04,1,1.33), 
                   labels =  sppName_labs, 
                   label_x = c(-0.06,-0.04,0.032), 
                   label_y = c(0.96,1,1), 
                   label_size = 8, label_fontface = 'plain')
ggsave(filename = 'habitat_area_stack.png',
       height = 180, width = 100, units = 'mm',
       dpi = 600)


# Each niche dimension on a separate axis - vector plot
vector_df <- full_join(hab_comb_df, climate_df) %>% 
  filter(gcm == 'hadgem2es') %>% 
  select(gcm, rcp, species, decade_num, habitat_med, climate_med) %>% 
  group_by(gcm, rcp, species) %>% 
  mutate(start_hab = habitat_med[decade_num == 1],
         start_clim = climate_med[decade_num == 1]) %>% 
  # mutate(start_hab = habitat_med[decade == 'Hist.'],
  #        start_clim = climate_med[decade == 'Hist.']) %>% 
  rowwise() %>% 
  mutate(hab_diff = (habitat_med - start_hab)/start_hab,
         clim_diff = (climate_med - start_clim)/start_clim) 

vector_df_pts <- vector_df %>% 
  filter(decade_num !=4)

library(ggforce)

vector_plot <- 
  ggplot(vector_df) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_path(aes(x = clim_diff, y = hab_diff, group = interaction(rcp, species)),
            arrow = arrow(length = unit(0.15, "cm")), size = 0.6, color = 'grey20') +
  geom_path(aes(x = clim_diff, y = hab_diff, color = interaction(rcp, species)), 
            arrow = arrow(length = unit(0.15, "cm")), size = 0.5) +
  geom_point(data = vector_df_pts,
             aes(x = clim_diff, y = hab_diff, fill = interaction(rcp, species)),
             shape = 21, size = 0.7) +
  facet_zoom(xy = species == 'mart', xlim = c(-1,1), ylim = c(-1,1)) +
  scale_shape_manual(values = c(21, 24)) + 
  scale_fill_manual('Species; scenario', 
                    values = fig_cols, 
                    labels = c('BBWO; Warm','BBWO; Hot',
                               'MART; Warm','MART; Hot',
                               'TAHU; Warm','TAHU; Hot')) +
  scale_color_manual('Species; scenario',
                     values = fig_cols,
                     labels = c('BBWO; Warm','BBWO; Hot',
                                'MART; Warm','MART; Hot',
                                'TAHU; Warm','TAHU; Hot')) +
  # scale_x_continuous(expand = c(0,0), breaks = c(-1,-0.5,0,0.5,1)) +
  # scale_y_continuous(expand = c(0,0)) +
  theme_bw() %+replace%
  theme(
    axis.text = element_text(color="black",size=6),
    axis.title = element_text(size=8),
    axis.ticks = element_line(color="black",size=0.25),
    # legend
    legend.position ="bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.height = unit(0.25, 'cm'),
    legend.key.width = unit(0.25, 'cm')) +
  labs(x = 'Climate change', y = 'Vegetation change')
vector_plot

ggsave(filename = 'vector_plot.png',
       height = 75, width = 100, units = 'mm',
       dpi = 600)




### OLD -------------------------------------------------------------------------
# Using facets to create a broken axis... not quite working
clim_quad_plot <- function(spp){
  fake <- expand.grid('gcm' = unique(climate_df$gcm),
                      'rcp' = unique(climate_df$rcp),
                      'decade' = c(2009,2011),
                      'hist' = '1')
  
  c_df <- filter(climate_df, species == spp) %>% 
    mutate(hist = ifelse(decade == 2010, '1','2')) %>% 
    full_join(., fake)
  
  h_df <- filter(habitat_df, species == spp) %>% 
    mutate(hist = ifelse(decade == 2010, '1','2')) %>% 
    full_join(., fake)
  
  gmc_plot <- map(c('canesm2','hadgem2es'), function(gcm_i) {
    c_df_i <- filter(c_df, gcm == gcm_i)
    h_df_i <- filter(h_df, gcm == gcm_i)
    
    p <- 
      ggplot() +
      geom_line(data = c_df_i, aes(x = decade, y = climate_med, color = 'Climate-only')) +
      geom_line(data = h_df_i, aes(x = decade, y = habitat_med, color = 'Vegetation-only')) +
      geom_line(data = h_df_i, aes(x = decade, y = combined_med_med, color = 'Joint habitat')) +
      geom_pointrange(data = c_df_i, aes(x = decade, y = climate_med, ymin = climate_low, ymax = climate_high, color = 'Climate-only')) +
      geom_pointrange(data = h_df_i, aes(x = decade, y = habitat_med, ymin = habitat_low, ymax = habitat_high, color = 'Vegetation-only')) +
      geom_pointrange(data = h_df_i, aes(x = decade, y = combined_med_med, ymin = combined_med_low, ymax = combined_med_high, color = 'Joint habitat')) +
      facet_grid(rcp~hist, 
                 labeller = labeller(rcp = rcp_labs),
                 scales = 'free_x', space = 'free_x') +
      scale_color_manual(values = c('#35978f','#bf812d','#762a83')) +
      scale_y_continuous(breaks = seq(0,250000, 50000), labels = seq(0,250,50),
                         sec.axis = sec_axis(~ ./279488, name = 'Proportion of study area')) +
      scale_x_continuous(breaks = c(2010,2030,2050,2070), 
                         labels = c('Hist.','2030','2050','2070')) +
      labs(y = 'Habitat area (ha)', x = 'Year') +
      theme_bw(base_size = 14) +
      theme(strip.text.x = element_blank(),
            strip.background =  element_rect('transparent'),
            strip.text = element_text(face = 'bold'))
    return(p)
  })
  
  
}
combined_plot




# The full 'combined' niche model 
combined_plot <- 
  ggplot(habitat_df) +
  geom_line(aes(x = decade, 
                y = combined_med_med, 
                group = as.character(rcp),
                color = interaction(rcp, species))) +
  geom_pointrange(aes(x = decade, 
                      y = combined_med_med, 
                      ymin = combined_med_low, 
                      ymax = combined_med_high,
                      color = interaction(rcp, species))) +
  facet_grid(species~gcm, labeller = labeller(gcm = gcm_labs, species = spp_labs)) +
  scale_color_manual('RCP; Species', values = fig_cols) +
  #scale_x_continuous(breaks = seq(2010,2100,10)) +
  scale_y_continuous(breaks = seq(0,125000, 25000), labels = seq(0,125,25),
                     sec.axis = sec_axis(~ ./279488, name = 'Proportion of study area')) +
  labs(y = 'Habitat area (ha)', x = 'Year') +
  theme_bw(base_size = 14) +
  theme(strip.background =  element_rect('transparent'),
        strip.text = element_text(face = 'bold', size = 14))

ggsave(combined_plot, filename = 'land5_sims20_combined.png', width = 9, height = 7, dpi = 300)


# Plotting niche dimensions separately 
model_cols <- c('#35978f','#762a83','#bf812d')

# Just plot RCP 85
habitat_85 <- filter(habitat_df, rcp == '85')
climate_85 <- filter(climate_df, rcp == '85')

habitat_v_climate <- 
  ggplot() +
  geom_ribbon(data = climate_85, aes(x = decade, ymin = climate_low, ymax = climate_high), fill = '#35978f', alpha = 0.4) +
  geom_line(data = climate_85, aes(x = decade, y = climate_med), color = '#35978f', size = 1.2) +
  geom_ribbon(data = habitat_85, aes(x = decade, ymin = habitat_low, ymax = habitat_high), fill = '#bf812d', alpha = 0.4) +
  geom_line(data = habitat_85, aes(x = decade, y = habitat_med), color = '#bf812d', size = 1.2) +
  geom_ribbon(data = habitat_85, aes(x = decade, ymin = combined_low_med, ymax = combined_high_med), fill = '#762a83', alpha = 0.4) +
  geom_line(data = habitat_85, aes(x = decade, y = combined_med_med), color = '#762a83', size = 1.2) +
  facet_grid(species~gcm, labeller = labeller(gcm = gcm_labs, species = spp_labs)) +
  #scale_x_continuous(breaks = seq(2020,2100,10)) +
  scale_y_continuous(breaks = seq(0,250000, 50000), labels = seq(0,250,50),
                     sec.axis = sec_axis(~ ./279488, name = 'Proportion of study area')) +
  labs(y = 'Habitat area (1000s ha)', x = 'Year') +
  theme_bw(base_size = 14) +
  theme(strip.background =  element_rect('transparent'),
        strip.text = element_text(face = 'bold', size = 14)) 

ggsave(habitat_v_climate, filename = 'habit_v_climate.png', width = 8, height = 8, dpi = 300)



ggsave(vector_plot, filename = 'vector_plot.png', width = 8, height = 6, dpi = 300)



# Version with overlapping RCPs
combined_plot <- 
  ggplot() +
  geom_line(data = climate_df,
            aes(x = decade, 
                y = climate_med, 
                color = 'Climate-only',
                group = rcp), 
            position = position_dodge(width = 0.3), size = 0.75) +
  geom_line(data = habitat_df,
            aes(x = decade, 
                y = habitat_med, 
                color = 'Vegetation-only',
                group = rcp), 
            position = position_dodge(width = 0.3), size = 0.75) +
  geom_line(data = habitat_df,
            aes(x = decade, 
                y = combined_med_med, 
                color = 'Joint habitat',      
                group = rcp), 
            position = position_dodge(width = 0.3), size = 0.75) +
  geom_pointrange(data = climate_df,
                  aes(x = decade, 
                      y = climate_med, 
                      ymin = climate_low,
                      ymax = climate_high,
                      color = 'Climate-only',
                      fill = 'Climate-only',
                      shape = rcp), 
                  position = position_dodge(width = 0.3), size = 0.75) +
  geom_pointrange(data = habitat_df,
                  aes(x = decade, 
                      y = habitat_med, 
                      ymin = habitat_low, 
                      ymax = habitat_high,
                      color = 'Vegetation-only',
                      fill = 'Vegetation-only',
                      shape = rcp), 
                  position = position_dodge(width = 0.3), size = 0.75) +
  geom_pointrange(data = habitat_df,
                  aes(x = decade, 
                      y = combined_med_med, 
                      ymin = combined_med_low, 
                      ymax = combined_med_high,
                      color = 'Joint habitat',      
                      fill = 'Joint habitat',
                      shape = rcp), 
                  position = position_dodge(width = 0.3), size = 0.75) +
  facet_grid(species~gcm, labeller = labeller(gcm = gcm_labs, species = spp_labs)) +
  scale_color_manual(values = c('#35978f','#bf812d','#762a83')) +
  scale_fill_manual(values = c('white','white','white'), guide = F) +
  scale_shape_manual(values = c(21,19)) +
  scale_y_continuous(breaks = seq(0,250000, 50000), labels = seq(0,250,50),
                     sec.axis = sec_axis(~ ./279488, name = 'Proportion of study area')) +
  scale_x_discrete(breaks = c('2010','2030','2050','2070','2080'), 
                   labels = c('Hist.','2030','2050','2070','2080')) +
  labs(y = 'Habitat area (ha)', x = 'Year') +
  theme_bw(base_size = 14) +
  theme(strip.background =  element_rect('transparent'),
        strip.text = element_text(face = 'bold', size = 14))
combined_plot


# OLD VECTOR PLOT

combined_sum <- full_join(hab_comb_df, climate_df) %>% 
  filter(gcm == 'hadgem2es') %>% 
  mutate(period = ifelse(decade < 2050, 'early', 'late')) 
combined_mean <- combined_sum %>% 
  group_by(gcm, rcp, species, period) %>% 
  summarise(across(everything(), ~ max(.x) - mean(.x)))

vector_plot <- 
  ggplot(combined_sum) +
  geom_vline(xintercept = 150000) +
  geom_hline(yintercept = 75000) +
  geom_point(aes(x = climate_med, y = habitat_med, fill = interaction(rcp, species), shape = period), size = 3, alpha = 0.9) +
  # geom_point(data = combined_mean, aes(x = climate_med, y = habitat_med, fill = interaction(rcp, species), shape = period), 
  #            size = 6, alpha = 0.6) +
  # geom_point(data = combined_mean, aes(x = climate_med, y = habitat_med, color = interaction(rcp, species), shape = period), 
  #            size = 6, stroke = 1) +
  geom_path(data = combined_mean, 
            aes(x = climate_med, y = habitat_med, group = interaction(rcp, species)), 
            arrow = arrow(length = unit(0.5, "cm")), size = 1.5, color = 'grey20') +
  geom_path(data = combined_mean, 
            aes(x = climate_med, y = habitat_med, color = interaction(rcp, species)), 
            arrow = arrow(length = unit(0.5, "cm")), size = 1) +
  #geom_abline(intercept = 0, slope = 0.5) +
  
  #geom_step(aes(x = climate_med, y = habitat_med, color = interaction(rcp, species), shape = gcm), size = 1) +
  scale_shape_manual(values = c(21, 24)) + 
  scale_y_continuous(breaks = seq(0,250000, 50000), labels = seq(0,250,50)) +
  scale_x_continuous(breaks = seq(0,250000, 50000), labels = seq(0,250,50)) +
  scale_fill_manual('RCP; Species', values = fig_cols) +
  scale_color_manual('RCP; Species', values = fig_cols) +
  theme_bw(base_size = 14) +
  labs(x = 'Climate niche (ha)', y = 'Habitat niche (ha)')

























