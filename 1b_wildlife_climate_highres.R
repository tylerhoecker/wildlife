library(tidyverse)
library(sf)
library(raster)
library(biomod2)
library(tmap)

setwd('/Volumes/thexternal/iLand_gye/wildlife_pc/biomod2ing/')

# Load occurrence data ----------------------------------------------------------
occ_records <- readRDS('../vertnet_data/vertnet_occ_data_thin_1km.rds')
# Save species names as a shortcut for later
species <- names(occ_records)

# Load historical climate data -------------------------------------------------
# Data downloaded from WordClim project: https://worldclim.org/data/index.html 

# Use these bioclim predictors
#bioclim_vars <- c(1,5,6,12,13,14)
bioclim_vars <- seq(1,19,1)

# These variables need to be dividied

# Load data that has already been cropped to the bounding box of the occurrence points. Re-do if points get updated.
current_stack <- stack('bioclim_data/wc21_cropped/bioclim_NA_crop.tif') 
names(current_stack) <- paste0('bio_',bioclim_vars)

# current_coarse <- aggregate(current_stack, fact = 10, fun = mean, na.rm = T)
# writeRaster(current_coarse, 'bioclim_data/wc21_cropped/bioclim_NA_coarse.tif', overwrite = T)
current_coarse <- stack('bioclim_data/wc21_cropped/bioclim_NA_coarse.tif') 
names(current_coarse) <- paste0('bio_',bioclim_vars)

# To run again:

# current_stack <- paste0('bioclim_data/wc2.1_30s_bio/wc2.1_30s_bio_', bioclim_vars,'.tif') %>%
#   stack()
# names(current_stack) <- paste0('bio_',bioclim_vars)

# Crop this global product to the extent of occurrence data for these species
# occ_extents <- map(unname(occ_records), extent)
# max_extents <- do.call(raster::merge, occ_extents)
# current_stack <- crop(current_stack, max_extents)
# writeRaster(current_stack, 'bioclim_data/wc21_cropped/bioclim_NA_crop.tif', overwrite = T)

# For quicker testing - predict only on the study area!
# landscapes <- c('grte','synp','nynp','wynp','sgye')
# ls_outlines <- map(landscapes, function(x){
#   read_sf(paste0('/Volumes/thexternal/iLand_gye/stand_grids/',x,'_polygon/',x,'_polygon.shp')) %>%
#     st_transform(., crs = st_crs(current_stack)) %>%
#     nngeo::st_remove_holes()}) %>%
#   set_names(landscapes) %>%
#   bind_rows(.id = 'landscape')
# 
# ls_bb <- st_bbox(st_buffer(ls_outlines, 0.1)) %>% st_as_sfc()
# current_gye <- crop(current_stack, st_buffer(ls_outlines, 0.1))
# writeRaster(current_gye, 'bioclim_data/wc21_cropped/bioclim_gye.tif', overwrite = T)
current_gye <- stack('bioclim_data/wc21_cropped/bioclim_gye.tif') 
names(current_gye) <- paste0('bio_',bioclim_vars)

# Set up SDM modeling data
# .x = occ_records[[3]]
# .y = names(occ_records)[[3]]

# Create the models -----------------------------------------------------------
# This is slow, so load from RData
# biomod_formatted_dfs <-
#   map2(occ_records, names(occ_records), function(.x, .y){
#     # Format the data per BIOMOD requirements
#     BIOMOD_FormatingData(resp.var = rep(1, dim(.x)[1]),
#                          resp.xy = st_coordinates(.x),
#                          expl.var = current_stack,
#                          resp.name = .y,
#                          PA.nb.rep = 1,
#                          PA.nb.absences = 5000,
#                          #PA.dist.min = 0.1,
#                          #PA.dist.max = 1,
#                          PA.strategy = 'random')})
# saveRDS(biomod_formatted_dfs, 'biomod_formatted_dfs_1km.R')
biomod_formatted_dfs <- readRDS('biomod_formatted_dfs_1km.R')

n_runs <- 10
# Already ran 10 models -----------
biomod_models <-
  map(biomod_formatted_dfs, function(biomod_df){

    # Explicity define Maxent tuning parameters
    # Commented out variables note unchanged defaults
    biomod_opts <-
      BIOMOD_ModelingOptions(
        MAXENT.Phillips = list(
          memory_allocated = 2000,
          # background_data_dir = 'default',
          # maximumbackground = 'default',
          # maximumiterations = 200,
          # visible = FALSE,
          # linear = TRUE,
          # quadratic = TRUE,
          product = FALSE,
          # threshold = TRUE,
          # hinge = TRUE,
          # lq2lqptthreshold = 80,
          # l2lqthreshold = 10,
          # hingethreshold = 15,
          # beta_threshold = -1,
          # beta_categorical = -1,
          # beta_lqp = -1,
          # beta_hinge = -1,
          betamultiplier = 2
          #defaultprevalence = 0.5
          )
    )

    # Create models
    biomod_mods <-
      BIOMOD_Modeling(
        biomod_df,
        models = c("MAXENT.Phillips"), #,"GLM","RF","GBM"
        models.options = biomod_opts,
        NbRunEval = n_runs,
        DataSplit = 70,
        VarImport = 10,
        models.eval.meth = c("TSS","ROC"),
        SaveObj = TRUE,
        rescal.all.models = FALSE,
        do.full.models = FALSE,
        modeling.id = 'maxentOnly')})
# ------------

biomod_models <- map(species, function(sp){
  grep('.*maxentOnly.models.out', 
       list.files(paste0('/Volumes/thexternal/iLand_gye/wildlife_pc/biomod2ing/',sp,'/')), 
       value = T)
})

# Load the models into the global environment
# This does not work in a function because they need to be in the global... apparently
for(i in 1:length(species)){
  load(paste0(species[i],'/',biomod_models[i]))
  BIOMOD_LoadModels(get(biomod_models[[i]]))
}

# Don't run --------------------------------------------------------------------
map2_df(biomod_models, species, function(mod_obj, sp_name){
  the_model <- get(mod_obj)
  get_evaluations(the_model, as.data.frame = T) %>%
    # If using multiple model types
    # separate(col = Model.name, into = c('Model_type','Run'), sep = '_') %>%
    # dplyr::select(-Evaluating.data, -Run) %>%
    # group_by(Model_type, Eval.metric) %>%
    group_by(Eval.metric) %>% 
    summarise(across(where(is.numeric),
                     list(mean = ~mean(.x, na.rm = T), se = ~(sd(.x, na.rm = T)/sqrt(n_runs))))) %>% 
    mutate(species = sp_name) 
  }) %>% 
  write_csv(paste0('eval_stats_allSpecies.csv'))


# Variable importance and response curves... models loaded from disk
walk2(biomod_models, species, function(mod_obj, sp_name){

  var_import <- get_variables_importance(the_model, as.data.frame = T) %>%
    as.data.frame.table() %>%
    #mutate(import_z = scale(Freq)) %>%
    group_by(Var1) %>%
    summarise(mean_import = mean(Freq, na.rm = T),
              sd = sd(Freq, na.rm = T)) %>%
    arrange(desc(mean_import))

  write_csv(var_import, paste0('var_import_',sp_name,'.csv'))

  # Produce response curve data and do some adjustments for plotting purposes
  response_df <-
    response.plot2(
      models = get_built_models(the_model),
      Data = get_formal_data(the_model, 'expl.var'),
      show.variables = get_formal_data(the_model,'expl.var.names'),
      do.bivariate = FALSE,
      fixed.var.metric = 'mean',
      data_species = get_formal_data(the_model, 'resp.var'), # Calculate median only over presences... not sure...
      plot = F) %>%
    as_tibble() %>%
    full_join(var_import, by = c('expl.name' = 'Var1')) %>%
    mutate(expl.name = factor(expl.name)) %>%
    group_by(pred.name, expl.name) %>%
    mutate(sd_pred = sd(pred.val)) %>%
    ungroup() %>%
    #mutate(pred.val = ifelse(sd_pred <0.0001, NA, pred.val)) %>%
    separate(pred.name, into = c('trash1','trash2','run','model_type'))

  # Add column of labels for the bioclim variables that make sense
  bio_labs <- c('Ann.T', 'T Range','Isothermality','Seasonality','Max T Warmest M','Min T Coldest M','T Ann. Range','Mean T Wettest Q',
                'Mean T Driest Q', 'Mean T Warmest Q', 'Mean T Coldest Q', 'Ann. Precip', 'Precip Wettest M', 'Precip Driest M',
                'Precip Season', 'Precip Wettest Q', 'Precip Driest Q', 'Precip Warmest Q','Precip Coldest Q')
  names(bio_labs) <- as.character(unique(response_df$expl.name))

  # Rename biomod variables and order the dataframe by variable importance
  response_df <- response_df %>%
    ungroup() %>%
    mutate(bio_vars = recode(expl.name, !!!bio_labs),
           bio_vars = fct_reorder(bio_vars, mean_import, .desc = T))

  ggplot(response_df) +
    geom_line(aes(x = expl.val, y = pred.val, group = run), alpha = 0.7) +
    facet_wrap(~bio_vars, scales = 'free_x') +
    labs(x = 'Predictor variable value (units vary)', y = 'Estimated probability') +
    ggtitle(sp_name) +
    theme_bw(base_size = 11)
  ggsave(filename = paste0('response_curves_',sp_name,'.png'), width = 10, height = 8)
})

# Make variable importance and evaluation stats figures
vimp_csvs <- paste0('var_import_',species,'.csv') %>% 
  set_names(species)
vimp_df <- map_df(vimp_csvs, read_csv, .id = 'species') %>% 
  mutate(se2 = 1.96*sd/sqrt(3))

# Add column of labels for the bioclim variables that make sense
bio_labs <- c('Ann T', 'T Range','Isothermality','Seasonality','Max T Warmest M','Min T Coldest M','T Ann Range','Mean T Wettest Q',
              'Mean T Driest Q', 'Mean T Warmest Q', 'Mean T Coldest Q', 'Ann. Precip', 'Precip Wettest M', 'Precip Driest M',
              'Precip Season', 'Precip Wettest Q', 'Precip Driest Q', 'Precip Warmest Q','Precip Coldest Q')
names(bio_labs) <- as.character(unique(vimp_df$Var1))

vimp_high <- vimp_df %>%  
  group_by(Var1) %>% 
  mutate(avg_iv = mean(mean_import)) %>% 
  ungroup() %>% 
  # group_by(species) %>%
  # filter(mean_import > median(mean_import)) 
  mutate(Var1 = recode(Var1, !!!bio_labs),
         Var1 = fct_reorder(Var1, avg_iv, .desc = F))

dodge <- position_dodge(width=0.9)
spp_cols <- c('#807dba','#41ab5d','#f16913') #c('#f16913','#7f2704','#41ab5d','#00441b','#807dba','#3f007d')
spp_labs <- c('bbwo' = 'P. arcticus', 'mart' = 'M. americana', 'tahu' = 'T. hudsonicus')

ggplot(vimp_high) +
  geom_col(aes(x = Var1, y = mean_import, fill = species), color = 'black', size = 0.35,
           position = dodge) +
  # geom_errorbar(aes(x = Var1, 
  #                   ymin = mean_import - se2, 
  #                   ymax = mean_import + se2,
  #                   color = species),
  #               position = dodge, width = 0.25) +
  scale_color_manual('Species', values = spp_cols, labels = spp_labs) +
  scale_fill_manual('Species', values = spp_cols, labels = spp_labs) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(y = 'Mean variable importance') +
  #facet_wrap(~species, ncol = 1) +
  coord_flip() +
  theme_bw() %+replace%
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black",size=6),
    axis.title = element_text(size=8),
    axis.ticks = element_line(color="black",size=0.25),
    axis.title.y = element_blank(),
    # legend
    legend.position="bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.25, 'cm'))

ggsave(filename = 'variable_importance.png',
       height = 120, width = 100, units = 'mm',
       dpi = 600)

# Create ensemble models and project them on contemporary/historical data
walk2(biomod_models, species, function(mod_obj, sp_name){

  the_model <- get(mod_obj)
  model_names <- BIOMOD_LoadModels(the_model)

  biomod_em <-
    BIOMOD_EnsembleModeling(
      modeling.output = the_model,
      chosen.models = model_names,
      em.by = 'all',
      eval.metric = c('ROC'),
      eval.metric.quality.threshold = 0.8,
      prob.mean = T,
      prob.ci = T,
      models.eval.meth = c('TSS','ROC'))

  em_cutoff <- get_evaluations(biomod_em, as.data.frame = T)

  biomod_proj_na <-
    BIOMOD_Projection(
      modeling.output = the_model,
      selected.models = 'all',
      new.env = stack(current_coarse),
      proj.name = 'hist_na',
      build.clamping.mask = FALSE)

  biomod_proj_gye <-
    BIOMOD_Projection(
      modeling.output = the_model,
      selected.models = 'all',
      new.env = stack(current_gye),
      proj.name = 'hist_gye',
      build.clamping.mask = FALSE)

  biomod_em_proj_na <-
    BIOMOD_EnsembleForecasting(
      EM.output = biomod_em,
      projection.output = biomod_proj_na)

  biomod_em_proj_gye <-
    BIOMOD_EnsembleForecasting(
      EM.output = biomod_em,
      projection.output = biomod_proj_gye)

  })

# Map the historical model projections
map_colors <- list('Purples','Greens','Oranges') %>%
  set_names(species)
map_labs <- list('Picoides arcticus','Martes americana','Tamiasciurus hudsonicus') %>%
  set_names(species)
walk(species, function(sp){

  na_dist <- raster(paste0(sp, '/proj_hist_na/proj_hist_na_',sp,'_ensemble.grd'))
  gye_dist <- raster(paste0(sp, '/proj_hist_gye/proj_hist_gye_',sp,'_ensemble.grd'))

  north_america_map <-
    tm_shape(na_dist/10) +
    tm_raster(title = 'ROR',
              palette = map_colors[[sp]], breaks = seq(0,100,10)) +
    tm_shape(occ_records[[sp]]) +
    tm_bubbles(size = 0.03, col = '#d9d9d9', border.col = 'black') +
    tm_shape(ls_bb) +
    tm_polygons(col = '#feb24c', border.col = 'black', alpha = 0.6) +
    tm_scale_bar(breaks = c(0,1000,2000), text.size = 0.8) +
    tm_layout(title = paste0(map_labs[[sp]], ' (n = ', dim(occ_records[[sp]])[1], ')'),
              title.position = c('left', 'top'),
              title.size = 1,
              legend.position = c('right', 'top'),
              legend.title.size = 0.9,
              legend.text.size = 0.8,
              bg.color = '#6baed6')

  gye_map <-
    tm_shape(gye_dist/10) +
    tm_raster(palette = map_colors[[sp]], breaks = seq(0,100,10)) +
    tm_shape(occ_records[[sp]]) +
    tm_bubbles(size = 0.03, col = '#d9d9d9', border.col = 'black') +
    tm_shape(ls_outlines) +
    tm_borders(col = 'black') +
    tm_scale_bar(breaks = c(0,10,20), text.size = 0.8) +
    tm_layout(legend.show = FALSE)

  tmap_save(north_america_map,
            insets_tm = gye_map,
            insets_vp = grid::viewport(0.102, 0.48, width = 1, height = 0.89),
            filename = paste0('proj_hist_map_',sp,'.png'), dpi=600)
})

# ------------------------------------------------------------------------------

# Access previously created ensemble models
biomod_ems <- map(species, function(sp){
  grep('.*maxentOnlyensemble.models.out', 
       list.files(paste0('/Volumes/thexternal/iLand_gye/wildlife_pc/biomod2ing/',sp,'/')), 
       value = T)})

# Now load the ensemble models (again, loop instead of functio)
for(i in 1:length(species)){
  load(paste0(species[i],'/',biomod_ems[i]))
  BIOMOD_LoadModels(get(biomod_ems[[i]]))
}


# Get the optimal cutoff (threshold) values for the ensemble models, save as a small species-length object
em_mean_cutoffs <- map(biomod_ems, function(em){
  get_evaluations(get(em), as.data.frame = T) %>% 
    filter(Eval.metric == 'ROC') %>%
    dplyr::select(Model.name, Cutoff) %>% 
    pivot_wider(names_from = 'Model.name', values_from = 'Cutoff') %>% 
    rename(mean = EMmeanByROC_mergedAlgo_mergedRun_mergedData, 
          lower = EMciInfByROC_mergedAlgo_mergedRun_mergedData, 
          upper =  EMciSupByROC_mergedAlgo_mergedRun_mergedData)}) %>% 
  set_names(species)

saveRDS(em_cutoffs, 'cutoff_values.Rds')

# Predict future climate models ------------------------------------------------
# TESTING
# landscape = c('synp')
# gcm= c('canesm2')
# rcp = c('45')
# period = c('2030')

future_biomod_fn <- function(landscape, gcm, rcp, period){
  
  # Avoid re-running things that are done --------------------------------------
  if(
    file.exists(paste0('future_rasters/',landscape,'_',gcm,'_',rcp,'_',period,'.tif'))
    ){
    print(paste0(landscape,'_',gcm,'_',rcp,'_',period, ' already complete'))
    return(NULL) 
  }
  # ----------------------------------------------------------------------------
  
  # Load future climate data  --------------------------------------------------
  # Same bio variables, but future downscaled data from: http://ccafs-climate.org/data_spatial_downscaling/  
  future_stack <- paste0('bioclim_data/ccafs_iland/',
                         landscape,'_',gcm,'_',rcp,'_',
                         period,'_bio_',bioclim_vars,'.tif') %>% 
    stack() 
  
  names(future_stack) <- names(current_stack)
  
  # The CCAFS data provide temperature in tenths of degrees (200 = 20.0 C)
  # CHECK 3, 4 
  temp_vars <- paste0('bio_',c(1,2,4,5,6,7,8,9,10,11))
  
  for(var in temp_vars){
    future_stack[[var]] <- future_stack[[var]]/10
  }
    
  # mod_objs = biomod_models[[2]]
  # sp_name = species[[2]]
  
  biomod_futures <- 
    pmap(list(biomod_models, biomod_ems, species), function(mod_obj, em_obj, sp_name){
      
      the_model <- get(mod_obj)
      the_em <- get(em_obj)
      #model_names <- BIOMOD_LoadModels(the_model)
      proj_name <- paste0(landscape,'_',gcm,'_',rcp,'_',period,'_future')
 
      biomod_proj_fut <- 
        BIOMOD_Projection(
          modeling.output = the_model,
          selected.models = 'all',
          new.env = future_stack, 
          proj.name = proj_name,
          build.clamping.mask = FALSE)
      
      biomod_em_proj_na <- 
        BIOMOD_EnsembleForecasting(
          EM.output = the_em,
          projection.output = biomod_proj_fut)
      
      fut_proj_rast <- 
        raster(paste0(sp_name,'/proj_',proj_name,'/proj_',proj_name,'_',sp_name,'_ensemble.grd'))
  
     # fut_map <- 
     #    tm_shape(fut_proj_rast/10) +
     #    tm_raster(title = 'ROR',
     #              palette = map_colors[[sp_name]], breaks = seq(0,100,10))
     # 
     # tmap_save(fut_map, filename = paste0(sp_name,'/',proj_name,'map_',sp_name,'.png'), dpi=600)
      
      return(fut_proj_rast)
    }) %>% 
    set_names(species)
  
  
  env_stack <- stack(biomod_futures) 
  
  writeRaster(env_stack, paste0('future_rasters/',landscape,'_',gcm,'_',rcp,'_',period,'.tif'), 
              #bylayer = T, 
              suffix = names(occ_records),
              overwrite= T)
 }

scenarios <- expand.grid('landscape' = c('grte','synp','nynp','wynp','sgye'), # ,'synp','nynp','wynp','sgye'
                         'gcm'= c('canesm2','hadgem2es'), # 'canesm2','hadgem2cc','hadgem2es'
                         'rcp' = c('45','85'),
                         'period' = c(2030,2050,2070)) 

pwalk(scenarios, future_biomod_fn)







