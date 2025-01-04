library(ggplot2)
library(dplyr)
library(plotly)
library(raster)
library(rnaturalearth)
library(sf)
library(scales)
library(nmfspalette)
library(terra)
library(here)

#function to calculate anomaly, prep map data for dashboard to display
prep_map_data <- function(indicator, raster_ann, raster_clim){ #, min_x, max_x, min_y, max_y){
  
  #indicator <- "SST"
  target <- raster_ann
  names(target) <- c("")
  raster_df <- target%>%
    rasterToPoints(.) %>% #convert to points   
    as.data.frame() %>% #convert to dataframe
    mutate(x_disp = ifelse(x >180, -x + 360, x)) %>%
    mutate(labels = paste0(signif(x_disp, 3), ", ", signif(y, 3),
                           "<br>", indicator,
                           ": ", signif(layer, 3))) %>%
    mutate(ID = indicator)
  
  anom <- raster_ann-raster_clim
  raster_anom_df <- anom %>%
    rasterToPoints(.) %>% #convert to points   
    as.data.frame() %>% #convert to dataframe
    mutate(x_disp = ifelse(x >180, -x + 360, x)) %>% #recenter at 180
    mutate(layer_disp = ifelse(layer > 0, paste0("+ ", layer), layer)) %>% #format hovertext
    mutate(layer_disp = ifelse(layer == 0, 0, layer)) %>%
    mutate(labels = paste0(signif(x_disp, 3), ", ", signif(y, 3),
                           "<br>", indicator,
                           " Anomaly: ", signif(layer, 3))) %>%
    mutate(ID = paste0(indicator, "_anom"))
  
  raster_df['x'] <-  - raster_df['x_disp']  #raster_df['x_disp'] -180  
  raster_anom_df['x'] <- - raster_df['x_disp'] #raster_anom_df['x_disp'] -180  
  #combine datasets
  raster_df <- bind_rows(raster_df, raster_anom_df)

  #return data
  return(list(target, anom, raster_df))
  
}

#read in raster data
sst_2023 <- raster(here("Indicator_Dashboard/PreProcessing", "sst-MH-2023-mean.nc"))
sst_clim <- raster(here("Indicator_Dashboard/PreProcessing", "sst-MH-1985-2009-mean.nc"))
chl_2023 <- raster(here("Indicator_Dashboard/PreProcessing", "chl-MH-2023-mean.nc"))
chl_clim <- raster(here("Indicator_Dashboard/PreProcessing", "chl-MH-1998-2009-mean.nc"))

# get lat coordinates from raw data
# Different data sources have different bounds
satellite_bbox <- bbox(sst_2023)
s_min_x <- satellite_bbox[1,1] 
s_max_x <- satellite_bbox[1,2]
s_min_y <- satellite_bbox[2,1]
s_max_y <- satellite_bbox[2,2]

#set target resolution in degrees
target_res <- 0.05

# Workaround for missing information
srs <- sst_2023@srs
chl_2023@srs <- srs

#get prepped data 
sst_out <- prep_map_data("SST", sst_2023, sst_clim)
chl_out <- prep_map_data("Chl", chl_2023, chl_clim)

#quick tests - annual maps
#plot(tatd_out[[1]])
plot(sst_out[[1]])
plot(chl_out[[1]])
plot(log(chl_out[[1]]))


#quick tests - anomaly maps
plot(sst_out[[2]])
plot(chl_out[[2]])
plot(log(chl_out[[2]]))


#combine data and write to file
raster_df_all <- bind_rows(sst_out[[3]], chl_out[[3]])
write.csv(raster_df_all, here("Indicator_Dashboard", "Data", "Dashboard_Map_Data_2023.csv"))
#write.csv(sst_out[[3]], here("Sea_Surface_Temperature", "SST_map_data_2023.csv"))
#write.csv(chl_out[[3]], here("Ocean_Color", "Chl_map_data_2023.csv"))

#get lon coordinates from cropped data
cropped_bbox <- bbox(sst_out[[1]])
#min_x <- cropped_bbox[1,1] 
#max_x <- cropped_bbox[1,2]# 180-cropped_bbox[1,2]
min_x <- cropped_bbox[1,1]-360
max_x <- cropped_bbox[1,2]-360
min_y <- cropped_bbox[2,1]
max_y <- cropped_bbox[2,2]

#read in coast data
world_coasts <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf") %>% st_make_valid()
world_countries <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") %>% st_make_valid()

#CRS to recenter for plotting - warnings here ok
#plot_crs <- "+proj=longlat +x_0=0 +y_0=0 +lat_0=0 +lon_0=180 +datum=WGS84 +no_defs"
plot_crs <- srs
#base maps
land_proj <- world_countries %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = min_y, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTIPOLYGON") #recast for plotting

coast_proj <- world_coasts %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = min_y, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTILINESTRING") #recast for plotting

#quick test
ggplot() + geom_sf(data = land_proj, color = "grey", fill = "grey") +
  geom_sf(data = coast_proj, color = "black")


#save data
saveRDS(land_proj, "Indicator_Dashboard/Data/rnatearth_MH_land.RData")
saveRDS(coast_proj, "Indicator_Dashboard/Data/rnatearth_MH_coast.RData")

##############################not run###################################
#base maps for ENSO
land_proj_enso <- world_countries %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = -5, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTIPOLYGON") #recast for plotting

coast_proj_enso <- world_coasts %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = min_x, ymin = -5, 
                       xmax = max_x, ymax = max_y)) %>% #crop to area
  st_cast(., "MULTILINESTRING") #recast for plotting

#base maps for Tropical Cyclones
land_proj_tcs <- world_countries %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = -90, ymin = -50, 
                       xmax = 100, ymax = 50)) %>% #crop to area
  st_cast(., "MULTIPOLYGON") #recast for plotting

coast_proj_tcs <- world_coasts %>% 
  st_break_antimeridian(lon_0 = 180) %>% #trim out dateline to avoid artifacts
  st_transform(crs = plot_crs) %>% #adjust CRS for Pacific
  st_crop(x = ., y = c(xmin = -90, ymin = -50, 
                       xmax = 100, ymax = 50)) %>% #crop to area
  st_cast(., "MULTILINESTRING") #recast for plotting

#quick test
ggplot() + geom_sf(data = land_proj, color = "grey", fill = "grey") +
  geom_sf(data = coast_proj, color = "black")

ggplot() + geom_sf(data = land_proj_enso, color = "grey", fill = "grey") +
  geom_sf(data = coast_proj_enso, color = "black")

ggplot() + geom_sf(data = land_proj_tcs, color = "grey", fill = "grey") +
  geom_sf(data = coast_proj_tcs, color = "black")

#save data
saveRDS(land_proj, "Indicator_Dashboard/Data/rnatearth_land.RData")
saveRDS(coast_proj, "Indicator_Dashboard/Data/rnatearth_coast.RData")
 
saveRDS(land_proj_enso, "Indicator_Dashboard/Data/rnatearth_enso_land.RData")
saveRDS(coast_proj_enso, "Indicator_Dashboard/Data/rnatearth_enso_coast.RData")

saveRDS(land_proj_tcs, "Indicator_Dashboard/Data/rnatearth_tcs_land.RData")
saveRDS(coast_proj_tcs, "Indicator_Dashboard/Data/rnatearth_tcs_coast.RData")
