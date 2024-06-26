
###                                                                    ###
# Go to section 'Calculating reference levels' for the final good method #
###                                                                    ###


if(Sys.info()[4] == "D01RI1700308") {
  wd <- ""
}else if(Sys.info()[4] == "S-JRCIPRAP320P") {
  wd <- ""
}else if(Sys.info()[4] %in% c("jeodpp-terminal-jd001-03", "jeodpp-terminal-03", "jeodpp-terminal-dev-12", 
                              "jeodpp-terminal-jd002-03", "jeodpp-terminal-jd004-03.cidsn.jrc.it",
                              "jeodpp-terminal-dev-jd002-12.cidsn.jrc.it")) {
  if(!dir.exists("/eos/jeodpp/home/users/rotllxa/Birds_Map_Indicators/")) 
    dir.create("/eos/jeodpp/home/users/rotllxa/Birds_Map_Indicators/")
  wd <- "/eos/jeodpp/home/users/rotllxa/Birds_Map_Indicators/"
  gbif_creds <- "/home/rotllxa/Documents/"
  WhereAmI <- "bdap"
}else if(Sys.info()[4] %in% c("MacBook-MacBook-Pro-de-Xavier.local",
                              "MacBookdeXavier.station")){
  if(!dir.exists("/Users/xavi_rp/Documents/D3_NRL_farmlandBirdsStudy/")) 
    dir.create("/Users/xavi_rp/Documents/D3_NRL_farmlandBirdsStudy/")
  wd <- "/Users/xavi_rp/Documents/D3_NRL_farmlandBirdsStudy/"
  WhereAmI <- "mac"
}else{
  wd <- ""
  gbif_creds <- "C:/Users/rotllxa/Documents/"
}

setwd(wd)

library(sf)
library(tidyverse)
library(data.table)
library(terra)


## Data ####
### Farmland birds indicator ####

if(WhereAmI == "mac"){
  list.files("./Study1_Richness/Map_Farmland_Bird_Indicators")
}else{
  list.files("./Map_Farmland_Bird_Indicators_Study1/")
}


if(WhereAmI == "mac"){
  fbi_eu_grid <- read_sf(dsn = "./Study1_Richness/Map_Farmland_Bird_Indicators/", layer = "fbi_eu_grid")
}else{
  fbi_eu_grid <- read_sf(dsn = "./Map_Farmland_Bird_Indicators_Study1/", layer = "fbi_eu_grid")
  
}
fbi_eu_grid

kk <- sample(1:nrow(fbi_eu_grid), 100, replace = TRUE)
fbi_eu_grid_kk <- fbi_eu_grid[kk, ]



### Biogeographical Regions (EEA - 2016) ####

if(WhereAmI == "mac"){
  biogeoregions <- read_sf(dsn = "/Users/xavi_rp/Documents/D5_FFGRCC/BiogeographicalRegions2016/eea_v_3035_1_mio_biogeo-regions_p_2016_v01_r00", 
                           layer = "BiogeoRegions2016")
}else{
  biogeoregions <- read_sf(dsn = "/eos/jeodpp/home/users/rotllxa/BiogeographicalRegions/eea_v_3035_1_mio_biogeo-regions_p_2016_v01_r00/", 
                           layer = "BiogeoRegions2016")
}

biogeoregions <- st_transform(biogeoregions, st_crs(fbi_eu_grid))
biogeoregions



## Extracting biogeographical regions  ####

i <- st_intersection(fbi_eu_grid_kk, biogeoregions)
i
# this creates new features for squares in fbi_eu_grid which partially belong to several bioregions
# For simplicity, we'll use the centroid of the 10km "cells"

ii <- st_intersection(fbi_eu_grid, biogeoregions)
ii

sum(is.na(ii$name))  # 
sum(!is.na(ii$name)) # 



## FB indicator centroids  ####

fbi_eu_grid_centroid <- st_centroid(fbi_eu_grid)
fbi_eu_grid_centroid


## Extracting biogeographical regions for centroids  ####

fbi_bioreg <- st_join(fbi_eu_grid_centroid, biogeoregions, left = TRUE)
fbi_bioreg

sum(is.na(fbi_bioreg$name))   #  3254 these centroids might be outside boundaries, etc. (6.74% of total cells)
sum(!is.na(fbi_bioreg$name))  # 44999


table(fbi_bioreg$name)  # Available points by biogeographical region

#       Alpine Bio-geographical Region      Atlantic Bio-geographical Region     Black Sea Bio-geographical Region 
#                       4166                                  7845                                   122 
#       Boreal Bio-geographical Region   Continental Bio-geographical Region  Macaronesian Bio-geographical Region 
#                       8706                                 13222                                   109 
#       Mediterranean Bio-geographical Region     Pannonian Bio-geographical Region       Steppic Bio-geographical Region 
#                       9099                                  1327                                   403 




## Averaging indicators and threshold for good condition by bioregion  ####

fbi_bioreg_dt <- data.table(fbi_bioreg)

GoodCond_thresholds <- fbi_bioreg_dt %>% 
  group_by(name) %>%
  summarize(N = n(),
            Mean_FBIMAP = mean(FBIMAP, na.rm = TRUE), Mean_FBCONSMAP = mean(FBCONSMAP, na.rm = TRUE),
            Max_FBIMAP = max(FBIMAP, na.rm = TRUE), Max_FBCONSMAP = max(FBCONSMAP, na.rm = TRUE),          # Max 
            Max60_FBIMAP = round((max(FBIMAP, na.rm = TRUE) * 0.6), 0), Max60_FBCONSMAP = round((max(FBCONSMAP, na.rm = TRUE) * 0.6), 0),     # threshold for good condition = 60%
            Min_FBIMAP = min(FBIMAP, na.rm = TRUE), Min_FBCONSMAP = min(FBCONSMAP, na.rm = TRUE))

GoodCond_thresholds
#   name                                      N Mean_FBIMAP Mean_FBCONSMAP Max_FBIMAP Max_FBCONSMAP Max60_FBIMAP Max60_FBCONSMAP
# 1 Alpine Bio-geographical Region         4166       13.5           12.8          31            26           19              16
# 2 Atlantic Bio-geographical Region       7845       19.5           15.9          33            27           20              16
# 3 Black Sea Bio-geographical Region       122       27.1           23.7          31            27           19              16
# 4 Boreal Bio-geographical Region         8706       11.7           12.0          24            23           14              14
# 5 Continental Bio-geographical Region   13222       22.8           20.7          32            29           19              17
# 6 Macaronesian Bio-geographical Region    109        7.11           6.17         12            10            7               6
# 7 Mediterranean Bio-geographical Region  9099       25.1           19.9          34            28           20              17
# 8 Pannonian Bio-geographical Region      1327       25.6           24.5          31            28           19              17
# 9 Steppic Bio-geographical Region         403       29.4           27.3          31            29           19              17
# 10 NA                                    3254       17.4           14.2          31            26           19              16



GoodCond_thresholds[, c("name", "Max60_FBIMAP")]
GoodCond_thresholds[, c("name", "Max60_FBCONSMAP")]



unique(fbi_bioreg_dt$name)
summary(fbi_bioreg_dt[name == "Alpine Bio-geographical Region", FBIMAP])
summary(fbi_bioreg_dt[name == "Mediterranean Bio-geographical Region", FBIMAP])
summary(fbi_bioreg_dt[name == "Continental Bio-geographical Region", FBIMAP])
summary(fbi_bioreg_dt[name == "Atlantic Bio-geographical Region", FBIMAP])

ggplot(fbi_bioreg_dt, aes(x = FBIMAP)) + 
  geom_bar() + 
  facet_wrap(~name)




summary(fbi_bioreg_dt[name == "Alpine Bio-geographical Region", FBCONSMAP])
summary(fbi_bioreg_dt[name == "Mediterranean Bio-geographical Region", FBCONSMAP])
summary(fbi_bioreg_dt[name == "Continental Bio-geographical Region", FBCONSMAP])
summary(fbi_bioreg_dt[name == "Atlantic Bio-geographical Region", FBCONSMAP])


ggplot(fbi_bioreg_dt, aes(x = FBCONSMAP)) + 
  geom_bar() + 
  facet_wrap(~name)



## Calculating good condition for each cell  ####

GoodCond_cells <- fbi_bioreg_dt %>%
  left_join(GoodCond_thresholds[c("name", "Max60_FBIMAP", "Max60_FBCONSMAP")], by = "name") %>%
  mutate(Condition_FBIMAP = if_else(FBIMAP < Max60_FBIMAP, "Not-Good", "Good")) %>%
  mutate(Condition_FBCONSMAP = if_else(FBCONSMAP < Max60_FBCONSMAP, "Not-Good", "Good")) %>%
  group_by(name) %>%
  summarize( N = n(), 
             Good_Cond_FBIMAP = sum(Condition_FBIMAP == "Good", na.rm = TRUE),
             Good_Cond_FBCONSMAP = sum(Condition_FBCONSMAP == "Good", na.rm = TRUE),
             Good_Cond_FBIMAP_perc = (sum(Condition_FBIMAP == "Good", na.rm = TRUE) * 100) / n(),
             Good_Cond_FBCONSMAP_perc = (sum(Condition_FBCONSMAP == "Good", na.rm = TRUE) * 100) / n())

GoodCond_cells
#    name                                      N Good_Cond_FBIMAP Good_Cond_FBCONSMAP Good_Cond_FBIMAP_perc Good_Cond_FBCONSMAP_perc
#  1 Alpine Bio-geographical Region         4166             1344                1575                  32.3                     37.8
#  2 Atlantic Bio-geographical Region       7845             3612                4581                  46.0                     58.4
#  3 Black Sea Bio-geographical Region       122              121                 121                  99.2                     99.2
#  4 Boreal Bio-geographical Region         8706             2892                2875                  33.2                     33.0
#  5 Continental Bio-geographical Region   13222            11600               11562                  87.7                     87.4
#  6 Macaronesian Bio-geographical Region    109               76                  65                  69.7                     59.6
#  7 Mediterranean Bio-geographical Region  9099             8185                7620                  90.0                     83.7
#  8 Pannonian Bio-geographical Region      1327             1327                1327                 100                      100  
#  9 Steppic Bio-geographical Region         403              403                 403                 100                      100  
#  10 NA                                    3254             1322                1379                  40.6                     42.4



## Mapping Good/Not-Good condition  ####

GoodCond_map <- fbi_bioreg %>%
  left_join(GoodCond_thresholds[c("name", "Max60_FBIMAP", "Max60_FBCONSMAP")], by = "name")  %>% 
  mutate(Condition_FBIMAP = if_else(FBIMAP < Max60_FBIMAP, "Not-Good", "Good")) %>%
  mutate(Condition_FBCONSMAP = if_else(FBCONSMAP < Max60_FBCONSMAP, "Not-Good", "Good"))
  
GoodCond_map
names(GoodCond_map)

st_write(GoodCond_map, "GoodCond_map.shp")
 
gg


## Condition Indicator (from Vallecillo et al, 2022. doi:10.2760/13048) ####

# Condition_Indicator = (V-VL) / (VH-VL)
# V: Observed value of the variable
# VL: Low condition value (lower reference level)
# VH: High condition value (upper reference level)

# But we might want to establish 60% of max value (or any other reference level) as 1, and minimum value as 0


fbi_bioreg_dt

#ConditionIndicator <- fbi_bioreg_dt %>%
#  left_join(GoodCond_thresholds[c("name", "Max_FBIMAP", "Max_FBCONSMAP", "Min_FBIMAP", "Min_FBCONSMAP")], by = "name") %>%
#  mutate(ConditionIndicator_FBIMAP = round((FBIMAP - Min_FBIMAP) / (Max_FBIMAP - Min_FBIMAP), 2)) %>%
#  mutate(ConditionIndicator_FBCONSMAP = round((FBCONSMAP - Min_FBCONSMAP) / (Max60_FBCONSMAP - Min_FBCONSMAP), 2)) #%>% head()


ConditionIndicator <- fbi_bioreg_dt %>%
  left_join(GoodCond_thresholds[c("name", "Max60_FBIMAP", "Max60_FBCONSMAP", "Min_FBIMAP", "Min_FBCONSMAP")], by = "name") %>%
  mutate(ConditionIndicator_FBIMAP = ifelse(FBIMAP >= Max60_FBIMAP, 
                                            1, 
                                            round((FBIMAP - Min_FBIMAP) / (Max60_FBIMAP - Min_FBIMAP), 2))) %>%
  mutate(ConditionIndicator_FBCONSMAP = ifelse(FBCONSMAP >= Max60_FBCONSMAP, 
                                               1, 
                                               round((FBCONSMAP - Min_FBCONSMAP) / (Max60_FBCONSMAP - Min_FBCONSMAP), 2))) #%>% head()

ConditionIndicator

# some checks
range(ConditionIndicator$ConditionIndicator_FBIMAP)
summary(ConditionIndicator$ConditionIndicator_FBIMAP)
ConditionIndicator[ConditionIndicator_FBIMAP == 1][1001]
ConditionIndicator[ConditionIndicator_FBIMAP != 1][1001]



# statistics by biogeographical region

ConditionIndicator %>% 
  group_by(name) %>%
  summarise(mean_indicator_FBIMAP = mean(ConditionIndicator_FBIMAP),
            sd_indicator_FBIMAP = sd(ConditionIndicator_FBIMAP),
            mean_indicator_FBCONSMAP = mean(ConditionIndicator_FBCONSMAP),
            sd_indicator_FBCONSMAP = sd(ConditionIndicator_FBCONSMAP))

#     name                                  mean_indicator_FBIMAP sd_indicator_FBIMAP mean_indicator_FBCONSMAP sd_indicator_FBCONSMAP
#  1 Alpine Bio-geographical Region                        0.617              0.343                     0.652                 0.339 
#  2 Atlantic Bio-geographical Region                      0.824              0.229                     0.827                 0.250 
#  3 Black Sea Bio-geographical Region                     0.992              0.0905                    0.992                 0.0905
#  4 Boreal Bio-geographical Region                        0.660              0.313                     0.672                 0.306 
#  5 Continental Bio-geographical Region                   0.965              0.113                     0.958                 0.133 
#  6 Macaronesian Bio-geographical Region                  0.777              0.376                     0.794                 0.292 
#  7 Mediterranean Bio-geographical Region                 0.987              0.0483                    0.976                 0.0684
#  8 Pannonian Bio-geographical Region                     1                  0                         1                     0     
#  9 Steppic Bio-geographical Region                       1                  0                         1                     0     
#  10 NA                                                   0.844              0.194                     0.830                 0.215 


## Mapping Condition Indicator  ####

CondIndic_map <- fbi_bioreg %>%
  left_join(GoodCond_thresholds[c("name", "Max60_FBIMAP", "Max60_FBCONSMAP", "Min_FBIMAP", "Min_FBCONSMAP")], by = "name") %>%
  mutate(ConditionIndicator_FBIMAP = ifelse(FBIMAP >= Max60_FBIMAP, 
                                            1, 
                                            round((FBIMAP - Min_FBIMAP) / (Max60_FBIMAP - Min_FBIMAP), 2))) %>%
  mutate(ConditionIndicator_FBCONSMAP = ifelse(FBCONSMAP >= Max60_FBCONSMAP, 
                                               1, 
                                               round((FBCONSMAP - Min_FBCONSMAP) / (Max60_FBCONSMAP - Min_FBCONSMAP), 2)))


CondIndic_map
range(CondIndic_map$ConditionIndicator_FBIMAP)
plot(CondIndic_map["FBIMAP"])

st_write(CondIndic_map, "CondIndic_map.shp")






## Crop type (Rega et al. 2020) ####

list.files("/eos/jeodpp/data/projects/REFOCUS/data/BIODIVERSITY/Rega/")

Crop_management_systems_dom50_def <- rast("/eos/jeodpp/data/projects/REFOCUS/data/BIODIVERSITY/Rega/Crop_management_systems_dom50_def.tif")
Crop_management_systems_dom50_def
Crop_management_systems_dom50_def_vals <- values(Crop_management_systems_dom50_def, dataframe = TRUE)
Crop_management_systems_dom50_def_unique <- unique(Crop_management_systems_dom50_def_vals$Cropmgmt50)
sort(levels(Crop_management_systems_dom50_def_unique))

rcl_df <- data.frame(is = sort(levels(Crop_management_systems_dom50_def_unique)),
                     becomes = sort(levels(Crop_management_systems_dom50_def_unique)))

rcl_df$becomes <- gsub(" - High", "", rcl_df$becomes)
rcl_df$becomes <- gsub(" - Low" , "", rcl_df$becomes)
rcl_df$becomes <- gsub(" - Medium" , "", rcl_df$becomes)
rcl_df$becomes <- gsub("Grasslands and meadows" , "Grasslands_meadows", rcl_df$becomes)
rcl_df$becomes <- gsub("Mixed systems with prevalence of arable crops" , "Mixed_Prevalence_Arable_Crops", rcl_df$becomes)
rcl_df$becomes <- gsub("Mixed systems with prevalence of grasslands" , "Mixed_Prevalence_Grasslands", rcl_df$becomes)
rcl_df$becomes <- gsub("Mixed systems with prevalence of permanent crops" , "Mixed_Prevalence_Permanent_Crops", rcl_df$becomes)
rcl_df$becomes <- gsub("No Data" , "No_Data", rcl_df$becomes)
rcl_df$becomes <- gsub("Non-agricultural areas" , "Non-agricultural", rcl_df$becomes)
rcl_df$becomes <- gsub("Specialist field crops - cereals" , "Specialist_Cereals", rcl_df$becomes)
rcl_df$becomes <- gsub("Specialist field crops - industrial crops" , "Specialist_Industrial_Crops", rcl_df$becomes)
rcl_df$becomes <- gsub("Specialist Forage crops" , "Specialist_Forage", rcl_df$becomes)
rcl_df$becomes <- gsub("Specialist fruits and citrus fruits" , "Specialist_Fruits_Citrus", rcl_df$becomes)
rcl_df$becomes <- gsub("Specialist Olives" , "Specialist_Olives", rcl_df$becomes)
rcl_df$becomes <- gsub("Specialist Vegetables, flowers and horticulture" , "Specialist_Horticulture", rcl_df$becomes)
rcl_df$becomes <- gsub("Specialist Vineyards" , "Specialist_Vineyards", rcl_df$becomes)
rcl_df

#rcl_df$becomes_num <- 1:nrow(rcl_df)
rcl_df
head(rcl_df)
View(rcl_df)


#Crop_systems <- classify(Crop_management_systems_dom50_def, 
#                         rcl = rcl_df[, c(1, 3)], 
#                         filename = "crop_systems_from_Rega.tif")

head(Crop_management_systems_dom50_def_vals)
tail(Crop_management_systems_dom50_def_vals)
length(Crop_management_systems_dom50_def_vals$Cropmgmt50)


Crop_management_systems_dom50_def_vals <- as.data.table(Crop_management_systems_dom50_def_vals)
Crop_management_systems_dom50_def_vals <- left_join(Crop_management_systems_dom50_def_vals, 
                                                    rcl_df, 
                                                    by = join_by(Cropmgmt50 == is))
head(Crop_management_systems_dom50_def_vals)
Crop_management_systems_dom50_def_vals

setnames(Crop_management_systems_dom50_def_vals, "becomes", "Crop_System")
head(Crop_management_systems_dom50_def_vals)

Crop_systems <- setValues(Crop_management_systems_dom50_def, 
                          Crop_management_systems_dom50_def_vals[, .SD, .SDcols = "Crop_System"],
                          keepnames = TRUE)
Crop_systems
names(Crop_systems) <- "Crop_System"

#writeRaster(Crop_systems, filename = "crop_systems_from_Rega.tif", overwrite = TRUE)


## 'Crop_Systems' is at 100m grid. It needs to be resampled to 'fbi_eu_grid' (10km)
fbi_eu_grid

Crop_systems <- project(Crop_systems, crs(fbi_eu_grid))
Crop_systems

writeRaster(Crop_systems, filename = "crop_systems_from_Rega.tif", overwrite = TRUE)
#Crop_systems <- rast("crop_systems_from_Rega.tif")



## aggregate to 10km

#Crop_systems_vals <- values(Crop_systems, dataframe = TRUE)  

rcl_df1 <- data.frame(rcl_df[, "becomes"])
rcl_df1 <- data.frame(rcl_df1[!duplicated(rcl_df1), ])
rcl_df1$code <- 1:nrow(rcl_df1)
names(rcl_df1) <- c("Crop_System", "Crop_System_Code")
head(rcl_df1)
View(rcl_df1)
rcl_df1

Crop_systems_dt <- as.data.table(values(Crop_systems, dataframe = TRUE))
Crop_systems_dt

Crop_systems_num <- left_join(Crop_systems_dt, rcl_df1, by = join_by(Crop_System)) 
Crop_systems_num

Crop_systems[["Crop_systems_num"]] <- Crop_systems_num$Crop_System_Code
Crop_systems

writeRaster(Crop_systems[["Crop_systems_num"]], filename = "crop_systems_num_from_Rega.tif", overwrite = TRUE)
#Crop_systems_num <- rast("crop_systems_num_from_Rega.tif")
Crop_systems_num
unique(values(Crop_systems_num))


mode_fun <- function(x){
  mode_res_1 <- table(x)
  mode_res_max <- max(mode_res_1, na.rm = TRUE)
  if(mode_res_max >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
    mode_res <- names(which.max(mode_res_1))
  }else{
    if(sum(mode_res_1[names(mode_res_1) %in% c("Mixed_Prevalence_Arable_Crops",
                                           "Specialist_Cereals",
                                           "Specialist_Forage",
                                           "Specialist_Horticulture",
                                           "Specialist_Industrial_Crops")]) >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
      mode_res <- "Mixed_Prevalence_Arable_Crops"
    }else if(sum(mode_res_1[names(mode_res_1) %in% c("Mixed_Prevalence_Permanent_Crops",
                                                     "Specialist_Fruits_Citrus",
                                                     "Specialist_Olives",
                                                     "Specialist_Vineyards")]) >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
      mode_res <- "Mixed_Prevalence_Permanent_Crops"
    }else if(sum(mode_res_1[names(mode_res_1) %in% c("Mixed_Prevalence_Grasslands",
                                                    "Grasslands_meadows")]) >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
      mode_res <- "Mixed_Prevalence_Grasslands"
    }else{
      mode_res <- "Non-agricultural"
    }
  }

  return(mode_res)
}

mode_fun_num <- function(x){
  if(all(is.na(x))){
    mode_res <- NA
  }else{
    mode_res_1 <- table(x)
    mode_res_max <- max(mode_res_1, na.rm = TRUE)
    if(mode_res_max >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
      mode_res <- as.numeric(names(which.max(mode_res_1)))
    }else{
      if(sum(mode_res_1[names(mode_res_1) %in% c(2,
                                                 7,
                                                 9,
                                                 12,
                                                 8)]) >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
        mode_res <- 2
      }else if(sum(mode_res_1[names(mode_res_1) %in% c(4,
                                                       10,
                                                       11,
                                                       13)]) >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
        mode_res <- 4
      }else if(sum(mode_res_1[names(mode_res_1) %in% c(3,
                                                       1)]) >= (sum(mode_res_1, na.rm = TRUE) * 0.5)){
        mode_res <- 3
      }else{
        mode_res <- 6
      }
    }
    
  }
  
  return(mode_res)

}


#
#Crop_systems
#cat_coords <-  c(3500000, 3800000, 1900000, 2300000)   # Catalonia (LAEA, m) (xmin, xmax, ymin, ymax)
#Crop_systems_cat <- crop(Crop_systems, ext(cat_coords))

#Crop_systems_num <- merge(as.data.table(Crop_systems_cat), rcl_df1, by = "Crop_System", all.x = TRUE)
#Crop_systems_num

#Crop_systems_cat[["Crop_systems_num"]] <- Crop_systems_num$Crop_System_Code
#Crop_systems_cat

#
t0 <- Sys.time()
Crop_systems_10km <- aggregate(Crop_systems_num,
                               fact = 100,
                               fun = mode_fun_num,
                               cores = 1) #, filename="", overwrite=FALSE, wopt=list())
Crop_systems_10km
Sys.time() - t0

writeRaster(Crop_systems_10km, filename = "crop_systems_from_Rega_10km.tif", overwrite = TRUE)
#Crop_systems_10km <- rast("crop_systems_from_Rega_10km.tif")



# checks

plot(Crop_systems_10km)
s <- sel(Crop_systems_10km)
s
plot(s)

s1 <- sel(s)
s1
plot(s1)

s2 <- crop(Crop_systems, s1)
plot(s2[["Crop_systems_num"]])
table(values(s2[["Crop_systems_num"]]))



## Aggregation to 1km

t0 <- Sys.time()
Crop_systems_1km <- aggregate(Crop_systems_num,
                              fact = 10,
                              fun = mode_fun_num,
                              cores = 10) #, filename="", overwrite=FALSE, wopt=list())
Crop_systems_1km
Sys.time() - t0

writeRaster(Crop_systems_1km, filename = "crop_systems_from_Rega_1km.tif", overwrite = TRUE)
#Crop_systems_1km <- rast("crop_systems_from_Rega_1km.tif")
plot(Crop_systems_1km)


rcl_df1_char <- data.frame(is = c(2,
                                  7,
                                  9,
                                  12,
                                  8,
                                  4,
                                  10,
                                  11,
                                  13,
                                  3,
                                  1,
                                  5), 
                           becomes = c("Mixed_Prevalence_Arable_Crops",
                                       "Specialist_Cereals",
                                       "Specialist_Forage",
                                       "Specialist_Horticulture",
                                       "Specialist_Industrial_Crops",
                                       "Mixed_Prevalence_Permanent_Crops",
                                       "Specialist_Fruits_Citrus",
                                       "Specialist_Olives",
                                       "Specialist_Vineyards",
                                       "Mixed_Prevalence_Grasslands",
                                       "Grasslands_meadows",
                                       "No_Data")) %>% arrange(is)
  

rcl_df1_char

Crop_systems_1km_vals <- data.table(values(Crop_systems_1km))
Crop_systems_1km_vals <- left_join(Crop_systems_1km_vals, rcl_df1_char,
                                   by = join_by(Crop_systems_num == is))
Crop_systems_1km_vals

Crop_systems_1km_char <- Crop_systems_1km
Crop_systems_1km_char[["code_char"]] <- Crop_systems_1km_vals$becomes
Crop_systems_1km_char

plot(Crop_systems_1km_char[["code_char"]])

writeRaster(Crop_systems_1km_char[["code_char"]], filename = "crop_systems_from_Rega_1km_char.tif", overwrite = TRUE)
#Crop_systems_1km_char <- rast("crop_systems_from_Rega_1km_char.tif")



## so far, nothing worked to define meaningful reference sites

## FB Change Map to define reference sites ####
# Reference sites defined as cells with positive or stable change and with a "baseline" sufficiently high
# As we don't have the separated values used to calculate the change (1980s vs 2010s), and we don't have the 
# baseline, we use the 2010s combined with the change to define the reference sites as follows:
# if:
#   change < 0  (negative tendency)                    --->   bad (not a reference site)
#   change = 0  (stable)         --->   FBI_2010s low  --->   Not a reference
#                                --->   FBI_2010s high --->   Reference
#   change > 0  (low positive)   --->   FBI_2010s low  --->   Not a reference
#                                --->   FBI_2010s high --->   Reference
#   change >> 0 (high poritive)  --->   FBI_2010s low  --->   Not a reference
#                                --->   FBI_2010s high --->   Reference ??
#  
#  FBI_2010s to be high or low can be defined by a high percentil (e.g. 90 or 95)


if(WhereAmI == "mac"){
  list.files("./Study2_Change/FBChange_maps_10km")
}else{
  list.files("./Map_Farmland_Bird_Indicators_Study2Change/FBChange_maps_10km")
}


if(WhereAmI == "mac"){
  fb_change_eu_grid <- read_sf("./Study2_Change/FBChange_maps_10km/fb_change_maps_10_km.shp")
}else{
  fb_change_eu_grid <- read_sf("./Map_Farmland_Bird_Indicators_Study2Change/FBChange_maps_10km/fb_change_maps_10_km.shp")
}

fb_change_eu_grid

# FBAGR       Indicator based on all farmland species according to the classification done in EBBA2
#             for Agricultural / grassland birds (Keller et al. 2020)
# FBI         Indicator based on farmland species included in the Farmland Bid Index
# FBCONSALL   Indicator based on all farmland species that are considered of conservation concern
#             (SPEC) in Europe by BirdLife International (2017)
# FBCONS10    Indicator based on farmland species that are considered of conservation concern
#             (SPEC) in Europe by BirdLife International (2017) and were used to produce FCONS 10-km map
# FBMEAN      Average of the previous indicators
#
# Changes between 1980s and 2010s. 

fb_change_eu_grid_vals <- data.table(fb_change_eu_grid)
fb_change_eu_grid_vals

summary(fb_change_eu_grid_vals$FBMEAN)
mean(fb_change_eu_grid_vals$FBMEAN)  # -4.641058
quantile(fb_change_eu_grid_vals$FBMEAN, seq(0, 1, 0.1))
#   0%       10%       20%       30%        40%       50%       60%       70%       80%       90%      100% 
# -95.80860 -27.06095 -19.06340 -13.45610  -8.61930  -4.10845   0.42140   5.26340  10.26380  16.97945  73.19500 

mean(fb_change_eu_grid_vals$FBI)  # -3.345352
quantile(fb_change_eu_grid_vals$FBI, seq(0, 1, 0.1))
#   0%       10%       20%       30%        40%       50%       60%       70%       80%       90%      100% 
# -99.16000 -23.86820 -16.29200 -10.49200  -6.30800  -2.14400   0.78400   5.00800   9.78400  16.56735  70.00000

ggplot(fb_change_eu_grid_vals, aes(x = FBI)) + 
  geom_histogram(bins = 10) +
  stat_bin(aes(y = after_stat(count), label = after_stat(count)), geom  = "text", bins = 10, vjust= - 0.5) 



# scatterplot: cells change (1980s:2010s) vs. fbi (2010s)

fbi_eu_grid_vals <- data.table(fbi_eu_grid)

sum(fbi_eu_grid_vals$CELLCODE %in% fb_change_eu_grid_vals$CELLCODE)  # 39050
sum(!fbi_eu_grid_vals$CELLCODE %in% fb_change_eu_grid_vals$CELLCODE) #  9203

sum(fbi_eu_grid_vals$CELLCODE %in% fbi_bioreg$CELLCODE)
sum(!fbi_eu_grid_vals$CELLCODE %in% fbi_bioreg$CELLCODE)

change_fbi <- left_join(fb_change_eu_grid_vals, 
                        #fbi_eu_grid_vals,
                        fbi_bioreg,   # including also biogeographical region 
                        by = "CELLCODE") %>% drop_na()

change_fbi
names(change_fbi)

ggplot(change_fbi, aes(x = FBIMAP, y = FBI)) + 
  geom_point() + 
  xlab("FBI 2010s") +
  ylab("FBI change (1980s:2010s)") 


perctl <- 90
perctl <- 95


change_fbi <- change_fbi %>%
  group_by(name) %>%  # by biogeographical region
  mutate(percentile = quantile(FBIMAP, perctl/100)) %>% 
  mutate(reference_site = ifelse((FBIMAP >= percentile & FBI >= 0), TRUE, FALSE)) #%>% View()

change_fbi
names(change_fbi)
unique(change_fbi$reference_site)


p_p95 <- ggplot(change_fbi, aes(x = FBIMAP, y = FBI, colour = reference_site)) + 
  geom_point(size = 1) + 
  xlab("FBI 2010s") +
  ylab("FBI change (1980s:2010s)") + 
  labs(title = "Reference sites",
       color = paste0("Change >= 0\nand percentile ", perctl)) +
  facet_wrap(~ short_name, nrow = 4, ncol = 3) 


#jpeg("Reference_sites_perc90.jpg", height = 15, width = 18, units = "cm", res = 300, pointsize = 8)
p_p90
#dev.off()

jpeg("Reference_sites_perc95.jpg", height = 15, width = 18, units = "cm", res = 300, pointsize = 8)
p_p95
dev.off()


# map

change_fbi
fb_change_eu_grid
names(fb_change_eu_grid)

fb_change_eu_grid_1 <- left_join(fb_change_eu_grid, 
                                 select(change_fbi, CELLCODE, reference_site),   
                                 by = "CELLCODE")

m_p95 <- ggplot() +
  #geom_sf(data = fb_change_eu_grid, mapping = aes(fill = FBI), color = NA)
  #geom_sf(data = fbi_eu_grid, mapping = aes(fill = FBIMAP), color = NA)
  geom_sf(data = fb_change_eu_grid_1,
          mapping = aes(fill = reference_site),
          color = NA) +
  #geom_sf(data = biogeoregions, fill = NA) +
  labs(title = "Reference sites",
       fill = paste0("Change >= 0\nand percentile ", perctl))


jpeg("Reference_sites_perc95_map.jpg", height = 15, width = 18, units = "cm", res = 300, pointsize = 8)
m_p95
dev.off()

#jpeg("Reference_sites_perc90_map.jpg", height = 15, width = 18, units = "cm", res = 300, pointsize = 8)
m_p90
#dev.off()



## Just to check cells in the highest percentile, but with negative tendency of change
change_fbi <- left_join(fb_change_eu_grid_vals, 
                        #fbi_eu_grid_vals,
                        fbi_bioreg,   # including also biogeographical region 
                        by = "CELLCODE") %>% drop_na()

change_fbi <- change_fbi %>%
  group_by(name) %>%
  mutate(percentile = quantile(FBIMAP, perctl/100)) %>% 
  mutate(reference_site = ifelse((FBIMAP >= percentile & FBI < 0), TRUE, FALSE)) #%>% View()

fb_change_eu_grid_1 <- left_join(fb_change_eu_grid, 
                                 select(change_fbi, CELLCODE, reference_site),   
                                 by = "CELLCODE")

ggplot() +
  geom_sf(data = fb_change_eu_grid_1,
          mapping = aes(fill = reference_site),
          color = NA) +
  labs(title = "Reference sites",
       fill = paste0("Change < 0\nand percentile ", perctl))
##


## Calculating reference levels ####
# By biogeographical regions
# Upper reference = 0.6 * max(FBIMAP)
change_fbi
names(change_fbi)
range(change_fbi$FBIMAP)  # 3 - 34
unique(change_fbi$reference_site)

GoodCond_thresholds_RefSites <- change_fbi %>%
  filter(reference_site == TRUE) %>%
  group_by(name) %>%
  summarize(N = n(),
            Mean_FBIMAP = mean(FBIMAP, na.rm = TRUE), #Mean_FBCONSMAP = mean(FBCONSMAP, na.rm = TRUE),
            Max_FBIMAP = max(FBIMAP, na.rm = TRUE), #Max_FBCONSMAP = max(FBCONSMAP, na.rm = TRUE),          # Max 
            Max60_FBIMAP = round((max(FBIMAP, na.rm = TRUE) * 0.6), 0), 
            Max80_FBIMAP = round((max(FBIMAP, na.rm = TRUE) * 0.8), 0), 
            Max70_FBIMAP = round((max(FBIMAP, na.rm = TRUE) * 0.7), 0), 
            Min_FBIMAP_RefSites = min(FBIMAP, na.rm = TRUE))

  
GoodCond_thresholds_RefSites
#   name                                      N  Mean_FBIMAP  Max_FBIMAP  Max60_FBIMAP  Min_FBIMAP_RefSites
# 1 Alpine Bio-geographical Region          141         24.9          31            19          24
# 2 Atlantic Bio-geographical Region        325         30.8          33            20          30
# 3 Black Sea Bio-geographical Region        17         30.1          31            19          30
# 4 Boreal Bio-geographical Region           62         21.7          23            14          21
# 5 Continental Bio-geographical Region     426         29.8          32            19          29
# 6 Mediterranean Bio-geographical Region   196         31.7          34            20          31
# 7 Pannonian Bio-geographical Region       157         28.2          31            19          28
# 8 Steppic Bio-geographical Region          35         31            31            19          31


GoodCond_thresholds_RefSites_MinBGRegion <- fbi_bioreg %>%
  data.table() %>%
  group_by(name) %>%
  summarize(Min_FBIMAP = min(FBIMAP, na.rm = TRUE))
  
GoodCond_thresholds_RefSites <- GoodCond_thresholds_RefSites %>%
  left_join(GoodCond_thresholds_RefSites_MinBGRegion, by = "name")

GoodCond_thresholds_RefSites

## Mapping good condition 

names(fb_change_eu_grid_1)
names(fbi_bioreg)
plot(fb_change_eu_grid_1["FBI"])
plot(fbi_bioreg["FBIMAP"])

GoodCond_RefSites_map_60 <- fbi_eu_grid %>%
  st_join(select(fbi_bioreg, CELLCODE, name), by = "CELLCODE") %>% #names()
  left_join(GoodCond_thresholds_RefSites[c("name", "Max60_FBIMAP")], by = "name")  %>% 
  mutate(Condition_FBIMAP = if_else(FBIMAP < Max60_FBIMAP, "Not-Good", "Good")) 

GoodCond_RefSites_map_70 <- fbi_eu_grid %>%
  st_join(select(fbi_bioreg, CELLCODE, name), by = "CELLCODE") %>% #names()
  left_join(GoodCond_thresholds_RefSites[c("name", "Max70_FBIMAP")], by = "name")  %>% 
  mutate(Condition_FBIMAP = if_else(FBIMAP < Max70_FBIMAP, "Not-Good", "Good")) 

GoodCond_RefSites_map_80 <- fbi_eu_grid %>%
  st_join(select(fbi_bioreg, CELLCODE, name), by = "CELLCODE") %>% #names()
  left_join(GoodCond_thresholds_RefSites[c("name", "Max80_FBIMAP")], by = "name")  %>% 
  mutate(Condition_FBIMAP = if_else(FBIMAP < Max80_FBIMAP, "Not-Good", "Good")) 


GoodCond_RefSites_map_60
names(GoodCond_RefSites_map_60)
unique(GoodCond_RefSites_map_60$Condition_FBIMAP)


p_60 <- ggplot(data = GoodCond_RefSites_map_60) +
  geom_sf(mapping = aes(color = Condition_FBIMAP, fill = Condition_FBIMAP)) + 
  labs(title = "Threshold for good cond.: Max * 0.6") +
  theme(plot.title = element_text(hjust = 0.5))

p_70 <- ggplot(data = GoodCond_RefSites_map_70) +
  geom_sf(mapping = aes(color = Condition_FBIMAP, fill = Condition_FBIMAP)) +
  labs(title = "Threshold for good cond.: Max * 0.7") +
  theme(plot.title = element_text(hjust = 0.5))

p_80 <- ggplot(data = GoodCond_RefSites_map_80) +
  geom_sf(mapping = aes(color = Condition_FBIMAP, fill = Condition_FBIMAP)) +
  labs(title = "Threshold for good cond.: Max * 0.8") +
  theme(plot.title = element_text(hjust = 0.5))


jpeg("GoodCond_RefSites_map.jpg", height = 30, width = 20, units = "cm", res = 300, pointsize = 8)
gridExtra::grid.arrange(p_60, p_70, p_80)
dev.off()



## Condition Indicator (from Vallecillo et al, 2022. doi:10.2760/13048) ####

# Condition_Indicator = (V-VL) / (VH-VL)
# V: Observed value of the variable
# VL: Low condition value (lower reference level)
# VH: High condition value (upper reference level)

# But we establish 60% of max value (or any other reference level) as 1, and minimum value as 0

names(GoodCond_RefSites_map_60)
names(fbi_bioreg_dt)
names(fbi_eu_grid)

 
ConditionIndicator_60 <- fbi_eu_grid %>%
  st_join(select(fbi_bioreg, CELLCODE, name), by = "CELLCODE") %>% #names()
  left_join(GoodCond_thresholds_RefSites[c("name", "Max_FBIMAP","Max60_FBIMAP", "Min_FBIMAP")], by = "name")  %>% #View()
  mutate(ConditionIndicator_FBIMAP = ifelse(FBIMAP >= Max60_FBIMAP, 
                                            1, 
                                            #round((FBIMAP - Min_FBIMAP) / (Max_FBIMAP - Min_FBIMAP), 2))) #%>% names()
                                            round((FBIMAP - Min_FBIMAP) / (Max60_FBIMAP - Min_FBIMAP), 2))) #%>% names()

ConditionIndicator_60
names(ConditionIndicator_60)
range(ConditionIndicator_60$ConditionIndicator_FBIMAP, na.rm = TRUE)
sort(unique(ConditionIndicator_60$ConditionIndicator_FBIMAP, na.rm = TRUE))
View(select(ConditionIndicator_60, name, FBIMAP, Max60_FBIMAP, Max_FBIMAP, Min_FBIMAP, ConditionIndicator_FBIMAP))  
View(ConditionIndicator_60)


p60_indicator_adjusted <- ggplot(data = ConditionIndicator_60) +
  geom_sf(mapping = aes(fill = ConditionIndicator_FBIMAP, color = ConditionIndicator_FBIMAP)) + 
  #scale_colour_gradientn(colors = terrain.colors(20)) +
  #scale_fill_gradientn(colors = terrain.colors(20)) +
  viridis::scale_fill_viridis(discrete = FALSE, direction = -1, name = "Condition indicator\n(adjusted)") +
  viridis::scale_colour_viridis(discrete = FALSE, direction = -1, name = "Condition indicator\n(adjusted)") +
  labs(title = "Condition Indicator",
       subtitle = "Adjusted scale (0 - 0.6 adjusted to 0 - 1)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

jpeg("CondIndicator_RefSites_map_adjusted.jpg", height = 15, width = 20, units = "cm", res = 300, pointsize = 8)
p60_indicator_adjusted
dev.off()



## not adjusting scale
ConditionIndicator_60_1 <- fbi_eu_grid %>%
  st_join(select(fbi_bioreg, CELLCODE, name), by = "CELLCODE") %>% #names()
  left_join(GoodCond_thresholds_RefSites[c("name", "Max_FBIMAP","Max60_FBIMAP", "Min_FBIMAP")], by = "name")  %>% #View()
  mutate(ConditionIndicator_FBIMAP = round((FBIMAP - Min_FBIMAP) / (Max_FBIMAP - Min_FBIMAP), 2)) #%>% names()

range(ConditionIndicator_60_1$ConditionIndicator_FBIMAP, na.rm = TRUE)
sum(ConditionIndicator_60_1$ConditionIndicator_FBIMAP > 1, na.rm = TRUE)
View(ConditionIndicator_60_1[ConditionIndicator_60_1$ConditionIndicator_FBIMAP > 1 &
                               !is.na(ConditionIndicator_60_1$ConditionIndicator_FBIMAP), ])

ConditionIndicator_60_1[ConditionIndicator_60_1$ConditionIndicator_FBIMAP > 1 &
                          !is.na(ConditionIndicator_60_1$ConditionIndicator_FBIMAP), "ConditionIndicator_FBIMAP"] <- 1
range(ConditionIndicator_60_1$ConditionIndicator_FBIMAP, na.rm = TRUE)


p60_indicator <- ggplot(data = ConditionIndicator_60_1) +
  geom_sf(mapping = aes(fill = ConditionIndicator_FBIMAP, color = ConditionIndicator_FBIMAP)) + 
  #scale_colour_gradientn(colors = terrain.colors(20)) +
  #scale_fill_gradientn(colors = terrain.colors(20)) +
  viridis::scale_fill_viridis(discrete = FALSE, direction = -1, name = "Condition indicator\n(not adjusted)") +
  viridis::scale_colour_viridis(discrete = FALSE, direction = -1, name = "Condition indicator\n(not adjusted)") +
  labs(title = "Condition Indicator",
       subtitle = "(Not adjusted)") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

p60_indicator

jpeg("CondIndicator_RefSites_map_NotAdjusted.jpg", height = 15, width = 20, units = "cm", res = 300, pointsize = 8)
p60_indicator
dev.off()
