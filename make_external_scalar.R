# Alberto Rovellini
# 3/7/2023
# Code to write external scaling input for plankton growth rates
# Take ROMS-derived monthly indices from NEP10K ROMS-NPZ
# Calculate 2000-2012 climatologies
# Calculate 2013-2016 climatology
# Get ratio 
# Pack to NetCDF

# See https://confluence.csiro.au/display/Atlantis/External+Mortality%2C+Growth+and+Recruitment+Scaling

# this exercise is hinged on the following mapping from NPZ to Atlantis GOA:
# time-averaged small coastal copepod concentration <- Mesozooplankton
# time-averaged detritus concentration <- Detritus (split lab, ref, car as 40:40:20)
# time-averaged euphausiid concentration <- Euphausiids
# time-averaged large microzooplankton concentration <- Microzooplankton
# time-averaged small microzooplankton concentration <- Microzooplankton
# time-averaged neocalanus spp. concentration <- Mesozooplankton
# time-averaged ammonia concentration <- NH3
# time-averaged nitrate concentration <- NO3
# time-averaged large phytoplankton concentration <- Diatoms
# time-averaged small phytoplankton concentration <- Picophytoplankton

library(tidync)
library(dplyr)
library(tidyr)
library(ncdf4)
library(RNetCDF)
library(lubridate)
library(ggplot2)

options(digits=3)

# # use biomass or production?
# useprod <- FALSE
# 
# # read integrated quantities
# npz <- read.csv('../../../ROMS_to_Index/Data/output_from_loon/nep_sum_hind.csv')
# 
# if(useprod){
#   npz <- npz %>%
#     filter(depthclass == 'All', varname %in% c('prod_Cop','prod_Eup','prod_MZL','prod_MZS','prod_PhL','prod_PhS'))
#   npz$varname <- gsub('prod_', '', npz$varname)
# } else {
#   npz <- npz %>%
#     filter(depthclass == 'All', varname %in% c('Cop','Eup','MZL','MZS','PhL','PhS'))
# }
# 
# key <- data.frame('varname'=c('Cop','Eup','MZL','MZS','PhL','PhS'),
#                   'Code'=c('ZM',
#                          'EUP',
#                          'ZS',
#                          'ZS',
#                          'PL',
#                          'PS'))
# # join Atlantis groups
# npz <- npz %>%
#   left_join(key, by = 'varname')
# 
# npz_atlantis <- npz %>%
#   group_by(date, NMFS_AREA, Code) %>%
#   summarise(Biomass = sum(value)) %>%
#   ungroup() %>%
#   mutate(year = year(date), month = month(date)) %>%
#   rowwise() %>%
#   mutate(season = case_when(
#     month %in% 1:3 ~ 1,
#     month %in% 4:6 ~ 2,
#     month %in% 7:9 ~ 3,
#     month %in% 10:12 ~ 4
#   ))
# 
# # 2000-2012 climatologies
# # do it by season and decide what to use
# clim_2000_2012 <- npz_atlantis %>%
#   filter(year %in% 2000:2012) %>%
#   group_by(season, NMFS_AREA, Code) %>%
#   summarise(climatology = mean(Biomass)) %>%
#   ungroup()
# 
# # 2013-2016 climatology
# clim_2013_2016 <- npz_atlantis %>%
#   filter(year %in% 2013:2016) %>%
#   group_by(season, NMFS_AREA, Code) %>%
#   summarise(climatology = mean(Biomass)) %>%
#   ungroup()
# 
# # scalars
# scalars <- clim_2000_2012 %>%
#   left_join(clim_2013_2016, by = c('season','NMFS_AREA','Code')) %>%
#   mutate(scalar = climatology.y / climatology.x) %>%
#   filter(NMFS_AREA != 'All')
# 
# # plot to compare with realized output from Atlantis runs - how do these compare and do we need to tune these scalars
# plotscalar <- scalars %>%
#   group_by(season, Code) %>%
#   summarize(scalar = mean(scalar)) %>% # average across NMFS areas
#   ungroup() %>%
#   ggplot(aes(x = season, y = scalar))+
#   geom_point()+
#   geom_hline(yintercept = 1, color = 'red')+
#   theme_bw()+
#   facet_wrap(~Code)
# plotscalar
# 
# # scale the scalar further for the purpose of tuning
# # mod <- 3
# # scalars <- scalars %>%
# #   rowwise() %>%
# #   mutate(scalar = ifelse(scalar>0, scalar * mod, scalar / mod)) %>%
# #   ungroup()
#   
# # now we need to map NMFS areas to Atlantis boxes
# # for BC boxes we use values from Area 650
# box_nmfs <- rbind(data.frame('NMFS_AREA' = '610', 'box_id' = 0:22),
#                  data.frame('NMFS_AREA' = '620', 'box_id' = 23:37),
#                  data.frame('NMFS_AREA' = '630', 'box_id' = 38:60),
#                  data.frame('NMFS_AREA' = '640', 'box_id' = 61:75),
#                  data.frame('NMFS_AREA' = '650', 'box_id' = 76:108))
# 
# # add boxes to scalars
# scalars <- scalars %>%
#   full_join(box_nmfs, by = 'NMFS_AREA')
# 
# nc_name <- '../data/scalar_heatwave_from_NPZ_x3.nc'
# 
# nc_file <- create.nc(nc_name)
# 
# nbox <- 109
# nlayer <- 7
# t_units <- "seconds since 1990-01-01 00:00:00 -9" # for consistency with the init.nc files
# seconds_timestep <- 43200 * 2 * (365/4) # this is a time step per quarter
# this_geometry <- "GOA_WGS84_V4_final.bgm"
# 
# # time array
# time_array <- seq(0, seconds_timestep * (40*4), seconds_timestep) # in quarters
# ntime <- length(time_array)
# 
# # scalar array
# # dummy array with 1 in all positions, this should result in no changes in the run
# # scalar_array <- array(data = 1, dim = c(nlayer, nbox, ntime))
# 
# # and one with 30 years at 1 then 10 years at 0.8
# # scalar_array <- array(data = c(rep(matrix(1, nrow = nlayer, ncol = nbox), 10),
# #                                rep(matrix(0.8, nrow = nlayer, ncol = nbox), 5)),
# #                       dim = c(nlayer, nbox, ntime))
# 
# # for heatwave runs
# # define dimensions
# dim.def.nc(nc_file, "t", unlim=TRUE)
# dim.def.nc(nc_file, "b", nbox) # manual 
# dim.def.nc(nc_file, "z", nlayer) # manual
# 
# # define time variable
# var.def.nc(nc_file, "t", "NC_DOUBLE", "t")
# 
# # define attributes
# att.put.nc(nc_file, "t", "units", "NC_CHAR", t_units)
# att.put.nc(nc_file, "t", "dt", "NC_DOUBLE", seconds_timestep)
# att.put.nc(nc_file, "NC_GLOBAL", "title", "NC_CHAR", "External scalar for plankton growth")
# att.put.nc(nc_file, "NC_GLOBAL", "geometry", "NC_CHAR", this_geometry)
# att.put.nc(nc_file, "NC_GLOBAL", "parameters", "NC_CHAR", " ")
# 
# # enter time variable in the netcdf
# var.put.nc(nc_file, "t", time_array)
# 
# 
# pack_scalars <- function(fg){
#   
#   this_scalar <- scalars %>%
#     filter(Code == fg) %>%
#     mutate(scalar = as.numeric(round(scalar, digits = 3))) %>%
#     pull(scalar) %>%
#     rep(., each = nlayer)
#   
#   scalar_array <- array(data = c(rep(matrix(1, nrow = nlayer, ncol = nbox), 31*4), # by quarter
#                                     rep(this_scalar, 5), 
#                                     rep(matrix(1, nrow = nlayer, ncol = nbox), 5*4)), # by quarter
#                            dim = c(nlayer, nbox, ntime))
#   
#   var.def.nc(nc_file, paste0(fg, "growth_ff"), "NC_DOUBLE", c("z","b","t"))
#   
#   att.put.nc(nc_file, paste0(fg, "growth_ff"), "_FillValue", "NC_DOUBLE", 1)
#   
#   var.put.nc(nc_file, paste0(fg, "growth_ff"), scalar_array)
#   
# }
# 
# purrr::map(c('PL', 'PS', 'ZM', 'ZL', 'EUP'), pack_scalars)
# 
# close.nc(nc_file)

# Dummy scalars for trial runs --------------------------------------------

nc_name <- '../data/scalar_heatwave_PL_ZM_EUP_03_Hindcast.nc'

nc_file <- create.nc(nc_name)

nbox <- 109
nlayer <- 7
t_units <- "seconds since 1990-01-01 00:00:00 -9" # for consistency with the init.nc files
seconds_timestep <- 43200 * 2 * 365
this_geometry <- "GOA_WGS84_V4_final.bgm"

# time array
time_array <- seq(0, seconds_timestep * 70, seconds_timestep)
ntime <- length(time_array)

# scalar array
# dummy array with 1 in all positions, this should result in no changes in the run
# scalar_array <- array(data = 1, dim = c(nlayer, nbox, ntime))

# and one with 30 years at 1 then 10 years at 0.8
# scalar_array <- array(data = c(rep(matrix(1, nrow = nlayer, ncol = nbox), 10),
#                                rep(matrix(0.8, nrow = nlayer, ncol = nbox), 5)),
#                       dim = c(nlayer, nbox, ntime))

# for heatwave runs
# ZM
scalar_array <- array(data = c(rep(matrix(1, nrow = nlayer, ncol = nbox), 62),
                               rep(matrix(0.3, nrow = nlayer, ncol = nbox), 4),
                               rep(matrix(1, nrow = nlayer, ncol = nbox), 4)),
                      dim = c(nlayer, nbox, ntime))

# and PL
# scalar_array_PL <- array(data = c(rep(matrix(1, nrow = nlayer, ncol = nbox), 31),
#                                   rep(matrix(0.5, nrow = nlayer, ncol = nbox), 20)),
#                                   #rep(matrix(1, nrow = nlayer, ncol = nbox), 5)),
#                          dim = c(nlayer, nbox, ntime))

dim.def.nc(nc_file, "t", unlim=TRUE)
dim.def.nc(nc_file, "b", nbox) # manual
dim.def.nc(nc_file, "z", nlayer) # manual I am not sure about this, do we need the sediment layer in the fluxes??

var.def.nc(nc_file, "t", "NC_DOUBLE", "t")
#var.def.nc(nc_file, "PSgrowth_ff", "NC_DOUBLE", c("z","b","t"))
var.def.nc(nc_file, "PLgrowth_ff", "NC_DOUBLE", c("z","b","t"))
#var.def.nc(nc_file, "ZSgrowth_ff", "NC_DOUBLE", c("z","b","t"))
var.def.nc(nc_file, "ZMgrowth_ff", "NC_DOUBLE", c("z","b","t"))
#var.def.nc(nc_file, "ZLgrowth_ff", "NC_DOUBLE", c("z","b","t"))
var.def.nc(nc_file, "EUPgrowth_ff", "NC_DOUBLE", c("z","b","t"))

#att.put.nc(nc_file, "PSgrowth_ff", "_FillValue", "NC_DOUBLE", 1)
att.put.nc(nc_file, "PLgrowth_ff", "_FillValue", "NC_DOUBLE", 1)
#att.put.nc(nc_file, "ZSgrowth_ff", "_FillValue", "NC_DOUBLE", 1)
att.put.nc(nc_file, "ZMgrowth_ff", "_FillValue", "NC_DOUBLE", 1)
#att.put.nc(nc_file, "ZLgrowth_ff", "_FillValue", "NC_DOUBLE", 1)
att.put.nc(nc_file, "EUPgrowth_ff", "_FillValue", "NC_DOUBLE", 1)
att.put.nc(nc_file, "t", "units", "NC_CHAR", t_units)
att.put.nc(nc_file, "t", "dt", "NC_DOUBLE", seconds_timestep)
att.put.nc(nc_file, "NC_GLOBAL", "title", "NC_CHAR", "External scalar for plankton growth")
att.put.nc(nc_file, "NC_GLOBAL", "geometry", "NC_CHAR", this_geometry)
att.put.nc(nc_file, "NC_GLOBAL", "parameters", "NC_CHAR", " ")

var.put.nc(nc_file, "t", time_array)
#var.put.nc(nc_file, "PSgrowth_ff", scalar_array)
var.put.nc(nc_file, "PLgrowth_ff", scalar_array)
#var.put.nc(nc_file, "ZSgrowth_ff", scalar_array)
var.put.nc(nc_file, "ZMgrowth_ff", scalar_array)
#var.put.nc(nc_file, "ZLgrowth_ff", scalar_array)
var.put.nc(nc_file, "EUPgrowth_ff", scalar_array)


close.nc(nc_file)
