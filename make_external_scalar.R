# Alberto Rovellini
# 3/7/2023
# Code to write external scaling input for plankton growth rates

# See https://confluence.csiro.au/display/Atlantis/External+Mortality%2C+Growth+and+Recruitment+Scaling

library(tidync)
library(ncdf4)
library(RNetCDF)

getwd()

nc_name <- '../data/scalar_test_08.nc'

nc_file <- create.nc(nc_name)

nbox <- 109
nlayer <- 7
t_units <- "seconds since 1990-01-01 00:00:00 -9" # for consistency with the init.nc files
seconds_timestep <- 43200 * 2 * 365
this_geometry <- "GOA_WGS84_V4_final.bgm"

# time array
time_array <- seq(0, seconds_timestep * 14, seconds_timestep)
ntime <- length(time_array)

# scalar array
# dummy array with 1 in all positions, this should result in no changes in the run
# scalar_array <- array(data = 1, dim = c(nlayer, nbox, ntime))

# and one with 30 years at 1 then 10 years at 0.8
scalar_array <- array(data = c(rep(matrix(1, nrow = nlayer, ncol = nbox), 10),
                               rep(matrix(0.8, nrow = nlayer, ncol = nbox), 5)),
                      dim = c(nlayer, nbox, ntime))

dim.def.nc(nc_file, "t", unlim=TRUE)
dim.def.nc(nc_file, "b", nbox) # manual 
dim.def.nc(nc_file, "z", nlayer) # manual I am not sure about this, do we need the sediment layer in the fluxes??

var.def.nc(nc_file, "t", "NC_DOUBLE", "t")
var.def.nc(nc_file, "ZMgrowth_ff", "NC_DOUBLE", c("z","b","t"))

att.put.nc(nc_file, "ZMgrowth_ff", "_FillValue", "NC_DOUBLE", 1)
att.put.nc(nc_file, "t", "units", "NC_CHAR", t_units)
att.put.nc(nc_file, "t", "dt", "NC_DOUBLE", seconds_timestep)
att.put.nc(nc_file, "NC_GLOBAL", "title", "NC_CHAR", "Expernal scalar for zooplankton growth")
att.put.nc(nc_file, "NC_GLOBAL", "geometry", "NC_CHAR", this_geometry)
att.put.nc(nc_file, "NC_GLOBAL", "parameters", "NC_CHAR", " ")

var.put.nc(nc_file, "t", time_array)
var.put.nc(nc_file, "ZMgrowth_ff", scalar_array)

close.nc(nc_file)