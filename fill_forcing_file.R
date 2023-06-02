
getwd()

file.create('fill_force.txt')

y <- 70
spinup <- 40 + 5 # accounting for 1991-1995

# hydro
for (i in 0:(y-1)){
  if(i < (spinup)) {
    cat(paste0('hd',i, '.name ../hydro_forcings/hydro/goa_hydro_1999.nc', '\n'), file = 'fill_force.txt', append = T)
  } else {
    cat(paste0('hd',i ,'.name ../hydro_forcings/hydro/goa_hydro_', (1996 + i - spinup), '.nc', '\n'), 
        file = 'fill_force.txt', append = T)
  }
}

# temp
for (i in 0:(y-1)){
  if(i < (spinup)) {
    cat(paste0('Temperature',i, '.name ../hydro_forcings/temp/goa_temp_1999.nc', '\n'), file = 'fill_force.txt', append = T)
  } else {
    cat(paste0('Temperature',i ,'.name ../hydro_forcings/temp/goa_temp_', (1996 + i - spinup), '.nc', '\n'), 
        file = 'fill_force.txt', append = T)
  }
}

# salt
for (i in 0:(y-1)){
  if(i < (spinup)) {
    cat(paste0('Salinity',i, '.name ../hydro_forcings/salt/goa_salt_1999.nc', '\n'), file = 'fill_force.txt', append = T)
  } else {
    cat(paste0('Salinity',i ,'.name ../hydro_forcings/salt/goa_salt_', (1996 + i - spinup), '.nc', '\n'), 
        file = 'fill_force.txt', append = T)
  }
}
