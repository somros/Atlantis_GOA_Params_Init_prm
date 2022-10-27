# view the init.nc that we produced

library(shinyrAtlantis)
library(tidyverse)

bgm.file <-'../data/GOA_WGS84_V4_final.bgm'
nc.file <- '../data/GOA_cb_new.nc'

obj <- make.sh.init.object(bgm.file, nc.file)
sh.init(obj)
