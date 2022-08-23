#' Build NC initial conditions files for Atlantis GOA
#' 
#' @author Alberto Rovellini
#' @date February 2022
#' 

# Setup -------------------------------------------------------------------

# List of packages for session
.packages = c("Rcpp","shinyrAtlantis","devtools","tidyverse","stringi","ncdf4", "data.table", "rbgm", "sf")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst], dependencies = TRUE)

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

#devtools::install_github("shanearichards/shinyrAtlantis")
#Forked by J. Porobic 
#devtools::install_github("Atlantis-Ecosystem-Model/shinyrAtlantis")

#NOTE: shinyRAtlantis as of 2/4/2022 does not include grouptypes JELLIES and SPONGES. If you have these groups, you need to add them to the AttributeTemplate.csv file in the package data

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

#' clean up the space
rm(list=ls())
select <- dplyr::select

data_path <- "../data"
results_path <- "../output"

setwd(data_path)

# Vertebrates -------------------------------------------------------------

grp.file   <- "GOA_Groups.csv"
bgm.file   <- "GOA_WGS84_V4_final.bgm"
cum.depths <- c(0,30,100,200,500,1000,4000)
csv.name   <- "GOAtemplate"

func.groups <- read_csv(grp.file) %>% 
  select(Code,Name,LongName) 

# Used initially
#creates empty csv files for editing to then produce NetCDF files
csv.name   <- "GOAtemplate"
make.init.csv(grp.file, bgm.file, cum.depths, csv.name)

# Fixing a couple of issues with horiz.csv and init.csv that create problems in init.nc
# 1. The NetDCF is written out without long-name attribute for Rugosity. This tracer seems to be the only one where name=longname in init.csv. Fix that.

init.data <- read.csv('GOAtemplate_init.csv')
init.data[init.data$name=='Rugosity',]$long_name <- 'Rugosity (R)'

# 2. the tracer pH does not feature in init.csv. Let's add it - copying atts from CalCurrent
pH.frame <- data.frame("pH", TRUE, 2, "[ z b ]", "tracer", "log10H", "pH of water", 0, 1, 1, 1, 1, 0, 0, 8.1, NA, NA, NA, NA, NA, NA, NA, "constant", "uniform", 8.1, 8.1)
colnames(pH.frame)<-colnames(init.data)

init.data1 <- rbind(init.data[1:which(init.data$name=='Temp'),],
                   pH.frame,
                   init.data[(which(init.data$name=='Temp')+1):nrow(init.data),])

write.csv(init.data1,'GOAtemplate_init.csv',row.names = F)
rm(init.data, init.data1)

#takes life history created in GOA_Bioparam_table.gsheet and estimates numbers by age

life.history <- read_csv("life_history_parameters.csv")
mammals.w.age <-  read_csv("marine_mammal_w_age.csv") # for marine mammals weight at age are known
seasonal.distribution <- read_csv("seasonal_distribution.csv") 
# bird mammal consumption is used to estimate mum and calculated in mammal_x.xls
bird.mammal.consumption <- read_csv("bird_mammal_consumption.csv")

migratory.proportion <- read_csv("migiobox_proportion.csv")
#these are groups that have all or a proportion of the population outside the model on January 1st
#see Table Migration parameters in Tech memo

mig.groups <- c("WHT","WHG","WHH","PIN","BDF","BDI","BSF","SCH","SCM","SCO","SPI","SSO","HAK") # "SHP" and "EUL" migrate at other times of the year, i.e. they are in the model on Jan 1

salmon.codes <- NA # this was used in AMPS

#make.init.nc(bgm.file, cum.depths, init.file, horiz.file, nc.file)
init.file  <- "GOAtemplate_init.csv"
horiz.file <- "GOAtemplate_horiz.csv"
nc.file    <- "GOAtemplate.nc"

horiz.data <- read_csv(horiz.file)
init.data <- read_csv(init.file)

#clearance_parameters  (Cmax = AC * W^BC)
ac_clearance  <-  0.3
bc_clearance  <-  -0.3 +1
number_1 <- 0.8 # need to track this down
number_2 <-3.3 # see manual from Fish Bioenergetics 3.0
#Horne et al. 2010 Tech Memo, near equation 13: 
#"We considered the consumption rates from this equation to represent daily averages, and we assumed that individuals generally operate at about 30% 
#of their potential maximum. Thus we multiplied the resulting Gmax by three to obtain theoretical maxima. As with fish groups, we assumed a growth 
#efficiency of 10% for calculating the maximum growth rate “g.” 

age.structured.sp <- life.history %>% distinct(Code) %>% .$Code

mammals.sp <- c("KWT","WHT","KWR","DOL","WHG","WHH","WHB","SSL","PIN")
birds.sp <- c("BDF","BDI","BSF","BSI")

get_age_numbers <- function(eachspcode){
  
  print(eachspcode)
  
  this.life.history <- life.history %>% 
    filter(Code==eachspcode)
  
  this.longname <- func.groups %>% filter(Code==eachspcode) %>% .$LongName
  this.name <- func.groups %>% filter(Code==eachspcode) %>% .$Name
  
   
  if(eachspcode %in% mig.groups){
    
    #mig_prop is the proportion that stays in the model
    mig.prop <- migratory.proportion %>% 
      filter(atlantis_fg == eachspcode) %>% 
      mutate(mig_prop= 1-migiobox) %>% 
      .$mig_prop
    
    out.prop = 1- mig.prop
    
  } else {
    
    out.prop = 0
    mig.prop = 1
    
  }
  
  if(eachspcode %in% mammals.sp){
    
    wet.weight <- mammals.w.age %>%
      dplyr::select(contains(eachspcode)) %>% 
      drop_na() %>%
      setNames(c("wet_weight_kg","actual_age")) %>% 
      mutate(Code = eachspcode,
             age_class=1:nrow(.),
             wet_weight_g = wet_weight_kg * 1000) %>% 
      dplyr::select(Code,age_class,actual_age,wet_weight_g)
    
    wet.weight.recruitment <- tibble(Code = eachspcode,
                                     age_class = 0,
                                     actual_age = 0,
                                     wet_weight_g = (this.life.history$a_FUNC*((this.life.history$Linf_FUNC*(1-exp(-this.life.history$k_FUNC*(min(wet.weight$actual_age)/2)))))^this.life.history$b_FUNC))
    
    
    nums_age <- bind_rows(wet.weight.recruitment,wet.weight)
    
    this.consumption <- bird.mammal.consumption %>% 
      filter(atlantis_fg==eachspcode) %>% 
      .$`Individual Consumption (Pred Mass^-1 * day^-1)`
    
  } else if(eachspcode %in% birds.sp) {
    
    nums_age <- tibble(Code = eachspcode,
                       age_class=0:this.life.history$num_age_class,
                       actual_age=this.life.history$ypa_FUNC*age_class,
                       wet_weight_g = this.life.history$Average_body_size_g)
    
    this.consumption <- bird.mammal.consumption %>% 
      filter(atlantis_fg==eachspcode) %>% 
      .$`Individual Consumption (Pred Mass^-1 * day^-1)`
    
    
  } else {
    
    wet.weight <- tibble(Code = eachspcode,
                       age_class=1:this.life.history$num_age_class, 
                       actual_age=this.life.history$ypa_FUNC*age_class,
                       wet_weight_g = this.life.history$a_FUNC*((this.life.history$Linf_FUNC*(1-exp(-this.life.history$k_FUNC*(actual_age)))))^this.life.history$b_FUNC)
    
    wet.weight.recruitment <- tibble(Code = eachspcode,
                                     age_class = 0,
                                     actual_age = 0,
                                     wet_weight_g = this.life.history$a_FUNC*((this.life.history$Linf_FUNC*(1-exp(-this.life.history$k_FUNC*(min(wet.weight$actual_age)/2)))))^this.life.history$b_FUNC)
    
                                     
          #wet_weight_g = (this.life.history$a_FUNC*((this.life.history$Linf_FUNC*(1-exp(-this.life.history$k_FUNC*(this.life.history$arecruit_FUNC/360)))))^this.life.history$b_FUNC))
    
     nums_age <- bind_rows(wet.weight.recruitment,wet.weight)
    
  }
  
  if(eachspcode %in% mammals.sp) {
   # for marine mammals w at age is known, found in marine_mammals_w.xls
    nums_age_table <- nums_age %>% 
      mutate(struct_N_mg = wet_weight_g*this.life.history$Ntoall/this.life.history$drytowet*this.life.history$sNtN*1000,
             resN_mg = wet_weight_g*this.life.history$Ntoall/this.life.history$drytowet*this.life.history$rNtN*1000,
             mum_growth_rate = (struct_N_mg + resN_mg)*this.consumption*number_2,
             clearance = mum_growth_rate/10,
             exponential_decay = exp(-this.life.history$M_FUNC*(age_class-1)*this.life.history$ypa_FUNC),
             proportional_age_distribution = exponential_decay/ sum(exponential_decay),
             avg_ind_w = wet_weight_g * proportional_age_distribution,
             total_numbers = this.life.history$Biomass_mt*1e+06/sum(avg_ind_w),
             numbers_at_age = round(proportional_age_distribution*total_numbers,0),
             #prop_age_j = if_else(age_class<this.life.history$mat_FUNC,1,0), proportional_age_dist_juv = prop_age_j/sum(prop_age_j), # potential issue here
             #prop_age_a = if_else(age_class>this.life.history$mat_FUNC,1,0), proportional_age_dist_a = prop_age_a/sum(prop_age_a),
             tot_biomass_age = numbers_at_age * wet_weight_g,
             prop_biomass_age = tot_biomass_age / sum(tot_biomass_age),
             numbers_mig = numbers_at_age * mig.prop, 
             nums_out_mig = numbers_at_age * out.prop)
    
  } else if (eachspcode %in% birds.sp){
    
    #birds have one weigth,they fledge at full size
    nums_age_table <- nums_age %>% 
      mutate(struct_N_mg = wet_weight_g*this.life.history$Ntoall/this.life.history$drytowet*this.life.history$sNtN*1000,
             resN_mg = wet_weight_g*this.life.history$Ntoall/this.life.history$drytowet*this.life.history$rNtN*1000,
             mum_growth_rate = (struct_N_mg + resN_mg)*this.consumption*number_2,
             clearance = mum_growth_rate/10,
             exponential_decay = exp(-this.life.history$M_FUNC*(age_class-1)*this.life.history$ypa_FUNC),
             proportional_age_distribution = exponential_decay/ sum(exponential_decay),
             avg_ind_w = this.life.history$Average_body_size_g,
             total_numbers = this.life.history$Biomass_mt*1e+06/avg_ind_w,
             numbers_at_age = round(proportional_age_distribution*total_numbers,0),
             #prop_age_j = if_else(age_class<this.life.history$mat_FUNC,1,0), proportional_age_dist_juv = prop_age_j/sum(prop_age_j), # potential issue here
             #prop_age_a = if_else(age_class>this.life.history$mat_FUNC,1,0), proportional_age_dist_a = prop_age_a/sum(prop_age_a),
             tot_biomass_age = numbers_at_age * wet_weight_g,
             prop_biomass_age = tot_biomass_age / sum(tot_biomass_age),
             numbers_mig = numbers_at_age * mig.prop,
             nums_out_mig = numbers_at_age * out.prop)
  } else {
    
    # for everything else average weight is calculated
    nums_age_table <- nums_age %>% 
      mutate(struct_N_mg = wet_weight_g*this.life.history$Ntoall/this.life.history$drytowet*this.life.history$sNtN*1000,
             resN_mg = wet_weight_g*this.life.history$Ntoall/this.life.history$drytowet*this.life.history$rNtN*1000,
             mum_growth_rate = ac_clearance *((struct_N_mg + resN_mg)^bc_clearance)*number_1 * number_2,
             clearance = mum_growth_rate/10,
             exponential_decay = exp(-this.life.history$M_FUNC*(age_class-1)*this.life.history$ypa_FUNC),
             proportional_age_distribution = exponential_decay/ sum(exponential_decay),
             avg_ind_w = wet_weight_g * proportional_age_distribution,
             total_numbers = ceiling(this.life.history$Biomass_mt*1e+06/sum(avg_ind_w)),
             numbers_at_age = ceiling(proportional_age_distribution*total_numbers),
             #prop_age_j = if_else(age_class<this.life.history$mat_FUNC,1,0), proportional_age_dist_juv = prop_age_j/sum(prop_age_j), # potential issue here
             #prop_age_a = if_else(age_class>=this.life.history$mat_FUNC,1,0), proportional_age_dist_a = prop_age_a/sum(prop_age_a),
             tot_biomass_age = numbers_at_age * wet_weight_g,
             prop_biomass_age = tot_biomass_age / sum(tot_biomass_age),
             numbers_mig = ceiling(numbers_at_age * mig.prop),
             nums_out_mig = ceiling(numbers_at_age * out.prop))
  }
  
  print(nums_age_table)
  
}

num.biomass.age <- lapply(age.structured.sp, get_age_numbers) %>% 
  bind_rows %>% 
  left_join(func.groups, by="Code")

num.biomass.age %>% # skipping age 0 for now it seems - check that we are not double-counting the age 0 animals since our biomass estimates already include them!!
  filter(age_class!=0) %>% 
  group_by(Code) %>% 
  summarise(tot_nums=sum(numbers_at_age)) %>% 
  write_csv("tot_nums_groups.csv")

write_csv(num.biomass.age, "nums_age_functional_groups.csv")

#make nums, strN, and resN for each group, to be added to horiz.data
nc_vector <- function(eachspcode, num.biomass.age) {
  
  print(eachspcode)
  
  #get life history for this group
  this.life.history <- life.history %>% 
    filter(Code==eachspcode)
  #atlantis name and long name
  this.longname <- func.groups %>% filter(Code==eachspcode) %>% .$LongName
  this.name <- func.groups %>% filter(Code==eachspcode) %>% .$Name
  
  print(this.longname)
  #weight at age, nums, resN and strN
  nums.age.table <- num.biomass.age %>% 
    filter(Code == eachspcode)
  
  print(this.life.history$Distribution)
  
  #get salmon distribution, uses adult distributions # ALBI: this does not happen for the GOA yet, salmon groups set to NA
  if(eachspcode %in% salmon.codes){
    
    this.distribution <- salmon.distribution %>% 
      dplyr::select(contains(this.name)) %>% 
      dplyr::select(ends_with("_S1")) %>% 
      setNames("distribution")
    
  } else {
    #if not salmon, distribution can be seasonal or monthly
    if(this.life.history$Distribution == "Seasonal") {
      
      this.distribution.j <- seasonal.distribution %>% 
        dplyr::select(contains(this.name)) %>% 
        dplyr::select(ends_with("J_S1")) %>% # how do we handle juvenile vs adult distributions?
        setNames("distribution")
      
      this.distribution.a <- seasonal.distribution %>% 
        dplyr::select(contains(this.name)) %>% 
        dplyr::select(ends_with("A_S1")) %>% # how do we handle juvenile vs adult distributions?
        setNames("distribution")
      
    }
    
  }
  
  #nums at age as numeric
  nums.age <- nums.age.table %>% 
    filter(age_class!=0) %>% 
    dplyr::select(numbers_at_age) %>% 
    .$numbers_at_age
  
  boxes <- paste("box",0:108,sep="") # could set this soft to box no from bgm
  
  num.age.names <- 1:length(nums.age)
  
  print(paste("num age classes is",nums.age))
  
  #check if species code starts at totally or partially outside of the model domain on January 1st 
  #multiply nums,strN, and resN by residency if so
  if(eachspcode %in% mig.groups){
    
    mig.prop <- migratory.proportion %>% 
      filter(atlantis_fg == eachspcode) %>% 
      mutate(mig_prop= 1-migiobox) %>% 
      .$mig_prop
    
    out.prop = 1- mig.prop
  }
  
  #make vector for each age class
  calc_age <- function(eachage) {
    
    print(paste("Analyzing age",eachage))
    #get value
    
    # set this.distribution to Juvenile or Adult as appropriate based on the age at maturity of the species and the age of the group here
    if(eachage < life.history[life.history$Code==eachspcode,]$mat_FUNC){
      this.distribution <- this.distribution.j
    } else {
      this.distribution <- this.distribution.a
    }
    
    this.num.age  <- nums.age.table %>% 
      filter(age_class==eachage) %>% 
     .$numbers_at_age
    
   #get name
    this.age.name <- paste(paste(this.name,eachage,sep=""),"Nums",sep="_")
    
 this.strN.age <- nums.age.table %>% 
      filter(age_class==eachage) %>% 
      .$struct_N_mg
 
 this.strN.name <- paste(paste(this.name,eachage,sep=""),"StructN",sep="_")
 
 
 this.resN.age <- nums.age.table %>% 
   filter(age_class==eachage) %>% 
   .$resN_mg
 
 this.resN.name <- paste(paste(this.name,eachage,sep=""),"ResN",sep="_")
 
    
    if(exists("mig.prop")){
      this.num.age <-  this.num.age * mig.prop
      
        }
    
 init.data <<- init.data %>% 
   mutate(fill.value=replace(fill.value, name==this.strN.name, this.strN.age)) %>% 
   mutate(fill.value=replace(fill.value, name==this.resN.name, this.resN.age)) 
 
    test.dist <- sum(this.distribution$distribution)
    
    if(test.dist == 0 | this.num.age == 0) {
      
      age.distribution <- tibble(value = 0,
                                 box = boxes) 
      
      this.age.distribution <- age.distribution %>% pivot_wider(names_from=box, values_from=value) %>% 
        mutate(variable = this.age.name)
      
   
      this.sp.distribution <- this.age.distribution %>% 
        dplyr::select(variable,everything())
      
    } else {
      
      norm.distribution <- this.distribution %>% 
        mutate(distribution = if_else(is.na(distribution),0,distribution), 
               distribution_norm = distribution/ sum(distribution))
      
      this.age.distribution <- norm.distribution %>% 
        mutate(value = this.num.age * distribution_norm,
               variable = this.age.name,
               box = boxes) %>% 
        dplyr::select(variable, value, box)
      
      this.sp.distribution <- this.age.distribution%>% # this step rearranges the order of the boxes- does it matter?
        pivot_wider(names_from=box,values_from=value)
      
    }
    
    print(this.sp.distribution)
    
    this.rnsn.distribution <- matrix(0,2,109) %>% # we are setting these all to 0 and then using the fillvalue as ResN and StrucN are constant throught 
      as_tibble() %>% 
      setNames(c(paste("box",0:108,sep=""))) %>% 
      mutate(variable = c(this.strN.name, this.resN.name)) %>% 
      select(variable,box0:box108) %>%
      bind_rows(this.sp.distribution)
      
     
    return(this.rnsn.distribution)
    
  }
  
  age.list <- 1:length(num.age.names)
  age.nums <- lapply(age.list,calc_age) %>% 
    bind_rows
  
  col.order <- c("variable",boxes)
  
  age.vector  <- age.nums[c(col.order)] 
  
  return(age.vector)
  
}

#this final vector should have 1653 rows

num.biomass.frame <- lapply(age.structured.sp, nc_vector, num.biomass.age = num.biomass.age) %>% 
  bind_rows %>% 
  arrange(variable) %>% 
  dplyr::rename(Variable=variable)

#this can be used if it is the first time you create init files
init.file  <- "GOA_init.csv"
horiz.file <- "GOA_horiz.csv"
nc.file    <- "GOA_cb_new.nc"

#add Nums to empty horiz.csv file created above
horiz.data %>% 
  filter(!grepl("_Nums",Variable)) %>% 
  mutate_if(is.numeric, round, 5) %>% 
  bind_rows(num.biomass.frame) %>%
  write_csv(horiz.file)

init.data %>% 
  write_csv(init.file)

# AR: here using a modified version of make.init.nc() to solve the issue of zeroes being packed as cdf vectors instead of '_'
source(file = '../code/modified_funs.R')

make.init.nc.ar(bgm.file,cum.depths,init.file,horiz.file,nc.file) # check if values for RN and SN are 0 or _

# Biomass pools -----------------------------------------------------------

# will need the BGM for box areas here
atlantis.bgm <- read_bgm('../data/GOA_WGS84_V4_final.bgm')
atlantis.box <- atlantis.bgm %>% box_sf()

# Plankton and nutrients:
# values calculated in init_calculator.R at https://github.com/somros/plankton_nutrients_init
# local address is C:/Users/Alberto Rovellini/Documents/GOA/SDM/Plankton_and_nutrients/outputs/init/

plankton.files <- list.files('C:/Users/Alberto Rovellini/Documents/GOA/SDM/Plankton_and_nutrients/outputs/init/', full.names = T)

outnc <- nc_open(nc.file, write=TRUE) # open .nc file

for(i in 1:length(plankton.files)){
  this.file <- plankton.files[i]
  this.var <- gsub('.txt','',gsub('.*init/','',this.file))
  this.data <- read.table(this.file,sep=',')
  this.data.horiz <- t(this.data[,-ncol(this.data)])
  
  ncvar_put(outnc, varid = this.var, vals = this.data.horiz)
}

nc_close(outnc)

# Invertebrates that we have SDMs for
# We read in virgin biomass and spatial distributions
# read in spatial distributions of invertebrates
invert.distrib <- read.csv('../data/seasonal_distribution_inverts.csv')

# read in virgin biomass
# As of 2/8/2022 this is mainly from biomass estimates per km2 from Aydin et al. (2007) multiplied by the model area, so to be treated as ballpark estimates
# Our own SDMs are not reliable enough for benthic invertebrates to come up with our own biomass estimates
invert.biomass <- read.csv('../data/inverts_virgin_biomass.csv')

# select species for which we have info to do this
dist.species <- unique(gsub('_A_S.','',colnames(invert.distrib)))
biom.species <- invert.biomass$name
all.species <- intersect(dist.species,biom.species) # these are the species we will do this for

# let's drop the plankton from here, as we already get raw concentration from NPZ
all.species <- setdiff(all.species, c('Euphausiids','Macrozooplankton','Mesozooplankton','Microzooplankton'))

# make and write concentrations per box - in mg N m-2 for benthos and in mg N m-3 for CEP and PWN, for these assuming that at t0
# all the biomass in the box is concentrated in the first 1 m of depth from the bottom - it should get redistributed by the
# vertical params in biol.prm at t0 anyway
# do we need to put infauna in the sediment? 

outnc <- nc_open(nc.file, write=TRUE) # open .nc file

# we need to get the thickness of the bottom layer of each box, to divide the concentration per m2

make.lyr <- function(depth){
  lyrvec <- c(cum.depths[cum.depths<depth],depth)
  lyrframe <- data.frame(lyr = length(lyrvec):1, lower = lyrvec)
  return(lyrframe)
}

depth.key <- atlantis.box %>%
  st_set_geometry(NULL) %>%
  mutate(botz = -botz) %>%
  mutate(lyr = purrr::map(botz, make.lyr)) %>%
  unnest(cols = c(lyr)) %>%
  group_by(box_id) %>%
  mutate(dz = lower-lag(lower, default = 0)) %>%
  ungroup() %>%
  filter(lyr==1) %>%
  select(box_id,dz)

for(i in 1:length(all.species)){
  
  this.species <- all.species[i]
  this.dist <- invert.distrib %>% select(contains(paste(this.species,'A_S1',sep='_')))
  this.unit <- init.data[grepl(this.species,init.data$name),]$units
  
  this.init <- this.dist %>%
    set_names('Prop') %>%
    mutate(Wet.weight.t = Prop*invert.biomass[invert.biomass$name==this.species,]$biomass_mt,
           mgN = Wet.weight.t/20/5.7*1000000000, #t to mg 10^9
           .bx0 = 0:108) %>%
    left_join((atlantis.box %>% st_set_geometry(NULL) %>% select(.bx0,area)), by = '.bx0') %>%
    mutate(mgNm2 = round(mgN/area,digits=5)) %>%
    select(mgNm2) %>%
    t()
  
  # have them all on the bottom, even CEP and PWN that are mg N m-3. Need to look into how these guys get distributed with vert
  # If BC, BD, BO, put them in the sediment instead
  
  infauna <- c('Benthic_carnivores','Meiobenthos','Deposit_feeders')
  
  if(this.unit[1]=='mg N m-3') {
    if(this.species %in% infauna) {
      this.init <- rbind(matrix(0,nrow=6,ncol=ncol(this.init)), # sediment layer is 1 m so this is fine
                         this.init)
    } else { # attach depth key and divide by layer thickness so that we have conc per m3 
      
      this.init <- depth.key %>%
        mutate(Nm2 = this.init[1,],
               Nm3 = Nm2/dz) %>%
        mutate(Nm3 = replace(Nm3, is.nan(Nm3), 0)) %>%
        pull(Nm3) %>%
        rbind(matrix(0,nrow=6,ncol=ncol(this.init)))
    }
  }
  
  ncvar_put(outnc, varid = paste(this.species,'N',sep='_'), vals = this.init)

}

# note that we are getting some pretty high values here for some of this benthos. 
# this is entirely reliant on old ecopath estimates, so it may not be sensible

nc_close(outnc)

# Habitat types -----------------------------------------------------------

# uses the habitat cover values calculated at Bottom_cover, currently hosted at https://github.com/somros/GOA_bottom_cover

cover <- read.csv('../../../Bottom_cover/abiotic_habitat.csv')

reef <- cover %>% filter(atlantis_class=='Reef') %>% mutate(cover = round(cover, digits = 3)) %>% select(cover) %>% t()
flat <- cover %>% filter(atlantis_class=='Sand') %>% mutate(cover = round(cover, digits = 3)) %>% select(cover) %>% t()
soft <- cover %>% filter(atlantis_class=='Soft') %>% mutate(cover = round(cover, digits = 3)) %>% select(cover) %>% t()

outnc <- nc_open(nc.file, write=TRUE) # open .nc file

ncvar_put(outnc, varid = 'reef', vals = reef)
ncvar_put(outnc, varid = 'flat', vals = flat)
ncvar_put(outnc, varid = 'soft', vals = soft)

nc_close(outnc)

# Tracers ----------------------------------------------------------

# Some important tracers are initialized with default values in the wc and sed in the shinyRAtlantis package.
# However, according to the log file they may limit the initial distributions of several groups. 
# These tracers are: Temp, salt, Oxygen, pH, Si, Det_Si
# possibly more (e.g. Chl-a and DON) - keep an eye on these
# I am not sure why the code won't just let the fillvalues do their part, but instead of changing the fill values and the function, I am fixing them here

# Let's set them to values that make more sense:
# Temp = 6 C (as the fillvalue of 15 will kill off some of my critters with the thermal windows we got from Aquamaps)
# salt = 33 psu
# Oxygen = 8000 mg O2 m-3
# pH = 8.1
# Si = 2194 mg Si m-3 https://www.sciencedirect.com/science/article/abs/pii/030442039190021N
# Alternatively Si = 336 mg Si m-3 (from 12 micromolar for Silicates in Jan 1990 in <100 m waters on GOA shelf from COBALT https://gulf-of-alaska.portal.aoos.org/#map)
# Det_Si = 2194 mg Si m-3 https://www.sciencedirect.com/science/article/abs/pii/030442039190021N

outnc <- nc_open(nc.file, write=TRUE) # open .nc file

tracers <- data.frame('tracer'=c('Temp','salt','Oxygen','pH','Si','Det_Si'), 'fillval'=c(6,33,8000,8.1,336,336))

for (i in 1:nrow(tracers)){
  this_tracer <- tracers[i,]$tracer
  this_table <- ncvar_get(outnc, this_tracer)
  this_table[which(this_table != 0)] <- tracers[i,]$fillval
  this_table[which(this_table == 0)] <- NA # and set missing layers to NA so that they get packed as '_', it will not matter but good for consistency
  ncvar_put(outnc, varid = this_tracer, vals = this_table)
}

nc_close(outnc)

# DayLight ----------------------------------------------------------------
# As of code version v6645 and later we need to add the Tracer DayLight. Initialize it as the same as Light
require(RNetCDF)

outnc <- open.nc(nc.file, write=TRUE) # open .nc file

light <- var.get.nc(outnc, "Light")

# define the new variable
var.def.nc(outnc, "DayLight", "NC_DOUBLE", c("z","b","t"))

# prepare all attributes for the new variable
att.put.nc(outnc, "DayLight", "units", "NC_CHAR", "missing")
att.put.nc(outnc, "DayLight", "_FillValue", "NC_DOUBLE", 0)
att.put.nc(outnc, "DayLight", "long_name", "NC_CHAR", "Daylight intensity on the surface of a cell")
att.put.nc(outnc, "DayLight", "bmtype", "NC_CHAR", "tracer")
att.put.nc(outnc, "DayLight", "dtype", "NC_INT", 0)
att.put.nc(outnc, "DayLight", "sumtype", "NC_INT", 1)
att.put.nc(outnc, "DayLight", "inwc", "NC_INT", 0)
att.put.nc(outnc, "DayLight", "insed", "NC_INT", 0)
att.put.nc(outnc, "DayLight", "dissol", "NC_INT", 1)
att.put.nc(outnc, "DayLight", "decay", "NC_DOUBLE", 0)
att.put.nc(outnc, "DayLight", "partic", "NC_INT", 0)
att.put.nc(outnc, "DayLight", "fill.value", "NC_DOUBLE", 0)

# add variable to the NC file
var.put.nc(ncfile = outnc, variable = "DayLight", data = light, count = c(NA, NA, NA))

# close the file
close.nc(outnc)

# Bacteria ----------------------------------------------------------------

# BB and PB are all packed as 0 by the make code. This means that their fillvalues do not stick (and who knows about those fillvalues anyway)
# Here I coerce them to 0.1 for PB in all layers but sed, and to 0.1 in sediment for BB

outnc <- nc_open(nc.file, write=TRUE) # open .nc file

bacts <- c('Pelagic_bacteria_N','Benthic_bacteria_N')

for(i in 1:length(bacts)){
  this_bact <- bacts[i]
  this_table <- ncvar_get(outnc, this_bact)
  if(this_bact=='Pelagic_bacteria_N'){
    this_table1 <- rbind(matrix(0.1,nrow = (nrow(this_table)-1), ncol = ncol(this_table)),rep(0,109))
  } else {
    this_table1 <- rbind(matrix(0,nrow = (nrow(this_table)-1), ncol = ncol(this_table)),rep(0.1,109))
  }
  ncvar_put(outnc, varid = this_bact, vals = this_table1)
}

nc_close(outnc)
