#' Code to set up biology prm for Atlantis GOA model
#' 
#' Reads in values and writes prm
#' @author Alberto Rovellini, adapted from Hem Nalini Morzaria Luna's code for AMPS
#' @keywords Gulf of Alaska, Atlantis
#' 
#' January 2022
#' 
#' 
rm(list=ls())

#devtools::install_github("tidyverse/googlesheets4")

.packages = c("cowplot","RNetCDF", "data.table", "tidyverse","magrittr","ggplot2",
              "readxl", "httr","gsubfn","parallel")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

#functional groups
grp.file <- read_csv("../data/GOA_Groups.csv")
#last three slots are assumed to be dedtritus groups

# nums and w at age, estimated in PS_initial_conditions.nc
num.biomass.age <- read_csv("../data/Nums_age_functional_groups.csv") # UPDATE AFTER NOTES AND BC
life.history <- read_csv("../data/life_history_parameters.csv") # UPDATE AFTER DOING NOTES
#distributions
seasonal.distribution <- read_csv("../data/seasonal_distribution.csv") # UPDATE AFTER CANADA REWORK AND ISAAC"S SALMON
invert.distribution <- read_csv("../data/seasonal_distribution_inverts.csv") # UPDATE AFTER CANADA REWORK

prey.matrix <- read_csv("../data/goa_pprey_matrix.csv")

#nc file set up has spatial distributions inverts
goa.horiz <- read_csv("../data/GOAtemplate_horiz.csv")

#distributions for these invert groups are read separately
invert.groups <- c("Octopus", "Squid", "Crab_tanner", "Crab_king", "Crab_other", "Shrimp_pandalid", "Shrimp_other", "Epibenthic_carn", 
                   "Epibenthic_graze", "Corals", "Sponges", "Filter_feeders", "Bivalves", "Benthic_carnivores", "Meiobenthos", 
                   "Deposit_feeders", "Macroalgae", "Euphausiids", "Macrozooplankton", "Mesozooplankton", "Microzoplankton", 
                   "Jellyfish", "Gelatinous_other", "Peropods", "Diatoms", "Picophytoplankton")

invert.codes <- c("OCT", "SQD", "TAN", "KIN", "CRO", "PAN", "PWN", "EBC", "BG", "COR", "SPG", "BFF", "BIV", "BC", "BO", "BD", "MA", 
                  "EUP", "ZL", "ZM", "ZS", "JEL", "GEL", "PTE", "PL", "PS")

# non.coastal.polygons <- c(18,19,20,32,43)

# bgm.info <- read_csv("BGM_boxinfo.csv") %>% 
#   dplyr::rename(BOX_ID=box_id) %>% 
#   mutate(coastal_polygon = if_else(BOX_ID %in% non.coastal.polygons, 0,1)) 

#name bioparamfile
#this copies the header template
file.copy("../data/GOAbioparamheader.txt", "../data/GOAbioparam_test.prm", overwrite = TRUE)

goa.header <- "../data/GOAbioparam_test.prm"

#excel_sheets lists all sheets
flags.vector <- read_xlsx("../data/GOA_biological_parameters_fg.xlsx",sheet = "flags") %>% 
  dplyr::select(-Comment) %>% 
  mutate(flag = paste(Flag, Value,Description, sep = " ")) %>% 
  .$flag

#flags in prm
cat("## Flags and switches\n \n",file=goa.header, append=TRUE)
cat(flags.vector, sep = "\n \n", file=goa.header, append=TRUE)
cat("\n \n",file=goa.header, append=TRUE)

#template with info on each parameter, it is a vector, with units and description
# ALBI: most likely here we will have issues, e.g. mismatch in lengths, missing params, etc.
prm.template <- read_xlsx("../data/GOA_biological_parameters_fg.xlsx",sheet = "template") %>% 
  mutate(index = 1:nrow(.))

prm.rows <- prm.template$index

#this is data not in vectors
#The result of using IFERROR is that NA's are read as characters, and so are the non-NA's. See how this propagates, if we are
#only querying and writing to text file it should hardly matter, but if there is any calculation involved we need to go
# back and rewrite formulas so that NA's are NA's
prm.parameters <- read_xlsx("../data/GOA_biological_parameters_fg.xlsx",sheet = "parameters")

migratory.groups <- prm.parameters %>% 
  dplyr::select(atlantis_fg, flagXXXMigrate) %>% 
  filter(!is.na(flagXXXMigrate)) %>% 
  filter(flagXXXMigrate>0) %>% .$atlantis_fg


write.prm <- function(eachrow){
  
  # eachrow <- 97 # for testing
  
  print(paste("Analyzing row from template",eachrow))
  #get info on parameter to analyze, from sheet template
  this.row <- prm.template[eachrow,]
  
  print(this.row)
  
  this.description <- paste("#",this.row$note)
  
  this.parameter <- this.row$parameter
  
  #for parameters included in spreadsheet
  if(this.row$spreadsheet==1){
    
    print(paste(this.parameter,"parameter estimated in spreadsheet"))
    
    if(this.row$vector==0) {
      
      print(paste(this.parameter,"this is not a vector"))
      
      this.data <- prm.parameters %>% 
        dplyr::select(atlantis_fg,this.parameter) %>% 
        setNames(c("fg_group","parameter")) %>% 
        mutate(parameter = as.numeric(parameter)) %>% 
        filter(!is.na(parameter)) 
      
      group.list <- this.data$fg_group %>% as.list()
      
      param.vector <-  rapply(group.list, function (x) gsub("XXX",x,this.parameter), how = "replace") %>% 
        unlist() %>% 
        as.data.frame() %>% 
        setNames("param_name") %>% 
        bind_cols(this.data) %>% 
        dplyr::select(-fg_group) %>% 
        mutate(variable = paste(param_name,parameter, sep = "   ")) %>% 
        dplyr::select(variable) %>% .$variable
      
      cat(this.description, sep="\n",file=goa.header, append=TRUE)
      cat("\n",file=goa.header, append=TRUE)
      cat(param.vector, sep="\n",file=goa.header, append=TRUE)
      cat("\n",file=goa.header, append=TRUE)
      
    } else if (this.row$vector==1){
      
      print(paste(this.parameter,"this is a vector"))
      
      # ALBI: if value1 is NA, none of the vector will be read. This makes sure we do not read empty inverts params or similar
      vector.data <- read_xlsx("../data/GOA_biological_parameters_fg.xlsx",sheet = this.parameter) %>% 
        filter(value1!="NA") %>%  # because of this it should be able to handle character "NA"
        filter(!is.na(value1)) %>% 
        dplyr::select( `Long Name`, index, atlantis_fg, contains("value"))
      
      vector.n <- ncol(vector.data)-3
      
      group.codes <- vector.data %>% 
        dplyr::select(atlantis_fg) %>%
        .$atlantis_fg
      
      cat(this.description, file=goa.header, append=TRUE)
      cat("\n",file=goa.header, append=TRUE)
      
      for(eachgroup in group.codes){
        
        positions <- c(4:ncol(vector.data))
        
        vector.dat <-  vector.data %>% 
          filter(atlantis_fg==eachgroup) %>% 
          select(positions) %>% 
          as.numeric(.) %>% 
          round(.,8)
        
        prm.vector <- vector.dat[!is.na(vector.dat)]
        
        vector.n <- length(prm.vector)
        
        prm.vector.coll <- prm.vector %>% paste(.,collapse =" ")
        
        vector.name <- gsub("XXX",eachgroup,this.parameter) %>% 
          paste(.,vector.n, sep="  ")
        
        cat(trimws(paste0(vector.name, "\n \n", prm.vector.coll,"\n \n", sep = ""),"l"), file=goa.header, append=TRUE)
        cat("\n",file=goa.header, append=TRUE)
      }
      
      
    }
    
  } else if(this.row$spreadsheet==0) {
    
    print(paste(this.parameter,"parameter not estimated in spreadsheet"))
    
    if(this.parameter == "KMIGa_XXX" | this.parameter == "KMIGa_XXXrn" |  this.parameter ==   "KMIGa_XXXsn" | this.parameter == "C_XXX" | this.parameter == "mum_XXX"){
      
      if(this.parameter == "KMIGa_XXX"){
        
        vector.data <- num.biomass.age %>% 
          dplyr::select(Code,age_class,nums_out_mig) %>% 
          setNames(c("Code","age_class","param")) %>% 
          mutate(param=if_else(Code %in% migratory.groups, param, 0))
        
      } else if(this.parameter == "KMIGa_XXXrn") {
        
        vector.data <- num.biomass.age %>% 
          dplyr::select(Code,age_class,resN_mg) %>% 
          setNames(c("Code","age_class","param")) %>% 
          mutate(param=if_else(Code %in% migratory.groups, param, 0))
        
      } else if(this.parameter == "KMIGa_XXXsn") {
        
        vector.data <- num.biomass.age %>% 
          dplyr::select(Code,age_class,struct_N_mg) %>% 
          setNames(c("Code","age_class","param")) %>% 
          mutate(param=if_else(Code %in% migratory.groups, param, 0))
        
      } else if(this.parameter == "C_XXX") {
        
        vector.data <- num.biomass.age %>% 
          dplyr::select(Code,age_class,clearance) %>% 
          setNames(c("Code","age_class","param"))
        
      } else if(this.parameter == "mum_XXX") {
        
        vector.data <- num.biomass.age %>% 
          dplyr::select(Code,age_class,mum_growth_rate) %>% 
          setNames(c("Code","age_class","param"))        
      }
      
      group.code <- num.biomass.age %>% 
        distinct(Code) %>%
        .$Code
      
      cat(this.description, file=goa.header, append=TRUE)
      cat("\n",file=goa.header, append=TRUE)
      
      for(eachgroup in group.code){
        
        prm.vector <- vector.data %>% 
          filter(Code==eachgroup, age_class!=0) %>% 
          dplyr::select(-Code) %>% 
          spread(age_class, param) %>% 
          round(.,8) %>% 
          as.character() %>% 
          paste(.,collapse =" ")
        
        vector.n <- vector.data %>% 
          filter(Code==eachgroup, age_class!=0) %>% 
          .$age_class %>% max()
        
        vector.name <- gsub("XXX",eachgroup,this.parameter) %>% 
          paste(.,vector.n, sep="  ")
        
        cat(trimws(paste0(vector.name, "\n \n", prm.vector,"\n \n", sep = ""),"l"), file=goa.header, append=TRUE)
        cat("\n",file=goa.header, append=TRUE)
      }
      
      
      
    }
    
    if(this.parameter == "FXXX_S#" | this.parameter == "FXXX_S#juv"){
      
      #spatial distribution is for vertebrates and invertebrates that move
      group.code <- grp.file %>%
        filter(!Code %in% c("BB","PB","MA","DL","DF","DC")) %>% # so I guess we do not give distributions for detritus but we do for plankton?
        distinct(Code) %>%
        .$Code
      
      cat(this.description, file=goa.header, append=TRUE)
      cat("\n",file=goa.header, append=TRUE)
      
      for(eachgroup in group.code){
        
        # eachgroup <- group.code[70] # for testing
        
        this.name <- grp.file %>% filter(Code==eachgroup) %>% .$`Name`
        this.longname <- grp.file %>% filter(Code==eachgroup) %>% .$`LongName`
        print(this.name)
        
        if(this.name %in% invert.groups){
          
          this.distribution <- invert.distribution %>% 
            dplyr::select(contains(this.name)) %>% 
            setNames(c("S_1","S_2","S_3","S_4"))
          
        } else {
          
          this.life.history <- life.history %>% 
            filter(Code==eachgroup)
          
          if(this.life.history$Distribution == "Seasonal"){
            
            if(this.parameter=="FXXX_S#juv"){
              
              this.distribution <- seasonal.distribution %>% 
                dplyr::select(contains(this.name)) %>%
                dplyr::select(contains('_J_'))
              
            } else {
              
              this.distribution <- seasonal.distribution %>% 
                dplyr::select(contains(this.name)) %>%
                dplyr::select(contains('_A_'))
              
            }
            
            
          }
        }
        
        # name.vector <- names(this.distribution)  %>% 
        #   gsub(paste(this.longname,"_",sep=""),"",.)
        # 
        # this.distribution <- this.distribution %>% setNames(name.vector)
        
        dist.col <- 1:ncol(this.distribution)
        
        for(eachcol in dist.col) {
          
          #eachcol <- 1 for testing
          
          prm.data <- this.distribution[,eachcol] %>% 
            mutate(index = 1:nrow(.)) %>% 
            spread(index,1)
          
          prm.vector <- prm.data %>% 
            as.numeric %>% 
            round(.,8) %>% 
            paste(.,collapse =" ") 
          
          vector.n <- ncol(prm.data)
          
          vector.name <- gsub("XXX",eachgroup,this.parameter) %>% 
            gsub("#",eachcol,.) %>% 
            paste(.,vector.n, collapse="   ")
          
          cat(trimws(paste(vector.name, "\n \n", prm.vector,"\n \n", sep = ""),"l"), file=goa.header, append=TRUE)
          cat("\n",file=goa.header, append=TRUE)
        }
        
      }
      
    }
    
    if(this.parameter == "pPREY") {
      
      cat(this.description, file=goa.header, append=TRUE)
      cat("\n",file=goa.header, append=TRUE)
      
      
      for(eachrow in 1:nrow(prey.matrix)) {
        
        pprey.name <- prey.matrix %>% 
          slice(eachrow) %>% 
          dplyr::select(name) %>% .$name
        
        sp.row <- prey.matrix %>% 
          slice(eachrow) %>% 
          dplyr::select(-name) %>% 
          as.vector %>% unlist
        
        cat(pprey.name, file=goa.header, append=TRUE)
        cat("\n", file=goa.header, append=TRUE)
        cat(sp.row, file=goa.header, append=TRUE)
        cat("\n \n", file=goa.header, append=TRUE)
        
      }
    }
    
  }
  
}




#cannot be paralellized, the parameters come out jumbled
lapply(prm.rows,write.prm)



