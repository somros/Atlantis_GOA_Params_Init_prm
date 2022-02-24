# debugging make.init.csv

def.att.file <- system.file("extdata", "AttributeTemplate.csv", package = "shinyrAtlantis")
df.atts <- read.csv(file = def.att.file, header = TRUE, stringsAsFactors = FALSE)


# read in group data from group csv file
df.grp <- read.csv(file = grp.file, header = TRUE, stringsAsFactors = FALSE)
# make sure GroupType column title exists
col.titles <- names(df.grp)
col.InvertType <- which(col.titles == "InvertType")
if (!length(col.InvertType) == 0) {
  names(df.grp)[col.InvertType] <- "GroupType"
}
df.grp$GroupType <- as.character(df.grp$GroupType)
# find epibenthos groups
epi.grps.def <- c("SED_EP_FF", "SED_EP_OTHER", "EP_OTHER", "MOB_EP_OTHER",
                  "MICROPHTYBENTHOS", "PHYTOBEN", "SEAGRASS", "CORAL")
epi.grps <- df.grp$Name[df.grp$GroupType %in% epi.grps.def]
df.grp <- df.grp %>% mutate(isEpiGrp = GroupType %in% epi.grps.def)

# find cover groups. These groups need _cover in boxTracers
cover.grps <- df.grp$Name[df.grp$IsCover == 1]

# set up flags for groups that need multiple N values (e.g. _N1, _N2, ...)
df.grp <- df.grp %>% mutate(multiN =
                              (IsCover == 1       & (NumCohorts > 1)) |
                              (GroupType == "PWN" & (NumCohorts > 1)) |
                              (GroupType == "CEP" & (NumCohorts > 1)))

# groups with nums, structural and reserve N values
sr.grps <- c("FISH", "BIRD", "SHARK", "MAMMAL", "REPTILE", "FISH_INVERT")
# set up flags for groups that need _Nums, _ResN, _StructN
df.grp <- df.grp %>% mutate(needsNums = GroupType %in% sr.grps)

# set up a flag for groups that need light adaptation
light.adpn.grps <- c("DINOFLAG", "MICROPHTYBENTHOS", "SM_PHY",
                     "MED_PHY", "LG_PHY")
df.grp <- df.grp %>% mutate(needsLight = GroupType %in% light.adpn.grps)

# create data frame for invert biological variables (ignores multiple stocks)
Variable  <- NULL # variable name
long_name <- NULL # long name
att.index <- NULL # corresponding row of df.atts
for (grp in 1:length(df.grp$Name)) {
  if (!df.grp$needsNums[grp]) {
    if (!df.grp$multiN[grp]) { # single group
      grp <- which(df.grp$Name=="Jellyfish")
      Variable  <- c(Variable, paste(df.grp$Name[grp], "_N", sep = ""))
      indx <- which(df.atts$name==paste(df.grp$GroupType[grp], "_N", sep = ""))
      long_name <- c(long_name, paste(df.grp$Name[grp],
                                      df.atts$long_name[indx], sep = " "))
      att.index <- c(att.index, indx)
    } else { # multiple groups
      for (j in 1:df.grp$NumCohorts[grp]) {
        Variable <- c(Variable, paste(df.grp$Name[grp], "_N", as.character(j),
                                      sep = ""))
        indx <- which(df.atts$name==paste(df.grp$GroupType[grp], "_N",
                                          sep = ""))
        long_name <- c(long_name, paste(df.grp$Name[grp], "cohort",
                                        as.character(j), df.atts$long_name[indx], sep = " "))
        att.index <- c(att.index, indx)
      }
    }
    if (df.grp$IsCover[grp]) { # single cover group
      Variable <- c(Variable, paste(df.grp$Name[grp], "_Cover", sep = ""))
      indx <- which(df.atts$name == "Cover")
      long_name <- c(long_name, paste("Percent cover by",
                                      df.grp$Name[grp], sep = " "))
      att.index <- c(att.index, indx)
    }
    if (df.grp$IsSiliconDep[grp]) { # single silicon group
      Variable <- c(Variable, paste(df.grp$Name[grp], "_S", sep = ""))
      indx <- which(df.atts$name == "Si3D")
      long_name <- c(long_name, paste(df.grp$Name[grp],
                                      "Silicon", sep = " "))
      att.index <- c(att.index, indx)
    }
    if (df.grp$needsLight[grp]) { # single light adaptation group
      Variable <- c(Variable, paste("Light_Adaptn_", df.grp$Code[grp],
                                    sep = ""))
      indx <- which(df.atts$name == "Light3D")
      long_name <- c(long_name, paste("Light adaption of",
                                      df.grp$Name[grp], sep = " "))
      att.index <- c(att.index, indx)
    }
  }
}
df.invert <- data.frame(Variable, long_name, att.index,
                        stringsAsFactors = FALSE)

# create data frame of vertebrate default variables (ignores multiple stocks)
Variable  <- NULL # variable name
long_name <- NULL # Long name
att.index <- NULL # corresponding row of df.atts
for (grp in 1:length(df.grp$Name)) {
  if (df.grp$needsNums[grp]) {
    Variable <- c(Variable, paste(df.grp$Name[grp], "_N", sep = ""))
    indx <- which(df.atts$name==paste(df.grp$GroupType[grp], "_N", sep = ""))
    long_name <- c(long_name, paste(df.grp$Name[grp],
                                    df.atts$long_name[indx], sep = " "))
    att.index <- c(att.index, indx)
    for (j in 1:df.grp$NumCohorts[grp]) {
      Variable <- c(Variable, paste(df.grp$Name[grp], as.character(j),
                                    "_Nums", sep = ""))
      indx <- which(df.atts$name=="Nums3D")
      long_name <- c(long_name, paste("Numbers of", df.grp$Name[grp], "cohort",
                                      as.character(j), sep = " "))
      att.index <- c(att.index, indx)
    }
    for (j in 1:df.grp$NumCohorts[grp]) {
      Variable <- c(Variable, paste(df.grp$Name[grp], as.character(j),
                                    "_StructN", sep = ""))
      indx <- which(df.atts$name=="StructN3D")
      long_name <- c(long_name, paste("Individual structural N for",
                                      df.grp$Name[grp], "cohort", as.character(j), sep = " "))
      att.index <- c(att.index, indx)
    }
    for (j in 1:df.grp$NumCohorts[grp]) {
      Variable <- c(Variable, paste(df.grp$Name[grp], as.character(j),
                                    "_ResN", sep = ""))
      indx <- which(df.atts$name=="ResN3D")
      long_name <- c(long_name, paste("Individual reserve N for",
                                      df.grp$Name[grp], "cohort", as.character(j), sep = " "))
      att.index <- c(att.index, indx)
    }
  }
}

df.vert <- data.frame(Variable, long_name, att.index, stringsAsFactors = FALSE)

df.return <- rbind(df.invert, df.vert)