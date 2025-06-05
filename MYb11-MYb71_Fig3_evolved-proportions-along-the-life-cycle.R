########################################
#            MYb11 - MYb71             #
#          Evolved proportions         #
#     along the biphasic life cycle    #
#            Nancy Obeng               #
########################################

#### = Set R dependencies = ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Specify libraries
library(plyr); library(ggplot2); library(multcomp); library(ggpubr)

# Source files
source("my_theme.R")

#### = Import data = ####

# Import proportion data across assays
fileNameCol <- "./data/MYb11-MYb71_colonization_meanProps.csv"
fileNameGrowth <- "./data/MYb11-MYb71_growth_meanProps.csv"
# in both cases proportion means proportion of MYb11 in co-culture

# Proportions calculated as median proportion across biol. replicates
dataCol <- read.csv(fileNameCol, header = T)
dataGrowth <- read.csv(fileNameGrowth, header = T)
dataGrowth <- dataGrowth[dataGrowth$strain == "MYb11",]
dataGrowth <- subset(dataGrowth, select=-strain)

# Decompose releas assay data as well as different plates
dataCol$assay <- as.character(dataCol$assay)

# 1. Paste proportions of releas into the proportions column
releas <- subset(dataCol[!is.na(dataCol$propReleas), ],)
releasWorm <- subset(releas, select = - propReleas)
releasOut <- subset(releas, select = - prop)

names(releasOut)[names(releasOut) == "propReleas"]  <- "prop"
releasWorm[releasWorm$assay == "releas", 'assay'] <- "Short-term persistence"
releasOut[releasOut$assay == "releas", 'assay'] <- "Release"

dataCol <- subset(dataCol[is.na(dataCol$propReleas), ], select = - propReleas)
dataCol <- rbind(dataCol, releasWorm, releasOut)

rm(releas, releasWorm, releasOut) # Clean up

# 2. Paste persistence data from different environments
dataCol[dataCol$assay == "col.pe" & dataCol$lawn == "b", 'assay'] <- "Colonization (Ad)"
dataCol[dataCol$assay == "col.pe" & dataCol$lawn == "h", 'assay'] <- "Persistence (Ad; +hkOP50)"
dataCol[dataCol$assay == "col.pe" & dataCol$lawn == "e", 'assay'] <- "Long-term persistence" # prev. Persistence (Ad; empty)"

dataCol <- subset(dataCol, select = -lawn)

# 3. Rename assays for better graph readibility
dataCol[dataCol$assay == "col.l4", 'assay'] <- "Colonization"
dataCol[dataCol$assay == "col.2h", 'assay'] <- "Early colonization"
dataGrowth$assay <- as.character(dataGrowth$assay)
dataGrowth[dataGrowth$assay == "growth_t0h", 'assay'] <- "Growth (0h)"
dataGrowth[dataGrowth$assay == "growth_24h", 'assay'] <- "Growth (24h)"
dataGrowth[dataGrowth$assay == "growth_hCy", 'assay'] <- "Growth (72h)"
dataGrowth[dataGrowth$assay == "growth_fCy", 'assay'] <- "Growth (168h)"

# 4. Add growth data
data <- rbind(dataCol, dataGrowth)

# Set focal phenotypes & adjust factor level order
data <- data[!data$assay%in%c("Persistence (Ad; +hkOP50)", "Colonization (Ad)"), ]

data$assay <- factor(data$assay, levels = c("Growth (0h)", "Growth (24h)", "Growth (72h)", "Growth (168h)", 
                                            "Early colonization", "Colonization",
                                            "Short-term persistence", 
                                            "Long-term persistence", "Release"))

data$treatE <- toupper(data$treatE)
data$treatE <- factor(data$treatE, levels = rev(toupper(c("anc", "bi", "mono"))))

data$strain <- rep("MYb11", length(data$assay))

dataMYb71 <- data
dataMYb71$prop <- 1 - dataMYb71$prop
dataMYb71$strain <- rep("MYb71", length(dataMYb71$assay))

dataWithStrain <- rbind(data, dataMYb71)


#### = Plotting = ####

# Possible EE cycles
cycs <- c(4, 10)

# Plotting color scheme
colsMYb <- c("#B9006A", "#139150", "grey27")

# Function to stacked bar graphs
plotStackedBarGraphs <- function(data, cycle){
  # Set cycle to exclude
  ifelse(cycle == 4, cyc <- 10, cyc <- 4)
  
  data <- data[data$cycle!= cyc, ]
  data$treatE <- factor(data$treatE, levels = rev(c("ANC", "BI", "MONO")))
  
  ggplot(data, aes(fill=strain, y=treatE, x=prop)) + 
    geom_bar(position="stack", stat="identity")+
    facet_grid(~assay)+
    scale_fill_manual(values=colsMYb)+
    theme_bw()+
    theme(
      axis.text.x=element_text(hjust=1, angle=45)
    )
  
  ggsave(paste("./plots/evoProportions_stacked-bar_",cycle,".svg", sep=""), height=2, width=14,dpi=300)
  
  return()
}

# Loop through cycles assayed
for (i in 1:length(cycs)){ plotStackedBarGraphs(dataWithStrain, cycs[i]) } 

# Clean up
rm(data, dataCol, dataGrowth, dataMYb71, dataWithStrain, colsMYb, currentDirectory, cycs, fileNameCol, 
   fileNameGrowth, i, my_theme, plotCircles, plotStackedBarGraphs)
