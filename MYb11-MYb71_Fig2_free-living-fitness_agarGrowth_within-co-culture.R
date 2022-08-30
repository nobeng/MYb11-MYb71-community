########################################
#            MYb11 - MYb71             #
#     Free-living fitness (NGM agar)   #
#            Nancy Obeng               #
########################################

#### Set dependencies + open data ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

library(plyr); library(stringr); library(ggplot2); library(dplyr); library(multcomp)

source("my_theme.R")
#colsTreat <- c("grey25", rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))
#colsMYb <- c("#FC2F62", "#189E00", "grey27")
colsMYb <- c("#B9006A", "#139150", "grey27")

#### 1. Open data ####
data <- read.csv("./data/MYb11-MYb71_growthAgar_absData_co-culture.csv", header =T)

#### 2. Plotting #### 
# Proper axis labels
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)}

# focus data
data <- data[data$treatE == "anc",]

#### Plot growth over time ####

# order factor levels
data$tp <- factor(data$tp, levels=c("t0h", "24h", "hCy", "fCy"))
data$strain <- factor(data$strain, levels = c("MYb11", "MYb71", "co-culture"))
data$pop <- as.factor(data$pop)


# Add time points in hours to dataframe
timepoints <- c("t0h", "24h", "hCy", "fCy")
hours <- c(0, 24, 72, 168)

data$hours <- rep(1, length(data$tp))
for(n in 1:length(timepoints)){data[data$tp==timepoints[n], 'hours'] <- hours[n]  }

# Leave 20% space on top of max. data point for significance stars later
top <- max(data$cfu)*2

# Plot data
ggplot(data, aes(x=hours, y=cfu, fill = as.factor(strain), group = interaction(strain, hours)))+
  stat_summary(fun.y = median, geom = "line", aes(group=strain, col = strain), lty=2, size=1)+
  geom_boxplot(size = 0.3)+
  geom_jitter(position=position_dodge(width=15),aes(group=strain), size=1.3, alpha = 0.5)+
  geom_point(aes(x=1, y=top), alpha=0)+
  scale_y_log10(labels=fancy_scientific)+
  labs(x = "Time (h)", y = "CFU/plate")+ 
  scale_fill_manual(values = colsMYb)+
  scale_color_manual(values = colsMYb)+
  my_theme()

# Save plot to file
absPlots <- paste(currentDirectory, "/plots", sep="")
ifelse(dir.exists(absPlots) == F, 
       dir.create(absPlots), "Directory for absolute data already exists.")
fileExtensions <- c(".png", ".svg")
for(i in 1:length(fileExtensions)){
  ggsave(paste(absPlots, "/cfuAgar~time_ANCESTORS_within-co-culture", fileExtensions[i], sep=""),
         dpi = 300, height = 3.2, width = 5.5)  }

#### Stats: abundances MYb11 vs MYb71 per time point ####
# Prepare output file name
fileName <- paste("stats/", "cfuAgar_3d_ANCESTORS_paired-t-test_within-co-culture",".txt", sep="")

# Adjust factor levels of strains to test mono against co-culture abundances
data$strain <- factor(data$strain, levels = c("MYb11", "MYb71"))

sink(fileName)
print("t = 0h")
t.test(log10(data[data$strain == "MYb11" & data$tp == "t0h", 'cfu']),
       log10(data[data$strain == "MYb71" & data$tp == "t0h", 'cfu']),
       paired = T)

print("t = 24h")
t.test(log10(data[data$strain == "MYb11" & data$tp == "24h", 'cfu']),
       log10(data[data$strain == "MYb71" & data$tp == "24h", 'cfu']),
       paired = T)

print("t = hCy")
t.test(log10(data[data$strain == "MYb11" & data$tp == "hCy", 'cfu']),
       log10(data[data$strain == "MYb71" & data$tp == "hCy", 'cfu']),
       paired = T)

print("t = fCy")
t.test(log10(data[data$strain == "MYb11" & data$tp == "fCy", 'cfu']),
       log10(data[data$strain == "MYb71" & data$tp == "fCy", 'cfu']),
       paired = T)
sink()

# Final clean up
rm(data, absPlots, colsMYb, currentDirectory, fileExtensions, hours ,i, n, timepoints, top, 
   fancy_scientific, my_theme)
