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
data <- read.csv("./data/MYb11-MYb71_growthAgar_absData.csv", header =T)

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

# order factor levels
data$tp <- factor(data$tp, levels=c("t0h", "24h", "hCy", "fCy"))
data$strain <- factor(data$strain, levels = c("MYb11", "MYb71", "co-culture"))
data$pop <- as.factors(data$pop)

# Leave 20% space on top of max. data point for significance stars later
top <- max(data$cfu)*2

# Fig.1: Growth at 3d
p <- ggplot(data[data$tp == "hCy",], aes(x = strain, y = cfu, fill = strain)) +
  geom_boxplot(size = 0.7)+
  geom_jitter(size = 2.5, alpha = 0.5, width = 0.1)+
  geom_point(aes(x=1, y=top), alpha=0)+
  scale_fill_manual(values = colsMYb) +
  scale_y_log10(labels = fancy_scientific)+
  xlab("Bacteria") + ylab("CFU/plate")+
  my_theme(d)+
  theme(axis.text.x = element_text(hjust=1, angle=45));p

# Save plot to file
# Check presence of directory for rel.loads
absPlots <- paste(currentDirectory, "/plots", sep="")
ifelse(dir.exists(absPlots) == F, 
       dir.create(absPlots), "Directory for absolute data already exists.")

# Define file name
fileExtensions <- c(".png", ".svg")
for(i in 1:length(fileExtensions)){
  fileName <- paste(absPlots, "/cfuAgar.3d~treatE_ANCESTORS", fileExtensions[i], sep = "") 
  
  # Save plot to directory with rel. loads
  ggsave(fileName, p, height=3.2, width=2.7,dpi=300) }

#### 3. Check difference between mono and co-culture (Dunnet tests) ####
# Prepare output file name
fileName <- paste("stats/", "cfuAgar_3d_ANCESTORS_ANOVA+Dunnet",".txt", sep="")

# Adjust factor levels of strains to test mono against co-culture abundances
data$strain <- factor(data$strain, levels = c("co-culture", "MYb11", "MYb71"))

# Run linear regression (LM with ANOVA)
m <- lm(log10(cfu) ~ strain*pop, data=data[data$tp == "hCy",])

sink(fileName)
print(as.data.frame(anova(m))) # ANOVA output
print("")
print(summary(glht(m, mcp(strain="Tukey")))) # Post-hoc output: Tukey
print("")
print(summary(glht(m, mcp(strain="Dunnett")))) # Post-doc output: Dunnett
sink()

# Clean up
rm(m, p)

#### 4. Plot for supplement: growth over time ####
# Add time points in hours to dataframe
timepoints <- c("t0h", "24h", "hCy", "fCy")
hours <- c(0, 24, 72, 168)

data$hours <- rep(1, length(data$tp))
for(n in 1:length(timepoints)){data[data$tp==timepoints[n], 'hours'] <- hours[n]  }

# Plot data
ggplot(data, aes(x=hours, y=cfu, fill = as.factor(strain), group = interaction(strain, hours)))+
  stat_summary(fun.y = median, geom = "line", aes(group=strain, col = strain), lty=2, size=1)+
  geom_boxplot(size = 0.3)+
  geom_jitter(position=position_dodge(width=15),aes(group=strain), size=1.3, alpha = 0.5)+
  scale_y_log10(labels=fancy_scientific)+
  labs(x = "Time (h)", y = "CFU/plate")+ 
  scale_fill_manual(values = colsMYb)+
  scale_color_manual(values = colsMYb)+
  my_theme()

# Save plot to file
for(i in 1:length(fileExtensions)){
  ggsave(paste(absPlots, "/cfuAgar~time_ANCESTORS", fileExtensions[i], sep=""),
         dpi = 300, height = 3.2, width = 5.5)  }

# Final clean up
rm(data, absPlots, colsMYb, currentDirectory, fileExtensions, fileName, hours ,i, n, timepoints, top, 
   fancy_scientific, my_theme)
