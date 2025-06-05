########################################
#            MYb11 - MYb71             #
#   Within host fitness (CFU/L4 worm)  #
#            Nancy Obeng               #
########################################

#### Set dependencies + open data ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Specify libraries
library(plyr); library(ggplot2); library(lawstat); library(lme4); library(emmeans); library(multcomp)

# Source files
source("my_theme.R")

# Specify assay name
assayName <- "col.l4"

# Open file of absolute cfus in colonization assay
data <- read.csv("./data/MYb11-MYb71_col.l4_absData.csv", header = T)
data <- subset(data, select=-id) # remove id

# Extract EE cycles and (possible) strains and lawn types
cycs <- c(4, 10)
treats <- c("bi", "mono")
strains <- unique(data$strain)
lawns <- unique(data$lawn)

#### Function definitions ####

# Averaging technical replicates
averageReplicates <- function(data){
  ifelse(assayName == "releas",
         data <- ddply(data, .(assay, treatE, cycle, strain, repO, lawn), summarize, 
                       cfu_w = mean(cfu_w), 
                       cfu_releas = mean(cfu_releas)), 
         data <- ddply(data, .(assay, treatE, cycle, strain, repO, lawn), summarize, 
                           cfu_w = mean(cfu_w)))
  return(data) }

# Plotting absolute load
#colsMYb <- c("#FC2F62", "#189E00","grey27")
colsMYb <- c("#B9006A", "#139150", "grey27")
colsTreat <- c("grey25", rgb(0, 158,115, maxColorValue = 255), rgb(230, 159,0, maxColorValue = 255))

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

# Plot load
plotLoad <- function(data, cycle, cols, strains, dimW, dimH, lawns){
  # Set cycle and subset respective cycle
  data <- data[data$cycle == 0, ]
  data <- data[!is.na(data$cfu_w),]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data$cfu_w)*2
  
  # Order factor levels
  data$strain <- factor(data$strain, levels = c("MYb11", "MYb71", "co-culture"))
  
  # Print plot
  p <- ggplot(data, aes(x = strain, y = cfu_w, fill = strain)) +
    geom_boxplot(size = 0.7)+
    geom_jitter(size = 2.5, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_fill_manual(values = colsMYb) +
    scale_y_log10(labels = fancy_scientific)+
    xlab("Bacteria") + ylab("CFU/worm")+
    my_theme(d)+
    theme(axis.text.x = element_text(hjust=1, angle=45));p
  
  # Check presence of directory for rel.loads
  absPlots <- paste(currentDirectory, "/plots", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    fileName <- paste(absPlots, "/", assayName, "_abs.load~treatE_ANCESTORS", fileExtensions[i], sep = "")
  
    # Save plot to directory with rel. loads
    ggsave(fileName, p, height=3.2, width=2.7,dpi=300) }
  
  return(p)}

#### = 1. Average replicates + plot  = ####
# Average technical replicates
data <- averageReplicates(data)

data <- data[data$strain!= "OP50",]

# Plot 
for (l in 1:length(lawns)){ # Loop through lawn types
  for (i in 1:length(cycs)){ # Loop through evol. cycles
    plotLoad(data, cycs[i], colsMYb, "ancestors", 3, 2.5, lawns[l]) }}


#### 2. Check difference between mono and co-culture (Dunnett tests) ####
# Collect data
data <- data[data$cycle == 0, ]
data <- data[!is.na(data$cfu_w),]

# Prepare output file name
fileName <- paste("stats/", "col.l4_ANCESTORS_ANOVA+Dunnet",".txt", sep="")

# Adjust factor levels of strains to test mono against co-culture abundances
data$strain <- factor(data$strain, levels = c("co-culture", "MYb11", "MYb71"))
data$strain <- factor(data$strain, levels = c("MYb11", "MYb71", "co-culture"))

# Run linear regression (LM with ANOVA)
m <- lm(log10(cfu_w) ~ strain*repO, data=data)

sink(fileName)
print(as.data.frame(anova(m))) # ANOVA output
print("")
print(summary(glht(m, mcp(strain="Tukey")))) # Post-hoc output: Tukey
print("")
print(summary(glht(m, mcp(strain="Dunnett")))) # Post-doc output: Dunnett
sink()


# Directly compare MYb11 and MYb71
t.test(log10(data[data$strain == "MYb11", 'cfu_w']), log10(data[data$strain == "MYb71", 'cfu_w']))

# Clean up
rm(m, data, assayName, colsMYb, colsTreat, currentDirectory, cycs, fileName, i, l, 
   lawns, strains, treats, averageReplicates, fancy_scientific, my_theme, plotLoad)
