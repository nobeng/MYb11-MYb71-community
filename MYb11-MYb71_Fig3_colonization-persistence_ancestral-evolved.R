########################################
#            MYb11 - MYb71             #
#   Abundances in co-culture along     #
#       the biphasic life cycle        #
#            Nancy Obeng               #
########################################

#### Set dependencies + open data ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Specify libraries
library(plyr); library(ggplot2); library(lawstat); library(lme4); library(emmeans); library(dplyr)

# Source files
source("my_theme.R")

# Specify assay name
assayNames <- c("col.2h","col.l4", "releas","col.pe") 

# Extract EE cycles and treatments
cycs <- c(4, 10)
treats <- c("bi", "mono")

colsMYb <- c("#B9006A", "#139150", "grey27")
colsEvol <- c("black", "grey27", "grey60","#B9006A","#FF5F90","#139150", "#1AC96F")

#### = Function definitions = ####

# Averaging technical replicates
averageReplicates <- function(data){
  ifelse(assayName == "releas",
         data <- ddply(data, .(assay, treatE, cycle, strain, repO, lawn), summarize, 
                       cfu_w = mean(cfu_w), 
                       cfu_releas = mean(cfu_releas)), 
         data <- ddply(data, .(assay, treatE, cycle, strain, repO, lawn), summarize, 
                       cfu_w = mean(cfu_w)))
  
  # Average when multiple counts are present
  
  return(data) }

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
plotLoad <- function(data, cycle, cols, dimW, dimH, lawns, treatment){
  # Set cycle, evol. treatment of focus and subset respective cycle and lawn
  ifelse(cycle == 10, cyc <- 4, cyc <- 10)
  ifelse(treatment == "bi", excludedTreatE <- "mono", excludedTreatE <- "bi")
  data <- data[data$cycle != cyc & data$lawn == lawns, ]
  
  # Collect totals of co-cultures
  dataSum <- data %>%
    group_by(assay, treatE, cycle, repO,lawn) %>% 
    summarise(cfu_w = sum(cfu_w))
  dataSum$strain <- rep("co-culture", length(dataSum$assay))
  
  # Rearrange data frame to attach co-culture sums
  strain <- data$strain
  data <- subset(data, select=-strain)
  data <- cbind(data, strain)
  data <- rbind(data, dataSum)
  
  data$strain <- factor(data$strain, levels = c("MYb11", "MYb71", "co-culture"))
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data$cfu_w)*2
  
  # Add sample identifier for proper color scheme
  data$sample <- paste(data$strain, data$treatE, sep="")
  
  # Plot absolute counts of species within co-culture (anc vs evolved)
  p <- ggplot(data[data$treatE != excludedTreatE,], aes(x = strain, y = cfu_w)) +
    geom_boxplot(aes(fill = sample), size = 0.7)+
    geom_point(aes(fill = "black", group=treatE), size = 2, alpha = 0.5, 
                position=position_dodge(width=0.75))+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_y_log10(labels = fancy_scientific)+
    xlab("Evol. treatment") + ylab("CFU/worm")+
    scale_fill_manual(values = colsEvol) +
    my_theme(d)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  # Check presence of directory for plots
  absPlots <- paste(currentDirectory, "/plots/fig.3", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){

    fileName1 <- paste(absPlots, "/", assayName, "_abs.load~treatE_C", cycle,"_",
                       toupper(treatment),"_components-co-culture",
                       "_lawn", toupper(lawns), fileExtensions[i], sep = "")    
    
    # Save plot to directory with rel. loads
    ggsave(fileName1, p, height=dimH, width=dimW,dpi=300) }
  
  return("Plots saved to file")}

plotReleas <- function(data, cycle, cols, dimW, dimH, lawns, treatment){
  # Set cycle, evol. treatment of focus and subset respective cycle and lawn
  ifelse(cycle == 10, cyc <- 4, cyc <- 10)
  ifelse(treatment == "bi", excludedTreatE <- "mono", excludedTreatE <- "bi")
  data <- data[data$cycle != cyc & data$lawn == lawns, ]
  
  # Collect totals of co-cultures
  dataSum <- data %>%
    group_by(assay, treatE, cycle, repO,lawn) %>% 
    summarise(cfu_w = sum(cfu_w), cfu_releas = sum(cfu_releas))
  dataSum$strain <- rep("co-culture", length(dataSum$assay))
  
  # Rearrange data frame to attach co-culture sums
  strain <- data$strain
  data <- subset(data, select=-strain)
  data <- cbind(data, strain)
  data <- rbind(data, dataSum)
  
  data$strain <- factor(data$strain, levels = c("MYb11", "MYb71", "co-culture"))
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data$cfu_w)*2
  
  # Add sample identifier for proper color scheme
  data$sample <- paste(data$strain, data$treatE, sep="")
  
  # Plot absolute counts of species within co-culture (anc vs evolved)
  p1 <- ggplot(data[data$treatE != excludedTreatE,], aes(x = strain, y = cfu_w)) +
    geom_boxplot(aes(fill = sample), size = 0.7)+
    geom_point(aes(fill = "black", group=treatE), size = 2, alpha = 0.5, 
               position=position_dodge(width=0.75))+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_y_log10(labels = fancy_scientific)+
    xlab("Evol. treatment") + ylab("CFU/worm")+
    scale_fill_manual(values = colsEvol) +
    my_theme(d)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  # Release
  # Where there is only one strain, set other strain to log10=0 (i.e. single place holder cell)
  data[data$cfu_releas <= 0, 'cfu_releas'] <- 1
  p2 <- ggplot(data[data$treatE != excludedTreatE,], aes(x = strain, y = cfu_releas)) +
    geom_boxplot(aes(fill = sample), size = 0.7)+
    geom_point(aes(fill = "black", group=treatE), size = 2, alpha = 0.5, 
               position=position_dodge(width=0.75))+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_y_log10(labels = fancy_scientific)+
    xlab("Evol. treatment") + ylab("CFU released/worm")+
    scale_fill_manual(values = colsEvol) +
    my_theme(d)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  # Check presence of directory for plots
  absPlots <- paste(currentDirectory, "/plots/fig.3", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    
    fileName1 <- paste(absPlots, "/", assayName, "_abs.load~treatE_C", cycle,"_",
                       toupper(treatment),"_components-co-culture",
                       "_lawn", toupper(lawns), fileExtensions[i], sep = "")    
    fileName2 <- paste(absPlots, "/", assayName, "_abs.releas~treatE_C", cycle,"_",
                       toupper(treatment),"_components-co-culture",
                       "_lawn", toupper(lawns), fileExtensions[i], sep = "")   
    
    # Save plot to directory with rel. loads
    ggsave(fileName1, p1, height=dimH, width=dimW,dpi=300)
    ggsave(fileName2, p2, height=dimH, width=dimW,dpi=300)}
  
  return("Plots saved to file")}

#### Plot biological replicate data ####
for(a in 1:length(assayNames)){
  
  # Subset assay
  assayName <- assayNames[a]
  
  # Open file of absolute cfus in colonization assay
  data <- read.csv(paste("./data/", assayName, "_absDataMIX.csv", sep = ""), header = T)
  data <- subset(data, select=-id) # remove id
  
  # Extract strains and lawn types
  strains <- unique(data$strain)
  lawns <- unique(data$lawn)
  
  # Average technical replicates
  data <- averageReplicates(data)
  
  # Plot & save to file
  for (l in 1:length(lawns)){ # Loop through lawn types
    for (i in 1:length(cycs)){ # Loop through evol. cycles
      for (t in 1:length(treats)){
      ifelse(assayName == "releas",
             plotReleas(data, cycs[i], colsMYb, 3.2, 3.2, lawns[l], treats[t]),
             plotLoad(data, cycs[i], colsMYb, 3.2, 3.2, lawns[l], treats[t])) # sizing for ancestral plots
             }}}}

# Clean up
rm(data, a, assayName, assayNames, colsEvol, colsMYb, currentDirectory, cycs, i, l, lawns,
   strains, t, treats, averageReplicates, fancy_scientific, my_theme, plotLoad, plotReleas)