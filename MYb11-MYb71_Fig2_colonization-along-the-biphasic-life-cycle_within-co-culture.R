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
assayNames <- c("col.2h", "col.l4", "releas", "col.pe")

# Extract EE cycles and treatments
cycs <- c(4, 10)
treats <- c("bi", "mono")

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

# Plotting absolute load
#colsMYb <- c("grey27", "#FC2F62", "#189E00")
colsMYb <- c("grey27","#B9006A", "#139150")
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
plotLoad <- function(data, cycle, cols, dimW, dimH, lawns){
  # Set cycle and subset respective cycle and lawn
  ifelse(cycle == 10, cyc <- 4, cyc <- 10)
  data <- data[data$cycle != cyc & data$lawn == lawns, ]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data$cfu_w)*2
  
  # Print plot of total co-culture
  dataSum <- data %>%
             group_by(treatE, cycle, repO,lawn) %>% 
             summarise(cfu_w = sum(cfu_w))
  
  dataSum$strain <- rep("Co-culture (MYb11 + MYb71)", length(dataSum$treatE))
  
  # Plot total load of co-culture
  p1 <- ggplot(dataSum, aes(x = toupper(treatE), y = cfu_w)) +
    geom_boxplot()+
    geom_point(aes(fill = "black"), size = 2, alpha = 0.5)+
    geom_point(aes(x=1, y=top), alpha=0)+
    xlab("Evo. life cycle") + ylab("CFU/worm")+
    scale_y_log10(labels = fancy_scientific)+
    scale_fill_manual(values = colsMYb) +
    facet_wrap(~strain)+
    my_theme(d); p1
  
  # Plot load of species in co-culture in ancestral mix
  p2 <- ggplot(data[data$treatE == "anc", ], aes(x = strain, y = cfu_w)) +
    geom_boxplot(aes(fill = strain), size = 0.7)+
    geom_jitter(aes(fill = "black"), size = 2, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_y_log10(labels = fancy_scientific)+
    xlab("Bacteria") + ylab("CFU/worm")+
    scale_fill_manual(values = colsMYb) +
    my_theme(d);p2
  
  # Check for differences in strains within host accordingly:paired test & save to file
  # Prepare file name
  dirStatsOut <- "./stats"
  fileName <- paste(dirStatsOut, "/composition_co-culture_C10_pairedTests_",assayName,  "_lawn", toupper(lawns),".txt", sep="")
  
  testData <- data[data$treatE == "anc", ]
  
  # Where there is only one strain, set other strain to log10=0 (i.e. single place holder cell)
  ifelse(assayName == "col.2h", 
         testData <- rbind(testData,
                           data.frame(assay="col.2h",
                                      treatE="anc",
                                      cycle=0,
                                      strain="MYb71",
                                      repO=6,
                                      lawn="B",
                                      cfu_w=1)),
         0)
  
  # Run tests and save to file alongside
  sink(fileName)
  
  print("Paired t-test")
  print(t.test(log10(testData[testData$strain == "MYb11", 'cfu_w']), 
               log10(testData[testData$strain == "MYb71", 'cfu_w']),
               paired = T, alternative = "two.sided"))
  print("")
  print("Paired Mann-U-Whitney")
  print(wilcox.test(log10(testData[testData$strain == "MYb11", 'cfu_w']), 
                    log10(testData[testData$strain == "MYb71", 'cfu_w']),
                    paired = T, alternative = "two.sided"))
  # Close file
  sink()
  
  # Old version
  p3 <- ggplot(data, aes(x = toupper(treatE), y = cfu_w)) +
    geom_boxplot(aes(fill = strain), size = 0.7)+
    geom_jitter(aes(fill = "black"), size = 2, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_y_log10(labels = fancy_scientific)+
    xlab("Life cycle") + ylab("CFU/worm")+
    facet_wrap(~strain)+
    scale_fill_manual(values = colsMYb) +
    my_theme(d);p3
  
  # Check presence of directory for plots
  absPlots <- paste(currentDirectory, "/plots", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){

    fileName1 <- paste(absPlots, "/", assayName, "_abs.load~treatE_C", cycle, "_sumCoculture",
                       "_lawn", toupper(lawns), fileExtensions[i], sep = "")    
    fileName2 <- paste(absPlots, "/", assayName, "_abs.load~treatE_C", cycle, "_withinCoculture_ANC",
                      "_lawn", toupper(lawns), fileExtensions[i], sep = "")
    
    # Save plot to directory with rel. loads
    ggsave(fileName1, p1, height=dimH, width=dimW,dpi=300)
    ggsave(fileName2, p2, height=dimH, width=dimW,dpi=300) }
  
  return("Plots saved to file")}

plotReleas <- function(data, cycle, cols, dimW, dimH, lawns){
  # Set cycle and subset respective cycle and lawn
  ifelse(cycle == 10, cyc <- 4, cyc <- 10)
  data <- data[data$cycle != cyc & data$lawn == lawns, ]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(data$cfu_releas)*2
  
  # Print plot of total co-culture
  dataSum <- data %>%
    group_by(treatE, cycle, repO,lawn) %>% 
    summarise(cfu_releas = sum(cfu_releas))
  
  dataSum$strain <- rep("Co-culture (MYb11 + MYb71)", length(dataSum$treatE))
  
  # Where there is only one strain, set other strain to log10=0 (i.e. single place holder cell)
  data[data$cfu_releas <= 0, 'cfu_releas'] <- 1
  
  # Plot total load of co-culture
  p1 <- ggplot(dataSum, aes(x = toupper(treatE), y = cfu_releas)) +
    geom_boxplot(size = 0.7)+
    geom_point(aes(fill = "black"), size = 2, alpha = 0.5)+
    geom_point(aes(x=1, y=top), alpha=0)+
    xlab("Evo. life cycle") + ylab("CFU/worm")+
    scale_y_log10(labels = fancy_scientific)+
    scale_fill_manual(values = colsMYb) +
    facet_wrap(~strain)+
    my_theme(d); p1
  
  # Release in co-culture of ancestral mix
  p2 <- ggplot(data[data$treatE == "anc", ], aes(x = strain, y = cfu_releas)) +
    geom_boxplot(aes(fill = strain), size =0.7)+
    geom_jitter(aes(fill = "black"), size = 2, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_y_log10(labels = fancy_scientific)+
    xlab("Bacteria") + ylab("CFU released/worm")+
    scale_fill_manual(values = colsMYb) +
    my_theme(d);p2
  
  # Check for differences in strains within host accordingly:paired test & save to file
  # Prepare file name
  dirStatsOut <- "./stats"
  fileName <- paste(dirStatsOut, "/composition_co-culture_C10_pairedTests_",assayName,  "_releas_lawn", toupper(lawns),".txt", sep="")
  
  testData <- data[data$treatE == "anc", ]
  
  # Where there is only one strain, set other strain to log10=0 (i.e. single place holder cell)
  testData[testData$cfu_releas <= 0, 'cfu_releas'] <- 1
  
  # Run tests and save to file alongside
  sink(fileName)
  
  print("Paired t-test")
  print(t.test(log10(testData[testData$strain == "MYb11", 'cfu_releas']), 
               log10(testData[testData$strain == "MYb71", 'cfu_releas']),
               paired = T, alternative = "two.sided"))
  print("")
  print("Paired Mann-U-Whitney")
  print(wilcox.test(log10(testData[testData$strain == "MYb11", 'cfu_releas']), 
                    log10(testData[testData$strain == "MYb71", 'cfu_releas']),
                    paired = T, alternative = "two.sided"))
  # Close file
  sink()
  
  # Old version
  p3 <- ggplot(data, aes(x = strain, y = cfu_releas)) +
    geom_boxplot(aes(fill = strain), size = 0.7)+
    geom_jitter(aes(fill = "black", shape = as.factor(data$repO)), size = 2, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    scale_y_log10(labels = fancy_scientific)+
    xlab("Life cycle") + ylab("CFU released/worm")+
    facet_wrap(~toupper(treatE))+
    scale_fill_manual(values = colsMYb) +
    scale_shape_manual(values = c(7, 21:25))+
    my_theme(d)+
    theme(axis.text.x = element_text(angle = 45, hjust=1));p3
  
  # Check presence of directory for loads
  absPlots <- paste(currentDirectory, "/plots", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".png", ".svg")
  for(i in 1:length(fileExtensions)){
    fileName1 <- paste(absPlots, "/", assayName, "_releas~treatE_C",cycle, "_sumCoculture",
                       "_lawn", toupper(lawns), fileExtensions[i], sep = "")
    fileName2 <- paste(absPlots, "/", assayName, "_releas~treatE_C",cycle, "_withinCoculture",
                      "_lawn", toupper(lawns), fileExtensions[i], sep = "")
    
    # Save plot to directory with rel. loads
    ggsave(fileName1, p1, height=dimH, width=dimW,dpi=300)
    ggsave(fileName2, p2, height=dimH, width=dimW,dpi=300) }
  
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
      plotLoad(data, cycs[i], colsMYb, 2.7, 3.2, lawns[l]) # sizing for ancestral plots
      
      # Also plot CFU released, when appropriate
      if(assayName == "releas") {
        plotReleas(data, cycs[i], colsMYb, 2.7, 3.2, lawns[l]) # sizing for ancestral plots
      }}}
  
}

#### Plot species proportions in co-culture ###

dataAllAssays <- data.frame()

for (a in 1:length(assayNames)){
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
  
  # Adjust column names of release assay
  ifelse(assayName == "releas",
          # release assay: format subassays of short-term persistence and release
          {dataShortPersist <- subset(data, select = -c(cfu_releas))
           dataShortPersist$assay <- rep("shortPersist", length(dataShortPersist$assay))
           dataReleas <-subset(data, select = -c(cfu_w)) 
           colnames(dataReleas)[colnames(dataReleas) == "cfu_releas"] <- "cfu_w"
           
           dataAllAssays <- rbind(dataAllAssays, dataShortPersist)
           dataAllAssays <- rbind(dataAllAssays, dataReleas)}, 
           
          # all other phenotypes
          dataAllAssays <- rbind(dataAllAssays, data))
}

# Clean up 
rm(data, dataReleas, dataShortPersist)

# Focus on ancestral
dataAllAssays <- dataAllAssays[dataAllAssays$treatE == "anc", ]
dataAllAssays$assaySpec <- paste(dataAllAssays$assay, dataAllAssays$lawn, sep=".")

# Exclude negative entries
dataAllAssays <- dataAllAssays[dataAllAssays$cfu_w >= 0, ]

# Rename phenotypes for better graph readibility
oldNames <- c("col.pe.b", "col.pe.h", "col.pe.e", "releas.b", "shortPersist.b",
              "col.l4.b","col.2h.b")

newNames <- c("Colonization (Ad)", "Persistence (Ad; +hkOP50)", "Persistence (Ad; empty)",
              "Short-term persistence (L4)", "Release (L4)", "Colonization (L4)",
              "Early colonization")

for(i in 1:length(oldNames)){
  dataAllAssays$assaySpec[dataAllAssays$assaySpec == oldNames[i]] <- newNames[i]}

dataAllAssays$assaySpec <- factor(dataAllAssays$assaySpec, levels =
                                 c("Early colonization","Colonization (L4)", "Colonization (Ad)", 
                                   "Short-term persistence (L4)", "Persistence (Ad; +hkOP50)", 
                                   "Persistence (Ad; empty)", "Release (L4)"))

# Calculate total co-culture and proportions
dataSum <- dataAllAssays %>%
           group_by(assay, treatE, cycle,repO,lawn) %>%
           summarise(total_cfu_w = sum(cfu_w))

dataAllAssays <- merge(dataAllAssays, dataSum, by = c("assay", "treatE", "cycle","repO","lawn"))
dataAllAssays$prop <- dataAllAssays$cfu_w / dataAllAssays$total_cfu_w

means.sem <- ddply(dataAllAssays, c("assaySpec", "strain"), summarise,
                   mean=mean(prop), sem=sd(prop)/sqrt(length(prop)))
means.sem <- transform(means.sem, lower=mean-sem, upper=mean+sem)
colnames(means.sem)[colnames(means.sem) == "mean"] <- "prop"
means.sem$assaySpec <- as.factor(means.sem$assaySpec)

# Plot proportions
ggplot(dataAllAssays, aes(x = strain, y = prop, fill =  strain))+
  geom_boxplot()+
  #stat_summary(fun.y = "mean",geom = "bar")+
  geom_point(alpha = 0.5)+
  geom_errorbar(data=means.sem, aes(x = strain, ymax=upper,  ymin=lower), width=0.15) + 
  labs(x = "Bacterial species", y = "Porportion in co-culture")+ 
  facet_wrap(~assaySpec)+
  scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))+
  scale_fill_manual(values = colsMYb[2:3])+
  my_theme()

# Clean up
rm(a, assayName,i, l, plotLoad, plotReleas)