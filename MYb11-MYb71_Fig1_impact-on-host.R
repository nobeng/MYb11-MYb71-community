########################################
#            MYb11 - MYb71             #
#       Impact on C. elegans MY316     #
#            Nancy Obeng               #
########################################

#### = 1. Dependencies, directories and variables = ####
library(plyr); library(ggplot2); library(gridExtra); library(multcomp)

# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

# Source files
source("my_theme.R")

# Variables
colsTreat <- c("black","#8DB500", "#FFB733")
#colsMYb <- c("#FC2F62", "grey27", "#189E00", "blue")
colsMYb <- c("#B9006A", "grey27", "#139150")
assay <- assayName <- "wormSize"
cycs <- c(4, 10)

#### = 2. Function definitions = ####

# Function to plot worm dimensions as a function of treatE
plot_wormDimension.treatE <- function(dataset, cycle, cols, strains, dims, labelX, plotH, plotW){
  # 1. Set cycle to exclude from data set
  ifelse(cycle == 4, cyc <- 10, ifelse(cycle == 10, cyc <- 4, 0))
  dataset <- dataset[dataset$cycle != cyc,]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(dataset[,dims])*1.05
  
  # 2. Print plot
  p <- ggplot(dataset, aes(x = toupper(treatE), y = dataset[, dims]))+
    geom_boxplot(aes(fill = strain), size = 0.7) +
    geom_jitter(size = 2, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    xlab("Life Cycle") + ylab(labelX)+
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    facet_grid(~strain)+
    my_theme(dataset)
  
  # 3. Save plot to file
  # Check presence of directory
  absPlots <- paste(currentDirectory, "/plots", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".svg", ".png")
  for(i in 1:length(fileExtensions)){
    fileName <- paste(absPlots, "/", assayName, "_", dims,"~treatE_C", cycle, "_", strains,
                      fileExtensions[i], sep = "")
    # Save plot
    ggsave(fileName, p, height=plotH, width=plotW,dpi=300) }
  
  # Check for differences between treatE
  # Prepare output file name
  fileName <- paste("stats/monoVSco_",dims, "_compTREATE_ANOVA+Dunnet",".txt", sep="")
  
  # Adjust factor levels of strains to test mono against co-culture abundances
  dataset$strain <- factor(dataset$strain, levels = c("co-culture", "MYb11", "MYb71"))
  dataset$treatE <- factor(dataset$treatE, levels = c("anc", "mono", "bi"))
  
  # Run linear regression (LM with ANOVA)
  m <- lm(dataset[, dims] ~ treatE*repO, data=dataset)
  
  sink(fileName)
  print(as.data.frame(anova(m))) # ANOVA output
  print("")
  print(summary(glht(m, mcp(treatE="Tukey")))) # Post-hoc output: Tukey
  print("")
  print(summary(glht(m, mcp(treatE="Dunnett")))) # Post-doc output: Dunnett
  sink()
  
  return(p)  }

# Function to plot worm dimensions of ancestors
plot_wormDimension.strainANC <- function(dataset, cycle, cols, strains, dims, labelX, plotH, plotW){
  # 1. Set cycle to exclude from data set
  ifelse(cycle == 4, cyc <- 10, ifelse(cycle == 10, cyc <- 4, 0))
  dataset <- dataset[dataset$cycle != cyc,]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(dataset[,dims])*1.05
  
  dataset$strain <- factor(dataset$strain, levels = c("MYb11", "MYb71", "co-culture"))
  
  # 2. Print plot
  p <- ggplot(dataset[dataset$treatE == "anc", ], aes(x = strain, y = dataset[, dims]))+
    geom_boxplot(aes(fill = strain), size = 0.7) +
    geom_jitter(size = 2.5, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    xlab("Bacteria") + ylab(labelX)+
    scale_fill_manual(values = cols[c(1, 3, 2)]) +
    scale_color_manual(values = cols[c(1, 3, 2)]) +
    #facet_grid(~strain)+
    my_theme(dataset)+
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
  
  # 3. Save plot to file
  # Check presence of directory
  absPlots <- paste(currentDirectory, "/plots", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".svg", ".png")
  for(i in 1:length(fileExtensions)){
    fileName <- paste(absPlots, "/", assayName, "_", dims,"~strain_ANCESTORS_", strains, fileExtensions[i], sep = "")
    # Save plot
    ggsave(fileName, p, height=plotH, width=plotW,dpi=300) }
  
  # Check for differences between either monoculture and the co-culture within each treatE
  # Prepare output file name
  fileName <- paste("stats/monoVSco_", dims, "_compANC_ANOVA+Dunnet",".txt", sep="")
  
  # Adjust factor levels of strains to test mono against co-culture abundances
  dataset$strain <- factor(dataset$strain, levels = c("co-culture", "MYb11", "MYb71"))
  
  # Run linear regression (LM with ANOVA)
  m <- lm(dataset[, dims] ~ strain*repO, data=dataset)
  
  sink(fileName)
  print(as.data.frame(anova(m))) # ANOVA output
  print("")
  print(summary(glht(m, mcp(strain="Tukey")))) # Post-hoc output: Tukey
  print("")
  print(summary(glht(m, mcp(strain="Dunnett")))) # Post-doc output: Dunnett
  sink()
  
  return(p)  }

# Function to plot # washed worms detected as a function of treatE
plot_nWorms.treatE <- function(dataset, cycle, cols, strains,plotH, plotW){
  # 1. Set cycle to exclude from data set
  ifelse(cycle == 4, cyc <- 10, ifelse(cycle == 10, cyc <- 4, 0))
  dataset <- dataset[dataset$cycle != cyc,]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(dataset$nWorms)*1.05
  
  # 2. Print plot
  p <- ggplot(dataset, aes(x = toupper(treatE), y = nWorms))+
    geom_boxplot(aes(fill = strain), size = 0.7) +
    geom_jitter(size = 2, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    xlab("Life Cycle") + ylab("# worms/150?l")+
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    facet_grid(~strain)+
    my_theme(dataset)
  
  # 3. Save plot to file
  # Check presence of directory
  absPlots <- paste(currentDirectory, "/plots", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".svg", ".png")
  for(i in 1:length(fileExtensions)){
    fileName <- paste(absPlots, "/", assayName, "_nWorms~treatE_C", cycle, "_", strains, fileExtensions[i], sep = "")
    
    # Save plot
    ggsave(fileName, p, height=plotH, width=plotW,dpi=300) }
  
  # Check for differences between treatE
  # Prepare output file name
  fileName <- paste("stats/monoVSco_nWorms_compTREATE_ANOVA+Dunnet",".txt", sep="")
  
  # Adjust factor levels of strains to test mono against co-culture abundances
  dataset$treatE <- factor(dataset$treatE, levels = c("anc", "mono", "bi"))
  
  # Run linear regression (LM with ANOVA)
  m <- lm(nWorms ~ treatE*repO, data=dataset)
  
  sink(fileName)
  print(as.data.frame(anova(m))) # ANOVA output
  print("")
  print(summary(glht(m, mcp(treatE="Tukey")))) # Post-hoc output: Tukey
  print("")
  print(summary(glht(m, mcp(treatE="Dunnett")))) # Post-doc output: Dunnett
  sink()
  
  
  return(p)  }

# Function to plot # washed worms of ancestors
plot_nWorms.strainANC <- function(dataset, cycle, cols, strains,plotH, plotW){
  # 1. Set cycle to exclude from data set
  ifelse(cycle == 4, cyc <- 10, ifelse(cycle == 10, cyc <- 4, 0))
  dataset <- dataset[dataset$cycle != cyc,]
  
  # Leave 20% space on top of max. data point for significance stars later
  top <- max(dataset$nWorms)*1.05
  
  dataset$strain <- factor(dataset$strain, levels = c("MYb11", "MYb71", "co-culture"))
  
  # 2. Print plot
  p <- ggplot(dataset, aes(x = strain, y = nWorms))+
    geom_boxplot(aes(fill = strain), size = 0.7) +
    geom_jitter(size = 2, alpha = 0.5, width = 0.1)+
    geom_point(aes(x=1, y=top), alpha=0)+
    xlab("Bacteria") + ylab("# worms/150?l")+
    scale_fill_manual(values = cols[c(1,3,2)]) +
    scale_color_manual(values = cols[c(1,3,2)]) +
    my_theme(dataset)+
    theme(axis.text.x = element_text(hjust = 1, angle = 45))
  
  # 3. Save plot to file
  # Check presence of directory
  absPlots <- paste(currentDirectory, "/plots", sep="")
  ifelse(dir.exists(absPlots) == F, 
         dir.create(absPlots), "Directory for absolute data already exists.")
  
  # Define file name
  fileExtensions <- c(".svg", ".png")
  for(i in 1:length(fileExtensions)){
    fileName <- paste(absPlots, "/", assayName, "_nWorms~strain_ANC_", strains, fileExtensions[i], sep = "")
    
    # Save plot
    ggsave(fileName, p, height=plotH, width=plotW,dpi=300) }
  
  # Check for differences between either monoculture and the co-culture within each treatE
  # Prepare output file name
  fileName <- paste("stats/monoVSco_nWorms_compANC_ANOVA+Dunnet",".txt", sep="")
  
  # Adjust factor levels of strains to test mono against co-culture abundances
  dataset$strain <- factor(dataset$strain, levels = c("co-culture", "MYb11", "MYb71"))
  
  # Run linear regression (LM with ANOVA)
  m <- lm(nWorms ~ strain*repO, data=dataset)
  
  sink(fileName)
  print(as.data.frame(anova(m))) # ANOVA output
  print("")
  print(summary(glht(m, mcp(strain="Tukey")))) # Post-hoc output: Tukey
  print("")
  print(summary(glht(m, mcp(strain="Dunnett")))) # Post-doc output: Dunnett
  sink()
  
  return(p)  }

# Subset data from a specific evol. cycle
subsetCycle <- function(data, cycle){
  ifelse(cycle == 10,
         cyc <- 4, 
         cyc <- 10)
  
  return(data[data$cycle != cyc, ])
}

# Check parametric assumptions in prep. for regressions
checkParametricAssumptions <- function(data, dims){
  # Subset data
  d <- data
  pheno <- d[,dims]
  d <- subset(d, select = c(strain, treatE, cycle, repO))
  d <- cbind(d, pheno)
  
  # Test for normality
  normalTest <- shapiro.test(d$pheno)$p.value
  
  # Test for homogeneity of variance
  homoTest <- levene.test(d$pheno, group = data$treatE)$p.value
  
  return(ifelse(normalTest > 0.05 & homoTest > 0.05, T, F)) }

# Run linear and generalized linear models
runRegression <- function(data, parametric, strains, dims){
  # Subset data
  d <- data
  pheno <- d[,dims]
  d <- subset(d, select = c(strain, treatE, cycle, repO))
  d <- cbind(d, pheno)
  
  # Only consider comparison of evolved changes
  cyc <- unique(data$cycle)[unique(data$cycle) != 0]
  
  # Prepare output file name
  fileNameTukey <- paste("stats out/", assayName,"/", dims, "_Tukey_C",cyc,"_",strains,".txt", sep = "")
  
  # Choose model depending on data structure
  ifelse(parametric == T, 
         # Run linear regression (LM with ANOVA)
         {
           # Check LM & save to file
           ifelse(strains == "MYb", 
                  m <- lm(pheno ~ treatE+strain+repO, data = d),
                  m <- lm(pheno ~ treatE+repO, data = d))
           
           fileName <- paste("stats out/", assayName,"/", dims, "_ANOVA_C", cyc,"_",strains,".txt", sep = "")
           
           write.csv(as.data.frame(anova(m)), fileName, sep = ",", row.names = T)
           
           # Tukey post-hoc comparison of treatE & save to file
           marginal <- emmeans(m, ~ treatE)
           write.csv(pairs(marginal), fileNameTukey, row.names = F)
           
         }, 
         # Run GLM with Gamma distribution
         {
           # Check GLM with Gamma family
           ifelse(strains == "MYb", 
                  m <- glm(pheno ~ treatE+strain+repO, family = "Gamma", data = data),
                  m <- glm(pheno ~ treatE+repO, family = "Gamma", data = data))
           
           fileName <- paste("stats out/", assayName,"/", dims, "_GLM_C", cyc,"_",strains,".txt", sep = "")
           
           sink(fileName)
           print(summary(m))
           sink()
           
           # Tukey post-hoc comparison of treatE & save to file
           marginal <-emmeans(m, ~ treatE)
           write.csv(pairs(marginal), fileNameTukey, row.names = F)
           
         })
  
  return("Regression analysis done and saved.")
}

#### = 3. Import data = ####
dataMeans <- read.csv("./data/MYb11-MYb71_sizesMY316_absData.NO.CIRCULARITY.csv", header = T)

#### = 4. Plot data = ####
wormDimensions <- c("avg.area","avg.minor.major" )#,"avg.major", "avg.minor")
labelsX <- c(expression("Area (mean mm"^"2"*")"),
             "Minor/Major axis length (mean in mm)")
             #"Major axis length (mean in mm)",
             #"Minor axis length (mean in mm)",              )

c <- 2 # focus on evolved cycle 10 material and ancestor
  for (d in 1:length(wormDimensions)){ # worm dimensions measured
    # main figure
    plot_wormDimension.treatE(dataMeans[dataMeans$strain == "co-culture",], cycs[c], colsMYb, "co-culture", wormDimensions[d],labelsX[d], 3.2, 2.5)
    plot_nWorms.treatE(dataMeans[dataMeans$strain == "co-culture",], cycs[c], colsMYb, "co-culture", 3.2, 2.5)
    
    plot_wormDimension.strainANC(dataMeans[dataMeans$treatE == "anc",],
                              cycs[c], colsMYb, "ANC_comparisons", wormDimensions[d],labelsX[d], 3.2, 2.5)
    plot_nWorms.strainANC(dataMeans[dataMeans$treatE == "anc",],
                       cycs[c], colsMYb, "ANC_comparisons", 3.2, 2.5)
    
    # supplement
    plot_wormDimension.treatE(dataMeans,cycs[c], colsMYb, "MYb", wormDimensions[d],labelsX[d], 3, 7)
    plot_nWorms.treatE(dataMeans, cycs[c], colsMYb, "MYb", 3, 7) }#}


#### = Clean up = ####
rm(dataMeans,  assayName, c, colsMYb, colsTreat,currentDirectory, cycs, d, labelsX, wormDimensions, 
   my_theme, plot_nWorms.treatE, plot_wormDimension.treatE, assay, 
   plot_wormDimension.strainANC, plot_nWorms.strainANC, subsetCycle, runRegression, 
   checkParametricAssumptions)

   