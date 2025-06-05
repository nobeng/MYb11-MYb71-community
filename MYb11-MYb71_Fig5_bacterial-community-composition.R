### ===    Colonization of evolving worms   === ###
###                  Nancy Obeng                ###
###################################################
# Clean working directory
rm(list=ls())

# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

library(plyr); library(ggplot2); library(gridExtra); library(dplyr); library(ggh4x)
  
# Source files
source("my_theme.R")

### == Load and organize data === ####

# Import data
evoAll <- read.csv("data/MYb11-MYb71_colonization-worms_proportions_evo-anc_absData.csv", header = T)

# Calculate species proportions
evoAll$propMYb11 <- evoAll$cfu_MYb11 / (evoAll$cfu_MYb11 + evoAll$cfu_MYb71)

### === Plotting == ###
plotH <- 3.2
plotW <- 7

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
} # proper y axis

evoAll <- evoAll %>%
          group_by(popID, treatment, cycle) %>%
          summarize(propMYb11 = mean(propMYb11))

evoAll$treatment <- as.factor(evoAll$treatment)

# # Leave 20% space on top of max. data point for significance stars later
top <- 1.1

ggplot(data = evoAll[evoAll$cycle %in% c(1, 10), ],aes(x = as.factor(cycle), y = propMYb11, fill = as.factor(cycle))) +
  geom_boxplot(size = 0.7, outlier.shape = NA) +
  geom_jitter(size = 2.5, alpha = 0.5, width = 0.1)+
  geom_point(aes(x=1, y=top), alpha=0)+
  facet_grid(~treatment) +
  xlab("") + ylab("MYb11 proportion\nin worm population")+
  scale_fill_manual(values = c("lightgrey", "grey30")) +
  my_theme(dataset)+
  force_panelsizes(rows = unit(2.15, "in"),
                   cols = unit(1.3, "in"))

# Save to files
fileExtensions <- c(".svg", ".png")
for(i in 1:length(fileExtensions)){
  fileName <- paste("plots/colonization_wormPop_co-culture", fileExtensions[i], sep = "")
  
  # Save plot
  ggsave(fileName, height=plotH, width=plotW, dpi=300) }


### === Statistics == ###
wilcox.test(pull(evoAll[evoAll$treatment == "biphasic" & evoAll$cycle == 1, 'propMYb11']),
            pull(evoAll[evoAll$treatment == "biphasic" & evoAll$cycle == 10, "propMYb11"]))

wilcox.test(pull(evoAll[evoAll$treatment == "monophasic" & evoAll$cycle == 1, 'propMYb11']),
            pull(evoAll[evoAll$treatment == "monophasic" & evoAll$cycle == 10, "propMYb11"]))
  