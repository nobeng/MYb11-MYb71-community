########################################
#            MYb11 - MYb71             #
#     competitions (colonization,      # 
#   short-term persistence)            #
#     Anna Czerwinski/ Nancy Obeng     #
########################################

#### Set dependencies, prepare environment ####
# Set working directory to source location
currentDirectory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(currentDirectory)

library(plyr); library(stringr); library(ggplot2); library(dplyr); library(multcomp)

source("my_theme.R")
colsMYb <- c("#B9006A","red4", "#139150",  "#139151")

# Proper axis labels
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format-*
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)}

#### Open & curate data ####
data <- read.csv("data/MYb11-MYb71_anc-evo-props_Col_STP_absData.csv.csv", header =T)

# order factor levels
data$assay <- factor(data$assay, levels=c("col.l4", "STP"))
data$strain <- factor(data$strain, levels = c("MYb11","wspE", "MYb71"))
data$treatP <- factor(data$treatP, levels = c("10_90", "50_50","90_10","99_1"))

# Add time points in hours to data frame
Assay <- c("col.l4", "STP")
names <- c("Colonization", "Short-term persistence")

data$names <- rep(1, length(data$assay))
for(n in 1:length(Assay)){data[data$assay==Assay[n], 'names'] <- names[n]  }
data$names <- factor(data$names, levels = c("Colonization","Short-term persistence"))

# Exclude negative entries
data <- data[data$cfu >= 0, ]
# Where there is only one strain, set other strain to log10=0 (i.e. single place holder cell)
data[data$cfu< 1, 'cfu'] <- 1
data$cfu <-as.numeric(data$cfu)

# === Focus on adaptive advantage in colonization and short-term persistence ====

# Focus on 1:1 competitions
#data_focus <- data[data$treatP == "50_50", ]
data_focus <- data

# Calculate proportion of MYb11/wspE in mix

# 1. Get total load per competition
data_focus_load <- data_focus %>%
  group_by(assay, cycle, treatP, Rep, names) %>%
  summarise(load_cfu = sum(cfu))

data_focus <- merge(data_focus, data_focus_load)

data_focus$proportion <- data_focus$cfu / data_focus$load_cfu
data_focus <- data_focus[data_focus$strain != "MYb71", ]

# Plot proportions
top <- max(data_focus$proportion)*2 # Leave 20% space on top of max. data point for significance stars later

# Main figure 3
ggplot(data_focus[data_focus$treatP =="50_50", ], aes(x=strain, y=proportion))+
  geom_abline(slope=0, intercept=10, color="red")+
  geom_boxplot(aes(fill = strain),size = 0.3)+
  geom_point(aes(fill=strain),position = position_dodge(width=0.7), size=1.3, alpha = 0.5)+
  geom_abline(slope=0, intercept=0.5, lty=2)+
  facet_grid(~assay)+
  labs(x = "MYb11 isolate", y = "Proportion in worm")+ 
  scale_fill_manual(values = colsMYb)+
  scale_color_manual(values = colsMYb)+
  my_theme()

 # Save plot to file
fileExtensions <- c(".png", ".svg")

# Statistical comparison
d <- data_focus[data_focus$treatP =="50_50" & data_focus$assay == "col.l4", ]
print(wilcox.test(d[d$strain == "MYb11", 'proportion'],d[d$strain == "wspE", 'proportion'], alternative = "two.sided", paired= TRUE))

d <- data_focus[data_focus$treatP =="50_50" & data_focus$assay == "STP", ]
print(wilcox.test(d[d$strain == "MYb11", 'proportion'],d[d$strain == "wspE", 'proportion'], alternative = "two.sided",paired= TRUE))

for(i in 1:length(fileExtensions)){
  ggsave(paste("plots","/props~MYb11_anc-vs-evo",fileExtensions[i], sep=""),
         dpi = 300, height = 3.2, width = 6)  }
