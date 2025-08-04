#BiocManager::install("coRdon")
#devtools::install_github("jlw-ecoevo/gRodon2")
library("gRodon")
library("Biostrings")

# MYb71
genome <- readDNAStringSet("data/MYb71.ffn.gz") # bakta
highly_expressed <- grepl("ribosomal protein",names(genome),ignore.case = T)
predictGrowth(genome, highly_expressed)

# MYb11
genome <- readDNAStringSet("data/MYb11.ffn.gz") # bakta
highly_expressed <- grepl("ribosomal protein",names(genome),ignore.case = T)
predictGrowth(genome, highly_expressed)
