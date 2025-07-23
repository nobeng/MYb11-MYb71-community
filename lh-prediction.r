library(data.table)
library(stringr)
library(ggplot2)
library(ggrepel)

csr.dt <- fread("data/csr-cub.csv")
dat <- fread("data/lht.csv.gz")
taxonomy.dt <- fread("data/taxonomy.csv")
colors.tax <- fread("data/colors_taxa.csv")

taxonomy.dt <- merge(taxonomy.dt, colors.tax, by.x="genus", by.y="group", all.x=T)
taxonomy.dt[is.na(color), color:=setdiff(BacArena::colpal3,color)[1:.N]]
dat[,group:=taxonomy.dt$genus[match(org, taxonomy.dt$id)]]
csr.dt[,group:=taxonomy.dt$genus[match(org, taxonomy.dt$id)]]

# Subfigure 1 (max growth)
ggplot(dat, aes(x=log(2)/grodon_d, y=0)) + 
  labs(color="") + geom_vline(xintercept=log(2)/5, linetype="dashed", color = "red") + 
  geom_point(data=dat[!org%in%c("MYb11","MYb71")], size=5, color="gray", alpha=0.5) + 
  geom_point(data=dat[org%in%c("MYb11","MYb71")], size=6, aes(color=group)) + 
  xlab("Maximal growth rate (codon usage)") + ylab("") +
  geom_label_repel(data=dat[org%in%c("MYb11","MYb71")], size=6, nudge_y = 0.1, aes(label=org, color=group), max.overlaps = 100, show.legend = FALSE) + 
  scale_color_manual(values=setNames(taxonomy.dt$color,taxonomy.dt$genus)) + 
  labs(color="") + 
  scale_y_continuous(limits = c(-0.1,0.15)) + 
  scale_x_continuous(limits = c(0,1.2)) +
  theme_minimal(base_size = 14) +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position="none")
ggsave("plots/cembio-max-growth.pdf", height=2, width=5)


# Subfigure 2 (CSR triangle)
coord.lst <- list(c=c(-1,0), s=c(0,1), r=c(1,0))
csr.coord <- lapply(csr.dt$csr, function(csr) Reduce("+",coord.lst[unlist(str_split(csr,","))]))
csr.dt <- cbind(csr.dt, data.table(x=sapply(csr.coord, function(x) x[1]), y=sapply(csr.coord, function(x) x[2])))
csr.dt[grepl(",",csr), `:=`(x=x/2, y=y/2)]
ggplot(csr.dt) + 
  geom_jitter(data=csr.dt[!org%in%c("MYb11","MYb71")], width=0.05,height=0.05,aes(x=x,y=y), size=2, alpha=0.5) + 
  geom_point(data=csr.dt[org%in%c("MYb11","MYb71")], aes(x=x,y=y), size=3) + 
  geom_polygon(data=data.frame(x=c(-1,1, 0), y=c(0,0,1)), aes(x=x,y=y), alpha=0.1, fill="blue") + 
  annotate(geom="text", x=c(-1.1,0,1.1),y=c(0,1.05,0), label=c("C","S","R"), color="red", size=10) + 
  geom_label_repel(data=csr.dt[org%in%c("MYb11","MYb71")], size=6, nudge_y = -0.1, nudge_x = 0.4, aes(x=x,y=y,label=org, color=group), max.overlaps = 100, show.legend = FALSE) + 
  scale_color_manual(values=setNames(taxonomy.dt$color,taxonomy.dt$genus)) + 
  xlab("") + ylab("") + theme_minimal(base_size=14) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")
ggsave("plots/cembio-csr.pdf", height=3.5, width=3.7)
