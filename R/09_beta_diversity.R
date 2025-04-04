# SETUP ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(patchwork)
sampletypecolors <- c("darkblue","#fca311")
set.seed(666)


# data ####
ps <- readRDS("./Output/ps_cleaned_w_tree.RDS")
ps_ra <- ps %>% 
  transform_sample_counts(function(x){x/sum(x)})
ord <- ordinate(ps_ra,method = "NMDS",distance = 'unifrac')

ord_plot <- 
plot_ordination(ps,ord,color="sample_type") + 
  stat_ellipse() +
  theme_bw() + 
  scale_color_manual(values = sampletypecolors,labels=c("Endophyte","Epiphyte")) +
  theme(axis.title = element_text(face='bold',size=16),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        legend.position = 'none') +
  labs(color="Sample type")
saveRDS(ord_plot,"./Output/figs/ordination_plot_sampletype.RDS")  
ggsave("./Output/figs/ordination_plot_sampletype.png",dpi=300,height = 6,width = 6)
# plot_ordination(ps,ord,color="distance_from_edge")

# Permutational ANOVA
mod <- adonis2(ps_ra@otu_table ~ ps_ra@sam_data$sample_type + ps_ra@sam_data$distance_from_edge)
mod
# OUTCOME ~ PREDICTORS




alpha_plot <- readRDS("./Output/figs/alpha_plot.RDS")

diversity_plot <- 
alpha_plot / ord_plot +
  plot_layout(guides='collect') +
  plot_annotation(tag_levels = "A")
saveRDS(diversity_plot,"./Output/figs/Figure_3.RDS")
ggsave("./Output/figs/Figure_3.png",dpi=300,height = 6,width = 8)
