# setup ####
library(phyloseq)
library(tidyverse)

'%ni%' <- Negate('%in%')

# load data ####
ps <- readRDS("./Output/ps_cleaned_w_tree.RDS")

ps@sam_data$sample_type %>% table
# relax parameters in cleaning step to avoid losing so many samples???

# for each tree where there are both epi and endo...
# calculate proportion of endo taxa also found in epi sample

prop_nested <- c()
for(i in ps@sam_data$tree){
  x <- ps %>% subset_samples(tree == i)
  x <- x %>% subset_taxa(taxa_sums(x) > 0)
  if(nsamples(x) < 2){prop_nested[i] <- NA;next}
  endo <- x %>% subset_samples(sample_type == "endophyte")
  endo <- endo %>% subset_taxa(taxa_sums(endo) > 0)
  epi <- x %>% subset_samples(sample_type == "epiphyte")
  epi <- epi %>% subset_taxa(taxa_sums(epi) > 0)
  prop_nested[i] <- sum(taxa_names(endo) %in% taxa_names(epi)) / ntaxa(endo)
}
dev.off()
prop_nested %>% plot

nested_df <- 
microbiome::meta(ps) %>% 
  full_join(data.frame(nestedness=prop_nested,tree=names(prop_nested))) %>% 
  select(-c(sample,sample_type,fwd_filepath,rev_filepath)) 
nested_df %>% 
  ggplot(aes(x=distance_from_edge,y=nestedness)) +
  geom_point() +
  theme_minimal() +
  labs(x="Distance from edge (m)",y="Proportion of endophyte taxa\nnested within epipythe community") +
  theme(axis.title = element_text(face='bold',size=14),
        axis.text = element_text(face='bold',size=12)) +
  geom_smooth(method="lm",se=FALSE,color='gray')
ggsave("./Output/Figs/community_nestedness_vs_edge.png",height = 6,width = 6,dpi=300)
nested_df %>% 
  glm(data=.,formula=nestedness~distance_from_edge) %>% 
  summary()
prop_nested %>% summary
