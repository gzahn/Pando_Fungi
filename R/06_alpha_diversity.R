# Alpha diversity

# SETUP ####
library(phyloseq)
library(vegan)
library(tidyverse)
library(broom)
library(ggmap)
library(zahntools)
library(patchwork)
options(scipen = 999)
sampletypecolors <- c("darkblue","#fca311")
# functions ####

add_alphadiv_measures <- function(physeq){
  # build species-level physeq
  ps_species <- 
    tax_glom(physeq,"Species")
  
  alpha_asv <- estimate_richness(physeq,measures = c("Observed","Shannon","Simpson"))
  alpha_spp <- estimate_richness(ps_species,measures = c("Observed","Shannon","Simpson"))
  
  # add to sample metadata
  physeq@sam_data$asv_richness <- alpha_asv$Observed
  physeq@sam_data$asv_shannon <- alpha_asv$Shannon
  physeq@sam_data$asv_simpson <- alpha_asv$Observed
  
  # build species-level physeq
  ps_species <- 
    tax_glom(physeq,"Species")
  
  alpha_asv <- estimate_richness(physeq,measures = c("Observed","Shannon","Simpson"))
  alpha_spp <- estimate_richness(ps_species,measures = c("Observed","Shannon","Simpson"))
  
  # add to sample metadata
  physeq@sam_data$asv_richness <- alpha_asv$Observed
  physeq@sam_data$asv_shannon <- alpha_asv$Shannon
  physeq@sam_data$asv_simpson <- alpha_asv$Observed
  
  physeq@sam_data$spp_richness <- alpha_spp$Observed
  physeq@sam_data$spp_shannon <- alpha_spp$Shannon
  physeq@sam_data$spp_simpson <- alpha_spp$Observed
  
  return(physeq)
}

ra <- function(x){x/sum(x)}

'%ni%' <- Negate('%in%')

## data ####
# asv-level physeq
ps <- readRDS("./Output/ps_cleaned_w_tree.RDS")
# species-level physeq
ps_species <- 
  tax_glom(ps,"Species")


# richness
# shannon div
# at different taxonomic levels
# difference between endo/epi?

# CALCULATE ALPHA-DIV ####

# add to sample metadata
ps <- add_alphadiv_measures(ps)

# extract sample data
sam <- ps@sam_data %>% as("data.frame")

# check correlation of vars
glm(data=sam,formula = distance_from_edge ~ sample_type) %>% 
  summary

m1 <- glm(data=sam,formula = asv_richness ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "asv_richness")
m2 <- glm(data=sam,formula = asv_shannon ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "asv_shannon")
m3 <- glm(data=sam,formula = asv_simpson ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "asv_simpson")
m4 <- glm(data=sam,formula = spp_richness ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "spp_richness")
m5 <- glm(data=sam,formula = spp_shannon ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "spp_shannon")
m6 <- glm(data=sam,formula = spp_simpson ~ sample_type * distance_from_edge) %>% 
  tidy() %>% mutate(outcome = "spp_simpson")

mod_tab <- list(m1,m2,m3,m4,m5,m6) %>% 
  reduce(full_join) %>% 
  filter(term != "(Intercept)") %>% 
  mutate(term = term %>% str_remove("sample_type"))
saveRDS(mod_tab,"./Output/alpha_div_models.RDS")


mod_tab %>% 
  kableExtra::kable() %>% 
  kableExtra::kable_classic() %>% 
  kableExtra::column_spec(5,bold = if_else(mod_tab$p.value < 0.05,TRUE,FALSE))

sam_long <- 
sam %>% 
  mutate(across(all_of(c("asv_richness","asv_shannon","asv_simpson",
                         "spp_richness","spp_shannon","spp_simpson")),scale)) %>% 
  pivot_longer(c(asv_richness,asv_shannon,asv_simpson,
                 spp_richness,spp_shannon,spp_simpson))


  
# plot alpha div
ps_melted <- 
  psmelt(ps)
ps_melted %>% names

df <- 
ps_melted %>% 
  pivot_longer(c(spp_richness,spp_shannon),names_to = "alpha_div",values_to = "measure") %>% 
  dplyr::select(alpha_div,measure,Sample,sample_type) %>% 
  unique.data.frame()

df$Sample <- factor(df$Sample,levels= df %>% 
                      filter(alpha_div == "spp_richness") %>% 
                      arrange(measure) %>% 
                      pluck("Sample"))

alpha_plot <- 
df %>% 
  mutate(alpha_div = case_when(alpha_div == "spp_richness" ~ "Richness",
                               alpha_div == "spp_shannon" ~ "Shannon"),
         sample_type = sample_type %>% str_to_sentence) %>% 
  ggplot(aes(x=Sample,y=measure,color=sample_type)) +
  geom_point(size=3,alpha=.75) +
  facet_wrap(~alpha_div,scales='free') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.title = element_text(face='bold',size=16),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=16),
        strip.text = element_text(face='bold',size=16),
        strip.background = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(color="Sample type",x="Samples",y="Diversity value") +
  scale_color_manual(values=sampletypecolors)
saveRDS(alpha_plot,"./Output/figs/alpha_plot.RDS")
ggsave("./Output/figs/alpha_plot.png",dpi=300,height = 6,width = 6)

ps_species@tax_table[,1] <- ps_species@tax_table[,1] %>% str_remove(".__")
ps_species@tax_table[,2] <- ps_species@tax_table[,2] %>% str_remove(".__")
ps_species@tax_table[,3] <- ps_species@tax_table[,3] %>% str_remove(".__")
ps_species@tax_table[,4] <- ps_species@tax_table[,4] %>% str_remove(".__")
ps_species@tax_table[,5] <- ps_species@tax_table[,5] %>% str_remove(".__")
ps_species@tax_table[,6] <- ps_species@tax_table[,6] %>% str_remove(".__")
ps_species@tax_table[,7] <- ps_species@tax_table[,7] %>% str_remove(".__")

relabund <- 
ps_species %>% 
  transform_sample_counts(ra)

class_ra <- 
ps_species %>% 
  tax_glom("Class") %>% 
  transform_sample_counts(ra)

class_mat <- 
class_ra@otu_table %>%
  as("matrix")
class_mat[is.nan(class_mat)] <- 0
# class_ra@otu_table %>% View
class_ra@otu_table[,which(class_mat %>% colMeans() >= 0.01)] %>% colnames()

lowabund_asvs <- (class_ra@tax_table %>% rownames()) %ni% (class_ra@otu_table[,which(class_mat %>% colMeans() >= 0.01)] %>% colnames())

class_ra@tax_table[lowabund_asvs,"Class"] <- "Other"


class_ra@tax_table[,3] <- factor(class_ra@tax_table[,3],levels = c("Leotiomycetes","Dothideomycetes","Agaricomycetes" ,"Tremellomycetes","Other"))

# reorder samples by dist from edge

class_ra@sam_data$distance_from_edge <- 
  factor(class_ra@sam_data$distance_from_edge,
         levels = class_ra@sam_data %>% 
           as("data.frame") %>% 
           arrange(desc(distance_from_edge)) %>% 
           pluck("seq_coast_tube_id"))


# plot taxonomic breakdown
class_ra %>% 
  plot_bar2(fill="Class") +
  # scale_x_discrete(limits = class_ra@sam_data %>% 
  #                    as("data.frame") %>% 
  #                    arrange(desc(distance_from_edge)) %>% 
  #                    pluck("seq_coast_tube_id")) +
  facet_wrap(~sample_type,scales = 'free') +
  scale_fill_manual(values = c("#A00102","#6098a3","#c97724","#297a18","gray6"),
                    labels = c("Leotiomycetes","Dothideomycetes","Agaricomycetes" ,"Tremellomycetes","Other")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face='bold',size=20),
        legend.title = element_text(face='bold',size=20),
        legend.text = element_text(face='bold',size=18),
        axis.title = element_text(face='bold',size=20)) +
  labs(y="Relative abundance")



ggsave("./Output/figs/Figure_2.png",dpi=500,width = 12,height = 8)
ps_species@sam_data$sample_type <- factor(ps_species@sam_data$sample_type, levels  = c("epiphyte","endophyte"))



# make fig 1 - map

# Load google maps API key from .Renviron 
ggmap::register_google(key = Sys.getenv("APIKEY")) # Key kept private

# customization
mapstyle <- rjson::fromJSON(file = "./R/mapstyle.json") %>% # from JSON file exported from snazzymaps.com
  googlemap_json_to_string(.)

# data
path <- "./Data/pando_sample_metadata_filled.csv"
dat <- read_csv(path)

latlon <- 
  dat %>% 
  select(tree,lon,lat) %>% 
  filter(!is.na(lon) & !is.na(lat))

# BUILD MAP ####
area <- 
  ggmap::get_googlemap(center = c(lon = mean(latlon$lon), 
                                  lat = mean(latlon$lat)),
                       zoom = 16,
                       scale = 2,
                       style=mapstyle)

sam
map_epi <- 
  ggmap::ggmap(area) +
  geom_point(data=sam %>% filter(sample_type == "epiphyte"),aes(x=lon,y=lat,text=tree,fill=asv_shannon),
             size=4,alpha=.7,shape=24) +
  scale_fill_viridis_c(option = 'rocket',begin = .2) +
  labs(fill="Shannon\ndiversity",x="Longitude",y="Latitude",title = "Epiphyte") +
  theme(legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=16),
        plot.title = element_text(face='bold',size=20,hjust = .5))
map_epi

map_endo <- 
  ggmap::ggmap(area) +
  geom_point(data=sam %>% filter(sample_type == "endophyte"),aes(x=lon,y=lat,text=tree,fill=asv_shannon),
             size=4,alpha=.7,shape=24) +
  scale_fill_viridis_c(option = 'rocket',begin = .2,end = .6) +
  labs(fill="Shannon\ndiversity",x="Longitude",y="Latitude",title = "Endophyte") +
  theme(legend.position = 'none') +
  theme(legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=14),
        axis.title = element_text(face='bold',size=16),
        plot.title = element_text(face='bold',size=20,hjust = .5))


map_endo + map_epi +
  patchwork::plot_layout(guides='collect')
ggsave("./Output/figs/Figure_1.png",dpi=500,width = 12,height = 6)

ggmap::ggmap(area) +
  geom_point(data=latlon,aes(x=lon,y=lat,text=tree),color='red',size=2,alpha=1,shape=4)














library(corncob)
da_analysis <- differentialTest(formula = ~ sample_type, #abundance
                 phi.formula = ~ 1, #dispersion
                 formula_null = ~ 1, #mean
                 phi.formula_null = ~ 1,
                 test = "Wald", boot = FALSE,
                 data = ps_species,
                 fdr_cutoff = 0.05,
                 full_output = TRUE)

plot(da_analysis)



# da_analysis$significant_taxa[c(1,4)]
# bbd1 <- bbdml(formula = GCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACATTGCGCCCTTTGGCATTCCGAAGGGCATGCCTGTTCGAGCGTCATTACAACCACTCAAGCACTCGCTTGGCCTTGGGGCACCCGGCGTCGGGGCCCTCAAAACCAGCGGCGGTGCTCGTCAGCTCTACGCGTAGTAATACTCCTCGCGTCTGGCCCTGGCGAGCCACCGGCCGGCAACCCCCCACACTTCTCAGGTTGACCTCGGATCAGGTAGGGATA ~ sample_type,
#       phi.formula = ~ 1, #dispersion
#       formula_null = ~ 1, #mean
#       phi.formula_null = ~ 1,
#       test = "Wald", boot = FALSE,
#       data = ps_species,
#       fdr_cutoff = 0.05)
# plot(bbd1,color = "sample_type")
# 
# bbd2 <- bbdml(formula = GCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCTTTGAACGCACCTTGCGCTCCTTGGTATTCCGAGGAGCATGCCCGTTTGAGTGTCTTGATATCATCAAACGCCCCTGCCTTTTCTGGTGGGTCGGCGATTTGGACTTGGACGTGCTGCTGGGCGCTCTGCGCTCGGCTCGTCTTGAATGCATTAGCTGGATCTCCCTCCGGCCTTGGACCTTATGGCGTGATAAGATCACCCGTCGAGCAGTCCTTGGTCAACGTTGGGTT ~ sample_type,
#               phi.formula = ~ 1, #dispersion
#               formula_null = ~ 1, #mean
#               phi.formula_null = ~ 1,
#               test = "Wald", boot = FALSE,
#               data = ps_species,
#               fdr_cutoff = 0.05)
# plot(bbd2,color = "sample_type")

library(fungaltraits); packageVersion("fungaltraits")
# download traits metadata
traits_meta <- read_csv("https://github.com/traitecoevo/fungaltraits/releases/download/v0.0.3/funtothefun.csv")

# download FungalTraits database
traits_db <- fungaltraits::fungal_traits()
# match taxa at genus level
genera <- ps_species@tax_table[,6] %>% str_remove("^g__")
species <- ps_species@tax_table[,7] %>% str_remove("^s__")
fungal_traits <- 
  data.frame(Genus=genera) %>% 
  mutate(species=paste(Genus,species,sep="_")) %>% 
  left_join(traits_db,by=c("species","Genus"),multiple='all')

fungal_traits[which(fungal_traits$Genus == "Cytospora"),]

fungal_traits[!is.na(fungal_traits$guild_fg),]
