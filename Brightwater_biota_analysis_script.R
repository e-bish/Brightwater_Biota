# load libraries
library(tidyverse)
library(janitor)
library(vegan)
library(indicspecies)
library(viridis)
library(RColorBrewer)

# load data
import_100ft <- read_csv("raw_data/import_100ft.csv")
import_300ft <- read_csv("raw_data/import_300ft.csv")
import_600ft <- read_csv("raw_data/import_600ft.csv")
import_600ft_ref <- read_csv("raw_data/import_600ft_ref.csv")
species_metadata <- read_csv("raw_data/species_metadata.csv")

# set seed for permutations
set.seed(2025)

##### prepare dataframes for analysis ####

#combine dataframes
combined_import <- bind_rows(import_100ft, import_300ft, import_600ft, import_600ft_ref) %>% #combine imported data into a single dataframe
  select(!contains("..."))
  
#extract QC rows 
QC_data <- combined_import %>% #these rows are repeats of data we already have so we can set them aside for now
  filter(grepl("QC", Notes))

#tidy data
combined_tidy <- combined_import %>% 
  filter(!is.na(Year)) %>% #remove empty rows
  anti_join(QC_data) %>% #remove QC rows
  filter(!Year == 2017) %>% #remove for this analysis because there wasn't enough replication
  rename(Depth_ft = Depth) %>% 
  mutate(Replicate = as.factor(Replicate), 
         Year = factor(Year, labels = c("2", "10"))) %>% #recode as deployment intervals
  mutate(Depth = case_when(Depth_ft == 100 ~ 30,
                           Depth_ft == 300 ~ 90, 
                           TRUE ~ 200),
         Site = ifelse(grepl("R", Grid.ID), paste0(Depth, "_ref"), Depth), .after = Depth_ft) %>% 
  mutate(Site = factor(Site, levels = c("30", "90", "200", "200_ref"))) %>% 
  full_join(species_metadata, by = "SpeciesID") %>% 
  filter(!is.na(Year))

#create motile and nonmotile dataframes for separate analysis
motile_tidy <- combined_tidy %>% 
  filter(Mobility == "Motile") 

nonmotile_tidy <- combined_tidy %>% 
  filter(Mobility == "Nonmotile") 

remove_algae <- combined_tidy %>% 
  filter(is.na(Recode_for_MS)) %>% 
  mutate(Grid.Area.cm = Grid.Area.cm - Org.Area.cm) %>% 
  select(PhotoID, Grid.Area.cm)

nonmotile_tidy<- nonmotile_tidy %>% 
  rows_update(remove_algae, by = "PhotoID") %>% 
  filter(!is.na(Recode_for_MS))

#### format cover and count data for tables ####

rep_grid_area_sum <- nonmotile_tidy %>% 
  group_by(Year, Site, Replicate) %>% 
  summarize(total_grid_area = sum(Grid.Area.cm)) 

nonmotile_table <- nonmotile_tidy %>% 
  group_by(Site, Year, Replicate, Phylum, Recode_for_MS, MS_Species_Group) %>% 
  summarize(group_area = sum(Org.Area.cm)) %>% 
  right_join(rep_grid_area_sum) %>% 
  mutate(group_prop = round(group_area/total_grid_area, 4)*100) %>% 
  group_by(Site, Year, Phylum, MS_Species_Group) %>% 
  summarize(mean = round(mean(group_prop), 2), sd = round(sd(group_prop), 2)) %>% 
  ungroup() %>% 
  replace_na(list(sd = 0)) %>% 
  mutate(mean_sd = paste0(mean, " (", sd, ")")) %>% 
  select(!c(mean, sd)) %>% 
  mutate(site_year = paste(Year, Site, sep = "_")) %>% 
  select(!Year:Site) %>% 
  pivot_wider(names_from = site_year, values_from = mean_sd) %>% 
  arrange(Phylum, MS_Species_Group) %>% 
  replace(is.na(.), "-")

# write_csv(nonmotile_table, file = "outputs/nonmotile_table.csv")

motile_table <- motile_tidy %>% 
  group_by(Site, Year, Replicate, Phylum, MS_Species_Group) %>% 
  summarize(Org_sum = sum(Species.Count)) %>% 
  group_by(Site, Year, Phylum, MS_Species_Group) %>% 
  summarize(mean = round(mean(Org_sum), 2), sd = round(sd(Org_sum), 2)) %>% 
  ungroup() %>% 
  mutate(mean_sd = paste0(mean, " (", sd, ")")) %>% 
  select(!c(mean, sd)) %>% 
  mutate(site_year = paste(Year, Site, sep = "_")) %>% 
  select(!Year:Site) %>% 
  pivot_wider(names_from = site_year, values_from = mean_sd) %>% 
  arrange(Phylum, MS_Species_Group) %>% 
  replace(is.na(.), "-")

# write_csv(motile_table, file = "outputs/motile_table.csv")

## summary plots based on % cover 
#calculate nonmotile percent live cover
perc_cover_df <- nonmotile_tidy %>% 
  group_by(Year, Depth, Site, Replicate) %>% 
  summarize(grid.area.sum = sum(Grid.Area.cm), #sum the total area across all grids on a given plate
            live.area.covered = sum(Org.Area.cm)) %>% 
  ungroup() %>% 
  mutate(perc.live.cover = round((live.area.covered/grid.area.sum)*100, 2)) #divide the area covered by the total area to get % cover

depth_colors <- c("#fde725","#4de3cc","#440154", "#fc8961")

#### Figure 3 ####
perc_cover_df %>%
  mutate(Site = factor(Site, labels = c("Shallow outfall", "Mid outfall", "Deep outfall", "Deep reference"))) %>% 
  ggplot(aes(x = Year, y = perc.live.cover, fill = Site, shape = Site)) +
  geom_point(size = 4, alpha = 0.85) +
  theme_classic() +
  labs(y = "Percent live cover",
       x = "Years deployed",
       fill = "Site") +
  scale_fill_manual(values = depth_colors) +
  scale_shape_manual(values = c(21, 23, 24, 25)) + 
  theme(text = element_text(size = 14), 
        strip.background = element_blank()) 

# ggsave("outputs/figure3.tiff", width = 8, height = 6, dpi = 300)

#### summary of nonmotile phyla ####

## summary plot for number of species
spp_richness_df <- nonmotile_tidy %>% 
  group_by(Year, Site, Replicate) %>% 
  distinct(Species_Group) %>% 
  summarize(group_richness = n()) %>% 
  ungroup()

spp_richness_df %>% 
  group_by(Year, Site) %>% 
  summarize(mean_sr = mean(group_richness),
            sd_sr = sd(group_richness)) %>% 
  mutate(sd_sr = ifelse(is.na(sd_sr), 0, sd_sr)) %>% 
  ggplot(aes(x = factor(Site, labels = c("Shallow outfall", 
                                         "Mid outfall",
                                         "Deep outfall",
                                         "Deep reference")), 
             y = mean_sr, fill = Year)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = mean_sr - sd_sr, ymax = mean_sr + sd_sr), 
                width=.1, position = position_dodge(width = .9)) +
  scale_fill_manual(values = c("#f3e0af", "#1e9b8a")) +
  theme_classic() +
  labs(y = "Average number of distinct 
       organism categories", 
       x = "Site", 
       fill = "Years\ndeployed")

# ggsave("outputs/figure4.tiff", width = 8, height = 6, dpi = 300)

#phylum proportion values 
area_sum <- nonmotile_tidy %>% 
  group_by(Year, Site) %>% 
  mutate(total_live_area = sum(Org.Area.cm)) %>% 
  dplyr::select(Year, Site, total_live_area) %>% 
  distinct()

#breakdown by organism category
group_perc_live_cover <- nonmotile_tidy %>% 
  group_by(Year, Site, Recode_for_MS) %>% 
  summarize(group_area = sum(Org.Area.cm)) %>% 
  right_join(area_sum) %>% 
  mutate(group_prop = round(group_area/total_live_area, 4)*100)

#### univariate stats ####
shapiro.test(x = perc_cover_df$perc.live.cover)

#test for differences across all sites
perc_cover_mod <- aov(perc.live.cover ~ Site + Year + Site*Year, data = perc_cover_df)
summary(perc_cover_mod)

#model diagnostics
aov.res <- residuals(perc_cover_mod)

hist(perc_cover_mod$residuals)
qqnorm(perc_cover_mod$residuals)
qqline(perc_cover_mod$residuals)

#### Figure 5 ####
## plot proportions of nonmotile taxa
phylum_colors <- c(
  "#44627a", "#6b8ba4", "#8fb9a8",
  "#d6c690", "#eea98d", "#fff4b3", "#e9e3da"
)

phylum_colors1 <- c(
  "#6b8ba4", "#5c3ea0",  "#2f77bd",
  "#3fb26b",  "#d6c690","#ffd64d", "#f4a259"
)

phylum_colors2 <- c( "#fc8961","#1e9b8a",
                     "#f3e0af",
                     "#4de3cc","#6b8ba4",
                     "#fde725","lightgrey")

nonmotile_tidy %>% 
  mutate(Group = ifelse(Phylum %in% c("Unknown", "Tunicata", "Chordata", "Brachiopoda"), "Other", Phylum)) %>% 
  ggplot(aes(x = Year, y = Org.Area.cm, fill = Group)) +
  geom_col(position = "fill") +
  theme_classic() + 
  facet_wrap(~factor(Site, 
                     labels = c("Shallow\noutfall", "Mid\noutfall", "Deep\noutfall", "Deep\nreference"))) + 
  # theme(text = element_text(size = 14)) +
  labs(y = "Proportion of total live area", x = "Years deployed", fill = "Phylum") +
  scale_fill_manual(values = phylum_colors2) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, vjust = 1),
        strip.background = element_blank()) 

# ggsave("outputs/figure5.tiff", width = 8, height = 6, dpi = 300)

#### nonmotile multivariate analysis ####

#convert the data to wide format
nonmotile_tidy_w <- nonmotile_tidy %>% 
  group_by(Year, Depth, Site, Replicate, Recode_for_MS) %>% 
  summarize(sum_area = sum(Org.Area.cm)) %>% 
  filter(!is.na(Recode_for_MS)) %>% 
  ungroup() %>% 
  arrange(Recode_for_MS) %>% 
  pivot_wider(names_from = Recode_for_MS, values_from = sum_area, values_fill = 0) %>% 
  clean_names() %>% 
  arrange(year, site, replicate)

#prep abundance and metadata matrices

abund <- nonmotile_tidy_w %>%
  select(!contains("unidentified")) %>% #not including unidentified orgs
  select(amphipod_tube:urticina_sp) 

meta <- nonmotile_tidy_w %>%
  select(year:replicate)

#conduct data transformations and/or relativizations 
#not removing rare species because we're looking for differences in diversity

#function for calculating the coefficient of variation
CV <- function(x) { 100 * sd(x) / mean(x) }

CV(x = rowSums(abund)) #low (<50), relativizing by row will not make a big difference
CV(x = colSums(abund)) #high (100+), will have an effect to relativize between species

#### nonmotile NMDS #### 
nonmotile.nmds <- metaMDS(abund, distance = "bray", 
                          autotransform = FALSE,
                          engine = "monoMDS",
                          k = 3,
                          weakties = TRUE,
                          model = "global",
                          maxit = 400,
                          try = 40,
                          trymax = 100, 
                          trace = FALSE)
#post-hoc test
stressplot(object = nonmotile.nmds, lwd = 5)
#stress values less than 0.2 but greater than 0.10 are decent but shouldn't be relied on for details

#extract points and bind to metadata
nmds_points <- data.frame(nonmotile.nmds$points)
nmds_points <- bind_cols(meta, nmds_points) %>% 
  mutate(depth = factor(depth, labels = c("30", "90", "200")),
         year = factor(year))

#### Figure 6 ####
#nmds plot
ggplot(data=nmds_points,
       aes(x=MDS1, y=MDS2,
           fill= factor(site,
                        labels = c("Shallow outfall", "Mid outfall", 
                                   "Deep outfall", "Deep reference")), 

           shape= year)) + 
  geom_point(color = "black", size = 4) +
  # stat_ellipse(aes(group = depth, color = depth), 
  #              linetype = "dashed", show.legend = FALSE) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  labs(fill = "Site", shape = "Year") + 
  scale_shape_manual(values = c(21,24,22)) +
  scale_fill_manual(values = depth_colors) + 
  scale_color_manual(values = depth_colors) +
  guides(fill=guide_legend(override.aes=list(color=c(depth_colors)))) +
  annotate("text", x = -1, y = 1.4, 
           label = paste("Stress = ", round(nonmotile.nmds$stress, 3))) +
  theme(text = element_text(size = 14))

# ggsave("outputs/figure6.tiff", width = 8, height = 6, dpi = 300)

#### Nonmotile multi stats ####
## test for differences between depths and years

adonis2(abund ~ depth:year, method = "bray",
        data = meta, permutations = 9999, by = "margin")
#significant interaction term
#don't have enough replication to do pairwise tests

permutest(betadisper(vegdist(abund, method = "bray"), group = meta$depth, type = "median"), pairwise = TRUE)
#not significant

permutest(betadisper(vegdist(abund, method = "bray"), group = meta$year, type = "median"), pairwise = TRUE)
#not significant

## see if any species are specific to a particular depth
# #see how many species are unique to each depth
# unique_by_depth <- combined_tidy %>% 
#   group_by(Depth) %>% 
#   distinct(Recode_for_MS) %>% 
#   mutate(present = 1) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = Depth, values_from = present, values_fill = 0) %>% 
#   mutate(depths_present = rowSums(select_if(., is.numeric))) %>% 
#   arrange(depths_present)
# 
# unique_group_by_depth <- combined_tidy %>% 
#   group_by(Depth) %>% 
#   distinct(Recode_for_MS) %>% 
#   mutate(present = 1) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = Depth, values_from = present, values_fill = 0) %>% 
#   mutate(depths_present = rowSums(select_if(., is.numeric))) %>% 
#   arrange(depths_present)


abund_10 <- nonmotile_tidy_w %>%
  select(!contains("unidentified")) %>% #not including unidentified orgs
  filter(year == "10") %>% 
  select(amphipod_tube:urticina_sp) 

meta_10 <- nonmotile_tidy_w %>%
  filter(year == "10") %>% 
  select(year:replicate)

ISA_depth <- multipatt(x = abund_10, cluster = meta_10$site, duleg = TRUE)
summary(ISA_depth)


## SIMPER 
SIMPER_depth <- simper(comm = abund_10, group = meta_10$site, permutations = 9999)
summary(SIMPER_depth)

SIMPER_test <- keep(SIMPER_depth, cumsum <0.7)

#### motile analysis ####
motile_tidy_w <- motile_tidy %>% 
  arrange(Recode_for_MS) %>% 
  select(!MS_Species_Group) %>% 
  pivot_wider(names_from = Recode_for_MS, values_from = Species_count, values_fill = 0) %>% 
  clean_names()

#prep final raw dataframes
motile.abund <- motile_tidy_w %>% 
  select(!year:phylum)

motile.meta <- motile_tidy_w %>% 
  select(year:phylum)

CV(x = rowSums(motile.abund)) #high, indicates there would be high impact of relativizing between samples
CV(x = colSums(motile.abund)) #high, indicates there would be high impact of relativizing between species

motile.abund_rel <- wisconsin(motile.abund) #double standardization: standardizes by row total and max of columns

# motile.nmds <- metaMDS(motile.abund_rel,
#                        distance = "bray", 
#                        autotransform = FALSE,
#                        engine = "monoMDS",
#                        k = 3,
#                        weakties = TRUE,
#                        model = "global",
#                        maxit = 400,
#                        try = 40,
#                        trymax = 100, 
#                        trace = FALSE)
#post-hoc test
# stressplot(object = motile.nmds, lwd = 5)
# 
# #tidyverse
# motile.nmds_points <- data.frame(motile.nmds$points)
# motile.nmds_points <- bind_cols(motile.meta, motile.nmds_points) %>% 
#   mutate(depth = factor(depth, labels = c("30", "90", "200", "200 ref")),
#          year = factor(year))

#### Figure 5 ####
#motile nmds
# ggplot(data=motile.nmds_points,
#        aes(x=MDS1, y=MDS2,
#            fill= depth,
#            shape= year)) + 
#   geom_point(color = "black", size = 3) +
#   stat_ellipse(aes(group = depth, color = depth), 
#                linetype = "dashed", show.legend = FALSE) +
#   theme(axis.line = element_blank(), 
#         axis.ticks = element_blank(),
#         axis.text =  element_blank()) +
#   theme_classic() +
#   labs(fill = "Site", shape = "Year") + 
#   scale_shape_manual(values = c(21, 24, 22)) +
#   scale_fill_manual(values = depth_colors) + 
#   scale_color_manual(values = depth_colors) +
#   guides(fill=guide_legend(override.aes=list(color=c(depth_colors)))) +
#   annotate("text", x = 2.1, y = -2.7, 
#            label = paste("Stress = ", round(motile.nmds$stress, 3))) +
#   theme(text = element_text(size = 14))

# ggsave("figures/figure5.tiff", width = 8, height = 6, dpi = 300)

#### motile multi analysis ####

## compare depths at the outfall outfall
motile.abund_outfall <- motile_tidy_w %>%
  filter(!depth == "600ft_ref") %>% 
  select(amphipod_shrimp:sea_urchin)

motile.abund_outfall_rel <- wisconsin(motile.abund_outfall)

motile.meta_outfall <- motile_tidy_w %>% 
  filter(!depth == "600ft_ref") %>% 
  select(year:replicate)

adonis2(motile.abund_outfall_rel ~ depth:year, method = "bray",
        data = motile.meta_outfall, permutations = 9999, by = "margin")
#interaction is significant

## compare the outfall and reference 200 m sites
motile.abund_ref <- motile_tidy_w %>%
  filter(depth == 600 | depth == "600_ref") %>%
  select(amphipod_shrimp:sea_urchin)

motile.abund_ref_rel <- wisconsin(motile.abund_ref)

motile.meta_ref <- motile_tidy_w %>% 
  filter(depth == 600 | depth == "600_ref") %>%
  select(year:replicate)

adonis2(motile.abund_ref_rel ~ depth:year, method = "bray",
        data = motile.meta_ref, permutations = 9999, by = "margin")

#interaction is not significant; test for main effects 
adonis2(motile.abund_ref_rel ~ depth + year, method = "bray",
        data = motile.meta_ref, permutations = 9999, by = "margin")
#no significant effect

## SIMPER comparing year 2 to year 10 across all sites
motile.SIMPER_year <- simper(comm = motile.abund_rel, group = motile.meta$year, permutations = 9999)

summary(motile.SIMPER_year)$'2_10'


