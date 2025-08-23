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
combined_import <- bind_rows(import_100ft, import_300ft, import_600ft, import_600ft_ref) #combine imported data into a single dataframe

#extract QC rows 
QC_data <- combined_import %>% #these rows are repeats of data we already have so we can set them aside for now
  filter(grepl("QC", Notes))

#tidy data
combined_tidy <- combined_import %>% 
  filter(!is.na(Year)) %>% #remove empty rows
  anti_join(QC_data) %>% #remove QC rows
  filter(!Year == 2017) %>% #remove for this analysis because there wasn't enough replication
  mutate(Depth = as.character(Depth), 
         Replicate = as.factor(Replicate), 
         Year = factor(Year, labels = c("2", "10"))) %>% #recode as deployment intervals
  mutate(Depth = ifelse(grepl("R", Grid.ID), paste0(Depth, "_ref"), Depth)) %>% #recode reference plates
  full_join(species_metadata, by = "SpeciesID")

#create motile and nonmotile dataframes for separate analysis
motile_tidy <- combined_tidy %>% filter(Mobility == "Motile")
nonmotile_tidy <- combined_tidy %>% filter(Mobility == "Nonmotile")

#### summary of nonmotile phyla ####
#phylum proportion values 
area_sum <- nonmotile_tidy %>% 
  group_by(Year, Depth) %>% 
  mutate(total_live_area = sum(Org.Area.cm)) %>% 
  dplyr::select(Year, Depth, total_live_area) %>% 
  distinct()

#breakdown by organism category
group_perc_live_cover <- nonmotile_tidy %>% 
  group_by(Year, Depth, Species_Group) %>% 
  summarize(group_area = sum(Org.Area.cm)) %>% 
  right_join(area_sum) %>% 
  mutate(group_prop = round(group_area/total_live_area, 4)*100)

#focus on dominant groups
top_10group_perc <- group_perc_live_cover %>% 
  filter(group_prop > 10) 

#plot proportions of nonmotile taxa
nonmotile_tidy %>% 
  ggplot(aes(x = factor(Depth, 
                        # labels = c("30 m\noutfall", "90 m\noutfall", "200 m\noutfall", "200 m\nreference")), 
             labels = c("30", "90", "200", "200\nref")), 
             y = Org.Area.cm, fill = Genera)) +
  geom_col(position = "fill") +
  theme_classic() + 
  facet_wrap(~Year) + 
  # theme(text = element_text(size = 14)) +
  labs(y = "Proportion of total live area", x = "Site", fill = "Phylum") +
  scale_fill_manual(values = rev(brewer.pal(n = 11, "Spectral"))) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0, vjust = 1))

# ggsave("figures/figure2.png", width = 8, height = 6, dpi = 300)

#### summary plots based on % cover ####
#calculate nonmotile percent live cover
perc_cover_df <- nonmotile_tidy %>% 
  group_by(Year, Depth, Replicate) %>% 
  summarize(grid.area.sum = sum(Grid.Area.cm), #sum the total area across all grids on a given plate
            live.area.covered = sum(Org.Area.cm)) %>% 
  ungroup() %>% #we're done summarizing data so we can remove the groupings
  mutate(perc.live.cover = round((live.area.covered/grid.area.sum)*100, 2)) #divide the area covered by the total area to get % cover

depth_colors <- c("#fde725","#21918c", "#440154", "#fc8961")

## Figure 2
perc_cover_df %>%
  ggplot(aes(x = Year, y = perc.live.cover, fill = Depth)) +
  geom_point(size = 3, shape = 21) +
  theme_classic() +
  labs(y = "Percent live cover",
       x = "Pipe material deployment time (years)",
       fill = "Site") +
  scale_fill_manual(values = depth_colors,
                    labels = c("30", "90", "200", "200 ref")) +
  theme(text = element_text(size = 14))

# ggsave("figures/figure3.png", width = 8, height = 6, dpi = 300)

#### univariate stats ####
shapiro.test(x = perc_cover_df$perc.live.cover)
shapiro.test(x = log(perc_cover_df$perc.live.cover))

#test the outfall sites only 
perc_cover_mod <- aov(log(perc.live.cover) ~ Depth + Year + Depth*Year, data = filter(perc_cover_df, !Depth == "600_ref"))
summary(perc_cover_mod)

#post-hoc test
perc_cover_mod2 <- aov(log(perc.live.cover) ~ Depth + Year, data = filter(perc_cover_df, !Depth == "600_ref"))
TukeyHSD(perc_cover_mod2)

#model diagnostics
aov.res <- residuals(perc_cover_mod)

hist(perc_cover_mod$residuals)
qqnorm(perc_cover_mod$residuals)
qqline(perc_cover_mod$residuals)

#test for differences between the -200 m outfall and reference sites
compare_tidy <- perc_cover_df %>% 
  filter(Depth %in% c("600", "600_ref")) 

shapiro.test(x = compare_tidy$perc.live.cover)

summary(aov(perc.live.cover ~ Depth + Year + Depth*Year, data = compare_tidy))

compare_live_wo_5 <- compare_tidy %>% 
  filter(!Year == 5)

summary(aov(perc.live.cover ~ Depth + Year + Depth*Year, data = compare_live_wo_5))

#### nonmotile multivariate analysis ####

#convert the data to wide format
nonmotile_tidy_w <- nonmotile_tidy %>% 
  group_by(Year, Depth, Replicate, Species_Group) %>% 
  summarize(sum_area = sum(Area.Covered.cm)) %>% 
  filter(!is.na(Species_Group)) %>% 
  ungroup() %>% 
  arrange(Species_Group) %>% 
  pivot_wider(names_from = Species_Group, values_from = sum_area, values_fill = 0) %>% 
  clean_names() %>% 
  arrange(year, depth, replicate)

#prep abundance and metadata matrices

abund <- nonmotile_tidy_w %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta <- nonmotile_tidy_w %>%
  select(year:replicate)

#conduct data transformations and/or relativizations 
#not removing rare species because we're looking for differences in diversity

#function for calculating the coefficient of variation
CV <- function(x) { 100 * sd(x) / mean(x) }

CV(x = rowSums(abund)) #low (<50), relativizing by row will not make a big difference
CV(x = colSums(abund)) #medium (50-100), could have an effect to relativize between species

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
  mutate(depth = factor(depth, labels = c("30", "90", "200", "200 ref")),
         year = factor(year))

#nmds plot
ggplot(data=nmds_points,
       aes(x=MDS1, y=MDS2,
           fill= depth,
           shape= year)) + 
  geom_point(color = "black", size = 3) +
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
#warning message isn't a problem

# ggsave("figures/figure5.png", width = 8, height = 6, dpi = 300)

#### Nonmotile multi stats ####
## compare depths at the pipe
abund_outfall <- nonmotile_tidy_w %>%
  filter(!depth == "600_ref") %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta_outfall <- nonmotile_tidy_w %>%
  filter(!depth == "600_ref") %>%
  select(year:replicate)

adonis2(abund_outfall ~ depth:year, method = "bray",
        data = meta_outfall, permutations = 9999, by = "margin")
#significant interaction term

permutest(betadisper(vegdist(abund_outfall, method = "bray"), group = meta_outfall$depth, type = "median"))
#significant

year.beta <- betadisper(vegdist(abund_outfall, method = "bray"), group = meta_outfall$depth, type = "median")
permutest(year.beta, pairwise = TRUE)

permutest(betadisper(vegdist(abund_outfall, method = "bray"), group = meta_outfall$year, type = "median"))
#not significant

## compare the -200 m outfall and reference sites
abund_ref <- nonmotile_tidy_w %>%
  filter(depth == 600 | depth == "600_ref") %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta_ref <- nonmotile_tidy_w %>%
  filter(depth == 600 | depth == "600_ref") %>%
  select(year:replicate)

adonis2(abund_ref ~ depth:year, method = "bray",
        data = meta_ref, permutations = 9999, by = "margin")
#significant interaction term

permutest(betadisper(vegdist(abund_ref, method = "bray"), group = meta_ref$depth, type = "median"))
permutest(betadisper(vegdist(abund_ref, method = "bray"), group = meta_ref$year, type = "median"))
#neither are significant

## see if any species are specific to a particular depth
#see how many species are unique to each depth
unique_by_depth <- combined_tidy %>% 
  group_by(Depth) %>% 
  distinct(SpeciesID) %>% 
  mutate(present = 1) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Depth, values_from = present, values_fill = 0) %>% 
  mutate(depths_present = rowSums(select_if(., is.numeric))) %>% 
  arrange(depths_present)

unique_group_by_depth <- combined_tidy %>% 
  group_by(Depth) %>% 
  distinct(Species_Group) %>% 
  mutate(present = 1) %>% 
  ungroup() %>% 
  pivot_wider(names_from = Depth, values_from = present, values_fill = 0) %>% 
  mutate(depths_present = rowSums(select_if(., is.numeric))) %>% 
  arrange(depths_present)

ISA_depth <- multipatt(x = abund, cluster = meta$depth, duleg = TRUE)
summary(ISA_depth)

nonmotile_tidy %>% 
  group_by(Depth) %>%
  distinct(SpeciesID)

## SIMPER 
SIMPER_depth <- simper(comm = abund, group = meta$depth, permutations = 9999)
summary(SIMPER_depth)

#see which species are contributing to differences in years
SIMPER_year <- simper(comm = abund, group = meta$year, permutations = 9999)

summary(SIMPER_year)$'2_10' %>%
  round(3) %>%
  head()

#### motile analysis ####
motile_tidy_w <- motile_tidy %>% 
  group_by(Year, Depth, Replicate, Species_Group) %>% 
  summarize(sum_area = sum(Species.Count)) %>% 
  filter(!is.na(Species_Group)) %>% 
  ungroup() %>% 
  arrange(Species_Group) %>% 
  pivot_wider(names_from = Species_Group, values_from = sum_area, values_fill = 0) %>% 
  clean_names()

#prep final raw dataframes
motile.abund <- motile_tidy_w %>% 
  select(amphipod_shrimp:sea_urchin)

motile.meta <- motile_tidy_w %>% 
  select(year:replicate)

CV(x = rowSums(motile.abund)) #high, indicates there would be high impact of relativizing between samples
CV(x = colSums(motile.abund)) #high, indicates there would be high impact of relativizing between species

motile.abund_rel <- wisconsin(motile.abund) #double standardization: standardizes by row total and max of columns

motile.nmds <- metaMDS(motile.abund_rel,
                       distance = "bray", 
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
stressplot(object = motile.nmds, lwd = 5)

#tidyverse
motile.nmds_points <- data.frame(motile.nmds$points)
motile.nmds_points <- bind_cols(motile.meta, motile.nmds_points) %>% 
  mutate(depth = factor(depth, labels = c("30", "90", "200", "200 ref")),
         year = factor(year))

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

# ggsave("figures/figure6.png", width = 8, height = 6, dpi = 300)

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

## compare motile dispersion
permutest(betadisper(vegdist(motile.abund_outfall_rel, method = "bray"), group = motile.meta_outfall$depth, type = "median"))
permutest(betadisper(vegdist(motile.abund_outfall_rel, method = "bray"), group = motile.meta_outfall$year, type = "median"))
permutest(betadisper(vegdist(motile.abund_ref_rel, method = "bray"), group = motile.meta_ref$depth, type = "median"))
permutest(betadisper(vegdist(motile.abund_ref_rel, method = "bray"), group = motile.meta_ref$year, type = "median"))
#none are significant

## SIMPER comparing year 2 to year 10 across all sites
motile.SIMPER_year <- simper(comm = motile.abund_rel, group = motile.meta$year, permutations = 9999)

summary(motile.SIMPER_year)$'2_10'


