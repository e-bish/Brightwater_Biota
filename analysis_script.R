# load libraries
library(tidyverse)
library(janitor)
library(vegan)
library(ggordiplots)
library(indicspecies)
library(viridis)
library(RColorBrewer)

# load data
import_100ft <- read_csv("raw_data/import_100ft.csv")
import_300ft <- read_csv("raw_data/import_300ft.csv")
import_600ft <- read_csv("raw_data/import_600ft.csv")
import_600ft_ref <- read_csv("raw_data/import_600ft_ref.csv")
species_metadata <- read_csv("raw_data/species_metadata.csv")

set.seed(2025)

# prep dataframes
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
  mutate(Depth = as.character(Depth), 
         Replicate = as.factor(Replicate), 
         Year = factor(Year, labels = c("2", "5", "10"))) %>% 
  mutate(Depth = ifelse(grepl("R", Grid.ID), paste0(Depth, "_ref"), Depth)) %>% 
  full_join(species_metadata, by = "SpeciesID")

#species identified in 2014 report but not the current effort at the outfall pipe sites
not_found <- combined_tidy %>% filter(is.na(Depth))

#remove rows with empty depth
combined_tidy <- combined_tidy %>% 
  anti_join(not_found)

motile_tidy <- combined_tidy %>% filter(Mobility == "Motile")
nonmotile_tidy <- combined_tidy %>% filter(Mobility == "Nonmotile")

## summary plots
#percent cover
perc_cover_df <- nonmotile_tidy %>% 
  group_by(Year, Depth, Replicate) %>% #this step tells R that we're interested in summarizing values across unique combinations of year, depth, and replicate
  summarize(grid.area.sum = sum(Grid.Area.cm), #sum the total area across all grids on a given tile
            area.covered.sum = sum(Area.Covered.cm), #sum the area covered across all grids on a given tile
            live.area.covered = sum(Org.Area.cm)) %>% 
  ungroup() %>% #we're done summarizing data so we can remove the groupings
  mutate(perc.cover = round((area.covered.sum/grid.area.sum)*100, 2),
         perc.live.cover = round((live.area.covered/grid.area.sum)*100, 2)) #divide the area covered by the total area to get % cover

depth_colors <- c("#fde725","#21918c", "#440154", "#fc8961")

perc_cover_df %>%
  filter(!Depth == "600_ref") %>% 
  ggplot(aes(x = Year, y = perc.live.cover, fill = Depth)) +
  geom_point(size = 3, shape = 21) +
  theme_classic() +
  labs(y = "Percent live cover", x = "Pipe material deployment time (years)", color = "Depth (ft)") +
  scale_fill_manual(values = depth_colors[1:3])

#compare cover

perc_cover_df %>% 
  filter(Depth %in% c("600", "600_ref")) %>% 
  group_by(Year, Depth) %>% 
  summarize(mean = mean(perc.live.cover), sd = sd(perc.live.cover)) %>% 
  ggplot(aes(x = Year, y = mean, fill = factor(Depth, labels = c("Outfall pipe", "Reference")))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.95) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.1, position = position_dodge(width = .9)) +
  theme_classic() +
  labs(y = "Average percent live cover at -600 ft MLLW", x = "Pipe material deployment time (years)", fill = "Site") +
  scale_fill_manual(values = depth_colors[3:4])

#phyla
nonmotile_tidy %>% 
  filter(!Depth == "600_ref") %>% 
  ggplot(aes(x = Depth, y = Org.Area.cm, fill = Genera)) +
  geom_col(position = "fill") +
  theme_classic() + 
  facet_wrap(~Year) + 
  # theme(text = element_text(size = 14)) +
  labs(y = "Proportion of total live area", x = "Depth (ft)", fill = "Phylum") +
  scale_fill_manual(values = rev(brewer.pal(n = 11, "Spectral")))

#nonmotile analysis

nonmotile_tidy_w <- nonmotile_tidy %>% 
  group_by(Year, Depth, Replicate, Species_Group) %>% 
  summarize(sum_area = sum(Area.Covered.cm)) %>% 
  filter(!is.na(Species_Group)) %>% 
  ungroup() %>% 
  arrange(Species_Group) %>% 
  pivot_wider(names_from = Species_Group, values_from = sum_area, values_fill = 0) %>% 
  clean_names() %>% 
  arrange(year, depth, replicate)

#prep final raw dataframes

abund <- nonmotile_tidy_w %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta <- nonmotile_tidy_w %>%
  select(year:replicate)

#conduct data transformations and/or relativizations
#not removing rare species because we're looking for differences in diversity

CV <- function(x) { 100 * sd(x) / mean(x) }

CV(x = rowSums(abund2)) #low, relativizing will not make a big difference
CV(x = colSums(abund2)) #high, could have an effect to relativize between species


#### NMDS #### 
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

#base R
# plot(nonmotile.nmds, "sites")
# with(abund, 
#      points(nonmotile.nmds,
#             display = "sites"))
# ordihull(nonmotile.nmds,
#          meta$depth, display = "sites", draw = c("polygon"))

#tidyverse
nmds_points <- data.frame(nonmotile.nmds$points)
nmds_points <- bind_cols(meta, nmds_points)
# nmds_depth_hulls <- nmds_points %>% 
#   group_by(depth) %>% 
#   slice(chull(MDS1, MDS2))

# nmds_year_hulls <- nmds_points %>% 
#   group_by(year) %>% 
#   slice(chull(MDS1, MDS2))

# stress_df <- data.frame(stress = nonmotile.nmds$stress, x = 2, y = 1)
#stress values less than 0.2 but greater than 0.10 are decent but shouldn't be relied on for details

ggplot(data=nmds_points,
       aes(x=MDS1, y=MDS2,
           fill= factor(meta$depth, labels = c("100", "300", "600", "600 reference")),
           shape= meta$year)) + 
  geom_point(color = "black", size = 3) +
  stat_ellipse(aes(group = meta$depth, color = meta$depth), 
               linetype = "dashed", show.legend = FALSE) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  labs(fill = "Depth", shape = "Year") + 
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_manual(values = depth_colors) + 
  scale_color_manual(values = depth_colors) +
  guides(fill=guide_legend(override.aes=list(color=c(depth_colors))))

#### alt interaction plot code. Leaving here for reference but decided to use betadisper instead because 
#the documentation seems to match up with Anderson et al. 2017 ###

nmds_points <- nmds_points %>% 
  mutate(int = paste(year, depth, sep = "_"))

## using envfit
centroids <- envfit(nonmotile.nmds,nmds_points[,"int"])[["factors"]][["centroids"]]

nonmotile.int.nmds <- metaMDS(centroids, distance = "euc", 
                              autotransform = FALSE,
                              engine = "monoMDS",
                              k = 3,
                              weakties = TRUE,
                              model = "global",
                              maxit = 400,
                              try = 40,
                              trymax = 100)

nmds_int_points <- data.frame(nonmotile.int.nmds$points) %>% 
  as.data.frame() %>% 
  mutate(depth = rep(c("100", "300", "600", "600_ref"), times = 3),
         year = rep(c("10", "2", "5"), each = 4)) %>% 
  mutate(year = factor(year, levels = c("2", "5", "10")))


ggplot(data=nmds_int_points,
       aes(x=MDS1, y=MDS2, 
           # label = rownames(centroids), 
           color = depth, 
           shape = year)) +
  geom_point(size = 3) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic()

## using betadisper
abund_dist <- vegdist(abund, method = "bray")
centroids <- betadisper(abund_dist, group = nmds_points$int, type = "centroid")$centroids

centroids_points <- data.frame(centroids) %>% 
  as.data.frame() %>% 
  select(!PCoA3:PCoA31) %>% 
  mutate(depth = rep(c("100", "300", "600", "600_ref"), times = 3),
         year = rep(c("10", "2", "5"), each = 4)) %>% 
  mutate(year = factor(year, levels = c("2", "5", "10")))


nonmotile.int.nmds <- metaMDS(centroids, distance = "euc", 
                              autotransform = FALSE,
                              engine = "monoMDS",
                              k = 3,
                              weakties = TRUE,
                              model = "global",
                              maxit = 400,
                              try = 40,
                              trymax = 100)

# nonmotile.int.nmds$stress

nonmotile.int.points <- nonmotile.int.nmds$points

ggplot(data=nonmotile.int.points ,
       aes(x=MDS1, y=MDS2, 
           # label = rownames(centroids), 
           fill = centroids_points$depth, 
           shape = centroids_points$year)) +
  geom_point(size = 3) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  labs(fill = "Depth", shape = "Year") + 
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_manual(values = depth_colors) + 
  scale_color_manual(values = depth_colors) +
  guides(fill=guide_legend(override.aes=list(color=c(depth_colors)))) 

#### Nonmotile permanova ####
## compare 600_ref separately
abund_out <- nonmotile_tidy_w %>%
  filter(!depth == "600_ref") %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta_out <- nonmotile_tidy_w %>%
  filter(!depth == "600_ref") %>%
  select(year:replicate)


abund_comp <- nonmotile_tidy_w %>%
  filter(depth == 600 | depth == "600_ref") %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta_comp <- nonmotile_tidy_w %>%
  filter(depth == 600 | depth == "600_ref") %>%
  select(year:replicate)

#make the distance matrix separately because the argument for it in pairwise Adonis isn't working
nonmotile_dist_out <- vegdist(abund_out, method = "bray")
adonis2(nonmotile_dist_out ~ depth:year, data = meta_out, permutations = 9999, by = "margin")

nonmotile_dist_comp <- vegdist(abund_comp, method = "bray")
adonis2(nonmotile_dist_comp ~ depth:year, data = meta_comp, permutations = 9999, by = "margin")

## see if any species are specific to a particular depth
ISA_depth <- multipatt(x = abund, cluster = meta$depth, duleg = TRUE)
summary(ISA_depth)

ISA_year <- multipatt(x = abund, cluster = meta$year, duleg = TRUE)
summary(ISA_year)

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

motile.nmds <- metaMDS(motile.abund_rel, distance = "bray")
#post-hoc test
stressplot(object = motile.nmds, lwd = 5)

#base R
plot(motile.nmds, "sites")
with(motile.abund_rel, 
     points(motile.nmds,
            display = "sites"))
ordihull(motile.nmds,
         motile.meta$depth, display = "sites", draw = c("polygon"))

#tidyverse
motile.nmds_points <- data.frame(motile.nmds$points)
motile.nmds_points <- bind_cols(motile.meta, motile.nmds_points)


calculate_hull <- function(df) {
  df[chull(df$MDS1, df$MDS2), ]
}

motile.nmds_depth_hulls <- motile.nmds_points %>% 
  group_by(factor(depth)) %>% 
  group_split() %>%
  map_dfr(~ calculate_hull(.), .id = "hull_no") 

motile.nmds_year_hulls <- motile.nmds_points %>% 
  group_by(factor(year)) %>% 
  group_split() %>%
  map_dfr(~ calculate_hull(.), .id = "hull_no") 

stress_df <- data.frame(stress = motile.nmds$stress, x = 2, y = 1)
#stress values less than 0.2 but greater than 0.10 are decent but shouldn't be relied on for details

ggplot(data=motile.nmds_points,
       aes(x=MDS1, y=MDS2,
           fill= factor(motile.meta$depth, labels = c("100", "300", "600", "600 reference")),
           shape= motile.meta$year)) + 
  geom_point(color = "black", size = 3) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  labs(fill = "Depth", shape = "Year") + 
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_manual(values = depth_colors) + 
  scale_color_manual(values = depth_colors) + 
  guides(fill=guide_legend(override.aes=list(color=c(depth_colors))))
  # geom_polygon(data = motile.nmds_depth_hulls,
  #              aes(x = MDS1, y = MDS2, 
  #              fill = hull_no, group = depth), 
  #            show.legend = FALSE, alpha = 0.75, inherit.aes = FALSE)


# anosim(abund,
#   meta$depth,
#   permutations = 999,
#   distance = "bray",
#   strata = NULL,
#   parallel = getOption("mc.cores")
# )
# 
# anosim(abund,
#        meta$year,
#        permutations = 999,
#        distance = "bray",
#        strata = NULL,
#        parallel = getOption("mc.cores")
# )
#The test statistic (R) ranges from -1 (all lowest ranks are between groups – an unusual situation) to +1 (all lowest ranks are within groups).

