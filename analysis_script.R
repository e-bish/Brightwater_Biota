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
  mutate(Depth = as.character(Depth), 
         Replicate = as.factor(Replicate), 
         Year = factor(Year, labels = c("2", "5", "10"))) %>% #recode as deployment intervals
  mutate(Depth = ifelse(grepl("R", Grid.ID), paste0(Depth, "_ref"), Depth)) %>% #recode reference plates
  full_join(species_metadata, by = "SpeciesID")

#create motile and nonmotile dataframes for separate analysis
motile_tidy <- combined_tidy %>% filter(Mobility == "Motile")
nonmotile_tidy <- combined_tidy %>% filter(Mobility == "Nonmotile")

#### summary of nonmotile phyla ####
nonmotile_tidy %>% 
  filter(!Depth == "600_ref") %>% 
  ggplot(aes(x = factor(Depth, labels = c("30", "90", "200")), 
             y = Org.Area.cm, fill = Genera)) +
  geom_col(position = "fill") +
  theme_classic() + 
  facet_wrap(~Year) + 
  # theme(text = element_text(size = 14)) +
  labs(y = "Proportion of total live area", x = "Depth (m)", fill = "Phylum") +
  scale_fill_manual(values = rev(brewer.pal(n = 11, "Spectral")))

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

## Figure 1 
perc_cover_df %>%
  filter(!Depth == "600_ref") %>% 
  ggplot(aes(x = Year, y = perc.live.cover, fill = Depth)) +
  geom_point(size = 3, shape = 21) +
  theme_classic() +
  labs(y = "Percent live cover", 
       x = "Pipe material deployment time (years)", 
       color = "Depth (ft)") +
  scale_fill_manual(values = depth_colors[1:3], 
                    labels = c("30 m", "90 m", "200 m"))

# ggsave("figures/figure3.png", width = 8, height = 6, dpi = 300)

#compare cover to reference location
perc_cover_df %>% 
  filter(Depth %in% c("600", "600_ref")) %>% 
  group_by(Year, Depth) %>% 
  summarize(mean = mean(perc.live.cover), sd = sd(perc.live.cover)) %>% 
  ggplot(aes(x = Year, y = mean, fill = factor(Depth, labels = c("Outfall pipe", "Reference")))) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.95) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width=.1, position = position_dodge(width = .9)) +
  theme_classic() +
  labs(y = "Average percent live cover at -200 m MLLW", x = "Pipe material deployment time (years)", fill = "Site") +
  scale_fill_manual(values = depth_colors[3:4])

# ggsave("figures/figure4.png", width = 8, height = 6, dpi = 300)

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
  mutate(depth = factor(depth, labels = c("30 m pipe", "90 m pipe", "200 m pipe", "200 m reference")),
         year = factor(year))

#nmds plot
ggplot(data=nmds_points,
       aes(x=MDS1, y=MDS2,
           fill= depth,
           shape= year)) + 
  geom_point(color = "black", size = 3) +
  stat_ellipse(aes(group = depth, color = depth), 
               linetype = "dashed", show.legend = FALSE) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  labs(fill = "Site", shape = "Year") + 
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_manual(values = depth_colors) + 
  scale_color_manual(values = depth_colors) +
  guides(fill=guide_legend(override.aes=list(color=c(depth_colors)))) +
  annotate("text", x = -1.2, y = 1.4, 
           label = paste("Stress = ", round(nonmotile.nmds$stress, 3)))
#warning message isn't a problem

# ggsave("figures/figure5.png", width = 8, height = 6, dpi = 300)

#### Nonmotile multi stats ####
## compare depths at the pipe
abund_pipe <- nonmotile_tidy_w %>%
  filter(!depth == "600_ref") %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta_pipe <- nonmotile_tidy_w %>%
  filter(!depth == "600_ref") %>%
  select(year:replicate)

adonis2(abund_pipe ~ depth:year, method = "bray",
        data = meta_pipe, permutations = 9999, by = "margin")
#significant interaction term

## compare the -200 m pipe and reference sites
abund_ref <- nonmotile_tidy_w %>%
  filter(depth == 600 | depth == "600_ref") %>%
  select(amphipod_tube:ulva_sp) #not including unidentified orgs

meta_ref <- nonmotile_tidy_w %>%
  filter(depth == 600 | depth == "600_ref") %>%
  select(year:replicate)

adonis2(abund_ref ~ depth:year, method = "bray",
        data = meta_ref, permutations = 9999, by = "margin")
#significant interaction term

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
  mutate(depth = factor(depth, labels = c("30 m pipe", "90 m pipe", "200 m pipe", "200 m reference")),
         year = factor(year))

#motile nmds
ggplot(data=motile.nmds_points,
       aes(x=MDS1, y=MDS2,
           fill= depth,
           shape= year)) + 
  geom_point(color = "black", size = 3) +
  stat_ellipse(aes(group = depth, color = depth), 
               linetype = "dashed", show.legend = FALSE) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  labs(fill = "Site", shape = "Year") + 
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_manual(values = depth_colors) + 
  scale_color_manual(values = depth_colors) +
  guides(fill=guide_legend(override.aes=list(color=c(depth_colors)))) +
  annotate("text", x = 2.1, y = -2.7, 
           label = paste("Stress = ", round(motile.nmds$stress, 3)))

# ggsave("figures/figure6.png", width = 8, height = 6, dpi = 300)

