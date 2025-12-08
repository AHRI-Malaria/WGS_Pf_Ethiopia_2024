# Install necessary packages
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(dplyr)
library(sf)
library(ggforce)
library(ggpattern) 
library(rnaturalearth)
library(geodata)
library(readr)
library(ggrepel)
library(ggnewscale)

# Load Ethiopia outline (national borders)
ethiopia <- ne_countries(scale = "medium", country = "Ethiopia", returnclass = "sf")

# Get Ethiopia admin level 1 (regions)
ethiopia_admin1 <- geodata::gadm("GADM", country = "ETH", level = 1)

ethiopia_admin1_sf <- st_as_sf(ethiopia_admin1)  # Convert to sf for ggplot2

# Create a world background (to make the outside grey)
world <- ne_countries(scale = "medium", returnclass = "sf")

districts <- data.frame(
  site = c("Mizan", "Gambella", "Gondar","Arba Minch", "Metehera Fentale"),
  lon = c(35.583, 34.583, 37.467, 37.5573,39.9170),
  lat = c(6.15, 8.25, 12.6, 6.0350,8.8995)
)
samples <- read.csv("../input_files/Metadata_Popgen_Ethiopia.csv", sep=",", header=TRUE)
sample_counts <- samples %>%
  count(site) %>%
  rename(site = site, sample_count = n)

districts <- left_join(districts, sample_counts, by = "site")
districts$sample_count[is.na(districts$sample_count)] <- 0  # Replace NAs with 0

# Count samples per site-study
study_counts <- samples %>%
  count(site, Study) %>%
  group_by(site) %>%
  mutate(
    total = sum(n),
    frac = n / total,
    angle = 2 * pi * frac,
    start = cumsum(lag(angle, default = 0)),
    end = cumsum(angle)
  ) %>%
  left_join(districts, by = c("site"))


# Plot map site names, sample count
ggplot() +
  geom_sf(data = world, fill = "grey90", color = NA) +
  geom_sf(data = ethiopia_admin1_sf, fill = "white", color = "black", size = 0.8) +
  geom_point(data = districts, aes(x = lon, y = lat, size = sample_count, color = site)) +
  scale_size_continuous(
    name = "Sample count",      
    range = c(2, 8),
    breaks = c(10, 20, 30,40),
    limits = c(10, 40)
  ) +
  scale_color_discrete(name = "Site") + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  coord_sf(xlim = c(32, 48), ylim = c(3, 15)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 1.5, size = 14)
  )
ggsave("figs/map_sites_samplesize.pdf", width = 6.5, height = 5.5, dpi = 300)


# Plot map with studies
ggplot() +
  geom_sf(data = world, fill = "grey90", color = NA) +
  geom_sf(data = ethiopia_admin1_sf, fill = "white", color = "black", size = 0.8) +
  geom_arc_bar(
    data = study_counts,
    aes(
      x0 = lon, y0 = lat,
      r0 = 0, r = 0.4,       # Adjust radius size here
      start = start, end = end,
      fill =Study
    ),
    color = "black", size = 0.2
  ) +
  scale_fill_manual(
    name = "Study",  # Optional: sets legend title
    values = c(
      "TES" = "darkgreen",  
      "Leucodepletion" = "lightblue"  
    )
  )+
  coord_sf(xlim = c(32, 48), ylim = c(3, 15)) +
  theme_void() +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )                   


ggsave("figs/map_studies.pdf", width = 8, height = 6, dpi = 300)



# Load Africa countries
africa <- ne_countries(scale = "medium", continent = "Africa", returnclass = "sf")

# Get Ethiopia admin level 1 (regions) - to highlight Ethiopia
ethiopia_admin1 <- geodata::gadm("GADM", country = "ETH", level = 1)
ethiopia_admin1_sf <- st_as_sf(ethiopia_admin1)

# Your 5 main Ethiopian districts to highlight with specific colors
main_district_names <- c("Mizan", "Gambella", "Gondar", "Arba Minch", "Metehera Fentale")

# Load the comprehensive database
all_samples <- read.csv("../input_files/Popgen_db_samples_072025_filt_FINAL.csv", sep = ",", header = TRUE)

# Get all sites from the comprehensive database
all_sites <- all_samples %>%
  filter(!is.na(Admin.level.1.latitude) & !is.na(Admin.level.1.longitude)) %>%
  distinct(District, Admin.level.1.latitude, Admin.level.1.longitude) %>%
  rename(lon = Admin.level.1.longitude, lat = Admin.level.1.latitude)

# Count samples per site in comprehensive database
all_sample_counts <- all_samples %>%
  filter(!is.na(Admin.level.1.latitude) & !is.na(Admin.level.1.longitude)) %>%
  count(District) %>%
  rename(sample_count = n)

# Merge coordinates with sample counts
all_sites <- left_join(all_sites, all_sample_counts, by = "District")

# Mark which sites are your main 5 Ethiopian districts
all_sites <- all_sites %>%
  mutate(main_site = District %in% main_district_names)

# Filter out invalid coordinates (Mali entries with dates)
# Keep only numeric coordinates within valid lat/lon ranges
all_sites <- all_sites %>%
  filter(
    is.numeric(lon) & is.numeric(lat) &
      lon >= -180 & lon <= 180 &
      lat >= -90 & lat <= 90
  )

# Define which countries to highlight (black borders, white fill)

highlight_countries <- unique(all_samples$Country)

country_colors <- c(
  "Ethiopia" = "#E41A1C", 
  "Kenya" = "#377EB8", 
  "Tanzania" = "#4DAF4A",
  "Uganda" = "#984EA3",   
  "Sudan" = "#FF7F00",  
  "Ghana" = "#FFFF33", 
  "Nigeria" = "#A65628",
  "Mali" = "#eb88d4"
)

# Create separate layers for highlighted vs other countries
africa_highlighted <- africa %>% filter(name %in% highlight_countries)
africa_other <- africa %>% filter(!name %in% highlight_countries)


# Define East and West Africa bounds
xlim_range <- c(-20, 52)  # West to East
ylim_range <- c(-5, 20)   # South to North

# Plot map with site names and sample counts
# Plot map
ggplot() +
  # Non-highlighted countries in grey
  geom_sf(data = africa_other, fill = "grey90", color = "darkgrey", size = 0.3) +
  # Highlighted countries in white with colored borders
  geom_sf(data = africa_highlighted, aes(color = name), fill = "white", linewidth = 0.8) +
  scale_color_manual(values = country_colors, name = "Countries") +
  # Ethiopia admin regions
  geom_sf(data = ethiopia_admin1_sf, fill = NA, color = "grey", size = 0.3) +
  # Create a new color scale for the points (so it doesn't conflict with country borders)
  ggnewscale::new_scale_color() +
  # All non-main sites in grey
  geom_point(data = all_sites %>% filter(!main_site), 
             aes(x = lon, y = lat, size = sample_count), 
             color = "grey50", alpha = 0.8) +
  # Main 5 Ethiopian sites highlighted with original colors
  geom_point(data = all_sites %>% filter(main_site), 
             aes(x = lon, y = lat, size = sample_count, color = District)) +
  scale_color_discrete(name = "Study sites") + 
  scale_size_continuous(
    name = "Sample count",      
    range = c(2, 8),
    breaks = c(10, 20, 30, 40, 100, 200),
    limits = c(1, NA)
  ) +
  guides(color = guide_legend(override.aes = list(size = 6))) +
  coord_sf(xlim = xlim_range, ylim = ylim_range) +
  theme_void() +
  theme(
    legend.position = "right"
  )


ggsave("figs/map_all_samples.pdf", width = 13, height = 6, dpi = 300)

