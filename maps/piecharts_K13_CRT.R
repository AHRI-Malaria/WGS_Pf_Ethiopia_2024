library(dplyr)
library(ggplot2)
library(sf)
library(scatterpie)
library(ggforce)
library(RColorBrewer)
library(geodata)
library(rnaturalearth)
library(rnaturalearthdata)
library(readxl)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Load map data
world <- ne_countries(scale = "medium", returnclass = "sf")
ethiopia_admin1 <- geodata::gadm("GADM", country = "ETH", level = 1) |> st_as_sf()


# Create long format datatframe
# python kelch13_convert.py
R622I <-read.csv("../input_files/R622I_longformat.csv")
met <- read.csv("../input_files/Popgen_db_samples_072025_filt_FINAL.csv")

met_selected <- met %>%
  select(wgs_id, Country,District, Year)

R622I_met <- left_join(R622I, met_selected, by = c("Sample" = "wgs_id"))

R622I_counts <- R622I_met %>%
  filter(Country=="Ethiopia") %>%
  group_by(District, Genotype, Year) %>%
  summarise(Count = n(), .groups = "drop")%>%
  filter(District != "NA")

R622I_percentage <- R622I_counts %>%
  group_by(District, Year) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

districts<-met %>%
  filter(Country == "Ethiopia") %>%
  select(District, District.longitude, District.latitude) %>%
  distinct() %>%
  na.omit()

genotype_map_data <- R622I_percentage %>%
  left_join(districts, by = c("District" = "District"))


sample_totals <- genotype_map_data %>%
  group_by(District) %>%
  summarise(total_samples = sum(Count), .groups = "drop")

# Add total sample info back to the data
genotype_map_data <- genotype_map_data %>%
  left_join(sample_totals, by = "District")


# Pivot to wide format for scatterpie
pie_data <- genotype_map_data %>%
  select(District, District.longitude, District.latitude, Genotype, Percentage, Year, total_samples) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Percentage, values_fill = 0)

# --- Define colors
colors_genotype <- c("0/0" = "#5c965c", "0/1" = "#ff7f0e", "1/1" = "#1f77b4")
year_colors <- c("2013" = "#9b61ba", "2015" = "#6873b0", "2017" = "grey", "2018" = "#c4897e", "2024" = "#d4ce7b")

# --- Define scaling for pie radius
scale_pie_size <- function(n) sqrt(n) * 0.09
pie_data <- pie_data %>%
  mutate(radius = scale_pie_size(total_samples))

# --- Create map
base_map <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "white") +
  geom_sf(data = ethiopia_admin1, fill = "white", color = "black", size = 0.3) +
  coord_sf(xlim = c(32, 48), ylim = c(3, 15), expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

# --- Add outer ring (colored by Year)
base_map <- base_map +
  geom_circle(
    data = pie_data,
    aes(x0 = District.longitude, y0 = District.latitude, r = radius * 1.1, color = as.factor(Year)),  # slightly larger than pie
    linewidth = 3, fill = NA
  )

# --- Add pies
base_map <- base_map +
  geom_scatterpie(
    data = pie_data,
    aes(x = District.longitude, y = District.latitude, r = radius),
    cols = names(colors_genotype),
    color = "white", # border between slices
    alpha = 1
  ) +
  scale_fill_manual(values = colors_genotype, name = "PfK13 R622I") +
  scale_color_manual(values = year_colors, name = "Year") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    color = guide_legend(order = 2),
    size = guide_legend(order = 3)
  )

pie_legend <- data.frame(
  District.longitude = c(48.5, 48.5),   # place just outside map
  District.latitude = c(13, 10),
  total_samples = c(20, 30)
) %>%
  mutate(radius = scale_pie_size(total_samples))

# --- Add legend for pie size
base_map +
  geom_circle(
    data = pie_legend,
    aes(x0 = District.longitude, y0 = District.latitude, r = radius),
    fill = "grey90", color = "grey50"
  ) +
  geom_text(
    data = pie_legend,
    aes(x = District.longitude + radius + 0.3, y = District.latitude, label = total_samples),
    hjust = 0, size = 3
  ) +
  annotate("text", x = 10, y = 14, label = "Sample count", fontface = "bold", size = 3)

# --- Final map
map<-base_map+
  geom_text(data = districts, aes(x = District.longitude, y = District.latitude, label = District), size = 3, hjust = -0.5)

map

ggsave("figs/map_R622I.pdf",base_map, width = 8, height = 6, dpi = 300)




# # # # # # A675V
A675V <-read.csv("../input_files/A675V_longformat.csv")
met <- read.csv("../input_files/Popgen_db_samples_072025_filt_FINAL.csv")

met_selected <- met %>%
  select(wgs_id, Country,District, Year)

A675V_met <- left_join(A675V, met_selected, by = c("Sample" = "wgs_id"))

A675V_counts <- A675V_met %>%
  filter(Country=="Ethiopia") %>%
  group_by(District, Genotype, Year) %>%
  summarise(Count = n(), .groups = "drop")%>%
  filter(District != "NA")


A675V_percentage <- A675V_counts %>%
  group_by(District, Year) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

districts<-met %>%
  filter(Country == "Ethiopia") %>%
  select(District, District.longitude, District.latitude) %>%
  distinct() %>%
  na.omit()

genotype_map_data <- A675V_percentage %>%
  left_join(districts, by = c("District" = "District"))

sample_totals <- genotype_map_data %>%
  group_by(District) %>%
  summarise(total_samples = sum(Count), .groups = "drop")

# Add total sample info back to the data
genotype_map_data <- genotype_map_data %>%
  left_join(sample_totals, by = "District")


# Pivot to wide format for scatterpie
pie_data <- genotype_map_data %>%
  select(District, District.longitude, District.latitude, Genotype, Percentage, Year, total_samples) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Percentage, values_fill = 0)

# --- Define colors
colors_genotype <- c("0/0" = "#5c965c", "1/1" = "#1f77b4")
year_colors <- c("2013" = "#9b61ba", "2015" = "#6873b0", "2017" = "grey", "2018" = "#c4897e", "2024" = "#d4ce7b")

# --- Define scaling for pie radius
scale_pie_size <- function(n) sqrt(n) * 0.09
pie_data <- pie_data %>%
  mutate(radius = scale_pie_size(total_samples))

# --- Create map
base_map <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "white") +
  geom_sf(data = ethiopia_admin1, fill = "white", color = "black", size = 0.3) +
  coord_sf(xlim = c(32, 48), ylim = c(3, 15), expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

# --- Add outer ring (colored by Year)
base_map <- base_map +
  geom_circle(
    data = pie_data,
    aes(x0 = District.longitude, y0 = District.latitude, r = radius * 1.1, color = as.factor(Year)),  # slightly larger than pie
    linewidth = 3, fill = NA
  )

# --- Add pies
base_map <- base_map +
  geom_scatterpie(
    data = pie_data,
    aes(x = District.longitude, y = District.latitude, r = radius),
    cols = names(colors_genotype),
    color = "white", # border between slices
    alpha = 1
  ) +
  scale_fill_manual(values = colors_genotype, name = "PfK13 A675V") +
  scale_color_manual(values = year_colors, name = "Year") +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    color = guide_legend(order = 2),
    size = guide_legend(order = 3)
  )

pie_legend <- data.frame(
  District.longitude = c(48.5, 48.5),   # place just outside map
  District.latitude = c(13, 10),
  total_samples = c(20, 30)
) %>%
  mutate(radius = scale_pie_size(total_samples))

# --- Add legend for pie size
base_map +
  geom_circle(
    data = pie_legend,
    aes(x0 = District.longitude, y0 = District.latitude, r = radius),
    fill = "grey90", color = "grey50"
  ) +
  geom_text(
    data = pie_legend,
    aes(x = District.longitude + radius + 0.3, y = District.latitude, label = total_samples),
    hjust = 0, size = 3
  ) +
  annotate("text", x = 10, y = 14, label = "Sample count", fontface = "bold", size = 3)

# --- Final map
map<-base_map+
  geom_text(data = districts, aes(x = District.longitude, y = District.latitude, label = District), size = 3, hjust = -0.5)

map

ggsave("figs/map_A675V.pdf",base_map, width = 8, height = 6, dpi = 300)


# # # # # # Combined


