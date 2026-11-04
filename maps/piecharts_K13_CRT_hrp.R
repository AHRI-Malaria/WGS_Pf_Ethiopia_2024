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

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("~/Google Drive/My Drive/Research_Fellow/Popgen/WGS_Pf_Ethiopia_2024/maps")

# --- Load map data
world <- ne_countries(scale = "medium", returnclass = "sf")
ethiopia_admin1 <- geodata::gadm("GADM", country = "ETH", level = 1) |> st_as_sf()


# Create long format dataframe
# python kelch13_convert.py
R622I <-read.csv("../input_files/R622I_longformat.csv")
met <- read.csv("../input_files/Popgen_db_samples_072025_filt_FINAL.csv")

met_selected <- met %>%
  dplyr::select(wgs_id, Country,District, Year)

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
  dplyr::select(District, District.longitude, District.latitude) %>%
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
  dplyr::select(District, District.longitude, District.latitude, Genotype, Percentage, Year, total_samples) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Percentage, values_fill = 0)

# --- Define colors
colors_genotype <- c("0/0" = "#5c965c", "0/1" = "#ff7f0e", "1/1" = "#1f77b4")
year_colors <- c("2013" = "#9b61ba", "2015" = "#6873b0", "2017" = "grey", "2018" = "#c4897e", "2019" = "#e8a04a", "2024" = "#d4ce7b")

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


# # # # # # Combined R622I + A675V with Eritrea
R622I_eri <- read.csv("../input_files/R622I_longformat.csv")
A675V_eri <- read.csv("../input_files/A675V_longformat.csv")
met_eth <- read.csv("../input_files/Popgen_db_samples_072025_filt_FINAL.csv")
met_eri <- read.csv("../input_files/Popgen_db_samples_012026_filt_FINAL_fws_Eritrea_filtered.csv")

# Prepare Ethiopia metadata
met_eth_sel <- met_eth %>%
  filter(Country == "Ethiopia") %>%
  dplyr::select(wgs_id, Country, District, District.latitude, District.longitude, Year)

# Prepare Eritrea metadata: use Admin.level.1 as District
met_eri_sel <- met_eri %>%
  filter(Country == "Eritrea") %>%
  dplyr::select(wgs_id, Country, District = Admin.level.1,
         District.latitude = Admin.level.1.latitude,
         District.longitude = Admin.level.1.longitude, Year)

# Combine metadata
met_combined <- bind_rows(met_eth_sel, met_eri_sel)

# Normalize genotype notation (pipe to slash)
normalize_genotype <- function(df) {
  df %>% mutate(Genotype = case_when(
    Genotype == "1|1" ~ "1/1",
    Genotype == "0|1" ~ "0/1",
    Genotype == "1|0" ~ "0/1",
    Genotype == "0|0" ~ "0/0",
    TRUE ~ Genotype
  ))
}
R622I_eri <- normalize_genotype(R622I_eri)
A675V_eri <- normalize_genotype(A675V_eri)

# Deduplicate: if a sample has multiple rows (multiple positions), keep most severe genotype
genotype_rank <- function(g) case_when(g == "1/1" ~ 3, g == "0/1" ~ 2, g == "0/0" ~ 1, TRUE ~ 0)

R622I_eri_dedup <- R622I_eri %>%
  mutate(rank = genotype_rank(Genotype)) %>%
  group_by(Sample) %>%
  slice_max(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Sample, Genotype)

A675V_eri_dedup <- A675V_eri %>%
  mutate(rank = genotype_rank(Genotype)) %>%
  group_by(Sample) %>%
  slice_max(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(Sample, Genotype)

R622I_eri_met <- left_join(R622I_eri_dedup, met_combined, by = c("Sample" = "wgs_id"))
A675V_eri_met <- left_join(A675V_eri_dedup, met_combined, by = c("Sample" = "wgs_id"))

# Classify R622I
R622I_eri_class <- R622I_eri_met %>%
  filter(Country %in% c("Ethiopia", "Eritrea")) %>%
  mutate(Classification = case_when(
    Genotype == "1/1" ~ "R622I",
    Genotype == "0/1" ~ "WT+R622I",
    Genotype == "0/0" ~ "WT_R622I"
  )) %>%
  dplyr::select(Sample, District, Year, R622I_status = Classification)

# Classify A675V
A675V_eri_class <- A675V_eri_met %>%
  filter(Country %in% c("Ethiopia", "Eritrea")) %>%
  mutate(Classification = case_when(
    Genotype == "1/1" ~ "A675V",
    Genotype == "0/0" ~ "WT_A675V"
  )) %>%
  dplyr::select(Sample, A675V_status = Classification)

# Combine classifications
combined_eri_class <- R622I_eri_class %>%
  left_join(A675V_eri_class, by = "Sample") %>%
  mutate(Final_Classification = case_when(
    R622I_status == "R622I" & A675V_status == "A675V" ~ "R622I+A675V",
    R622I_status == "R622I" ~ "R622I",
    R622I_status == "WT+R622I" ~ "Mixed R622I + WT",
    A675V_status == "A675V" ~ "A675V",
    TRUE ~ "WT"
  )) %>%
  filter(!is.na(District) & District != "NA")

combined_eri_class <- combined_eri_class %>%
  mutate(Final_Classification = factor(Final_Classification,
                                       levels = c("WT", "R622I", "Mixed R622I + WT", "A675V")))

# Count and calculate percentages
combined_eri_counts <- combined_eri_class %>%
  group_by(District, Final_Classification, Year) %>%
  summarise(Count = n(), .groups = "drop")

combined_eri_pct <- combined_eri_counts %>%
  group_by(District, Year) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

districts_combined <- met_combined %>%
  dplyr::select(District, District.longitude, District.latitude) %>%
  distinct() %>%
  na.omit()

combined_eri_map <- combined_eri_pct %>%
  left_join(districts_combined, by = "District")

sample_totals_eri <- combined_eri_map %>%
  group_by(District) %>%
  summarise(total_samples = sum(Count), .groups = "drop")

combined_eri_map <- combined_eri_map %>%
  left_join(sample_totals_eri, by = "District")

pie_data_eri <- combined_eri_map %>%
  dplyr::select(District, District.longitude, District.latitude, Final_Classification, Percentage, Year, total_samples) %>%
  tidyr::pivot_wider(names_from = Final_Classification, values_from = Percentage, values_fill = 0)

colors_combined_eri <- c(
  "WT" = "#6cd669",
  "R622I" = "red",
  "Mixed R622I + WT" = "#ff7f0e",
  "A675V" = "yellow"
)

scale_pie_size <- function(n) sqrt(n) * 0.07
pie_data_eri <- pie_data_eri %>%
  mutate(radius = scale_pie_size(total_samples)) %>%
  mutate(Year = ifelse(Year == "2017/2018", 2017.5, as.numeric(Year)))

# Load Eritrea admin boundaries
eritrea_admin1 <- geodata::gadm("GADM", country = "ERI", level = 1) |> st_as_sf()

# Country outlines for thick borders
eth_eri_outline <- world %>% filter(name %in% c("Ethiopia", "Eritrea"))

# Create map with both countries
base_map_eri <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "white") +
  geom_sf(data = ethiopia_admin1, fill = "white", color = "black", linewidth = 0.2) +
  geom_sf(data = eritrea_admin1, fill = "white", color = "black", linewidth = 0.2) +
  geom_sf(data = eth_eri_outline, fill = NA, color = "black", linewidth = 1.2) +
  coord_sf(xlim = c(32, 48), ylim = c(3, 18), expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  geom_circle(
    data = pie_data_eri,
    aes(x0 = District.longitude, y0 = District.latitude, r = radius * 1.1, color = Year),
    linewidth = 2, fill = NA
  ) +
  geom_scatterpie(
    data = pie_data_eri,
    aes(x = District.longitude, y = District.latitude, r = radius),
    cols = names(colors_combined_eri)[names(colors_combined_eri) %in% names(pie_data_eri)],
    color = "white",
    size = 0.2,
    alpha = 1
  ) +
  scale_fill_manual(values = colors_combined_eri, name = "PfK13 Genotype") +
  scale_color_gradient(
    low = "#98d8ed",
    high = "#09449c",
    name = "Collection year",
    breaks = c(2013, 2016, 2019, 2022, 2024)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    color = guide_colorbar(order = 2, barwidth = 1, barheight = 8)
  )

base_map_eri

ggsave("figs/map_kelch13_with_Eritrea.pdf", base_map_eri, width = 8, height = 6, dpi = 300)


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
year_colors <- c("2013" = "#9b61ba", "2015" = "#6873b0", "2017" = "grey", "2018" = "#c4897e", "2019" = "#e8a04a", "2024" = "#d4ce7b")

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
# Merge both mutation datasets
R622I_classification <- R622I_met %>%
  filter(Country == "Ethiopia") %>%
  mutate(Classification = case_when(
    Genotype %in% c("1/1") ~ "R622I",
    Genotype %in% c("0/1") ~ "WT+R622I",
    Genotype == "0/0" ~ "WT_R622I"
  )) %>%
  select(Sample, District, Year, R622I_status = Classification)

A675V_classification <- A675V_met %>%
  filter(Country == "Ethiopia") %>%
  mutate(Classification = case_when(
    Genotype == "1/1" ~ "A675V",
    Genotype == "0/0" ~ "WT_A675V"
  )) %>%
  select(Sample, District, Year, A675V_status = Classification)

# Combine classifications
combined_classification <- R622I_classification %>%
  left_join(A675V_classification, by = c("Sample", "District", "Year")) %>%
  mutate(Final_Classification = case_when(
    R622I_status == "R622I" & A675V_status == "A675V" ~ "R622I+A675V",
    R622I_status == "R622I" ~ "R622I",
    R622I_status == "WT+R622I" ~ "Mixed R622I + WT",
    A675V_status == "A675V" ~ "A675V",
    TRUE ~ "WT"
  )) %>%
  filter(District != "NA")

combined_classification <- combined_classification %>%
  mutate(Final_Classification = factor(Final_Classification, 
                                       levels = c("WT", "R622I", "Mixed R622I + WT", "A675V")))

# Count samples by classification
combined_counts <- combined_classification %>%
  group_by(District, Final_Classification, Year) %>%
  summarise(Count = n(), .groups = "drop")

# Calculate percentages
combined_percentage <- combined_counts %>%
  group_by(District, Year) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# Join with district coordinates
combined_map_data <- combined_percentage %>%
  left_join(districts, by = "District")

# Get sample totals
sample_totals_combined <- combined_map_data %>%
  group_by(District) %>%
  summarise(total_samples = sum(Count), .groups = "drop")

combined_map_data <- combined_map_data %>%
  left_join(sample_totals_combined, by = "District")

# Pivot to wide format for scatterpie
pie_data_combined <- combined_map_data %>%
  select(District, District.longitude, District.latitude, Final_Classification, Percentage, Year, total_samples) %>%
  tidyr::pivot_wider(names_from = Final_Classification, values_from = Percentage, values_fill = 0)

# Define colors for combined classification
colors_combined <- c(
  "WT" = "#6cd669", 
  "R622I" = "red",
  "Mixed R622I + WT" = "#ff7f0e",
  "A675V" = "yellow"
)

# Scale pie size
pie_data_combined <- pie_data_combined %>%
  mutate(radius = scale_pie_size(total_samples)) %>%
  mutate(Year = ifelse(Year == "2017/2018", 2017.5, as.numeric(Year)))

# Create map
base_map_combined <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "white") +
  geom_sf(data = ethiopia_admin1, fill = "white", color = "black", size = 0.3) +
  coord_sf(xlim = c(32, 48), ylim = c(3, 15), expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  # Add outer ring (colored by Year)
  geom_circle(
    data = pie_data_combined,
    aes(x0 = District.longitude, y0 = District.latitude, r = radius * 1.1, color = Year),
    linewidth = 2,  # Changed to thinner
    fill = NA
  ) +
  # Add pies
  geom_scatterpie(
    data = pie_data_combined,
    aes(x = District.longitude, y = District.latitude, r = radius),
    cols = names(colors_combined)[names(colors_combined) %in% names(pie_data_combined)],
    color = "white",
    size = 0.2,  # Thin white borders between slices
    alpha = 1
  ) +
  scale_fill_manual(values = colors_combined, name = "PfK13 Genotype") +
  scale_color_gradient(
    low = "#98d8ed",
    high = "#09449c",
    name = "Collection year",
    breaks = c(2013, 2016, 2019, 2022, 2024)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    color = guide_colorbar(order = 2, barwidth = 1, barheight = 8)  # horizontal bar
  )

base_map_combined


ggsave("figs/map_kelch13.pdf",base_map_combined, width = 8, height = 6, dpi = 300)





# CRT K76T
# python crt_convert.py
K76T <-read.csv("../input_files/K76T_longformat.csv")
met <- read.csv("../input_files/Popgen_db_samples_012026_filt_FINAL_fws.csv")

K76T <- K76T %>%
  mutate(Genotype = case_when(
    Genotype == "1|1" ~ "1/1",
    Genotype == "0|1" ~ "0/1",
    Genotype == "0|0" ~ "0/0",
    TRUE ~ Genotype 
  ))

met_selected <- met %>%
  select(wgs_id, Country,District, Year)

K76T_met <- left_join(K76T, met_selected, by = c("Sample" = "wgs_id"))

K76T_counts <- K76T_met %>%
  filter(Country=="Ethiopia") %>%
  group_by(District, Genotype, Year) %>%
  summarise(Count = n(), .groups = "drop")%>%
  filter(District != "NA")

K76T_percentage <- K76T_counts %>%
  group_by(District, Year) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

districts<-met %>%
  filter(Country == "Ethiopia") %>%
  select(District, Admin.level.1.longitude, Admin.level.1.latitude) %>%
  distinct() %>%
  na.omit()

genotype_map_data <- K76T_percentage %>%
  left_join(districts, by = c("District" = "District"))


sample_totals <- genotype_map_data %>%
  group_by(District) %>%
  summarise(total_samples = sum(Count), .groups = "drop")

# Add total sample info back to the data
genotype_map_data <- genotype_map_data %>%
  left_join(sample_totals, by = "District")


# Pivot to wide format for scatterpie
pie_data <- genotype_map_data %>%
  select(District, Admin.level.1.longitude, Admin.level.1.latitude, Genotype, Percentage, Year, total_samples) %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = Percentage, values_fill = 0)

# --- Define colors
colors_genotype <- c("0/0" = "#6cd669", "0/1" = "#ff7f0e", "1/1" = "red")

# --- Define scaling for pie radius
scale_pie_size <- function(n) sqrt(n) * 0.09
pie_data <- pie_data %>%
  mutate(radius = scale_pie_size(total_samples))%>%
  mutate(Year = ifelse(Year == "2017/2018", 2017.5, as.numeric(Year)))

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
    aes(x0 = Admin.level.1.longitude, y0 = Admin.level.1.latitude, r = radius * 1.1, color = Year),  # slightly larger than pie
    linewidth = 3, fill = NA
  )

# --- Add pies
base_map <- base_map +
  geom_scatterpie(
    data = pie_data,
    aes(x = Admin.level.1.longitude, y = Admin.level.1.latitude, r = radius),
    cols = names(colors_genotype),
    color = "white", # border between slices
    alpha = 1
  ) +
  scale_fill_manual(values = colors_genotype, name = "PfCRT K76T") +
  scale_color_gradient(
    low = "#98d8ed",
    high = "#09449c",
    name = "Collection year",
    breaks = c(2013, 2016, 2019, 2022, 2024)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    color = guide_colorbar(order = 2, barwidth = 1, barheight = 8)  # horizontal bar
  )

pie_legend <- data.frame(
  Admin.level.1.longitude = c(48.5, 48.5),   # place just outside map
  Admin.level.1.latitude = c(13, 10),
  total_samples = c(20, 30)
) %>%
  mutate(radius = scale_pie_size(total_samples))

# --- Add legend for pie size
base_map +
  geom_circle(
    data = pie_legend,
    aes(x0 = Admin.level.1.longitude, y0 = Admin.level.1.latitude, r = radius),
    fill = "grey90", color = "grey50"
  ) 


base_map

ggsave("figs/map_K76T.pdf",base_map, width = 8, height = 6, dpi = 300)


# --- K76T frequency 2024: Gambella vs other districts
# 76T present = 0/1 or 1/1 (mixed infections counted as carrying the mutant)
K76T_2024 <- K76T_met %>%
  filter(Country == "Ethiopia", Year == "2024",
         District %in% c("Gambella", "Arba Minch", "Mizan", "Metehera Fentale", "Gondar")) %>%
  mutate(K76T_present = Genotype %in% c("0/1", "1/1"),
         Group = ifelse(District == "Gambella", "Gambella", "Other"))

# Per district
K76T_2024 %>%
  group_by(District) %>%
  summarise(n = n(), n_76T = sum(K76T_present), freq_76T = round(100 * mean(K76T_present), 1))

# Gambella vs other districts
K76T_2024 %>%
  group_by(Group) %>%
  summarise(n = n(), n_76T = sum(K76T_present), freq_76T = round(100 * mean(K76T_present), 1))

K76T_tab <- table(K76T_2024$Group, K76T_2024$K76T_present)
K76T_tab
fisher.test(K76T_tab)







# HRP
hrp <-read.csv("../input_files/sample_meta_with_hrp_status_updated.csv")
met <- read.csv("../input_files/Popgen_db_samples_072025_filt_FINAL.csv")

# Join hrp with metadata
hrp_met <- left_join(hrp, met_selected, by = c("Sample" = "wgs_id"))

# If site column already has district info, you may use that instead
# Classify hrp_status
hrp_classification <- hrp_met %>%
  filter(Country == "Ethiopia") %>%
  mutate(HRP_status = case_when(
    hrp_status == "wild" ~ "WT",
    hrp_status == "hrp2_only" ~ "hrp2 deletion",
    hrp_status == "hrp3_only" ~ "hrp3 deletion",
    hrp_status == "mutant" ~ "hrp2+hrp3 deletion",
    TRUE ~ hrp_status
  )) %>%
  filter(!is.na(District) & District != "NA")

# Count and calculate percentages
hrp_counts <- hrp_classification %>%
  group_by(District, HRP_status, Year) %>%
  summarise(Count = n(), .groups = "drop")

hrp_percentage <- hrp_counts %>%
  group_by(District, Year) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

# Join with district coordinates
hrp_map_data <- hrp_percentage %>%
  left_join(districts, by = "District")


# Get sample totals
sample_totals_hrp <- hrp_map_data %>%
  group_by(District) %>%
  summarise(total_samples = sum(Count), .groups = "drop")

hrp_map_data <- hrp_map_data %>%
  left_join(sample_totals_hrp, by = "District")

# Pivot to wide format for scatterpie
pie_data_hrp <- hrp_map_data %>%
  select(District, District.longitude, District.latitude, HRP_status, Percentage, Year, total_samples) %>%
  tidyr::pivot_wider(names_from = HRP_status, values_from = Percentage, values_fill = 0)

# Define colors
colors_hrp <- c(
  "WT"                 = "#6cd669",
  "hrp2 deletion"      = "#ff7f0e",
  "hrp3 deletion"      = "yellow",
  "hrp2+hrp3 deletion" = "red"
)

# Scale and clean year
pie_data_hrp <- pie_data_hrp %>%
  mutate(radius = scale_pie_size(total_samples)) %>%
  mutate(Year = ifelse(Year == "2017/2018", 2017.5, as.numeric(Year)))

# Create map
base_map_hrp <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "white") +
  geom_sf(data = ethiopia_admin1, fill = "white", color = "black", size = 0.3) +
  coord_sf(xlim = c(32, 48), ylim = c(3, 15), expand = FALSE) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  geom_circle(
  data = pie_data_hrp,
  aes(x0 = District.longitude, y0 = District.latitude, r = radius * 1.1, color = Year),
  linewidth = 2, fill = NA
  ) +
  geom_scatterpie(
    data = pie_data_hrp,
    aes(x = District.longitude, y = District.latitude, r = radius),
    cols = names(colors_hrp)[names(colors_hrp) %in% names(pie_data_hrp)],
    color = "white",
    size = 0.2,
    alpha = 1
  ) +
  scale_fill_manual(values = colors_hrp, name = "HRP deletion status") +
  scale_color_continuous(
    "#09449c",
    name = "Collection year",
    breaks = c(2024)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 5)),
    color = guide_colorbar(order = 2, barwidth = 1, barheight = 8)
  )

base_map_hrp


ggsave("figs/map_hrp.pdf",base_map_hrp, width = 8, height = 6, dpi = 300)
