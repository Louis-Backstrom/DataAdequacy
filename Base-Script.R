#### 0. Initialization ---------------------------------------------------------

library(tidyverse)
library(sf)
library(auk)
library(readxl)
library(lubridate)
library(ggpubr)
library(units)
library(performance)

source("Data-Paths.R") # this will read in all the directories where the datasets are stored locally (not on this repo)

`%notin%` <- Negate(`%in%`)
sf_use_s2(TRUE)

#### 1. Set Up -----------------------------------------------------------------

taxonomy_fixer <- read_xlsx("Data/TaxonomyFixerMaster.xlsx", sheet = 2)
wlab <- read_xlsx(wlab_v4pt2_path)
aus <- read_sf(aus_path) %>% 
  st_make_valid()
trait <- read_xlsx("Data/TraitMaster.xlsx")

# cell_size <- 1 # now set by the variable cell_size script

all_species <- taxonomy_fixer %>% 
  filter(CLASS %in% c("Island", "Normal", "Shorebird")) %>% 
  pull(OUT) %>% 
  unique()

grd <- st_make_grid(round(st_bbox(aus), 0), cellsize = cell_size, crs = 4326) %>% 
  st_sf() %>% 
  rowid_to_column(var = "cell_id") %>% 
  mutate(full_area = as.numeric(st_area(.)))

grd_clip <- st_intersection(grd, aus) %>% 
  st_make_valid() %>% 
  mutate(clip_area = as.numeric(st_area(.))) %>% 
  st_drop_geometry()

grd <- left_join(grd, grd_clip, by = c("cell_id", "full_area")) %>% 
  mutate(prop_land = replace_na(clip_area / full_area, 0)) %>% 
  select(-full_area, -clip_area, -cell_id) %>% 
  filter(prop_land > 0) %>% 
  st_make_valid() %>% 
  rowid_to_column(var = "cell_id")

rm(grd_clip)


##### 1.1. eBird ----------------------------------------------------------------

if (file.exists(paste0("Data/", cell_size, "/ebird_cleaned.txt")) == TRUE) {
  rebuild_ebird <- askYesNo(msg = "Do you want to rebuild the eBird Dataset?")
} else {
  rebuild_ebird <- TRUE
}

if (rebuild_ebird == TRUE) {
  
  ebird_default <- read_ebd(ebd_path) %>% 
    select("checklist_id", "all_species_reported", "scientific_name", "observation_date", "latitude", "longitude", "exotic_code") %>% 
    filter(exotic_code %notin% c("X"),
           all_species_reported == TRUE) %>%
    left_join(., taxonomy_fixer, by = c("scientific_name" = "IN")) %>% 
    rename("species_binomial" = "OUT") %>% 
    filter(CLASS %in% c("Island", "Normal", "Shorebird")) %>% 
    select(-all_species_reported, -scientific_name, -exotic_code, -CLASS, -WLAB, -Population)
  
  ebird_sensitive <- read_ebd(ebd_sensitive_path) %>% 
    select("checklist_id", "all_species_reported", "scientific_name", "observation_date", "latitude", "longitude", "exotic_code") %>% 
    filter(exotic_code %notin% c("X"),
           all_species_reported == TRUE) %>%
    left_join(., taxonomy_fixer, by = c("scientific_name" = "IN")) %>% 
    rename("species_binomial" = "OUT") %>% 
    filter(CLASS %in% c("Island", "Normal", "Shorebird")) %>% 
    select(-all_species_reported, -scientific_name, -exotic_code, -CLASS, -WLAB, -Population)
  
  ebird <- bind_rows(ebird_default, ebird_sensitive) %>% 
    mutate(year = year(observation_date))
  
  loc <- ebird %>% 
    select(latitude, longitude) %>% 
    distinct() %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = F) %>% 
    mutate(in_aus = as.logical(lengths(st_intersects(., aus)))) %>% 
    st_drop_geometry()
  
  ebird <- ebird %>% 
    left_join(., loc) %>% 
    filter(in_aus == TRUE) %>% 
    select(-in_aus)
  
  rm(ebird_default, ebird_sensitive, loc)
  
  write_tsv(ebird, paste0("Data/", cell_size, "/ebird_cleaned.txt"))
  
} else {
  
  ebird <- read_tsv(paste0("Data/", cell_size, "/ebird_cleaned.txt"))
  
} 

##### 1.2. Birdata --------------------------------------------------------------

if (file.exists(paste0("Data/", cell_size, "/birdata_cleaned.txt")) == TRUE) {
  rebuild_birdata <- askYesNo(msg = "Do you want to rebuild the Birdata Dataset?")
} else {
  rebuild_birdata <- TRUE
}

if (rebuild_birdata == TRUE) {
  
  birdata <- read_csv(birdata_path) %>% 
    select("checklist_id" = "Survey.ID", "all_species_reported" = "All.Species.Recorded",
           "scientific_name" = "Scientific.Name", "observation_date" = "Start.Date",
           "latitude" = "Latitude", "longitude" = "Longitude") %>% 
    filter(all_species_reported == "Yes") %>% 
    left_join(., taxonomy_fixer, by = c("scientific_name" = "IN")) %>% 
    rename("species_binomial" = "OUT") %>% 
    filter(CLASS %in% c("Island", "Normal", "Shorebird")) %>% 
    select(-all_species_reported, -scientific_name, -CLASS, -WLAB, -Population) %>% 
    mutate(year = year(observation_date))
  
  loc <- birdata %>% 
    select(latitude, longitude) %>% 
    distinct() %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = F) %>% 
    mutate(in_aus = as.logical(lengths(st_intersects(., aus)))) %>% 
    st_drop_geometry()
  
  birdata <- birdata %>% 
    left_join(., loc) %>% 
    filter(in_aus == TRUE) %>% 
    select(-in_aus)
  
  rm(loc)
  
  write_tsv(birdata, paste0("Data/", cell_size, "/birdata_cleaned.txt"))
  
} else {
  
  birdata <- read_tsv(paste0("Data/", cell_size, "/birdata_cleaned.txt"))
  
} 

##### 1.3. ABG ------------------------------------------------------------------

if (file.exists(paste0("Data/", cell_size, "/abg_poly_cleaned.gpkg")) == TRUE) {
  rebuild_abg <- askYesNo(msg = "Do you want to rebuild the ABG Maps?")
} else {
  rebuild_abg <- TRUE
}

if (rebuild_abg == TRUE) {
  
  sf_use_s2(FALSE)
  
  abg <- read_sf(abg_path) %>% 
    left_join(., filter(wlab, TaxonLevel == "sp"), by = c("sp_id" = "SpID")) %>%
    filter(rnge %in% c(1, 5, 6, 7)) %>% 
    left_join(., taxonomy_fixer, by = c("TaxonScientificName" = "IN")) %>% 
    rename("species_binomial" = "OUT") %>% 
    filter(CLASS %in% c("Island", "Normal", "Shorebird")) %>% 
    select(species_binomial) %>% 
    drop_na(species_binomial) %>% 
    group_by(species_binomial) %>% 
    summarise() %>% 
    st_cast() %>% 
    st_transform(crs = 4326) %>% 
    st_intersection(., aus) %>% 
    st_make_valid()
  
  abg_grd <- st_intersects(abg, grd) %>% 
    as.data.frame() %>% 
    rename("sp_id" = "row.id", "cell_id" = "col.id") %>%
    left_join(., st_drop_geometry(rowid_to_column(abg, var = "sp_id")), by = "sp_id") %>% 
    select(cell_id, species_binomial)
  
  abg_grd <- expand_grid(species_binomial = all_species,
                         cell_id = grd$cell_id) %>%
    left_join(., mutate(abg_grd, in_abg = TRUE)) %>% 
    mutate(in_abg = replace_na(in_abg, FALSE))
  
  abg_grd_sf <- full_join(grd, abg_grd, by = "cell_id") %>% 
    filter(in_abg == TRUE)
  
  sf_use_s2(TRUE)
  
  write_sf(abg, paste0("Data/", cell_size, "/abg_poly_cleaned.gpkg"))
  write_sf(abg_grd_sf, paste0("Data/", cell_size, "/abg_grid_cleaned.gpkg"))
  write_tsv(abg_grd, paste0("Data/", cell_size, "/abg_cleaned.gpkg"))
  
} else {
  
  abg <- read_sf(paste0("Data/", cell_size, "/abg_poly_cleaned.gpkg"))
  abg_grd_sf <- read_sf(paste0("Data/", cell_size, "/abg_grid_cleaned.gpkg"))
  abg_grd <- read_tsv(paste0("Data/", cell_size, "/abg_cleaned.gpkg"))
  
}

##### 1.4. BLI ------------------------------------------------------------------

if (file.exists(paste0("Data/", cell_size, "/bli_poly_cleaned.gpkg")) == TRUE) {
  rebuild_bli <- askYesNo(msg = "Do you want to rebuild the BLI Maps?")
} else {
  rebuild_bli <- TRUE
}

if (rebuild_bli == TRUE) {
  
  sf_use_s2(FALSE)
  
  bli <- read_sf(bli_path) %>% 
    filter(origin %in% c(1, 2, 3, 5),
           presence == 1) %>% 
    left_join(., taxonomy_fixer, by = c("binomial" = "IN")) %>% 
    rename("species_binomial" = "OUT") %>% 
    filter(CLASS %in% c("Island", "Normal", "Shorebird")) %>% 
    select(species_binomial) %>% 
    drop_na(species_binomial) %>% 
    group_by(species_binomial) %>% 
    summarise() %>% 
    st_cast() %>% 
    st_transform(crs = 4326) %>% 
    st_intersection(., aus) %>% 
    st_make_valid()
  
  bli_grd <- st_intersects(bli, grd) %>% 
    as.data.frame() %>% 
    rename("sp_id" = "row.id", "cell_id" = "col.id") %>%
    left_join(., st_drop_geometry(rowid_to_column(bli, var = "sp_id")), by = "sp_id") %>% 
    select(cell_id, species_binomial)
  
  bli_grd <- expand_grid(species_binomial = all_species,
                         cell_id = grd$cell_id) %>%
    left_join(., mutate(bli_grd, in_bli = TRUE)) %>% 
    mutate(in_bli = replace_na(in_bli, FALSE))
  
  bli_grd_sf <- full_join(grd, bli_grd, by = "cell_id") %>% 
    filter(in_bli == TRUE)
  
  sf_use_s2(TRUE)
  
  write_sf(bli, paste0("Data/", cell_size, "/bli_poly_cleaned.gpkg"))
  write_sf(bli_grd_sf, paste0("Data/", cell_size, "/bli_grid_cleaned.gpkg"))
  write_tsv(bli_grd, paste0("Data/", cell_size, "/bli_cleaned.gpkg"))
  
} else {
  
  bli <- read_sf(paste0("Data/", cell_size, "/bli_poly_cleaned.gpkg"))
  bli_grd_sf <- read_sf(paste0("Data/", cell_size, "/bli_grid_cleaned.gpkg"))
  bli_grd <- read_tsv(paste0("Data/", cell_size, "/bli_cleaned.gpkg"))
  
}

#### 2. Preparation ------------------------------------------------------------

add_grd <- left_join(abg_grd, bli_grd) %>% 
  mutate(in_either = in_bli + in_abg != 0,
         in_both = in_bli + in_abg == 2) %>% 
  left_join(., grd) %>% 
  select(cell_id, prop_land, species_binomial, in_abg, in_bli, in_either, in_both)

ebird_loc_grd <- ebird %>% 
  select(latitude, longitude) %>% 
  distinct() %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>% 
  st_join(., grd) %>% 
  st_drop_geometry() %>% 
  select(latitude, longitude, cell_id)

ebird <- left_join(ebird, ebird_loc_grd, by = c("latitude", "longitude")) %>% 
  mutate(year = if_else(year < 1989, 1989, year)) %>% 
  filter(year < 2023)

ebird_grd <- ebird %>% 
  group_by(cell_id, species_binomial, year) %>% 
  summarise(n_obs.ebird = n())

ebird_year_grd <- expand_grid(
  cell_id = unique(add_grd$cell_id),
  year = range(ebird$year)[1]:range(ebird$year)[2],
  species_binomial = all_species)

ebird_checklist_grd <- ebird %>% 
  group_by(cell_id, year) %>% 
  summarise(n_checklists.ebird = n_distinct(checklist_id))

ebird_full_grd <- left_join(ebird_year_grd, ebird_grd, by = c("cell_id", "year", "species_binomial")) %>% 
  mutate(n_obs.ebird = replace_na(n_obs.ebird, 0)) %>% 
  left_join(., ebird_checklist_grd, by = c("cell_id", "year")) %>% 
  mutate(n_checklists.ebird = replace_na(n_checklists.ebird, 0))

birdata_loc_grd <- birdata %>% 
  select(latitude, longitude) %>% 
  distinct() %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE) %>% 
  st_join(., grd) %>% 
  st_drop_geometry() %>% 
  select(latitude, longitude, cell_id)

birdata <- left_join(birdata, birdata_loc_grd, by = c("latitude", "longitude")) %>% 
  mutate(year = if_else(year < 1989, 1989, year)) %>% 
  filter(year < 2023)

birdata_grd <- birdata %>% 
  group_by(cell_id, species_binomial, year) %>% 
  summarise(n_obs.birdata = n())

birdata_year_grd <- expand_grid(
  cell_id = unique(add_grd$cell_id),
  year = range(birdata$year)[1]:range(birdata$year)[2],
  species_binomial = all_species)

birdata_checklist_grd <- birdata %>% 
  group_by(cell_id, year) %>% 
  summarise(n_checklists.birdata = n_distinct(checklist_id))

birdata_full_grd <- left_join(birdata_year_grd, birdata_grd, by = c("cell_id", "year", "species_binomial")) %>% 
  mutate(n_obs.birdata = replace_na(n_obs.birdata, 0)) %>% 
  left_join(., birdata_checklist_grd, by = c("cell_id", "year")) %>% 
  mutate(n_checklists.birdata = replace_na(n_checklists.birdata, 0))

full_grd <- left_join(ebird_full_grd, birdata_full_grd) %>% 
  mutate(n_obs.combined = n_obs.ebird + n_obs.birdata,
         n_checklists.combined = n_checklists.ebird + n_checklists.birdata) %>% 
  left_join(., add_grd, by = c("cell_id", "species_binomial"))

full_grd_cumulative <- full_grd %>% 
  group_by(species_binomial, cell_id) %>% 
  arrange(year) %>% 
  mutate(across(str_subset(colnames(full_grd), "\\."), ~cumsum(.x))) %>% 
  filter(year > 1989)

full_grd <- full_grd %>% 
  filter(year > 1989)

write_tsv(full_grd, paste0("Data/", cell_size, "/full_grid.txt"))
write_tsv(full_grd_cumulative, paste0("Data/", cell_size, "/full_grid_cumulative.txt"))

rm(list=ls()[! ls() %in% c("cell_size", "aus", "full_grd", "grd", "all_species", 
                           "trait", "%notin%", "full_grd_cumulative")])

save.image()
gc()

#### 3. Plots ------------------------------------------------------------------

if (length(str_subset(list.files(paste0("Plots/Species/", cell_size, "/")), "png")) != 0) {
  rebuild_plots <- askYesNo(msg = "Do you want to rebuild the Map/Record Plots?")
} else {
  rebuild_plots <- TRUE
}

if (rebuild_plots == TRUE) {
  
  sf_use_s2(FALSE)
  
  plot_grd <- full_grd_cumulative %>% 
    filter(year == 2022) %>% 
    mutate(in_ebird = n_obs.ebird > 0,
           in_birdata = n_obs.birdata > 0) %>% 
    select(cell_id, species_binomial, in_ebird, in_birdata, in_abg, in_bli) %>% 
    left_join(grd, ., by = "cell_id")
  
  for (species in all_species) {
    abg_species <- plot_grd %>% 
      filter(species_binomial == species,
             in_abg == TRUE)
    bli_species <- plot_grd %>% 
      filter(species_binomial == species,
             in_bli == TRUE)
    ebird_species <- plot_grd %>% 
      filter(species_binomial == species,
             in_ebird == TRUE) %>% 
      st_buffer(dist = -set_units(0.2, degrees))
    birdata_species <- plot_grd %>% 
      filter(species_binomial == species,
             in_birdata == TRUE) %>% 
      st_buffer(dist = -set_units(0.2, degrees))
    
    plot <- ggplot() +
      geom_sf(data = aus) +
      geom_sf(data = abg_species, fill = "blue", alpha = 0.3) +
      geom_sf(data = bli_species, fill = "red", alpha = 0.3) +
      geom_sf(data = ebird_species, fill = "green", alpha = 0.3) +
      geom_sf(data = birdata_species, fill = "blue", alpha = 0.3) +
      theme_classic()
    
    ggsave(paste0("Plots/Species/", cell_size, "/", species, ".png"), plot, width = 7.5, height = 4.5)
  }
  
  sf_use_s2(TRUE)
  
}

#### 4. Metrics ----------------------------------------------------------------

hoover <- function (x, distribution = NULL) {
  if (is.null(distribution)) {
    return (0.5 * sum(abs(x - mean(x))) / sum(x))
  }
  
  if (length(x) != length(distribution)) {
    stop ("Vector lengths are not equal")
  }
  
  if (isTRUE(all.equal(sum(x), sum(distribution))) == FALSE) {
    new_x <- x * sum(distribution) / sum(x)
    warning ("Vector sums are not equal. Applying transformation to make them equal: this may have unexpected results!")
    return (c(0.5 * sum(abs(new_x - distribution)) / sum(new_x), sum(distribution) / sum(x)))
  }
  
  return (0.5 * sum(abs(x - distribution)) / sum(x))
}

##### 4.1. Mean Inventory Completeness -----------------------------------------

mic <- full_grd %>%
  filter(in_either == TRUE) %>%
  group_by(cell_id, year) %>%
  summarise(expected_diversity = n_distinct(species_binomial),
            observed_diversity.ebird = sum(n_obs.ebird > 0, na.rm = TRUE),
            observed_diversity.birdata = sum(n_obs.birdata > 0, na.rm = TRUE),
            observed_diversity.combined = sum(n_obs.combined > 0, na.rm = TRUE)) %>%
  mutate(diversity_completeness.ebird = observed_diversity.ebird / expected_diversity,
         diversity_completeness.birdata = observed_diversity.birdata / expected_diversity,
         diversity_completeness.combined = observed_diversity.combined / expected_diversity) %>%
  left_join(full_grd, .) %>%
  filter(in_either == TRUE) %>%
  group_by(species_binomial, year) %>%
  summarise(mean_inventory_completeness.ebird = weighted.mean(x = diversity_completeness.ebird, w = prop_land, na.rm = T),
            mean_inventory_completeness.birdata = weighted.mean(x = diversity_completeness.birdata, w = prop_land, na.rm = T),
            mean_inventory_completeness.combined = weighted.mean(x = diversity_completeness.combined, w = prop_land, na.rm = T))

mic_cumulative <- full_grd_cumulative %>%
  filter(in_either == TRUE) %>%
  group_by(cell_id, year) %>%
  summarise(expected_diversity = n_distinct(species_binomial),
            observed_diversity.ebird = sum(n_obs.ebird > 0, na.rm = TRUE),
            observed_diversity.birdata = sum(n_obs.birdata > 0, na.rm = TRUE),
            observed_diversity.combined = sum(n_obs.combined > 0, na.rm = TRUE)) %>%
  mutate(diversity_completeness.ebird = observed_diversity.ebird / expected_diversity,
         diversity_completeness.birdata = observed_diversity.birdata / expected_diversity,
         diversity_completeness.combined = observed_diversity.combined / expected_diversity) %>%
  left_join(full_grd_cumulative, .) %>%
  filter(in_either == TRUE) %>%
  group_by(species_binomial, year) %>%
  summarise(mean_inventory_completeness.ebird = weighted.mean(x = diversity_completeness.ebird, w = prop_land, na.rm = T),
            mean_inventory_completeness.birdata = weighted.mean(x = diversity_completeness.birdata, w = prop_land, na.rm = T),
            mean_inventory_completeness.combined = weighted.mean(x = diversity_completeness.combined, w = prop_land, na.rm = T))

##### 4.2. Total Range Completeness --------------------------------------------

trc <- full_grd %>% 
  filter(in_either == TRUE) %>% 
  mutate(detected.ebird = replace_na(n_obs.ebird > 0, FALSE),
         detected.birdata = replace_na(n_obs.birdata > 0, FALSE),
         detected.combined = replace_na(n_obs.combined > 0, FALSE)) %>% 
  group_by(species_binomial, year) %>% 
  summarise(total_range_completeness.ebird = weighted.mean(detected.ebird, prop_land, na.rm = T),
            total_range_completeness.birdata = weighted.mean(detected.birdata, prop_land, na.rm = T),
            total_range_completeness.combined = weighted.mean(detected.combined, prop_land, na.rm = T))

trc_cumulative <- full_grd_cumulative %>% 
  filter(in_either == TRUE) %>% 
  mutate(detected.ebird = replace_na(n_obs.ebird > 0, FALSE),
         detected.birdata = replace_na(n_obs.birdata > 0, FALSE),
         detected.combined = replace_na(n_obs.combined > 0, FALSE)) %>% 
  group_by(species_binomial, year) %>% 
  summarise(total_range_completeness.ebird = weighted.mean(detected.ebird, prop_land, na.rm = T),
            total_range_completeness.birdata = weighted.mean(detected.birdata, prop_land, na.rm = T),
            total_range_completeness.combined = weighted.mean(detected.combined, prop_land, na.rm = T))

##### 4.3. Total Number of Records ---------------------------------------------

tnr <- full_grd %>% 
  group_by(species_binomial, year) %>% 
  summarise(total_number_records.ebird = sum(n_obs.ebird),
            total_number_records.birdata = sum(n_obs.birdata),
            total_number_records.combined = sum(n_obs.combined))

tnr_cumulative <- full_grd_cumulative %>% 
  group_by(species_binomial, year) %>% 
  summarise(total_number_records.ebird = sum(n_obs.ebird),
            total_number_records.birdata = sum(n_obs.birdata),
            total_number_records.combined = sum(n_obs.combined))

##### 4.4. Checklist Spatial Bias ----------------------------------------------

csb <- full_grd %>% 
  filter(in_either == TRUE) %>% 
  group_by(species_binomial, year) %>% 
  summarise(checklist_spatial_bias.ebird = 1-suppressWarnings(hoover(x = n_checklists.ebird, distribution = prop_land)[1]),
            checklist_spatial_bias.birdata = 1-suppressWarnings(hoover(x = n_checklists.birdata, distribution = prop_land)[1]),
            checklist_spatial_bias.combined = 1-suppressWarnings(hoover(x = n_checklists.combined, distribution = prop_land)[1]))

csb_cumulative <- full_grd_cumulative %>% 
  filter(in_either == TRUE) %>% 
  group_by(species_binomial, year) %>% 
  summarise(checklist_spatial_bias.ebird = 1-suppressWarnings(hoover(x = n_checklists.ebird, distribution = prop_land)[1]),
            checklist_spatial_bias.birdata = 1-suppressWarnings(hoover(x = n_checklists.birdata, distribution = prop_land)[1]),
            checklist_spatial_bias.combined = 1-suppressWarnings(hoover(x = n_checklists.combined, distribution = prop_land)[1]))

#### 5. Plots ------------------------------------------------------------------

##### 5.1. Mean Inventory Completeness -----------------------------------------

mic_year_plot <- mic %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(mean_inventory_completeness.ebird),
            q5.ebird = quantile(mean_inventory_completeness.ebird, 0.05),
            q95.ebird = quantile(mean_inventory_completeness.ebird, 0.95),
            median.birdata = median(mean_inventory_completeness.birdata),
            q5.birdata = quantile(mean_inventory_completeness.birdata, 0.05),
            q95.birdata = quantile(mean_inventory_completeness.birdata, 0.95),
            median.combined = median(mean_inventory_completeness.combined),
            q5.combined = quantile(mean_inventory_completeness.combined, 0.05),
            q95.combined = quantile(mean_inventory_completeness.combined, 0.95)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Mean Inventory Completeness\n(How adequately surveyed is this species' range?)") + 
  lims(y = c(0, 1))

mic_cumulative_plot <- mic_cumulative %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(mean_inventory_completeness.ebird),
            q5.ebird = quantile(mean_inventory_completeness.ebird, 0.05),
            q95.ebird = quantile(mean_inventory_completeness.ebird, 0.95),
            median.birdata = median(mean_inventory_completeness.birdata),
            q5.birdata = quantile(mean_inventory_completeness.birdata, 0.05),
            q95.birdata = quantile(mean_inventory_completeness.birdata, 0.95),
            median.combined = median(mean_inventory_completeness.combined),
            q5.combined = quantile(mean_inventory_completeness.combined, 0.05),
            q95.combined = quantile(mean_inventory_completeness.combined, 0.95)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Mean Inventory Completeness\n(How adequately surveyed is this species' range?)") + 
  lims(y = c(0, 1))

mic_final_plot <- mic_cumulative %>% 
  filter(year == 2022) %>% 
  select(-year) %>% 
  pivot_longer(!species_binomial, names_to = "dataset", values_to = "mean_inventory_completeness") %>% 
  mutate(dataset = as_factor(str_replace(dataset, "mean_inventory_completeness\\.", ""))) %>% 
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = dataset, y = mean_inventory_completeness, fill = dataset), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("ebird" = "green", "birdata" = "blue", "combined" = "black")) +
  scale_x_discrete(labels = c("ebird" = "eBird", "birdata" = "Birdata", "combined" = "Combined")) +
  labs(x = "Dataset", y = "Mean Inventory Completeness\n(How adequately surveyed is this species' range?)") +
  theme(legend.position = "none") +
  lims(y = c(0, 1))

##### 5.2. Total Range Completeness --------------------------------------------

trc_year_plot <- trc %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(total_range_completeness.ebird),
            q5.ebird = quantile(total_range_completeness.ebird, 0.05),
            q95.ebird = quantile(total_range_completeness.ebird, 0.95),
            median.birdata = median(total_range_completeness.birdata),
            q5.birdata = quantile(total_range_completeness.birdata, 0.05),
            q95.birdata = quantile(total_range_completeness.birdata, 0.95),
            median.combined = median(total_range_completeness.combined),
            q5.combined = quantile(total_range_completeness.combined, 0.05),
            q95.combined = quantile(total_range_completeness.combined, 0.95)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Total Range Completeness\n(How well surveyed is this species across its range?)") +
  lims(y = c(0, 1))

trc_cumulative_plot <- trc_cumulative %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(total_range_completeness.ebird),
            q5.ebird = quantile(total_range_completeness.ebird, 0.05),
            q95.ebird = quantile(total_range_completeness.ebird, 0.95),
            median.birdata = median(total_range_completeness.birdata),
            q5.birdata = quantile(total_range_completeness.birdata, 0.05),
            q95.birdata = quantile(total_range_completeness.birdata, 0.95),
            median.combined = median(total_range_completeness.combined),
            q5.combined = quantile(total_range_completeness.combined, 0.05),
            q95.combined = quantile(total_range_completeness.combined, 0.95)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Total Range Completeness\n(How well surveyed is this species across its range?)") +
  lims(y = c(0, 1))

trc_final_plot <- trc_cumulative %>% 
  filter(year == 2022) %>% 
  select(-year) %>% 
  pivot_longer(!species_binomial, names_to = "dataset", values_to = "total_range_completeness") %>% 
  mutate(dataset = as_factor(str_replace(dataset, "total_range_completeness\\.", ""))) %>% 
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = dataset, y = total_range_completeness, fill = dataset), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("ebird" = "green", "birdata" = "blue", "combined" = "black")) +
  scale_x_discrete(labels = c("ebird" = "eBird", "birdata" = "Birdata", "combined" = "Combined")) +
  labs(x = "Dataset", y = "Total Range Completeness\n(How well surveyed is this species across its range?)") +
  theme(legend.position = "none") +
  lims(y = c(0, 1))

##### 5.3. Total Number of Records ----------------------------------------------

tnr_year_plot <- tnr %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(total_number_records.ebird),
            q5.ebird = quantile(total_number_records.ebird, 0.05),
            q95.ebird = quantile(total_number_records.ebird, 0.95),
            median.birdata = median(total_number_records.birdata),
            q5.birdata = quantile(total_number_records.birdata, 0.05),
            q95.birdata = quantile(total_number_records.birdata, 0.95),
            median.combined = median(total_number_records.combined),
            q5.combined = quantile(total_number_records.combined, 0.05),
            q95.combined = quantile(total_number_records.combined, 0.95)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Total Number of Records\n(How much data do we have on this species?)") +
  scale_y_log10(labels = scales::comma) 

tnr_cumulative_plot <- tnr_cumulative %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(total_number_records.ebird),
            q5.ebird = quantile(total_number_records.ebird, 0.05),
            q95.ebird = quantile(total_number_records.ebird, 0.95),
            median.birdata = median(total_number_records.birdata),
            q5.birdata = quantile(total_number_records.birdata, 0.05),
            q95.birdata = quantile(total_number_records.birdata, 0.95),
            median.combined = median(total_number_records.combined),
            q5.combined = quantile(total_number_records.combined, 0.05),
            q95.combined = quantile(total_number_records.combined, 0.95)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Total Number of Records\n(How much data do we have on this species?)") +
  scale_y_log10(labels = scales::comma) 

tnr_final_plot <- tnr_cumulative %>% 
  filter(year == 2022) %>% 
  select(-year) %>% 
  pivot_longer(!species_binomial, names_to = "dataset", values_to = "total_number_records") %>% 
  mutate(dataset = as_factor(str_replace(dataset, "total_number_records\\.", "")),
         total_number_records = if_else(total_number_records == 0, 1, total_number_records)) %>% # to deal with log-transformation issues!
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = dataset, y = total_number_records, fill = dataset), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("ebird" = "green", "birdata" = "blue", "combined" = "black")) +
  scale_x_discrete(labels = c("ebird" = "eBird", "birdata" = "Birdata", "combined" = "Combined")) +
  labs(x = "Dataset", y = "Total Number of Records\n(How much data do we have on this species?)") +
  theme(legend.position = "none") +
  scale_y_log10(labels = scales::comma)

##### 5.4. Checklist Spatial Bias ----------------------------------------------

csb_year_plot <- csb %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(checklist_spatial_bias.ebird, na.rm = T),
            q5.ebird = quantile(checklist_spatial_bias.ebird, 0.05, na.rm = T),
            q95.ebird = quantile(checklist_spatial_bias.ebird, 0.95, na.rm = T),
            median.birdata = median(checklist_spatial_bias.birdata, na.rm = T),
            q5.birdata = quantile(checklist_spatial_bias.birdata, 0.05, na.rm = T),
            q95.birdata = quantile(checklist_spatial_bias.birdata, 0.95, na.rm = T),
            median.combined = median(checklist_spatial_bias.combined, na.rm = T),
            q5.combined = quantile(checklist_spatial_bias.combined, 0.05, na.rm = T),
            q95.combined = quantile(checklist_spatial_bias.combined, 0.95, na.rm = T)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Checklist Spatial Bias\n(How biased are the data we have for this species?)") +
  lims(y = c(0, 1))

csb_cumulative_plot <- csb_cumulative %>% 
  group_by(year) %>% 
  summarise(median.ebird = median(checklist_spatial_bias.ebird, na.rm = T),
            q5.ebird = quantile(checklist_spatial_bias.ebird, 0.05, na.rm = T),
            q95.ebird = quantile(checklist_spatial_bias.ebird, 0.95, na.rm = T),
            median.birdata = median(checklist_spatial_bias.birdata, na.rm = T),
            q5.birdata = quantile(checklist_spatial_bias.birdata, 0.05, na.rm = T),
            q95.birdata = quantile(checklist_spatial_bias.birdata, 0.95, na.rm = T),
            median.combined = median(checklist_spatial_bias.combined, na.rm = T),
            q5.combined = quantile(checklist_spatial_bias.combined, 0.05, na.rm = T),
            q95.combined = quantile(checklist_spatial_bias.combined, 0.95, na.rm = T)) %>% 
  ggplot(aes(x = year)) +
  theme_classic() +
  geom_ribbon(aes(ymin = q5.ebird, ymax = q95.ebird), alpha = 0.1, fill = "green") +
  geom_line(aes(y = median.ebird), linewidth = 1, colour = "green") +
  geom_ribbon(aes(ymin = q5.birdata, ymax = q95.birdata), alpha = 0.1, fill = "blue") +
  geom_line(aes(y = median.birdata), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(ymin = q5.combined, ymax = q95.combined), alpha = 0.1, fill = "black") +
  geom_line(aes(y = median.combined), linewidth = 1, colour = "black") +
  labs(x = "Year", y = "Checklist Spatial Bias\n(How biased are the data we have for this species?)") +
  lims(y = c(0, 1))

csb_final_plot <- csb_cumulative %>% 
  filter(year == 2022) %>% 
  select(-year) %>% 
  pivot_longer(!species_binomial, names_to = "dataset", values_to = "checklist_spatial_bias") %>% 
  mutate(dataset = as_factor(str_replace(dataset, "checklist_spatial_bias\\.", ""))) %>% 
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = dataset, y = checklist_spatial_bias, fill = dataset), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("ebird" = "green", "birdata" = "blue", "combined" = "black")) +
  scale_x_discrete(labels = c("ebird" = "eBird", "birdata" = "Birdata", "combined" = "Combined")) +
  labs(x = "Dataset", y = "Checklist Spatial Bias\n(How biased are the data we have for this species?)") +
  theme(legend.position = "none") +
  lims(y = c(0, 1))

##### 5.7. Combined ------------------------------------------------------------

combined_year_plot <- ggarrange(
  mic_year_plot, trc_year_plot, tnr_year_plot, csb_year_plot,
  ncol = 2, nrow = 2, align = "hv"
)

ggsave(plot = combined_year_plot, paste0("Plots/", cell_size, "/_year.png"), height = 10, width = 10)

combined_cumulative_plot <- ggarrange(
  mic_cumulative_plot,   trc_cumulative_plot, tnr_cumulative_plot, csb_cumulative_plot,
  ncol = 2, nrow = 2, align = "hv"
)

ggsave(plot = combined_cumulative_plot, paste0("Plots/", cell_size, "/_cumulative.png"), height = 10, width = 10)

combined_final_plot <- ggarrange(
  mic_final_plot, trc_final_plot, tnr_final_plot, csb_final_plot,
  ncol = 2, nrow = 2, align = "hv"
)

ggsave(plot = combined_final_plot, paste0("Plots/", cell_size, "/_final.png"), height = 10, width = 10)

#### 6. Clean Up ---------------------------------------------------------------

combined <- reduce(list(mic, trc, tnr, csb), left_join)
combined_cumulative <- reduce(
  list(mic_cumulative, trc_cumulative, tnr_cumulative, csb_cumulative), left_join)

write_tsv(combined, paste0("Data/", cell_size, "/combined.txt"))
write_tsv(combined_cumulative, paste0("Data/", cell_size, "/combined_cumulative.txt"))

save.image()
gc()

#### 7. Threat Analysis --------------------------------------------------------

trait_combined <- left_join(combined, trait, by = "species_binomial") %>% 
  mutate(threat_status = factor(threat_status, levels = c("LC/NA", "NT", "VU", "EN", "CR")))
trait_cumulative <- left_join(combined_cumulative, trait, by = "species_binomial") %>% 
  mutate(threat_status = factor(threat_status, levels = c("LC/NA", "NT", "VU", "EN", "CR")))

mic_threat_plot <- trait_cumulative %>% 
  filter(year == 2022) %>% 
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = threat_status, y = mean_inventory_completeness.combined, fill = threat_status), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("LC/NA" = "#60c659", "NT" = "#cce226", "VU" = "#f9e814", "EN" = "#fc7f3f", "CR" = "#d81e05")) +
  labs(x = "Threat Status (APAB 2020)", y = "Mean Inventory Completeness\n(How adequately surveyed is this species' range?)") +
  theme(legend.position = "none") +
  lims(y = c(0, 1))

trc_threat_plot <- trait_cumulative %>% 
  filter(year == 2022) %>% 
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = threat_status, y = total_range_completeness.combined, fill = threat_status), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("LC/NA" = "#60c659", "NT" = "#cce226", "VU" = "#f9e814", "EN" = "#fc7f3f", "CR" = "#d81e05")) +
  labs(x = "Threat Status (APAB 2020)", y = "Total Range Completeness\n(How well surveyed is this species across its range?)") +
  theme(legend.position = "none") +
  lims(y = c(0, 1))

tnr_threat_plot <- trait_cumulative %>% 
  filter(year == 2022) %>% 
  mutate(total_number_records.combined = if_else(total_number_records.combined == 0, 1, total_number_records.combined)) %>% # to deal with log-transformation issues!
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = threat_status, y = total_number_records.combined, fill = threat_status), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("LC/NA" = "#60c659", "NT" = "#cce226", "VU" = "#f9e814", "EN" = "#fc7f3f", "CR" = "#d81e05")) +
  labs(x = "Threat Status (APAB 2020)", y = "Total Number of Records\n(How much data do we have on this species?)") +
  theme(legend.position = "none") +
  scale_y_log10(labels = scales::comma)

csb_threat_plot <- trait_cumulative %>% 
  filter(year == 2022) %>% 
  ggplot() +
  theme_classic() +
  geom_violin(aes(x = threat_status, y = checklist_spatial_bias.combined, fill = threat_status), alpha = 0.5, linewidth = 1, draw_quantiles = 0.5) +
  scale_fill_manual(values = c("LC/NA" = "#60c659", "NT" = "#cce226", "VU" = "#f9e814", "EN" = "#fc7f3f", "CR" = "#d81e05")) +
  labs(x = "Threat Status (APAB 2020)", y = "Checklist Spatial Bias\n(How biased are the data we have for this species?)") +
  theme(legend.position = "none") +
  lims(y = c(0, 1))

combined_threat_plot <- ggarrange(
  mic_threat_plot, trc_threat_plot, tnr_threat_plot, csb_threat_plot,
  ncol = 2, nrow = 2, align = "hv"
)

ggsave(plot = combined_threat_plot, paste0("Plots/", cell_size, "/_threat.png"), height = 10, width = 10)

#### 8. Trait Modeling ---------------------------------------------------------

model_data <- combined_cumulative %>% 
  filter(year == 2022) %>% 
  left_join(., trait) %>% 
  rowid_to_column("species_id")

mic_model <- lm(
  data = model_data,
  formula = mean_inventory_completeness.combined ~ threat_status + 
    log10(species_uniqueness) + log10(body_mass) + log10(average_count) + 
    log10(range_size) + log10(density) + log10(human_density)) 

mic_summary <- summary(mic_model)
mic_check <- check_model(mic_model)
mic_residuals <- mic_model$residuals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(species_id = 1, mic_residual = 2) %>% 
  mutate(species_id = as.integer(species_id))

trc_model <- lm(
  data = model_data,
  formula = total_range_completeness.combined ~ threat_status + 
    log10(species_uniqueness) + log10(body_mass) + log10(average_count) + 
    log10(range_size) + log10(density) + log10(human_density)) 

trc_summary <- summary(trc_model)
trc_check <- check_model(trc_model)
trc_residuals <- trc_model$residuals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(species_id = 1, trc_residual = 2) %>% 
  mutate(species_id = as.integer(species_id))

tnr_model <- lm(
  data = model_data,
  formula = log10(total_number_records.combined) ~ threat_status + 
    log10(species_uniqueness) + log10(body_mass) + log10(average_count) + 
    log10(range_size) + log10(density) + log10(human_density)) 

tnr_summary <- summary(tnr_model)
tnr_check <- check_model(tnr_model)
tnr_residuals <- tnr_model$residuals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(species_id = 1, tnr_residual = 2) %>% 
  mutate(species_id = as.integer(species_id))

csb_model <- lm(
  data = model_data,
  formula = checklist_spatial_bias.combined ~ threat_status + 
    log10(species_uniqueness) + log10(body_mass) + log10(average_count) + 
    log10(range_size) + log10(density) + log10(human_density)) 

csb_summary <- summary(csb_model)
csb_check <- check_model(csb_model)
csb_residuals <- csb_model$residuals %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  rename(species_id = 1, csb_residual = 2) %>% 
  mutate(species_id = as.integer(species_id))

residual_combined <- reduce(
  list(model_data, mic_residuals, trc_residuals, tnr_residuals, csb_residuals), left_join)

coefficients_combined <- bind_rows(
  rownames_to_column(mutate(as.data.frame(mic_summary$coefficients), Metric = "MIC"), var = "Coefficient"),
  rownames_to_column(mutate(as.data.frame(trc_summary$coefficients), Metric = "TRC"), var = "Coefficient"),
  rownames_to_column(mutate(as.data.frame(tnr_summary$coefficients), Metric = "TNR"), var = "Coefficient"),
  rownames_to_column(mutate(as.data.frame(csb_summary$coefficients), Metric = "CSB"), var = "Coefficient")
) %>% 
  mutate(Metric = factor(Metric, levels = c("CSB", "TNR", "TRC", "MIC")))

coefficients_plot <- ggplot(data = coefficients_combined, aes(x = Estimate, y = Coefficient, group = Metric, colour = Metric)) +
  geom_point(position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(xmin = Estimate - qnorm(0.975) * `Std. Error`, xmax = Estimate + qnorm(0.975) * `Std. Error`),
                position = position_dodge(width = 0.75)) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_discrete(
    breaks = sort(unique(coefficients_combined$Coefficient), decreasing = TRUE),
    labels = c("Threat Status: VU", "Threat Status: NT", "Threat Status: LC", 
               "Threat Status: EN", "Species Uniqueness", "Range Size", 
               "Human Density", "Density", "Body Mass", "Average Count", "Intercept")
  )

ggsave(plot = coefficients_plot, paste0("Plots/", cell_size, "/_coefficients.png"), height = 5, width = 5)  
ggsave(plot = coefficients_plot, paste0("Plots/", cell_size, "/_coefficients.svg"), height = 5, width = 5)  
