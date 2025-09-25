archivo="C:/Users/jabra/Dropbox/Posdoc/Hybrids/SampleSheet/Sacc718_8.csv"
archivo="C:/Users/jabra/Dropbox/Posdoc/Hybrids/SampleSheet/Sacc718_8_build.csv" #Correcciones detectadas durante el análisis

df=read.csv(archivo)
tousedf_set1 = df[df$STRAIN_SETS == 1,]
head(tousedf_set1)


#### Parragraph ####

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(glue)
  library(tidyr)
})

# ====== Limpieza y preparación ======
tousedf_set1_clean <- tousedf_set1 %>%
  mutate(
    GenomSpecies = toupper(str_trim(GenomSpecies)),
    Region       = str_trim(Region),
    Location     = str_trim(Location),
    OBJ_ID       = str_trim(OBJ_ID)
  )

species_expected <- c("SACE", "SAPA", "HYBRID")

tousedf_species_filtered <- tousedf_set1_clean %>%
  filter(GenomSpecies %in% species_expected)

# Helper para porcentajes seguros (1 decimal)
safe_pct <- function(num, den) {
  ifelse(den > 0, round(100 * num / den, 1), 0)
}

# ====== Totales globales ======
total_isolates <- nrow(tousedf_species_filtered)

species_counts_summary <- tousedf_species_filtered %>%
  count(GenomSpecies) %>%
  complete(GenomSpecies = species_expected, fill = list(n = 0)) %>%
  mutate(percent_global = safe_pct(n, total_isolates))

n_sace   <- species_counts_summary %>% filter(GenomSpecies == "SACE")   %>% pull(n)
p_sace   <- species_counts_summary %>% filter(GenomSpecies == "SACE")   %>% pull(percent_global)

n_sapa   <- species_counts_summary %>% filter(GenomSpecies == "SAPA")   %>% pull(n)
p_sapa   <- species_counts_summary %>% filter(GenomSpecies == "SAPA")   %>% pull(percent_global)

n_hybrid <- species_counts_summary %>% filter(GenomSpecies == "HYBRID") %>% pull(n)
p_hybrid <- species_counts_summary %>% filter(GenomSpecies == "HYBRID") %>% pull(percent_global)

# ====== Totales por región (para West_I y South_central) ======
region_totals <- tousedf_species_filtered %>%
  count(Region, name = "n_region_total")

# West_I → híbridos
west_region_total <- region_totals %>%
  filter(Region == "West_I") %>%
  pull(n_region_total) %>% { if (length(.) == 0) 0 else . }

west_hybrids_count <- tousedf_species_filtered %>%
  filter(Region == "West_I", GenomSpecies == "HYBRID") %>%
  nrow()

west_hybrids_pct_global <- safe_pct(west_hybrids_count, total_isolates)      # % sobre total global
west_hybrids_pct_region <- safe_pct(west_hybrids_count, west_region_total)   # % dentro de West_I

# South_central → S. cerevisiae (SACE)
south_region_total <- region_totals %>%
  filter(Region == "South_central") %>%
  pull(n_region_total) %>% { if (length(.) == 0) 0 else . }

south_sace_count <- tousedf_species_filtered %>%
  filter(Region == "South_central", GenomSpecies == "SACE") %>%
  nrow()

south_sace_pct_global <- safe_pct(south_sace_count, total_isolates)        # % sobre total global
south_sace_pct_region <- safe_pct(south_sace_count, south_region_total)    # % dentro de South_central

# ====== Locations y crioviales con >1 especie ======
multi_species_locations <- tousedf_species_filtered %>%
  group_by(Location) %>%
  summarize(n_species = n_distinct(GenomSpecies), .groups = "drop") %>%
  filter(!is.na(Location), n_species >= 2) %>%
  nrow()
p_multi_species_locations <- safe_pct(multi_species_locations, n_distinct(tousedf_species_filtered$Location)-1) #-1 para quitar el 'empty'

multi_species_criovials <- tousedf_species_filtered %>%
  group_by(OBJ_ID) %>%
  summarize(n_species = n_distinct(GenomSpecies), .groups = "drop") %>%
  filter(!is.na(OBJ_ID), n_species >= 2) %>%
  nrow()
p_multi_species_criovials <- safe_pct(multi_species_criovials, n_distinct(tousedf_species_filtered$OBJ_ID)-1)#-1 para quitar el 'empty'

# ====== Párrafo final ======
paragraph <- glue(
  "In total, there are {total_isolates} sequenced isolates from fermentation tanks. \n",
  "Out of them, {n_sace} ({p_sace}%) are S. cerevisiae, ",
  "{n_sapa} ({p_sapa}%) are S. paradoxus and ",
  "{n_hybrid} ({p_hybrid}%) are hybrids. \n",
  "The three species are not evenly distributed across the country. \n",
  "In West_I, {west_hybrids_count} out of {west_region_total} sequences ",
  "are hybrids (",
  "{west_hybrids_pct_region}% within West_I), while in South_central, ",
  "{south_sace_count} out of {south_region_total} ({south_sace_pct_region}%) sequences ",
  "are S. cerevisiae  (FIGURE 1). \n",
  "Differential distribution of the three species remains after controlling for the sampling effort ",
  "in the number of tanks, and number of locations (Supp table 1). \n",
  "We found {multi_species_locations} locations ({p_multi_species_locations}%) and ",
  "{multi_species_criovials} criovial-samples ({p_multi_species_criovials}%) in which more than one species was found, ",
  "suggesting that the fermentation system is important for interactions between the species\n"
)

cat(paragraph, "\n")


#### Tabla por región region_summary ####
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
})

# --- Punto de partida: tousedf_species_filtered del script anterior ---

# Conteo por región y especie
region_species_counts <- tousedf_species_filtered %>%
  count(Region, GenomSpecies, name = "n_species") %>%
  pivot_wider(
    names_from = GenomSpecies,
    values_from = n_species,
    values_fill = 0
  )

# Total de muestras por región
region_totals <- tousedf_species_filtered %>%
  count(Region, name = "Number_of_samples")

# Número de locations únicos por región
region_locations <- tousedf_species_filtered %>%
  group_by(Region) %>%
  summarize(Number_of_locations = n_distinct(Location), .groups = "drop")

# Unir todo
region_summary <- region_totals %>%
  left_join(region_locations, by = "Region") %>%
  left_join(region_species_counts, by = "Region") %>%
  mutate(
    perc_SACE   = round(100 * SACE   / Number_of_samples, 1),
    perc_SAPA   = round(100 * SAPA   / Number_of_samples, 1),
    perc_HYBRID = round(100 * HYBRID / Number_of_samples, 1)
  ) %>%
  arrange(Region)

# Mostrar tabla
region_summary
#####Tabla por OBJ_ID o por Location repetido multi_species_locations_table multi_sequence_objids  #####
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# Asegura que por cada OBJ_ID existan las tres especies como columnas (rellena con 0)
objid_species_counts <- tousedf_species_filtered %>%
  count(OBJ_ID, GenomSpecies, name = "n_especie") %>%
  group_by(OBJ_ID) %>%
  tidyr::complete(
    GenomSpecies = c("SACE", "SAPA", "HYBRID"),
    fill = list(n_especie = 0)
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = GenomSpecies,
    values_from = n_especie
  )

# Total de secuencias por OBJ_ID
objid_total_sequences <- tousedf_species_filtered %>%
  count(OBJ_ID, name = "seqs")

# Location y Region por OBJ_ID
objid_loc_region <- tousedf_species_filtered %>%
  group_by(OBJ_ID) %>%
  summarize(
    location = paste(sort(unique(Location)), collapse = " | "),
    region   = paste(sort(unique(Region)),   collapse = " | "),
    .groups  = "drop"
  )

# Unir todo, renombrar Hybrid -> Hyb, calcular multi_species correctamente y ordenar
multi_sequence_objids <- objid_total_sequences %>%
  left_join(objid_species_counts, by = "OBJ_ID") %>%
  left_join(objid_loc_region,     by = "OBJ_ID") %>%
  rename(Hyb = HYBRID) %>%
  mutate(
    multi_species = rowSums(across(c(SACE, SAPA, Hyb), ~ .x > 0)) > 1
  ) %>%
  select(OBJ_ID, seqs, SACE, SAPA, Hyb, multi_species, location, region) %>%
  arrange(desc(multi_species), desc(seqs))

# Mostrar la tabla
multi_sequence_objids

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# --- Conteos por Location y especie ---
location_species_counts <- tousedf_species_filtered %>%
  count(Location, GenomSpecies, name = "n_especie") %>%
  group_by(Location) %>%
  tidyr::complete(
    GenomSpecies = c("SACE", "SAPA", "HYBRID"),
    fill = list(n_especie = 0)
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = GenomSpecies,
    values_from = n_especie
  )

# --- Total de secuencias por Location ---
location_total_sequences <- tousedf_species_filtered %>%
  count(Location, name = "seqs")

# --- Region asociada a cada Location (si hubiera varias, concatenar) ---
location_region <- tousedf_species_filtered %>%
  group_by(Location) %>%
  summarize(
    region = paste(sort(unique(Region)), collapse = " | "),
    .groups = "drop"
  )

# --- Unir todo ---
multi_species_locations_table <- location_total_sequences %>%
  left_join(location_species_counts, by = "Location") %>%
  left_join(location_region,         by = "Location") %>%
  rename(Hyb = HYBRID) %>%
  mutate(
    multi_species = rowSums(across(c(SACE, SAPA, Hyb), ~ .x > 0)) > 1
  ) %>%
  select(Location, seqs, SACE, SAPA, Hyb, multi_species, region) %>%
  arrange(desc(multi_species), desc(seqs))

# Mostrar tabla
multi_species_locations_table

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# --- Conteos por Location y especie ---
location_species_counts <- tousedf_species_filtered %>%
  count(Location, GenomSpecies, name = "n_especie") %>%
  group_by(Location) %>%
  tidyr::complete(
    GenomSpecies = c("SACE", "SAPA", "HYBRID"),
    fill = list(n_especie = 0)
  ) %>%
  ungroup() %>%
  pivot_wider(
    names_from  = GenomSpecies,
    values_from = n_especie
  )

# --- Total de secuencias por Location ---
location_total_sequences <- tousedf_species_filtered %>%
  count(Location, name = "seqs")

# --- Region asociada a cada Location (si hubiera varias, concatenar) ---
location_region <- tousedf_species_filtered %>%
  group_by(Location) %>%
  summarize(
    region = paste(sort(unique(Region)), collapse = " | "),
    .groups = "drop"
  )

# --- Unir todo ---
multi_species_locations_table <- location_total_sequences %>%
  left_join(location_species_counts, by = "Location") %>%
  left_join(location_region,         by = "Location") %>%
  rename(Hyb = HYBRID) %>%
  mutate(
    multi_species = rowSums(across(c(SACE, SAPA, Hyb), ~ .x > 0)) > 1
  ) %>%
  select(Location, seqs, SACE, SAPA, Hyb, multi_species, region) %>%
  arrange(desc(multi_species), desc(seqs))

# Mostrar tabla








##### Bootstrap sampling! #####

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})

set.seed(123)  # reproducibilidad
n_reps <- 50   # número de repeticiones bootstrap

# ===== Helper para calcular "region_summary" a partir de un subconjunto =====
make_region_summary <- function(df_subset) {
  region_species_counts <- df_subset %>%
    count(Region, GenomSpecies, name = "n_species") %>%
    pivot_wider(
      names_from = GenomSpecies,
      values_from = n_species,
      values_fill = 0
    )
  
  region_totals <- df_subset %>%
    count(Region, name = "Number_of_samples")
  
  region_locations <- df_subset %>%
    group_by(Region) %>%
    summarize(Number_of_locations = n_distinct(Location), .groups = "drop")
  
  region_summary <- region_totals %>%
    left_join(region_locations, by = "Region") %>%
    left_join(region_species_counts, by = "Region") %>%
    mutate(
      perc_SACE   = 100 * SACE   / Number_of_samples,
      perc_SAPA   = 100 * SAPA   / Number_of_samples,
      perc_HYBRID = 100 * HYBRID / Number_of_samples
    )
  
  return(region_summary)
}

# ===== Bootstrap a nivel de OBJ_ID =====
boot_objid <- map_dfr(1:n_reps, function(i) {
  # tomar 1 secuencia al azar por OBJ_ID
  sampled_df <- tousedf_species_filtered %>%
    group_by(OBJ_ID) %>%
    sample_n(1) %>%
    ungroup()
  
  make_region_summary(sampled_df) %>%
    mutate(rep = i)
})

summary_objid_boot <- boot_objid %>%
  group_by(Region) %>%
  summarize(
    across(c(Number_of_samples, Number_of_locations, SACE, SAPA, HYBRID,
             perc_SACE, perc_SAPA, perc_HYBRID),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

# ===== Bootstrap a nivel de Location =====
boot_location <- map_dfr(1:n_reps, function(i) {
  # tomar 1 secuencia al azar por Location
  sampled_df <- tousedf_species_filtered %>%
    group_by(Location) %>%
    sample_n(1) %>%
    ungroup()
  
  make_region_summary(sampled_df) %>%
    mutate(rep = i)
})

summary_location_boot <- boot_location %>%
  group_by(Region) %>%
  summarize(
    across(c(Number_of_samples, Number_of_locations, SACE, SAPA, HYBRID,
             perc_SACE, perc_SAPA, perc_HYBRID),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

View(summary_objid_boot)
View(summary_location_boot)


# Gráficas ##########
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})

# ---------- Helpers ----------
species_levels <- c("SACE", "SAPA", "Hyb")

make_pie <- function(data_long, facet_col = NULL, title_text = "") {
  # Espera columnas: Region (opcional), Species, Count, Percent
  plot_data <- data_long %>%
    mutate(
      Species = factor(Species, levels = species_levels),
      Count   = ifelse(is.na(Count), 0, Count),
      Percent = ifelse(is.na(Percent), 0, Percent),
      label   = paste0(Species, "\n n=", Count, " | ", sprintf("%.1f", Percent), "%")
    )
  
  base <- ggplot(plot_data, aes(x = "", y = Count, fill = Species)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    geom_text(aes(label = ifelse(Count > 0, label, "")), ## Comentar o quitar comentario de esta linea para poner o quitar los textos de los pays
              position = position_stack(vjust = 0.5), size = 3) +## Comentar o quitar comentario de esta linea para poner o quitar los textos de los pays
    labs(x = NULL, y = NULL, title = title_text, fill = "Species") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  if (!is.null(facet_col)) {
    base <- base + facet_wrap(reformulate(facet_col), scales = "free", ncol = 3)
  }
  base
}

# ---------- Página 1: Todos los datos ----------
overall_long <- tousedf_species_filtered %>%
  mutate(GenomSpecies = ifelse(GenomSpecies == "HYBRID", "Hyb", GenomSpecies)) %>%
  count(GenomSpecies, name = "Count") %>%
  complete(GenomSpecies = species_levels, fill = list(Count = 0)) %>%
  rename(Species = GenomSpecies) %>%
  mutate(Percent = 100 * Count / sum(Count))

plot_page1 <- make_pie(overall_long, title_text = "Set1: species composition")

# ---------- Página 2: region_summary (un subplot por región) ----------
# Se asume que region_summary tiene columnas: Region, Number_of_samples, Number_of_locations, SACE, SAPA, HYBRID, perc_SACE, perc_SAPA, perc_HYBRID
region_summary_long <- region_summary %>%
  mutate(Hyb = HYBRID) %>%
  select(Region, SACE, SAPA, Hyb, perc_SACE, perc_SAPA, perc_HYBRID) %>%
  pivot_longer(cols = c(SACE, SAPA, Hyb), names_to = "Species", values_to = "Count") %>%
  left_join(
    region_summary %>%
      transmute(
        Region,
        p_SACE = perc_SACE,
        p_SAPA = perc_SAPA,
        p_Hyb  = perc_HYBRID
      ) %>%
      pivot_longer(cols = c(p_SACE, p_SAPA, p_Hyb),
                   names_pattern = "p_(.*)", names_to = "Species", values_to = "Percent"),
    by = c("Region", "Species")
  ) %>%
  mutate(Species = factor(Species, levels = species_levels))

plot_page2 <- make_pie(region_summary_long, facet_col = "Region",
                       title_text = "Per region (region_summary): species composition")

# ---------- Página 3: summary_objid_boot (medias; un subplot por región) ----------
# Se asume que summary_objid_boot tiene columnas con *_mean y perc_*_mean
objid_boot_long <- summary_objid_boot %>%
  transmute(
    Region,
    SACE = SACE_mean,
    SAPA = SAPA_mean,
    Hyb  = HYBRID_mean,
    perc_SACE = perc_SACE_mean,
    perc_SAPA = perc_SAPA_mean,
    perc_Hyb  = perc_HYBRID_mean
  ) %>%
  pivot_longer(cols = c(SACE, SAPA, Hyb), names_to = "Species", values_to = "Count") %>%
  left_join(
    summary_objid_boot %>%
      transmute(
        Region,
        p_SACE = perc_SACE_mean,
        p_SAPA = perc_SAPA_mean,
        p_Hyb  = perc_HYBRID_mean
      ) %>%
      pivot_longer(cols = c(p_SACE, p_SAPA, p_Hyb),
                   names_pattern = "p_(.*)", names_to = "Species", values_to = "Percent"),
    by = c("Region", "Species")
  ) %>%
  mutate(
    Species = factor(Species, levels = species_levels),
    Count   = replace_na(Count, 0),
    Percent = replace_na(Percent, 0)
  )

plot_page3 <- make_pie(objid_boot_long, facet_col = "Region",
                       title_text = "Per region (bootstrap by OBJ_ID): mean species composition")

# ---------- Página 4: summary_location_boot (medias; un subplot por región) ----------
location_boot_long <- summary_location_boot %>%
  transmute(
    Region,
    SACE = SACE_mean,
    SAPA = SAPA_mean,
    Hyb  = HYBRID_mean,
    perc_SACE = perc_SACE_mean,
    perc_SAPA = perc_SAPA_mean,
    perc_Hyb  = perc_HYBRID_mean
  ) %>%
  pivot_longer(cols = c(SACE, SAPA, Hyb), names_to = "Species", values_to = "Count") %>%
  left_join(
    summary_location_boot %>%
      transmute(
        Region,
        p_SACE = perc_SACE_mean,
        p_SAPA = perc_SAPA_mean,
        p_Hyb  = perc_HYBRID_mean
      ) %>%
      pivot_longer(cols = c(p_SACE, p_SAPA, p_Hyb),
                   names_pattern = "p_(.*)", names_to = "Species", values_to = "Percent"),
    by = c("Region", "Species")
  ) %>%
  mutate(
    Species = factor(Species, levels = species_levels),
    Count   = replace_na(Count, 0),
    Percent = replace_na(Percent, 0)
  )

plot_page4 <- make_pie(location_boot_long, facet_col = "Region",
                       title_text = "Per region (bootstrap by Location): mean species composition")

# ---------- Exportar a PDF (4 páginas) ----------
pdf("Species_pies_Fig1.pdf", width = 10, height = 7)
print(plot_page1)  # Página 1
print(plot_page2)  # Página 2
print(plot_page3)  # Página 3
print(plot_page4)  # Página 4
dev.off()

#### Guardar CSVs de tablas relevantes  ######
# ====== Guardar tablas en CSV ======
# 1) Resumen global por región
write.csv(region_summary, "region_summary.csv", row.names = FALSE)
# 2) Tabla de crioviales con más de una secuencia
write.csv(multi_sequence_objids, "multi_sequence_objids.csv", row.names = FALSE)
# 3) Tabla de locations con más de una especie
write.csv(multi_species_locations_table, "multi_species_locations_table.csv", row.names = FALSE)
# 4) Bootstraps
write.csv(summary_objid_boot, "summary_objid_boot.csv", row.names = FALSE)
write.csv(summary_location_boot, "summary_location_boot.csv", row.names = FALSE)

