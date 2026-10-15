# ============================================================
# Figure 3&4. Detection across monitoring techniques
# ============================================================
#
# Run this script from the repository root.

# Packages
library(tidyverse)
library(ggparliament)

# Run from the repository root.
dir.create("outputs", showWarnings = FALSE)

mol <- read_csv("data/molecular_occ.csv", show_col_types = FALSE)
seine <- read_csv("data/seine_2024.csv", show_col_types = FALSE)
trawl <- read_csv("data/trawl_2024.csv", show_col_types = FALSE)

# Assign molecular samples to major bay systems
mol <- mol %>%
  mutate(
    area = case_when(
      str_starts(sample_ID_2, "AB") ~ "Aransas Bay System",
      str_starts(sample_ID_2, "BB") ~ "Upper Laguna Madre System",
      str_starts(sample_ID_2, "CB") ~ "Aransas Bay System",
      str_starts(sample_ID_2, "KE") ~ "Matagorda Bay System",
      str_starts(sample_ID_2, "AY") ~ "San Antonio Bay System",
      str_starts(sample_ID_2, "CA") ~ "Matagorda Bay System",
      str_starts(sample_ID_2, "TP") ~ "Matagorda Bay System",
      str_starts(sample_ID_2, "CR") ~ "Aransas Bay System",
      str_starts(sample_ID_2, "SC") ~ "Aransas Bay System",
      TRUE ~ NA_character_
    )
  )

# Collapse any technical replicate information to one result per environmental sample
mol_sample <- mol %>%
  group_by(sample_ID_2, sample_ID, area, year) %>%
  summarise(
    qPCR = ifelse(all(is.na(qPCR)), NA_integer_, as.integer(any(qPCR == 1, na.rm = TRUE))),
    ddPCR = ifelse(all(is.na(ddPCR)), NA_integer_, as.integer(any(ddPCR == 1, na.rm = TRUE))),
    metabarcoding = ifelse(
      all(is.na(metabarcoding)),
      NA_integer_,
      as.integer(any(metabarcoding == 1, na.rm = TRUE))
    ),
    .groups = "drop"
  )

mol_long <- mol_sample %>%
  pivot_longer(
    cols = c(qPCR, ddPCR, metabarcoding),
    names_to = "method",
    values_to = "occurrence"
  ) %>%
  filter(!is.na(occurrence)) %>%
  mutate(
    method = recode(
      method,
      qPCR = "qPCR",
      ddPCR = "ddPCR",
      metabarcoding = "Metabarcoding"
    ),
    status = case_when(
      occurrence == 0 ~ "Not detected",
      occurrence == 1 & year == 2024 ~ "Detected in 2024",
      occurrence == 1 & year == 2025 ~ "Detected in 2025"
    )
  )

# Fisheries-independent comparison is restricted to February and August 2024,
# matching the species-specific molecular sampling periods.
traditional_long <- bind_rows(
  seine %>% mutate(gear = "Seine"),
  trawl %>% mutate(gear = "Trawl")
) %>%
  filter(MONTH %in% c(2, 8)) %>%
  mutate(
    method = "Traditional fisheries",
    year = 2024,
    occurrence = as.integer(OCC),
    status = case_when(
      occurrence == 1 ~ "Detected in 2024",
      occurrence == 0 ~ "Not detected"
    )
  ) %>%
  dplyr::select(method, occurrence, status, gear, area, year)

plot_df <- bind_rows(
  mol_long %>% dplyr::select(method, occurrence, status, area, year),
  traditional_long %>% dplyr::select(method, occurrence, status, area, year)
) %>%
  mutate(
    method = factor(
      method,
      levels = c("Traditional fisheries", "Metabarcoding", "qPCR", "ddPCR")
    ),
    status = factor(
      status,
      levels = c("Not detected", "Detected in 2024", "Detected in 2025")
    )
  )


# Figure 3: semicircle detection plot

# Count samples in each detection-status category. Each sample will become
# one point (a "seat") in the corresponding method's semicircle.
parl_counts <- plot_df %>%
  group_by(method, status) %>%
  summarise(seats = n(), .groups = "drop")

# Calculate the totals used in each facet heading
method_summary <- plot_df %>%
  group_by(method) %>%
  summarise(
    n = n(),
    detections = sum(occurrence == 1, na.rm = TRUE),
    detection_rate = detections / n,
    .groups = "drop"
  ) %>%
  mutate(
    method_label = paste0(
      method, "\n", detections, "/", n,
      " (", round(100 * detection_rate, 1), "%)"
    ),
    method_label = factor(method_label, levels = method_label[order(method)])
  )

# ggparliament::parliament_data() creates one semicircle at a time, so this
# helper converts the counts for one detection method into seat coordinates.
make_parliament_semicircle <- function(df_method) {
  method_name <- unique(df_method$method)
  
  # Rename the status columns to the names expected by ggparliament.
  tmp <- df_method %>%
    ungroup() %>%
    filter(seats > 0) %>%
    mutate(
      party_short = as.character(status),
      party_long = as.character(status),
      seats = as.integer(seats)
    ) %>%
    dplyr::select(party_short, party_long, seats)
  
  # Use fewer rows for methods with fewer samples so their points remain
  # large enough to see. This changes only the layout, not the counts.
  n_seats <- sum(tmp$seats)
  n_rows <- case_when(
    n_seats <= 50 ~ 4,
    n_seats <= 150 ~ 6,
    TRUE ~ 8
  )
  
  ggparliament::parliament_data(
    election_data = tmp,
    party_seats = tmp$seats,
    parl_rows = n_rows,
    type = "semicircle"
  ) %>%
    mutate(method = method_name)
}

# Build one semicircle per method, combine their coordinates, and attach the formatted labels
parl_df <- parl_counts %>%
  group_split(method) %>%
  map_dfr(make_parliament_semicircle) %>%
  left_join(method_summary, by = "method") %>%
  mutate(
    method_label = factor(
      method_label,
      levels = levels(method_summary$method_label)
    )
  )

# Plot one point per sample, colored by its detection status, with one facet
# for each monitoring or molecular method.
Figure3 <- ggplot(parl_df, aes(x = x, y = y, colour = party_short)) +
  geom_parliament_seats(size = 1.9) +
  facet_wrap(~method_label, ncol = 4) +
  scale_colour_manual(
    values = c(
      "Detected in 2024" = "#0072B2",
      "Detected in 2025" = "#F5891D",
      "Not detected" = "grey78"
    ),
    name = NULL
  ) +
  coord_equal() +
  theme_void(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

Figure3

ggsave(
  "outputs/Figure3_detection_across_techniques.pdf",
  Figure3,
  width = 12,
  height = 3,
  units = "in"
)

# Figure 4: detection-rate heat map by bay system

area_summary <- plot_df %>%
  filter(!is.na(area)) %>%
  group_by(method, area) %>%
  summarise(
    n = n(),
    detections = sum(occurrence == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    detection_rate = detections / n,
    label = paste0(
      detections, "/", n,
      "\n(", scales::percent(detection_rate, accuracy = 0.1), ")"
    ),
    method = factor(
      method,
      levels = c("Traditional fisheries", "Metabarcoding", "qPCR", "ddPCR")
    ),
    area = factor(
      area,
      levels = c(
        "Upper Laguna Madre System",
        "Aransas Bay System",
        "San Antonio Bay System",
        "Matagorda Bay System"
      )
    )
  )

Figure4 <- ggplot(area_summary, aes(x = method, y = area, fill = detection_rate)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(
    aes(
      label = label,
      color = ifelse(detection_rate >= 0.15, "white", "black")
    ),
    size = 4
  ) +
  scale_color_identity() +
  scale_fill_gradient(
    low = "grey95",
    high = "grey30",
    labels = scales::percent_format(),
    name = "Detection\nrate"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_blank()
  )

Figure4

ggsave(
  "outputs/Figure4_detection_heatmap.pdf",
  Figure4,
  width = 8,
  height = 4,
  units = "in"
)
