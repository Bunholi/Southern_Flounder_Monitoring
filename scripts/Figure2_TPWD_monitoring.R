# ============================================================
# Figure 2. TPWD fisheries-independent monitoring
# ============================================================
#
# This script reproduces the four panels of Figure 2:
#   A. Bag-seine annual CPUE trend
#   B. Bay-trawl annual CPUE trend
#   C. Gill-net annual CPUE trend
#   D. Detection probability by sampling gear
#
# Run this script from the repository root.

# Packages

library(dplyr)
library(emmeans)
library(ggplot2)
library(mgcv)
library(patchwork)
library(scales)
library(stringr)


# Import and prepare TPWD datasets

area_labels <- c(
  "3" = "Matagorda Bay System",
  "4" = "San Antonio Bay System",
  "5" = "Aransas Bay System",
  "7" = "Upper Laguna Madre System"
)

prepare_tpwd_data <- function(path) {
  read.csv(path) %>%
    mutate(
      YEAR = as.numeric(YEAR),
      MONTH = as.numeric(MONTH),
      MAJOR_AREA_CODE = as.numeric(MAJOR_AREA_CODE),
      TOTAL_SAMPLED = as.numeric(TOTAL_SAMPLED),
      area = recode(as.character(MAJOR_AREA_CODE), !!!area_labels),
      area = factor(area, levels = unname(area_labels)),
      MONTH_F = factor(MONTH)
    )
}

df.s <- prepare_tpwd_data("data/TPWD_BS_SouthernFlounder_CPUE.csv")
df.t <- prepare_tpwd_data("data/TPWD_BT_SouthernFlounder_CPUE.csv")
df.g <- prepare_tpwd_data("data/TPWD_GN_SouthernFlounder_CPUE.csv")



# Color palette and plotting theme

pal <- c(
  "Matagorda Bay System" = "#566E3D",
  "San Antonio Bay System" = "#03256C",
  "Aransas Bay System" = "#7FCDBB",
  "Upper Laguna Madre System" = "#2C7FB8"
)

theme_flounder <- theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "black", size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12)
  )


# Bag-seine CPUE model

gam_year.s <- gam(
  CATCH_PER_HECTARE ~ s(YEAR, k = 10) + area +
    s(MONTH, bs = "cc", k = 10),
  data = df.s,
  family = tw(link = "log"),
  method = "REML"
)

newdat_year.s <- expand.grid(
  YEAR = seq(min(df.s$YEAR, na.rm = TRUE),
             max(df.s$YEAR, na.rm = TRUE), by = 0.1),
  MONTH = 1:12,
  area = levels(df.s$area)
)

p.s <- predict(gam_year.s, newdata = newdat_year.s,
               type = "link", se.fit = TRUE)

pred_year.s <- newdat_year.s %>%
  mutate(
    fit = exp(p.s$fit),
    lwr = exp(p.s$fit - 2 * p.s$se.fit),
    upr = exp(p.s$fit + 2 * p.s$se.fit)
  ) %>%
  group_by(YEAR, area) %>%
  summarise(
    predicted = mean(fit, na.rm = TRUE),
    conf.low = mean(lwr, na.rm = TRUE),
    conf.high = mean(upr, na.rm = TRUE),
    .groups = "drop"
  )

seine.yr <- ggplot(
  pred_year.s,
  aes(x = YEAR, y = predicted, color = area, fill = area)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal, name = "Bay system") +
  scale_fill_manual(values = pal, name = "Bay system") +
  guides(color = "none", fill = "none") +
  labs(
    title = "Bag seine",
    x = "Year",
    y = "CPUH (catch per hectare)"
  ) +
  theme_flounder


# Bay-trawl CPUE model

gam_year.t <- gam(
  CATCH_PER_HOUR ~ s(YEAR, k = 10) + area +
    s(MONTH, bs = "cc", k = 10),
  data = df.t,
  family = tw(link = "log"),
  method = "REML"
)

newdat_year.t <- expand.grid(
  YEAR = seq(min(df.t$YEAR, na.rm = TRUE),
             max(df.t$YEAR, na.rm = TRUE), by = 0.1),
  MONTH = 1:12,
  area = levels(df.t$area)
)

p.t <- predict(gam_year.t, newdata = newdat_year.t,
               type = "link", se.fit = TRUE)

pred_year.t <- newdat_year.t %>%
  mutate(
    fit = exp(p.t$fit),
    lwr = exp(p.t$fit - 2 * p.t$se.fit),
    upr = exp(p.t$fit + 2 * p.t$se.fit)
  ) %>%
  group_by(YEAR, area) %>%
  summarise(
    predicted = mean(fit, na.rm = TRUE),
    conf.low = mean(lwr, na.rm = TRUE),
    conf.high = mean(upr, na.rm = TRUE),
    .groups = "drop"
  )

trawl.yr <- ggplot(
  pred_year.t,
  aes(x = YEAR, y = predicted, color = area, fill = area)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal, name = "Bay system") +
  scale_fill_manual(values = pal, name = "Bay system") +
  guides(color = "none", fill = "none") +
  labs(
    title = "Bay trawl",
    x = "Year",
    y = "CPUE (catch per hour)"
  ) +
  theme_flounder


# Gill-net CPUE model

gam_year.g <- gam(
  CATCH_PER_HOUR ~ s(YEAR, k = 10) + area + MONTH_F,
  data = df.g,
  family = tw(link = "log"),
  method = "REML"
)

newdat_year.g <- expand.grid(
  YEAR = seq(min(df.g$YEAR, na.rm = TRUE),
             max(df.g$YEAR, na.rm = TRUE), by = 0.1),
  MONTH_F = levels(df.g$MONTH_F),
  area = levels(df.g$area)
)

p.g <- predict(gam_year.g, newdata = newdat_year.g,
               type = "link", se.fit = TRUE)

pred_year.g <- newdat_year.g %>%
  mutate(
    fit = exp(p.g$fit),
    lwr = exp(p.g$fit - 2 * p.g$se.fit),
    upr = exp(p.g$fit + 2 * p.g$se.fit)
  ) %>%
  group_by(YEAR, area) %>%
  summarise(
    predicted = mean(fit, na.rm = TRUE),
    conf.low = mean(lwr, na.rm = TRUE),
    conf.high = mean(upr, na.rm = TRUE),
    .groups = "drop"
  )

gillnet.yr <- ggplot(
  pred_year.g,
  aes(x = YEAR, y = predicted, color = area, fill = area)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal, name = "Bay system") +
  scale_fill_manual(values = pal, name = "Bay system") +
  guides(color = "none", fill = "none") +
  labs(
    title = "Gill net",
    x = "Year",
    y = "CPUE (catch per hour)"
  ) +
  theme_flounder


# Detection-probability analysis

df.clean.s <- df.s %>%
  filter(
    START_SALINITY_NUM <= 60,
    START_TURBIDITY_NUM <= 300
  ) %>%
  mutate(gear = "Seine")

df.clean.t <- df.t %>%
  filter(
    START_SALINITY_NUM <= 60,
    START_TURBIDITY_NUM <= 300
  ) %>%
  mutate(gear = "Trawl")

df.clean.g <- df.g %>%
  filter(
    START_SALINITY_NUM <= 60,
    START_TURBIDITY_NUM <= 300
  ) %>%
  mutate(gear = "Gillnet")

df.all <- bind_rows(df.clean.s, df.clean.t, df.clean.g) %>%
  mutate(
    area = str_squish(as.character(area)),
    area = factor(area, levels = names(pal)),
    gear = factor(gear, levels = c("Seine", "Trawl", "Gillnet")),
    presence = as.integer(TOTAL_SAMPLED > 0)
  )

#model
detection_model <- gam(
  presence ~ gear + area + s(MONTH) + s(YEAR),
  data = df.all,
  family = binomial(link = "logit"),
  method = "REML"
)

# Figure 2d 
det_summary <- df.all %>%
  group_by(gear, area) %>%
  summarise(
    n = n(),
    detections = sum(presence == 1, na.rm = TRUE),
    p = mean(presence == 1, na.rm = TRUE),
    se = sqrt(p * (1 - p) / n),
    lwr = pmax(0, p - 1.96 * se),
    upr = pmin(1, p + 1.96 * se),
    .groups = "drop"
  ) %>%
  mutate(area = factor(area, levels = names(pal)))

gear_type <- ggplot(det_summary, aes(x = gear, y = p, fill = area)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  scale_fill_manual(values = pal, name = "System") +
  labs(
    x = "Gear type",
    y = "Detection probability"
  ) +
  theme_classic(base_size = 14)

gear_type

# Assemble and save Figure 2

Figure2 <- (
  (seine.yr | trawl.yr) /
    (gillnet.yr | gear_type)
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.tag = element_text(face = "bold", size = 13)
  )

Figure2

ggsave(
  filename = "/Users/ingridbunholi/Downloads/Southern_Flounder_Monitoring_GitHub/data/Figure2_TPWD_monitoring.pdf",
  plot = Figure2,
  width = 12,
  height = 8,
  units = "in"
)


# Model outputs for manuscript tables

summary(gam_year.s)
summary(gam_year.t)
summary(gam_year.g)
summary(detection_model)

nobs(gam_year.s)
nobs(gam_year.t)
nobs(gam_year.g)
nobs(detection_model)

emmeans(detection_model, pairwise ~ gear, type = "response")
