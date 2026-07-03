library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)

target_species <- c(
  "Yellowfin Bream (Acanthopagrus australis)"
)

# Set your video frame rate
fps <- 60

# 60-second interval size in frames
interval_seconds <- 60
interval_frames <- fps * interval_seconds

# 1. MeanCount using 60-second intervals
meancount <- MxMd %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_"),
    TimeBin_60s = floor(Frame / interval_frames)
  ) %>%
  group_by(CameraID, Estuary, Habitat, Site, Tide, Camera, Species, TimeBin_60s) %>%
  summarise(
    interval_count = n(),
    interval_mean_size = mean(`Size (cm)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(CameraID, Estuary, Habitat, Site, Tide, Camera, Species) %>%
  summarise(
    MeanCount = mean(interval_count, na.rm = TRUE),
    MeanCount_mean_size = mean(interval_mean_size, na.rm = TRUE),
    n_intervals = n(),
    .groups = "drop"
  )

# 2. MIDS abundance and mean size per deployment
mids <- MxMd2 %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_")
  ) %>%
  group_by(CameraID, Estuary, Habitat, Site, Tide, Camera, Species) %>%
  summarise(
    MIDS = n(),
    MIDS_mean_size = mean(`Size (cm)`, na.rm = TRUE),
    .groups = "drop"
  )
comparison_df <- meancount %>%
  inner_join(
    mids,
    by = c(
      "CameraID",
      "Estuary",
      "Habitat",
      "Site",
      "Tide",
      "Camera",
      "Species"
    )
  )

wilcox.test(
  comparison_df$MIDS,
  comparison_df$MeanCount,
  paired = TRUE,
  exact = FALSE
)
comparison_df %>%
  summarise(
    Median_MeanCount = median(MeanCount),
    Median_MIDS = median(MIDS),
    Median_Difference = median(MIDS - MeanCount)
  )
comparison_long <- comparison_df %>%
  select(CameraID, MeanCount, MIDS) %>%
  pivot_longer(
    cols = c(MeanCount, MIDS),
    names_to = "Method",
    values_to = "Abundance"
  )

comparison_long %>%
  wilcox_test(
    Abundance ~ Method,
    paired = TRUE
  )

comparison_long %>%
  wilcox_effsize(
    Abundance ~ Method,
    paired = TRUE
  )

summary_data <- comparison_long %>%
  group_by(Method) %>%
  summarise(
    median_abundance = median(Abundance),
    q1 = quantile(Abundance, 0.25),
    q3 = quantile(Abundance, 0.75)
  )


median_iqr <- function(x) {
  data.frame(
    y = median(x),
    ymin = quantile(x, 0.25),
    ymax = quantile(x, 0.75)
  )
}

ggplot(comparison_long, aes(Method, Abundance)) +

  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.7,
    colour = "grey30"
  ) +

  stat_summary(
    fun.data = median_iqr,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.8
  ) +

  stat_summary(
    fun = median,
    geom = "point",
    shape = 18,     # diamond
    size = 4,
    colour = "orange"
  ) +

  theme_classic()


# 3. Join
comparison <- mids %>%
  full_join(
    meancount,
    by = c(
      "CameraID", "Estuary", "Habitat",
      "Site", "Tide", "Camera", "Species"
    )
  ) %>%
  mutate(
    abundance_difference = MIDS_n - MeanCount,
    size_difference = MIDS_mean_size - MeanCount_mean_size,
    abundance_ratio = MIDS_n / MeanCount
  )

comparison

# 4. Abundance long format
abundance_plot <- comparison %>%
  select(
    CameraID, Estuary, Habitat, Site, Tide, Camera, Species,
    MIDS_n, MeanCount
  ) %>%
  pivot_longer(
    cols = c(MIDS_n, MeanCount),
    names_to = "Method",
    values_to = "Abundance"
  ) %>%
  mutate(
    Method = recode(
      Method,
      MIDS_n = "MIDS",
      MeanCount = "MeanCount"
    ),
    Method = factor(Method, levels = c("MeanCount", "MIDS"))
  ) %>%
  filter(!is.na(Abundance))

# 5. Mean size long format
size_plot <- comparison %>%
  select(
    CameraID, Estuary, Habitat, Site, Tide, Camera, Species,
    MIDS_mean_size, MeanCount_mean_size
  ) %>%
  pivot_longer(
    cols = c(MIDS_mean_size, MeanCount_mean_size),
    names_to = "Method",
    values_to = "Mean_size"
  ) %>%
  mutate(
    Method = recode(
      Method,
      MIDS_mean_size = "MIDS",
      MeanCount_mean_size = "MeanCount"
    ),
    Method = factor(Method, levels = c("MeanCount", "MIDS"))
  ) %>%
  filter(!is.na(Mean_size))

# 6. Abundance plot
ggplot(abundance_plot,
       aes(x = Method,
           y = Abundance,
           fill = Method)) +

  geom_jitter(
    width = 0.15,
    alpha = 0.5
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    shape = 18,
    colour = "orange"
  ) +

  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.2
  ) +

  facet_grid(
    Estuary ~ Habitat,
    scales = "free_y"
  ) +

  labs(
    x = NULL,
    y = "Abundance",
    fill = "Method"
  ) +

  theme_classic()

# 7. Mean size plot
ggplot(size_plot,
       aes(x = Method,
           y = Mean_size,
           fill = Method)) +

  geom_jitter(
    width = 0.15,
    alpha = 0.5
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    shape = 18,
    colour = "orange"
  ) +

  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.2
  ) +

  facet_grid(
    Estuary ~ Habitat,
    scales = "free_y"
  ) +

  labs(
    x = NULL,
    y = "Mean length per 60-second interval (cm)",
    fill = "Method"
  ) +

  theme_classic()

# 8. Overall abundance test
abundance_test_overall <- abundance_plot %>%
  group_by(Species) %>%
  wilcox_test(
    Abundance ~ Method,
    paired = TRUE
  )

abundance_test_overall

# 9. Overall abundance effect size
abundance_effect_overall <- abundance_plot %>%
  group_by(Species) %>%
  wilcox_effsize(
    Abundance ~ Method,
    paired = TRUE
  )

abundance_effect_overall

# 10. Abundance tests by Estuary x Habitat
abundance_tests_by_group <- abundance_plot %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_test(
    Abundance ~ Method,
    paired = TRUE
  ) %>%
  adjust_pvalue(method = "BH")

abundance_tests_by_group

# 11. Abundance effect sizes by Estuary x Habitat
abundance_effects_by_group <- abundance_plot %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_effsize(
    Abundance ~ Method,
    paired = TRUE
  )

abundance_effects_by_group

# 12. Summary table
summary_table <- comparison %>%
  group_by(Species, Estuary, Habitat) %>%
  summarise(
    n_deployments = n(),
    median_MIDS_n = median(MIDS_n, na.rm = TRUE),
    median_MeanCount = median(MeanCount, na.rm = TRUE),
    mean_MIDS_n = mean(MIDS_n, na.rm = TRUE),
    mean_MeanCount = mean(MeanCount, na.rm = TRUE),
    median_ratio = median(abundance_ratio, na.rm = TRUE),
    mean_ratio = mean(abundance_ratio, na.rm = TRUE),
    median_MIDS_size = median(MIDS_mean_size, na.rm = TRUE),
    median_MeanCount_size = median(MeanCount_mean_size, na.rm = TRUE),
    .groups = "drop"
  )

summary_table
target_species <- c(
  "Yellowfin Bream (Acanthopagrus australis)"
)

# MeanCount abundance and mean size per deployment
meancount <- MxMd %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_")
  ) %>%
  group_by(CameraID, Estuary, Habitat, Species, Frame) %>%
  summarise(
    frame_count = n(),
    frame_mean_size = mean(`Size (cm)`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(CameraID, Estuary, Habitat, Species) %>%
  summarise(
    MeanCount = mean(frame_count, na.rm = TRUE),
    MeanCount_mean_size = mean(frame_mean_size, na.rm = TRUE),
    n_frames = n(),
    .groups = "drop"
  )

# MIDS abundance and mean size per deployment
mids <- MxMd2 %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_")
  ) %>%
  group_by(CameraID, Estuary, Habitat, Species) %>%
  summarise(
    MIDS_n = n(),
    MIDS_mean_size = mean(`Size (cm)`, na.rm = TRUE),
    .groups = "drop"
  )

# Join MIDS and MeanCount
comparison <- mids %>%
  full_join(
    meancount,
    by = c("CameraID", "Estuary", "Habitat", "Species")
  ) %>%
  mutate(
    abundance_difference = MIDS_n - MeanCount,
    size_difference = MIDS_mean_size - MeanCount_mean_size
  )

comparison

# Abundance long format
abundance_plot <- comparison %>%
  select(CameraID, Estuary, Habitat, Species, MIDS_n, MeanCount) %>%
  pivot_longer(
    cols = c(MIDS_n, MeanCount),
    names_to = "Method",
    values_to = "Abundance"
  ) %>%
  mutate(
    Method = recode(
      Method,
      MIDS_n = "MIDS",
      MeanCount = "MeanCount"
    ),
    Method = factor(Method, levels = c("MeanCount", "MIDS"))
  ) %>%
  filter(!is.na(Abundance))

# Mean size long format
size_plot <- comparison %>%
  select(
    CameraID,
    Estuary,
    Habitat,
    Species,
    MIDS_mean_size,
    MeanCount_mean_size
  ) %>%
  pivot_longer(
    cols = c(MIDS_mean_size, MeanCount_mean_size),
    names_to = "Method",
    values_to = "Mean_size"
  ) %>%
  mutate(
    Method = recode(
      Method,
      MIDS_mean_size = "MIDS",
      MeanCount_mean_size = "MeanCount"
    ),
    Method = factor(Method, levels = c("MeanCount", "MIDS"))
  ) %>%
  filter(!is.na(Mean_size))

# Abundance plot
ggplot(abundance_plot,
       aes(x = Method,
           y = Abundance,
           fill = Method)) +

  geom_jitter(
    width = 0.15,
    alpha = 0.5
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    shape = 18,
    colour = "orange"
  ) +

  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.2
  ) +

  facet_grid(
    Estuary ~ Habitat,
    scales = "free_y"
  ) +

  labs(
    x = NULL,
    y = "Abundance",
    fill = "Method"
  ) +

  theme_classic()

# Mean size plot
ggplot(size_plot,
       aes(x = Method,
           y = Mean_size,
           fill = Method)) +

  geom_jitter(
    width = 0.15,
    alpha = 0.5
  ) +

  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    shape = 18,
    colour = "orange"
  ) +

  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    width = 0.2
  ) +

  facet_grid(
    Estuary ~ Habitat,
    scales = "free_y"
  ) +

  labs(
    x = NULL,
    y = "Mean length per deployment (cm)",
    fill = "Method"
  ) +

  theme_classic()

# Overall abundance Wilcoxon test
abundance_test_overall <- abundance_plot %>%
  group_by(Species) %>%
  wilcox_test(
    Abundance ~ Method,
    paired = TRUE
  )

abundance_test_overall

# Overall abundance effect size
abundance_effect_overall <- abundance_plot %>%
  group_by(Species) %>%
  wilcox_effsize(
    Abundance ~ Method,
    paired = TRUE
  )

abundance_effect_overall

# Abundance Wilcoxon tests by Estuary x Habitat
abundance_tests_by_group <- abundance_plot %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_test(
    Abundance ~ Method,
    paired = TRUE
  ) %>%
  adjust_pvalue(method = "BH")

abundance_tests_by_group

# Abundance effect sizes by Estuary x Habitat
abundance_effects_by_group <- abundance_plot %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_effsize(
    Abundance ~ Method,
    paired = TRUE
  )

abundance_effects_by_group

# Overall mean size Wilcoxon test
size_test_overall <- size_plot %>%
  group_by(Species) %>%
  wilcox_test(
    Mean_size ~ Method,
    paired = TRUE
  )

size_test_overall

# Overall mean size effect size
size_effect_overall <- size_plot %>%
  group_by(Species) %>%
  wilcox_effsize(
    Mean_size ~ Method,
    paired = TRUE
  )

size_effect_overall

# Mean size Wilcoxon tests by Estuary x Habitat
size_tests_by_group <- size_plot %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_test(
    Mean_size ~ Method,
    paired = TRUE
  ) %>%
  adjust_pvalue(method = "BH")

size_tests_by_group

# Mean size effect sizes by Estuary x Habitat
size_effects_by_group <- size_plot %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_effsize(
    Mean_size ~ Method,
    paired = TRUE
  )

size_effects_by_group

# Summary table
summary_table <- comparison %>%
  group_by(Species, Estuary, Habitat) %>%
  summarise(
    n_deployments = n(),
    median_MIDS_n = median(MIDS_n, na.rm = TRUE),
    median_MeanCount = median(MeanCount, na.rm = TRUE),
    mean_MIDS_n = mean(MIDS_n, na.rm = TRUE),
    mean_MeanCount = mean(MeanCount, na.rm = TRUE),
    median_MIDS_size = median(MIDS_mean_size, na.rm = TRUE),
    median_MeanCount_size = median(MeanCount_mean_size, na.rm = TRUE),
    mean_MIDS_size = mean(MIDS_mean_size, na.rm = TRUE),
    mean_MeanCount_size = mean(MeanCount_mean_size, na.rm = TRUE),
    .groups = "drop"
  )

summary_table

#frequency
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)

target_species <- c(
  "Yellowfin Bream (Acanthopagrus australis)"
)

# 1. MeanCount mean length per frame
meancount_lengths <- MxMd %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_")
  ) %>%
  group_by(CameraID, Estuary, Habitat, Site, Tide, Camera, Species, Frame) %>%
  summarise(
    `Size (cm)` = mean(`Size (cm)`, na.rm = TRUE),
    frame_count = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Method = "MeanCount"
  ) %>%
  select(
    CameraID, Estuary, Habitat, Site, Tide, Camera,
    Species, Frame, Method, `Size (cm)`, frame_count
  ) %>%
  filter(!is.na(`Size (cm)`))

# 2. MIDS length observations
mids_lengths <- MxMd2 %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_"),
    Method = "MIDS"
  ) %>%
  select(
    CameraID, Estuary, Habitat, Site, Tide, Camera,
    Species, Method, `Size (cm)`
  ) %>%
  filter(!is.na(`Size (cm)`))

# 3. Combine
length_freq_data <- bind_rows(
  meancount_lengths,
  mids_lengths
) %>%
  mutate(
    Method = factor(Method, levels = c("MeanCount", "MIDS"))
  )

# Check sample sizes overall
length_freq_data %>%
  count(Method)

# Check sample sizes by Estuary x Habitat
length_freq_data %>%
  count(Estuary, Habitat, Method)

# 4. Create length bins
length_freq_data <- length_freq_data %>%
  mutate(
    Length_bin = cut(
      `Size (cm)`,
      breaks = seq(
        floor(min(`Size (cm)`, na.rm = TRUE)),
        ceiling(max(`Size (cm)`, na.rm = TRUE)) + 2,
        by = 2
      ),
      include.lowest = TRUE
    )
  )

# 5. Summarise length-frequency counts
length_freq_summary <- length_freq_data %>%
  group_by(Species, Estuary, Habitat, Method, Length_bin) %>%
  summarise(
    n = n(),
    .groups = "drop"
  )

length_freq_summary

# 6. Length-frequency plot by Estuary x Habitat
ggplot(length_freq_data,
       aes(x = `Size (cm)`,
           colour = Method)) +

  geom_histogram(
    binwidth = 2,
    fill = NA,
    linewidth = 1,
    position = "identity"
  ) +

  scale_colour_manual(
    values = c(
      "MeanCount" = "black",
      "MIDS" = "orange"
    ) ) +

  labs(
    x = "Fish length (cm)",
    y = "Number of observations",
    colour = "Method"
  ) +

  theme_classic()

# 7. Overall Wilcoxon test
length_wilcox_overall <- length_freq_data %>%
  group_by(Species) %>%
  wilcox_test(
    `Size (cm)` ~ Method
  )

length_wilcox_overall

# 8. Overall Fligner-Killeen test
length_fligner_overall <- length_freq_data %>%
  group_by(Species) %>%
  summarise(
    fligner_statistic = fligner.test(`Size (cm)` ~ Method)$statistic,
    fligner_p = fligner.test(`Size (cm)` ~ Method)$p.value,
    .groups = "drop"
  )

length_fligner_overall

# 9. Overall KS test
length_ks_overall <- length_freq_data %>%
  group_by(Species) %>%
  summarise(
    ks_statistic = ks.test(
      `Size (cm)`[Method == "MeanCount"],
      `Size (cm)`[Method == "MIDS"]
    )$statistic,
    ks_p = ks.test(
      `Size (cm)`[Method == "MeanCount"],
      `Size (cm)`[Method == "MIDS"]
    )$p.value,
    .groups = "drop"
  )

length_ks_overall

# 10. Wilcoxon tests by Estuary x Habitat
length_wilcox_by_group <- length_freq_data %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_test(
    `Size (cm)` ~ Method
  ) %>%
  adjust_pvalue(method = "BH")

length_wilcox_by_group

# 11. Effect sizes by Estuary x Habitat
length_effects_by_group <- length_freq_data %>%
  group_by(Species, Estuary, Habitat) %>%
  wilcox_effsize(
    `Size (cm)` ~ Method
  )

length_effects_by_group

# 12. Fligner-Killeen tests by Estuary x Habitat
length_fligner_by_group <- length_freq_data %>%
  group_by(Species, Estuary, Habitat) %>%
  summarise(
    fligner_statistic = fligner.test(`Size (cm)` ~ Method)$statistic,
    fligner_p = fligner.test(`Size (cm)` ~ Method)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    fligner_p_adj = p.adjust(fligner_p, method = "BH")
  )

length_fligner_by_group

# 13. KS tests by Estuary x Habitat
length_ks_by_group <- length_freq_data %>%
  group_by(Species, Estuary, Habitat) %>%
  summarise(
    ks_statistic = ks.test(
      `Size (cm)`[Method == "MeanCount"],
      `Size (cm)`[Method == "MIDS"]
    )$statistic,
    ks_p = ks.test(
      `Size (cm)`[Method == "MeanCount"],
      `Size (cm)`[Method == "MIDS"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    ks_p_adj = p.adjust(ks_p, method = "BH")
  )

length_ks_by_group

# 14. Summary table
length_summary <- length_freq_data %>%
  group_by(Species, Estuary, Habitat, Method) %>%
  summarise(
    n = n(),
    mean_length = mean(`Size (cm)`, na.rm = TRUE),
    median_length = median(`Size (cm)`, na.rm = TRUE),
    sd_length = sd(`Size (cm)`, na.rm = TRUE),
    min_length = min(`Size (cm)`, na.rm = TRUE),
    max_length = max(`Size (cm)`, na.rm = TRUE),
    .groups = "drop"
  )

length_summary

# 1. MeanCount length composition: pooled individual lengths from 60-s intervals
meancount_lengths <- MxMd %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_"),
    TimeBin_60s = floor(Frame / interval_frames)
  ) %>%
  group_by(CameraID, Estuary, Habitat, Site, Tide, Camera, Species, TimeBin_60s) %>%
  mutate(
    # keeps one sampled frame per 60-s interval
    sampled_frame = min(Frame, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(Frame == sampled_frame) %>%
  mutate(Method = "MeanCount") %>%
  select(
    CameraID, Estuary, Habitat, Site, Tide, Camera,
    Species, Frame, TimeBin_60s, Method, `Size (cm)`
  ) %>%
  filter(!is.na(`Size (cm)`))

# 2. MIDS length observations
mids_lengths <- MxMd2 %>%
  filter(Species %in% target_species) %>%
  mutate(
    CameraID = paste(Estuary, Habitat, Site, Tide, Camera, sep = "_"),
    Method = "MIDS"
  ) %>%
  select(
    CameraID, Estuary, Habitat, Site, Tide, Camera,
    Species, Method, `Size (cm)`
  ) %>%
  filter(!is.na(`Size (cm)`))

# 3. Combine
length_freq_data <- bind_rows(
  meancount_lengths,
  mids_lengths
) %>%
  mutate(
    Method = factor(Method, levels = c("MeanCount", "MIDS"))
  )

# Check sample sizes
length_freq_data %>%
  count(Method)
length_freq_data <- length_freq_data %>%
  mutate(
    Length_bin = cut(
      `Size (cm)`,
      breaks = seq(
        floor(min(`Size (cm)`, na.rm = TRUE)),
        ceiling(max(`Size (cm)`, na.rm = TRUE)) + 2,
        by = 2
      ),
      include.lowest = TRUE
    )
  )

length_freq_summary <- length_freq_data %>%
  group_by(Species, Method, Length_bin) %>%
  summarise(n = n(), .groups = "drop")
library(ggplot2)
library(dplyr)

ggplot(length_freq_data,
       aes(x = `Size (cm)`, fill = Method)) +
  geom_histogram(
    binwidth = 2,
    position = "identity",
    alpha = 1,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "MeanCount" = "darkgray",
      "MIDS" = "#F4A300"
    )
  ) +
  labs(
    x = "Fish length (cm)",
    y = "Number of observations",
    fill = "Method"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = "right"
  )
