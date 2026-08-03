# Mahogany-plantations-and-woody-plant-diversity-Fiji
Plot-level data from Nukurua and Nausori Highlands

# =============================================================================
# Woody Plant Composition and Richness in Mahogany Plantations and Adjacent
# Native Forests in Fiji
#
# Authors: Iliesa Koroi, Ole R. Vetaas, Tamara Osborne-Naikatini
#
# This script reproduces all analyses, figures, tables, and appendices
# reported in the manuscript. Run sections sequentially; all objects
# created in earlier sections are required by later ones.
#
# Working directory must contain the eight CSV data files:
#   NKR_native.csv, NKR_stocked.csv, NKR_unstocked.csv, NKR_mixed.csv
#   NHL_native.csv, NHL_stocked.csv, NHL_unstocked.csv, NHL_mixed.csv
#
# Output files written to the working directory:
#   Figures 2-5 (.png and .tiff)
#   Table1_richness_abundance.csv
#   Table2_PERMANOVA_betadisper_CCA.csv
#   AppendixS2_GLM_results.csv
#   wilcoxon_results_all.csv
#   df_stats_summary.csv
#   seedling_frequency_by_group.csv
#   tree_frequency_by_group.csv
#   rank_abundance_seedlings_all.csv
#   sac_individual_based_all.csv
#   dca_scores_combined_intersite.csv
#   dca_scores_nkr_within_site.csv
#   dca_scores_nhl_within_site.csv
#
# R version: 4.4.1
# Key packages: vegan 2.7-3, ggplot2, cowplot, ggpubr, gt
# =============================================================================
rm(list = ls())

setwd("C:/Users/Administrator/Desktop/Rwf")

# ---- 0. Packages -------------------------------------------------------------
library(tidyverse)   # dplyr, tidyr, ggplot2, stringr, readr
library(vegan)       # adonis2, betadisper, cca, decorana, specaccum, vegdist
library(cowplot)     # plot_grid, ggdraw, draw_plot, draw_label, get_legend
library(ggpubr)      # stat_compare_means

# Close any stray graphics devices from previous runs
while (!is.null(dev.list())) dev.off()


# =============================================================================
# SECTION 0: Load and combine raw data
# =============================================================================

read_site <- function(site_prefix) {
  groups <- c("native", "stocked", "unstocked", "mixed")
  out <- lapply(groups, function(g) {
    fname <- paste0(site_prefix, "_", g, ".csv")
    # Force DBH to character at read time: some files store blanks or text
    # in the DBH column, which causes type conflicts when binding rows.
    df <- read.csv(fname, stringsAsFactors = FALSE,
                   colClasses = c(DBH = "character"))
    df$Group <- str_to_title(g)
    df
  })
  dplyr::bind_rows(out)
}

NKR_all <- read_site("NKR")
NHL_all <- read_site("NHL")

# Report non-numeric DBH rows before converting (surfaces data-entry errors)
report_bad_dbh <- function(df, site_label) {
  bad <- df %>%
    filter(is.na(suppressWarnings(as.numeric(DBH))) & !is.na(DBH) & DBH != "")
  if (nrow(bad) > 0) {
    cat(sprintf("\n%s: %d row(s) with non-numeric DBH:\n", site_label, nrow(bad)))
    print(bad %>% select(Group, Transect, Plot, Species, DBH))
  } else {
    cat(sprintf("%s: all DBH values numeric or blank. OK.\n", site_label))
  }
}
report_bad_dbh(NKR_all, "NKR")
report_bad_dbh(NHL_all, "NHL")

# Build unique PlotID: Plot (1-5) repeats across Transect (1-3) and Group,
# so a composite key is required to distinguish physical plots.
NKR_all <- NKR_all %>% mutate(PlotID = paste(Group, Transect, Plot, sep = "_"))
NHL_all <- NHL_all %>% mutate(PlotID = paste(Group, Transect, Plot, sep = "_"))

# Set factor levels for consistent ordering throughout
group_levels <- c("Native", "Stocked", "Unstocked", "Mixed")
NKR_all$Group <- factor(NKR_all$Group, levels = group_levels)
NHL_all$Group <- factor(NHL_all$Group, levels = group_levels)

# Convert DBH to numeric (any bad values become NA with a warning)
NKR_all$DBH <- as.numeric(NKR_all$DBH)
NHL_all$DBH <- as.numeric(NHL_all$DBH)


# =============================================================================
# SECTION 0b: Verify plot counts
# =============================================================================

check_plot_counts <- function(df, site_label) {
  counts <- df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(n_plots = dplyr::n_distinct(PlotID), .groups = "drop")
  cat(sprintf("\n%s plot counts:\n", site_label))
  print(counts)
  flagged <- counts %>% dplyr::filter(n_plots < 15)
  if (nrow(flagged) > 0) {
    cat("WARNING: groups missing plots:\n"); print(flagged)
  } else {
    cat("OK: all groups have 15 plots.\n")
  }
}

check_plot_counts(NKR_all, "NKR")
check_plot_counts(NHL_all, "NHL")


# =============================================================================
# SECTION 1: Per-plot species richness and seedling abundance
# =============================================================================

compute_plot_richness <- function(data_all) {
  
  # Total woody species per plot: union of trees + seedlings + climbers
  woody_long <- data_all %>%
    pivot_longer(cols = c(Species, Seedlings, Climbers),
                 names_to = "Layer", values_to = "WoodyTaxon") %>%
    filter(!is.na(WoodyTaxon), WoodyTaxon != "")
  
  woody <- woody_long %>%
    group_by(PlotID, Group) %>%
    summarise(TotalWoodySpecies = n_distinct(WoodyTaxon), .groups = "drop")
  
  # Seedling species richness per plot
  seedlings <- data_all %>%
    filter(!is.na(Seedlings), Seedlings != "") %>%
    group_by(PlotID, Group) %>%
    summarise(SeedlingRichness = n_distinct(Seedlings), .groups = "drop")
  
  # Seedling abundance per plot (sum of Abundance column)
  seedling_abundance <- data_all %>%
    filter(!is.na(Seedlings), Seedlings != "") %>%
    mutate(Abundance = as.numeric(Abundance)) %>%
    group_by(PlotID, Group) %>%
    summarise(SeedlingAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  list(woody             = woody,
       seedlings         = seedlings,
       seedling_abundance = seedling_abundance)
}

NKR_richness <- compute_plot_richness(NKR_all)
NHL_richness <- compute_plot_richness(NHL_all)

# Helper: standard error
se <- function(x) sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))

# Summary table (pre-zero-fill, for inspection only)
summarise_richness <- function(plot_counts, value_col) {
  plot_counts %>%
    group_by(Group) %>%
    summarise(Mean = mean(.data[[value_col]], na.rm = TRUE),
              SE   = se(.data[[value_col]]),
              n    = n(),
              .groups = "drop")
}

cat("\n=== NKR total woody: mean +/- SE ===\n")
print(summarise_richness(NKR_richness$woody, "TotalWoodySpecies"))
cat("\n=== NKR seedling richness: mean +/- SE ===\n")
print(summarise_richness(NKR_richness$seedlings, "SeedlingRichness"))
cat("\n=== NKR seedling abundance: mean +/- SE ===\n")
print(summarise_richness(NKR_richness$seedling_abundance, "SeedlingAbundance"))
cat("\n=== NHL total woody: mean +/- SE ===\n")
print(summarise_richness(NHL_richness$woody, "TotalWoodySpecies"))
cat("\n=== NHL seedling richness: mean +/- SE ===\n")
print(summarise_richness(NHL_richness$seedlings, "SeedlingRichness"))
cat("\n=== NHL seedling abundance: mean +/- SE ===\n")
print(summarise_richness(NHL_richness$seedling_abundance, "SeedlingAbundance"))


# =============================================================================
# SECTION 1b: Zero-fill missing plots
# =============================================================================
# Plots with zero seedlings recorded drop out of the richness/abundance
# groupby above. Zero-filling reinstates them as true ecological zeros
# before Wilcoxon tests and GLMs, so all n=15 per group.

zero_fill_metric <- function(all_data, metric_df, value_col) {
  all_plot_ids <- all_data %>% distinct(PlotID, Group)
  metric_df %>%
    right_join(all_plot_ids, by = c("PlotID", "Group")) %>%
    mutate(!!value_col := ifelse(is.na(.data[[value_col]]), 0,
                                 .data[[value_col]]))
}

zero_fill_all <- function(all_data, richness_list) {
  lapply(names(richness_list), function(nm) {
    df        <- richness_list[[nm]]
    value_col <- setdiff(names(df), c("PlotID", "Group"))[1]
    zero_fill_metric(all_data, df, value_col)
  }) %>% setNames(names(richness_list))
}

NKR_richness <- zero_fill_all(NKR_all, NKR_richness)
NHL_richness <- zero_fill_all(NHL_all, NHL_richness)

cat("\n=== Post-zero-fill NKR total woody (used in all tests) ===\n")
print(summarise_richness(NKR_richness$woody, "TotalWoodySpecies"))
cat("\n=== Post-zero-fill NHL total woody ===\n")
print(summarise_richness(NHL_richness$woody, "TotalWoodySpecies"))


# =============================================================================
# SECTION 2: Wilcoxon rank-sum pairwise tests
# =============================================================================
# Pairwise comparisons: each plantation stand type vs adjacent native,
# and among plantation types. No multiple-comparison correction applied
# (see Methods for rationale).

run_wilcoxon_pairs <- function(plot_counts, value_col, pairs,
                               site_label, metric_label) {
  results <- lapply(pairs, function(pr) {
    sub <- plot_counts %>% filter(Group %in% pr)
    if (n_distinct(sub$Group) < 2) {
      return(data.frame(Site = site_label, Metric = metric_label,
                        Comparison = paste(pr, collapse = " vs "),
                        W = NA, n1 = NA, n2 = NA, P = NA, Note = "Missing group"))
    }
    wt          <- wilcox.test(sub[[value_col]] ~ droplevels(sub$Group),
                               exact = FALSE)
    n_per_group <- table(droplevels(sub$Group))
    data.frame(Site       = site_label,
               Metric     = metric_label,
               Comparison = paste(pr, collapse = " vs "),
               W          = unname(wt$statistic),
               n1         = n_per_group[1],
               n2         = n_per_group[2],
               P          = wt$p.value,
               Note       = "")
  })
  dplyr::bind_rows(results)
}

pair_list <- list(
  c("Mixed",     "Native"),
  c("Unstocked", "Native"),
  c("Stocked",   "Native"),
  c("Mixed",     "Stocked"),
  c("Stocked",   "Unstocked")
)

NKR_woody_wilcox        <- run_wilcoxon_pairs(NKR_richness$woody,
                                              "TotalWoodySpecies",  pair_list, "NKR", "Total woody")
NKR_seedling_wilcox     <- run_wilcoxon_pairs(NKR_richness$seedlings,
                                              "SeedlingRichness",   pair_list, "NKR", "Seedling richness")
NKR_seedling_abund_wilcox <- run_wilcoxon_pairs(NKR_richness$seedling_abundance,
                                                "SeedlingAbundance",  pair_list, "NKR", "Seedling abundance")
NHL_woody_wilcox        <- run_wilcoxon_pairs(NHL_richness$woody,
                                              "TotalWoodySpecies",  pair_list, "NHL", "Total woody")
NHL_seedling_wilcox     <- run_wilcoxon_pairs(NHL_richness$seedlings,
                                              "SeedlingRichness",   pair_list, "NHL", "Seedling richness")
NHL_seedling_abund_wilcox <- run_wilcoxon_pairs(NHL_richness$seedling_abundance,
                                                "SeedlingAbundance",  pair_list, "NHL", "Seedling abundance")

wilcoxon_results_all <- dplyr::bind_rows(
  NKR_woody_wilcox,    NKR_seedling_wilcox,    NKR_seedling_abund_wilcox,
  NHL_woody_wilcox,    NHL_seedling_wilcox,    NHL_seedling_abund_wilcox
)

cat("\n=== Wilcoxon rank-sum results (all comparisons) ===\n")
print(wilcoxon_results_all, row.names = FALSE)
write.csv(wilcoxon_results_all, "wilcoxon_results_all.csv", row.names = FALSE)


# =============================================================================
# SECTION 3: Figure 2 — boxplots with Wilcoxon significance brackets
# =============================================================================

forest_cols <- c("Native"    = "#005F73",
                 "Stocked"   = "#E29578",
                 "Unstocked" = "#FFB703",
                 "Mixed"     = "#007A33")

box_comparisons <- list(c("Mixed",     "Native"),
                        c("Stocked",   "Unstocked"),
                        c("Mixed",     "Stocked"),
                        c("Native",    "Unstocked"))

plot_richness_boxplot <- function(plot_counts, value_col, comparisons,
                                  title_label, y_label) {
  ggplot(plot_counts, aes(x = Group, y = .data[[value_col]], fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 20,
                 size = 3, colour = "black") +
    scale_fill_manual(values = forest_cols) +
    stat_compare_means(comparisons = comparisons,
                       method      = "wilcox.test",
                       label       = "p.signif",
                       hide.ns     = TRUE) +
    labs(title = title_label, x = "Forest type", y = y_label) +
    theme_classic(base_size = 14) +
    theme(axis.title        = element_text(face = "bold"),
          axis.text         = element_text(colour = "black"),
          legend.position   = "none",
          plot.title        = element_text(face = "bold", hjust = 0.5))
}

p_nkr_woody <- plot_richness_boxplot(
  NKR_richness$woody, "TotalWoodySpecies", box_comparisons,
  "Total woody species richness — NKR", "Total woody species richness")

p_nhl_woody <- plot_richness_boxplot(
  NHL_richness$woody, "TotalWoodySpecies", box_comparisons,
  "Total woody species richness — NHL", "Total woody species richness")

p_nkr_seed <- plot_richness_boxplot(
  NKR_richness$seedlings, "SeedlingRichness", box_comparisons,
  "Seedling richness — NKR", "Seedling richness")

p_nhl_seed <- plot_richness_boxplot(
  NHL_richness$seedlings, "SeedlingRichness", box_comparisons,
  "Seedling richness — NHL", "Seedling richness")

figure2 <- plot_grid(p_nkr_woody, p_nhl_woody,
                     p_nkr_seed,  p_nhl_seed,
                     ncol = 2, labels = c("a)", "b)", "c)", "d)"),
                     label_size = 12)

ggsave("Figure2_boxplots.png",  figure2, width = 12, height = 10, dpi = 300)
ggsave("Figure2_boxplots.tiff", figure2, width = 12, height = 10, dpi = 300,
       compression = "lzw")
cat("\nFigure 2 saved.\n")


# =============================================================================
# SECTION 4: PERMANOVA
# =============================================================================

build_presence_matrix <- function(data_all) {
  data_all %>%
    filter(!is.na(Species), Species != "") %>%
    select(PlotID, Group, Species) %>%
    distinct() %>%
    mutate(Presence = 1) %>%
    pivot_wider(names_from  = Species,
                values_from = Presence,
                values_fill = 0)
}

NKR_presence <- build_presence_matrix(NKR_all)
NHL_presence <- build_presence_matrix(NHL_all)

run_permanova <- function(presence_df, site_label) {
  species_mat <- presence_df %>% select(-PlotID, -Group)
  env         <- presence_df %>% select(PlotID, Group)
  set.seed(123)
  perm <- adonis2(species_mat ~ Group, data = env,
                  method = "bray", permutations = 999)
  print(perm)
  df_model <- perm$Df[1];   df_resid <- perm$Df[2]
  Fval     <- perm$F[1];    R2val    <- perm$R2[1]
  Pval     <- perm$`Pr(>F)`[1]
  cat(sprintf("\n%s PERMANOVA: F(%d,%d) = %.2f, R2 = %.3f, P = %.3f\n\n",
              site_label, df_model, df_resid, Fval, R2val, Pval))
  list(model    = perm, df_model = df_model, df_resid = df_resid,
       F = Fval, R2 = R2val, P = Pval)
}

NKR_permanova <- run_permanova(NKR_presence, "NKR")
NHL_permanova <- run_permanova(NHL_presence, "NHL")


# =============================================================================
# SECTION 5: Betadisper (test for homogeneity of multivariate dispersions)
# =============================================================================

run_betadisper <- function(presence_df, site_label) {
  species_mat <- presence_df %>% select(-PlotID, -Group)
  bd          <- betadisper(vegdist(species_mat, method = "bray"),
                            presence_df$Group)
  bd_anova <- anova(bd)
  print(bd_anova)
  set.seed(123)
  bd_perm <- permutest(bd, permutations = 999)
  print(bd_perm)
  df_model <- bd_anova$Df[1];        df_resid <- bd_anova$Df[2]
  Fval     <- bd_anova$`F value`[1]; Pval     <- bd_anova$`Pr(>F)`[1]
  cat(sprintf("\n%s betadisper: F(%d,%d) = %.2f, P = %.3f\n\n",
              site_label, df_model, df_resid, Fval, Pval))
  list(anova    = bd_anova, permutest = bd_perm,
       df_model = df_model, df_resid  = df_resid,
       F = Fval, P = Pval)
}

NKR_betadisper <- run_betadisper(NKR_presence, "NKR")
NHL_betadisper <- run_betadisper(NHL_presence, "NHL")


# =============================================================================
# SECTION 6: Canonical Correspondence Analysis (CCA)
# =============================================================================

run_cca <- function(presence_df, site_label) {
  species_mat <- presence_df %>% select(-PlotID, -Group)
  env         <- presence_df %>% select(PlotID, Group)
  cca_model   <- cca(species_mat ~ Group, data = env)
  set.seed(123)
  cca_test  <- anova.cca(cca_model, permutations = 999)
  print(cca_test)
  df_model  <- cca_test$Df[1];  df_resid <- cca_test$Df[2]
  Fval      <- cca_test$F[1];   Pval     <- cca_test$`Pr(>F)`[1]
  var_expl  <- cca_model$CCA$tot.chi / cca_model$tot.chi
  cat(sprintf("\n%s CCA: F(%d,%d) = %.3f, P = %.3f, variance explained = %.1f%%\n\n",
              site_label, df_model, df_resid, Fval, Pval, var_expl * 100))
  list(model       = cca_model, test     = cca_test,
       df_model    = df_model,  df_resid = df_resid,
       F = Fval, P = Pval, var_explained = var_expl)
}

NKR_cca <- run_cca(NKR_presence, "NKR")
NHL_cca <- run_cca(NHL_presence, "NHL")


# =============================================================================
# SECTION 7: Consolidated statistics summary table (Table 2 source)
# =============================================================================

stats_summary <- data.frame(
  Test     = c("PERMANOVA", "PERMANOVA", "betadisper", "betadisper", "CCA", "CCA"),
  Site     = c("NKR", "NHL", "NKR", "NHL", "NKR", "NHL"),
  Df_model = c(NKR_permanova$df_model, NHL_permanova$df_model,
               NKR_betadisper$df_model, NHL_betadisper$df_model,
               NKR_cca$df_model, NHL_cca$df_model),
  Df_resid = c(NKR_permanova$df_resid, NHL_permanova$df_resid,
               NKR_betadisper$df_resid, NHL_betadisper$df_resid,
               NKR_cca$df_resid, NHL_cca$df_resid),
  F        = round(c(NKR_permanova$F, NHL_permanova$F,
                     NKR_betadisper$F, NHL_betadisper$F,
                     NKR_cca$F, NHL_cca$F), 3),
  P        = round(c(NKR_permanova$P, NHL_permanova$P,
                     NKR_betadisper$P, NHL_betadisper$P,
                     NKR_cca$P, NHL_cca$P), 4)
)

cat("\n=== Consolidated PERMANOVA / betadisper / CCA summary ===\n")
print(stats_summary, row.names = FALSE)
write.csv(stats_summary, "df_stats_summary.csv", row.names = FALSE)


# =============================================================================
# SECTION 8: Appendix S3 / Figure 3a — Species frequency bar charts
# =============================================================================
# Frequency (%) = proportion of plots (out of 15) in which a species
# was recorded at least once — standard ecological frequency definition.

compute_frequency <- function(data_all, name_col, site_label, life_stage) {
  plots_per_group <- data_all %>%
    distinct(PlotID, Group) %>%
    count(Group, name = "TotalPlots")
  
  data_all %>%
    filter(!is.na(.data[[name_col]]), .data[[name_col]] != "") %>%
    distinct(PlotID, Group, .data[[name_col]]) %>%
    count(Group, .data[[name_col]], name = "PlotsPresent") %>%
    left_join(plots_per_group, by = "Group") %>%
    mutate(FreqPct    = 100 * PlotsPresent / TotalPlots,
           Site       = site_label,
           LifeStage  = life_stage) %>%
    rename(SpeciesName = !!name_col)
}

NKR_seedling_freq <- compute_frequency(NKR_all, "Seedlings", "NKR", "Seedling")
NHL_seedling_freq <- compute_frequency(NHL_all, "Seedlings", "NHL", "Seedling")
NKR_tree_freq     <- compute_frequency(NKR_all, "Species",   "NKR", "Tree")
NHL_tree_freq     <- compute_frequency(NHL_all, "Species",   "NHL", "Tree")

seedling_freq_all <- bind_rows(NKR_seedling_freq, NHL_seedling_freq)
tree_freq_all     <- bind_rows(NKR_tree_freq, NHL_tree_freq)

# Order species by descending mean frequency within each site
order_species <- function(freq_df) {
  freq_df %>%
    group_by(Site, SpeciesName) %>%
    mutate(MeanFreq = mean(FreqPct)) %>%
    ungroup() %>%
    arrange(Site, desc(MeanFreq))
}
seedling_freq_all <- order_species(seedling_freq_all)
tree_freq_all     <- order_species(tree_freq_all)

# Top-N trim (no-op here: NKR has 34, NHL has 40 seedling spp — both <= 40)
top_n_per_site <- function(freq_df, n) {
  freq_df %>%
    distinct(Site, SpeciesName, MeanFreq) %>%
    arrange(Site, desc(MeanFreq), SpeciesName) %>%
    group_by(Site) %>%
    slice_head(n = n) %>%
    ungroup() %>%
    select(Site, SpeciesName) %>%
    inner_join(freq_df, by = c("Site", "SpeciesName"))
}
seedling_freq_top40 <- top_n_per_site(seedling_freq_all, n = 40)

# Greyscale palette consistent across Figures 3-5
group_cols_grey <- c("Mixed"     = "grey75",
                     "Stocked"   = "grey15",
                     "Native"    = "grey95",
                     "Unstocked" = "grey45")

build_freq_panel <- function(site_label, panel_letter) {
  seed_site  <- seedling_freq_top40 %>% filter(Site == site_label)
  tree_site  <- tree_freq_all       %>% filter(Site == site_label)
  
  seed_order <- seed_site %>% distinct(SpeciesName, MeanFreq) %>%
    arrange(desc(MeanFreq)) %>% pull(SpeciesName)
  tree_order <- tree_site %>% distinct(SpeciesName, MeanFreq) %>%
    arrange(desc(MeanFreq)) %>% pull(SpeciesName)
  
  seed_site$SpeciesName <- factor(seed_site$SpeciesName, levels = seed_order)
  tree_site$SpeciesName <- factor(tree_site$SpeciesName, levels = tree_order)
  
  main_plot <- ggplot(seed_site,
                      aes(x = SpeciesName, y = FreqPct, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = group_cols_grey, name = "Forest type",
                      breaks = c("Mixed", "Native", "Stocked", "Unstocked")) +
    labs(x = "Species", y = "Frequency (%)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x           = element_text(angle = 60, hjust = 1,
                                               size = 7, colour = "black"),
          axis.text.y           = element_text(size = 7, colour = "black"),
          axis.title            = element_text(face = "bold", size = 9),
          panel.grid            = element_blank(),
          panel.border          = element_blank(),
          axis.line             = element_line(colour = "black", linewidth = 0.3),
          legend.position.inside = c(0.90, 0.88),
          legend.position       = "inside",
          legend.background     = element_rect(fill = "white", colour = NA),
          legend.title          = element_text(face = "bold", size = 8),
          legend.text           = element_text(size = 7),
          legend.key.size       = unit(0.35, "cm"),
          legend.key            = element_rect(fill = "white")) +
    guides(fill = guide_legend(nrow = 2, ncol = 2)) +
    annotate("text", x = -Inf, y = Inf, label = paste0(panel_letter, ")"),
             hjust = -0.3, vjust = 1.3, size = 5, fontface = "bold")
  
  inset_plot <- ggplot(tree_site,
                       aes(x = SpeciesName, y = FreqPct, fill = Group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = group_cols_grey, guide = "none") +
    scale_y_continuous(breaks = seq(0, 100, 25)) +
    labs(x = NULL, y = "Frequency (%)") +
    theme_bw(base_size = 8) +
    theme(axis.text.x    = element_text(angle = 60, hjust = 1,
                                        size = 5, colour = "black"),
          axis.text.y    = element_text(size = 6, colour = "black"),
          axis.title.y   = element_text(size = 7),
          panel.grid     = element_blank(),
          panel.border   = element_blank(),
          axis.line      = element_line(colour = "black", linewidth = 0.25),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin    = margin(3, 3, 3, 3))
  
  ggdraw(main_plot) +
    draw_plot(inset_plot, x = 0.28, y = 0.52, width = 0.45, height = 0.46)
}

panel_a <- build_freq_panel("NKR", "a")
panel_b <- build_freq_panel("NHL", "b")
figure3a <- plot_grid(panel_a, panel_b, ncol = 1, rel_heights = c(1, 1))

ggsave("Figure3a_species_frequency.png",  figure3a,
       width = 16, height = 14, dpi = 300)
ggsave("Figure3a_species_frequency.tiff", figure3a,
       width = 16, height = 14, dpi = 300, compression = "lzw")

write.csv(seedling_freq_all, "seedling_frequency_by_group.csv", row.names = FALSE)
write.csv(tree_freq_all,     "tree_frequency_by_group.csv",     row.names = FALSE)
cat("\nFigure 3a (species frequency) and CSVs saved.\n")


# =============================================================================
# SECTION 9: Figure 3b — DCA ordination (three panels)
# =============================================================================
# Panel (a): combined inter-site ordination (NKR vs NHL)
# Panel (b): within-site NKR, coloured by forest type
# Panel (c): within-site NHL, coloured by forest type

# Build combined presence matrix for inter-site ordination
NKR_presence_tagged <- NKR_presence %>% mutate(Site = "NKR")
NHL_presence_tagged <- NHL_presence %>% mutate(Site = "NHL")

combined_presence <- bind_rows(NKR_presence_tagged, NHL_presence_tagged)
species_cols <- setdiff(colnames(combined_presence), c("PlotID", "Group", "Site"))
combined_presence[species_cols] <- lapply(combined_presence[species_cols],
                                          function(x) ifelse(is.na(x), 0, x))

cat(sprintf("\nCombined matrix: %d plots, %d species.\n",
            nrow(combined_presence), length(species_cols)))

run_dca <- function(presence_df, exclude_cols = c("PlotID", "Group", "Site")) {
  species_mat <- presence_df %>% select(-any_of(exclude_cols))
  species_mat <- species_mat[, colSums(species_mat) > 0, drop = FALSE]
  decorana(species_mat)
}

dca_combined <- run_dca(combined_presence)
dca_nkr      <- run_dca(NKR_presence)
dca_nhl      <- run_dca(NHL_presence)

print(dca_combined); print(dca_nkr); print(dca_nhl)

extract_dca_scores <- function(dca_model, presence_df, group_col) {
  scores_df           <- as.data.frame(scores(dca_model, display = "sites",
                                              choices = 1:2))
  scores_df[[group_col]] <- presence_df[[group_col]]
  scores_df
}

scores_combined <- extract_dca_scores(dca_combined, combined_presence, "Site")
scores_nkr      <- extract_dca_scores(dca_nkr, NKR_presence, "Group")
scores_nhl      <- extract_dca_scores(dca_nhl, NHL_presence, "Group")

# Colour and shape schemes (greyscale for print)
site_cols   <- c("NKR" = "grey20", "NHL" = "grey70")
site_shapes <- c("NKR" = 16,       "NHL" = 17)

group_shapes <- c("Mixed"     = 21, "Stocked"   = 22,
                  "Native"    = 23, "Unstocked" = 24)

dca_theme <- theme_bw(base_size = 10) +
  theme(panel.grid     = element_blank(),
        panel.border   = element_blank(),
        axis.line      = element_line(colour = "black", linewidth = 0.3),
        axis.text      = element_text(colour = "black", size = 8),
        axis.title     = element_text(face = "bold", size = 9),
        legend.title   = element_text(face = "bold", size = 8),
        legend.text    = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))

panel_dca_a <- ggplot(scores_combined,
                      aes(x = DCA1, y = DCA2, colour = Site, shape = Site)) +
  geom_point(size = 2.5, stroke = 1) +
  scale_colour_manual(values = site_cols) +
  scale_shape_manual(values = site_shapes) +
  labs(x = "DCA1", y = "DCA2") + dca_theme +
  annotate("text", x = -Inf, y = Inf, label = "a)", hjust = -0.3, vjust = 1.3,
           size = 5, fontface = "bold")

panel_dca_b <- ggplot(scores_nkr,
                      aes(x = DCA1, y = DCA2, fill = Group, shape = Group)) +
  geom_point(size = 2.5, colour = "black", stroke = 0.4) +
  scale_fill_manual(values = group_cols_grey) +
  scale_shape_manual(values = group_shapes) +
  labs(x = "DCA1", y = "DCA2") + dca_theme +
  theme(legend.position = "none") +
  annotate("text", x = -Inf, y = Inf, label = "b)", hjust = -0.3, vjust = 1.3,
           size = 5, fontface = "bold")

panel_dca_c <- ggplot(scores_nhl,
                      aes(x = DCA1, y = DCA2, fill = Group, shape = Group)) +
  geom_point(size = 2.5, colour = "black", stroke = 0.4) +
  scale_fill_manual(values = group_cols_grey) +
  scale_shape_manual(values = group_shapes) +
  labs(x = "DCA1", y = "DCA2") + dca_theme +
  guides(fill = guide_legend(override.aes = list(shape = group_shapes))) +
  annotate("text", x = -Inf, y = Inf, label = "c)", hjust = -0.3, vjust = 1.3,
           size = 5, fontface = "bold")

shared_dca_legend <- get_legend(panel_dca_c)
panel_dca_c_no_leg <- panel_dca_c + theme(legend.position = "none")

bottom_row <- plot_grid(panel_dca_b, panel_dca_c_no_leg, shared_dca_legend,
                        ncol = 3, rel_widths = c(3, 3, 1))
figure3b <- plot_grid(panel_dca_a, bottom_row, ncol = 1, rel_heights = c(1, 1))

ggsave("Figure3b_DCA_ordination.png",  figure3b,
       width = 12, height = 10, dpi = 300)
ggsave("Figure3b_DCA_ordination.tiff", figure3b,
       width = 12, height = 10, dpi = 300, compression = "lzw")

write.csv(scores_combined, "dca_scores_combined_intersite.csv",  row.names = TRUE)
write.csv(scores_nkr,      "dca_scores_nkr_within_site.csv",     row.names = TRUE)
write.csv(scores_nhl,      "dca_scores_nhl_within_site.csv",     row.names = TRUE)
cat("\nFigure 3b (DCA) saved.\n")


# =============================================================================
# SECTION 10: Figure 4 — Rank-abundance curves (RACs) for seedlings
# =============================================================================
# log10(abundance + 1) transformation following Magurran (2004).

compute_rank_abundance <- function(data_all, site_label) {
  data_all %>%
    filter(!is.na(Seedlings), Seedlings != "") %>%
    mutate(Abundance = as.numeric(Abundance)) %>%
    group_by(Group, Seedlings) %>%
    summarise(TotalAbundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Group) %>%
    arrange(Group, desc(TotalAbundance)) %>%
    mutate(Rank           = row_number(),
           Log10Abundance = log10(TotalAbundance + 1)) %>%
    ungroup() %>%
    mutate(Site = site_label)
}

rac_nkr <- compute_rank_abundance(NKR_all, "NKR")
rac_nhl <- compute_rank_abundance(NHL_all, "NHL")
rac_all <- bind_rows(rac_nkr, rac_nhl)

cat("\n=== RAC species count per forest type ===\n")
print(rac_all %>%
        group_by(Site, Group) %>%
        summarise(MaxRank = max(Rank), .groups = "drop"))
write.csv(rac_all, "rank_abundance_seedlings_all.csv", row.names = FALSE)

rac_theme <- theme_bw(base_size = 10) +
  theme(panel.grid    = element_blank(),
        panel.border  = element_blank(),
        axis.line     = element_line(colour = "black", linewidth = 0.3),
        axis.text     = element_text(colour = "black", size = 8),
        axis.title    = element_text(face = "bold", size = 9),
        legend.title  = element_text(face = "bold", size = 8),
        legend.text   = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"))

build_rac_panel <- function(site_label, panel_letter) {
  rac_all %>%
    filter(Site == site_label) %>%
    ggplot(aes(x = Rank, y = Log10Abundance,
               colour = Group, shape = Group, group = Group)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2, fill = "white", stroke = 0.8) +
    scale_colour_manual(values = group_cols_grey) +
    scale_shape_manual(values  = group_shapes) +
    labs(x = "Species rank",
         y = expression(log[10]*"(density + 1)")) +
    rac_theme +
    annotate("text", x = -Inf, y = Inf, label = paste0(panel_letter, ")"),
             hjust = -0.3, vjust = 1.3, size = 5, fontface = "bold")
}

panel_rac_a <- build_rac_panel("NKR", "a") + theme(legend.position = "none")
panel_rac_b <- build_rac_panel("NHL", "b")
shared_rac_legend  <- get_legend(panel_rac_b)
panel_rac_b_no_leg <- panel_rac_b + theme(legend.position = "none")

figure4 <- plot_grid(panel_rac_a, panel_rac_b_no_leg, shared_rac_legend,
                     ncol = 3, rel_widths = c(3, 3, 1))

ggsave("Figure4_rank_abundance.png",  figure4, width = 11, height = 5, dpi = 300)
ggsave("Figure4_rank_abundance.tiff", figure4, width = 11, height = 5, dpi = 300,
       compression = "lzw")
cat("\nFigure 4 (RACs) saved.\n")


# =============================================================================
# SECTION 11: Figure 5 — Individual-based species accumulation curves (SACs)
# =============================================================================
# specaccum(method = "rarefaction") accumulates over individuals (not plots),
# so it requires a plot-by-species matrix of individual counts.
#
# Trees:    one row per stem in the raw data → abundance = row count per plot
# Seedlings: Abundance column = individual count; use sum(Abundance) per plot

build_tree_abund_matrix <- function(data_all, group_filter) {
  data_all %>%
    filter(!is.na(Species), Species != "", Group == group_filter) %>%
    count(PlotID, Species, name = "Count") %>%
    pivot_wider(names_from  = Species,
                values_from = Count,
                values_fill = 0) %>%
    select(-PlotID)
}

build_seedling_abund_matrix <- function(data_all, group_filter) {
  data_all %>%
    filter(!is.na(Seedlings), Seedlings != "", Group == group_filter) %>%
    mutate(Abundance = as.numeric(Abundance)) %>%
    group_by(PlotID, Seedlings) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from  = Seedlings,
                values_from = Total,
                values_fill = 0) %>%
    select(-PlotID)
}

run_sac_all <- function(data_all, site_label) {
  results <- list()
  for (life_stage in c("Tree", "Seedling")) {
    for (ft in group_levels) {
      mat <- if (life_stage == "Tree") {
        build_tree_abund_matrix(data_all, ft)
      } else {
        build_seedling_abund_matrix(data_all, ft)
      }
      if (nrow(mat) < 2 || sum(mat) == 0) {
        cat(sprintf("Skipping %s %s %s: insufficient data.\n",
                    site_label, life_stage, ft))
        next
      }
      sac <- specaccum(mat, method = "rarefaction")
      results[[paste(site_label, life_stage, ft)]] <- data.frame(
        Site       = site_label,
        LifeStage  = life_stage,
        Group      = ft,
        Individuals = sac$individuals,
        Richness   = sac$richness,
        SD         = sac$sd
      )
    }
  }
  bind_rows(results)
}

sac_nkr <- run_sac_all(NKR_all, "NKR")
sac_nhl <- run_sac_all(NHL_all, "NHL")
sac_all <- bind_rows(sac_nkr, sac_nhl)
write.csv(sac_all, "sac_individual_based_all.csv", row.names = FALSE)

sac_theme <- theme_bw(base_size = 10) +
  theme(panel.grid      = element_blank(),
        panel.border    = element_blank(),
        axis.line       = element_line(colour = "black", linewidth = 0.3),
        axis.text       = element_text(colour = "black", size = 8),
        axis.title      = element_text(face = "bold", size = 9),
        strip.background = element_rect(fill = "grey95", colour = "black",
                                        linewidth = 0.3),
        strip.text      = element_text(face = "bold", size = 8),
        legend.title    = element_text(face = "bold", size = 8),
        legend.text     = element_text(size = 7),
        legend.key.size = unit(0.35, "cm"),
        legend.direction = "horizontal")

build_sac_panel <- function(site_label, panel_letter = NULL) {
  site_data <- sac_all %>%
    filter(Site == site_label) %>%
    mutate(LifeStage = factor(LifeStage, levels = c("Tree", "Seedling")))
  
  p <- ggplot(site_data,
              aes(x = Individuals, y = Richness,
                  colour = Group, fill = Group)) +
    geom_ribbon(aes(ymin = pmax(Richness - SD, 0),
                    ymax = Richness + SD),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ LifeStage, scales = "free_x") +
    scale_colour_manual(values = group_cols_grey) +
    scale_fill_manual(values   = group_cols_grey) +
    labs(x = "Number of individuals", y = NULL) +
    sac_theme +
    guides(colour = guide_legend(nrow = 1),
           fill   = guide_legend(nrow = 1))
  
  if (!is.null(panel_letter)) {
    p <- p + annotate("text", x = -Inf, y = Inf,
                      label = paste0(panel_letter, ")"),
                      hjust = -0.3, vjust = 1.3,
                      size = 5, fontface = "bold")
  }
  p
}

panel_sac_a <- build_sac_panel("NKR", "a") + theme(legend.position = "none")
panel_sac_b <- build_sac_panel("NHL", "b")
shared_sac_legend  <- get_legend(panel_sac_b)
panel_sac_b_no_leg <- panel_sac_b + theme(legend.position = "none")

y_label <- ggdraw() +
  draw_label("Cumulative number of species", angle = 90, size = 9,
             fontface = "bold", x = 0.5, y = 0.5, hjust = 0.5, vjust = 0.5)

panels_stack      <- plot_grid(panel_sac_a, panel_sac_b_no_leg,
                               ncol = 1, rel_heights = c(1, 1))
panels_with_label <- plot_grid(y_label, panels_stack,
                               ncol = 2, rel_widths = c(0.06, 1))
figure5           <- plot_grid(panels_with_label, shared_sac_legend,
                               ncol = 1, rel_heights = c(8, 1))

ggsave("Figure5_SAC.png",  figure5, width = 9, height = 11, dpi = 300)
ggsave("Figure5_SAC.tiff", figure5, width = 9, height = 11, dpi = 300,
       compression = "lzw")
cat("\nFigure 5 (SACs) saved.\n")


# =============================================================================
# SECTION 12: Table 1 — Richness and abundance with significance symbols
# =============================================================================

sig_symbol <- function(p) {
  if (is.na(p))    return("")
  if (p < 0.001)   return("***")
  if (p < 0.01)    return("**")
  if (p < 0.05)    return("*")
  return("")
}

get_sig_vs_native <- function(site_label, metric_label) {
  wilcoxon_results_all %>%
    filter(Site == site_label, Metric == metric_label,
           grepl("vs Native$", Comparison)) %>%
    mutate(Group = sub(" vs Native$", "", Comparison),
           Sig   = sapply(P, sig_symbol)) %>%
    select(Group, Sig)
}

format_cell <- function(summary_df, sig_df) {
  summary_df %>%
    left_join(sig_df, by = "Group") %>%
    mutate(Sig  = ifelse(is.na(Sig), "", Sig),
           Cell = sprintf("%.2f \u00b1 %.2f%s", Mean, SE, Sig)) %>%
    select(Group, Cell)
}

build_row <- function(summary_df, sig_df, row_label) {
  cells <- format_cell(summary_df, sig_df)
  ordered <- cells %>%
    mutate(Group = factor(Group,
                          levels = c("Stocked", "Unstocked", "Mixed", "Native"))) %>%
    arrange(Group)
  setNames(c(row_label, ordered$Cell),
           c("Metric", "Stocked", "Unstocked", "Mixed", "Native"))
}

nkr_woody_sum <- summarise_richness(NKR_richness$woody,             "TotalWoodySpecies")
nhl_woody_sum <- summarise_richness(NHL_richness$woody,             "TotalWoodySpecies")
nkr_seedR_sum <- summarise_richness(NKR_richness$seedlings,         "SeedlingRichness")
nhl_seedR_sum <- summarise_richness(NHL_richness$seedlings,         "SeedlingRichness")
nkr_seedA_sum <- summarise_richness(NKR_richness$seedling_abundance,"SeedlingAbundance")
nhl_seedA_sum <- summarise_richness(NHL_richness$seedling_abundance,"SeedlingAbundance")

table1_rows <- bind_rows(
  build_row(nkr_woody_sum, get_sig_vs_native("NKR","Total woody"),        "Total woody (NKR)"),
  build_row(nhl_woody_sum, get_sig_vs_native("NHL","Total woody"),        "Total woody (NHL)"),
  build_row(nkr_seedR_sum, get_sig_vs_native("NKR","Seedling richness"), "Seedling richness (NKR)"),
  build_row(nhl_seedR_sum, get_sig_vs_native("NHL","Seedling richness"), "Seedling richness (NHL)"),
  build_row(nkr_seedA_sum, get_sig_vs_native("NKR","Seedling abundance"),"Seedling abundance (NKR)"),
  build_row(nhl_seedA_sum, get_sig_vs_native("NHL","Seedling abundance"),"Seedling abundance (NHL)")
)

cat("\n=== Table 1 ===\n")
print(table1_rows, row.names = FALSE)
write.csv(table1_rows, "Table1_richness_abundance.csv", row.names = FALSE)


# =============================================================================
# SECTION 13: Table 2 — PERMANOVA / betadisper / CCA summary
# =============================================================================

format_effect <- function(test_name, R2 = NA, var_expl_pct = NA) {
  if (test_name == "PERMANOVA") return(sprintf("R\u00b2 = %.3f", R2))
  if (test_name == "betadisper") return("\u2014")
  if (test_name == "CCA") return(sprintf("%.1f%% of total variance",
                                         var_expl_pct * 100))
}

table2_rows <- bind_rows(
  data.frame(Test = "PERMANOVA", Site = "NKR",
             Df = sprintf("%d, %d", NKR_permanova$df_model, NKR_permanova$df_resid),
             F  = sprintf("%.2f", NKR_permanova$F),
             P  = sprintf("%.3f", NKR_permanova$P),
             Effect = format_effect("PERMANOVA", R2 = NKR_permanova$R2)),
  data.frame(Test = "PERMANOVA", Site = "NHL",
             Df = sprintf("%d, %d", NHL_permanova$df_model, NHL_permanova$df_resid),
             F  = sprintf("%.2f", NHL_permanova$F),
             P  = sprintf("%.3f", NHL_permanova$P),
             Effect = format_effect("PERMANOVA", R2 = NHL_permanova$R2)),
  data.frame(Test = "betadisper", Site = "NKR",
             Df = sprintf("%d, %d", NKR_betadisper$df_model, NKR_betadisper$df_resid),
             F  = sprintf("%.2f", NKR_betadisper$F),
             P  = sprintf("%.3f", NKR_betadisper$P),
             Effect = format_effect("betadisper")),
  data.frame(Test = "betadisper", Site = "NHL",
             Df = sprintf("%d, %d", NHL_betadisper$df_model, NHL_betadisper$df_resid),
             F  = sprintf("%.2f", NHL_betadisper$F),
             P  = sprintf("%.3f", NHL_betadisper$P),
             Effect = format_effect("betadisper")),
  data.frame(Test = "CCA", Site = "NKR",
             Df = sprintf("%d, %d", NKR_cca$df_model, NKR_cca$df_resid),
             F  = sprintf("%.3f", NKR_cca$F),
             P  = sprintf("%.3f", NKR_cca$P),
             Effect = format_effect("CCA", var_expl_pct = NKR_cca$var_explained)),
  data.frame(Test = "CCA", Site = "NHL",
             Df = sprintf("%d, %d", NHL_cca$df_model, NHL_cca$df_resid),
             F  = sprintf("%.3f", NHL_cca$F),
             P  = sprintf("%.3f", NHL_cca$P),
             Effect = format_effect("CCA", var_expl_pct = NHL_cca$var_explained))
)
colnames(table2_rows) <- c("Test", "Site", "Df", "F", "P", "Effect size")

cat("\n=== Table 2 ===\n")
print(table2_rows, row.names = FALSE)
write.csv(table2_rows, "Table2_PERMANOVA_betadisper_CCA.csv", row.names = FALSE)


# =============================================================================
# SECTION 14: Appendix S2 — GLM results (Poisson + QuasiPoisson)
# =============================================================================
# GLMs model richness as a Poisson count with forest type as predictor and
# adjacent secondary native forest as the reference level. Rate ratios
# (= exp[coefficient]) give the multiplicative difference in expected richness
# relative to native forest for each stand type.
# QuasiPoisson is fitted alongside to correct SEs where residual deviance/df > 1
# (confirmed overdispersion in NHL total woody richness; phi = 1.75).

build_richness_for_glm <- function(data_all) {
  all_plots <- data_all %>% distinct(PlotID, Group)
  
  woody_long <- data_all %>%
    pivot_longer(cols = c(Species, Seedlings, Climbers),
                 names_to  = "Layer",
                 values_to = "WoodyTaxon") %>%
    filter(!is.na(WoodyTaxon), WoodyTaxon != "")
  
  woody <- woody_long %>%
    group_by(PlotID, Group) %>%
    summarise(TotalWoody = n_distinct(WoodyTaxon), .groups = "drop")
  
  seedR_raw <- data_all %>%
    filter(!is.na(Seedlings), Seedlings != "") %>%
    group_by(PlotID, Group) %>%
    summarise(SeedlingRichness = n_distinct(Seedlings), .groups = "drop")
  
  seedR <- all_plots %>%
    left_join(seedR_raw, by = c("PlotID", "Group")) %>%
    mutate(SeedlingRichness = ifelse(is.na(SeedlingRichness), 0L,
                                     as.integer(SeedlingRichness)))
  
  woody %>% full_join(seedR, by = c("PlotID", "Group"))
}

nkr_glm_data <- build_richness_for_glm(NKR_all)
nhl_glm_data <- build_richness_for_glm(NHL_all)
nkr_glm_data$Group <- relevel(factor(nkr_glm_data$Group), ref = "Native")
nhl_glm_data$Group <- relevel(factor(nhl_glm_data$Group), ref = "Native")

fit_glm <- function(data, response_col, family_fn, site_label, metric_label,
                    model_type) {
  formula_str <- paste(response_col, "~ Group")
  model <- glm(as.formula(formula_str), data = data, family = family_fn)
  cat(sprintf("\n=== %s GLM: %s ~ forest type at %s ===\n",
              model_type, metric_label, site_label))
  if (model_type == "QuasiPoisson") {
    s <- summary(model)
    cat(sprintf("Dispersion parameter (phi): %.3f\n", s$dispersion))
  }
  print(summary(model))
  model
}

nkr_woody_poi   <- fit_glm(nkr_glm_data, "TotalWoody",       poisson(link="log"),
                           "NKR", "Total woody",       "Poisson")
nkr_seedR_poi   <- fit_glm(nkr_glm_data, "SeedlingRichness", poisson(link="log"),
                           "NKR", "Seedling richness", "Poisson")
nhl_woody_poi   <- fit_glm(nhl_glm_data, "TotalWoody",       poisson(link="log"),
                           "NHL", "Total woody",       "Poisson")
nhl_seedR_poi   <- fit_glm(nhl_glm_data, "SeedlingRichness", poisson(link="log"),
                           "NHL", "Seedling richness", "Poisson")

nkr_woody_qp    <- fit_glm(nkr_glm_data, "TotalWoody",       quasipoisson(link="log"),
                           "NKR", "Total woody",       "QuasiPoisson")
nkr_seedR_qp    <- fit_glm(nkr_glm_data, "SeedlingRichness", quasipoisson(link="log"),
                           "NKR", "Seedling richness", "QuasiPoisson")
nhl_woody_qp    <- fit_glm(nhl_glm_data, "TotalWoody",       quasipoisson(link="log"),
                           "NHL", "Total woody",       "QuasiPoisson")
nhl_seedR_qp    <- fit_glm(nhl_glm_data, "SeedlingRichness", quasipoisson(link="log"),
                           "NHL", "Seedling richness", "QuasiPoisson")

# Goodness-of-fit check
cat("\n=== Goodness-of-fit: residual deviance / df (Poisson models) ===\n")
for (info in list(
  list(nkr_woody_poi, "NKR total woody"),
  list(nkr_seedR_poi, "NKR seedling richness"),
  list(nhl_woody_poi, "NHL total woody"),
  list(nhl_seedR_poi, "NHL seedling richness")
)) {
  m     <- info[[1]]; label <- info[[2]]
  ratio <- m$deviance / m$df.residual
  flag  <- if (ratio > 2) " ** OVERDISPERSED" else if (ratio < 0.5) " (underdispersed)" else ""
  cat(sprintf("  %-30s deviance/df = %.2f%s\n", label, ratio, flag))
}

# Extract clean coefficient table
extract_glm_table <- function(model, site_label, metric_label, model_type) {
  s       <- summary(model)
  coef_df <- as.data.frame(s$coefficients)
  coef_df$Term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  coef_df$Term <- gsub("^Group", "", coef_df$Term)
  coef_df$Term[coef_df$Term == "(Intercept)"] <- "Native (reference)"
  coef_df$RateRatio <- exp(coef_df$Estimate)
  p_col          <- grep("Pr\\(", names(coef_df), value = TRUE)
  coef_df$P      <- coef_df[[p_col]]
  coef_df$Sig    <- cut(coef_df$P,
                        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                        labels = c("***", "**", "*", ""))
  data.frame(
    Site      = site_label,
    Metric    = metric_label,
    Model     = model_type,
    Term      = coef_df$Term,
    Estimate  = round(coef_df$Estimate, 4),
    SE        = round(coef_df$`Std. Error`, 4),
    RateRatio = round(coef_df$RateRatio, 3),
    Statistic = round(coef_df[[grep("value", names(coef_df))[1]]], 3),
    P         = signif(coef_df$P, 3),
    Sig       = as.character(coef_df$Sig),
    stringsAsFactors = FALSE
  )
}

appendix_s2 <- bind_rows(
  extract_glm_table(nkr_woody_poi, "NKR", "Total woody richness", "Poisson"),
  extract_glm_table(nkr_seedR_poi, "NKR", "Seedling richness",    "Poisson"),
  extract_glm_table(nhl_woody_poi, "NHL", "Total woody richness", "Poisson"),
  extract_glm_table(nhl_seedR_poi, "NHL", "Seedling richness",    "Poisson"),
  extract_glm_table(nkr_woody_qp,  "NKR", "Total woody richness", "QuasiPoisson"),
  extract_glm_table(nkr_seedR_qp,  "NKR", "Seedling richness",    "QuasiPoisson"),
  extract_glm_table(nhl_woody_qp,  "NHL", "Total woody richness", "QuasiPoisson"),
  extract_glm_table(nhl_seedR_qp,  "NHL", "Seedling richness",    "QuasiPoisson")
)

cat("\n=== Appendix S2: GLM results ===\n")
print(appendix_s2, row.names = FALSE)
write.csv(appendix_s2, "AppendixS2_GLM_results.csv", row.names = FALSE)
cat("\nAppendix S2 CSV saved.\n")


# =============================================================================
# SECTION 15: Session information
# =============================================================================

cat("\n=== Session information ===\n")
print(sessionInfo())

cat("\n=== All outputs written. Files produced: ===\n")
outputs <- c(
  "Figure2_boxplots.png / .tiff",
  "Figure3a_species_frequency.png / .tiff",
  "Figure3b_DCA_ordination.png / .tiff",
  "Figure4_rank_abundance.png / .tiff",
  "Figure5_SAC.png / .tiff",
  "Table1_richness_abundance.csv",
  "Table2_PERMANOVA_betadisper_CCA.csv",
  "AppendixS2_GLM_results.csv",
  "wilcoxon_results_all.csv",
  "df_stats_summary.csv",
  "seedling_frequency_by_group.csv",
  "tree_frequency_by_group.csv",
  "rank_abundance_seedlings_all.csv",
  "sac_individual_based_all.csv",
  "dca_scores_combined_intersite.csv",
  "dca_scores_nkr_within_site.csv",
  "dca_scores_nhl_within_site.csv"
)
cat(paste0("  ", outputs, "\n"), sep = "")



#####################Seedling detection probability#####################
library(ggplot2)

# --- Plot and subsample geometry -------------------------------------------
plot_area    <- 10 * 10                      # m^2, full 10 m x 10 m plot
quadrat_area <- 0.5 * 0.5                    # m^2, one nested quadrat
n_quadrats   <- 4
sample_area  <- quadrat_area * n_quadrats    # 1.0 m^2 total seedling subsample

p_actual <- sample_area / plot_area          # 0.01 -> 1% of plot area sampled

# --- Detection probability under complete spatial randomness ---------------
# If a species' N individuals are independently, uniformly distributed
# across the 100 m^2 plot, the probability that >=1 individual falls inside
# the sampled area is P(detect) = 1 - (1 - p)^N.
# This is a best-case (upper-bound) estimate: patchy/clumped distributions
# will detect worse than this, not better.
detect_prob <- function(N, p) 1 - (1 - p)^N

N_seq <- 1:20
df <- data.frame(
  N = N_seq,
  p = p_actual,
  scenario = "Actual sampling fraction (1%)"
)
df$P_detect <- detect_prob(df$N, df$p)

ggplot(df, aes(N, P_detect)) +
  geom_line(linewidth = 1, colour = "steelblue") +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed", colour = "grey50") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Individuals of a species present in the 100 m² plot (N)",
       y = "Probability of detecting the species\nin the seedling subsample",
       title = "Detection probability under complete spatial randomness",
       subtitle = "Seedling subsample = 1.0 m² of a 100 m² plot (1% sampling fraction)") +
  theme_minimal()

# --- Abundance needed to reach 50% / 90% detection probability -------------
N_for_p <- function(target, p) log(1 - target) / log(1 - p)

N_for_p(0.5, p_actual)   # individuals needed for a coin-flip chance of detection
N_for_p(0.9, p_actual)   # individuals needed for near-certain detection



#############################################

##overdispersion_check <- function(mod) {
# =============================================================================
# SECTION 2: Fit the four Poisson GLMs (forest type as predictor,
# Native as reference level, log link)
# =============================================================================

mod_nkr_total <- glm(TotalWoodySpecies ~ Group, data = NKR_richness$woody,
                     family = poisson(link = "log"))
mod_nkr_seed  <- glm(SeedlingRichness ~ Group, data = NKR_richness$seedlings,
                     family = poisson(link = "log"))
mod_nhl_total <- glm(TotalWoodySpecies ~ Group, data = NHL_richness$woody,
                     family = poisson(link = "log"))
mod_nhl_seed  <- glm(SeedlingRichness ~ Group, data = NHL_richness$seedlings,
                     family = poisson(link = "log"))

# =============================================================================
# SECTION 3: Overdispersion diagnostics
# =============================================================================

overdispersion_check <- function(mod) {
  dev_df  <- deviance(mod) / df.residual(mod)
  pearson <- sum(residuals(mod, type = "pearson")^2) / df.residual(mod)
  data.frame(deviance_over_df = dev_df, pearson_dispersion = pearson)
}

results <- lapply(list(NKR_total = mod_nkr_total,
                       NKR_seed  = mod_nkr_seed,
                       NHL_total = mod_nhl_total,
                       NHL_seed  = mod_nhl_seed),
                  overdispersion_check)

do.call(rbind, results)

# Direct quasipoisson dispersion check for the one model your text flags
# as overdispersed (should land near phi = 1.75)
mod_nhl_total_qp <- glm(TotalWoodySpecies ~ Group, data = NHL_richness$woody,
                        family = quasipoisson(link = "log"))
summary(mod_nhl_total_qp)$dispersion
