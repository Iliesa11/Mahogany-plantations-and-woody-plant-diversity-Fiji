# =============================================================================
# Woody plant composition and richness in Swietenia macrophylla plantations
# and secondary native forest, Viti Levu, Fiji — REVISED ANALYSIS SCRIPT
#
# Authors : Iliesa Koroi, Ole R. Vetaas, Tamara Osborne-Naikatini
# Target  : Pacific Science (University of Hawai'i Press)
# Version : 2026-10-26. Supersedes "Corrected mahogany vegetation composition
#           and richness analysis.R" (1,248 lines; New Forests submission).
#
# -----------------------------------------------------------------------------
# CHANGE LOG — errors in the previous script and how they are corrected
# (R1/R2 = New Forests reviewer IDs; D-codes = data audit of 30 Sep 2026)
# -----------------------------------------------------------------------------
#  1. Taxon names counted as typed: case/whitespace variants ('Wadamu' vs
#     'wadamu', 'Bau ') and spelling variants ('wabosocu'/'wabosucu') became
#     separate "species". -> Names normalised and mapped to botanical names via
#     taxon_lookup.csv; the script stops if any name is unmapped.          [D1]
#  2. 'Dead' and 'standing dead' stems counted as species in richness,
#     frequency, PERMANOVA, CCA and DCA. -> Removed.                        [D1]
#  3. Herbs (Zingiber, Ageratum, Stachytarpheta, Saccharum edule, Musa, ...),
#     sedges and epiphytic ferns counted as woody, contrary to the Methods.
#     -> Removed via habit_class in taxon_lookup.csv.                       [D2]
#  4. Climbers included in total woody richness, but none were entered for
#     NHL native plots, biasing NHL native richness downwards. -> Climbers
#     excluded at both sites (author decision, 30 Sep 2026).               [D10]
#  5. Composition analyses (PERMANOVA, betadisper, CCA, DCA) used the TREE
#     layer only, and S. macrophylla was NOT removed, although the Methods
#     state that mahogany was excluded from all ordination matrices.
#     -> Trees + seedlings; every composition analysis run with and without
#     mahogany and labelled accordingly.                          [R2-7, W4]
#  6. Sites analysed separately; similar results presented as "cross-site
#     validation". -> Joint models, site x forest type, interaction test.
#                                                             [R2-5, R2-m7]
#  7. Plots treated as independent (Wilcoxon, GLM, free permutations) despite
#     5 plots per transect. -> Transect random intercept in GLMMs; PERMANOVA
#     tests forest type against the among-transect mean square, permuting
#     whole transects within site.                        [R2-3, R1 L226]
#  8. No multiple-comparison correction; only 5 of 6 pairwise contrasts run
#     (mixed vs unstocked missing). -> All six contrasts per site from the
#     GLMM, Tukey-adjusted (emmeans).                  [R2-8, R1 L273, L413]
#  9. Figure 2 brackets came from ggpubr::stat_compare_means on a DIFFERENT
#     set of 4 comparisons (stocked vs native omitted), so figure and table
#     disagreed. -> Compact letters from the same adjusted contrasts.
#                                                          [R1 L386, Fig. 4]
# 10. Figure 2 used colour; all other figures greyscale. -> One greyscale
#     scheme; native forest drawn with a dark outline / dashed line so it is
#     visible (R1 Figs 3 and 5).
# 11. SAC ribbons were +/- 1 SD, but the caption says 95 % confidence
#     envelopes. -> +/- 1.96 SE; sample-based (plot) curves added.
#                                                         [R1 L447, R2-m5]
# 12. RAC y-axis labelled "density" but plotted summed abundance.
#     -> Seedling density per m2 computed explicitly.
# 13. Frequency figure counted 'Dead' and used local names. -> Botanical
#     names; dead stems excluded; denominator = plots actually sampled.
# 14. Detection-probability plot used N = 1-20, but the 50 % and 90 %
#     thresholds lie near 69 and 229 individuals. -> Range extended.
# 15. Duplicate GLM/dispersion code at the end of the script with a hard-coded
#     expected value (phi = 1.75). -> Removed; DHARMa diagnostics instead.
# 16. NHL native transect 2 crossed planted Pinus caribaea. -> Excluded from
#     the primary analysis; retained in sensitivity analysis S3.            [D9]
# 17. setwd() to a hard-coded Windows path, rm(list = ls()), no package
#     check/install block. -> Relative paths via data_dir; install block.
# 18. CCA and DCA duplicated the PERMANOVA information; CCA inertia was
#     over-interpreted. -> One NMDS for display.               [R1 L248, R2-6]
# 19. DATA: NHL_mixed.csv contained, besides the mixed-stand records, copies
#     of 100 % of NHL native, 95 % of NHL stocked and 84 % of NHL unstocked
#     rows (657 rows; 165 unique). This inflated NHL mixed richness (~18 vs
#     ~8 taxa per plot) and created the site x forest type interaction for
#     woody richness. -> Replaced with the 165 unique rows (30 Sep 2026);
#     the script now stops if > 20 identical rows are shared between files.
#     Forest-type main effects from the additive model are reported for
#     responses without an interaction (T3b).
# 20. FIGURES (2 Oct 2026, after checking the first exports): Fig. 3 letters
#     now follow the reported contrasts (averaged over sites where there is no
#     interaction); Fig. 4 NMDS uses step-across dissimilarities for display
#     because plots sharing no taxa made the ordination degenerate; Fig. 5
#     axis title and panel order corrected, SNF line heavier; Fig. 2 whole
#     number axis and (a)/(b) labels; Fig. S1 corrected spellings.
# 21. (10 Oct 2026) Fig. 4 NMDS drawn as one panel per site (eight ellipses in
#     one panel were unreadable); data_dir found automatically, with a clear
#     error naming the folders searched.
# 22. (10 Oct 2026) Section 4b added for reviewer risks: tree-only richness
#     (T9), stem-density covariate (T10), sensitivity without the three
#     pine-only and one garden plot (S13), 3-D NMDS stress (S14), software
#     citations (S15). No earlier result changes.
# 23. (10 Oct 2026) Author decision: tree-only richness is primary. Fig. 3
#     panel (a) is now tree taxa per plot with letters from T9; pooled woody
#     richness moves to Online Resource Fig. S4 (letters from T3b).
# 24. (10 Oct 2026) S16: the first-plot-removed (S10) and all-120-plot (S3)
#     checks repeated for tree-only richness. No earlier result changes.
# 25. (10 Oct 2026) Author decision: Fig. 4 shows axes 1-2 of a 3-D NMDS
#     (trymax 500, up to 3 restarts from previous.best); S14 reports the
#     same 3-D solution. Fig. 3 panels start at zero.
# 26. (11 Oct 2026) S17: category-free check of tree richness in plantation
#     plots against mahogany share of basal area (FHCL criteria unavailable).
# 27. (11 Oct 2026) S17b: the same share effect adjusted for log(stems).
# 28. (13 Oct 2026) Target journal now Pacific Science: figures exported at
#     140 mm (two-column width); each figure plus caption must fit 197 mm.
#     Supplementary figures renumbered in citation order: S1 frequency,
#     S2 pooled woody richness, S3 seedling rank-abundance, S4 detection.
#     Main figures renumbered in citation order: 2 accumulation, 3 species
#     overlap, 4 richness, 5 NMDS (tables likewise: old Table 5 is Table 2).
# 29. (26 Oct 2026) At 140 mm the one-row legends of the accumulation and NMDS
#     figures were cut off ("Mixe"): legends now in two rows. The stress
#     line under the NMDS removed (it repeats the caption).
# 30. (26 Oct 2026) data_dir set to the author's Pacific Science folder; the
#     script checks that the folder and all nine input files exist.
# 31. (26 Oct 2026) Fig. S4 y-axis title was clipped at 90 mm height: now on
#     two lines. Fig. S1 labels Podocarpus (kuasi) at genus level, as in the
#     Methods.
#
# Not errors in the code, but they must be reported correctly in the Methods:
#   * Trees were recorded from DBH >= 10 cm in practice (1 of 1,911 stems
#     < 10 cm), not >= 1 cm as stated.                                     [D11]
#   * Basal area exceeds 80 m2/ha in most groups: check whether DBH is girth
#     or plot size differs before reporting BA. The script warns.
#
# STATUS: written without access to R in the drafting environment; the
# logic was cross-checked against an independent Python implementation, but
# expect minor fixes on first run. Report the values this script produces.
# -----------------------------------------------------------------------------
#
# Design (as recorded in the data files)
#   2 sites (NKR = Nukurua, NHL = Nausori Highlands)
#   x 4 forest types (native = secondary native forest, stocked, unstocked, mixed)
#   x 3 transects per forest type per site (24 transects)
#   x 5 plots (10 m x 10 m) per transect = 120 plots
#   Seedlings (< 50 cm): four 0.5 m x 0.5 m quadrats per plot (1 m2 per plot)
#
# Analytical decisions (see the reviewer response matrix)
#   * Response = trees (DBH >= 1 cm) + seedlings. Climbers are EXCLUDED at both
#     sites because climber records were not entered for NHL native plots.
#   * "Dead" / "standing dead" stems, herbs, epiphytic ferns and unresolved
#     non-woody records are removed via taxon_lookup.csv (column habit_class).
#   * Replicate structure: plots nested in transects; transect enters every
#     model as a random intercept, and permutation tests move whole transects.
#   * Sites analysed jointly (site x forest type), replacing per-site tests.
#   * Pairwise contrasts: all six per site, Tukey-adjusted (emmeans).
#   * Composition: PERMANOVA (split-plot, transects permuted within site) +
#     PERMDISP + one NMDS for display. CCA and DCA are dropped.
#   * Every composition analysis is run for the whole community and for the
#     associated community (S. macrophylla excluded).
#   * NHL secondary native forest transect 2 (planted Pinus caribaea) is
#     excluded from the primary analysis: 115 plots on 23 transects. S3 repeats
#     the GLMM contrasts with all 120 plots.
#   * Tree layer: stems were recorded from DBH >= 10 cm in practice (1 of 1,911
#     stems < 10 cm); the Methods must say so.
#
# Inputs (working directory = folder holding this script; see PATHS)
#   NKR_native.csv NKR_stocked.csv NKR_unstocked.csv NKR_mixed.csv
#   NHL_native.csv NHL_stocked.csv NHL_unstocked.csv NHL_mixed.csv
#   taxon_lookup.csv   (local name -> botanical name, habit class, origin)
#
# Outputs: ./outputs/ (tables as .csv, figures as .tiff 600 dpi and .eps)
# =============================================================================

# ---- 0. Packages: check, install if missing, load ----------------------------
required_pkgs <- c(
  "tidyverse",  # dplyr, tidyr, readr, stringr, purrr, ggplot2, forcats
  "glmmTMB",    # GLMMs (Poisson, negative binomial, generalised Poisson)
  "DHARMa",     # simulation-based residual diagnostics for GLMMs
  "emmeans",    # estimated marginal means and adjusted pairwise contrasts
  "vegan",      # adonis2, betadisper, metaMDS, specaccum, rarefy
    "cowplot",    # multi-panel figures
  "ggupset",    # UpSet-style intersection plot in ggplot2
  "multcompView" # compact letter displays from adjusted pairwise P values
)
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                      logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  # On Windows without Rtools, never compile from source: when CRAN's source
  # version is newer than its binary (as for glmmTMB), take the binary.
  if (.Platform$OS.type == "windows") {
    options(install.packages.compile.from.source = "never")
    pkg_type <- "binary"
  } else {
    pkg_type <- getOption("pkgType")
  }
  tryCatch(
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org", type = pkg_type),
    error = function(e) stop("Package installation failed: ", conditionMessage(e))
  )
}
still_missing <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                       logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  stop("These packages could not be installed: ", paste(still_missing, collapse = ", "),
       "\nFor glmmTMB try: install.packages('glmmTMB', repos = ",
       "c('https://glmmtmb.r-universe.dev', 'https://cloud.r-project.org'))")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))
# multcomp attaches MASS, whose select() masks dplyr::select(); restore dplyr's
select <- dplyr::select
filter <- dplyr::filter
# emmeans and glmmTMB both export functions used below; namespace explicitly
# where a clash is possible.

set.seed(20260930)  # reproducible permutations, NMDS starts and SACs

# ---- PATHS -------------------------------------------------------------------
# Folder holding the nine input CSV files (author's working directory). If the
# files move, change this one line. The script stops with a clear message if
# the folder or any input file is missing.
data_dir <- "C:/Users/Administrator/Desktop/PhD Chapter Drafts/Chapter 4_Vegetation paper/Pacific Science"
# Use forward slashes in this path: R reads a backslash as an escape character.
required_files <- c("NKR_native.csv", "NKR_stocked.csv", "NKR_unstocked.csv", "NKR_mixed.csv",
                    "NHL_native.csv", "NHL_stocked.csv", "NHL_unstocked.csv", "NHL_mixed.csv",
                    "taxon_lookup.csv")
if (!dir.exists(data_dir)) {
  stop("data_dir does not exist: ", data_dir,
       "\nCheck the folder path (use forward slashes) and edit the data_dir line.")
}
missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
if (length(missing_files) > 0) {
  stop("These input files are missing from ", data_dir, ": ",
       paste(missing_files, collapse = ", "),
       "\nNote the file names: NHL_native.csv (with underscore), not NHLnative.csv.")
}
message("Reading data from: ", normalizePath(data_dir))
out_dir  <- file.path(data_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sites        <- c("NKR", "NHL")
forest_types <- c("native", "stocked", "unstocked", "mixed")
ft_labels    <- c(native = "Secondary native forest", stocked = "Stocked",
                  unstocked = "Unstocked", mixed = "Mixed")
mahogany     <- "Swietenia macrophylla"

# Author decision (30 Sep 2026): NHL secondary native forest transect 2 crossed
# planted Pinus caribaea on customary land. It is not secondary native forest,
# so it is EXCLUDED from the primary analysis and retained only in a
# sensitivity analysis (S3). Mahogany in NKR SNF is naturalised (author field
# observation) and is retained; NKR SNF is described as secondary forest with
# naturalised S. macrophylla.
exclude_units <- tibble::tibble(site = "NHL", ft = "native", transect = 2L)
is_excluded <- function(d) {
  paste(d$site, d$ft, d$transect) %in%
    paste(exclude_units$site, exclude_units$ft, exclude_units$transect)
}
plot_area_ha <- 0.01            # 10 m x 10 m
seed_area_m2 <- 1.0             # 4 x 0.25 m2 quadrats per plot

# Greyscale scheme (JAFPac prints in greyscale). Native carries a dark outline
# and a dashed line so it remains visible on white (R1, Figs 3 and 5).
ft_fill     <- c(native = "grey95", stocked = "grey15",
                 unstocked = "grey45", mixed = "grey75")
ft_linetype <- c(native = "dashed", stocked = "solid",
                 unstocked = "dotdash", mixed = "longdash")

# =============================================================================
# 1. Read, validate and clean
# =============================================================================
read_group <- function(site, ft) {
  f <- file.path(data_dir, sprintf("%s_%s.csv", site, ft))
  if (!file.exists(f)) stop("Input file not found: ", f)
  readr::read_csv(f, col_types = readr::cols(.default = "c"),
                  show_col_types = FALSE) |>
    mutate(site = site, ft = ft)
}

raw <- tidyr::expand_grid(site = sites, ft = forest_types) |>
  purrr::pmap(read_group) |>
  bind_rows()

required_cols <- c("Transect", "Plot", "Species", "DBH", "Seedlings", "Abundance")
if (!all(required_cols %in% names(raw))) {
  stop("Missing columns: ", paste(setdiff(required_cols, names(raw)), collapse = ", "))
}

# Guard against records copied between forest-type files (change log 19).
# Identical tree-only rows occur by chance between stands (same transect and
# plot number, a common species such as mahogany, and a common DBH such as
# 12 cm): in the verified NKR files 5-20 rows per file (<= 10 %) match this way.
# Copied files match wholesale (the original NHL_mixed.csv held 84-100 % of the
# other NHL files' rows). The check therefore uses the SHARE of each file's
# distinct rows found in another file at the same site, and stops above 25 %.
file_rows <- raw |>
  dplyr::select(site, ft, all_of(required_cols)) |>
  distinct()
cross_file_dups <- file_rows |>
  group_by(dplyr::pick(site, all_of(required_cols))) |>
  mutate(shared = n_distinct(ft) > 1) |>
  ungroup() |>
  group_by(site, ft) |>
  summarise(rows = n(),
            rows_also_in_another_stand = sum(shared),
            share_shared = mean(shared),
            shared_with_seedling_entry = sum(shared & !is.na(Seedlings)),
            .groups = "drop")
print(cross_file_dups)
if (any(cross_file_dups$share_shared > 0.25)) {
  stop("More than 25 % of a file's rows also occur in another forest-type file ",
       "at the same site; check for copied records (see cross_file_dups).")
}

lookup <- readr::read_csv(file.path(data_dir, "taxon_lookup.csv"),
                          show_col_types = FALSE)

# Normalise local names: trim, lower case (fixes 'Wadamu'/'wadamu', 'Bau ')
norm_name <- function(x) stringr::str_to_lower(stringr::str_squish(x))

long <- bind_rows(
  raw |> filter(!is.na(Species), Species != "") |>
    transmute(site, ft, transect = as.integer(Transect), plot = as.integer(Plot),
              layer = "tree", local_key = norm_name(Species),
              n = 1, dbh_cm = suppressWarnings(as.numeric(DBH))),
  raw |> filter(!is.na(Seedlings), Seedlings != "") |>
    transmute(site, ft, transect = as.integer(Transect), plot = as.integer(Plot),
              layer = "seedling", local_key = norm_name(Seedlings),
              n = coalesce(suppressWarnings(as.numeric(Abundance)), 1),
              dbh_cm = NA_real_)
)

# Stop if any local name is absent from the lookup (no silent losses)
unmapped <- setdiff(unique(long$local_key), lookup$local_key)
if (length(unmapped) > 0) {
  stop("Local names missing from taxon_lookup.csv: ", paste(unmapped, collapse = "; "))
}

long <- long |>
  left_join(lookup, by = "local_key") |>
  filter(habit_class %in% c("woody", "woody_habit_monocot_or_treefern")) |>
  mutate(ft = factor(ft, levels = forest_types),
         site = factor(site, levels = sites),
         transect_id = interaction(site, ft, transect, drop = TRUE, sep = "_"),
         plot_id = paste(transect_id, plot, sep = "_"),
         exotic = origin_sourced == "exotic")

# Log what was removed, for the Methods text
removed <- raw |>
  pivot_longer(c(Species, Seedlings), names_to = "layer_col", values_to = "raw_name") |>
  filter(!is.na(raw_name), raw_name != "") |>
  mutate(local_key = norm_name(raw_name)) |>
  left_join(lookup, by = "local_key") |>
  filter(!habit_class %in% c("woody", "woody_habit_monocot_or_treefern")) |>
  count(habit_class, local_key, name = "records")
readr::write_csv(removed, file.path(out_dir, "S0_removed_records.csv"))

# Complete plot frame (keeps plots with zero records)
plots_all <- tidyr::expand_grid(site = factor(sites, levels = sites),
                                ft = factor(forest_types, levels = forest_types),
                                transect = 1:3, plot = 1:5) |>
  mutate(transect_id = interaction(site, ft, transect, drop = TRUE, sep = "_"),
         plot_id = paste(transect_id, plot, sep = "_"))
stopifnot(nrow(plots_all) == 120)

# Check that every plot in the frame exists in the data
absent <- setdiff(plots_all$plot_id, unique(long$plot_id))
if (length(absent) > 0) warning("Plots with no woody records: ", paste(absent, collapse = ", "))

# Primary analysis set: excluded transect(s) removed (see exclude_units)
long_all <- long
plots    <- plots_all |> filter(!is_excluded(pick(site, ft, transect))) |>
  mutate(transect_id = droplevels(transect_id))
long     <- long_all  |> filter(!is_excluded(pick(site, ft, transect))) |>
  mutate(transect_id = droplevels(transect_id))
message("Primary analysis: ", nrow(plots), " plots on ",
        n_distinct(plots$transect_id), " transects (",
        nrow(plots_all) - nrow(plots), " plots excluded).")
n_plots_group <- count(plots, site, ft, name = "n_plots")

# =============================================================================
# 2. Plot-level metrics and stand attributes
# =============================================================================
compute_plot_metrics <- function(long, plots) long |>
  group_by(plot_id) |>
  summarise(
    woody_rich       = n_distinct(botanical),
    woody_rich_noMah = n_distinct(botanical[botanical != mahogany]),
    tree_rich        = n_distinct(botanical[layer == "tree"]),
    tree_rich_noMah  = n_distinct(botanical[layer == "tree" & botanical != mahogany]),
    # numbers of individuals, for the stem-density covariate (section 4b)
    n_stems          = sum(layer == "tree"),
    n_stems_noMah    = sum(layer == "tree" & botanical != mahogany),
    n_ind            = sum(n),                                   # trees + seedlings
    n_ind_noMah      = sum(n[botanical != mahogany]),
    seed_rich        = n_distinct(botanical[layer == "seedling"]),
    seed_ab          = sum(n[layer == "seedling"]),
    stems_ha         = sum(layer == "tree") / plot_area_ha,
    ba_m2ha          = sum(pi * (dbh_cm[layer == "tree"] / 200)^2, na.rm = TRUE) / plot_area_ha,
    mah_ba_share     = sum(pi * (dbh_cm[layer == "tree" & botanical == mahogany] / 200)^2, na.rm = TRUE) /
                       sum(pi * (dbh_cm[layer == "tree"] / 200)^2, na.rm = TRUE),
    exotic_stem_share = mean(exotic[layer == "tree"]),
    .groups = "drop"
  ) |>
  right_join(plots, by = "plot_id") |>
  mutate(across(c(woody_rich, woody_rich_noMah, tree_rich, tree_rich_noMah, seed_rich,
                  seed_ab, stems_ha, ba_m2ha, n_stems, n_stems_noMah, n_ind, n_ind_noMah),
                ~ coalesce(.x, 0)))

plot_metrics     <- compute_plot_metrics(long, plots)          # primary
plot_metrics_all <- compute_plot_metrics(long_all, plots_all)  # sensitivity (all 120 plots)

# Plausibility check on basal area. Values > 80 m2/ha are implausible for
# tropical plantations and secondary forest; they suggest DBH was recorded as
# girth (GBH) or that plot area differs from 0.01 ha. VERIFY before reporting.
ba_check <- plot_metrics |> group_by(site, ft) |>
  summarise(ba_mean = mean(ba_m2ha), .groups = "drop")
if (any(ba_check$ba_mean > 80)) {
  warning("Mean basal area exceeds 80 m2/ha in: ",
          paste(ba_check$site[ba_check$ba_mean > 80], ba_check$ft[ba_check$ba_mean > 80], collapse = "; "),
          ". Check DBH units (diameter vs girth) and plot dimensions.")
}

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
descriptives <- plot_metrics |>
  group_by(site, ft) |>
  summarise(across(c(woody_rich, woody_rich_noMah, tree_rich, seed_rich, seed_ab,
                     stems_ha, ba_m2ha, mah_ba_share, exotic_stem_share),
                   list(mean = ~ mean(.x, na.rm = TRUE), se = se)),
            .groups = "drop")
readr::write_csv(descriptives, file.path(out_dir, "T1_descriptives_by_site_foresttype.csv"))

# =============================================================================
# 3. GLMMs: richness and seedling abundance ~ site x forest type + (1 | transect)
# =============================================================================
# Fit one GLMM. The formula is built inside this function and the family is
# passed as an evaluated object, so every refit below is done explicitly with
# fit_glmm() rather than update(): update() re-evaluates the stored call in the
# global environment, where `family` is base R's family() function (this caused
# "no applicable method for 'family' applied to an object of class NULL").
fit_glmm <- function(response, family, dat = plot_metrics,
                     rhs = "site * ft + (1 | transect_id)") {
  f <- as.formula(paste(response, "~", rhs))
  tryCatch(glmmTMB::glmmTMB(f, data = dat, family = family),
           error = function(e) { warning(response, ": ", conditionMessage(e)); NULL })
}

# A fit is usable only if the Hessian is positive definite and all fixed
# effects are finite (genpois can fail silently with "NA/NaN function evaluation")
fit_ok <- function(m) {
  !is.null(m) && isTRUE(m$sdr$pdHess) && all(is.finite(glmmTMB::fixef(m)$cond))
}

# Richness counts are frequently under-dispersed; Poisson is then conservative.
# Generalised Poisson (genpois) accommodates both under- and over-dispersion.
check_dispersion <- function(m, label) {
  sim <- DHARMa::simulateResiduals(m, n = 1000, plot = FALSE)
  disp <- DHARMa::testDispersion(sim, plot = FALSE)
  zi   <- DHARMa::testZeroInflation(sim, plot = FALSE)
  tibble(model = label, dispersion_ratio = unname(disp$statistic),
         dispersion_p = disp$p.value, zero_infl_p = zi$p.value)
}

responses <- tribble(
  ~response,           ~family,
  "woody_rich",        "poisson",
  "woody_rich_noMah",  "poisson",
  "seed_rich",         "poisson",
  "seed_ab",           "nbinom2"
)

model_family <- list()   # records the family finally used for each response
models <- purrr::map2(responses$response, responses$family, function(r, fam) {
  fam_obj <- switch(fam, poisson = poisson(), nbinom2 = glmmTMB::nbinom2())
  m <- fit_glmm(r, fam_obj)
  if (!fit_ok(m)) { warning(r, ": primary fit not usable"); return(NULL) }
  d <- check_dispersion(m, r)
  # Refit with generalised Poisson if Poisson residuals are clearly mis-dispersed,
  # but keep it only if that fit is numerically stable
  if (fam == "poisson" && d$dispersion_p < 0.05) {
    message(r, ": Poisson dispersion p = ", signif(d$dispersion_p, 3),
            " (ratio ", round(d$dispersion_ratio, 2), ") -> trying generalised Poisson")
    m2 <- fit_glmm(r, glmmTMB::genpois())
    if (fit_ok(m2)) {
      m <- m2; fam_obj <- glmmTMB::genpois()
    } else {
      message(r, ": generalised Poisson fit unstable; keeping Poisson. ",
              "If dispersion ratio < 1 (under-dispersion), Poisson is conservative.")
    }
  }
  model_family[[r]] <<- fam_obj
  m
}) |> purrr::set_names(responses$response)

diagnostics <- purrr::imap_dfr(purrr::compact(models), check_dispersion)
readr::write_csv(diagnostics, file.path(out_dir, "S1_glmm_diagnostics.csv"))

# Interaction test: likelihood-ratio test of site x forest type
lrt_interaction <- purrr::imap_dfr(purrr::compact(models), function(m, r) {
  m_add <- fit_glmm(r, model_family[[r]], rhs = "site + ft + (1 | transect_id)")
  if (!fit_ok(m_add)) return(tibble(response = r, note = "additive model failed"))
  a <- anova(m_add, m)
  tibble(response = r,
         family = model_family[[r]]$family,
         chisq = a$Chisq[2], df = a$`Chi Df`[2], p = a$`Pr(>Chisq)`[2],
         transect_sd = sqrt(as.numeric(glmmTMB::VarCorr(m)$cond$transect_id)))
})
readr::write_csv(lrt_interaction, file.path(out_dir, "T2_glmm_interaction_LRT.csv"))

# All six pairwise contrasts within each site, Tukey-adjusted, on response scale
contrasts_tbl <- purrr::imap_dfr(purrr::compact(models), function(m, r) {
  emm <- emmeans::emmeans(m, ~ ft | site, type = "response")
  bind_rows(
    # emmeans names the back-transformed mean 'rate' (Poisson) or 'response'
    # (genpois, nbinom2); the latter clashed with the label column below and
    # was overwritten. Standardise it to 'estimate' first.
    as_tibble(summary(emm)) |>
      rename_with(~ "estimate", .cols = any_of(c("rate", "response", "emmean"))) |>
      mutate(response = r, what = "estimated mean"),
    as_tibble(summary(emmeans::contrast(emm, method = "pairwise", adjust = "tukey"),
                      infer = TRUE)) |> mutate(response = r, what = "contrast")
  )
})
readr::write_csv(contrasts_tbl, file.path(out_dir, "T3_glmm_emmeans_contrasts.csv"))

# Compact letter displays for Fig. 3 (replace asterisks; R1 Fig. 4 comment).
# Built directly from the Tukey-adjusted contrasts above with multcompView,
# because current emmeans no longer supplies a cld() method (the previous
# multcomp::cld() call failed with "no applicable method for 'cld'").
# Forest types sharing a letter do not differ at alpha = 0.05 (within site).
cld_tbl <- contrasts_tbl |>
  filter(what == "contrast") |>
  mutate(pair = stringr::str_replace(as.character(contrast), "\\s+[-/]\\s+", "-")) |>
  group_by(response, site) |>
  group_modify(function(d, key) {
    pv <- setNames(d$p.value, d$pair)
    lt <- multcompView::multcompLetters(pv, threshold = 0.05)$Letters
    tibble(ft = factor(names(lt), levels = forest_types), .group = unname(lt))
  }) |>
  ungroup()
readr::write_csv(cld_tbl, file.path(out_dir, "S2_compact_letters.csv"))

# Main effect of forest type from the ADDITIVE model (site + forest type).
# Report these, not the per-site contrasts, for any response whose site x
# forest type interaction is not significant in T2 (woody richness with and
# without mahogany, and seedling density, in the corrected data). Forest-type
# LRT: additive model vs site-only model (df = 3). Contrasts are averaged over
# sites and Tukey-adjusted.
additive_tbl <- purrr::imap_dfr(purrr::compact(models), function(m, r) {
  m_add  <- fit_glmm(r, model_family[[r]], rhs = "site + ft + (1 | transect_id)")
  m_site <- fit_glmm(r, model_family[[r]], rhs = "site + (1 | transect_id)")
  if (!fit_ok(m_add) || !fit_ok(m_site)) {
    return(tibble(response = r, note = "additive or site-only model failed"))
  }
  a   <- anova(m_site, m_add)
  emm <- emmeans::emmeans(m_add, ~ ft, type = "response")
  bind_rows(
    tibble(response = r, what = "forest-type LRT (additive vs site-only)",
           chisq = a$Chisq[2], df = a$`Chi Df`[2], p.value = a$`Pr(>Chisq)`[2]),
    as_tibble(summary(emm)) |>
      rename_with(~ "estimate", .cols = any_of(c("rate", "response", "emmean"))) |>
      mutate(response = r, what = "estimated mean (averaged over sites)"),
    as_tibble(summary(emmeans::contrast(emm, method = "pairwise", adjust = "tukey"),
                      infer = TRUE)) |>
      mutate(response = r, what = "contrast (averaged over sites)")
  )
})
readr::write_csv(additive_tbl, file.path(out_dir, "T3b_glmm_additive_foresttype.csv"))

# Robustness (S10): drop the first plot of each transect, which was placed for
# access, and refit the woody-richness models (interaction and additive forms)
plot_metrics_no1 <- dplyr::filter(plot_metrics, plot != 1)
s10 <- purrr::map_dfr(c("woody_rich", "woody_rich_noMah"), function(r) {
  purrr::map_dfr(c(interaction = "site * ft + (1 | transect_id)",
                   additive    = "site + ft + (1 | transect_id)"), function(rhs) {
    m <- fit_glmm(r, model_family[[r]], dat = plot_metrics_no1, rhs = rhs)
    if (!fit_ok(m)) return(tibble(response = r, model = rhs, note = "fit failed"))
    spec <- if (grepl("\\*", rhs)) ~ ft | site else ~ ft
    emmeans::emmeans(m, spec, type = "response") |>
      emmeans::contrast(method = "pairwise", adjust = "tukey") |>
      summary(infer = TRUE) |> as_tibble() |> mutate(response = r, model = rhs)
  })
})
readr::write_csv(s10, file.path(out_dir, "S10_robustness_without_plot1.csv"))

# Sensitivity: all 120 plots, i.e. INCLUDING the excluded NHL native transect 2
# (planted Pinus caribaea). Same model family as the primary fit.
sens_out <- purrr::imap_dfr(purrr::compact(models), function(m, r) {
  m_all <- fit_glmm(r, model_family[[r]], dat = plot_metrics_all)
  if (!fit_ok(m_all)) return(tibble(response = r, note = "sensitivity fit failed"))
  emmeans::emmeans(m_all, ~ ft | site, type = "response") |>
    emmeans::contrast(method = "pairwise", adjust = "tukey") |>
    summary(infer = TRUE) |> as_tibble() |> mutate(response = r)
})
readr::write_csv(sens_out, file.path(out_dir, "S3_sensitivity_including_NHL_native_T2.csv"))

# =============================================================================
# 4. Composition: split-plot PERMANOVA, PERMDISP, NMDS
# =============================================================================
make_pa <- function(dat, drop_mahogany = FALSE) {
  d <- if (drop_mahogany) filter(dat, botanical != mahogany) else dat
  m <- d |> distinct(plot_id, botanical) |> mutate(v = 1L) |>
    pivot_wider(names_from = botanical, values_from = v, values_fill = 0L)
  plots |> select(plot_id, site, ft, transect, transect_id) |>
    left_join(m, by = "plot_id") |>
    mutate(across(-c(plot_id, site, ft, transect, transect_id), ~ coalesce(.x, 0L)))
}

# Split-plot pseudo-F: forest type and site x forest type are tested against the
# among-transect (within site x forest type) mean square, and whole transects
# are permuted among forest types within each site (R2 #3, #5; R1 L327).
ss_within <- function(d2, g) {
  sum(vapply(split(seq_along(g), g, drop = TRUE),
             function(i) sum(d2[i, i]) / (2 * length(i)), numeric(1)))
}
splitplot_permanova <- function(D, meta, nperm = 4999) {
  d2 <- as.matrix(D)^2; N <- nrow(d2); sst <- sum(d2) / (2 * N)
  parts <- function(ft_vec) {
    cell <- paste(meta$site, ft_vec); tr <- paste(cell, meta$tr_key)
    ss_site <- sst - ss_within(d2, meta$site)
    ss_ft   <- sst - ss_within(d2, ft_vec)
    ss_cell <- sst - ss_within(d2, cell)
    c(site = ss_site, ft = ss_ft, int = ss_cell - ss_site - ss_ft,
      tr = ss_within(d2, cell) - ss_within(d2, tr), res = ss_within(d2, tr))
  }
  df_ft <- nlevels(factor(meta$ft)) - 1
  df_tr <- dplyr::n_distinct(meta$tr_key) - dplyr::n_distinct(paste(meta$site, meta$ft))
  obs <- parts(meta$ft)
  F_obs <- c(ft  = unname((obs["ft"]  / df_ft) / (obs["tr"] / df_tr)),
             int = unname((obs["int"] / df_ft) / (obs["tr"] / df_tr)))
  # transect-level label table, permuted within site
  tr_tab <- distinct(meta, site, tr_key, ft)
  F_perm <- replicate(nperm, {
    shuffled <- tr_tab |> group_by(site) |> mutate(ft_new = sample(ft)) |> ungroup()
    ft_new <- shuffled$ft_new[match(meta$tr_key, shuffled$tr_key)]
    p <- parts(ft_new)
    c((p["ft"] / df_ft) / (p["tr"] / df_tr), (p["int"] / df_ft) / (p["tr"] / df_tr))
  })
  tibble(term = c("site", "forest type", "site x forest type", "transect (within cell)", "residual"),
         SS = unname(obs), R2 = unname(obs) / sst,
         pseudo_F = c(NA, F_obs["ft"], F_obs["int"], NA, NA),
         P_perm = c(NA, (sum(F_perm[1, ] >= F_obs["ft"]) + 1) / (nperm + 1),
                    (sum(F_perm[2, ] >= F_obs["int"]) + 1) / (nperm + 1), NA, NA))
}

composition_results <- purrr::map_dfr(c(FALSE, TRUE), function(drop_mah) {
  pa <- make_pa(long, drop_mah)
  X  <- pa |> select(-c(plot_id, site, ft, transect, transect_id)) |> as.matrix()
  keep <- rowSums(X) > 0                                  # empty plots cannot be compared
  D  <- vegan::vegdist(X[keep, ], method = "bray", binary = TRUE)  # = Sorensen (R1 L250)
  meta <- pa[keep, ] |> transmute(site, ft, tr_key = as.character(transect_id))
  splitplot_permanova(D, meta) |>
    mutate(community = if (drop_mah) "associated (S. macrophylla excluded)" else "whole community",
           n_plots = sum(keep))
})
readr::write_csv(composition_results, file.path(out_dir, "T4_permanova_splitplot.csv"))

# Associated-community Sorensen matrix, used for PERMDISP and NMDS
pa_assoc <- make_pa(long, drop_mahogany = TRUE)
X_assoc  <- pa_assoc |> select(-c(plot_id, site, ft, transect, transect_id)) |> as.matrix()
keep_a   <- rowSums(X_assoc) > 0
D_assoc  <- vegan::vegdist(X_assoc[keep_a, ], method = "bray", binary = TRUE)
env_a    <- pa_assoc[keep_a, c("site", "ft", "transect_id")]
# Note: vegan::adonis2 cannot run this whole-transect restricted design here
# ("Design must be balanced if permuting 'strata'"), because NHL secondary
# native forest retains 2 transects rather than 3 after excluding the planted
# pine transect. The split-plot test in T4 handles the unbalanced design and is
# the reported test; adonis2 is therefore not run.

# PERMDISP (homogeneity of dispersion) by site x forest type
grp <- interaction(env_a$site, env_a$ft, sep = "_")
bd  <- vegan::betadisper(D_assoc, grp, type = "centroid")
bd_test <- vegan::permutest(bd, permutations = 4999)
# tapply() returns a 1-d array, which tibble rejects as a matrix column;
# summarise the distances with dplyr instead
permdisp_tbl <- tibble(group = as.character(grp), dist = as.numeric(bd$distances)) |>
  group_by(group) |>
  summarise(mean_dist = mean(dist), se_dist = se(dist), n = n(), .groups = "drop") |>
  mutate(F = bd_test$tab$F[1], df1 = bd_test$tab$Df[1], df2 = bd_test$tab$Df[2],
         P = bd_test$tab$`Pr(>F)`[1])
readr::write_csv(permdisp_tbl, file.path(out_dir, "T5_permdisp.csv"))

# Pairwise dispersion differences (S11), Tukey HSD on distances to centroid
s11 <- as.data.frame(TukeyHSD(bd)$group) |>
  tibble::rownames_to_column("pair") |> as_tibble()
readr::write_csv(s11, file.path(out_dir, "S11_permdisp_tukey.csv"))

# NMDS for display only (Sorensen, associated community, both sites)
# Several plots share no taxa with most other plots once mahogany is excluded
# (e.g. a plot holding only Ixora pelagica, or only Pinus caribaea), so their
# Sorensen dissimilarity to nearly every plot is 1. Such ties at the maximum
# make NMDS degenerate: one plot is thrown far from the rest and the others
# collapse (the first export had NMDS1 running to 2.3). For DISPLAY ONLY the
# saturated dissimilarities are replaced by extended (step-across)
# dissimilarities (De'ath 1999; vegan::stepacross). Tests (T4, T5) still use
# the original Sorensen matrix.
n_saturated <- sum(as.matrix(D_assoc)[upper.tri(as.matrix(D_assoc))] >= 1)
message("Plot pairs with no shared taxa (Sorensen = 1): ", n_saturated, " of ",
        length(D_assoc))
low_connect <- tibble(plot_id = pa_assoc$plot_id[keep_a],
                      share_plots_with_shared_taxa =
                        rowMeans(as.matrix(D_assoc) < 1) - 1 / sum(keep_a)) |>
  arrange(share_plots_with_shared_taxa)
readr::write_csv(low_connect, file.path(out_dir, "S12_plot_connectivity.csv"))
D_display <- tryCatch(
  vegan::stepacross(D_assoc, path = "extended", toolong = 1, trace = FALSE),
  error = function(e) { warning("stepacross failed: ", conditionMessage(e)); D_assoc })
nmds <- vegan::metaMDS(D_display, k = 2, trymax = 200, autotransform = FALSE, trace = FALSE)
if (!isTRUE(nmds$converged)) warning("2-D NMDS did not converge; report stress and try more starts.")

# Change 25 (10 Oct 2026): Fig. 4 shows axes 1 and 2 of a THREE-dimensional
# NMDS (2-D stress 0.24; 3-D stress 0.18 in S14). More random starts, then up
# to three restarts from the best solution, to reach a converged solution.
nmds_fig <- vegan::metaMDS(D_display, k = 3, trymax = 500, autotransform = FALSE, trace = FALSE)
restarts <- 0
while (!isTRUE(nmds_fig$converged) && restarts < 3) {
  restarts <- restarts + 1
  nmds_fig <- vegan::metaMDS(D_display, k = 3, trymax = 500, autotransform = FALSE,
                             trace = FALSE, previous.best = nmds_fig)
}
if (!isTRUE(nmds_fig$converged)) {
  warning("3-D NMDS did not converge after ", restarts, " restarts; report this in the caption.")
}
message("Fig. 4 NMDS: k = 3, stress ", round(nmds_fig$stress, 3),
        ", converged = ", isTRUE(nmds_fig$converged), ", restarts = ", restarts)
nmds_df <- as_tibble(vegan::scores(nmds_fig, display = "sites")[, 1:2]) |>
  bind_cols(env_a) |> mutate(ft = factor(ft, levels = forest_types))

# One ordination of all plots, shown in one panel per site so that the four
# forest types can be told apart (eight overlapping ellipses in a single panel
# were unreadable). Axes are shared between panels, so positions are comparable.
nmds_linetype <- c(native = "dashed", stocked = "solid", unstocked = "dotdash", mixed = "dotted")
fig_nmds <- ggplot(nmds_df |> mutate(site = factor(site, levels = sites)),
                   aes(NMDS1, NMDS2)) +
  stat_ellipse(aes(linetype = ft, group = ft), level = 0.68, colour = "grey20") +
  geom_point(aes(fill = ft), shape = 21, size = 2.4, colour = "black", stroke = 0.4) +
  scale_fill_manual(values = ft_fill, labels = ft_labels, name = "Forest type") +
  scale_linetype_manual(values = nmds_linetype, labels = ft_labels, name = "Forest type") +
  facet_wrap(~ site, labeller = as_labeller(c(NKR = "(a) Nukurua", NHL = "(b) Nausori Highlands"))) +
  # Stress is reported in the caption, not on the figure (Pacific Science:
  # do not repeat caption information on the illustration).
  guides(fill = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  coord_equal() + theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.key.width = grid::unit(10, "mm"))
# If one or two points still sit far from the rest, list them from
# S12_plot_connectivity.csv and name them in the caption rather than dropping
# them silently.

# =============================================================================
# 4b. Reviewer-risk analyses (10 Oct 2026): tree-only richness (T9), the
#     stem-density covariate (T10), atypical-plot sensitivity (S13) and a
#     three-dimensional NMDS (S14)
# =============================================================================

# ---- helper: choose the family exactly as for the primary models -------------
# Poisson first; if DHARMa shows mis-dispersion (P < 0.05), try generalised
# Poisson and keep it only if the fit is stable.
choose_family <- function(r, dat = plot_metrics, rhs = "site * ft + (1 | transect_id)") {
  m <- fit_glmm(r, poisson(), dat = dat, rhs = rhs)
  if (!fit_ok(m)) stop(r, ": Poisson fit not usable")
  d <- check_dispersion(m, r)
  fam <- poisson()
  if (d$dispersion_p < 0.05) {
    m2 <- fit_glmm(r, glmmTMB::genpois(), dat = dat, rhs = rhs)
    if (fit_ok(m2)) fam <- glmmTMB::genpois()
  }
  list(family = fam, dispersion = d)
}

# ---- helper: interaction test, then additive or within-site contrasts --------
# Same decision rule as the primary analysis: no interaction (P >= 0.05) ->
# forest-type LRT in the additive model and Tukey contrasts averaged over sites;
# interaction -> Tukey contrasts within each site.
richness_tests <- function(r, fam, dat = plot_metrics, label = r) {
  m_int  <- fit_glmm(r, fam, dat = dat)
  m_add  <- fit_glmm(r, fam, dat = dat, rhs = "site + ft + (1 | transect_id)")
  m_site <- fit_glmm(r, fam, dat = dat, rhs = "site + (1 | transect_id)")
  if (!all(vapply(list(m_int, m_add, m_site), fit_ok, logical(1)))) {
    return(tibble(response = label, note = "model failed"))
  }
  a_int <- anova(m_add, m_int)
  a_ft  <- anova(m_site, m_add)
  p_int <- a_int$`Pr(>Chisq)`[2]
  head_rows <- tibble(
    response = label, family = fam$family,
    what = c("site x forest type LRT", "forest-type LRT (additive vs site-only)"),
    chisq = c(a_int$Chisq[2], a_ft$Chisq[2]), df = c(a_int$`Chi Df`[2], a_ft$`Chi Df`[2]),
    p.value = c(p_int, a_ft$`Pr(>Chisq)`[2]))
  if (p_int >= 0.05) {
    emm <- emmeans::emmeans(m_add, ~ ft, type = "response")
    scope <- "averaged over sites"
  } else {
    emm <- emmeans::emmeans(m_int, ~ ft | site, type = "response")
    scope <- "within site"
  }
  bind_rows(
    head_rows,
    as_tibble(summary(emm)) |>
      rename_with(~ "estimate", .cols = any_of(c("rate", "response", "emmean"))) |>
      mutate(response = label, what = paste("estimated mean,", scope)),
    as_tibble(summary(emmeans::contrast(emm, method = "pairwise", adjust = "tukey"),
                      infer = TRUE)) |>
      mutate(response = label, what = paste("contrast,", scope)))
}

# ---- T9: tree-only richness (trees >= 10 cm DBH, 100 m2 per plot) ------------
tree_responses <- c("tree_rich", "tree_rich_noMah")
tree_family <- purrr::map(tree_responses, choose_family) |> purrr::set_names(tree_responses)
t9 <- purrr::imap_dfr(tree_family, ~ richness_tests(.y, .x$family) |>
                        mutate(dispersion_ratio_poisson = .x$dispersion$dispersion_ratio,
                               dispersion_p_poisson     = .x$dispersion$dispersion_p))
readr::write_csv(t9, file.path(out_dir, "T9_tree_richness.csv"))

# ---- T10: does stem density account for the forest-type difference? ---------
# Richness is modelled with and without log(number of individuals) as a
# covariate, on the same plots. Plots with no individuals (richness 0) are
# dropped because log(0) is undefined; the number dropped is reported.
# NOTE: density may itself be a consequence of forest type (a mediator), so the
# covariate model shows how much of the forest-type difference is ASSOCIATED
# with density, not how much is caused by it.
density_pairs <- tribble(
  ~response,           ~count_var,      ~family_from,
  "woody_rich",        "n_ind",         "woody_rich",
  "woody_rich_noMah",  "n_ind_noMah",   "woody_rich_noMah",
  "tree_rich",         "n_stems",       "tree_rich",
  "tree_rich_noMah",   "n_stems_noMah", "tree_rich_noMah")
family_lookup <- c(model_family, purrr::map(tree_family, "family"))

ratio_vs_snf <- function(m) {
  emmeans::emmeans(m, ~ ft, type = "response") |>
    emmeans::contrast(method = "trt.vs.ctrl", ref = 1, adjust = "dunnett") |>  # ref = native
    summary(infer = TRUE) |> as_tibble() |>
    transmute(contrast = as.character(contrast), ratio, lower = asymp.LCL,
              upper = asymp.UCL, p.value)
}

t10 <- purrr::pmap_dfr(density_pairs, function(response, count_var, family_from) {
  fam <- family_lookup[[family_from]]
  dat <- plot_metrics |>
    mutate(n_cov = .data[[count_var]]) |>
    filter(n_cov > 0) |>
    mutate(log_n = log(n_cov))
  n_dropped <- nrow(plot_metrics) - nrow(dat)
  f <- function(rhs) fit_glmm(response, fam, dat = dat, rhs = rhs)
  m0  <- f("site + ft + (1 | transect_id)");         m0s <- f("site + (1 | transect_id)")
  m1  <- f("site + ft + log_n + (1 | transect_id)"); m1s <- f("site + log_n + (1 | transect_id)")
  if (!all(vapply(list(m0, m0s, m1, m1s), fit_ok, logical(1)))) {
    return(tibble(response = response, note = "model failed"))
  }
  lrt <- function(a, b) { x <- anova(a, b); c(chisq = x$Chisq[2], p = x$`Pr(>Chisq)`[2]) }
  l0 <- lrt(m0s, m0); l1 <- lrt(m1s, m1)
  b  <- glmmTMB::fixef(m1)$cond["log_n"]
  # Wald 95 % CI from the conditional-model covariance matrix
  se <- sqrt(vcov(m1)$cond["log_n", "log_n"])
  ci <- unname(b) + c(-1.96, 1.96) * se
  # How strongly is log(individuals) tied to forest type? R2 of a linear model
  # of log_n on site + forest type (high R2 = density and forest type confounded)
  r2_ft <- summary(lm(log_n ~ site + ft, data = dat))$r.squared
  left_join(ratio_vs_snf(m0) |> rename_with(~ paste0(.x, "_without_cov"), -contrast),
            ratio_vs_snf(m1) |> rename_with(~ paste0(.x, "_with_cov"), -contrast),
            by = "contrast") |>
    mutate(response = response, count_variable = count_var, family = fam$family,
           plots_used = nrow(dat), plots_dropped_zero = n_dropped,
           ft_LRT_chisq_without_cov = l0["chisq"], ft_LRT_p_without_cov = l0["p"],
           ft_LRT_chisq_with_cov = l1["chisq"], ft_LRT_p_with_cov = l1["p"],
           log_n_coef = unname(b), log_n_lower = unname(ci[1]), log_n_upper = unname(ci[2]),
           R2_logn_on_site_ft = r2_ft,
           # share of the excess richness (ratio - 1) that remains after adjustment
           attenuation_pct = 100 * (1 - (ratio_with_cov - 1) / (ratio_without_cov - 1)),
           .before = 1)
})
readr::write_csv(t10, file.path(out_dir, "T10_density_covariate.csv"))

# ---- S13: sensitivity without the atypical plots -----------------------------
# Three NHL unstocked plots hold only Pinus caribaea (apart from mahogany) and
# one NKR SNF plot holds only fruit trees typical of a garden (jackfruit,
# coconut, Delonix regia, mango). Repeat the main tests without them.
atypical_plots <- c("NHL_unstocked_1_1", "NHL_unstocked_1_2", "NHL_unstocked_3_4",
                    "NKR_native_2_5")
missing_ids <- setdiff(atypical_plots, plot_metrics$plot_id)
if (length(missing_ids) > 0) stop("Atypical plot IDs not found: ", paste(missing_ids, collapse = ", "))
pm_s13 <- filter(plot_metrics, !plot_id %in% atypical_plots)

s13_rich <- purrr::map_dfr(c("woody_rich", "woody_rich_noMah", "tree_rich", "tree_rich_noMah"),
  function(r) richness_tests(r, family_lookup[[r]], dat = pm_s13))
s13_comp <- purrr::map_dfr(c(FALSE, TRUE), function(drop_mah) {
  pa <- make_pa(long, drop_mah) |> filter(!plot_id %in% atypical_plots)
  X  <- pa |> select(-c(plot_id, site, ft, transect, transect_id)) |> as.matrix()
  keep <- rowSums(X) > 0
  D  <- vegan::vegdist(X[keep, ], method = "bray", binary = TRUE)
  meta <- pa[keep, ] |> transmute(site, ft, tr_key = as.character(transect_id))
  splitplot_permanova(D, meta) |>
    mutate(community = if (drop_mah) "associated (S. macrophylla excluded)" else "whole community",
           n_plots = sum(keep))
})
readr::write_csv(bind_rows(s13_rich |> mutate(part = "richness GLMMs"),
                           s13_comp |> mutate(part = "split-plot PERMANOVA")),
                 file.path(out_dir, "S13_sensitivity_without_atypical_plots.csv"))

# ---- S16: tree-only richness under the two earlier robustness checks ---------
# Tree taxa per plot became the primary richness response (author, 10 Oct
# 2026), so the two checks run earlier for pooled woody richness only are
# repeated for all four richness responses with the same decision rule as T9:
#   (a) without the first (access-point) plot of each transect (cf. S10);
#   (b) all 120 plots, i.e. including NHL native transect 2 (cf. S3).
# Pooled responses are included so the table can be checked against S10/S3.
s16_responses <- c("tree_rich", "tree_rich_noMah", "woody_rich", "woody_rich_noMah")
s16_sets <- list(
  "without first plot of each transect" = plot_metrics_no1,
  "all 120 plots (incl. NHL native transect 2)" = plot_metrics_all)
missing_cols <- setdiff(s16_responses, names(plot_metrics_all))
if (length(missing_cols) > 0) stop("S16: columns missing from plot_metrics_all: ",
                                   paste(missing_cols, collapse = ", "))
s16 <- purrr::imap_dfr(s16_sets, function(dat, set_name) {
  purrr::map_dfr(s16_responses, function(r) {
    tryCatch(
      richness_tests(r, family_lookup[[r]], dat = dat) |>
        mutate(data_set = set_name, n_plots = nrow(dat)),
      error = function(e) tibble(response = r, data_set = set_name,
                                 note = paste("failed:", conditionMessage(e))))
  })
})
readr::write_csv(s16, file.path(out_dir, "S16_tree_richness_robustness.csv"))

# ---- S17: category-free check - richness against mahogany share of basal area
# FHCL's criteria for the stocked / unstocked / mixed labels are not available,
# so we ask whether tree richness in plantation plots tracks a measured
# structural variable instead of the labels. Plantation plots only (SNF
# excluded); plots without trees (share undefined) are dropped and counted.
# Richness excluding mahogany is the main response, so the planted species
# itself cannot create the relationship.
pm_plant <- plot_metrics |>
  filter(ft != "native") |>
  mutate(mah_share = mah_ba_share, ft = droplevels(ft))
n_dropped_s17 <- sum(!is.finite(pm_plant$mah_share))
pm_plant <- filter(pm_plant, is.finite(mah_share))
s17 <- purrr::map_dfr(c("tree_rich_noMah", "tree_rich"), function(r) {
  tryCatch({
    fam <- family_lookup[[r]]
    m_share <- fit_glmm(r, fam, dat = pm_plant, rhs = "site + mah_share + (1 | transect_id)")
    m_label <- fit_glmm(r, fam, dat = pm_plant, rhs = "site + ft + (1 | transect_id)")
    m_site  <- fit_glmm(r, fam, dat = pm_plant, rhs = "site + (1 | transect_id)")
    if (!all(vapply(list(m_share, m_label, m_site), fit_ok, logical(1)))) stop("fit failed")
    a_share <- anova(m_site, m_share)
    a_label <- anova(m_site, m_label)
    b  <- glmmTMB::fixef(m_share)$cond[["mah_share"]]
    se <- sqrt(vcov(m_share)$cond["mah_share", "mah_share"])
    tibble(response = r, n_plots = nrow(pm_plant), n_dropped_no_trees = n_dropped_s17,
           share_coef = b, share_se = se,
           rate_ratio_per_0.1 = exp(0.1 * b),
           rr_lower = exp(0.1 * (b - 1.96 * se)), rr_upper = exp(0.1 * (b + 1.96 * se)),
           share_LRT_chisq = a_share$Chisq[2], share_LRT_df = a_share$`Chi Df`[2],
           share_LRT_p = a_share$`Pr(>Chisq)`[2],
           label_LRT_chisq = a_label$Chisq[2], label_LRT_df = a_label$`Chi Df`[2],
           label_LRT_p = a_label$`Pr(>Chisq)`[2],
           AIC_share = AIC(m_share), AIC_label = AIC(m_label), AIC_site_only = AIC(m_site))
  }, error = function(e) tibble(response = r, note = paste("failed:", conditionMessage(e))))
})
readr::write_csv(s17, file.path(out_dir, "S17_richness_vs_mahogany_share.csv"))

# ---- S17b: does the mahogany-share effect survive the number of stems? ------
# Mahogany and non-mahogany basal area share one 0.01-ha plot, so a high
# mahogany share may simply mean fewer other stems. Refit with log(stems) as a
# covariate: non-mahogany stems for richness excluding mahogany, all stems for
# all tree taxa. Plots with no counted stems are dropped (log(0) undefined).
s17b_pairs <- list(tree_rich_noMah = "n_stems_noMah", tree_rich = "n_stems")
s17b <- purrr::imap_dfr(s17b_pairs, function(n_col, r) {
  tryCatch({
    fam <- family_lookup[[r]]
    d <- pm_plant |>
      filter(.data[[n_col]] > 0) |>
      mutate(log_n = log(.data[[n_col]]))
    m_n       <- fit_glmm(r, fam, dat = d, rhs = "site + log_n + (1 | transect_id)")
    m_n_share <- fit_glmm(r, fam, dat = d, rhs = "site + log_n + mah_share + (1 | transect_id)")
    if (!all(vapply(list(m_n, m_n_share), fit_ok, logical(1)))) stop("fit failed")
    a  <- anova(m_n, m_n_share)
    b  <- glmmTMB::fixef(m_n_share)$cond[["mah_share"]]
    se <- sqrt(vcov(m_n_share)$cond["mah_share", "mah_share"])
    bn <- glmmTMB::fixef(m_n_share)$cond[["log_n"]]
    tibble(response = r, stem_count = n_col, n_plots = nrow(d),
           n_dropped_no_stems = nrow(pm_plant) - nrow(d),
           share_coef_adj = b, share_se_adj = se,
           rate_ratio_per_0.1_adj = exp(0.1 * b),
           rr_lower_adj = exp(0.1 * (b - 1.96 * se)), rr_upper_adj = exp(0.1 * (b + 1.96 * se)),
           share_LRT_chisq_adj = a$Chisq[2], share_LRT_p_adj = a$`Pr(>Chisq)`[2],
           log_n_coef = bn,
           cor_share_logn = cor(d$mah_share, d$log_n),
           AIC_logn = AIC(m_n), AIC_logn_share = AIC(m_n_share))
  }, error = function(e) tibble(response = r, note = paste("failed:", conditionMessage(e))))
})
readr::write_csv(s17b, file.path(out_dir, "S17b_mahogany_share_given_stems.csv"))

# ---- S14: three-dimensional NMDS (does stress fall to an acceptable level?) --
nmds3 <- nmds_fig  # the 3-D solution drawn in Fig. 4 (section 4, change 25)
s14 <- tibble(k = c(2, 3), stress = c(nmds$stress, nmds3$stress),
              converged = c(isTRUE(nmds$converged), isTRUE(nmds3$converged)))
readr::write_csv(s14, file.path(out_dir, "S14_nmds_stress_k2_k3.csv"))
message("NMDS stress: k = 2 ", round(nmds$stress, 3), "; k = 3 ", round(nmds3$stress, 3))

# ---- Software citations and versions for the reference list ------------------
sink(file.path(out_dir, "S15_software_citations.txt"))
print(R.version.string)
for (p in c("glmmTMB", "DHARMa", "emmeans", "vegan")) {
  cat("\n====", p, as.character(utils::packageVersion(p)), "====\n")
  print(utils::citation(p), style = "text")
}
cat("\n==== R ====\n"); print(utils::citation(), style = "text")
sink()

# =============================================================================
# 5. Species overlap, origin and life-form summaries (R1 major 1; W5)
# =============================================================================
membership <- long |> distinct(site, ft, botanical) |> mutate(present = TRUE)
overlap_tbl <- membership |>
  group_by(site) |>
  group_modify(function(d, key) {
    sets <- split(d$botanical, d$ft)
    purrr::map_dfr(forest_types, function(f) {
      others <- unique(unlist(sets[setdiff(forest_types, f)]))
      tibble(ft = f, n_taxa = length(sets[[f]]),
             unique_to_type = length(setdiff(sets[[f]], others)),
             shared_with_native = if (f == "native") NA_integer_ else length(intersect(sets[[f]], sets$native)),
             jaccard_vs_native  = if (f == "native") NA_real_ else
               length(intersect(sets[[f]], sets$native)) / length(union(sets[[f]], sets$native)))
    }) |> add_row(ft = "shared by all four", n_taxa = length(Reduce(intersect, sets)))
  }) |> ungroup()
readr::write_csv(overlap_tbl, file.path(out_dir, "T6_species_overlap.csv"))

fig_upset <- membership |>
  mutate(ft = ft_labels[as.character(ft)]) |>
  group_by(site, botanical) |> summarise(types = list(ft), .groups = "drop") |>
  ggplot(aes(x = types)) +
  geom_bar(fill = "grey30") +
  ggupset::scale_x_upset(n_intersections = 15) +
  facet_wrap(~ site, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(NKR = "(a) Nukurua", NHL = "(b) Nausori Highlands"))) +
  scale_y_continuous(breaks = scales::breaks_width(2)) +   # whole numbers of taxa
  labs(x = NULL, y = "Number of taxa") + theme_classic(base_size = 10)

origin_tbl <- long |>
  group_by(site, ft) |>
  summarise(
    taxa_native     = n_distinct(botanical[stringr::str_starts(origin_sourced, "native")]),
    taxa_exotic     = n_distinct(botanical[exotic]),
    taxa_unresolved = n_distinct(botanical[origin_sourced %in% c("unknown", "check")]),
    stem_share_exotic     = mean(exotic[layer == "tree"]),
    seedling_share_exotic = sum(n[layer == "seedling" & exotic]) / sum(n[layer == "seedling"]),
    .groups = "drop")
readr::write_csv(origin_tbl, file.path(out_dir, "T7_origin_summary.csv"))
# origin_sourced: Fiji native/exotic status from POWO distribution data
# (accessed 30 Sep 2026; source per taxon in column origin_source). Early
# Polynesian introductions (Aleurites, Artocarpus altilis, Cocos, Syzygium
# malaccense) are counted as exotic, following POWO. The successional-guild
# column (guild_to_source) is still to be completed by the author.

# =============================================================================
# 6. Exotics and major taxa quantified (R2 #11; R1 L376, L539-548)
# =============================================================================
# Densities use the number of plots actually retained per group (NHL SNF = 10)
species_tbl <- long |>
  left_join(n_plots_group, by = c("site", "ft")) |>
  group_by(site, ft, botanical, origin_sourced, n_plots) |>
  summarise(stems_ha     = sum(layer == "tree") / (first(n_plots) * plot_area_ha),
            ba_m2ha      = sum(pi * (dbh_cm[layer == "tree"] / 200)^2, na.rm = TRUE) / (first(n_plots) * plot_area_ha),
            plot_freq    = n_distinct(plot_id),
            seedlings_m2 = sum(n[layer == "seedling"]) / (first(n_plots) * seed_area_m2),
            .groups = "drop") |>
  arrange(site, ft, desc(stems_ha))
readr::write_csv(species_tbl, file.path(out_dir, "S5_species_by_foresttype.csv"))

key_exotics <- c(mahogany, "Maesopsis eminii", "Spathodoe campanulata",
                 "Pinus caribaea", "Cordia alliodora", "Eucalyptus deglupta", "Clidemia hirta")
readr::write_csv(filter(species_tbl, botanical %in% key_exotics),
                 file.path(out_dir, "T8_key_exotics.csv"))

# =============================================================================
# 7. Species accumulation: sample-based and individual-based (R1 L447; R2 #9)
# =============================================================================
abund_mat <- function(dat, s, f) {
  dat |> filter(site == s, ft == f) |>
    group_by(plot_id, botanical) |> summarise(n = sum(n), .groups = "drop") |>
    pivot_wider(names_from = botanical, values_from = n, values_fill = 0) |>
    select(-plot_id) |> as.matrix()
}
sac_df <- tidyr::expand_grid(site = sites, ft = forest_types) |>
  purrr::pmap_dfr(function(site, ft) {
    M <- abund_mat(long, site, ft)
    sb <- vegan::specaccum((M > 0) * 1, method = "random", permutations = 999)
    ib_n <- unique(round(seq(1, sum(M), length.out = 50)))
    ib <- vegan::rarefy(colSums(M), sample = ib_n, se = TRUE)
    bind_rows(
      tibble(site, ft, type = "Sample-based (plots)", x = sb$sites, S = sb$richness, se = sb$sd),
      tibble(site, ft, type = "Individual-based", x = ib_n, S = ib[1, ], se = ib[2, ]))
  }) |> mutate(ft = factor(ft, levels = forest_types))
readr::write_csv(sac_df, file.path(out_dir, "S6_accumulation_curves.csv"))

# Sites in the same order as the other figures (NKR first); plots panel on the
# left so the axis title reads correctly (the first export had it reversed).
# SNF is drawn as a heavier line so it stands out from the plantation types
# (R1, Fig. 5), and mixed stands are dotted so they are not confused with it.
sac_plot <- sac_df |>
  mutate(site = factor(site, levels = sites),
         type = factor(type, levels = c("Sample-based (plots)", "Individual-based")))
sac_linetype <- c(native = "dashed", stocked = "solid", unstocked = "dotdash", mixed = "dotted")
fig_sac <- ggplot(sac_plot, aes(x, S, linetype = ft)) +
  geom_ribbon(aes(ymin = pmax(S - 1.96 * se, 0), ymax = S + 1.96 * se, group = ft),
              fill = "grey85", alpha = 0.35) +
  geom_line(aes(linewidth = ft), colour = "black") +
  scale_linetype_manual(values = sac_linetype, labels = ft_labels, name = "Forest type") +
  scale_linewidth_manual(values = c(native = 1.0, stocked = 0.5, unstocked = 0.5, mixed = 0.6),
                         labels = ft_labels, name = "Forest type") +
  facet_grid(site ~ type, scales = "free_x") +
  labs(x = "Plots (left) or individuals (right) sampled", y = "Expected number of taxa") +
  # Two-row legend so the last forest type is not cut off at 140 mm width.
  guides(linetype = guide_legend(nrow = 2), linewidth = guide_legend(nrow = 2)) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom", legend.key.width = grid::unit(12, "mm"))

# =============================================================================
# 8. Richness figure with compact letters (replaces asterisks)
# =============================================================================
# Change 23 (10 Oct 2026): tree-only richness is the primary richness
# response, so panel (a) is now tree taxa per plot (T9). Pooled woody richness
# (trees + seedlings) becomes a secondary result, drawn as Online Resource
# Figure S4 with its own letters (T3b).
panel_labels <- c(tree_rich  = "(a) Tree taxa per plot",
                  seed_rich  = "(b) Seedling taxa per plot",
                  seed_ab    = "(c) Seedlings per m\u00b2",
                  woody_rich = "Woody taxa per plot (trees and seedlings)")
rich_long <- plot_metrics |>
  select(site, ft, all_of(names(panel_labels))) |>
  pivot_longer(all_of(names(panel_labels)), names_to = "response") |>
  mutate(site = factor(site, levels = sites))

# Letters must come from the SAME contrasts that the text reports:
#   * no site x forest type interaction (P >= 0.05) -> contrasts averaged
#     over sites from the additive model (T3b for woody taxa and seedling
#     density, T9 for tree taxa); the same letters at both sites;
#   * interaction present -> contrasts within site (T3).
letters_from <- function(d) {
  pv <- setNames(d$p.value, stringr::str_replace(as.character(d$contrast), "\\s+[-/]\\s+", "-"))
  lt <- multcompView::multcompLetters(pv, threshold = 0.05)$Letters
  tibble(ft = factor(names(lt), levels = forest_types), .group = unname(lt))
}
no_interaction <- lrt_interaction |> filter(!is.na(p), p >= 0.05) |> pull(response)
cld_additive <- additive_tbl |>
  filter(what == "contrast (averaged over sites)", response %in% no_interaction) |>
  group_by(response) |>
  group_modify(~ letters_from(.x)) |>
  ungroup()
# Tree taxa: T9 interaction LRT decides the scope of the letters.
t9_p_int <- t9 |> filter(response == "tree_rich", what == "site x forest type LRT") |> pull(p.value)
if (length(t9_p_int) != 1 || is.na(t9_p_int)) stop("T9 interaction test for tree_rich not found")
if (t9_p_int < 0.05) warning("tree_rich has a site x forest type interaction: letters are averaged over sites but should be within site")
cld_tree <- t9 |>
  filter(response == "tree_rich", what == "contrast, averaged over sites") |>
  letters_from() |>
  mutate(response = "tree_rich")
cld_fig <- bind_rows(
  cld_tbl |> filter(!response %in% no_interaction) |> mutate(site = factor(site, levels = sites)),
  bind_rows(cld_additive, cld_tree) |> tidyr::crossing(site = factor(sites, levels = sites)))
readr::write_csv(cld_fig, file.path(out_dir, "S2b_compact_letters_as_plotted.csv"))

# Draw one richness figure for any set of responses, with letters where available.
draw_richness <- function(resp, ncol_lab = TRUE) {
  d <- rich_long |>
    filter(response %in% resp) |>
    mutate(response = factor(panel_labels[response], levels = unname(panel_labels[resp])))
  p <- ggplot(d, aes(ft, value, fill = ft)) +
    geom_boxplot(outlier.shape = 21, colour = "black", width = 0.6) +
    scale_fill_manual(values = ft_fill, labels = ft_labels, guide = "none") +
    scale_x_discrete(labels = c(native = "SNF", stocked = "Stocked", unstocked = "Unstocked", mixed = "Mixed")) +
    facet_grid(response ~ site, scales = "free_y") +
    expand_limits(y = 0) +  # every panel starts at zero (panel a began at 1)
    labs(x = NULL, y = NULL) + theme_classic(base_size = 10)
  lab_pos <- d |> group_by(response) |> summarise(y = max(value) * 1.08, .groups = "drop")
  lt <- cld_fig |>
    filter(response %in% resp) |>
    mutate(response = factor(panel_labels[response], levels = unname(panel_labels[resp]))) |>
    left_join(lab_pos, by = "response")
  if (nrow(lt) > 0) p <- p + geom_text(data = lt, aes(ft, y, label = .group), inherit.aes = FALSE, size = 3)
  p
}
fig_rich   <- draw_richness(c("tree_rich", "seed_rich", "seed_ab"))  # Fig. 3
fig_pooled <- draw_richness("woody_rich")                           # Supplemental Material 3, Fig. S2

# =============================================================================
# 9. Retained from the previous script (corrected): frequency, rank-abundance,
#    seedling detection probability
# =============================================================================

# ---- 9a. Plot frequency (% of plots occupied) by forest type -----------------
# Denominator = plots actually sampled per site x forest type (n_plots_group),
# so plots with no records still count (previous script: 'Dead' included).
freq_tbl <- long |>
  distinct(site, ft, layer, botanical, plot_id) |>
  count(site, ft, layer, botanical, name = "plots_present") |>
  left_join(n_plots_group, by = c("site", "ft")) |>
  mutate(freq_pct = 100 * plots_present / n_plots)
readr::write_csv(freq_tbl, file.path(out_dir, "S7_frequency_by_foresttype.csv"))

# Figure S1: the 20 most frequent taxa per site and layer (mean frequency
# across forest types), dodged bars by forest type
freq_complete <- freq_tbl |>
  select(site, ft, layer, botanical, freq_pct) |>
  tidyr::complete(tidyr::nesting(site, layer, botanical), ft, fill = list(freq_pct = 0))
top_ids <- freq_complete |>
  group_by(site, layer, botanical) |>
  summarise(mean_freq = mean(freq_pct), .groups = "drop") |>
  group_by(site, layer) |>
  slice_max(mean_freq, n = 20, with_ties = FALSE) |>
  ungroup()
top_taxa <- freq_complete |>
  inner_join(top_ids, by = c("site", "layer", "botanical")) |>
  mutate(layer = factor(layer, c("tree", "seedling"),
                        c("Trees (DBH \u2265 10 cm)", "Seedlings (< 50 cm)")),
         botanical = forcats::fct_reorder(botanical, mean_freq))

# Labels use the corrected spellings (name_for_publication); zero-frequency
# bars are dropped so they do not leave stubs at the axis.
pub_names <- lookup |> distinct(botanical, name_for_publication) |>
  filter(!is.na(botanical)) |> distinct(botanical, .keep_all = TRUE)
top_taxa <- top_taxa |>
  mutate(botanical = as.character(botanical)) |>
  left_join(pub_names, by = "botanical") |>
  mutate(label = coalesce(name_for_publication, as.character(botanical)),
         # Kuasi is reported at genus level in the paper (Methods; POWO). Veitchia
         # sessilifolia is not relabelled because the data hold a separate
         # Veitchia sp.; it does not appear among the plotted taxa.
         label = dplyr::recode(label, "Podocarpus neriifolius" = "Podocarpus sp."),
         label = forcats::fct_reorder(label, mean_freq)) |>
  filter(freq_pct > 0)
fig_freq <- ggplot(top_taxa, aes(freq_pct, label, fill = ft)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.75,
           colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = ft_fill, labels = ft_labels, name = "Forest type") +
  facet_grid(layer ~ site, scales = "free_y", space = "free_y") +
  labs(x = "Frequency (% of plots)", y = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(face = "italic"), legend.position = "bottom")

# ---- 9b. Seedling rank-abundance curves (density per m2) ---------------------
rac_tbl <- long |>
  filter(layer == "seedling") |>
  group_by(site, ft, botanical) |>
  summarise(n = sum(n), .groups = "drop") |>
  left_join(n_plots_group, by = c("site", "ft")) |>
  mutate(density_m2 = n / (n_plots * seed_area_m2)) |>
  group_by(site, ft) |>
  arrange(desc(density_m2), .by_group = TRUE) |>
  mutate(rank = row_number()) |>
  ungroup()
readr::write_csv(rac_tbl, file.path(out_dir, "S8_seedling_rank_abundance.csv"))

fig_rac <- ggplot(rac_tbl, aes(rank, log10(density_m2 + 1), linetype = ft, shape = ft)) +
  geom_line(colour = "black", linewidth = 0.5) +
  geom_point(aes(fill = ft), size = 1.8, colour = "black", stroke = 0.3) +
  scale_linetype_manual(values = ft_linetype, labels = ft_labels, name = "Forest type") +
  scale_fill_manual(values = ft_fill, labels = ft_labels, name = "Forest type") +
  scale_shape_manual(values = c(native = 23, stocked = 22, unstocked = 24, mixed = 21),
                     labels = ft_labels, name = "Forest type") +
  facet_wrap(~ site) +
  labs(x = "Species rank", y = expression(log[10]*"(seedlings m"^-2*" + 1)")) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")

# ---- 9c. Seedling detection probability (Methods / Online Resource) ----------
# Under complete spatial randomness, P(detect) = 1 - (1 - f)^N, where f is the
# sampled fraction of the plot (1 m2 of 100 m2 = 0.01) and N the number of
# individuals of a species in the plot. Clumped distributions detect worse,
# so these values are upper bounds.
f_sampled <- seed_area_m2 / (plot_area_ha * 1e4)
N_for_p   <- function(target, f) log(1 - target) / log(1 - f)
detect_thresholds <- tibble(target = c(0.5, 0.9),
                            N_required = N_for_p(c(0.5, 0.9), f_sampled))
readr::write_csv(detect_thresholds, file.path(out_dir, "S9_detection_thresholds.csv"))

detect_df <- tibble(N = 1:300, P_detect = 1 - (1 - f_sampled)^(1:300))
fig_detect <- ggplot(detect_df, aes(N, P_detect)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = c(0.5, 0.9), linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = detect_thresholds$N_required, linetype = "dotted") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Individuals of a species in the 100 m² plot (N)",
       y = "Probability of detection\nin the 1 m² seedling subsample") +
  theme_classic(base_size = 10)

# =============================================================================
# 10. Export figures (Pacific Science: separate TIFF/EPS files, >= 300 dpi at
#     final size; widths 67 mm one column or 140 mm two columns; figure plus
#     caption <= 197 mm tall). TIFF 600 dpi greyscale + EPS.
# =============================================================================
save_fig <- function(p, name, w = 140, h = 150) {
  tryCatch({
    ggsave(file.path(out_dir, paste0(name, ".tiff")), p, width = w, height = h, units = "mm",
           dpi = 600, compression = "lzw", bg = "white")
    ggsave(file.path(out_dir, paste0(name, ".eps")), p, width = w, height = h, units = "mm",
           device = grDevices::cairo_ps, bg = "white")
  }, error = function(e) warning("Could not save ", name, ": ", conditionMessage(e)))
}
save_fig(fig_rich,  "Figure_4_richness", h = 165)
save_fig(fig_nmds,  "Figure_5_NMDS", h = 110)
save_fig(fig_upset, "Figure_3_species_overlap", h = 160)
save_fig(fig_sac,   "Figure_2_accumulation", h = 140)
save_fig(fig_freq,  "Figure_S1_frequency", h = 230)
save_fig(fig_rac,   "Figure_S3_seedling_rank_abundance", h = 100)
save_fig(fig_detect, "Figure_S4_detection_probability", w = 120, h = 90)
save_fig(fig_pooled, "Figure_S2_pooled_woody_richness", h = 80)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
message("Done. Outputs written to ", normalizePath(out_dir))
