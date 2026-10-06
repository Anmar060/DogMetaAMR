##############################################################################
## Chapter 4 (resistome) | Puppy faecal resistome and ARG host context
## Chapter3_master_analysis_v3.R  (open-science reference script)
##
## One fresh R session, one output root. Run from the top with source().
## Every number reported in the Chapter 4 Results and Appendix is produced here,
## and Section 11 writes 00_reproducibility_audit.csv, which compares each
## reproduced value against the value reported in the chapter.
##
## Supersedes, and merges:
##   Chapter3_master_analysis_v2.R / Chapter3_master_analysis_1.R  (cohort models,
##     within-dog tests, FLUPOUL/NanoMotif contig context, CUPID two-barcode)
##   Chapter4_AMR_resistome_analysis.R  (depth diagnostics, Good's coverage,
##     saturation model, categories, core/accessory, OXA vs MEG_4223,
##     surveillance panel, PCA/PERMANOVA, global and class-level models, Spearman)
##   methylation_analysis_v2.R  (motif counts, regime overlap, MAG quality)
##   CUPID_methylation_analysis_v3_rank_preserved.R  (rank-preserving CUPID)
##
## MARKER SCOPE (read this before interpreting any output)
##   "all"     = every AMR++/MEGARes gene_id, including "RequiresSNPConfirmation"
##               markers. PRIMARY scope: nothing is excluded from the dataset.
##   "non_snp" = gene_ids NOT annotated "RequiresSNPConfirmation". SENSITIVITY
##               scope, because SNP verification did not reduce any count and the
##               SNP-dependent markers include multi-copy rRNA genes whose counts
##               may partly track bacterial abundance.
##   Earlier scripts (Chapter4_AMR_resistome_analysis.R section 2b.ii and the
##   consolidated master) silently overwrote `amr` with the non-SNP subset, so
##   every analysis after that point dropped 91% of hits. Here no object is ever
##   overwritten: each analysis takes its scope explicitly and every cohort-level
##   output carries a `scope` column with both scopes.
##
## TERMINOLOGY: counts are ARG detections (gene hits), not phenotypic resistance.
##
## Required packages: tidyverse, janitor, readxl, glmmTMB, broom.mixed, vegan,
## patchwork, scales.
##############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(readxl)
  library(glmmTMB)
  library(broom.mixed)
  library(vegan)
  library(patchwork)
  library(scales)
})

# 0 | CONFIGURATION ------------------------------------------------------------

DATA_DIR <- path.expand("~/Documents/2110_dog_AMR/exploration")
OUT_DIR  <- file.path(DATA_DIR, "chapter4_thesis_outputs_v3")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TAB_DIR  <- file.path(OUT_DIR, "tables")
RDS_DIR  <- file.path(OUT_DIR, "rds")
invisible(lapply(c(OUT_DIR, FIG_DIR, TAB_DIR, RDS_DIR),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

MIN_READS      <- 50000L   # cohort QC threshold (same as Chapter 3)
SEED           <- 1L
BARCODES       <- c("barcode03", "barcode14", "barcode36", "barcode47")
CUPID_BARCODES <- c("barcode36", "barcode03")
SCOPES         <- c("all", "non_snp")

# Folders searched (recursively) for inputs whose location has varied between
# earlier scripts. The first match is used and recorded in 00_input_manifest.
SEARCH_ROOTS <- c(DATA_DIR,
                  file.path(DATA_DIR, "methylation_R_master"),
                  file.path(DATA_DIR, "contig_bins"),
                  file.path(DATA_DIR, "Chapter3", "Chapter3_analysis"))

find_input <- function(filename, required = TRUE) {
  for (root in SEARCH_ROOTS) {
    direct <- file.path(root, filename)
    if (file.exists(direct)) return(normalizePath(direct))
  }
  for (root in SEARCH_ROOTS[dir.exists(SEARCH_ROOTS)]) {
    hit <- list.files(root, pattern = paste0("^", gsub(".", "\\.", filename, fixed = TRUE), "$"),
                      recursive = TRUE, full.names = TRUE)
    if (length(hit)) return(normalizePath(hit[1]))
  }
  if (required) stop("Input not found under SEARCH_ROOTS: ", filename)
  NA_character_
}

inputs <- c(
  amr           = find_input("AMR_matrix_group.csv"),
  tech          = find_input("NanoStats_matrix_final.csv"),
  metadata      = find_input("litter_metadata.xlsx"),
  bracken       = find_input("merged_with_superkingdom_wave.csv"),
  runs          = find_input("methylation_run_manifest.tsv"),
  bins          = find_input("methylation_bin_master.tsv"),
  motifs        = find_input("methylation_motif_master_long.tsv"),
  inclusions    = find_input("methylation_include_contigs_master.tsv"),
  contamination = find_input("methylation_contamination_master.tsv"),
  abricate      = find_input("abricate_megares_master.tsv"),
  abricate_db   = find_input("abricate_database_versions.tsv", required = FALSE)
)
maps <- setNames(vapply(paste0(BARCODES, ".contig_bin.tsv"), find_input, ""), BARCODES)
cupid_files <- c("amrfinder.tsv", "ARG_read_methylation_master.tsv",
                 "plasmid_ARG_reads.tsv", "plasmid_host_top25.tsv",
                 "host_taxon_consensus.tsv", "plasmid_host_scores.tsv")
cupid_inputs <- unlist(lapply(CUPID_BARCODES, function(b)
  setNames(vapply(paste0(b, ".", cupid_files), find_input, "", required = FALSE),
           paste(b, c("amr", "overlap", "plasmids", "candidates", "consensus", "scores"),
                 sep = "_"))))
if (any(is.na(cupid_inputs[!grepl("_scores$", names(cupid_inputs))])))
  stop("Missing CUPID input(s): ",
       paste(names(cupid_inputs)[is.na(cupid_inputs) & !grepl("_scores$", names(cupid_inputs))],
             collapse = ", "))

# 0b | HELPERS -------------------------------------------------------------------

save_table <- function(x, stem) readr::write_csv(x, file.path(TAB_DIR, paste0(stem, ".csv")))
save_figure <- function(x, stem, width = 9, height = 6) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".tiff")), x, width = width, height = height,
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), x, width = width, height = height,
         bg = "white")
}
require_cols <- function(x, cols, label) {
  absent <- setdiff(cols, names(x))
  if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
  x
}
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  z <- tolower(trimws(as.character(x)))
  case_when(z %in% c("true", "t", "yes", "y", "1") ~ TRUE,
            z %in% c("false", "f", "no", "n", "0") ~ FALSE,
            TRUE ~ NA)
}
nonempty <- function(x) !is.na(x) & trimws(as.character(x)) != ""
quality_class <- function(completeness, contamination) {
  case_when(is.na(completeness) | is.na(contamination) ~ "QC unavailable",
            completeness >= 90 & contamination <= 5  ~ "High quality-like",
            completeness >= 50 & contamination <= 10 ~ "Medium quality-like",
            TRUE ~ "Low/problematic")
}
safe_signed_rank <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || all(x == 0)) return(NA_real_)
  tryCatch(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
           error = function(e) NA_real_)
}
fit_nb <- function(formula, data, label) {
  tryCatch(glmmTMB(formula, family = nbinom2(), data = data),
           error = function(e) {warning(label, ": ", conditionMessage(e)); NULL})
}
tidy_rate <- function(model, label) {
  if (is.null(model)) return(tibble(model = label, term = NA_character_))
  tryCatch(broom.mixed::tidy(model, effects = "fixed", component = "cond",
                             conf.int = TRUE, exponentiate = TRUE),
           error = function(e) broom.mixed::tidy(model, effects = "fixed",
                                                 component = "cond", exponentiate = TRUE)) |>
    mutate(model = label, .before = 1)
}
fit_stats <- function(fit, label) {
  if (is.null(fit)) return(tibble(model = label, n = NA_integer_, AIC = NA_real_,
                                    BIC = NA_real_, logLik = NA_real_))
  tibble(model = label, n = nobs(fit), AIC = AIC(fit), BIC = BIC(fit),
         logLik = as.numeric(logLik(fit)))
}

save_table(tibble(input = names(c(inputs, maps, cupid_inputs)),
                  path = unname(c(inputs, maps, cupid_inputs))), "00_input_manifest")

# 1 | INPUTS, SAMPLE QC AND SAMPLE TABLE ---------------------------------------

tech_all <- read_csv(inputs[["tech"]], show_col_types = FALSE) |>
  clean_names() |>
  require_cols(c("sample", "sample_id", "wave", "num_reads", "total_bases", "batch"),
               "NanoStats") |>
  mutate(across(c(sample, sample_id, batch), as.character),
         wave = as.integer(wave), num_reads = as.numeric(num_reads),
         total_bases = as.numeric(total_bases))
if (anyDuplicated(tech_all$sample)) stop("NanoStats has duplicate sample IDs")
sample_qc <- tech_all |>
  transmute(sample, sample_id, wave, num_reads, total_bases,
            included = !is.na(num_reads) & num_reads >= MIN_READS &
                       !is.na(total_bases) & total_bases > 0)
save_table(sample_qc, "01_cohort_sample_inclusion")
tech <- tech_all |> semi_join(filter(sample_qc, included), by = "sample")

# MEGARes v3 labels one gene (MEG_8700, SITABCD) with the misspelt category
# "Multi-compount"; it is recoded so that category totals sum to the cohort total.
amr_raw <- read_csv(inputs[["amr"]], show_col_types = FALSE) |>
  clean_names() |>
  require_cols(c("sample", "gene_id", "category", "subcategory", "family",
                 "gene_name", "group", "value"), "AMR matrix") |>
  mutate(sample = as.character(sample),
         value = suppressWarnings(as.numeric(value)),
         category_as_published = category,
         category = recode(category, "Multi-compount" = "Multi-compound"),
         snp_dependent = str_detect(coalesce(gene_name, ""), "RequiresSNPConfirmation"))
if (any(is.na(amr_raw$value)))
  stop("AMR matrix contains missing or non-numeric values; inspect before modelling")
save_table(amr_raw |> filter(category != category_as_published) |>
             distinct(gene_id, gene_name, category_as_published, category),
           "01_MEGARes_category_recodes")

# Genes present in the matrix but with zero hits in the retained samples
# (detected only in samples removed by the read filter) are excluded from every
# gene count, so "distinct ARGs detected" refers to the retained cohort only.
amr_cohort <- amr_raw |> semi_join(tech, by = "sample")
genes_zero_in_cohort <- amr_cohort |>
  group_by(gene_id, gene_name, snp_dependent) |>
  summarise(total = sum(value), .groups = "drop") |>
  filter(total == 0)
save_table(genes_zero_in_cohort, "01_genes_detected_only_in_excluded_samples")

amr_scope <- function(scope) {
  stopifnot(scope %in% SCOPES)
  if (scope == "all") amr_cohort else filter(amr_cohort, !snp_dependent)
}

meta <- read_excel(inputs[["metadata"]], sheet = "litter_metadata") |>
  clean_names() |>
  require_cols(c("sample", "q6_antibiotics", "q6_antibiotics_oral_previous",
                 "eats_raw", "q9_any_gi"), "litter_metadata") |>
  mutate(sample = as.character(sample))
if (anyDuplicated(meta$sample)) stop("Metadata has duplicate sample IDs")
meta <- meta |>
  transmute(sample, antibiotics = as_bool(q6_antibiotics),
            prev_antibiotics = as_bool(q6_antibiotics_oral_previous),
            eats_raw = as_bool(eats_raw), gi = as_bool(q9_any_gi))

sample_burden <- function(scope) {
  amr_scope(scope) |>
    group_by(sample) |>
    summarise(hits = sum(value),
              richness = n_distinct(gene_id[value > 0]),
              drug_hits = sum(value[group == "Drug_related"]),
              other_hits = sum(value[group != "Drug_related"]),
              singletons = sum(value == 1),
              .groups = "drop") |>
    rename_with(~paste0(.x, "_", scope), -sample)
}
samples <- tech |>
  left_join(meta, by = "sample") |>
  left_join(sample_burden("all"), by = "sample") |>
  left_join(sample_burden("non_snp"), by = "sample") |>
  mutate(across(matches("_(all|non_snp)$"), ~replace_na(.x, 0)),
         depth_gb = total_bases / 1e9,
         log10_total_bases = log10(total_bases),
         age = factor(wave, levels = 1:3, labels = c("4 months", "7 months", "12 months")))
for (s in SCOPES) {
  samples[[paste0("hits_per_gb_", s)]]      <- samples[[paste0("hits_", s)]] / samples$depth_gb
  samples[[paste0("hits_per_mreads_", s)]]  <- samples[[paste0("hits_", s)]] / samples$num_reads * 1e6
  samples[[paste0("goods_coverage_", s)]]   <- if_else(samples[[paste0("hits_", s)]] > 0,
    1 - samples[[paste0("singletons_", s)]] / samples[[paste0("hits_", s)]], NA_real_)
}
# Legacy column names from 02_cohort_sample_data.csv, retained only so old files can be
# compared: arg_hits / arg_richness / arg_hits_per_gb are NON-SNP-dependent values.
samples <- samples |>
  mutate(arg_hits = hits_non_snp, arg_richness = richness_non_snp,
         all_marker_hits = hits_all,
         arg_hits_per_gb = hits_per_gb_non_snp,
         all_marker_hits_per_gb = hits_per_gb_all)
save_table(samples, "01_cohort_sample_data")
saveRDS(samples, file.path(RDS_DIR, "cohort_samples.rds"))

cohort_summary <- samples |>
  summarise(n_samples = n(), n_dogs = n_distinct(sample_id),
            n_wave1 = sum(wave == 1), n_wave2 = sum(wave == 2), n_wave3 = sum(wave == 3),
            median_depth_gb = median(depth_gb), iqr_depth_gb = IQR(depth_gb),
            raw_fed = sum(eats_raw, na.rm = TRUE),
            antibiotics = sum(antibiotics, na.rm = TRUE),
            gi = sum(gi, na.rm = TRUE))
wave_summary <- samples |>
  group_by(age) |>
  summarise(n_samples = n(), n_dogs = n_distinct(sample_id),
            median_depth_gb = median(depth_gb),
            median_non_snp_hits_per_gb = median(hits_per_gb_non_snp),
            n_antibiotics = sum(antibiotics, na.rm = TRUE), .groups = "drop")
save_table(cohort_summary, "01_cohort_summary")
save_table(wave_summary, "01_wave_summary_tableA1")

# 2 | SEQUENCING-DEPTH SUFFICIENCY (both scopes) ------------------------------

bracken <- read_csv(inputs[["bracken"]], show_col_types = FALSE) |>
  clean_names() |>
  require_cols(c("sample", "superkingdom", "name", "fraction_total_reads"), "Bracken") |>
  mutate(sample = as.character(sample)) |>
  semi_join(tech, by = "sample")
kingdom_frac <- bracken |>
  group_by(sample, superkingdom) |>
  summarise(frac = sum(fraction_total_reads, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = superkingdom, values_from = frac, values_fill = 0)
samples_k <- samples |> left_join(kingdom_frac, by = "sample")

depth_diagnostics <- map_dfr(SCOPES, function(s) {
  h <- samples_k[[paste0("hits_", s)]]
  pg <- samples_k[[paste0("hits_per_gb_", s)]]
  pm <- samples_k[[paste0("hits_per_mreads_", s)]]
  rich <- samples_k[[paste0("richness_", s)]]
  good <- samples_k[[paste0("goods_coverage_", s)]]
  mm <- tryCatch(nls(rich ~ (Rmax * tb) / (K + tb),
                     data = tibble(rich = rich, tb = samples_k$total_bases),
                     start = list(Rmax = max(rich) * 1.5, K = median(samples_k$total_bases))),
                 error = function(e) NULL)
  cf <- if (is.null(mm)) c(Rmax = NA, K = NA) else coef(mm)
  tibble(scope = s,
         r_raw_vs_reads = cor(h, samples_k$num_reads),
         r_raw_vs_bases = cor(h, samples_k$total_bases),
         r_per_mreads_vs_reads = cor(pm, samples_k$num_reads),
         r_per_gb_vs_bases = cor(pg, samples_k$total_bases),
         r_per_gb_vs_reads = cor(pg, samples_k$num_reads),
         r_per_gb_vs_bacterial_fraction = if ("Bacteria" %in% names(samples_k))
           cor(pg, samples_k$Bacteria, use = "complete.obs") else NA_real_,
         mm_Rmax = unname(cf["Rmax"]),
         mm_depth_gb_90 = unname(0.90 * cf["K"] / 0.10 / 1e9),
         mm_depth_gb_95 = unname(0.95 * cf["K"] / 0.05 / 1e9),
         goods_median = median(good, na.rm = TRUE),
         goods_q1 = quantile(good, .25, na.rm = TRUE),
         goods_q3 = quantile(good, .75, na.rm = TRUE),
         goods_min = min(good, na.rm = TRUE),
         median_richness = median(rich))
})
save_table(depth_diagnostics, "02_depth_diagnostics_by_scope")

# 3 | RESISTOME COMPOSITION (scope: all markers unless stated) ----------------

n_samples <- n_distinct(samples$sample)
gene_table <- amr_cohort |>
  group_by(gene_id, category, subcategory, family, gene_name, group, snp_dependent) |>
  summarise(total_hits = sum(value), n_present = sum(value > 0), .groups = "drop") |>
  filter(total_hits > 0) |>
  mutate(prevalence = n_present / n_samples,
         compartment = case_when(prevalence >= .8 ~ "Core (>=80%)",
                                 prevalence >= .2 ~ "Intermediate (20-<80%)",
                                 TRUE ~ "Rare (<20%)"))
save_table(gene_table, "03_gene_prevalence_compartments")

snp_split <- gene_table |>
  group_by(snp_dependent) |>
  summarise(n_genes = n(), hits = sum(total_hits), .groups = "drop") |>
  mutate(pct_genes = 100 * n_genes / sum(n_genes), pct_hits = 100 * hits / sum(hits))
save_table(snp_split, "03_SNP_dependent_split_tableA4")

compartments <- gene_table |>
  group_by(snp_dependent, compartment) |>
  summarise(n_genes = n(), hits = sum(total_hits), .groups = "drop")
save_table(compartments, "03_compartments_by_marker_set_tableA3")
save_table(gene_table |> count(compartment, name = "n_genes") |>
             mutate(pct = 100 * n_genes / sum(n_genes)), "03_compartments_overall")
save_table(gene_table |> filter(compartment == "Core (>=80%)") |> arrange(desc(total_hits)),
           "03_core_genes")

category_by_sample <- amr_cohort |>
  group_by(sample, category) |>
  summarise(hits = sum(value), .groups = "drop") |>
  complete(sample = samples$sample, category, fill = list(hits = 0)) |>
  left_join(samples |> select(sample, depth_gb, age), by = "sample") |>
  mutate(hits_per_gb = hits / depth_gb)
category_summary <- category_by_sample |>
  group_by(category) |>
  summarise(total_hits = sum(hits), median_hits_per_gb = median(hits_per_gb),
            .groups = "drop") |>
  mutate(pct_hits = 100 * total_hits / sum(total_hits)) |>
  arrange(desc(total_hits))
save_table(category_summary, "03_category_totals_table4_4")

p_category <- category_by_sample |>
  group_by(age, category) |>
  summarise(median_per_gb = median(hits_per_gb), .groups = "drop") |>
  ggplot(aes(age, median_per_gb, fill = category)) +
  geom_col(position = position_dodge(width = .8)) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), labels = comma,
                     breaks = c(0, 1, 10, 100, 1000, 10000)) +
  labs(x = "Puppy age", y = "Median ARG hits per Gb (pseudo-log10)",
       fill = "Resistance category") +
  theme_bw(base_size = 11)
save_figure(p_category, "Fig_S_category_burden_by_age", 10, 6)

# Targeted comparison: OXA/Class D (non-SNP-dependent) vs MEG_4223/O23S (SNP-dependent).
# Both require the all-marker scope, because MEG_4223 is SNP-dependent.
target_amr <- bind_rows(
  amr_cohort |> filter(gene_id == "MEG_4223") |> mutate(target = "MEG_4223 / O23S"),
  amr_cohort |> filter(family == "Class_D_betalactamases") |>
    mutate(target = "OXA / Class D beta-lactamases"),
  amr_cohort |> filter(gene_id == "MEG_4741") |> mutate(target = "MEG_4741 (OXA) alone")
)
# Count distinct detected gene IDs across the retained cohort, rather than
# taking the maximum number detected in an individual sample.
target_gene_counts <- target_amr |>
  group_by(target) |>
  summarise(n_gene_ids = n_distinct(gene_id[value > 0]), .groups = "drop")
target_by_sample <- target_amr |>
  group_by(target, sample) |>
  summarise(hits = sum(value), n_gene_ids = n_distinct(gene_id[value > 0]), .groups = "drop") |>
  complete(target, sample = samples$sample, fill = list(hits = 0, n_gene_ids = 0L)) |>
  left_join(samples |> select(sample, depth_gb, age, eats_raw), by = "sample") |>
  mutate(hits_per_gb = hits / depth_gb, present = hits > 0)
target_summary <- target_by_sample |>
  group_by(target) |>
  summarise(total_hits = sum(hits),
            prevalence_pct = 100 * mean(present),
            median_hits_per_gb_where_detected = median(hits_per_gb[present]),
            median_hits_per_gb_all_samples = median(hits_per_gb), .groups = "drop") |>
  left_join(target_gene_counts, by = "target") |>
  relocate(n_gene_ids, .after = target)
save_table(target_summary, "03_OXA_vs_MEG4223_table4_5")
raw_target_tests <- target_by_sample |>
  filter(target != "MEG_4741 (OXA) alone", !is.na(eats_raw)) |>
  group_by(target, age) |>
  summarise(n_raw = sum(eats_raw), n_nonraw = sum(!eats_raw),
            p_raw = tryCatch(wilcox.test(hits_per_gb ~ eats_raw, exact = FALSE)$p.value,
                             error = function(e) NA_real_), .groups = "drop") |>
  mutate(p_raw_BH = p.adjust(p_raw, method = "BH"))
save_table(raw_target_tests, "03_OXA_MEG4223_raw_feeding_tests")

# Surveillance panel: descriptive prevalence and burden per wave (no hypothesis tests).
surveillance <- amr_cohort |>
  mutate(marker_set = case_when(
    family == "Class_D_betalactamases" ~ "OXA/Class D beta-lactamase",
    str_detect(gene_name, regex("^CTX", TRUE)) ~ "CTX-M",
    str_detect(gene_name, regex("^TEM", TRUE)) ~ "TEM",
    str_detect(gene_name, regex("^SHV", TRUE)) ~ "SHV",
    str_detect(gene_name, regex("^CMY", TRUE)) ~ "CMY (AmpC)",
    str_detect(gene_name, regex("^ACC", TRUE)) ~ "ACC (AmpC)",
    str_detect(gene_name, regex("^MCR", TRUE)) ~ "mcr",
    str_detect(gene_name, regex("^VAN", TRUE)) ~ "van",
    str_detect(gene_name, regex("^MEC", TRUE)) ~ "mec",
    str_detect(gene_name, regex("^(CFR|OPTRA|POXTA)", TRUE)) ~ "cfr/optrA/poxtA",
    str_detect(gene_name, regex("^NIM", TRUE)) ~ "nim",
    subcategory %in% c("Fluoroquinolones", "Sulfonamides", "Trimethoprim", "Phenicol") ~ subcategory,
    TRUE ~ NA_character_)) |>
  filter(!is.na(marker_set)) |>
  group_by(marker_set, sample, snp_dependent) |>
  summarise(hits = sum(value), .groups = "drop") |>
  left_join(samples |> select(sample, depth_gb, age), by = "sample") |>
  group_by(marker_set, snp_dependent, age) |>
  summarise(n_samples_present = sum(hits > 0),
            median_hits_per_gb_where_detected = median((hits / depth_gb)[hits > 0]),
            .groups = "drop")
save_table(surveillance, "03_surveillance_panel_by_wave")

# 4 | WITHIN-DOG PAIRED CHANGE (both scopes) -----------------------------------

interval_levels <- c("4 to 7 months", "7 to 12 months", "4 to 12 months")
make_pairs <- function(a, b, label, scope) {
  rate <- paste0("hits_per_gb_", scope)
  left <- samples |> filter(wave == a) |>
    transmute(sample_id, rate_a = .data[[rate]], abx_a = antibiotics)
  right <- samples |> filter(wave == b) |>
    transmute(sample_id, rate_b = .data[[rate]], abx_b = antibiotics)
  inner_join(left, right, by = "sample_id") |>
    mutate(interval = label, scope = scope,
           delta_log = log1p(rate_b) - log1p(rate_a),
           abx_either = case_when(abx_a %in% TRUE | abx_b %in% TRUE ~ TRUE,
                                  abx_a %in% FALSE & abx_b %in% FALSE ~ FALSE,
                                  TRUE ~ NA))
}
pairs <- map_dfr(SCOPES, function(s) bind_rows(
  make_pairs(1, 2, interval_levels[1], s),
  make_pairs(2, 3, interval_levels[2], s),
  make_pairs(1, 3, interval_levels[3], s))) |>
  mutate(interval = factor(interval, levels = interval_levels),
         marker_set = recode(scope, all = "All AMR++ markers",
                             non_snp = "Non-SNP-dependent markers"))
save_table(pairs, "04_within_dog_pairs")
paired_tests <- pairs |>
  group_by(scope, marker_set, interval) |>
  summarise(n_pairs = n(), median_delta_log = median(delta_log),
            p = safe_signed_rank(delta_log), .groups = "drop") |>
  group_by(scope) |>
  mutate(p_BH = p.adjust(p, method = "BH")) |>
  ungroup()
save_table(paired_tests, "04_within_dog_tests_tableA5")
abx_pair_test <- pairs |>
  filter(!is.na(abx_either)) |>
  group_by(scope, interval) |>
  summarise(n_exposed = sum(abx_either), n_unexposed = sum(!abx_either),
            p_abx = tryCatch(wilcox.test(delta_log ~ abx_either, exact = FALSE)$p.value,
                             error = function(e) NA_real_), .groups = "drop")
save_table(abx_pair_test, "04_antibiotic_pair_comparison_tableA6")

p_long <- ggplot(pairs, aes(interval, delta_log,
                            colour = factor(abx_either, levels = c(FALSE, TRUE)))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_jitter(width = .12, height = 0, alpha = .52, size = 1.5) +
  stat_summary(aes(group = interval), fun = median, geom = "point",
               shape = 95, size = 11, colour = "black") +
  facet_wrap(~ marker_set, ncol = 2) +
  scale_colour_manual(values = c("FALSE" = "#376C91", "TRUE" = "#C56E30"),
                      labels = c("No", "Yes"), na.value = "grey65") +
  labs(x = NULL, y = "Within-dog change in log(1 + ARG hits/Gb)",
       colour = "Antibiotics reported\nat either wave") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
save_figure(p_long, "Fig_4_3_longitudinal_ARG_change", 12, 5)

# 5 | MIXED-EFFECTS MODELS (both scopes) ---------------------------------------
# Primary model as in Table A2: negative binomial, depth offset in Gb (so the
# intercept is hits per Gb), dog identity as random intercept.

model_data <- samples |>
  filter(complete.cases(antibiotics, prev_antibiotics, eats_raw, gi, sample_id,
                        age, total_bases, batch)) |>
  mutate(across(c(antibiotics, prev_antibiotics, eats_raw, gi),
                ~factor(.x, levels = c(FALSE, TRUE))),
         batch = factor(batch))
covars <- "age + antibiotics + prev_antibiotics + eats_raw + gi"
model_specs <- tribble(
  ~outcome,       ~label,               ~adjust,
  "hits",         "total_burden",       "offset",
  "drug_hits",    "drug_related",       "offset",
  "other_hits",   "other_resistance",   "offset",
  "richness",     "richness",           "covariate",
  "hits",         "total_burden_batch", "offset_batch"
)
models <- list()
for (s in SCOPES) for (i in seq_len(nrow(model_specs))) {
  sp <- model_specs[i, ]
  y <- paste0(sp$outcome, "_", s)
  rhs <- switch(sp$adjust,
    offset       = paste(covars, "+ offset(log(depth_gb)) + (1 | sample_id)"),
    offset_batch = paste(covars, "+ batch + offset(log(depth_gb)) + (1 | sample_id)"),
    covariate    = paste(covars, "+ log10_total_bases + (1 | sample_id)"))
  nm <- paste(s, sp$label, sep = "__")
  models[[nm]] <- fit_nb(as.formula(paste(y, "~", rhs)), model_data, nm)
}
model_results <- imap_dfr(models, tidy_rate) |>
  separate(model, into = c("scope", "model"), sep = "__")
model_fit <- imap_dfr(models, fit_stats) |>
  separate(model, into = c("scope", "model"), sep = "__")
save_table(model_results, "05_mixed_model_rate_ratios_tableA2")
save_table(model_fit, "05_mixed_model_AIC_BIC")
saveRDS(models, file.path(RDS_DIR, "mixed_models.rds"))

# Class-level models (both scopes); classes present in fewer than 10 samples are skipped.
classes_to_model <- c("betalactams", "Tetracyclines", "MLS", "Aminoglycosides",
                      "Fluoroquinolones", "Sulfonamides", "Trimethoprim", "Phenicol",
                      "Metronidazole", "Oxazolidinone", "Glycopeptides",
                      "Cationic_antimicrobial_peptides", "Multi-drug_resistance",
                      "Drug_and_biocide_resistance", "Drug_and_biocide_and_metal_resistance")
class_models <- map_dfr(SCOPES, function(s) map_dfr(classes_to_model, function(cls) {
  dat <- amr_scope(s) |> filter(subcategory == cls) |>
    group_by(sample) |> summarise(class_hits = sum(value), .groups = "drop") |>
    right_join(model_data, by = "sample") |>
    mutate(class_hits = replace_na(class_hits, 0))
  n_present <- sum(dat$class_hits > 0)
  if (n_present < 10) return(tibble(scope = s, subcategory = cls, n_present = n_present,
                                    note = "skipped: present in <10 samples"))
  m <- fit_nb(as.formula(paste("class_hits ~", covars,
                               "+ offset(log(depth_gb)) + (1 | sample_id)")), dat, cls)
  tidy_rate(m, "class_model") |> mutate(scope = s, subcategory = cls,
                                        n_present = n_present, .before = 1)
}))
# BH correction across resistance classes within each scope and model term.
# Intercepts and rows without a model term retain NA adjusted p-values.
class_models <- class_models |>
  group_by(scope, term) |>
  mutate(p_BH = if_else(!is.na(term) & term != "(Intercept)",
                        p.adjust(p.value, method = "BH"), NA_real_)) |>
  ungroup()
save_table(class_models, "05_class_level_mixed_models")
save_table(class_models |> filter(term %in% c("age7 months", "age12 months",
                                              "antibioticsTRUE", "prev_antibioticsTRUE",
                                              "eats_rawTRUE", "giTRUE")),
           "05_class_level_key_terms")

# 6 | CLASS-LEVEL ORDINATION AND PERMANOVA (both scopes) -----------------------

perm_results <- list(); disp_results <- list()
for (s in SCOPES) {
  class_mat <- amr_scope(s) |>
    mutate(class_label = paste(category, subcategory, sep = ": ")) |>
    group_by(sample, class_label) |>
    summarise(hits = sum(value), .groups = "drop") |>
    complete(sample = samples$sample, class_label, fill = list(hits = 0)) |>
    left_join(samples |> select(sample, depth_gb), by = "sample") |>
    mutate(per_gb = hits / depth_gb) |>
    select(sample, class_label, per_gb) |>
    pivot_wider(names_from = class_label, values_from = per_gb, values_fill = 0) |>
    column_to_rownames("sample") |>
    as.matrix()
  class_mat <- class_mat[, colMeans(class_mat > 0) >= 0.05, drop = FALSE]
  profile <- log1p(class_mat)
  pca <- prcomp(profile, center = TRUE, scale. = TRUE)
  save_table(as_tibble(pca$x[, 1:3], rownames = "sample") |> mutate(scope = s),
             paste0("06_class_PCA_scores_", s))
  meta_perm <- samples[match(rownames(profile), samples$sample), ]
  d <- vegdist(profile, method = "bray")
  set.seed(SEED)
  overall <- adonis2(d ~ age + antibiotics + prev_antibiotics + eats_raw + gi,
                     data = meta_perm, permutations = 999, by = NULL)
  set.seed(SEED)
  terms_free <- adonis2(d ~ age + antibiotics + prev_antibiotics + eats_raw + gi,
                        data = meta_perm, permutations = 999, by = "terms")
  ctrl <- how(nperm = 999); setBlocks(ctrl) <- meta_perm$sample_id
  set.seed(SEED)
  margin_blocked <- adonis2(d ~ age + antibiotics + prev_antibiotics + eats_raw + gi,
                            data = meta_perm, permutations = ctrl, by = "margin")
  perm_results[[s]] <- bind_rows(
    as_tibble(overall, rownames = "term") |> mutate(test = "free, overall model"),
    as_tibble(terms_free, rownames = "term") |> mutate(test = "free, sequential terms"),
    as_tibble(margin_blocked, rownames = "term") |> mutate(test = "dog-blocked, marginal terms")) |>
    mutate(scope = s, .before = 1)
  set.seed(SEED)
  bd <- betadisper(d, meta_perm$age)
  pt <- permutest(bd, permutations = 999)
  disp_results[[s]] <- tibble(scope = s, F = pt$tab$F[1], p = pt$tab$`Pr(>F)`[1],
                              wave = names(bd$group.distances),
                              mean_distance_to_centroid = as.numeric(bd$group.distances))
}
save_table(bind_rows(perm_results), "06_PERMANOVA")
save_table(bind_rows(disp_results), "06_dispersion_betadisper")

# 7 | SAMPLE-LEVEL MICROBIOME-ARG CORRELATIONS (both scopes) -------------------

make_micro_mat <- function(rank) {
  bracken |>
    filter(superkingdom == "Bacteria") |>
    mutate(taxon = if (rank == "genus") word(name, 1) else name) |>
    group_by(sample, taxon) |>
    summarise(abund = sum(fraction_total_reads, na.rm = TRUE), .groups = "drop") |>
    pivot_wider(names_from = taxon, values_from = abund, values_fill = 0) |>
    column_to_rownames("sample") |> as.matrix()
}
mg <- make_micro_mat("genus"); ms <- make_micro_mat("species")
alpha <- tibble(sample = rownames(mg), shannon_genus = diversity(mg),
                richness_genus = specnumber(mg)) |>
  full_join(tibble(sample = rownames(ms), shannon_species = diversity(ms),
                   richness_species = specnumber(ms)), by = "sample")
spearman_results <- map_dfr(SCOPES, function(s) {
  dat <- samples_k |>
    left_join(alpha, by = "sample") |>
    mutate(arg_per_gb = .data[[paste0("hits_per_gb_", s)]],
           drug_per_gb = .data[[paste0("drug_hits_", s)]] / depth_gb,
           other_per_gb = .data[[paste0("other_hits_", s)]] / depth_gb,
           arg_richness = .data[[paste0("richness_", s)]])
  vars <- intersect(c("arg_per_gb", "drug_per_gb", "other_per_gb", "arg_richness",
                      "shannon_genus", "richness_genus", "shannon_species",
                      "richness_species", "Bacteria", "Viruses", "Fungi", "Archaea"),
                    names(dat))
  crossing(x = vars, y = vars) |> filter(x < y) |>
    pmap_dfr(function(x, y) {
      d <- dat |> select(all_of(c(x, y))) |> drop_na()
      ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "spearman"))
      tibble(x = x, y = y, n = nrow(d), rho = unname(ct$estimate), p = ct$p.value)
    }) |>
    mutate(q = p.adjust(p, method = "BH"), scope = s, .before = 1)
})
save_table(spearman_results, "07_spearman_ARG_microbiome_tableA8")

# 8 | ASSEMBLY / FLUPOUL-NANOMOTIF CONTEXT --------------------------------------

runs   <- read_tsv(inputs[["runs"]],   show_col_types = FALSE) |> clean_names()
bins   <- read_tsv(inputs[["bins"]],   show_col_types = FALSE) |> clean_names()
motifs <- read_tsv(inputs[["motifs"]], show_col_types = FALSE) |> clean_names()
require_cols(runs, c("barcode", "modification_dataset", "setting", "exit_status",
                     "n_motifs", "bin_motifs_exists"), "NanoMotif manifest")
require_cols(bins, c("barcode", "bin", "checkm_completeness", "checkm_contamination",
                     "checkm_marker_lineage"), "bin master")
require_cols(motifs, c("barcode", "bin", "motif", "mod_type", "mod_position",
                       "setting", "modification_dataset"), "motif master")
save_table(runs |> transmute(barcode, modification_dataset, setting, exit_status,
                             run_success = exit_status == 0 & as_bool(bin_motifs_exists),
                             n_motifs, n_bins_with_motifs), "08_nanomotif_runs")

# Motif records per barcode x regime x modification (Results 4.3.2 counts).
motif_counts <- motifs |>
  count(barcode, setting, modification_dataset, name = "n_motif_records") |>
  pivot_wider(names_from = c(setting, modification_dataset),
              values_from = n_motif_records, values_fill = 0)
save_table(motif_counts, "08_motif_records_by_regime_and_modification")
setting_overlap <- motifs |>
  distinct(barcode, modification_dataset, setting, bin, motif, mod_type, mod_position) |>
  mutate(motif_key = paste(bin, motif, mod_type, mod_position, sep = "|")) |>
  group_by(barcode, modification_dataset, motif_key) |>
  summarise(in_permissive = any(setting == "flupoul_like"),
            in_conservative = any(setting == "conservative"), .groups = "drop") |>
  group_by(barcode, modification_dataset) |>
  summarise(n_permissive = sum(in_permissive), n_conservative = sum(in_conservative),
            n_shared = sum(in_permissive & in_conservative),
            jaccard = n_shared / sum(in_permissive | in_conservative), .groups = "drop")
save_table(setting_overlap, "08_motif_regime_overlap_fig4_4")

bins <- bins |>
  mutate(barcode = as.character(barcode), bin = as.character(bin),
         checkm_completeness = as.numeric(checkm_completeness),
         checkm_contamination = as.numeric(checkm_contamination))
bin_qc <- bins |>
  distinct(barcode, bin, .keep_all = TRUE) |>
  transmute(barcode, bin, marker_lineage = checkm_marker_lineage,
            completeness = checkm_completeness, contamination = checkm_contamination,
            quality = quality_class(completeness, contamination))
save_table(bin_qc, "08_unique_bin_QC")
motif_quality <- motifs |>
  mutate(barcode = as.character(barcode), bin = as.character(bin)) |>
  left_join(bin_qc |> select(barcode, bin, quality), by = c("barcode", "bin")) |>
  count(setting, quality, name = "n_motif_records") |>
  group_by(setting) |>
  mutate(pct = 100 * n_motif_records / sum(n_motif_records)) |>
  ungroup()
save_table(motif_quality, "08_motif_records_by_bin_quality")

abricate <- read_tsv(inputs[["abricate"]], show_col_types = FALSE, name_repair = "minimal")
if (ncol(abricate) != 15L || !identical(names(abricate)[1:2], c("#FILE", "SEQUENCE")))
  stop("Unexpected ABRicate format; expected 15 columns starting '#FILE', 'SEQUENCE'")
names(abricate) <- c("source_file", "contig", "start", "end", "strand", "gene",
                     "coverage", "coverage_map", "gaps", "percent_coverage",
                     "percent_identity", "database", "accession", "product", "resistance")
abricate <- abricate |>
  mutate(barcode = str_extract(source_file, "barcode(?:03|14|36|47)"),
         percent_coverage = as.numeric(percent_coverage),
         percent_identity = as.numeric(percent_identity))
if (any(is.na(abricate$barcode))) stop("ABRicate source paths do not identify all barcodes")
save_table(abricate |> summarise(min_identity = min(percent_identity),
                                 min_coverage = min(percent_coverage),
                                 databases = paste(unique(database), collapse = ";")),
           "08_ABRicate_threshold_check")

contig_bins <- imap_dfr(maps, function(path, barcode_id) {
  read_tsv(path, col_names = c("contig", "original_bin"),
           col_types = cols(.default = col_character())) |>
    mutate(barcode = barcode_id, .before = 1)
})
arg_contigs <- abricate |> distinct(barcode, contig)
assigned <- read_tsv(inputs[["inclusions"]], show_col_types = FALSE) |>
  clean_names() |>
  require_cols(c("barcode", "contig", "bin", "assigned_bin", "confidence"),
               "NanoMotif inclusions") |>
  semi_join(arg_contigs, by = c("barcode", "contig")) |>
  filter(bin == "unbinned", confidence == "high_confidence", nonempty(assigned_bin)) |>
  group_by(barcode, contig) |>
  summarise(n_candidate_bins = n_distinct(assigned_bin),
            assigned_bin = if (n_distinct(assigned_bin) == 1L) first(assigned_bin)
                           else NA_character_, .groups = "drop")
arg_context <- arg_contigs |>
  left_join(contig_bins, by = c("barcode", "contig")) |>
  left_join(assigned, by = c("barcode", "contig")) |>
  left_join(bin_qc |> select(barcode, original_bin = bin, original_quality = quality),
            by = c("barcode", "original_bin")) |>
  left_join(bin_qc |> select(barcode, assigned_bin = bin, assigned_quality = quality),
            by = c("barcode", "assigned_bin")) |>
  mutate(context = case_when(
    nonempty(original_bin) & original_quality %in% c("High quality-like", "Medium quality-like") ~
      "Binned: credible quality",
    nonempty(original_bin) ~ "Binned: low quality or unresolved",
    nonempty(assigned_bin) & assigned_quality %in% c("High quality-like", "Medium quality-like") ~
      "Unbinned: credible candidate bin",
    TRUE ~ "Unbinned: no reliable candidate"))
if (anyDuplicated(arg_context[c("barcode", "contig")]))
  stop("ARG context join duplicated contigs")
save_table(arg_context, "08_ARG_contig_context")
save_table(arg_context |> count(barcode, context, name = "n_contigs"),
           "08_ARG_context_counts_fig4_5")

# 9 | CUPID READ-LEVEL HOST CANDIDATES (barcode36 and barcode03) ---------------

read_cupid <- function(b, what) read_tsv(cupid_inputs[[paste(b, what, sep = "_")]],
                                         show_col_types = FALSE)
plasmid <- map_dfr(CUPID_BARCODES, ~read_cupid(.x, "plasmids") |> mutate(barcode = .x, .before = 1))
top25   <- map_dfr(CUPID_BARCODES, ~read_cupid(.x, "candidates") |> mutate(barcode = .x, .before = 1))
cupid_amr <- map_dfr(CUPID_BARCODES, ~read_cupid(.x, "amr") |> mutate(barcode = .x, .before = 1))
cupid_overlap <- map_dfr(CUPID_BARCODES, ~read_cupid(.x, "overlap") |> mutate(barcode = .x, .before = 1))
require_cols(plasmid, c("read_id", "ARG_symbols", "AMRFinder_methods", "n_methylation_motifs"),
             "CUPID queries")
require_cols(top25, c("query_read_id", "rank", "host_taxid", "rss", "rmsd"), "CUPID top 25")

# Names and ranks resolved against the NCBI taxonomy (October 2026).
tax_lookup <- tribble(
  ~host_taxid, ~taxon, ~taxonomic_rank,
  "165179",  "Segatella copri",              "species",
  "2974251", "Segatella",                    "genus",
  "853",     "Faecalibacterium prausnitzii", "species",
  "437897",  "Megamonas funiformis",         "species",
  "816",     "Bacteroides",                  "genus",
  "135858",  "Catenibacterium",              "genus",
  "2478954", "Catenibacterium sp. co_0103",  "species",
  "2948645", "Burzaovirus (bacteriophage)",  "genus",
  "537007",  "Blautia hansenii DSM 20583",   "strain",
  "186803",  "Lachnospiraceae",              "family",
  "543",     "Enterobacteriaceae",           "family",
  "1843491", "Selenomonadaceae",             "family",
  "171549",  "Bacteroidales",                "order",
  "1239",    "Bacillota",                    "phylum",
  "131567",  "cellular organisms",           "root")
SEGATELLA_TAXIDS <- c("165179", "2974251")

queries <- plasmid |>
  transmute(barcode, query_read_id = read_id, ARG = ARG_symbols,
            method = AMRFinder_methods, n_motifs = n_methylation_motifs,
            query_label = paste0(barcode, " ", ARG, " (", n_motifs, " motifs)"))
candidates <- top25 |>
  mutate(host_taxid = as.character(host_taxid)) |>
  left_join(tax_lookup, by = "host_taxid") |>
  left_join(queries, by = c("barcode", "query_read_id")) |>
  mutate(taxon = coalesce(taxon, paste0("taxid ", host_taxid)),
         taxonomic_rank = coalesce(taxonomic_rank, "unresolved"))
if (any(is.na(candidates$ARG))) stop("CUPID candidates did not map to plasmid queries")
save_table(candidates, "09_CUPID_rank_preserved_candidates")

cupid_scope <- cupid_amr |>
  group_by(barcode) |>
  summarise(n_AMRFinder_hits = n(), n_ARG_positive_reads = n_distinct(`Contig id`),
            .groups = "drop") |>
  left_join(cupid_overlap |> group_by(barcode) |>
              summarise(n_with_4mC_5mC_profile = sum(as_bool(has_4mC_5mC_profile) %in% TRUE),
                        .groups = "drop"), by = "barcode") |>
  left_join(plasmid |> count(barcode, name = "n_plasmid_like_queries"), by = "barcode") |>
  left_join(top25 |> count(barcode, name = "n_candidate_rows"), by = "barcode")
save_table(cupid_scope, "09_CUPID_scope_tableA7")

query_summary <- candidates |>
  group_by(barcode, query_read_id, ARG, n_motifs) |>
  summarise(n_candidates = n(),
            rank1_taxon = taxon[which.min(rank)],
            top_taxon_support = max(table(host_taxid)) / n(),
            top_taxa = paste(names(which(table(taxon) == max(table(taxon)))), collapse = "; "),
            s_copri_species = mean(host_taxid == "165179"),
            segatella_lineage = mean(host_taxid %in% SEGATELLA_TAXIDS),
            rss_range = paste(round(min(rss), 3), round(max(rss), 3), sep = "-"),
            n_distinct_taxa = n_distinct(host_taxid), .groups = "drop")
save_table(query_summary, "09_CUPID_query_summary_fig4_6_4_7")
save_table(candidates |> count(barcode, ARG, taxon, taxonomic_rank, name = "n_top25"),
           "09_CUPID_exact_taxon_support")

# Background check: is Segatella over-represented in each top 25 relative to all
# candidate reads sharing at least one motif with the query? (hypergeometric,
# descriptive; the top 25 are not an independent random sample)
scores_available <- !is.na(cupid_inputs[paste0(CUPID_BARCODES, "_scores")])
if (all(scores_available)) {
  background <- map_dfr(CUPID_BARCODES, ~read_cupid(.x, "scores") |> mutate(barcode = .x, .before = 1)) |>
    mutate(host_taxid = as.character(host_taxid))
  # Top-25 composition is taken from the CUPID top-25 file itself (rank order as
  # written by CUPID), not re-derived, so ties in RSS cannot change membership.
  top_counts <- top25 |>
    mutate(host_taxid = as.character(host_taxid)) |>
    group_by(barcode, query_read_id) |>
    summarise(n_top = n(),
              top_segatella = sum(host_taxid %in% SEGATELLA_TAXIDS),
              top_s_copri = sum(host_taxid == "165179"), .groups = "drop")
  enrichment <- background |>
    group_by(barcode, query_read_id) |>
    summarise(n_background = n(),
              n_sharing_all_query_motifs = sum(n_shared_motifs == n_query_motifs),
              bg_segatella = sum(host_taxid %in% SEGATELLA_TAXIDS),
              bg_s_copri = sum(host_taxid == "165179"), .groups = "drop") |>
    inner_join(top_counts, by = c("barcode", "query_read_id")) |>
    mutate(expected_segatella = n_top * bg_segatella / n_background,
           p_segatella_enrichment = phyper(top_segatella - 1, bg_segatella,
                                           n_background - bg_segatella, n_top,
                                           lower.tail = FALSE),
           expected_s_copri = n_top * bg_s_copri / n_background,
           p_s_copri_enrichment = phyper(top_s_copri - 1, bg_s_copri,
                                         n_background - bg_s_copri, n_top,
                                         lower.tail = FALSE)) |>
    group_by(barcode) |>
    mutate(q_segatella_BH = p.adjust(p_segatella_enrichment, method = "BH"),
           q_s_copri_BH = p.adjust(p_s_copri_enrichment, method = "BH")) |>
    ungroup() |>
    left_join(queries |> select(barcode, query_read_id, ARG, n_motifs),
              by = c("barcode", "query_read_id"))
  save_table(enrichment, "09_CUPID_background_enrichment")
} else {
  warning("plasmid_host_scores.tsv not found for all CUPID barcodes; background check skipped")
}

# 10 | SESSION INFORMATION -----------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# 11 | REPRODUCIBILITY AUDIT -----------------------------------------------------
# Compares reproduced values with those reported in Chapter3_draft1_9.docx.
# `match` is TRUE when the reproduced value rounds to the reported value at the
# reported precision. Rows reported under one scope are checked against both.

dd <- depth_diagnostics
get_dd <- function(col, s) dd[[col]][dd$scope == s]
mfit <- function(s, m, col) model_fit[[col]][model_fit$scope == s & model_fit$model == m]
mrr <- function(s, m, term) model_results$estimate[model_results$scope == s &
                                                     model_results$model == m &
                                                     model_results$term == term]
pt <- function(s, iv, col) paired_tests[[col]][paired_tests$scope == s & paired_tests$interval == iv]
perm <- function(s, test, term, col) {
  x <- bind_rows(perm_results)
  x[[col]][x$scope == s & x$test == test & x$term == term]
}
sp <- function(s, x, y, col) spearman_results[[col]][spearman_results$scope == s &
                                                       spearman_results$x == x & spearman_results$y == y]
val1 <- function(x) if (length(x)) as.numeric(x[1]) else NA_real_
both <- function(metric, reported, digits, f) {
  tibble(metric = metric, reported = reported,
         reproduced_all = round(val1(f("all")), digits),
         reproduced_non_snp = round(val1(f("non_snp")), digits))
}
one <- function(metric, reported, digits, value) {
  tibble(metric = metric, reported = reported, reproduced_all = round(val1(value), digits),
         reproduced_non_snp = NA_real_)
}
cat_hits <- function(cat) category_summary$total_hits[category_summary$category == cat]
gene_n <- function(g) gene_table$n_present[gene_table$gene_id == g]
gene_hits <- function(g) gene_table$total_hits[gene_table$gene_id == g]
comp_n <- function(comp) sum(gene_table$compartment == comp)

audit <- bind_rows(
  one("samples", 191, 0, cohort_summary$n_samples),
  one("dogs", 73, 0, cohort_summary$n_dogs),
  one("median depth (Gb)", 2.16, 2, cohort_summary$median_depth_gb),
  one("IQR depth (Gb)", 2.03, 2, cohort_summary$iqr_depth_gb),
  one("raw-fed samples", 21, 0, cohort_summary$raw_fed),
  both("Table A1 wave 1 median hits/Gb (draft: non-SNP 529.0; all 7,766.8)", 529.0, 1,
       function(s) median(samples[[paste0("hits_per_gb_", s)]][samples$wave == 1])),
  one("antibiotic samples", 36, 0, cohort_summary$antibiotics),
  one("GI samples", 104, 0, cohort_summary$gi),
  both("r raw hits vs reads", 0.83, 2, function(s) get_dd("r_raw_vs_reads", s)),
  both("r raw hits vs bases", 0.84, 2, function(s) get_dd("r_raw_vs_bases", s)),
  both("r hits per million reads vs reads", 0.05, 2, function(s) get_dd("r_per_mreads_vs_reads", s)),
  both("r hits per Gb vs bases", 0.13, 2, function(s) get_dd("r_per_gb_vs_bases", s)),
  both("r hits per Gb vs bacterial fraction", 0.08, 2, function(s) get_dd("r_per_gb_vs_bacterial_fraction", s)),
  both("Michaelis-Menten Rmax", 167.73, 2, function(s) get_dd("mm_Rmax", s)),
  both("depth for 90% Rmax (Gb)", 5.95, 2, function(s) get_dd("mm_depth_gb_90", s)),
  both("depth for 95% Rmax (Gb)", 12.56, 2, function(s) get_dd("mm_depth_gb_95", s)),
  both("Good's coverage median", 0.998, 3, function(s) get_dd("goods_median", s)),
  one("distinct ARGs detected", 927, 0, nrow(gene_table)),
  one("SNP-dependent genes", 165, 0, sum(gene_table$snp_dependent)),
  one("non-SNP-dependent genes", 762, 0, sum(!gene_table$snp_dependent)),
  one("total hits", 3338527, 0, sum(gene_table$total_hits)),
  one("SNP-dependent hits", 3045971, 0, sum(gene_table$total_hits[gene_table$snp_dependent])),
  one("Drugs hits", 3286078, 0, cat_hits("Drugs")),
  one("Metals hits", 23935, 0, cat_hits("Metals")),
  one("Multi-compound hits", 19494, 0, cat_hits("Multi-compound")),
  one("Biocides hits", 9020, 0, cat_hits("Biocides")),
  one("Drugs median hits/Gb", 6529.2, 1,
      category_summary$median_hits_per_gb[category_summary$category == "Drugs"]),
  one("core genes", 25, 0, comp_n("Core (>=80%)")),
  one("intermediate genes", 223, 0, comp_n("Intermediate (20-<80%)")),
  one("rare genes", 679, 0, comp_n("Rare (<20%)")),
  one("MEG_3988 samples present", 190, 0, gene_n("MEG_3988")),
  one("MEG_3988 hits", 660361, 0, gene_hits("MEG_3988")),
  one("MEG_3992 hits", 276392, 0, gene_hits("MEG_3992")),
  one("MEG_7 samples present", 186, 0, gene_n("MEG_7")),
  one("MEG_4223 prevalence (%)", 72.8, 1,
      target_summary$prevalence_pct[target_summary$target == "MEG_4223 / O23S"]),
  one("MEG_4223 hits", 447749, 0,
      target_summary$total_hits[target_summary$target == "MEG_4223 / O23S"]),
  one("MEG_4223 median hits/Gb where detected", 1069, 0,
      target_summary$median_hits_per_gb_where_detected[target_summary$target == "MEG_4223 / O23S"]),
  one("OXA hits", 7451, 0,
      target_summary$total_hits[target_summary$target == "OXA / Class D beta-lactamases"]),
  one("OXA class median hits/Gb where detected (draft: 7.6)", 7.6, 1,
      target_summary$median_hits_per_gb_where_detected[target_summary$target == "OXA / Class D beta-lactamases"]),
  one("MEG_4741 median hits/Gb where detected", 7.6, 1,
      target_summary$median_hits_per_gb_where_detected[target_summary$target == "MEG_4741 (OXA) alone"]),
  both("paired 4 to 7 months median delta (all markers)", -0.254, 3,
       function(s) pt(s, "4 to 7 months", "median_delta_log")),
  both("paired 4 to 12 months median delta (all markers)", -0.226, 3,
       function(s) pt(s, "4 to 12 months", "median_delta_log")),
  both("paired 4 to 12 months BH p (all markers)", 0.0011, 4,
       function(s) pt(s, "4 to 12 months", "p_BH")),
  both("primary model AIC", 3045.5, 1, function(s) mfit(s, "total_burden", "AIC")),
  both("primary model BIC", 3074.8, 1, function(s) mfit(s, "total_burden", "BIC")),
  both("batch model AIC", 3032.1, 1, function(s) mfit(s, "total_burden_batch", "AIC")),
  both("batch model BIC", 3103.7, 1, function(s) mfit(s, "total_burden_batch", "BIC")),
  both("primary RR 7 months", 0.82, 2, function(s) mrr(s, "total_burden", "age7 months")),
  both("primary RR 12 months", 0.86, 2, function(s) mrr(s, "total_burden", "age12 months")),
  both("richness model depth RR", 1.88, 2, function(s) mrr(s, "richness", "log10_total_bases")),
  both("PERMANOVA free overall R2", 0.041, 3,
       function(s) perm(s, "free, overall model", "Model", "R2")),
  both("PERMANOVA free overall p", 0.155, 3,
       function(s) perm(s, "free, overall model", "Model", "Pr(>F)")),
  both("PERMANOVA blocked wave R2", 0.029, 3,
       function(s) perm(s, "dog-blocked, marginal terms", "age", "R2")),
  both("PERMANOVA blocked wave p", 0.015, 3,
       function(s) perm(s, "dog-blocked, marginal terms", "age", "Pr(>F)")),
  both("betadisper p", 0.064, 3,
       function(s) unique(disp_results[[s]]$p)),
  both("Spearman ARG richness vs genus Shannon rho", 0.226, 3,
       function(s) sp(s, "arg_richness", "shannon_genus", "rho")),
  both("Spearman ARG richness vs genus Shannon q", 0.0057, 4,
       function(s) sp(s, "arg_richness", "shannon_genus", "q")),
  one("CUPID barcode36 lnu(C) top-25 Segatella lineage", 21, 0,
      if (exists("enrichment")) enrichment$top_segatella[enrichment$barcode == "barcode36" &
                                                          grepl("lnu", enrichment$ARG)] else NA),
  one("CUPID barcode36 aph top-25 Segatella lineage", 22, 0,
      if (exists("enrichment")) enrichment$top_segatella[enrichment$barcode == "barcode36" &
                                                          grepl("aph", enrichment$ARG)] else NA)
) |>
  mutate(match_all = abs(reproduced_all - reported) < 1e-9,
         match_non_snp = abs(reproduced_non_snp - reported) < 1e-9)
save_table(audit, "00_reproducibility_audit")
print(audit, n = Inf, width = Inf)
message("Chapter 4 analysis finished. Outputs: ", OUT_DIR,
        "\nCheck tables/00_reproducibility_audit.csv first.")
