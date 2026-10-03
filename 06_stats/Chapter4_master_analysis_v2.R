# Chapter 4 | Puppy faecal resistome and ARG host context
# One fresh R session; one output root. Run from the top with source().
#
# Results questions:
# 1. Does within-dog ARG burden change with age or reported antibiotics?
# 2. Which assembled ARG contigs have credible MAG context (FLUPOUL/NanoMotif)?
# 3. How do host candidates differ between the high- and low-information
#    putative plasmid ARG reads (CUPID, barcode36 and barcode03)?
#
# AMR++ counts, ABRicate assemblies, and CUPID read hits are different assays.
# Never combine their row counts or call a methylation match a proven host.
# Existing inputs are not rewritten. A fresh run replaces files in the output root.
# Required packages: dplyr, tidyr, readr, readxl, stringr, ggplot2, purrr,
# tibble, janitor, glmmTMB, broom.mixed.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(janitor)
  library(glmmTMB)
  library(broom.mixed)
})

# 1 | CONFIGURATION AND SHARED QC ----------------------------------------------

DATA_DIR <- path.expand("~/Documents/2110_dog_AMR/exploration")
OUT_DIR <- file.path(DATA_DIR, "Chapter3/Final_analysis/Chapter3_analysis")
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDS_DIR <- file.path(OUT_DIR, "rds")
invisible(lapply(c(OUT_DIR, FIG_DIR, TAB_DIR, RDS_DIR),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

MIN_READS <- 50000L  # current Chapter 3 analysis threshold; verify Chapter 2 separately
BARCODES <- c("barcode03", "barcode14", "barcode36", "barcode47")
CUPID_BARCODES <- c("barcode36", "barcode03")
METHYL_DIR <- file.path(DATA_DIR, "Chapter3/Chapter3_analysis/FLUPOUL/FLUPOUL/methylation_R_master")
CUPID_DIR <- file.path(DATA_DIR, "Chapter3/Chapter3_analysis/CUPID")
CUPID_HOST_DIR <- file.path(CUPID_DIR)

inputs <- c(
  amr = file.path(DATA_DIR, "AMR_matrix_group.csv"),
  tech = file.path(DATA_DIR, "NanoStats_matrix_final.csv"),
  metadata = file.path(DATA_DIR, "litter_metadata.xlsx"),
  runs = file.path(METHYL_DIR, "methylation_run_manifest.tsv"),
  bins = file.path(METHYL_DIR, "methylation_bin_master.tsv"),
  motifs = file.path(METHYL_DIR, "methylation_motif_master_long.tsv"),
  inclusions = file.path(METHYL_DIR, "methylation_include_contigs_master.tsv"),
  contamination = file.path(METHYL_DIR, "methylation_contamination_master.tsv"),
  abricate = file.path(DATA_DIR, "abricate_megares_master.tsv")
)
# The same five read-level files are required for each CUPID barcode. Keeping
# this manifest in one place makes the low-depth replicate an actual input,
# rather than a number copied into a figure annotation.
cupid_paths <- function(barcode) c(
  amr = file.path(CUPID_DIR, paste0(barcode, ".amrfinder.tsv")),
  overlap = file.path(CUPID_DIR, paste0(barcode, ".ARG_read_methylation_master.tsv")),
  plasmids = file.path(CUPID_HOST_DIR, paste0(barcode, ".plasmid_ARG_reads.tsv")),
  candidates = file.path(CUPID_HOST_DIR, paste0(barcode, ".plasmid_host_top25.tsv")),
  consensus = file.path(CUPID_HOST_DIR, paste0(barcode, ".host_taxon_consensus.tsv"))
)
cupid_inputs <- unlist(lapply(CUPID_BARCODES, cupid_paths))
names(cupid_inputs) <- paste(rep(CUPID_BARCODES, each = 5L),
                             rep(names(cupid_paths(CUPID_BARCODES[1])),
                                 times = length(CUPID_BARCODES)), sep = "_")
inputs <- c(inputs, cupid_inputs)
maps <- setNames(file.path(DATA_DIR, "contig_bins",
                            paste0(BARCODES, ".contig_bin.tsv")), BARCODES)
missing <- c(inputs, maps)[!file.exists(c(inputs, maps))]
if (length(missing)) stop("Missing Chapter 3 input(s):\n", paste(missing, collapse = "\n"))

save_table <- function(x, stem) readr::write_csv(x, file.path(TAB_DIR, paste0(stem, ".csv")))
save_figure <- function(x, stem, width = 9, height = 6) {
  ggplot2::ggsave(file.path(FIG_DIR, paste0(stem, ".tiff")), x,
                  width = width, height = height, dpi = 600,
                  device = "tiff", compression = "lzw", bg = "white")
  ggplot2::ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), x,
                  width = width, height = height, bg = "white")
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
            completeness >= 90 & contamination <= 5 ~ "High quality-like",
            completeness >= 50 & contamination <= 10 ~ "Medium quality-like",
            TRUE ~ "Low/problematic")
}
safe_test <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || all(x == 0)) return(NA_real_)
  tryCatch(stats::wilcox.test(x, mu = 0, exact = FALSE)$p.value,
           error = function(e) NA_real_)
}
save_table(tibble(input = names(c(inputs, maps)), path = unname(c(inputs, maps)),
                  size_bytes = as.numeric(file.info(c(inputs, maps))$size)), "00_input_manifest")

# 2 | COHORT: BURDEN, PERSISTENCE, AND INFERENCE -------------------------------

tech_all <- read_csv(inputs[["tech"]], show_col_types = FALSE) |>
  clean_names() |>
  require_cols(c("sample", "sample_id", "wave", "num_reads", "total_bases"), "NanoStats") |>
  mutate(across(c(sample, sample_id), as.character))
tech_all <- tech_all |>
  mutate(wave = as.integer(wave), num_reads = as.numeric(num_reads),
         total_bases = as.numeric(total_bases))
if (anyDuplicated(tech_all$sample)) stop("NanoStats has duplicate sample IDs")
sample_qc <- tech_all |>
  transmute(sample, sample_id, wave, num_reads, total_bases,
            included = !is.na(num_reads) & num_reads >= MIN_READS &
                       !is.na(total_bases) & total_bases > 0)
save_table(sample_qc, "01_cohort_sample_inclusion")
tech <- tech_all |> semi_join(filter(sample_qc, included), by = "sample")

amr_raw <- read_csv(inputs[["amr"]], show_col_types = FALSE) |>
  clean_names() |>
  require_cols(c("sample", "sample_id", "wave", "gene_id", "gene_name", "value"),
               "AMR matrix") |>
  # read_csv() may infer purely numeric dog IDs as doubles. Keep all dog IDs
  # as characters so the gene-by-dog join agrees with the NanoStats cohort.
  mutate(sample = as.character(sample), sample_id = as.character(sample_id),
         wave = as.integer(wave), value = suppressWarnings(as.numeric(value)))
if (any(is.na(amr_raw$value)))
  stop("AMR matrix contains missing or nonnumeric values; inspect the input before modelling")
amr_raw <- amr_raw |>
  mutate(requires_snp_confirmation = str_detect(coalesce(gene_name, ""),
                                                 "RequiresSNPConfirmation"))
save_table(amr_raw |> semi_join(tech, by = "sample") |>
             group_by(requires_snp_confirmation) |>
             summarise(n_gene_ids = n_distinct(gene_id), hits = sum(value),
                       .groups = "drop"), "01_SNP_marker_compartments")
# The primary burden is the non-SNP-requiring compartment. This is a marker
# classification, not proof that every remaining hit is acquired or expressed.
amr <- amr_raw |>
  filter(!requires_snp_confirmation) |>
  semi_join(tech, by = "sample")
all_marker_sample <- amr_raw |>
  semi_join(tech, by = "sample") |>
  group_by(sample) |>
  summarise(all_marker_hits = sum(value), .groups = "drop")

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

amr_sample <- amr |>
  group_by(sample) |>
  summarise(arg_hits = sum(value), arg_richness = n_distinct(gene_id[value > 0]),
            .groups = "drop")
samples <- tech |>
  left_join(meta, by = "sample") |>
  left_join(amr_sample, by = "sample") |>
  left_join(all_marker_sample, by = "sample") |>
  mutate(arg_hits = replace_na(arg_hits, 0),
         all_marker_hits = replace_na(all_marker_hits, 0),
         arg_richness = replace_na(arg_richness, 0L),
         arg_hits_per_gb = arg_hits * 1e9 / total_bases,
         all_marker_hits_per_gb = all_marker_hits * 1e9 / total_bases,
         age = factor(wave, levels = 1:3, labels = c("4 months", "7 months", "12 months")))
save_table(samples, "02_cohort_sample_data")
saveRDS(samples, file.path(RDS_DIR, "cohort_samples.rds"))

gene_prevalence <- amr |>
  filter(value > 0) |>
  distinct(sample, gene_id, gene_name) |>
  count(gene_id, gene_name, name = "samples_positive") |>
  mutate(n_samples = nrow(samples), prevalence = samples_positive / n_samples,
         group = case_when(prevalence >= .8 ~ "Core (>=80%)",
                           prevalence >= .2 ~ "Intermediate (20–<80%)",
                           TRUE ~ "Rare (<20%)")) |>
  arrange(desc(samples_positive))
save_table(gene_prevalence, "02_gene_prevalence")

# Repeat detection is evaluated only among dogs observed at all three ages.
# The core refers to >=80% prevalence across samples; it does not imply that
# every core gene occurs in every dog at every wave.
three_wave_dogs <- samples |>
  count(sample_id, name = "n_waves") |>
  filter(n_waves == 3L) |>
  pull(sample_id)
core_gene_ids <- gene_prevalence |>
  filter(group == "Core (>=80%)") |>
  pull(gene_id)
core_dog_detection <- tidyr::expand_grid(sample_id = three_wave_dogs,
                                         gene_id = core_gene_ids) |>
  left_join(amr |>
              filter(value > 0, sample_id %in% three_wave_dogs,
                     gene_id %in% core_gene_ids) |>
              distinct(sample_id, gene_id, wave) |>
              count(sample_id, gene_id, name = "waves_detected"),
            by = c("sample_id", "gene_id")) |>
  mutate(waves_detected = replace_na(waves_detected, 0L))
save_table(core_dog_detection, "02_core_gene_three_wave_detections")
save_table(core_dog_detection |>
             summarise(n_dogs = length(three_wave_dogs),
                       n_core_genes = length(core_gene_ids),
                       n_gene_dog_pairs = n(),
                       n_detected_all_three = sum(waves_detected == 3L),
                       fraction_detected_all_three = mean(waves_detected == 3L)),
           "02_core_gene_persistence_summary")

# Pair by dog and explicit waves. Direct 4–12-month pairs are retained; a dog
# missing the 7-month sample is never mislabelled as a 4–7-month interval.
make_pairs <- function(a, b, label, rate_col, marker_set_label) {
  left <- samples |>
    filter(wave == a) |>
    transmute(sample_id, rate_a = .data[[rate_col]], abx_a = antibiotics)
  right <- samples |>
    filter(wave == b) |>
    transmute(sample_id, rate_b = .data[[rate_col]], abx_b = antibiotics)
  inner_join(left, right, by = "sample_id") |>
    mutate(interval = label, marker_set = marker_set_label,
           delta_log = log1p(rate_b) - log1p(rate_a),
           # NA exposure stays unknown, rather than being treated as FALSE.
           abx_either = case_when(abx_a %in% TRUE | abx_b %in% TRUE ~ TRUE,
                                  abx_a %in% FALSE & abx_b %in% FALSE ~ FALSE,
                                  TRUE ~ NA))
}
pairs <- map_dfr(list(
  c("Non-SNP-requiring markers", "arg_hits_per_gb"),
  c("All AMR++ markers", "all_marker_hits_per_gb")
), function(definition) {
  bind_rows(make_pairs(1, 2, "4 → 7 months", definition[2], definition[1]),
            make_pairs(2, 3, "7 → 12 months", definition[2], definition[1]),
            make_pairs(1, 3, "4 → 12 months", definition[2], definition[1]))
}) |>
  mutate(interval = factor(interval,
                           levels = c("4 → 7 months", "7 → 12 months", "4 → 12 months")),
         marker_set = factor(marker_set,
                             levels = c("All AMR++ markers", "Non-SNP-requiring markers")))
save_table(pairs, "02_within_dog_pairs")
paired_tests <- pairs |>
  group_by(marker_set, interval) |>
  summarise(n_pairs = n(), median_delta_log = median(delta_log),
            p_change = safe_test(delta_log), .groups = "drop") |>
  group_by(marker_set) |>
  mutate(p_BH = p.adjust(p_change, method = "BH")) |>
  ungroup()
save_table(paired_tests, "02_within_dog_tests")

# Exploratory antibiotic comparison of paired changes, with sample sizes.
abx_pair_test <- pairs |>
  filter(marker_set == "Non-SNP-requiring markers", !is.na(abx_either)) |>
  group_by(interval) |>
  summarise(n_exposed = sum(abx_either), n_unexposed = sum(!abx_either),
            p_abx = if (n_distinct(abx_either) < 2) NA_real_ else
              tryCatch(wilcox.test(delta_log ~ abx_either, exact = FALSE)$p.value,
                       error = function(e) NA_real_), .groups = "drop")
save_table(abx_pair_test, "02_antibiotic_pair_comparison")

model_data <- samples |>
  filter(complete.cases(antibiotics, prev_antibiotics, eats_raw, gi,
                        sample_id, age, total_bases)) |>
  mutate(across(c(antibiotics, prev_antibiotics, eats_raw, gi),
                ~factor(.x, levels = c(FALSE, TRUE))))
if (n_distinct(model_data$sample_id) < 3L) stop("Insufficient dogs for cohort model")
fit_model <- function(formula, data, label) {
  tryCatch(glmmTMB(formula, family = nbinom2(), data = data),
           error = function(e) {warning(label, ": ", conditionMessage(e)); NULL})
}
primary_model <- fit_model(
  arg_hits ~ age + antibiotics + prev_antibiotics + eats_raw + gi +
    offset(log(total_bases / 1e9)) + (1 | sample_id), model_data, "Primary model")
# Batch and age can be strongly confounded. This is sensitivity analysis,
# not evidence that either biological age or batch is causally responsible.
batch_model <- NULL
if ("batch" %in% names(model_data) && n_distinct(model_data$batch) > 1L) {
  batch_model <- fit_model(
    arg_hits ~ age + antibiotics + prev_antibiotics + eats_raw + gi +
      factor(batch) + offset(log(total_bases / 1e9)) + (1 | sample_id),
    model_data, "Batch sensitivity model")
}
tidy_rate <- function(model, label) {
  if (is.null(model)) return(tibble(model = label, term = NA_character_,
                                     estimate = NA_real_, conf.low = NA_real_,
                                     conf.high = NA_real_, p.value = NA_real_))
  tryCatch(broom.mixed::tidy(model, effects = "fixed", component = "cond",
                             conf.int = TRUE, exponentiate = TRUE) |>
             mutate(model = label, .before = 1),
           error = function(e) {
             warning("Could not summarise ", label, ": ", conditionMessage(e))
             tibble(model = label, term = NA_character_, estimate = NA_real_,
                    conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_)
           })
}
save_table(bind_rows(tidy_rate(primary_model, "primary"),
                     tidy_rate(batch_model, "batch_sensitivity")), "02_model_rate_ratios")
save_table(tibble(model = c("primary", "batch_sensitivity"),
                  n = nrow(model_data),
                  AIC = c(if (is.null(primary_model)) NA_real_ else AIC(primary_model),
                          if (is.null(batch_model)) NA_real_ else AIC(batch_model))),
           "02_model_status")
saveRDS(list(primary = primary_model, batch_sensitivity = batch_model),
        file.path(RDS_DIR, "cohort_models.rds"))

p_definition_sensitivity <- ggplot(pairs,
                 aes(interval, delta_log, colour = factor(abx_either,
                                                          levels = c(FALSE, TRUE)))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_jitter(width = .12, height = 0, alpha = .52, size = 1.5) +
  stat_summary(aes(group = interval), fun = median, geom = "point",
               shape = 95, size = 11, colour = "black") +
  facet_wrap(~ marker_set, ncol = 2) +
  scale_colour_manual(values = c("FALSE" = "#376C91", "TRUE" = "#C56E30"),
                      labels = c("No", "Yes"), na.value = "grey65") +
  labs(x = NULL, y = "Within-dog change in log(1 + ARG hits/Gb)",
       colour = "Antibiotics reported\nat either wave",
       title = "A  |  The age signal depends on the marker definition",
       subtitle = "Each point is a dog; black bars mark medians; SNP-requiring markers need separate interpretation") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
save_figure(p_definition_sensitivity, "Fig_S1_marker_definition_sensitivity", 12, 5)

p_long <- pairs |>
  filter(marker_set == "Non-SNP-requiring markers") |>
  ggplot(aes(interval, delta_log, colour = factor(abx_either,
                                                  levels = c(FALSE, TRUE)))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_jitter(width = .13, height = 0, alpha = .48, size = 1.5) +
  stat_summary(aes(group = interval), fun = median, geom = "point",
               shape = 95, size = 11, colour = "black") +
  scale_colour_manual(values = c("FALSE" = "#376C91", "TRUE" = "#C56E30"),
                      labels = c("No", "Yes"), na.value = "grey65") +
  labs(x = NULL, y = "Within-dog change in log(1 + ARG hits/Gb)",
       colour = "Antibiotics at\neither wave") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
save_figure(p_long, "Fig_A_longitudinal_ARG_change", 10, 5)

# A single supplementary line plot replaces the many individual ARG panels.
# These three names are exact non-SNP MEGARes gene names; MEG_4223 is a
# SNP-requiring 23S marker and must never be counted as an OXA beta-lactamase.
marker_names <- c("OXA", "TEM", "CTX")
marker_presence <- amr |>
  filter(value > 0, gene_name %in% marker_names) |>
  distinct(sample, gene_name) |>
  left_join(samples |> select(sample, age), by = "sample") |>
  count(age, gene_name, name = "n_positive") |>
  complete(age = factor(levels(samples$age), levels = levels(samples$age)),
           gene_name = marker_names,
           fill = list(n_positive = 0L)) |>
  left_join(samples |> count(age, name = "n_samples"), by = "age") |>
  mutate(prevalence = n_positive / n_samples,
         gene_name = factor(gene_name, levels = marker_names))
save_table(marker_presence, "02_targeted_marker_prevalence_by_age")
p_markers <- ggplot(marker_presence,
                    aes(age, prevalence, group = gene_name, colour = gene_name)) +
  geom_line(linewidth = .85) + geom_point(size = 2.3) +
  scale_colour_manual(values = c(OXA = "#286B80", TEM = "#D0783F", CTX = "#7755A5")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = NULL, y = "Samples with a detected marker", colour = "Gene family",
       title = "Surveillance marker prevalence across sampling ages",
       caption = "Cross-sectional proportions at each age; these lines do not imply persistence in the same dog.") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())
save_figure(p_markers, "Fig_S2_selected_marker_prevalence", 8, 4.5)

# 3 | FLUPOUL/NANOMOTIF: ASSEMBLY AND MAG CONTEXT ------------------------------

runs <- read_tsv(inputs[["runs"]], show_col_types = FALSE) |> clean_names()
bins <- read_tsv(inputs[["bins"]], show_col_types = FALSE) |> clean_names()
motifs <- read_tsv(inputs[["motifs"]], show_col_types = FALSE) |> clean_names()
require_cols(runs, c("barcode", "modification_dataset", "setting", "exit_status",
                     "n_motifs", "bin_motifs_exists"), "NanoMotif manifest")
require_cols(bins, c("barcode", "bin", "checkm_completeness",
                     "checkm_contamination", "checkm_marker_lineage"), "bin master")
require_cols(motifs, c("barcode", "bin", "motif", "mod_type", "mod_position"), "motif master")

run_summary <- runs |>
  transmute(barcode, modification_dataset, setting, exit_status,
            run_success = exit_status == 0 & as_bool(bin_motifs_exists),
            n_motifs, n_bins_with_motifs)
save_table(run_summary, "03_nanomotif_runs")

# A compact threshold-sensitivity audit, without treating each setting as an
# independent biological observation. Identities include their original bin.
setting_overlap <- motifs |>
  distinct(barcode, modification_dataset, setting, bin, motif,
           mod_type, mod_position) |>
  mutate(motif_key = paste(bin, motif, mod_type, mod_position, sep = "|")) |>
  group_by(barcode, modification_dataset, motif_key) |>
  summarise(in_flupoul = any(setting == "flupoul_like"),
            in_conservative = any(setting == "conservative"), .groups = "drop") |>
  group_by(barcode, modification_dataset) |>
  summarise(n_flupoul = sum(in_flupoul),
            n_conservative = sum(in_conservative),
            n_shared = sum(in_flupoul & in_conservative),
            jaccard = n_shared / sum(in_flupoul | in_conservative),
            .groups = "drop")
save_table(setting_overlap, "03_motif_setting_overlap")

bins <- bins |>
  mutate(barcode = as.character(barcode), bin = as.character(bin),
         checkm_completeness = as.numeric(checkm_completeness),
         checkm_contamination = as.numeric(checkm_contamination))
qc_disagreement <- bins |>
  group_by(barcode, bin) |>
  summarise(n_completeness = n_distinct(checkm_completeness, na.rm = TRUE),
            n_contamination = n_distinct(checkm_contamination, na.rm = TRUE),
            .groups = "drop") |>
  filter(n_completeness > 1 | n_contamination > 1)
save_table(qc_disagreement, "03_bin_QC_disagreements")
if (nrow(qc_disagreement)) stop("CheckM values disagree between runs for the same bin")
bin_qc <- bins |>
  distinct(barcode, bin, .keep_all = TRUE) |>
  transmute(barcode, bin, marker_lineage = checkm_marker_lineage,
            completeness = checkm_completeness,
            contamination = checkm_contamination,
            quality = quality_class(completeness, contamination))
save_table(bin_qc, "03_unique_bin_QC")

motif_counts <- motifs |>
  distinct(barcode, bin, motif, mod_type, mod_position) |>
  count(barcode, bin, name = "n_distinct_motif_calls")
# A motif call is evidence in the associated bin, not independent host proof.
bin_qc <- left_join(bin_qc, motif_counts, by = c("barcode", "bin")) |>
  mutate(n_distinct_motif_calls = replace_na(n_distinct_motif_calls, 0L))

# ABRicate's first header is literally '#FILE'. Rename source columns by
# validated position BEFORE janitor::clean_names(); otherwise `file` may resolve
# to base::file(), causing the previous script's str_extract() error.
abricate <- read_tsv(inputs[["abricate"]], show_col_types = FALSE,
                     name_repair = "minimal")
if (ncol(abricate) != 15L || !identical(names(abricate)[1:2], c("#FILE", "SEQUENCE")))
  stop("Unexpected ABRicate format; expected 15 columns starting '#FILE', 'SEQUENCE'")
names(abricate) <- c("source_file", "contig", "start", "end", "strand", "gene",
                     "coverage", "coverage_map", "gaps", "percent_coverage",
                     "percent_identity", "database", "accession", "product", "resistance")
abricate <- abricate |>
  mutate(barcode = str_extract(source_file, "barcode(?:03|14|36|47)"),
         percent_coverage = as.numeric(percent_coverage),
         percent_identity = as.numeric(percent_identity),
         start = as.numeric(start), end = as.numeric(end))
if (any(is.na(abricate$barcode))) stop("ABRicate source paths do not identify all barcodes")
if (any(abricate$database != "megares", na.rm = TRUE))
  warning("ABRicate includes a database other than MEGARes; inspect the master")

contig_bins <- imap_dfr(maps, function(path, barcode_id) {
  read_tsv(path, col_names = c("contig", "original_bin"),
           col_types = cols(.default = col_character())) |>
    mutate(barcode = barcode_id, .before = 1)
})
if (anyDuplicated(contig_bins[c("barcode", "contig")]))
  stop("An assembly contig maps to multiple original MetaBAT2 bins")

arg_contigs <- abricate |> distinct(barcode, contig)
# The inclusion master contains many rows with NA assignments. Filter to the
# ARG contigs before collapsing and never call max() on all-missing values.
include_arg <- read_tsv(inputs[["inclusions"]], show_col_types = FALSE) |>
  clean_names() |>
  require_cols(c("barcode", "contig", "bin", "assigned_bin", "confidence"),
               "NanoMotif inclusions") |>
  semi_join(arg_contigs, by = c("barcode", "contig")) |>
  filter(bin == "unbinned", confidence == "high_confidence", nonempty(assigned_bin)) |>
  distinct(barcode, contig, assigned_bin, modification_dataset, setting)
assigned <- include_arg |>
  group_by(barcode, contig) |>
  summarise(n_candidate_bins = n_distinct(assigned_bin),
            candidate_bins = paste(sort(unique(assigned_bin)), collapse = ";"),
            assigned_bin = if (n_distinct(assigned_bin) == 1L) first(assigned_bin)
                           else NA_character_, .groups = "drop")
save_table(assigned, "03_NanoMotif_high_confidence_ARG_candidates")

contam_arg <- read_tsv(inputs[["contamination"]], show_col_types = FALSE) |>
  clean_names() |>
  distinct(barcode, contig) |>
  mutate(methylation_inconsistency_flag = TRUE)

# One row per barcode-specific assembly contig; hit-level details remain below.
arg_context <- arg_contigs |>
  left_join(contig_bins, by = c("barcode", "contig")) |>
  left_join(assigned, by = c("barcode", "contig")) |>
  left_join(bin_qc |>
              rename(original_bin = bin, original_quality = quality,
                     original_lineage = marker_lineage,
                     original_completeness = completeness,
                     original_contamination = contamination,
                     original_motifs = n_distinct_motif_calls),
            by = c("barcode", "original_bin")) |>
  left_join(bin_qc |>
              rename(assigned_bin = bin, assigned_quality = quality,
                     assigned_lineage = marker_lineage,
                     assigned_completeness = completeness,
                     assigned_contamination = contamination,
                     assigned_motifs = n_distinct_motif_calls),
            by = c("barcode", "assigned_bin")) |>
  left_join(contam_arg, by = c("barcode", "contig")) |>
  mutate(methylation_inconsistency_flag = replace_na(methylation_inconsistency_flag, FALSE),
         context = case_when(
           nonempty(original_bin) & original_quality == "High quality-like" ~ "Binned: high quality-like",
           nonempty(original_bin) & original_quality == "Medium quality-like" ~ "Binned: medium quality-like",
           nonempty(original_bin) & original_quality == "Low/problematic" ~ "Binned: low/problematic",
           nonempty(original_bin) ~ "Binned: QC unavailable",
           !is.na(n_candidate_bins) & n_candidate_bins > 1L ~ "Unbinned: conflicting assignments",
           nonempty(assigned_bin) & assigned_quality == "High quality-like" ~ "Unbinned: candidate high quality-like",
           nonempty(assigned_bin) & assigned_quality == "Medium quality-like" ~ "Unbinned: candidate medium quality-like",
           nonempty(assigned_bin) & assigned_quality == "Low/problematic" ~ "Unbinned: candidate low/problematic",
           nonempty(assigned_bin) ~ "Unbinned: candidate QC unavailable",
           TRUE ~ "Unbinned: unresolved"))
if (nrow(arg_context) != nrow(arg_contigs) ||
    anyDuplicated(arg_context[c("barcode", "contig")]))
  stop("ARG context join changed the number of unique contigs")
save_table(arg_context, "03_ARG_contig_context")
save_table(abricate |> left_join(arg_context, by = c("barcode", "contig")),
           "03_ABRicate_hits_with_context")
saveRDS(arg_context, file.path(RDS_DIR, "ARG_contig_context.rds"))
save_table(arg_context |>
             count(barcode, context, name = "n_contigs"), "03_ARG_context_counts")

context_levels <- c("Binned: high quality-like", "Binned: medium quality-like",
                    "Binned: low/problematic", "Binned: QC unavailable",
                    "Unbinned: candidate high quality-like",
                    "Unbinned: candidate medium quality-like",
                    "Unbinned: candidate low/problematic",
                    "Unbinned: candidate QC unavailable",
                    "Unbinned: conflicting assignments", "Unbinned: unresolved")
context_colours <- c("#237C85", "#76B4B0", "#D99D63", "#B0B0B0",
                     "#2A8061", "#81B79E", "#C97C5C", "#CACACA",
                     "#78539A", "#E1E1E1")
p_assembly <- arg_context |>
  mutate(barcode = factor(barcode, levels = BARCODES),
         context = factor(context, levels = context_levels)) |>
  count(barcode, context, name = "n_contigs", .drop = FALSE) |>
  ggplot(aes(barcode, n_contigs, fill = context)) +
  geom_col(width = .75, colour = "white", linewidth = .15) +
  scale_fill_manual(values = setNames(context_colours, context_levels), drop = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, .08))) +
  labs(x = NULL, y = "Unique ARG-bearing contigs",
       fill = "Genomic context", title = "B  |  Assembly links ARGs to bins with uneven QC",
       subtitle = "Original bin membership precedes candidate methylation assignments") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 2))
save_figure(p_assembly, "Fig_S3_detailed_assembly_ARG_context", 10, 6.5)

# Four categories answer the main Results question without a 10-colour legend.
# A bin assignment is a genomic context candidate, not a validated ARG host.
arg_context_main <- arg_context |>
  mutate(context_main = case_when(
    !nonempty(original_bin) ~ "Not assembled into a bin",
    original_quality %in% c("High quality-like", "Medium quality-like") ~
      "Binned: high/medium quality-like",
    original_quality == "Low/problematic" ~ "Binned: low/problematic",
    TRUE ~ "Binned: CheckM QC unavailable")) |>
  count(context_main, name = "n_contigs") |>
  mutate(context_main = factor(context_main,
    levels = c("Binned: high/medium quality-like", "Binned: low/problematic",
               "Binned: CheckM QC unavailable", "Not assembled into a bin")))
save_table(arg_context_main, "03_ARG_context_main_figure")
p_assembly_main <- ggplot(arg_context_main,
                           aes(n_contigs, context_main, fill = context_main)) +
  geom_col(width = .67, show.legend = FALSE) +
  geom_text(aes(label = n_contigs), hjust = -.25, size = 3.5) +
  scale_fill_manual(values = c("Binned: high/medium quality-like" = "#287D78",
                               "Binned: low/problematic" = "#D19A61",
                               "Binned: CheckM QC unavailable" = "#AAAAAA",
                               "Not assembled into a bin" = "#735FA0")) +
  scale_x_continuous(limits = c(0, max(arg_context_main$n_contigs) * 1.18),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(limits = rev(levels(arg_context_main$context_main))) +
  labs(x = "ARG-bearing contigs", y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
save_figure(p_assembly_main, "Fig_B_ARG_context_summary", 8, 5)

# 4 | CUPID: READ-LEVEL CANDIDATES AND THEIR LIMITS ----------------------------

cupid <- setNames(lapply(CUPID_BARCODES, function(barcode) {
  paths <- cupid_paths(barcode)
  lapply(paths, function(path) read_tsv(path, show_col_types = FALSE))
}), CUPID_BARCODES)
cupid_hits <- imap_dfr(cupid, ~mutate(.x$amr, barcode = .y, .before = 1))
cupid_overlap <- imap_dfr(cupid, ~mutate(.x$overlap, barcode = .y, .before = 1))
plasmid <- imap_dfr(cupid, ~mutate(.x$plasmids, barcode = .y, .before = 1))
top25 <- imap_dfr(cupid, ~mutate(.x$candidates, barcode = .y, .before = 1))
cupid_consensus <- imap_dfr(cupid, ~mutate(.x$consensus, barcode = .y, .before = 1))
require_cols(cupid_hits, c("Contig id", "Element symbol", "Method"), "CUPID AMRFinder")
require_cols(cupid_overlap, c("read_id", "has_4mC_5mC_profile"),
             "CUPID methylation overlap")
require_cols(plasmid, c("read_id", "ARG_symbols", "AMRFinder_methods",
                        "n_methylation_motifs"), "CUPID queries")
require_cols(top25, c("query_read_id", "host_taxid", "rank", "rss", "rmsd"),
             "CUPID top-25 candidates")
require_cols(cupid_consensus, c("query_read_id", "host_taxid", "n_top25_reads"),
             "CUPID host consensus")

query_audit <- top25 |>
  group_by(barcode, query_read_id) |>
  summarise(n = n(), n_ranks = n_distinct(rank),
            first_rank = min(rank), last_rank = max(rank), .groups = "drop")
if (any(query_audit$n != 25L | query_audit$n_ranks != 25L |
        query_audit$first_rank != 1L | query_audit$last_rank != 25L))
  stop("CUPID queries do not each have exactly 25 unique candidate ranks")
consensus_audit <- cupid_consensus |>
  group_by(barcode, query_read_id) |>
  summarise(n = sum(n_top25_reads), .groups = "drop")
if (any(consensus_audit$n != 25L) || nrow(consensus_audit) != nrow(query_audit))
  stop("CUPID consensus counts do not sum to 25 for each query")
save_table(query_audit, "04_CUPID_top25_validation")

# These are the taxid/rank names supplied in the final host consensus master.
# Unknown taxids remain unresolved rather than being silently named a species.
tax_lookup <- tribble(
  ~host_taxid, ~taxon, ~taxonomic_rank,
  "165179", "Segatella copri", "species",
  "2974251", "Segatella", "genus",
  "853", "Faecalibacterium prausnitzii", "species",
  "543", "Enterobacteriaceae", "family",
  "131567", "cellular organisms", "root",
  "135858", "Catenibacterium", "genus",
  "1843491", "Selenomonadaceae", "family",
  "171549", "Bacteroidales", "order",
  "437897", "Megamonas funiformis", "species",
  "2948645", "taxid 2948645", "unresolved",
  "2478954", "taxid 2478954", "unresolved",
  "816", "Bacteroides", "genus",
  "186803", "Lachnospiraceae", "family",
  "1239", "Bacillota", "phylum",
  "537007", "Blautia hansenii DSM 20583", "no rank")
queries <- plasmid |>
  transmute(barcode, query_read_id = read_id, ARG = ARG_symbols,
            method = AMRFinder_methods, n_motifs = n_methylation_motifs,
            query_label = paste0(barcode, "  ", ARG, "  (", n_motifs, " motifs)"))
if (anyDuplicated(queries[c("barcode", "query_read_id")]))
  stop("Duplicate CUPID query IDs")
if (nrow(queries) != nrow(query_audit) ||
    nrow(anti_join(query_audit, queries, by = c("barcode", "query_read_id"))))
  stop("CUPID queries and candidate lists do not agree")
candidates <- top25 |>
  mutate(host_taxid = as.character(host_taxid)) |>
  left_join(tax_lookup, by = "host_taxid") |>
  left_join(queries, by = c("barcode", "query_read_id")) |>
  mutate(taxon = coalesce(taxon, paste0("taxid ", host_taxid)),
         taxonomic_rank = coalesce(taxonomic_rank, "unresolved"),
         taxon_rank_label = paste0(taxon, " [", taxonomic_rank, "; ", host_taxid, "]"))
if (any(is.na(candidates$ARG))) stop("CUPID candidates did not map to plasmid queries")
save_table(candidates, "04_CUPID_rank_preserved_candidates")
candidate_support <- candidates |>
  count(barcode, ARG, query_read_id, query_label, host_taxid, taxon,
        taxonomic_rank, taxon_rank_label,
        name = "n_top25") |>
  group_by(barcode, query_read_id) |>
  mutate(fraction_top25 = n_top25 / sum(n_top25)) |>
  ungroup()
save_table(candidate_support, "04_CUPID_exact_taxid_support")

score_separation <- candidates |>
  arrange(barcode, query_read_id, rank) |>
  group_by(barcode, query_read_id, ARG, n_motifs) |>
  summarise(n_candidates = n(), n_taxids = n_distinct(host_taxid),
            best_rss = first(rss), second_rss = nth(rss, 2),
            rss_gap = best_rss - second_rss,
            best_rmsd = first(rmsd), second_rmsd = nth(rmsd, 2),
            .groups = "drop")
save_table(score_separation, "04_CUPID_score_separation")
save_table(cupid_hits |> count(barcode, `Element symbol`, Method, name = "n_hits"),
           "04_CUPID_AMRFinder_methods")
cupid_scope <- imap_dfr(cupid, function(one, barcode) {
  if (!setequal(one$amr$`Contig id`, one$overlap$read_id))
    stop("AMRFinder and methylation-overlap read IDs disagree for ", barcode)
  tibble(barcode,
         n_AMRFinder_hits = nrow(one$amr),
         n_ARG_positive_reads = n_distinct(one$amr$`Contig id`),
         n_with_methylation_profile = sum(as_bool(one$overlap$has_4mC_5mC_profile) %in% TRUE),
         n_putative_plasmid_queries = nrow(one$plasmids),
         n_candidate_rows = nrow(one$candidates))
})
save_table(cupid_scope, "04_CUPID_scope")

candidate_categories <- candidates |>
  mutate(category = case_when(
    taxon == "Segatella copri" & taxonomic_rank == "species" ~ "S. copri (species)",
    taxon == "Segatella" & taxonomic_rank == "genus" ~ "Segatella (genus)",
    taxonomic_rank %in% c("species", "genus") ~ "Other species/genus",
    TRUE ~ "Higher rank/unresolved")) |>
  count(barcode, query_read_id, query_label, ARG, n_motifs, category,
        name = "n_top25") |>
  mutate(category = factor(category,
    levels = c("S. copri (species)", "Segatella (genus)",
               "Other species/genus", "Higher rank/unresolved")),
    query_label = factor(query_label, levels = rev(queries$query_label)))
category_audit <- candidate_categories |>
  group_by(barcode, query_read_id) |>
  summarise(n = sum(n_top25), .groups = "drop")
if (any(category_audit$n != 25L))
  stop("CUPID plot groups must each contain 25 candidates")
save_table(candidate_categories, "04_CUPID_main_figure_categories")
lineage_support <- candidate_categories |>
  group_by(barcode, query_read_id, ARG, n_motifs) |>
  summarise(n_segatella_species = sum(n_top25[category == "S. copri (species)"]),
            n_segatella_genus = sum(n_top25[category == "Segatella (genus)"]),
            n_top25 = sum(n_top25), .groups = "drop") |>
  mutate(n_lineage = n_segatella_species + n_segatella_genus,
         fraction_lineage = n_lineage / n_top25)
save_table(lineage_support, "04_CUPID_lineage_support")

p_cupid <- ggplot(candidate_categories,
                  aes(x = n_top25, y = query_label, fill = category)) +
  geom_col(width = .67, colour = "white", linewidth = .3) +
  geom_text(aes(label = if_else(n_top25 >= 4L, as.character(n_top25), "")),
            position = position_stack(vjust = .5), size = 3.1) +
  scale_fill_manual(values = c("S. copri (species)" = "#246E91",
                               "Segatella (genus)" = "#6BA9BE",
                               "Other species/genus" = "#D49A60",
                               "Higher rank/unresolved" = "#A5A5A5")) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25),
                     limits = c(0, 25), expand = expansion(mult = c(0, 0))) +
  labs(x = "Methylation-matched chromosome reads (top 25)", y = NULL,
       fill = "Kraken assignment") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(nrow = 2))
save_figure(p_cupid, "Fig_C_CUPID_candidate_taxa", 10, 5)

# Figures A and C are independent main-chapter candidates. The compact
# FLUPOUL bin-context chart B is optional supplementary material; the small
# context table is saved above if it reads more cleanly than a separate plot.

# Reproducible headline audit. No biological interpretation is hard-coded.
headline <- tibble(
  metric = c("cohort_samples", "cohort_dogs", "cohort_read_threshold",
             "ABRicate_hits", "ARG_contigs", "ARG_contigs_originally_binned",
             "ARG_contigs_originally_unbinned", "unbinned_ARG_contigs_high_confidence_unique_bin",
             "CUPID_putative_plasmid_queries", "CUPID_candidate_reads"),
  value = c(nrow(samples), n_distinct(samples$sample_id), MIN_READS,
            nrow(abricate), nrow(arg_context), sum(nonempty(arg_context$original_bin)),
            sum(!nonempty(arg_context$original_bin)),
            sum(!nonempty(arg_context$original_bin) & nonempty(arg_context$assigned_bin)),
            nrow(plasmid), nrow(top25)))
save_table(headline, "00_headline_counts")
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
message("Chapter 3 finished. Outputs: ", OUT_DIR)
