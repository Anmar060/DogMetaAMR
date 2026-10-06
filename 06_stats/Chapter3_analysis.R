##############################################################################
## Chapter 2 – Puppy gut microbiome: diversity, composition & GI upset
## Master analysis script (post-50k-read filter)
## Consolidates: postfilter_PCoA_PCA_GLMM_GI_AMR.R + microbiome_depth.R
##                + dog_feed_functional.R
## Last updated: May 2026
##############################################################################

# ── 0. PACKAGES ──────────────────────────────────────────────────────────────
# Run once to install:
# install.packages(c(
#   "tidyverse","janitor","readr","stringr","lubridate","skimr",
#   "broom","broom.mixed","performance",
#   "vegan","lme4","glmmTMB","compositions","zCompositions",
#   "ggrepel","pals","scales","colorspace","viridis","Polychrome",
#   "mgcv","emmeans","lmerTest"
# ))
# BiocManager::install(c("ANCOMBC","ALDEx2","Maaslin2","Nonpareil"))

library(tidyverse); library(janitor); library(skimr)
library(vegan);     library(glmmTMB); library(lme4);   library(lmerTest)
library(compositions); library(zCompositions)
library(ggplot2);   library(ggrepel)
library(pals);      library(colorspace); library(Polychrome)
library(broom);     library(broom.mixed); library(performance); library(emmeans)
library(ANCOMBC);   library(ALDEx2);  library(Maaslin2)
library(Nonpareil); library(mgcv); library(coin)

# Colour palette (25 distinct colours for litter plots)
c25 <- c(
  "dodgerblue2","#E31A1C","green4","#6A3D9A","#FF7F00","black","gold1",
  "skyblue2","#FB9A99","palegreen2","#CAB2D6","#FDBF6F","gray70","khaki2",
  "maroon","orchid1","deeppink1","blue1","steelblue4","darkturquoise",
  "green1","yellow4","yellow3","darkorange4","brown"
)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────
tech    <- read_csv("~/Documents/2110_dog_AMR/exploration/NanoStats_matrix_final.csv")      |> clean_names()
bracken <- read_csv("~/Documents/2110_dog_AMR/exploration/merged_with_superkingdom_wave.csv") |> clean_names()
amr     <- read_csv("~/Documents/2110_dog_AMR/exploration/AMR_matrix_group.csv")            |> clean_names()
meta    <- read_csv("~/Documents/2110_dog_AMR/exploration/litter_metadata.csv")       |> clean_names()
#diet-level file (~/Documents/2110_dog_AMR/dog_feed/dog_feed/dogfeed_diet.csv) and
#feed-level file (~/Documents/2110_dog_AMR/dog_feed/dog_feed/dogfeed_feed.csv) loaded in Section 6 only
diet    <- read_csv("~/Documents/2110_dog_AMR/dog_feed/dogfeed_diet.csv")       |> clean_names()
feed    <- read_csv("~/Documents/2110_dog_AMR/dog_feed/dogfeed_feed.csv")       |> clean_names()

# Fix known typo
amr <- amr |> mutate(category = recode(category, "Multi-compount" = "Multi-compound"))

# Harmonise types
for (df_name in c("tech","bracken","amr","meta")) {
  df <- get(df_name)
  df <- df |>
    mutate(sample = as.character(sample),
           wave      = as.double(wave))
  if ("litter_id" %in% names(df))
    df <- df |> mutate(litter_id = as.character(litter_id))
  assign(df_name, df)
}

# ── 2. DEPTH FILTER (≥50,000 reads) ──────────────────────────────────────────
keep_samples <- tech |> filter(num_reads >= 50000) |> pull(sample)

tech    <- tech    |> filter(sample %in% keep_samples)
bracken <- bracken |> filter(sample %in% keep_samples)
amr     <- amr     |> filter(sample %in% keep_samples)
meta    <- meta    |> filter(sample %in% keep_samples)

message(sprintf(
  "Retained %d samples (W1=%d, W2=%d, W3=%d)",
  length(keep_samples),
  sum(tech$wave==1), sum(tech$wave==2), sum(tech$wave==3)
))
# Expected: 191 samples (71 / 66 / 54)

# ── 3. BUILD MAIN SAMPLE TABLE ────────────────────────────────────────────────
samp <- meta |>
  transmute(
    sample, sample_id, wave,
    litter_id = coalesce(litter_id),
    gi            = as.logical(q9_any_gi),
    diarrhoea     = as.logical(q9a_diarrhoea),
    vomiting      = as.factor(q9b_vomiting),
    eats_raw      = as.logical(eats_raw),
    antibiotics   = as.factor(q6_antibiotics),
    scavenges     = as.logical(q4_any_scavenges),
    farm_contact, wild_contact, multi_abs,
    prev_antibiotics = as.factor(q6_antibiotics_oral_previous),
    abs_days      = coalesce(q6d_antibiotics_duration_days),
    days2_abs     = coalesce(q6_antibiotics_duration_before_sample),
    abs_oral      = as.logical(q6_antibiotics_oral),
    abs_topic     = as.logical(q6_antibiotics_topical),
    abs_injec     = as.logical(q6_antibiotics_injection)
  )

# Join sequencing tech stats
samp <- samp |>
  left_join(
    tech |> transmute(sample, sample_id, wave, batch, mean_read_length,
                      mean_read_quality, num_reads, total_bases),
    by = c("sample","sample_id","wave")
  )

# Summarise AMR per sample and join
amr_summ <- amr |>
  group_by(sample) |>
  summarise(
    amr_total      = sum(value, na.rm = TRUE),
    amr_richness   = n_distinct(gene_id[value > 0]),
    amr_drug_rel   = sum(value[group == "Drug_related"], na.rm = TRUE),
    .groups = "drop"
  )

samp <- samp |>
  left_join(amr_summ, by = "sample") |>
  mutate(
    amr_per_million_reads = (amr_total / num_reads) * 1e6,
    amr_per_gb            = (amr_total / total_bases) * 1e9,
    log_amr               = log1p(amr_total),
    log_totalbases        = log10(total_bases + 1),
    wave_f                = factor(wave)
  )

# Kingdom-level fractions per sample
kingdom_frac <- bracken |>
  group_by(sample, superkingdom) |>
  summarise(frac = sum(fraction_total_reads, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = superkingdom, values_from = frac, values_fill = 0)

samp <- samp |> left_join(kingdom_frac, by = "sample")

# ── 4. SEQUENCING DEPTH & NONPAREIL ──────────────────────────────────────────
###-------- 4a. AMR saturation curve (Michaelis-Menten) ------
nls_fit <- nls(
  amr_richness ~ (Rmax * total_bases) / (K + total_bases),
  data  = samp,
  start = list(Rmax = max(samp$amr_richness, na.rm = TRUE) * 1.5,
               K    = median(samp$total_bases, na.rm = TRUE))
)
# Rmax ≈ 168 unique AMR genes; 90% saturation at ~6 Gb, 95% at ~12.6 Gb
coef_nls  <- coef(nls_fit)
depth_90  <- (0.9  * coef_nls["K"]) / (1 - 0.9)  / 1e9   # ~5.95 Gb
depth_95  <- (0.95 * coef_nls["K"]) / (1 - 0.95) / 1e9   # ~12.56 Gb

###--------- 4b. Good's coverage per sample ----------
goods_coverage <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (!length(x)) return(NA_real_)
  1 - (sum(x == 1) / sum(x))
}

amr_wide <- amr |>
  group_by(sample, gene_id) |>
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = gene_id, values_from = value, values_fill = 0)

abund_list <- setNames(
  lapply(seq_len(nrow(amr_wide)), function(i) as.numeric(amr_wide[i, -1])),
  amr_wide$sample
)

cov_df <- tibble(
  sample       = names(abund_list),
  amr_coverage = vapply(abund_list, goods_coverage, numeric(1))
)
# Median Good's coverage ≈ 0.998; 4 samples below 0.90

samp <- samp |> left_join(cov_df, by = "sample")

###------ 4c. Nonpareil (run separately per wave; results summarised below) -----
# Source: microbiome_depth.R — run on Eddie HPC with Nonpareil .npo files
# Key results embedded as comments for record:
#   Wave 1: median coverage C ≈ higher, median Nd ≈ lower
#   Wave 3: median C ≈ lower, median Nd ≈ higher
#   Kruskal-Wallis: C and Nd differ significantly across all waves;
#                   LRstar differs between wave 1 and later waves
#   Wave 3 has ~45% lower total_bases than waves 1–2 (Kruskal-Wallis + pairwise Wilcox BH)
#   Conclusion: Wave 3 lower coverage partly technical (lower depth) partly biological
#               (genuine increase in community diversity across puppy development)

# ── 5. MICROBIOME COMPOSITION ─────────────────────────────────────────────────
###------ 5a. Build genus-level matrix (Bacteria only) --------
brack_bac <- bracken |>
  filter(superkingdom == "Bacteria") |>
  mutate(genus = word(name, 1)) |>
  group_by(sample, genus) |>
  summarise(abund = sum(fraction_total_reads, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = genus, values_from = abund, values_fill = 0)

micro_mat <- brack_bac |>
  column_to_rownames("sample") |>
  as.matrix()

###------ 5b. Alpha diversity (Shannon, Richness) on full micro_mat ---------
shannon  <- vegan::diversity(micro_mat, index = "shannon")
richness <- vegan::specnumber(micro_mat)

samp <- samp |>
  mutate(shannon  = shannon[sample],
         richness = richness[sample])

###------ 5c. CLR transformation (for PCA / PERMANOVA) -----------
# Prevalence filter: keep genera in ≥10% of samples
prev_filt <- colMeans(micro_mat > 0)
micro_filt <- micro_mat[, prev_filt >= 0.10]

# Sample filter: drop samples where >80% genera are zero
zero_prop  <- rowMeans(micro_filt == 0)
micro_filt <- micro_filt[zero_prop <= 0.80, ]
# One sample removed (617829_03); 190 samples × 44 genera remain

micro_no0 <- zCompositions::cmultRepl(
  micro_filt, label = 0, method = "CZM", z.delete = FALSE, z.warning = 1
)
micro_clr <- compositions::clr(micro_no0)

# Filtered metadata aligned to CLR matrix
samp_filt <- samp |> filter(sample %in% rownames(micro_clr))

###------ 5d. PCoA – Bray-Curtis (exploratory) ----------
bc_dist <- vegdist(micro_mat, method = "bray")
bc_ord  <- cmdscale(bc_dist, k = 2)
bc_df   <- as.data.frame(bc_ord) |>
  rownames_to_column("sample") |>
  left_join(samp, by = "sample")

###------ 5e. PCoA – Aitchison / CLR (inferential) ----------
aitch_dist <- dist(micro_clr, method = "euclidean")
aitch_ord  <- cmdscale(aitch_dist, k = 2, eig = TRUE)

ord_df <- as.data.frame(aitch_ord$points) |>
  rownames_to_column("sample") |>
  rename(PCo1 = V1, PCo2 = V2) |>
  left_join(samp, by = "sample")

###------ 5f. PCA on CLR (to extract PCs for GLMM) ----------
micro_pca <- prcomp(micro_clr, center = TRUE, scale. = FALSE)
# Variance explained: PC1≈11%, PC2≈9%, PC3≈8%, PC4≈7% ...
summary(micro_pca)$importance[, 1:6]

micro_scores <- micro_pca$x |>
  as.data.frame() |>
  rownames_to_column("sample") |>
  dplyr::select(sample, PC1, PC2, PC3)

samp_clr <- samp_filt |>
  left_join(micro_scores, by = "sample") |>
  mutate(
    batch_str  = as.character(batch),
    n_plates   = stringr::str_count(batch_str, "_") + 1
  )

###------ 5g. PERMANOVA (adonis2) -----------
meta_for_perm <- samp |>
  filter(sample %in% rownames(micro_clr)) |>
  arrange(match(sample, rownames(micro_clr)))

# Standard (free permutation)
adonis2(
  aitch_dist ~ gi + eats_raw + antibiotics + prev_antibiotics + scavenges + batch + wave,
  data = meta_for_perm, permutations = 999
)
# Model R² = 0.128, F = 1.31, p = 0.001

# Marginal (by = "margin")
adonis2(
  aitch_dist ~ gi + eats_raw + antibiotics + prev_antibiotics + scavenges + batch + wave,
  data = meta_for_perm, permutations = 999, by = "margin"
)
# eats_raw: R²=0.012, p=0.001; batch: R²=0.079, p=0.052; all others ns

# Block-permutation (within dog, accounting for repeated measures)
perm_block <- how(nperm = 999)
setBlocks(perm_block) <- meta_for_perm$sample_id

adonis2(
  aitch_dist ~ gi + eats_raw + antibiotics + prev_antibiotics + scavenges + batch + wave,
  data = meta_for_perm, permutations = perm_block
)
# Model R² = 0.128, p = 0.002; marginal: only batch p = 0.034

# Dispersion tests (betadisper)
disp_gi  <- betadisper(aitch_dist, meta_for_perm$gi)
anova(disp_gi); permutest(disp_gi)   # ns — no heterogeneity bias

disp_abx <- betadisper(aitch_dist, meta_for_perm$antibiotics)
anova(disp_abx)                       # ns

# ── 6. ALPHA DIVERSITY BY DIET: Alpha diversity GLMMs (diet and antibiotic associations)  ────────────────────────────────────────────────
# Source: dog_feed_functional.R
# Requires: dogfeed_feed.csv (particulate per-sample feed data)
# feed <- read_csv("dogfeed_feed.csv") |> clean_names()
## Requires: alpha_all (built below), feed_samp, bracken, tech, meta

### --------- 6a — Build abundance matrices at genus and species level ------------
make_micro_mat <- function(bracken, rank = c("genus","species")) {
  rank <- match.arg(rank)
  bracken %>%
    filter(superkingdom == "Bacteria") %>%
    mutate(taxon = case_when(
      rank == "genus"   ~ word(name, 1),
      rank == "species" ~ name
    )) %>%
    group_by(sample, taxon) %>%
    summarise(abund = sum(fraction_total_reads, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = taxon, values_from = abund, values_fill = 0) %>%
    column_to_rownames("sample") %>%
    as.matrix()
}

feed_micro_mat_genus   <- make_micro_mat(bracken, "genus")
feed_micro_mat_species <- make_micro_mat(bracken, "species")

###------------- 6b — Calculate alpha diversity -----------------
alpha_from_mat <- function(mat) {
  tibble(
    sample     = rownames(mat),
    shannon    = vegan::diversity(mat, index = "shannon"),
    richness   = vegan::specnumber(mat),
    simpson    = vegan::diversity(mat, index = "simpson"),
    invsimpson = vegan::diversity(mat, index = "invsimpson")
  )
}

alpha_genus   <- alpha_from_mat(feed_micro_mat_genus)   |> mutate(rank = "genus")
alpha_species <- alpha_from_mat(feed_micro_mat_species) |> mutate(rank = "species")
alpha_all <- alpha_bound |>
  left_join(
    samp |> dplyr::select(-any_of(c("shannon", "richness", "simpson", "invsimpson"))),
    by = "sample"
  )


### ---------- 6c — Final GLMMs (batch as random, sample_id as random) -------------
## Shannon ~ raw feeding (primary model)
m_shannon_genus <- glmmTMB(
  shannon ~ eats_raw + wave + antibiotics + prev_antibiotics + scavenges +
    (1 | batch) + (1 | sample_id),
  data   = alpha_all %>% filter(rank == "genus"),
  family = gaussian()
)

m_shannon_species <- glmmTMB(
  shannon ~ eats_raw + wave + antibiotics + prev_antibiotics + scavenges +
    (1 | batch) + (1 | sample_id),
  data   = alpha_all %>% filter(rank == "species"),
  family = gaussian()
)

## Richness ~ raw feeding + sequencing depth
m_richness_genus <- glmmTMB(
  richness ~ eats_raw + wave + antibiotics + prev_antibiotics + scavenges +
    log10(total_bases + 1) + (1 | sample_id),
  data   = alpha_all %>% filter(rank == "genus"),
  family = nbinom2()  # richness is a count 
)
summary(m_richness_genus)

m_richness_species <- glmmTMB(
  richness ~ eats_raw + wave + antibiotics + prev_antibiotics + scavenges +
    log10(total_bases + 1) + (1 | batch) + (1 | sample_id),
  data   = alpha_all %>% filter(rank == "species"),
  family = nbinom2()
)
summary(m_richness_species)

###-- 6d — Extract and save results -----------
coef_shannon <- bind_rows(
  broom.mixed::tidy(m_shannon_genus,   effects = "fixed") |> mutate(rank = "genus"),
  broom.mixed::tidy(m_shannon_species, effects = "fixed") |> mutate(rank = "species")
)

coef_richness <- bind_rows(
  broom.mixed::tidy(m_richness_genus,   effects = "fixed") |> mutate(rank = "genus"),
  broom.mixed::tidy(m_richness_species, effects = "fixed") |> mutate(rank = "species")
)

ci_shannon <- bind_rows(
  confint(m_shannon_genus,   parm = "eats_rawTRUE") |> as_tibble() |> mutate(rank = "genus"),
  confint(m_shannon_species, parm = "eats_rawTRUE") |> as_tibble() |> mutate(rank = "species")
) |> rename(ci_low = `2.5 %`, ci_high = `97.5 %`)

ci_richness <- bind_rows(
  confint(m_richness_genus,   parm = "eats_rawTRUE") |> as_tibble() |> mutate(rank = "genus"),
  confint(m_richness_species, parm = "eats_rawTRUE") |> as_tibble() |> mutate(rank = "species")
) |> rename(ci_low = `2.5 %`, ci_high = `97.5 %`)

glmm_results <- bind_rows(
  coef_shannon |> left_join(ci_shannon, by = "rank") |> mutate(outcome = "shannon"),
  coef_richness |> left_join(ci_richness, by = "rank") |> mutate(outcome = "richness")
)

write_csv(glmm_results, "alpha_diversity_GLMM_results.csv")

###------ 6e — Emmeans for raw feeding (Shannon genus) -------
emm_raw <- emmeans::emmeans(m_shannon_genus, ~ eats_raw)
write_csv(as_tibble(emm_raw), "alpha_diversity_emmeans_raw.csv")

###--- 6f — Model fit -----
performance::r2(m_shannon_genus)
AIC(m_shannon_genus, m_shannon_species)
AIC(m_richness_genus, m_richness_species)

broom.mixed::tidy(m_shannon_genus, effects = "fixed")
broom.mixed::tidy(m_shannon_species, effects = "fixed")

confint(m_shannon_species, parm = "eats_rawTRUE")

# Key finding (embedded for reference):
#   GLMMs (glmmTMB, batch as random, dog as random):
#   Shannon ~ eats_raw + wave + antibiotics + prev_antibiotics + scavenges
#     eats_rawTRUE: β=0.32, 95%CI [0.08, 0.56], genus level
#     eats_rawTRUE: β=0.38, 95%CI [0.12, 0.64], species level
#   Richness ~ eats_raw + ... + log10(total_bases):
#     genus level barely significant (p≈0.049); species ns
#   Conclusion: raw feeding associated with higher Shannon diversity;
#               richness driven primarily by sequencing depth

# emmeans (genus-level Shannon):
#   eats_raw FALSE: 1.20 (95%CI 0.96–1.43)
#   eats_raw TRUE:  1.52 (95%CI 1.21–1.83)

# ── 7. DIFFERENTIAL ABUNDANCE ─────────────────────────────────────────────────
###------ 7a. Helper functions -----
clean_tax <- function(x) x |> stringr::str_squish() |> make.names()
direction_lfc <- function(x, pos_lbl, neg_lbl) if_else(x > 0, pos_lbl, neg_lbl)

###------- 7b. Species-level count matrix (Bacteria, ≥10% prevalence, ≥50 total reads) -----
tax_counts_long <- bracken |>
  filter(superkingdom == "Bacteria") |>
  group_by(sample, species = name) |>
  summarise(count = sum(new_est_reads, na.rm = TRUE), .groups = "drop")

count_mat <- tax_counts_long |>
  pivot_wider(names_from = species, values_from = count, values_fill = 0) |>
  column_to_rownames("sample") |>
  as.matrix()

prev_sp  <- colMeans(count_mat > 0)
count_filt <- count_mat[, prev_sp >= 0.10 & colSums(count_mat) >= 50]
# 191 → 63 species retained

feat_df <- as.data.frame(count_filt)

meta_df2 <- samp |>
  filter(sample %in% rownames(feat_df)) |>
  distinct(sample, .keep_all = TRUE) |>
  arrange(match(sample, rownames(feat_df))) |>
  column_to_rownames("sample")

# Relative abundance (for MaAsLin2)
relab  <- sweep(count_filt, 1, rowSums(count_filt), "/")
relab[is.na(relab)] <- 0
relab2 <- relab; colnames(relab2) <- clean_tax(colnames(relab2))

meta_maas <- meta_df2 |>
  rownames_to_column("sample") |>
  mutate(
    log_totalbases   = log10(total_bases + 1),
    gi               = as.factor(gi),
    vomiting         = as.factor(vomiting),
    diarrhoea        = as.factor(diarrhoea),
    antibiotics      = as.factor(antibiotics),
    eats_raw         = as.factor(eats_raw),
    wave             = as.factor(wave),
    sample_id        = as.factor(sample_id)
  ) |>
  column_to_rownames("sample")
meta_maas <- meta_maas[rownames(relab2), , drop = FALSE]

###--------------- 7c. GI upset — ANCOM-BC2 ------------
ancom_gi <- ancombc2(
  data        = feat_df,
  meta_data   = meta_df2,
  fix_formula = "gi + antibiotics + eats_raw + wave + batch",
  taxa_are_rows = FALSE, p_adj_method = "BH",
  prv_cut = 0.10, lib_cut = 0, group = "gi",
  struc_zero = TRUE, neg_lb = TRUE
)
gi_ancom <- ancom_gi$res |>
  transmute(
    taxon, ancom_lfc = lfc_giTRUE, ancom_q = q_giTRUE,
    ancom_diff = diff_giTRUE,
    ancom_sig = q_giTRUE < 0.05 & diff_giTRUE
  )
# 27 taxa significant; 18 lower, 9 higher in GI upset (see Results §2.6.2)

###------------ 7d. GI upset — MaAsLin2 ------------
fit_gi <- Maaslin2(
  input_data     = as.data.frame(relab2),
  input_metadata = as.data.frame(meta_maas),
  output         = "maaslin2_gi",
  fixed_effects  = c("gi","antibiotics","eats_raw","wave","log_totalbases"),
  random_effects = c("sample_id"),
  normalization  = "NONE", transform = "LOG",
  analysis_method= "LM",   correction = "BH", max_significance = 0.05
)
gi_maas <- read_tsv("maaslin2_gi/all_results.tsv") |>
  filter(metadata == "gi") |>
  mutate(taxon = clean_tax(feature), maas_sig = qval < 0.05)
# 1 significant taxon: Faecalibacterium prausnitzii (coef = -0.92)

###----------- 7e. GI upset — ALDEx2 (two-group, no covariates; supplementary only) -------
ald_counts <- t(as.matrix(count_filt)); storage.mode(ald_counts) <- "numeric"
meta_al    <- meta_df2[colnames(ald_counts), , drop = FALSE]
set.seed(1)
ald <- ALDEx2::aldex.clr(ald_counts, ifelse(meta_al$gi, "GI","NoGI"),
                          mc.samples = 128, denom = "all", verbose = FALSE)
ald_res <- tibble(
  species = rownames(ALDEx2::aldex.ttest(ald)),
  p = ALDEx2::aldex.ttest(ald)$we.ep,
  q = p.adjust(ALDEx2::aldex.ttest(ald)$we.ep, "BH"),
  effect = ALDEx2::aldex.effect(ald)$effect
) |> arrange(q)
# Very conservative; no taxon passes q < 0.10
write_csv(ald_res, "aldex2_gi_results.csv")

###-------- 7f. Merge three-method comparison ------------
gi_comparison <- full_join(
  gi_ancom |> mutate(taxon = clean_tax(taxon)),
  gi_maas  |> dplyr::select(taxon, maas_coef = coef, maas_q = qval, maas_sig),
  by = "taxon"
) |>
  mutate(n_methods_sig = rowSums(across(c(ancom_sig, maas_sig), ~replace_na(.x, FALSE)))) |>
  arrange(desc(n_methods_sig), ancom_q, maas_q)
write_csv(gi_comparison, "GI_taxon_comparison_ANCOM_MaAsLin.csv")

###--------- 7g. Vomiting and diarrhoea (separate models; supplementary) -----------
run_maaslin <- function(outcome, outdir) {
  Maaslin2(
    input_data = as.data.frame(relab2), input_metadata = as.data.frame(meta_maas),
    output = outdir,
    fixed_effects = c(outcome,"antibiotics","eats_raw","wave","log_totalbases"),
    random_effects = "sample_id",
    normalization = "NONE", transform = "LOG", analysis_method = "LM",
    correction = "BH", max_significance = 0.05
  )
}
fit_vomiting  <- run_maaslin("vomiting",  "maaslin2_vomiting")
fit_diarrhoea <- run_maaslin("diarrhoea", "maaslin2_diarrhoea")

## ---- Section 7g.v2: MaAsLin2 — GI upset with sequencing depth as covariate ----
## Consolidates maaslin2_species_gi_bases_raw_GI from PCoA_PCA_GLMM_GI_AMR.R
## Requires: relab (species-level relative abundance matrix, rows = samples)
##           samp (metadata with gi, antibiotics, eats_raw, wave, log_totalbases, sample_id)

## Build relative abundance matrix (species level, post-filter samples)
relab_sp <- bracken |>
  filter(superkingdom == "Bacteria", confidence == 15) |>
  group_by(sample, name) |>
  summarise(abund = sum(fraction_total_reads, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = name, values_from = abund, values_fill = 0) |>
  column_to_rownames("sample") |>
  as.matrix()

## Build matched metadata (must have same row order as relab_sp)
meta_gi_depth <- samp |>
  dplyr::select(sample, gi, antibiotics, eats_raw, wave, log_totalbases, sample_id) |>
  filter(sample %in% rownames(relab_sp)) |>
  mutate(across(c(gi, antibiotics, eats_raw), as.logical),
         wave = as.factor(wave)) |>
  column_to_rownames("sample")

## Align row order
meta_gi_depth <- meta_gi_depth[rownames(relab_sp), , drop = FALSE]
stopifnot(identical(rownames(relab_sp), rownames(meta_gi_depth)))

##  Run MaAsLin2
fit_gi_depth <- Maaslin2(
  input_data     = as.data.frame(relab_sp),
  input_metadata = as.data.frame(meta_gi_depth),
  output         = "maaslin2_gi_depth",
  fixed_effects  = c("gi", "antibiotics", "eats_raw", "wave", "log_totalbases"),
  random_effects = c("sample_id"),
  normalization  = "NONE",
  transform      = "LOG",
  analysis_method = "LM",
  correction     = "BH",
  max_significance = 0.05
)

## 7gv.2 Read results back and save tidy summaries
gi_depth_all <- read_tsv("maaslin2_gi_depth/all_results.tsv")
gi_depth_sig <- read_tsv("maaslin2_gi_depth/significant_results.tsv")

## Depth-associated taxa specifically
depth_sig <- gi_depth_sig |>
  filter(metadata == "log_totalbases") |>
  arrange(desc(coef))

write_csv(gi_depth_all, "maaslin2_gi_depth_all_results.csv")
write_csv(depth_sig,    "maaslin2_depth_significant_taxa.csv")

message(sprintf(
  "Depth covariate MaAsLin2: %d taxa significant for log_totalbases (FDR<0.05)",
  nrow(depth_sig)
))

## 7g.v2. Quick summary for thesis section 2.7
gi_depth_sig |>
  count(metadata) |>
  arrange(desc(n)) |>
  print()

###---------- 7h. Wave (age) trend — MaAsLin2 + ANCOM-BC2 ------
meta_maas_wave <- meta_maas |>
  mutate(wave_num = as.numeric(as.character(wave)), wave = as.factor(wave))
meta_df2_wave  <- meta_df2  |>
  mutate(wave_num = as.numeric(as.character(wave)), wave = as.factor(wave))

fit_wave <- Maaslin2(
  input_data = as.data.frame(relab2), input_metadata = as.data.frame(meta_maas_wave),
  output = "maaslin2_wave_trend",
  fixed_effects = c("wave_num","gi","antibiotics","eats_raw","log_totalbases"),
  random_effects = "sample_id",
  normalization = "NONE", transform = "LOG", analysis_method = "LM",
  correction = "BH"
)
ancom_wave <- ancombc2(
  data = feat_df, meta_data = meta_df2_wave,
  fix_formula = "wave_num + gi + antibiotics + eats_raw + batch",
  group = NULL, taxa_are_rows = FALSE, p_adj_method = "BH",
  prv_cut = 0.10, struc_zero = FALSE, neg_lb = TRUE
)
# Save full results
write_csv(ancom_wave$res, "wave_trend_taxa.csv")
# Results stored in wave_trend_taxa.csv; 13 taxa decrease with wave in MaAsLin2 (FDR<0.05)

# Read back and save as a tidy CSV alongside the ANCOM output
wave_maaslin_sig <- read_tsv("maaslin2_wave_trend/significant_results.tsv") |>
  filter(metadata == "wave_num") |>
  arrange(qval)

write_csv(wave_maaslin_sig, "wave_trend_maaslin2_significant.csv")

# Quick summary
message(sprintf(
  "Wave trend MaAsLin2: %d significant taxa (FDR<0.05); %d decreasing, %d increasing",
  nrow(wave_maaslin_sig),
  sum(wave_maaslin_sig$coef < 0),
  sum(wave_maaslin_sig$coef > 0)
))

###------ 7i. Raw diet — MaAsLin2 + ANCOM-BC2-------
fit_raw <- Maaslin2(
  input_data = as.data.frame(relab2), input_metadata = as.data.frame(meta_maas),
  output = "maaslin2_raw_diet",
  fixed_effects = c("eats_raw","gi","antibiotics","wave","log_totalbases"),
  random_effects = "sample_id",
  normalization = "NONE", transform = "LOG", analysis_method = "LM",
  correction = "BH", max_significance = 0.05
)
ancom_raw <- ancombc2(
  data = feat_df, meta_data = meta_df2,
  fix_formula = "eats_raw + gi + antibiotics + wave + batch",
  group = "eats_raw", taxa_are_rows = FALSE, p_adj_method = "BH",
  prv_cut = 0.10, struc_zero = TRUE, neg_lb = TRUE
)

###------------ 7j. Antibiotic use — ANCOM-BC2 ---------
ancom_abx <- ancombc2(
  data = feat_df, meta_data = meta_df2,
  fix_formula = "antibiotics + gi + eats_raw + wave + batch",
  group = "antibiotics", taxa_are_rows = FALSE, p_adj_method = "BH",
  prv_cut = 0.10, struc_zero = TRUE, neg_lb = TRUE
)
# 14 taxa significant; 11 lower, 3 higher during antibiotic treatment

# ── 8. GLMMs: GI UPSET ────────────────────────────────────────────────────────
###------ 8a. Primary model (gi_final: GI ~ microbiome + metadata) ---------
gi_final <- glmmTMB(
  gi ~ shannon + PC1 + PC2 +
    eats_raw + antibiotics + prev_antibiotics + wave +
    (1 | sample_id) + (1 | batch),
  data   = samp_clr, family = binomial()
)
summary(gi_final)
# antibiotics p = 0.019 (OR ≈ 4.1); all microbiome terms ns
# Conclusion: microbiome composition/diversity does not predict GI beyond antibiotics

## 8b. Extended model (with depth + batch complexity)
gi_final_depth <- glmmTMB(
  gi ~ shannon + PC1 + PC2 +
    eats_raw + antibiotics + prev_antibiotics + wave +
    log_totalbases + n_plates +
    (1 | sample_id),
  data = samp_clr, family = binomial()
)
summary(gi_final_depth) # only recent antibiotics as significant

## 8c. Vomiting model (supplementary)
vomiting_final <- glmmTMB(
  vomiting ~ shannon + PC1 + PC2 +
    eats_raw + antibiotics + prev_antibiotics + wave +
    (1 | sample_id) + (1 | batch),
  data = samp_clr, family = binomial()
)
summary(vomiting_final)   # No significant predictors

## 8d. Diarrhoea model (supplementary)
diarr_final <- glmmTMB(
  diarrhoea ~ shannon + PC1 + PC2 +
    eats_raw + antibiotics + prev_antibiotics + wave +
    (1 | sample_id) + (1 | batch),
  data = samp_clr, family = binomial()
)
summary(diarr_final)   # antibiotics marginal (p = 0.064)

## 8e. Model comparison: with vs without PCs
gi_base <- glmmTMB(
  gi ~ eats_raw + antibiotics + scavenges + prev_antibiotics + wave +
    (1 | sample_id) + (1 | batch),
  data = samp_clr, family = binomial()
)
anova(gi_base, gi_final)   # p = 0.59 — PCs do not improve fit

# ── 9. GLMMs: AMR BURDEN ─────────────────────────────────────────────────────
###------ 9a. AMR total load ~ GI + wave (negative binomial, offset by depth) -----------
amr_gi_mod <- glmmTMB(
  amr_total ~ gi + antibiotics + eats_raw + wave_f +
    offset(log(total_bases)) + (1 | sample_id),
  data = samp, family = nbinom2()
)
summary(amr_gi_mod)
# Wave 2: exp(-0.325) ≈ 0.72 (~28% lower per Gb vs wave 1, p=0.0008)
# Wave 3: exp(-0.294) ≈ 0.75 (~25% lower, p=0.004); GI ns; antibiotics ns

###------ 9b. AMR richness (supplementary; convergence issues — interpret cautiously) ------
amr_rich_mod <- glmmTMB(
  amr_richness ~ gi + antibiotics + eats_raw + wave_f +
    offset(log(total_bases)) + (1 | sample_id),
  data = samp, family = nbinom2()
)
summary(amr_rich_mod)
# Wve 3: exp(0.274), p= 0.0394

# ── 10. ADDITIONAL PLOTS ──────────────────────────────────────────────────────
## Volcano plot — ANCOM-BC2 GI
volc_df <- ancom_gi$res |>
  transmute(taxon, lfc = lfc_giTRUE, q = q_giTRUE,
            neglog10q = -log10(q_giTRUE), sig = diff_giTRUE) |>
  filter(is.finite(neglog10q))

top_labels <- volc_df |> arrange(q) |> slice_head(n = 10)

ggplot(volc_df, aes(lfc, neglog10q)) +
  geom_point(aes(shape = sig), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_text_repel(data = top_labels, aes(label = taxon), size = 3) +
  theme_bw() +
  labs(x = "ANCOM-BC2 log fold-change (GI upset)", y = "-log10(FDR q-value)",
       shape = "Significant (FDR < 0.05)")

ggsave("ANCOMBC2_volcano_GI.png", width = 10, height = 7, dpi = 300)

## AMR saturation curve
grid_amr <- data.frame(
  total_bases = seq(min(samp$total_bases), max(samp$total_bases), length.out = 200)
)
grid_amr$pred <- predict(nls_fit, newdata = grid_amr)

ggplot(samp, aes(total_bases / 1e9, amr_richness)) +
  geom_point(alpha = 0.6) +
  geom_line(data = grid_amr, aes(total_bases / 1e9, pred), colour = "red") +
  geom_vline(xintercept = depth_90, linetype = "dashed",
             aes(colour = "90% saturation (~6 Gb)")) +
  theme_bw() +
  labs(x = "Sequencing depth (Gb)", y = "AMR gene richness",
       title = "AMR gene discovery: Michaelis-Menten saturation curve")

ggsave("AMR_saturation_curve.png", width = 9, height = 6, dpi = 300)

## PCoA coloured by GI upset
ggplot(ord_df, aes(PCo1, PCo2, colour = gi)) +
  geom_point(size = 2) +
  theme_bw() +
  labs(x = "PCoA1 (Aitchison)", y = "PCoA2", colour = "GI upset")

ggsave("PCoA_Aitchison_GI.png", width = 8, height = 6, dpi = 300)

## PCoA coloured by batch (for supplementary)
ggplot(ord_df, aes(PCo1, PCo2, colour = batch, shape = factor(wave))) +
  geom_point(size = 2) +
  theme_bw() +
  labs(x = "PCoA1", y = "PCoA2", colour = "Batch", shape = "Wave")

ggsave("PCoA_Aitchison_batch.png", width = 9, height = 6, dpi = 300)

# ── 11. MISSING DEDICATED MODELS: Raw Diet + Antibiotic ANCOM-BC2 & MaAsLin2 ──
# These three runs were absent from both the pre-filter and postfilter scripts.
# Run once; outputs are saved to disk and do not need re-running unless data changes.

###------ 11a. Raw diet — MaAsLin2 (eats_raw as outcome) ---------
# Uses the same relab2, meta_maas, and samp_clr objects from Section 7.
# meta_maas already has eats_raw as a factor column.

fit_raw_diet <- Maaslin2(
  input_data     = as.data.frame(relab2),
  input_metadata = as.data.frame(meta_maas),
  output         = "maaslin2_species_raw_diet",
  fixed_effects  = c("eats_raw", "gi", "antibiotics", "wave", "log_totalbases"),
  random_effects = c("sample_id"),
  normalization  = "NONE",
  transform      = "LOG",
  analysis_method= "LM",
  correction     = "BH",
  max_significance = 0.05
)

raw_diet_sig <- read_tsv("maaslin2_species_raw_diet/significant_results.tsv") |>
  filter(metadata == "eats_raw") |>
  arrange(qval)

raw_diet_all <- read_tsv("maaslin2_species_raw_diet/all_results.tsv") |>
  filter(metadata == "eats_raw") |>
  arrange(qval)

###------ 11b. Raw diet — ANCOM-BC2 (eats_raw as primary group) ---------
# Uses feat_df and meta_df2 from Section 7b.
# Note: meta_df2$eats_raw must be a factor/character column.

ancom_raw_diet <- ancombc2(
  data          = feat_df,
  meta_data     = meta_df2,
  fix_formula   = "eats_raw + gi + antibiotics + wave + batch",
  group         = "eats_raw",
  taxa_are_rows = FALSE,
  p_adj_method  = "BH",
  prv_cut       = 0.10,
  lib_cut       = 0,
  struc_zero    = TRUE,
  neg_lb        = TRUE
)

# Extract eats_rawTRUE results
raw_diet_ancom_res <- ancom_raw_diet$res |>
  transmute(
    taxon,
    ancom_lfc   = lfc_eats_rawTRUE,
    ancom_se    = se_eats_rawTRUE,
    ancom_q     = q_eats_rawTRUE,
    ancom_diff  = diff_eats_rawTRUE,
    ancom_sig   = q_eats_rawTRUE < 0.05 & diff_eats_rawTRUE,
    direction   = if_else(lfc_eats_rawTRUE > 0, "higher_in_raw", "lower_in_raw")
  ) |>
  arrange(ancom_q)

write_csv(raw_diet_ancom_res, "ancombc2_raw_diet_results.csv")

###------ 11c. Antibiotic use — ANCOM-BC2 (antibiotics as primary group) -----------
# antibiotics in meta_df2 is a factor (TRUE/FALSE for current antibiotics).

ancom_antibiotics <- ancombc2(
  data          = feat_df,
  meta_data     = meta_df2,
  fix_formula   = "antibiotics + gi + eats_raw + wave + batch",
  group         = "antibiotics",
  taxa_are_rows = FALSE,
  p_adj_method  = "BH",
  prv_cut       = 0.10,
  lib_cut       = 0,
  struc_zero    = TRUE,
  neg_lb        = TRUE
)

# Extract antibioticsTRUE results
abx_ancom_res <- ancom_antibiotics$res |>
  transmute(
    taxon,
    ancom_lfc   = lfc_antibioticsTRUE,
    ancom_se    = se_antibioticsTRUE,
    ancom_q     = q_antibioticsTRUE,
    ancom_diff  = diff_antibioticsTRUE,
    ancom_sig   = q_antibioticsTRUE < 0.05 & diff_antibioticsTRUE,
    direction   = if_else(lfc_antibioticsTRUE > 0,
                          "higher_during_antibiotics",
                          "lower_during_antibiotics")
  ) |>
  arrange(ancom_q)

write_csv(abx_ancom_res, "ancombc2_antibiotics_results.csv")

###-------- 11d. Quick summary prints --------
cat("\n=== Raw diet MaAsLin2 significant taxa ===\n")
print(raw_diet_sig |> dplyr::select(feature, coef, qval))

cat("\n=== Raw diet ANCOM-BC2 significant taxa ===\n")
print(raw_diet_ancom_res |> filter(ancom_sig))

cat("\n=== Antibiotic ANCOM-BC2 significant taxa ===\n")
print(abx_ancom_res |> filter(ancom_sig))

# ── 12. S. copri sensitivity analyses -------------------

###-------- 12a. Does the raw diet Shannon effect survive after controlling for S. copri? -------
# Add S. copri relative abundance as a covariate in the Shannon model.
# If eats_raw coefficient attenuates substantially, copri suppression is driving
# the diversity effect. If it barely changes, the diversity signal is independent.

samp_clr_copri <- samp_clr |>
  left_join(
    bracken |>
      filter(name == "Segatella copri") |>
      group_by(sample) |>
      summarise(copri_frac = sum(fraction_total_reads, na.rm = TRUE)),
    by = "sample"
  ) |>
  mutate(copri_frac = replace_na(copri_frac, 0))

# Model without S. copri (your existing result for reference)
m_raw_shannon_base <- glmmTMB(
  shannon ~ eats_raw + wave + antibiotics + prev_antibiotics +
    (1 | sample_id) + (1 | batch),
  data = samp_clr_copri, family = gaussian()
)
summary(m_raw_shannon_base)

# Model with S. copri fraction added
m_raw_shannon_copri <- glmmTMB(
  shannon ~ eats_raw + copri_frac + wave + antibiotics + prev_antibiotics +
    (1 | sample_id) + (1 | batch),
  data = samp_clr_copri, family = gaussian()
)
summary(m_raw_shannon_copri)
anova(m_raw_shannon_base, m_raw_shannon_copri) 

# Compare eats_rawTRUE coefficient between models
tidy(m_raw_shannon_base,  effects = "fixed") |> filter(term == "eats_rawTRUE")
tidy(m_raw_shannon_copri, effects = "fixed") |> filter(term == "eats_rawTRUE")

copri_comparison <- bind_rows(
  tidy(m_raw_shannon_base,  effects = "fixed") |> 
    filter(term == "eats_rawTRUE") |> mutate(model = "base"),
  tidy(m_raw_shannon_copri, effects = "fixed") |> 
    filter(term == "eats_rawTRUE") |> mutate(model = "copri_controlled")
)
write_csv(copri_comparison, "shannon_copri_sensitivity.csv")

###------------- 12b. Does raw diet ANCOM-BC signal survive on S. copri-excluded matrix? ---------
# Remove S. copri from the count matrix, re-normalise, re-run MaAsLin2.
# This directly tests whether enrichment of other taxa is real or compositional artefact.

relab2_no_copri <- relab2[, !grepl("Segatella.copri|Prevotella.copri", colnames(relab2))]
relab2_no_copri <- sweep(relab2_no_copri, 1, rowSums(relab2_no_copri), "/")

fit_raw_no_copri <- Maaslin2(
  input_data     = as.data.frame(relab2_no_copri),
  input_metadata = as.data.frame(meta_maas),
  output         = "maaslin2_raw_no_copri",
  fixed_effects  = c("eats_raw", "gi", "antibiotics", "wave", "log_totalbases"),
  random_effects = "sample_id",
  normalization  = "NONE", transform = "LOG",
  analysis_method = "LM", correction = "BH"
)
# If Bacteroides faecis, Parabacteroides distasonis, Phocaeicola vulgatus
# remain significant with similar coefficients → enrichment is genuine.
# If they disappear or flip → they were compositionally suppressed by S. copri.

###------ 12c. Kingdom-level diversity by raw diet ---------
# Viral and fungal alpha diversity per sample, then compare raw vs non-raw.

kingdom_alpha <- bracken |>
  filter(sample %in% samp$sample) |>
  group_by(sample, superkingdom) |>
  summarise(
    sk_shannon  = vegan::diversity(new_est_reads, index = "shannon"),
    sk_richness = vegan::specnumber(new_est_reads > 0),
    .groups = "drop"
  ) |>
  left_join(samp |> dplyr::select(sample, eats_raw, wave, antibiotics,
                                  sample_id, batch, total_bases),
            by = "sample")

# Model each kingdom separately
for (kingdom in c("Viruses", "Fungi", "Archaea")) {
  dat <- kingdom_alpha |>
    filter(superkingdom == kingdom, sk_richness > 0)
  
  if (nrow(dat) < 20) {
    message(kingdom, ": too few samples (n=", nrow(dat), ") — skip")
    next
  }
  
  m <- glmmTMB(
    sk_shannon ~ eats_raw + wave + log10(total_bases + 1) +
      (1 | sample_id) + (1 | batch),
    data = dat, family = gaussian()
  )
  cat("\n===", kingdom, "Shannon ~ eats_raw ===\n")
  print(tidy(m, effects = "fixed") |>
          dplyr::select(term, estimate, std.error, p.value))
}

# ── 13. Enterotype Shannon comparison ────────────────────────────────────────
# Requires: ent_alpha_all from enterotypes.R joined to postfilter samp table
# ent must have sample names cleaned back to underscore format (see enterotypes.R)

ent_genus <- ent_alpha_all |>
  filter(rank == "genus",
         enterotype %in% c("Bact_Phoc", "Firmicutes"),
         !is.na(shannon), !is.na(enterotype))

# Primary test
wilcox.test(shannon ~ enterotype, data = ent_genus, exact = FALSE)

# Effect size (rank-biserial correlation)
# install.packages("rstatix") if needed
rstatix::wilcox_effsize(shannon ~ enterotype, data = ent_genus)

# Descriptives
ent_genus |>
  group_by(enterotype) |>
  summarise(n = n(), median_shannon = median(shannon),
            IQR_low  = quantile(shannon, 0.25),
            IQR_high = quantile(shannon, 0.75))
#### ── 14. Enterotype GLMM ----------
# Requires ent_alpha_all from enterotypes.R 

library(glmmTMB)
library(broom.mixed)

# Binary outcome: Bact_Phoc (1) vs Firmicutes (0)
# Use genus-level rows only (one row per sample)
ent_genus <- ent_alpha_all |>
  filter(rank == "genus",
         enterotype %in% c("Bact_Phoc", "Firmicutes"),
         !is.na(enterotype))

ent_genus <- ent_genus |>
  mutate(
    is_bact_phoc = as.integer(enterotype == "Bact_Phoc"),
    wave_f       = factor(wave),
    eats_raw     = as.logical(eats_raw),
    antibiotics  = as.factor(antibiotics)
  )

# Primary model
m_enterotype <- glmmTMB(
  is_bact_phoc ~ eats_raw + wave_f + antibiotics + (1 | batch),
  data   = ent_genus,
  family = binomial()
)
summary(m_enterotype)
tidy(m_enterotype, effects = "fixed", exponentiate = TRUE, conf.int = TRUE)
# exponentiate=TRUE gives odds ratios directly

# Simpler model without antibiotics (likely too few events per cell)
m_enterotype_base <- glmmTMB(
  is_bact_phoc ~ eats_raw + wave_f + (1 | batch),
  data   = ent_genus,
  family = binomial()
)
summary(m_enterotype_base)

# Model comparison
anova(m_enterotype_base, m_enterotype)

# Descriptives to report alongside
ent_genus |>
  group_by(enterotype, eats_raw) |>
  summarise(n = n(), .groups = "drop") |>
  pivot_wider(names_from = eats_raw, values_from = n,
              names_prefix = "raw_")
##############################################################################
## END OF SCRIPT
##############################################################################

# GI upset GLMM (primary model)
summary(gi_final)

# Extended GI model with depth + batch complexity
summary(gi_final_depth)

# Vomiting and diarrhoea supplementary models
summary(vomiting_final)
summary(diarr_final)

# AMR wave model
summary(amr_gi_mod)

# Shannon ~ raw diet (base model without copri — your m_raw_shannon_base)
summary(m_raw_shannon_base)

# Wave trend ANCOM-BC2 object
summary(ancom_wave$res)   # or just: ancom_wave$res

summary(m_shannon_raw)
confint(m_shannon_raw, method = "Wald")
performance::r2(m_shannon_raw)
emmeans::emmeans(m_shannon_raw, ~ eats_raw)

summary(m_shannon_raw_sp)
confint(m_shannon_raw_sp, method = "Wald")
performance::r2(m_shannon_raw_sp)
