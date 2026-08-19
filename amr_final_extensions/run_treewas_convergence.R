#!/usr/bin/env Rscript

# Frozen discovery/validation phylogenetic convergence analysis for the
# K. pneumoniae-colistin project.
#
# The genotype panel and discovery/validation split are prepared without using
# phenotype associations. Discovery candidates must be significant under the
# treeWAS simulation-based Bonferroni gate, have evolutionary support from the
# simultaneous or subsequent test, and be resistance-directed. Only those frozen
# candidates are opened in the BioProject-disjoint validation cohort. A survivor
# remains an association, not a causal or novel resistance mechanism.

suppressPackageStartupMessages({
  library(ape)
  library(treeWAS)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_treewas_convergence.R PANEL_DIR OUT_DIR SEED")
}
panel_dir <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
seed <- as.integer(args[[3]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path, ...) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

read_genotypes <- function(path) {
  x <- read_tsv(path)
  if (!("sample_id" %in% names(x))) stop("Missing sample_id in ", path)
  ids <- as.character(x$sample_id)
  x$sample_id <- NULL
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  rownames(m) <- ids
  if (anyNA(m)) stop("NA genotype values in ", path)
  if (any(!(m %in% c(0, 1)))) stop("Non-binary genotype values in ", path)
  m
}

read_phenotypes <- function(path) {
  x <- read_tsv(path)
  if (!all(c("sample_id", "phenotype_binary") %in% names(x))) {
    stop("Phenotype file lacks required fields: ", path)
  }
  y <- as.numeric(x$phenotype_binary)
  names(y) <- as.character(x$sample_id)
  y
}

read_distance <- function(path) {
  d <- read_tsv(path, row.names = 1)
  m <- as.matrix(d)
  storage.mode(m) <- "numeric"
  if (nrow(m) != ncol(m)) stop("Distance matrix is not square")
  if (!identical(sort(rownames(m)), sort(colnames(m)))) stop("Distance IDs differ")
  m <- m[rownames(m), rownames(m), drop = FALSE]
  m <- (m + t(m)) / 2
  diag(m) <- 0
  if (any(!is.finite(m))) stop("Non-finite distance value")
  if (any(m < -1e-12)) stop("Negative distance value")
  m[m < 0] <- 0
  m
}

align_inputs <- function(geno, phen, dist) {
  ids <- rownames(geno)
  common <- ids[ids %in% names(phen) & ids %in% rownames(dist)]
  if (length(common) != length(ids)) {
    stop("Exact-ID mismatch: genotype=", length(ids), " common=", length(common))
  }
  list(
    geno = geno[common, , drop = FALSE],
    phen = phen[common],
    dist = dist[common, common, drop = FALSE]
  )
}

build_tree <- function(dist) {
  tree <- ape::bionj(as.dist(dist))
  tree$tip.label <- as.character(tree$tip.label)
  if (!setequal(tree$tip.label, rownames(dist))) stop("Tree tip mismatch")
  tree
}

# treeWAS 1.1.1 may return SNP.locus as a one-based numeric column index even
# when the input matrix has explicit pattern names. Normalize every reported
# locus against the exact genotype column order before any set operation.
normalize_locus_ids <- function(values, genotype_columns, context) {
  values <- as.character(values)
  out <- values
  numeric_positions <- which(grepl("^[0-9]+$", values))
  if (length(numeric_positions)) {
    indices <- suppressWarnings(as.integer(values[numeric_positions]))
    valid <- !is.na(indices) & indices >= 1L & indices <= length(genotype_columns)
    if (any(!valid)) {
      stop(context, ": treeWAS returned out-of-range numeric locus IDs: ",
           paste(values[numeric_positions][!valid], collapse = ","))
    }
    out[numeric_positions] <- genotype_columns[indices]
  }
  unknown <- setdiff(unique(out), genotype_columns)
  if (length(unknown)) {
    stop(context, ": unresolved treeWAS locus IDs: ", paste(unknown, collapse = ","))
  }
  unique(out)
}

extract_sig <- function(result, test_name, genotype_columns) {
  obj <- result[[test_name]]$sig.snps
  if (is.null(obj) || length(obj) == 0) {
    return(data.frame(pattern_id = character(), test = character(),
                      score = numeric(), p_value = numeric()))
  }
  if (is.atomic(obj) && !is.data.frame(obj)) {
    raw <- as.character(obj[!is.na(obj)])
    raw <- raw[raw != "No significant SNPs found."]
    if (!length(raw)) {
      return(data.frame(pattern_id = character(), test = character(),
                        score = numeric(), p_value = numeric()))
    }
    vals <- normalize_locus_ids(raw, genotype_columns,
                                paste0(test_name, " atomic sig.snps"))
    return(data.frame(pattern_id = vals, test = test_name,
                      score = NA_real_, p_value = NA_real_))
  }
  obj <- as.data.frame(obj, stringsAsFactors = FALSE)
  if (!("SNP.locus" %in% names(obj))) {
    return(data.frame(pattern_id = character(), test = character(),
                      score = numeric(), p_value = numeric()))
  }
  raw_ids <- as.character(obj$SNP.locus)
  normalized_per_row <- raw_ids
  numeric_positions <- which(grepl("^[0-9]+$", raw_ids))
  if (length(numeric_positions)) {
    indices <- suppressWarnings(as.integer(raw_ids[numeric_positions]))
    valid <- !is.na(indices) & indices >= 1L & indices <= length(genotype_columns)
    if (any(!valid)) {
      stop(test_name, ": treeWAS returned out-of-range numeric locus IDs: ",
           paste(raw_ids[numeric_positions][!valid], collapse = ","))
    }
    normalized_per_row[numeric_positions] <- genotype_columns[indices]
  }
  unknown <- setdiff(unique(normalized_per_row), genotype_columns)
  if (length(unknown)) {
    stop(test_name, ": unresolved treeWAS locus IDs: ", paste(unknown, collapse = ","))
  }
  keep <- !duplicated(normalized_per_row)
  data.frame(
    pattern_id = normalized_per_row[keep],
    test = test_name,
    score = if ("score" %in% names(obj)) as.numeric(obj$score[keep]) else NA_real_,
    p_value = if ("p.value" %in% names(obj)) as.numeric(obj$p.value[keep]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

combined_sig <- function(result, genotype_columns) {
  x <- result$treeWAS.combined$treeWAS.combined
  if (is.null(x)) return(character())
  normalize_locus_ids(x[!is.na(x)], genotype_columns, "treeWAS combined")
}

haldane_or <- function(a, b, c, d) {
  v <- as.numeric(c(a, b, c, d))
  if (any(v == 0)) v <- v + 0.5
  (v[[1]] * v[[4]]) / (v[[2]] * v[[3]])
}

contingency_table <- function(geno, phen, candidates = colnames(geno)) {
  candidates <- intersect(candidates, colnames(geno))
  if (length(candidates) == 0) return(data.frame())
  rows <- lapply(candidates, function(id) {
    g <- as.numeric(geno[, id])
    y <- as.numeric(phen[rownames(geno)])
    a <- sum(g == 1 & y == 1)
    b <- sum(g == 1 & y == 0)
    c <- sum(g == 0 & y == 1)
    d <- sum(g == 0 & y == 0)
    p <- tryCatch(
      fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
                  alternative = "greater")$p.value,
      error = function(e) NA_real_
    )
    data.frame(pattern_id = id, R_present = a, S_present = b,
               R_absent = c, S_absent = d,
               odds_ratio = haldane_or(a, b, c, d),
               fisher_one_sided_p = p, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$fisher_BH_q <- p.adjust(out$fisher_one_sided_p, method = "BH")
  out
}

run_treewas <- function(geno, phen, dist, label, seed_value) {
  tree <- build_tree(dist)
  if (!setequal(rownames(geno), tree$tip.label)) stop(label, ": tree/genotype mismatch")
  tree <- keep.tip(tree, rownames(geno))
  n_sim <- max(10000L, as.integer(ncol(geno) * 10L))
  chunk <- min(500L, ncol(geno))
  message(label, ": samples=", nrow(geno), " patterns=", ncol(geno),
          " simulated_null_loci=", n_sim)
  set.seed(seed_value)
  result <- treeWAS(
    snps = geno,
    phen = phen,
    tree = tree,
    phen.type = "discrete",
    n.snps.sim = n_sim,
    chunk.size = chunk,
    mem.lim = FALSE,
    test = c("terminal", "simultaneous", "subsequent"),
    snps.reconstruction = "parsimony",
    snps.sim.reconstruction = "parsimony",
    phen.reconstruction = "parsimony",
    p.value = 0.05,
    p.value.correct = "bonf",
    p.value.by = "count",
    plot.tree = FALSE,
    plot.manhattan = FALSE,
    plot.null.dist = FALSE,
    plot.dist = FALSE,
    seed = seed_value
  )
  saveRDS(result, file.path(out_dir, paste0(label, "_treewas_result.rds")))
  tests <- do.call(rbind, lapply(c("terminal", "simultaneous", "subsequent"),
                                function(x) extract_sig(result, x, colnames(geno))))
  write.csv(tests, file.path(out_dir, paste0(label, "_significant_by_test.csv")),
            row.names = FALSE)
  combined <- combined_sig(result, colnames(geno))
  writeLines(combined, file.path(out_dir, paste0(label, "_combined_significant.txt")))
  list(result = result, tests = tests, combined = combined,
       n_sim = n_sim, tree = tree)
}

manifest <- read.csv(file.path(panel_dir, "manifest.csv"), stringsAsFactors = FALSE)
meta <- read.csv(file.path(panel_dir, "pattern_metadata.csv"), stringsAsFactors = FALSE,
                 check.names = FALSE)
disc <- align_inputs(
  read_genotypes(file.path(panel_dir, "discovery_genotypes.tsv")),
  read_phenotypes(file.path(panel_dir, "discovery_phenotypes.tsv")),
  read_distance(file.path(panel_dir, "discovery_distance.tsv"))
)
val <- align_inputs(
  read_genotypes(file.path(panel_dir, "validation_genotypes.tsv")),
  read_phenotypes(file.path(panel_dir, "validation_phenotypes.tsv")),
  read_distance(file.path(panel_dir, "validation_distance.tsv"))
)

capture.output(sessionInfo(), file = file.path(out_dir, "R_SESSION_INFO.txt"))
writeLines(as.character(packageVersion("treeWAS")),
           file.path(out_dir, "TREEWAS_VERSION.txt"))

# Discovery is performed once on the phenotype-blind frozen panel.
disc_run <- run_treewas(disc$geno, disc$phen, disc$dist,
                        "discovery", seed)
disc_counts <- contingency_table(disc$geno, disc$phen, disc_run$combined)
evolutionary_ids <- unique(disc_run$tests$pattern_id[
  disc_run$tests$test %in% c("simultaneous", "subsequent")
])
frozen_ids <- intersect(disc_run$combined, evolutionary_ids)
if (nrow(disc_counts)) {
  frozen_ids <- intersect(frozen_ids,
                          disc_counts$pattern_id[disc_counts$odds_ratio > 1])
}
writeLines(frozen_ids, file.path(out_dir, "FROZEN_DISCOVERY_CONVERGENCE_CANDIDATES.txt"))

if (length(frozen_ids) == 0) {
  summary <- list(
    status = "NO_DISCOVERY_CONVERGENCE_CANDIDATE",
    n_discovery_samples = nrow(disc$geno),
    n_validation_samples = nrow(val$geno),
    n_input_patterns = ncol(disc$geno),
    n_discovery_combined_significant = length(disc_run$combined),
    n_discovery_evolutionary_supported_resistance_directed = 0,
    n_validation_replicated = 0,
    strict_candidates = character(),
    boundary = paste(
      "No frozen pattern passed the simulation-based Bonferroni discovery gate,",
      "an evolutionary treeWAS test, and a resistance-directed odds ratio.",
      "No new resistance marker is claimed."
    )
  )
  writeLines(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE),
             file.path(out_dir, "TREEWAS_CONVERGENCE_SUMMARY.json"))
  quit(status = 0)
}

# Validation is opened only after the candidate IDs are frozen.
val_ids <- intersect(frozen_ids, colnames(val$geno))
if (length(val_ids) != length(frozen_ids)) {
  stop("Frozen candidate missing from validation matrix")
}
val_geno <- val$geno[, val_ids, drop = FALSE]
val_run <- run_treewas(val_geno, val$phen, val$dist,
                       "validation_frozen", seed + 1L)
val_counts <- contingency_table(val_geno, val$phen, val_ids)

all_tests <- merge(
  data.frame(pattern_id = frozen_ids, stringsAsFactors = FALSE),
  disc_counts,
  by = "pattern_id", all.x = TRUE, suffixes = c("", "_discovery")
)
names(all_tests)[names(all_tests) %in% c("R_present", "S_present", "R_absent", "S_absent",
                                        "odds_ratio", "fisher_one_sided_p", "fisher_BH_q")] <-
  paste0(names(all_tests)[names(all_tests) %in% c("R_present", "S_present", "R_absent", "S_absent",
                                                  "odds_ratio", "fisher_one_sided_p", "fisher_BH_q")],
         "_discovery")
val_counts2 <- val_counts
names(val_counts2)[names(val_counts2) != "pattern_id"] <-
  paste0(names(val_counts2)[names(val_counts2) != "pattern_id"], "_validation")
all_tests <- merge(all_tests, val_counts2, by = "pattern_id", all.x = TRUE)
all_tests$discovery_combined <- all_tests$pattern_id %in% disc_run$combined
all_tests$discovery_simultaneous <- all_tests$pattern_id %in%
  disc_run$tests$pattern_id[disc_run$tests$test == "simultaneous"]
all_tests$discovery_subsequent <- all_tests$pattern_id %in%
  disc_run$tests$pattern_id[disc_run$tests$test == "subsequent"]
all_tests$validation_combined <- all_tests$pattern_id %in% val_run$combined
all_tests$validation_simultaneous <- all_tests$pattern_id %in%
  val_run$tests$pattern_id[val_run$tests$test == "simultaneous"]
all_tests$validation_subsequent <- all_tests$pattern_id %in%
  val_run$tests$pattern_id[val_run$tests$test == "subsequent"]
all_tests$strict_treewas_replication <- with(all_tests,
  validation_combined &
  (validation_simultaneous | validation_subsequent) &
  odds_ratio_validation > 1 &
  fisher_BH_q_validation <= 0.05
)
all_tests <- merge(all_tests, meta, by = "pattern_id", all.x = TRUE)
all_tests <- all_tests[order(!all_tests$strict_treewas_replication,
                             all_tests$fisher_BH_q_validation,
                             -all_tests$odds_ratio_validation), ]
write.csv(all_tests, file.path(out_dir, "DISCOVERY_VALIDATION_TREEWAS_EVIDENCE.csv"),
          row.names = FALSE)
strict <- all_tests[all_tests$strict_treewas_replication %in% TRUE, , drop = FALSE]
write.csv(strict, file.path(out_dir, "STRICT_TREEWAS_REPLICATED_PATTERNS.csv"),
          row.names = FALSE)

summary <- list(
  status = if (nrow(strict))
    "TREEWAS_REPLICATES_REQUIRE_KNOWN_MECHANISM_CONTEXT_AND_NOVELTY_AUDIT"
  else "NO_PATTERN_SURVIVED_DISCOVERY_AND_VALIDATION_TREEWAS_GATE",
  n_discovery_samples = nrow(disc$geno),
  n_validation_samples = nrow(val$geno),
  n_input_patterns = ncol(disc$geno),
  n_discovery_combined_significant = length(disc_run$combined),
  n_discovery_evolutionary_supported_resistance_directed = length(frozen_ids),
  n_validation_combined_significant = length(val_run$combined),
  n_strict_treewas_replicated = nrow(strict),
  strict_candidates = as.character(strict$pattern_id),
  boundary = paste(
    "A replicated treeWAS pattern is a phylogenetically supported association,",
    "not a causal or novel resistance mechanism. It must be integrated with",
    "known-mechanism, exact sequence/context, external-data and literature audits,",
    "and causality still requires laboratory validation."
  )
)
writeLines(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE),
           file.path(out_dir, "TREEWAS_CONVERGENCE_SUMMARY.json"))

report <- c(
  "# Frozen treeWAS discovery-validation convergence audit",
  "",
  paste0("- Discovery samples: **", nrow(disc$geno), "**"),
  paste0("- Validation samples: **", nrow(val$geno), "**"),
  paste0("- Input occurrence patterns: **", ncol(disc$geno), "**"),
  paste0("- Discovery combined significant: **", length(disc_run$combined), "**"),
  paste0("- Frozen evolutionary-supported, resistance-directed candidates: **", length(frozen_ids), "**"),
  paste0("- Complete discovery-validation treeWAS gate: **", nrow(strict), "**"),
  "",
  "## Claim boundary",
  "",
  summary$boundary
)
writeLines(report, file.path(out_dir, "TREEWAS_CONVERGENCE_REPORT.md"))
