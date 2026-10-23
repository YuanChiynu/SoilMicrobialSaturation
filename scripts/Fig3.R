DATA  <- "data"

GRID <- identical(Sys.getenv("GRID"), "1")
OUT  <- file.path("output", "Fig3")
POUT <- if (GRID) file.path(OUT, "grid") else OUT
STATS <- file.path("output", "Fig3", "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "Fig3")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(POUT,  showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(ggplot2)
  library(ggpubr)
  library(data.table)
})
if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr is required")

META      <- file.path(DATA, "data", "metagenome")
NMDS_SEED <- 20260915
PERM_SEED <- 20260901
PERM_N    <- 9999
CYCLES    <- c("carbon", "nitrogen", "phosphorus", "sulfur")

DENSITY_LEVELS   <- c("10%", "100%")
SAT_COL      <- c("10%" = "#ffd460", "100%" = "#2d4059")
TREAT_LEVELS <- c("CK", "C", "N", "pH", "T", "W")
TREAT_LABELS <- c("CK", "C", "N", "pH", "W", "D")

TREAT_SHAPES <- c(CK = 3, C = 15, N = 16, pH = 17, W = 18, D = 4)

NMDS_W <- 6.0; NMDS_H <- 6.0
BAR_W  <- 3.0; BAR_H  <- 6.0
COMP_W <- 2 * (NMDS_W + BAR_W)
COMP_H <- 2 * NMDS_H

.read_meta <- function() {
  gc <- data.table::fread(file.path(META, "reads_ra.csv.gz"), data.table = FALSE)
  gc <- as.data.frame(gc)[, c("sample", "KO", "RA")]
  names(gc) <- c("sample", "ko", "ra")
  gs <- aggregate(gc$ra, by = list(gc$sample, gc$ko), FUN = "sum")
  names(gs) <- c("sample", "ko", "gene_copies")
  hm <- tidyr::spread(gs, key = "sample", value = "gene_copies", fill = 0)
  KO <- read.csv(file.path(META, "KO_annotated.csv"), fileEncoding = "UTF-8-BOM")
  input <- merge(hm, KO, by.x = "ko", by.y = "KEGG_KO")

  rownames(input) <- make.unique(as.character(input$Gene_name))
  grp <- read.csv(file.path(META, "group.csv"), fileEncoding = "UTF-8-BOM")
  list(input = input, grp = grp)
}
META_DAT <- .read_meta()
INPUT    <- META_DAT$input
GRP      <- META_DAT$grp
stopifnot(all(GRP$sample %in% names(INPUT)))
stopifnot(all(c("sample", "satru", "treat") %in% names(GRP)))

SAMPLE_COLS <- names(INPUT)[names(INPUT) %in% GRP$sample]
stopifnot(length(SAMPLE_COLS) == nrow(GRP))
cat(sprintf("[Fig3] sample order = spread order (first 6: %s)\n",
            paste(head(SAMPLE_COLS, 6), collapse = ",")))

compute_cycle <- function(cyc) {
  th <- subset(INPUT, INPUT$cycle == cyc)
  if (nrow(th) < 3) stop("too few cycle genes: ", cyc)
  hd <- as.data.frame(t(th[, SAMPLE_COLS, drop = FALSE]))
  set.seed(NMDS_SEED)
  nm <- metaMDS(hd, autotransform = FALSE)
  site <- as.data.frame(nm$points)[, 1:2]
  names(site) <- c("NMDS1", "NMDS2")
  site$sample <- rownames(site)
  site <- merge(site, GRP, by = "sample", all.x = TRUE)
  site$satru <- factor(site$satru, levels = DENSITY_LEVELS)
  site$treat <- factor(site$treat, levels = TREAT_LEVELS, labels = TREAT_LABELS)
  stopifnot(!any(is.na(site$satru)), !any(is.na(site$treat)))

  dm <- as.matrix(vegdist(hd, method = "bray"))
  list(cycle = cyc, site = site, stress = nm$stress, dm = dm, n_gene = nrow(th))
}

compute_sim <- function(res) {
  g <- res$site

  rows <- do.call(rbind, lapply(DENSITY_LEVELS, function(s) {
    ck <- g$sample[g$satru == s & g$treat == "CK"]
    ot <- g$sample[g$satru == s & g$treat != "CK"]
    if (!length(ck) || !length(ot)) stop("CK / non-CK samples missing: ", s)
    dd <- res$dm[ck, ot, drop = FALSE]
    sim_ck <- 1 - rowMeans(dd)
    sim_ot <- 1 - colMeans(dd)
    stopifnot(abs(mean(c(sim_ck, sim_ot)) - mean(1 - dd)) < 1e-12)
    data.frame(satru = s,
               sample = c(ck, ot),
               role   = c(rep("CK", length(ck)), rep("non-CK", length(ot))),
               similarity = c(sim_ck, sim_ot),
               stringsAsFactors = FALSE)
  }))
  rows$satru <- factor(rows$satru, levels = DENSITY_LEVELS)

  tt <- t.test(similarity ~ satru, data = rows)
  agg <- do.call(rbind, lapply(DENSITY_LEVELS, function(s) {
    v <- rows$similarity[rows$satru == s]
    data.frame(satru = s, n = length(v), mean = mean(v), sd = sd(v),
               se = sd(v) / sqrt(length(v)), stringsAsFactors = FALSE)
  }))
  list(data = rows, agg = agg, tt = tt)
}

.p_star <- function(p) {
  ifelse(is.na(p), "ns",
    ifelse(p <= 0.001, "***",
      ifelse(p <= 0.01, "**",
        ifelse(p <= 0.05, "*", "ns"))))
}

compute_permanova <- function(res) {
  idx <- match(rownames(res$dm), GRP$sample)
  stopifnot(!any(is.na(idx)), nrow(res$dm) == nrow(GRP))
  pm <- data.frame(
    density   = factor(GRP$satru[idx], levels = DENSITY_LEVELS),
    treatment = factor(GRP$treat[idx], levels = TREAT_LEVELS, labels = TREAT_LABELS),
    stringsAsFactors = FALSE)
  stopifnot(!any(is.na(pm$density)), !any(is.na(pm$treatment)))
  set.seed(PERM_SEED)
  fit <- vegan::adonis2(res$dm ~ density * treatment, data = pm,
                        permutations = PERM_N, by = "terms")
  sm <- as.data.frame(fit)
  sm$term_raw <- trimws(rownames(sm))
  .pcol <- grep("Pr\\(>F\\)", names(sm), value = TRUE)[1]
  .fcol <- if ("F" %in% names(sm)) "F" else grep("^F$", names(sm), value = TRUE)[1]
  .r2col <- if ("R2" %in% names(sm)) "R2" else grep("R2", names(sm), value = TRUE)[1]
  stopifnot(!is.na(.pcol), !is.na(.fcol), !is.na(.r2col))

  resid_df <- sm$Df[sm$term_raw %in% c("Residual", "Residuals")][1]
  stopifnot(is.finite(resid_df), resid_df > 0)
  lab <- c("density"           = "RS (propagule density)",
           "treatment"         = "T (treatment)",
           "density:treatment" = "RS x T")
  keep <- sm$term_raw %in% names(lab)
  stopifnot(sum(keep) == 3)
  out <- data.frame(
    cycle = res$cycle, n_gene = res$n_gene, nmds_stress = res$stress,
    term_raw = sm$term_raw[keep], term = unname(lab[sm$term_raw[keep]]),
    df1 = sm$Df[keep], df2 = resid_df,
    F_value = sm[[.fcol]][keep], R2 = sm[[.r2col]][keep],
    p_value = sm[[.pcol]][keep],
    test = sprintf("vegan::adonis2(bray ~ density * treatment, permutations = %d, by = 'terms', set.seed(%d))",
                   PERM_N, PERM_SEED),
    stringsAsFactors = FALSE)
  out$symbol <- .p_star(out$p_value)
  rownames(out) <- NULL
  out
}

.ellipse_extent <- function(x, y, level = 0.95) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 5 || stats::sd(x) == 0 || stats::sd(y) == 0)
    return(list(x = range(x), y = range(y)))
  m <- cbind(x, y); cm <- colMeans(m); S <- stats::cov(m); n <- nrow(m)
  r <- sqrt(2 * stats::qf(level, 2, n - 1))
  th <- seq(0, 2 * pi, length.out = 200)
  eig <- eigen(S, symmetric = TRUE)
  E <- cbind(cos(th), sin(th)) %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
  E <- sweep(E * r, 2, cm, "+")
  list(x = range(E[, 1]), y = range(E[, 2]))
}

.panel_extent <- function(site) {
  ex <- range(site$NMDS1); ey <- range(site$NMDS2)
  for (g in levels(site$satru)) {
    d <- site[site$satru == g, ]
    if (nrow(d) < 5) next
    e <- .ellipse_extent(d$NMDS1, d$NMDS2)
    ex <- range(c(ex, e$x)); ey <- range(c(ey, e$y))
  }
  list(x = ex, y = ey)
}

.global_limits <- function(sites, pad = 0.04) {
  ext <- lapply(sites, .panel_extent)
  all <- range(c(unlist(lapply(ext, `[[`, "x")), unlist(lapply(ext, `[[`, "y"))))
  c0 <- mean(all)
  span <- diff(all) * (1 + pad)
  c(c0 - span / 2, c0 + span / 2, c0 - span / 2, c0 + span / 2)
}

.theme_panel <- theme(
  panel.background = element_rect(fill = "transparent", colour = NA),
  panel.grid       = element_blank(),
  panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.35),
  axis.line        = element_blank(),
  axis.ticks       = element_line(colour = "black", linewidth = 0.3),
  axis.ticks.length = unit(1.2, "pt"),
  axis.text        = element_text(size = 6, colour = "black"),
  axis.title       = element_text(size = 8, colour = "black"),
  legend.position  = "none",
  plot.margin      = margin(2, 2, 2, 2)
)

plot_nmds <- function(res, lims, perm = NULL) {
  g <- res$site
  p <- ggplot(g, aes(x = NMDS1, y = NMDS2)) +

    geom_point(aes(colour = satru, shape = treat), size = 1.5) +
    stat_ellipse(aes(colour = satru, fill = satru), level = 0.95, geom = "polygon",
                 alpha = 0.12, linewidth = 0.3, show.legend = FALSE) +
    scale_shape_manual(values = TREAT_SHAPES, name = "Treatments") +
    scale_colour_manual(values = SAT_COL, name = "Propagule density(%)") +
    scale_fill_manual(values = SAT_COL, name = "Propagule density(%)") +
    coord_fixed(ratio = 1, xlim = lims[1:2], ylim = lims[3:4], clip = "off") +

    annotate("text", x = lims[2] - 0.035 * diff(lims[1:2]),
                     y = lims[4] - 0.035 * diff(lims[3:4]),
             label = sprintf("Stress = %.3f", res$stress),
             size = 2.1, hjust = 1, vjust = 1) +
    labs(x = "NMDS 1", y = "NMDS 2") +
    theme_grey(base_size = 6) + .theme_panel

  if (GRID) {
    p <- p + theme(panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
                   panel.grid.minor = element_line(colour = "grey94", linewidth = 0.18))
  }

  if (!is.null(perm)) {
    .sym_of <- function(tm) perm$symbol[perm$term_raw == tm][1]
    p <- p + annotate("text",
             x = lims[1] + 0.035 * diff(lims[1:2]),
             y = lims[4] - 0.035 * diff(lims[3:4]),
             label = sprintf("RS%s  T%s", .sym_of("density"), .sym_of("treatment")),
             size = 2.1, hjust = 0, vjust = 1)
  }
  p
}

plot_sim <- function(sim) {
  stars <- .p_star(sim$tt$p.value)
  ggbarplot(sim$data, x = "satru", y = "similarity", add = "mean_se",
            color = "black", fill = "satru",
            position = position_dodge(0.8), width = 0.7) +
    scale_fill_manual(values = SAT_COL, name = "Propagule density(%)") +
    labs(x = "Propagule density(%)", y = "Bray-Curtis Similarity") +

    scale_y_continuous(limits = c(0, 1.10), breaks = seq(0, 1, 0.25)) +

    annotate("text", x = 1.5, y = 1.00, label = stars,
             size = 2.1, hjust = 0.5, vjust = 0) +
    theme_grey(base_size = 6) + .theme_panel
}

plist    <- list()
nres_all <- list()
test_all <- list()
sims_all <- list()
perm_all <- list()
for (i in seq_along(CYCLES)) {
  cyc <- CYCLES[i]
  res <- compute_cycle(cyc)
  sim <- compute_sim(res)
  pv  <- compute_permanova(res)
  nres_all[[cyc]] <- res
  sims_all[[cyc]] <- sim
  perm_all[[cyc]] <- pv
  cat(sprintf("[Fig3] (%s) %-4s genes=%-3d  stress=%.3f  similarity %.3f/%.3f  %s (t=%.2f, p=%.2e)\n",
              letters[i], cyc, res$n_gene, res$stress,
              sim$agg$mean[1], sim$agg$mean[2], .p_star(sim$tt$p.value),
              unname(sim$tt$statistic), sim$tt$p.value))
  cat(sprintf("        PERMANOVA(RS x T):  RS %s (F=%.2f, p=%.4g)   T %s (F=%.2f, p=%.4g)   RSxT %s (F=%.2f, p=%.4g)\n",
              pv$symbol[pv$term_raw == "density"],           pv$F_value[pv$term_raw == "density"],           pv$p_value[pv$term_raw == "density"],
              pv$symbol[pv$term_raw == "treatment"],         pv$F_value[pv$term_raw == "treatment"],         pv$p_value[pv$term_raw == "treatment"],
              pv$symbol[pv$term_raw == "density:treatment"], pv$F_value[pv$term_raw == "density:treatment"], pv$p_value[pv$term_raw == "density:treatment"]))
}

COMMON_LIMS <- .global_limits(lapply(nres_all, `[[`, "site"))
cat(sprintf("[Fig3] shared axis limits: x [%.3f, %.3f]  y [%.3f, %.3f] (span %.3f)\n",
            COMMON_LIMS[1], COMMON_LIMS[2], COMMON_LIMS[3], COMMON_LIMS[4],
            COMMON_LIMS[2] - COMMON_LIMS[1]))

for (i in seq_along(CYCLES)) {
  cyc <- CYCLES[i]
  res <- nres_all[[cyc]]
  sim <- sims_all[[cyc]]

  p_nmds <- plot_nmds(res, COMMON_LIMS, perm_all[[cyc]])
  p_sim  <- plot_sim(sim)

  if (!exists(".fig3_leg_nmds")) { .save_legend(p_nmds, "Fig3_NMDS_legend", POUT); .fig3_leg_nmds <- TRUE }
  if (!exists(".fig3_leg_sim"))  { .save_legend(p_sim,  "Fig3_similarity_legend", POUT); .fig3_leg_sim <- TRUE }

  ggsave(file.path(POUT, paste0("Fig3_", letters[i], "_NMDS.pdf")), p_nmds,
         width = NMDS_W, height = NMDS_H, units = "cm", useDingbats = FALSE)
  ggsave(file.path(POUT, paste0("Fig3_", letters[i], "_bar.pdf")), p_sim,
         width = BAR_W, height = BAR_H, units = "cm", useDingbats = FALSE)

  s10 <- res$site[res$site$satru == "10%", ]
  cov10   <- cov(s10$NMDS1, s10$NMDS2)
  d10     <- res$site[res$site$satru == "10%" & res$site$treat == "D", ]
  d_below <- as.integer(mean(d10$NMDS2) < mean(s10$NMDS2))

  test_all[[i]] <- data.frame(
    panel = paste0("(", letters[i], ")"), cycling = cyc,
    n_gene = res$n_gene, nmds_stress = res$stress,
    cov_sign_10 = sign(cov10), d10_below = d_below,
    n_10 = sim$agg$n[1], n_100 = sim$agg$n[2],
    mean_10 = sim$agg$mean[1], se_10 = sim$agg$se[1],
    mean_100 = sim$agg$mean[2], se_100 = sim$agg$se[2],
    t = unname(sim$tt$statistic), df = unname(sim$tt$parameter),
    p_value = sim$tt$p.value, stars = .p_star(sim$tt$p.value),
    test = paste0("Unpaired Welch two-sample t-test on sample-level units ",
                  "(n = 24 per density; bar = mean over all CK x non-CK pairs)"),
    stringsAsFactors = FALSE)

  plist[[i]] <- ggarrange(p_nmds, p_sim, ncol = 2, widths = c(NMDS_W, BAR_W),
                          labels = paste0("(", letters[i], ")"),
                          font.label = list(size = 7, face = "plain"))
}

ST <- do.call(rbind, test_all)

if (!GRID) {
  write.csv(ST, file.path(CACHE, "fig3_similarity_test.csv"), row.names = FALSE)
  if (dir.exists(STATS)) write.csv(ST, file.path(STATS, "fig3_similarity_test.csv"), row.names = FALSE)
}
cat("\n[Fig3] similarity test summary -> ", file.path(CACHE, "fig3_similarity_test.csv"),
    if (GRID) "(GRID mode: not rewritten; uses the delivered version)" else "", "\n")
print(ST[, c("panel", "cycling", "n_gene", "nmds_stress",
             "mean_10", "mean_100", "t", "df", "p_value", "stars")], row.names = FALSE)

PERM <- do.call(rbind, perm_all)
PERM$panel <- paste0("(", letters[match(PERM$cycle, CYCLES)], ")")
stopifnot(nrow(PERM) == 12, !any(is.na(PERM$p_value)),
          !any(is.na(PERM$F_value)), all(PERM$df2 == 36))
PERM_OUT <- PERM[, c("panel", "cycle", "n_gene", "term", "term_raw", "df1", "df2",
                     "F_value", "R2", "p_value", "symbol", "test")]
if (!GRID) write.csv(PERM_OUT, file.path(STATS, "fig3_permanova_ST.csv"), row.names = FALSE)
cat("\n[Fig3] PERMANOVA (RS x T) summary -> ", file.path(CACHE, "fig3_permanova_ST.csv"), "\n")
print(PERM_OUT[, c("panel", "cycle", "term", "df1", "df2", "F_value", "R2", "p_value", "symbol")],
      row.names = FALSE)

library(grid)
TOPPAD_PT <- 20

COMP_NAME <- if (GRID) "Fig3_10vs100_grid.pdf" else "Fig3_10vs100.pdf"
fig3 <- ggarrange(plotlist = plist, ncol = 2, nrow = 2)

fig3 <- gtable::gtable_add_rows(ggplot2::ggplotGrob(fig3), unit(TOPPAD_PT, "pt"), pos = 0)
ggsave(file.path(OUT, COMP_NAME), fig3,
       width = COMP_W, height = COMP_H + TOPPAD_PT / 28.3465,
       units = "cm", useDingbats = FALSE)
cat(sprintf("\n[Fig3] output -> %s (%.1f x %.2f cm, %.0f pt stamp band at top)\n",
            file.path(OUT, COMP_NAME), COMP_W,
            COMP_H + TOPPAD_PT / 28.3465, TOPPAD_PT))
cat("[Fig3] done\n")
