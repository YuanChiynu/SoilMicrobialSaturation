OUT   <- file.path("output", "Fig4")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "Fig4")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

library(igraph)
library(ggplot2)

n_sim <- 100
frac_grid <- seq(0, 0.8, by = 0.04)

natcon_adj <- function(adj) {
  N <- nrow(adj)
  if (N <= 1) return(0)
  evals <- eigen(adj, symmetric = TRUE, only.values = TRUE)$value
  m <- max(evals)
  lse <- m + log(sum(exp(evals - m))) - log(N)
  lse / (N - log(N))
}

nc_at_fractions <- function(ig, ord, frac_grid) {
  N <- vcount(ig)
  sapply(frac_grid, function(f) {
    k <- round(f * N)
    if (k <= 0) return(natcon_adj(as_adjacency_matrix(ig, sparse = FALSE)))
    k <- min(k, N - 1)
    tmp <- delete_vertices(ig, V(ig)$name[ord[seq_len(k)]])
    natcon_adj(as_adjacency_matrix(tmp, sparse = FALSE))
  })
}

norm_auc <- function(fr, nc) {
  nrm <- pmin(nc / nc[1], 1)
  auc <- sum(diff(fr) * (nrm[-1] + nrm[-length(nrm)]) / 2)
  auc / max(fr)
}

p50_of <- function(fr, nc) {
  nrm <- nc / nc[1]
  if (min(nrm, na.rm = TRUE) > 0.5) return(NA_real_)
  idx <- which(nrm <= 0.5)[1]
  if (idx <= 1) return(fr[1])
  x0 <- fr[idx - 1]; x1 <- fr[idx]
  y0 <- nrm[idx - 1]; y1 <- nrm[idx]
  x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)
}

safe_mean  <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else mean(x) }
safe_quant <- function(x, p) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else as.numeric(quantile(x, p)) }

eco_names <- c("Corpland", "Forest", "Grass")
eco_labs  <- c("Cropland", "Forest", "Grassland")
density_names <- c("10re", "40re", "70re", "100re")
density_labs  <- c("10%", "40%", "70%", "100%")
density_cols  <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")

.SKIP_COMPUTE <- !identical(Sys.getenv("RECOMPUTE"), "1") &&
  file.exists(file.path(CACHE, "robustness_curves.csv")) &&
  file.exists(file.path(CACHE, "auc_summary.csv"))

if (.SKIP_COMPUTE) {

  curves_df   <- read.csv(file.path(CACHE, "robustness_curves.csv"), stringsAsFactors = FALSE)
  auc_summary <- read.csv(file.path(CACHE, "auc_summary.csv"),       stringsAsFactors = FALSE)
  for (.f in c("robustness_curves.csv", "auc_summary.csv")) {
    .p <- file.path(CACHE, .f)
    cat(sprintf("[skip-compute] reading %-26s %6.1f KB  modified %s\n",
                .f, file.size(.p) / 1024, format(file.mtime(.p), "%Y-%m-%d %H:%M")))
  }
  cat("[skip-compute] random-attack simulation not rerun (set RECOMPUTE=1, about 1-3 h)\n")
} else {

set.seed(20260907)

auc_rows    <- list()
auc_sim_rows <- list()
curve_rows  <- list()

for (ei in seq_along(eco_names)) {
  eco <- eco_names[ei]
  gs <- list()
  for (si in seq_along(density_names)) {
    fn <- file.path(CACHE, paste0(eco, density_names[si], "_edge_screen0.0001_p0.05.csv"))
    ed <- read.csv(fn, header = TRUE)
    gs[[si]] <- graph_from_data_frame(ed[, c("source", "target")], directed = FALSE)
  }

  rand_nc <- lapply(seq_along(gs), function(si) {
    ig <- gs[[si]]; N <- vcount(ig)
    replicate(n_sim, { ord <- sample(N); nc_at_fractions(ig, ord, frac_grid) })
  })

  targ_nc <- lapply(seq_along(gs), function(si) {
    ig <- gs[[si]]
    ord <- order(degree(ig), betweenness(ig), decreasing = TRUE)
    nc_at_fractions(ig, ord, frac_grid)
  })

  for (si in seq_along(gs)) {
    mat <- rand_nc[[si]]
    nc_mean <- rowMeans(mat, na.rm = TRUE)
    nc_lo   <- apply(mat, 1, quantile, probs = 0.025, na.rm = TRUE)
    nc_hi   <- apply(mat, 1, quantile, probs = 0.975, na.rm = TRUE)
    auc_sim <- apply(mat, 2, function(nc) norm_auc(frac_grid, nc))
    p50_sim <- apply(mat, 2, function(nc) p50_of(frac_grid, nc))
    targ_auc <- norm_auc(frac_grid, targ_nc[[si]])
    targ_p50 <- p50_of(frac_grid, targ_nc[[si]])

    auc_sim_rows[[length(auc_sim_rows) + 1]] <- data.frame(
      ecosystem = eco_labs[ei], density = density_labs[si],
      sim = seq_len(n_sim), AUC = auc_sim, p50 = p50_sim,
      stringsAsFactors = FALSE)

    auc_rows[[length(auc_rows) + 1]] <- data.frame(
      ecosystem = eco_labs[ei], density = density_labs[si],
      random_AUC_mean = safe_mean(auc_sim),
      random_AUC_lo   = safe_quant(auc_sim, 0.025),
      random_AUC_hi   = safe_quant(auc_sim, 0.975),
      random_p50_mean = safe_mean(p50_sim),
      random_p50_lo   = safe_quant(p50_sim, 0.025),
      random_p50_hi   = safe_quant(p50_sim, 0.975),
      targeted_AUC    = targ_auc,
      targeted_p50    = targ_p50,
      stringsAsFactors = FALSE)

    curve_rows[[length(curve_rows) + 1]] <- data.frame(
      ecosystem = eco_labs[ei], density = density_labs[si], attack = "random",
      frac = frac_grid, nc_mean = nc_mean, nc_lo = nc_lo, nc_hi = nc_hi)
    curve_rows[[length(curve_rows) + 1]] <- data.frame(
      ecosystem = eco_labs[ei], density = density_labs[si], attack = "targeted",
      frac = frac_grid, nc_mean = targ_nc[[si]], nc_lo = NA, nc_hi = NA)
  }
}

auc_summary <- do.call(rbind, auc_rows)
auc_sim_df  <- do.call(rbind, auc_sim_rows)
curves_df   <- do.call(rbind, curve_rows)

kw_out <- list()
for (eco in eco_labs) {
  d_eco <- auc_sim_df[auc_sim_df$ecosystem == eco, ]
  d_eco$density <- factor(d_eco$density, levels = density_labs)
  kw_auc <- tryCatch(kruskal.test(AUC ~ density, data = d_eco), error = function(e) NULL)
  kw_p50 <- tryCatch(kruskal.test(p50 ~ density, data = d_eco), error = function(e) NULL)
  kw_out[[length(kw_out) + 1]] <- data.frame(
    ecosystem = eco, metric = "AUC",
    chi_sq = if (is.null(kw_auc)) NA_real_ else as.numeric(kw_auc$statistic),
    df     = if (is.null(kw_auc)) NA_real_ else as.numeric(kw_auc$parameter),
    p      = if (is.null(kw_auc)) NA_real_ else kw_auc$p.value)
  kw_out[[length(kw_out) + 1]] <- data.frame(
    ecosystem = eco, metric = "p50",
    chi_sq = if (is.null(kw_p50)) NA_real_ else as.numeric(kw_p50$statistic),
    df     = if (is.null(kw_p50)) NA_real_ else as.numeric(kw_p50$parameter),
    p      = if (is.null(kw_p50)) NA_real_ else kw_p50$p.value)
}
kw_stats <- do.call(rbind, kw_out)

auc_summary$ecosystem  <- factor(auc_summary$ecosystem, levels = eco_labs)
auc_summary$density <- factor(auc_summary$density, levels = density_labs)
write.csv(auc_summary, file = file.path(STATS, "auc_summary.csv"), row.names = FALSE)
write.csv(auc_sim_df, file = file.path(CACHE, "robustness_auc_sim.csv"), row.names = FALSE)
write.csv(kw_stats,   file = file.path(STATS, "robustness_kruskal.csv"), row.names = FALSE)
write.csv(curves_df,  file = file.path(CACHE, "robustness_curves.csv"), row.names = FALSE)

}

curves_df$density <- factor(curves_df$density, levels = density_labs)
curves_df$ecosystem  <- factor(curves_df$ecosystem, levels = eco_labs)

curves_df <- do.call(rbind, lapply(
  split(curves_df, list(curves_df$ecosystem, curves_df$density, curves_df$attack)),
  function(dd) {
    dd$nrm <- dd$nc_mean / dd$nc_mean[1]
    dd$lo  <- if (all(is.na(dd$nc_lo))) NA_real_ else dd$nc_lo / dd$nc_mean[1]
    dd$hi  <- if (all(is.na(dd$nc_hi))) NA_real_ else dd$nc_hi / dd$nc_mean[1]
    dd
  }))
rownames(curves_df) <- NULL

d_rand <- curves_df[curves_df$attack == "random", ]
d_rand$panel <- "Random attack"
d_targ <- curves_df[curves_df$attack == "targeted", ]
d_targ$panel <- "Sensitivity: degree-first attack"
d_comp <- rbind(d_rand, d_targ)
d_comp$panel <- factor(d_comp$panel, levels = c("Random attack", "Sensitivity: degree-first attack"))

.fit_quad_curve <- function(dd) {
  ok <- is.finite(dd$frac) & is.finite(dd$nrm)
  if (sum(ok) < 3) return(NULL)
  m <- suppressWarnings(lm(nrm ~ frac + I(frac^2), data = dd[ok, ]))
  data.frame(ecosystem = dd$ecosystem[ok][1], density = dd$density[ok][1],
             panel = dd$panel[ok][1],
             frac = dd$frac[ok], fit = as.numeric(predict(m)),
             stringsAsFactors = FALSE)
}
.fit_quad_stats <- function(dd) {
  ok <- is.finite(dd$frac) & is.finite(dd$nrm)
  if (sum(ok) < 3) return(NULL)
  m <- suppressWarnings(lm(nrm ~ frac + I(frac^2), data = dd[ok, ]))
  s <- summary(m)
  data.frame(ecosystem = dd$ecosystem[ok][1], density = dd$density[ok][1],
             panel = dd$panel[ok][1],
             b0 = unname(coef(m)[1]), b1 = unname(coef(m)[2]), b2 = unname(coef(m)[3]),
             r_sq = s$r.squared, adj_r_sq = s$adj.r.squared,
             stringsAsFactors = FALSE)
}
.fit_split  <- split(d_comp, list(d_comp$ecosystem, d_comp$density, d_comp$panel))
quad_df     <- do.call(rbind, lapply(.fit_split, .fit_quad_curve))
quad_stats  <- do.call(rbind, lapply(.fit_split, .fit_quad_stats))
write.csv(quad_stats, file = file.path(STATS, "robustness_quadratic_fit.csv"),
          row.names = FALSE)

fmt_p50 <- function(x) ifelse(is.na(x), ">0.8", sprintf("%.2f", x))

.k_sat <- as.integer(factor(as.character(auc_summary$density), levels = density_labs))
.annot_rows <- function(panel, p50_lab, auc_lab, x_p50, x_auc, y_tops) {
  data.frame(ecosystem = auc_summary$ecosystem, density = auc_summary$density,
             panel = panel,
             lab = c(sprintf("p50 = %s", p50_lab), sprintf("AUC = %.3f", auc_lab)),
             x = rep(c(x_p50, x_auc), each = length(.k_sat)),
             ypos = rep(y_tops[.k_sat], 2), stringsAsFactors = FALSE)
}
annot_df <- rbind(
  .annot_rows("Random attack",
              fmt_p50(auc_summary$random_p50_mean), auc_summary$random_AUC_mean,
              0.019, 0.265, c(0.496, 0.388, 0.280, 0.172)),
  .annot_rows("Sensitivity: degree-first attack",
              fmt_p50(auc_summary$targeted_p50), auc_summary$targeted_AUC,
              0.307, 0.546, c(1.328, 1.218, 1.108, 0.997)))
annot_df$panel <- factor(annot_df$panel, levels = c("Random attack", "Sensitivity: degree-first attack"))

.theme_rob <- function() {
  theme_bw(base_size = 8) +
    theme(panel.grid         = element_blank(),
          legend.position    = "top",
          legend.key         = element_blank(),
          legend.text        = element_text(size = 8),
          legend.title       = element_text(size = 8),
          strip.background   = element_blank(),
          strip.text.x       = element_text(size = 7, face = "bold", colour = "black",
                                            margin = margin(t = 2, r = 3, b = 3, l = 3)),
          strip.text.y.right = element_text(size = 7, face = "bold", colour = "black",
                                            margin = margin(t = 3, r = 2, b = 3, l = 4)),
          plot.margin        = margin(4, 4, 4, 4, "pt"))
}

p_comp <- ggplot(d_comp, aes(x = frac, y = nrm, color = density)) +

  geom_line(size = 1.05) +

  geom_line(data = quad_df, aes(x = frac, y = fit, color = density),
            linetype = "dashed", linewidth = 0.55, inherit.aes = FALSE) +
  scale_color_manual(name = "Propagule density(%)", values = density_cols) +

  scale_fill_manual(values = density_cols, guide = "none") +

  geom_text(data = annot_df, aes(x = x, y = ypos, label = lab, color = density),
            size = 2.46, hjust = 0, vjust = 1, show.legend = FALSE, inherit.aes = FALSE) +
  expand_limits(y = c(0, max(1.30, max(d_comp$hi, na.rm = TRUE) * 1.02))) +
  facet_grid(panel ~ ecosystem, scales = "free_y") +
  labs(x = "Proportions of removed nodes", y = "Normalized natural connectivity") +
  .theme_rob()

ggsave(file = file.path(OUT, "Roboutness_composite.pdf"), p_comp, width = 6.3, height = 3.7)

p_rand <- ggplot(d_rand, aes(x = frac, y = nrm, color = density)) +

  geom_line(size = 1.1) +

  geom_line(data = subset(quad_df, panel == "Random attack"),
            aes(x = frac, y = fit, color = density),
            linetype = "dashed", linewidth = 0.55, inherit.aes = FALSE) +
  scale_color_manual(values = density_cols) + scale_fill_manual(values = density_cols, guide = "none") +
  facet_wrap(~ ecosystem, ncol = 3) +
  labs(x = "Proportions of removed nodes", y = "Normalized natural connectivity") +
  .theme_rob()

p_targ <- ggplot(d_targ, aes(x = frac, y = nrm, color = density)) +
  geom_line(size = 1.1) +

  geom_line(data = subset(quad_df, panel == "Sensitivity: degree-first attack"),
            aes(x = frac, y = fit, color = density),
            linetype = "dashed", linewidth = 0.55, inherit.aes = FALSE) +
  scale_color_manual(values = density_cols) +
  facet_wrap(~ ecosystem, ncol = 3) +
  labs(x = "Proportions of removed nodes", y = "Normalized natural connectivity") +
  .theme_rob()

leg_p <- ggplot(curves_df[curves_df$attack == "random", ],
                aes(x = frac, y = nrm, color = density)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  scale_color_manual(values = density_cols, name = "Propagule density (%)") +
  theme_bw(base_size = 8) + theme(legend.position = "right")
if (exists(".save_legend")) .save_legend(leg_p, "Fig4_Robustness_legend", OUT)

cat("Robustness analysis done\n")
