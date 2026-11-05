pdf(NULL); dev.off()

DATA  <- "data"
OUT   <- file.path("output", "FigS15")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "FigS15")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

library(igraph)
library(ggplot2)
library(reshape2)

eco_names <- c("Corpland", "Forest", "Grass")
eco_labs  <- c("Cropland", "Forest", "Grassland")
density_names <- c("10re", "40re", "70re", "100re")
density_labs  <- c("10%", "40%", "70%", "100%")
density_cols  <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")

thr_vals  <- c(0.55, 0.65, 0.75, 0.85)
thr_labs  <- sprintf("%.2f", thr_vals)
eps       <- 1e-9

pos.mean <- function(v){ x <- v[v > 0]; if (!length(x)) 0 else mean(x) }
neg.mean <- function(v){ x <- v[v < 0]; if (!length(x)) 0 else mean(x) }

load_relabund_one <- function(eco, sat) {

  eco16 <- c(Corpland = "C", Forest = "F", Grass = "G")[[eco]]

  ot16 <- read.csv(file.path(DATA, "amplicon/prokaryote/OTU_rare_16S_all.csv"),
                   header = TRUE, row.names = 1, check.names = FALSE)
  ot16 <- as.data.frame(t(ot16))
  ot16$names <- rownames(ot16)
  cols16 <- setdiff(colnames(ot16), "names")

  colnames(ot16)[match(cols16, colnames(ot16))] <- paste0(cols16, "16s")
  cols16 <- paste0(cols16, "16s")
  grp16 <- read.csv(file.path(DATA, "amplicon/prokaryote/group.csv"), header = TRUE)
  m16 <- merge(ot16, grp16, by.x = "names", by.y = "sample")
  m16 <- m16[m16$ecosystem == eco16 & m16$group == sat, c("mental", cols16)]

  otit <- read.csv(file.path(DATA, "amplicon/fungi/OTU_rare_ITS.csv"),
                   header = TRUE, row.names = 1, check.names = FALSE)
  otit <- as.data.frame(t(otit))
  otit$names <- rownames(otit)
  colsit <- setdiff(colnames(otit), "names")
  colnames(otit)[match(colsit, colnames(otit))] <- paste0(colsit, "its")
  colsit <- paste0(colsit, "its")
  grpit <- read.csv(file.path(DATA, "amplicon/fungi/group.csv"), header = TRUE)
  mit <- merge(otit, grpit, by.x = "names", by.y = "sample")
  mit <- mit[mit$ecosystem == eco & mit$group == sat, c("mental", colsit)]

  common <- intersect(m16$mental, mit$mental)
  if (length(common) < length(union(m16$mental, mit$mental)))
    cat(sprintf("  [%s %s] mental intersect: %d of 16S %d / ITS %d\n",
                eco, sat, length(common), nrow(m16), nrow(mit)))
  stopifnot(length(common) >= 12)
  m16 <- m16[match(common, m16$mental), ]
  mit <- mit[match(common, mit$mental), ]
  X <- cbind(as.matrix(m16[, -1]), as.matrix(mit[, -1]))
  rownames(X) <- m16$mental

  otu <- t(X)
  otu1 <- otu; otu1[otu1 > 0] <- 1
  otu <- otu[which(rowSums(otu1) >= ncol(otu1) / 2), ]
  rowsum <- as.data.frame(rowSums(otu))
  sumtot <- colSums(rowsum)[1]
  RA <- rowsum / sumtot
  otu <- otu[which(RA[[1]] > 0.0001), ]
  otu <- t(otu)
  relabund <- otu / rowSums(otu)

  R <- suppressWarnings(cor(otu, method = "spearman"))
  R[is.na(R)] <- 0
  diag(R) <- 0
  list(relabund = relabund, R = R, otus = colnames(relabund))
}

cat("loading relabund for 12 subsets ...\n")
relabund_store <- list()
for (eco in eco_names) for (sat in density_names) {
  tag <- paste0(eco, sat)
  relabund_store[[tag]] <- load_relabund_one(eco, sat)
  cat("  ", tag, " OTUs(filtered) =", ncol(relabund_store[[tag]]$relabund), "\n")
}

tax_lookup <- data.frame()
for (eco in eco_names) for (sat in density_names) {
  nf <- file.path("output", "Fig4",
                  paste0(eco, sat, "_node_screen0.0001_p0.05.csv"))
  nd <- read.csv(nf, header = TRUE, stringsAsFactors = FALSE)
  tax_lookup <- rbind(tax_lookup, nd[, c("id", "kingdom", "phylum", "genus")])
}
tax_lookup <- tax_lookup[!duplicated(tax_lookup$id), ]

safe <- function(expr, alt = NA_real_, warn_is_na = TRUE) {

  if (warn_is_na) {
    tryCatch(expr, error = function(e) alt, warning = function(w) alt)
  } else {
    withCallingHandlers(
      tryCatch(expr, error = function(e) alt),
      warning = function(w) invokeRestart("muffleWarning"))
  }
}

metric_rows <- list()
for (ei in seq_along(eco_names)) {
  eco <- eco_names[ei]
  for (si in seq_along(density_names)) {
    sat <- density_names[si]
    tag <- paste0(eco, sat)
    ed_full <- read.csv(
      file.path("output", "Fig4", paste0(tag, "_edge_screen0.0001_p0.05.csv")),
      header = TRUE, stringsAsFactors = FALSE)
    rb <- relabund_store[[tag]]$relabund
    R  <- relabund_store[[tag]]$R

    for (thr in thr_vals) {
      ed <- ed_full[abs(ed_full$correlation) >= thr - eps, ]
      if (nrow(ed) == 0) {
        metric_rows[[length(metric_rows) + 1]] <- data.frame(
          ecosystem = eco_labs[ei], density = density_labs[si], threshold = thr,
          nodes = 0, edges = 0, pos_edges = 0, neg_edges = 0, neg_pos_ratio = NA,
          connectance = NA, average_degree = NA, average_path_length = NA,
          clustering_coefficient = NA, modularity = NA, degree_assortativity = NA,
          pos_cohesion_mean = NA, pos_cohesion_sd = NA,
          neg_cohesion_mean = NA, neg_cohesion_sd = NA)
        next
      }

      ed <- ed[ed$source != ed$target, ]
      ed$key <- paste(pmin(ed$source, ed$target), pmax(ed$source, ed$target), sep = "|")
      ed <- ed[!duplicated(ed$key), ]

      g <- graph_from_data_frame(ed[, c("source", "target", "correlation")], directed = FALSE)
      g <- simplify(g, edge.attr.comb = list(correlation = "first", weight = "first"))
      E(g)$weight <- abs(E(g)$correlation)

      pos_e <- sum(E(g)$correlation > 0); neg_e <- sum(E(g)$correlation < 0)
      N <- vcount(g); M <- ecount(g)

      set.seed(20261004 + ei * 1000 + si * 100 + match(thr, thr_vals))
      mod_louvain <- safe(igraph::modularity(
        g, membership(igraph::cluster_louvain(g, weights = E(g)$weight))))

      nodes_g <- V(g)$name
      cidx <- match(nodes_g, colnames(rb))
      if (anyNA(cidx)) stop("some nodes are absent from the relative-abundance matrix")
      sub_rb <- rb[, cidx, drop = FALSE]
      conn_pos_all <- apply(R, 2, pos.mean)
      conn_neg_all <- apply(R, 2, neg.mean)
      pos_coh <- as.vector(sub_rb %*% conn_pos_all[cidx])
      neg_coh <- as.vector(sub_rb %*% conn_neg_all[cidx])

      metric_rows[[length(metric_rows) + 1]] <- data.frame(
        ecosystem = eco_labs[ei], density = density_labs[si], threshold = thr,
        nodes = N, edges = M, pos_edges = pos_e, neg_edges = neg_e,
        neg_pos_ratio = ifelse(pos_e > 0, neg_e / pos_e, NA_real_),
        connectance = safe(edge_density(g, loops = FALSE)),
        average_degree = safe(mean(degree(g))),
        average_path_length = safe(average.path.length(g), warn_is_na = FALSE),
        clustering_coefficient = safe(transitivity(g)),
        modularity = mod_louvain,
        degree_assortativity = safe(assortativity_degree(g, directed = FALSE)),
        pos_cohesion_mean = mean(pos_coh), pos_cohesion_sd = sd(pos_coh),
        neg_cohesion_mean = mean(neg_coh), neg_cohesion_sd = sd(neg_coh))
      cat(sprintf("done %s | thr=%.2f | nodes=%d edges=%d\n", tag, thr, N, M))
    }
  }
}
metrics <- do.call(rbind, metric_rows)
metrics$ecosystem  <- factor(metrics$ecosystem,  levels = eco_labs)
metrics$density <- factor(metrics$density, levels = density_labs)
write.csv(metrics, file.path(STATS, "threshold_sensitivity_metrics.csv"), row.names = FALSE)

pat_rows <- list()
for (mtr in c("edges", "nodes", "average_degree", "connectance",
              "neg_pos_ratio", "modularity", "pos_cohesion_mean")) {
  for (eco in eco_labs) {
    for (thr in thr_vals) {
      v <- metrics[metrics$ecosystem == eco & metrics$threshold == thr, c("density", mtr)]
      v <- v[order(factor(v$density, levels = density_labs)), ]
      vals <- v[[mtr]]
      ok_monotonic <- if (any(is.na(vals))) NA else all(diff(vals) < 0)
      ok_10_gt_100 <- if (any(is.na(vals))) NA else vals[1] > vals[4]
      pat_rows[[length(pat_rows) + 1]] <- data.frame(
        metric = mtr, ecosystem = eco, threshold = thr,
        value_10 = vals[1], value_40 = vals[2], value_70 = vals[3], value_100 = vals[4],
        monotonic_10to100 = ok_monotonic, gt_10_vs_100 = ok_10_gt_100)
    }
  }
}
pattern <- do.call(rbind, pat_rows)
write.csv(pattern, file.path(STATS, "threshold_pattern_summary.csv"), row.names = FALSE)

n_pass <- sum(pattern$gt_10_vs_100 %in% TRUE)
n_tot  <- nrow(pattern)
cat(sprintf("ASSERT pattern '10%% > 100%%': %d / %d (metric × eco × thr)\n", n_pass, n_tot))
for (mtr in c("edges", "nodes", "average_degree", "connectance")) {
  sub <- pattern[pattern$metric == mtr, ]
  cat(sprintf("ASSERT %-16s : %d/%d thresholds×ecosystems keep 10%%>100%%; monotonic %d/%d\n",
              mtr, sum(sub$gt_10_vs_100 %in% TRUE), nrow(sub),
              sum(sub$monotonic_10to100 %in% TRUE), nrow(sub)))
}

metric_levels <- c("edges", "nodes", "average_degree", "connectance",
                   "neg_pos_ratio", "modularity", "pos_cohesion_mean", "neg_cohesion_mean")
metric_labs <- c(
  edges = "Number of edges", nodes = "Number of nodes",
  average_degree = "Average degree", connectance = "Connectance",
  neg_pos_ratio = "Neg:pos edge ratio", modularity = "Modularity",
  pos_cohesion_mean = "Positive cohesion", neg_cohesion_mean = "Negative cohesion")

long <- melt(metrics, id.vars = c("ecosystem", "density", "threshold"),
             measure.vars = metric_levels, variable.name = "metric")
long$metric <- factor(long$metric, levels = metric_levels, labels = metric_levels)
long$x_thr <- factor(sprintf("%.2f", long$threshold), levels = thr_labs)

.brk5 <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (!all(is.finite(r)) || diff(r) == 0) return(r[1])
  step <- signif(diff(r) / 4, 1)
  lo <- floor(r[1] / step) * step
  b <- lo + step * (0:4)
  while (max(b) < r[2] - 1e-9) b <- b + step
  b[1:5]
}

p <- ggplot(long, aes(x = x_thr, y = value,
                      color = density, group = density)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_color_manual(values = density_cols, name = "Propagule density (%)", drop = FALSE) +

  scale_y_continuous(breaks = .brk5) +
  facet_grid(metric ~ ecosystem, scales = "free_y",
             labeller = labeller(metric = setNames(metric_labs, metric_levels))) +
  labs(x = expression("Spearman |" * rho * "| threshold"),
       y = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(angle = 0, size = 6),
        legend.position = "none")

ggsave(file.path(OUT, "FigS15_composite.pdf"), p, width = 6.3, height = 9.4)
ggsave(file.path(CACHE, "FigS15_composite_p1.png"), p, width = 6.3, height = 9.4, dpi = 300)

leg_p <- ggplot(long, aes(x = x_thr, y = value,
                          color = density, group = density)) +
  geom_line(linewidth = 1) + geom_point(size = 2.6) +
  scale_color_manual(values = density_cols, name = "Propagule density (%)", drop = FALSE) +
  theme(legend.position = "bottom")
gt_leg <- ggplotGrob(leg_p)
idx_leg <- which(gt_leg$layout$name == "guide-box-bottom")
leg_grob <- gt_leg$grobs[[idx_leg]]
leg_pdf <- file.path(OUT, "FigS15_legend.pdf")
leg_png <- file.path(CACHE, "FigS15_legend.png")
pdf(leg_pdf, width = 14, height = 1.8); grid::grid.newpage(); grid::grid.draw(leg_grob); dev.off()
png(leg_png, width = 14, height = 1.8, units = "in", res = 300)
grid::grid.newpage(); grid::grid.draw(leg_grob); dev.off()
cat("[legend] FigS15_legend written via gtable guide-box-bottom\n")

cat("FigS15 done ->", file.path(OUT, "FigS15_composite.pdf"), "\n")
