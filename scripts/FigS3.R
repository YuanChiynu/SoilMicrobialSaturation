DATA  <- "data"

GRID  <- identical(Sys.getenv("GRID"), "1")
OUT   <- if (GRID) file.path("output", "FigS3", "grid") else file.path("output", "FigS3")
STATS <- file.path("output", "FigS3", "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "FigS3")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

library(vegan)
library(ggplot2)
library(ggpubr)

ECOS_LEVELS <- c("Croplands", "Forests", "Grasslands")
DENSITY_LEVELS  <- c("10%", "40%", "70%", "100%")
SAT_COLORS  <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")

compute_verify_pcoa <- function(otu_file, group_file, ecos_filter = NULL) {
  group <- read.csv(group_file, header = TRUE, stringsAsFactors = FALSE)

  group <- subset(group, treat %in% c("10%-D", "40%-D", "70%-D", "100%-D"))
  if (!is.null(ecos_filter)) group <- subset(group, ecosystem == ecos_filter)

  otu <- read.csv(otu_file, header = TRUE, row.names = 1, check.names = FALSE)
  otu <- as.data.frame(t(otu))

  common_samples <- intersect(rownames(otu), group$sample)
  otu <- otu[common_samples, , drop = FALSE]
  group <- group[group$sample %in% common_samples, ]

  otu <- otu[, colSums(otu, na.rm = TRUE) > 0, drop = FALSE]

  dist_bc <- vegdist(otu, method = "bray")
  pcoa <- cmdscale(dist_bc, k = 2, eig = TRUE)
  pc <- round(pcoa$eig / sum(pcoa$eig) * 100, digits = 2)

  pc12 <- as.data.frame(pcoa$points[, 1:2])
  colnames(pc12) <- c("V1", "V2")
  pc12$sample <- rownames(pc12)

  sample_site <- merge(pc12, group, by.x = "sample", by.y = "sample", all.x = TRUE)
  sample_site$treat <- factor(sample_site$treat,
                              levels = c("10%-D", "40%-D", "70%-D", "100%-D"))
  sample_site$Density <- factor(gsub("-D", "", sample_site$treat),
                                   levels = DENSITY_LEVELS)
  sample_site$Ecosystem <- factor(sample_site$ecosystem, levels = ECOS_LEVELS)
  sample_site$Density <- droplevels(sample_site$Density)
  sample_site$Ecosystem <- droplevels(sample_site$Ecosystem)
  stopifnot(!any(is.na(sample_site$Density)),
            !any(is.na(sample_site$Ecosystem)),
            nlevels(sample_site$Ecosystem) == 1)

  ord <- match(labels(dist_bc), sample_site$sample)
  stopifnot(!anyNA(ord), length(ord) == attr(dist_bc, "Size"))
  set.seed(20261008)
  sat <- factor(as.character(sample_site$Density[ord]), levels = DENSITY_LEVELS)
  ad  <- adonis2(dist_bc ~ sat, data = data.frame(sat = sat), permutations = 9999)
  bn  <- trimws(rownames(ad))
  perm <- data.frame(
    term  = "RS",
    df1   = ad$Df[1],
    df2   = ad$Df[which(bn %in% c("Residual", "Residuals"))],
    F_val = ad$F[1],
    R2    = ad$R2[1],
    p     = ad$`Pr(>F)`[1],
    stringsAsFactors = FALSE)
  stopifnot(is.finite(perm$df2), is.finite(perm$p))
  perm$symbol <- if (perm$p < 0.001) "***" else if (perm$p < 0.01) "**" else
                 if (perm$p < 0.05)  "*"   else "ns"

  list(site = sample_site, pc = pc, perm = perm)
}

.ellipse_extent <- function(x, y, level = 0.95) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 4) return(list(x = range(x), y = range(y)))
  m <- cbind(x, y); cm <- colMeans(m)
  S <- tryCatch(stats::cov(m), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S)) || det(S) <= 0)
    return(list(x = range(x), y = range(y)))
  r <- sqrt(2 * stats::qf(level, 2, n - 1))
  th <- seq(0, 2 * pi, length.out = 200)
  eig <- eigen(S, symmetric = TRUE)
  E <- cbind(cos(th), sin(th)) %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
  E <- sweep(E * r, 2, cm, "+")
  list(x = range(E[, 1]), y = range(E[, 2]))
}

.pcoa_own_lim <- function(site, level = 0.95) {
  xs <- numeric(0); ys <- numeric(0)
  for (g in levels(site$Density)) {
    d <- site[site$Density == g, ]
    e <- .ellipse_extent(d$V1, d$V2, level)
    xs <- c(xs, e$x, d$V1); ys <- c(ys, e$y, d$V2)
  }
  gx <- range(xs); gy <- range(ys)
  cx <- mean(gx); cy <- mean(gy)
  span <- max(diff(gx), diff(gy)) * 1.04
  c(cx - span / 2, cx + span / 2)
}

plot_verify_pcoa <- function(site, pc, out_pdf, lim, title, perm = NULL) {
  sample_site <- site
  p <- ggplot(sample_site, aes(x = V1, y = V2)) +
    geom_point(aes(color = Density), size = 1.6) +

    scale_color_manual(values = SAT_COLORS) +

    stat_ellipse(aes(group = Density, color = Density),
                 level = 0.95, linetype = 1, size = 0.6) +
    labs(x = paste0("PC1 ", pc[1], "%"),
         y = paste0("PC2 ", pc[2], "%"),
         color = "Density",
         title = title) +
    coord_fixed(ratio = 1, xlim = lim, ylim = lim, clip = "off") +

    theme(panel.grid.major = if (GRID) element_line(colour = "grey88", linewidth = 0.25) else element_blank(),
          panel.grid.minor = if (GRID) element_line(colour = "grey94", linewidth = 0.18) else element_blank(),
          panel.background = element_rect(color = NA, fill = 'transparent'),
          panel.border = element_rect(fill = NA, color = "black", size = 0.25),
          axis.ticks = element_line(color = "black", linewidth = 0.25),
          axis.text = element_text(size = 6, face = "plain"),
          axis.title = element_text(size = 8, face = "plain"),
          plot.title = element_text(size = 8, face = "plain", hjust = 0.5),
          legend.position = "none")

  if (!is.null(perm)) {
    .lab <- if (perm$symbol == "ns") "RS ns" else paste0("RS", perm$symbol)
    .sp  <- diff(lim)
    p <- p + annotate("text",
                      x = lim[1] + 0.03 * .sp, y = lim[2] - 0.03 * .sp,
                      label = .lab, hjust = 0, vjust = 1,
                      size = 6 / 2.845, colour = "black")
  }

  ggsave(out_pdf, p, width = 4.0, height = 4.0, units = "cm")
  return(p)
}

load_richness_data <- function(pairs, ecos_levels = ECOS_LEVELS) {
  alpha_list <- lapply(pairs, function(pr) {
    group <- read.csv(pr$group, header = TRUE, stringsAsFactors = FALSE)
    group <- subset(group, treat %in% c("10%-D", "40%-D", "70%-D", "100%-D"))
    otu <- read.csv(pr$otu, header = TRUE, row.names = 1, check.names = FALSE)
    otu <- as.data.frame(t(otu))
    common_samples <- intersect(rownames(otu), group$sample)
    otu <- otu[common_samples, , drop = FALSE]
    group <- group[group$sample %in% common_samples, ]
    stopifnot(nrow(otu) > 0, nrow(group) > 0)
    richness <- specnumber(otu)
    m <- match(names(richness), group$sample)
    data.frame(richness = as.numeric(richness), sample = names(richness),
               ecosystem = group$ecosystem[m], treat = group$treat[m],
               stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, alpha_list)
  d$treat <- factor(d$treat, levels = c("10%-D", "40%-D", "70%-D", "100%-D"))
  d$Density <- factor(gsub("-D", "", as.character(d$treat)), levels = DENSITY_LEVELS)
  d$Ecosystem <- factor(d$ecosystem, levels = ecos_levels)
  d <- droplevels(d)
  stopifnot(!any(is.na(d$Density)), !any(is.na(d$Ecosystem)))
  d
}

plot_verify_richness <- function(d, out_pdf, marker) {

  y0 <- max(d$richness)
  y_limit  <- c(0, y0 * 1.12)
  label_y  <- y0 * 1.04

  my_theme <- theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

  p <- ggbarplot(d, x = "Ecosystem", y = "richness",
                 add = "mean_se", color = "black", fill = "Density",
                 position = position_dodge(0.8)) +
    ylab("Richness") +
    xlab("Ecosystem") +
    stat_compare_means(aes(group = Density), label = "p.signif",
                       label.y = label_y, size = 2.1, method = "anova") +
    scale_fill_manual(values = setNames(SAT_COLORS, DENSITY_LEVELS)) +
    theme(panel.background = element_rect(fill = "transparent", colour = NA),
          panel.border = element_rect(fill = NA, colour = "black", size = 0.25),
          axis.ticks = element_line(color = "black", linewidth = 0.25),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 8, color = "black"),
          axis.title.x = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 6, color = "black",
                                     angle = 45, hjust = 1),
          legend.position = "none",
          legend.text = element_text(size = 8, color = "black")) +
    my_theme

  p <- p +
    geom_point(data = d,
               aes(x = Ecosystem, y = richness, fill = Density, group = Density),
               shape = 21, colour = "black", stroke = 0.12,
               size = 1.1, alpha = 1,
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               inherit.aes = FALSE, show.legend = FALSE)

  p <- p + ylim(y_limit)

  anova_rows <- list()
  for (ec in levels(d$Ecosystem)) {
    sub <- d[d$Ecosystem == ec, ]
    if (nrow(sub) < 3 || length(unique(sub$Density)) < 2) next
    mod <- suppressWarnings(aov(richness ~ Density, data = sub))
    s <- summary(mod)[[1]]
    anova_rows[[length(anova_rows) + 1]] <- data.frame(
      ecosystem = trimws(as.character(ec)),
      df1       = s$Df[1],
      df2       = s$Df[2],
      F_val     = s$`F value`[1],
      p_value   = s$`Pr(>F)`[1],
      n         = nrow(sub),
      stringsAsFactors = FALSE)
  }

  lsd_rows <- list()
  for (ec in levels(d$Ecosystem)) {
    sub <- d[d$Ecosystem == ec, ]
    if (nrow(sub) < 3 || length(unique(sub$Density)) < 2) next
    mod <- suppressWarnings(aov(richness ~ Density, data = sub))
    ls  <- tryCatch(agricolae::LSD.test(mod, "Density", alpha = 0.05),
                    error = function(e) NULL)
    if (is.null(ls) || is.null(ls$groups)) next
    g <- ls$groups
    lsd_rows[[length(lsd_rows) + 1]] <- data.frame(
      ecosystem  = trimws(as.character(ec)),
      Density = rownames(g),
      letter     = as.character(g$groups),
      ymax       = as.numeric(g$richness) + 0.03 * diff(range(d$richness)),
      stringsAsFactors = FALSE)
  }
  if (length(lsd_rows)) {
    lsd_df <- do.call(rbind, lsd_rows)
    lsd_df$Density <- factor(lsd_df$Density, levels = levels(d$Density))

    if (!GRID) write.csv(lsd_df, file.path(STATS, paste0("lsd_richness_", marker, ".csv")),
                         row.names = FALSE)
    p <- p + geom_text(data = lsd_df,
                       aes(x = ecosystem, y = ymax, label = letter, group = Density),
                       position = position_dodge(width = 0.8),
                       size = 2.1, vjust = -0.5, inherit.aes = FALSE)
  }

  ggsave(out_pdf, p, width = 4.0, height = 4.0, units = "cm")
  list(plot = p,
       anova = if (length(anova_rows)) do.call(rbind, anova_rows) else NULL)
}

GF_FILES <- list(
  `16S` = list(otu = file.path(DATA, "validation/prokaryote/Verify_rare_16S.csv"),
               group = file.path(DATA, "validation/prokaryote/group.csv")),
  ITS = list(otu = file.path(DATA, "validation/fungi/Verify_rare_ITS.csv"),
             group = file.path(DATA, "validation/fungi/group.csv")))
CRO_FILES <- list(
  `16S` = list(otu = file.path(DATA, "validation/prokaryote/Verify_rare_16S_Cropland.csv"),
               group = file.path(DATA, "validation/prokaryote/group_cropland.csv")),
  ITS = list(otu = file.path(DATA, "validation/fungi/Verify_rare_ITS_Cropland.csv"),
             group = file.path(DATA, "validation/fungi/group_cropland.csv")))

pcoa_res <- list(); p_pcoa <- list()
for (mk in c("16S", "ITS")) {
  for (ec in c("Grasslands", "Forests")) {
    res <- compute_verify_pcoa(GF_FILES[[mk]]$otu, GF_FILES[[mk]]$group,
                               ecos_filter = ec)
    cat(sprintf("[FigS3 %s %s] n=%d\n", mk, ec, nrow(res$site)))
    key <- paste(mk, ec)
    pcoa_res[[key]] <- res
    p_pcoa[[key]] <- plot_verify_pcoa(
      res$site, res$pc,
      file.path(OUT, sprintf("FigS3_%s_pcoa_%s.pdf", mk, tolower(ec))),
      .pcoa_own_lim(res$site), title = ec, perm = res$perm)
  }
  res <- compute_verify_pcoa(CRO_FILES[[mk]]$otu, CRO_FILES[[mk]]$group)
  cat(sprintf("[FigS3 %s Croplands] n=%d\n", mk, nrow(res$site)))
  key <- paste(mk, "Croplands")
  pcoa_res[[key]] <- res
  p_pcoa[[key]] <- plot_verify_pcoa(
    res$site, res$pc,
    file.path(OUT, sprintf("FigS3_%s_pcoa_cropland.pdf", mk)),
    .pcoa_own_lim(res$site), title = "Croplands", perm = res$perm)
}

perm_tab <- do.call(rbind, lapply(names(pcoa_res), function(k) {
  data.frame(marker    = sub(" .*", "", k),
             ecosystem = sub("^[^ ]+ ", "", k),
             n         = nrow(pcoa_res[[k]]$site),
             pcoa_res[[k]]$perm,
             stringsAsFactors = FALSE)
}))
perm_tab$q_BH <- p.adjust(perm_tab$p, method = "BH")
perm_tab <- perm_tab[, c("marker", "ecosystem", "term", "n",
                         "df1", "df2", "F_val", "R2", "p", "q_BH", "symbol")]
if (!GRID) write.csv(perm_tab, file.path(STATS, "figS3_pcoa_permanova.csv"), row.names = FALSE)
cat("\n[FigS3] PCoA PERMANOVA table (", nrow(perm_tab), " panels) -> ",
    file.path(STATS, "figS3_pcoa_permanova.csv"), "\n", sep = "")
print(perm_tab, row.names = FALSE)

rich_16S <- load_richness_data(list(GF_FILES[["16S"]], CRO_FILES[["16S"]]))
rich_ITS <- load_richness_data(list(GF_FILES[["ITS"]], CRO_FILES[["ITS"]]))
cat(sprintf("[FigS3 Richness] 16S n=%d (%s) / ITS n=%d (%s)\n",
            nrow(rich_16S), paste(table(rich_16S$Ecosystem), collapse = "/"),
            nrow(rich_ITS), paste(table(rich_ITS$Ecosystem), collapse = "/")))

r_rich_16S <- plot_verify_richness(rich_16S,
                                   file.path(OUT, "FigS3_16S_richness.pdf"), "16S")
r_rich_ITS <- plot_verify_richness(rich_ITS,
                                   file.path(OUT, "FigS3_ITS_richness.pdf"), "ITS")
p_rich_16S <- r_rich_16S$plot
p_rich_ITS <- r_rich_ITS$plot

anova_all <- list()
if (!is.null(r_rich_16S$anova)) { r_rich_16S$anova$marker <- "16S"; anova_all[[length(anova_all) + 1]] <- r_rich_16S$anova }
if (!is.null(r_rich_ITS$anova)) { r_rich_ITS$anova$marker <- "ITS"; anova_all[[length(anova_all) + 1]] <- r_rich_ITS$anova }
if (length(anova_all)) {
  anova_df <- do.call(rbind, anova_all)
  anova_df <- anova_df[, c("marker", "ecosystem", "df1", "df2", "F_val", "p_value", "n")]
  if (!GRID) write.csv(anova_df, file.path(STATS, "figS3_richness_anova.csv"), row.names = FALSE)
}

.save_legend(p_rich_16S, "FigS3_density_legend", OUT)

message("FigS3 outputs in ", OUT, ":")
print(list.files(path = OUT, pattern = "^FigS3_.*\\.pdf$", full.names = TRUE))
