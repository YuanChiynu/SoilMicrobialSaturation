OUT <- file.path("output", "FigS16")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TBL <- file.path("output", "TableS8_LEfSe")
TSTATS <- file.path("output", "TableS8_LEfSe", "stats")

suppressMessages({library(ggplot2); library(gtable); library(grid)})

N_TOP      <- 5
DENS       <- c("S10", "S40", "S70", "S100")
dens_labs  <- c("10%", "40%", "70%", "100%")
eco_levels <- c("Cropland", "Forest", "Grassland")
kin_levels <- c("16S", "ITS")
kin_labs   <- c("16S" = "Prokaryotes (16S)", "ITS" = "Fungi (ITS)")
pal_dens   <- c("10%" = "#8C2D04", "40%" = "#D94801", "70%" = "#F16913", "100%" = "#FDD0A2")

CLS_LEVELS <- c("r-selected (fast)", "K-selected (slow)",
                "Stress-tolerant", "Other / unassigned")
pal_cls <- c("r-selected (fast)"   = "#D1495B",
             "K-selected (slow)"   = "#00798C",
             "Stress-tolerant"     = "#EDAE49",
             "Other / unassigned"  = "#9AA0A6")

classify_strategy <- function(strategy, is_bac) {
  strategy <- as.character(strategy)
  out <- rep(CLS_LEVELS[4], length(strategy))
  ib <- as.logical(is_bac)
  if (length(ib) == 1L) ib <- rep(ib, length(strategy))
  ib[is.na(ib)] <- FALSE
  if (any(ib)) {
    s <- strategy[ib]; o <- rep(CLS_LEVELS[4], length(s))

    o[grepl("spore-forming|desiccation|low soil moisture", s, ignore.case = TRUE)] <- CLS_LEVELS[3]
    o[grepl("oligotrophic|slow-growing|K-selected|predatory|chemolithoautotrophic|ammonia-oxidizing|methanotrophic|sulfate-reducing",
            s, ignore.case = TRUE)] <- CLS_LEVELS[2]
    o[grepl("copiotrophic|r-selected", s, ignore.case = TRUE)] <- CLS_LEVELS[1]
    out[ib] <- o
  }
  if (any(!ib)) {
    s <- strategy[!ib]; o <- rep(CLS_LEVELS[4], length(s))

    o[grepl("saprotroph|decomposer|ruderal|mycoparasite", s, ignore.case = TRUE)] <- CLS_LEVELS[1]
    o[grepl("xerotolerant|stress-tolerant|drought-associated", s, ignore.case = TRUE)] <- CLS_LEVELS[3]
    o[grepl("mycorrhizal|biotroph", s, ignore.case = TRUE)] <- CLS_LEVELS[2]
    out[!ib] <- o
  }
  out
}

cat("== FigS16 (LEfSe box plots) ==\n")

res <- read.csv(file.path(TSTATS, "tableS8_lefse_resdiff_all.csv"),
                header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("kingdom", "level", "ecosystem", "factor", "taxon", "LDA", "P.adj",
                "kw_H", "kw_df", "kw_P") %in% names(res)))
res <- res[res$factor == "density" & res$level == "genus" &
           res$kingdom %in% kin_levels, , drop = FALSE]
stopifnot(nrow(res) > 0)

pick_top <- function(kk) {
  d <- res[res$kingdom == kk & !is.na(res$P.adj) & res$P.adj < 0.05, , drop = FALSE]
  if (nrow(d) == 0) return(character(0))
  n_eco <- tapply(d$ecosystem, d$taxon, function(x) length(unique(x)))
  lda_sum <- tapply(d$LDA, d$taxon, sum)
  keep <- names(n_eco)[n_eco >= 2]
  keep <- intersect(keep, names(lda_sum))
  keep <- keep[order(-lda_sum[keep])]
  head(keep, N_TOP)
}
top_by_kin <- lapply(kin_levels, pick_top); names(top_by_kin) <- kin_levels
for (kk in kin_levels) cat(sprintf("  [%s] selected genera: %d  (%s)\n", kk, length(top_by_kin[[kk]]),
                                   paste(top_by_kin[[kk]], collapse = " | ")))
stopifnot(all(lengths(top_by_kin) > 0))

sel <- do.call(rbind, lapply(kin_levels, function(kk) {
  if (!length(top_by_kin[[kk]])) return(NULL)
  data.frame(kingdom = kk, taxon = top_by_kin[[kk]], stringsAsFactors = FALSE)
}))
sel$key <- paste(sel$kingdom, sel$taxon, sep = "\r")

sel$strategy <- vapply(seq_len(nrow(sel)), function(i) {
  d <- res[res$kingdom == sel$kingdom[i] & res$taxon == sel$taxon[i], , drop = FALSE]
  if (nrow(d) == 0) return(NA_character_)
  as.character(d$strategy[1])
}, "")
sel$strategy_class <- classify_strategy(sel$strategy, sel$kingdom == "16S")
stopifnot(!anyNA(sel$strategy_class))
for (kk in kin_levels) {
  d <- sel[sel$kingdom == kk, , drop = FALSE]
  cat(sprintf("  [%s] strategies: %s\n", kk,
              paste(sprintf("%s=%s", d$taxon, d$strategy_class), collapse = " | ")))
}

val <- read.csv(file.path(TBL, "tableS8_lefse_values_long.csv"),
                header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
val <- val[val$factor == "density" & val$level == "genus", , drop = FALSE]
val$key <- paste(val$kingdom, val$taxon, sep = "\r")
val <- val[val$key %in% sel$key, , drop = FALSE]
stopifnot(nrow(val) > 0)
val <- val[, c("kingdom", "ecosystem", "taxon", "group", "sample", "ra")]
names(val)[names(val) == "group"] <- "density"

chk <- aggregate(list(n = val$sample),
                 by = list(kingdom = val$kingdom, ecosystem = val$ecosystem,
                           taxon = val$taxon, density = val$density),
                 FUN = length)
cat(sprintf("  [check] unique within-group n: %s (expected 24)\n", paste(sort(unique(chk$n)), collapse = ",")))
stopifnot(all(chk$n == 24))

panel_levels <- as.vector(t(outer(kin_levels, eco_levels, paste, sep = "\r")))
strip_lab <- setNames(
  sprintf("%s\n%s", kin_labs[vapply(panel_levels, function(z) strsplit(z, "\r")[[1]][1], "")],
          vapply(panel_levels, function(z) strsplit(z, "\r")[[1]][2], "")),
  panel_levels)
stopifnot(!any(is.na(strip_lab)))
val$panel <- factor(paste(val$kingdom, val$ecosystem, sep = "\r"),
                    levels = panel_levels, labels = unname(strip_lab))
stopifnot(!any(is.na(val$panel)))
val$taxon <- factor(val$taxon, levels = unique(sel$taxon))
val$density <- factor(val$density, levels = DENS, labels = dens_labs)
stopifnot(!any(is.na(val$density)))
cat(sprintf("  [check] panels: %d  (%s)\n", length(unique(val$panel)),
            paste(gsub("\n", " / ", levels(val$panel)), collapse = " | ")))

val$strategy_class <- sel$strategy_class[match(paste(val$kingdom, val$taxon, sep = "\r"),
                                               sel$key)]
stopifnot(!anyNA(val$strategy_class))
track <- unique(data.frame(panel = val$panel, taxon = val$taxon,
                           strategy_class = val$strategy_class, stringsAsFactors = FALSE))
ymax <- tapply(val$ra, val$panel, max, na.rm = TRUE)
track$ypos <- as.numeric(ymax[as.character(track$panel)]) * 1.115
track$strategy_class <- factor(track$strategy_class, levels = CLS_LEVELS)
stopifnot(nrow(track) == 30L)
cat(sprintf("  [check] strategy-track rows: %d; class distribution: %s\n", nrow(track),
            paste(names(table(track$strategy_class)), table(track$strategy_class),
                  sep = "=", collapse = " ")))

p <- ggplot(val, aes(x = taxon, y = ra, fill = density)) +
  geom_boxplot(width = 0.72, position = position_dodge(width = 0.78),
               outlier.size = 0.18, linewidth = 0.22, colour = "grey20") +

  geom_point(data = track, aes(x = taxon, y = ypos, colour = strategy_class),
             shape = 15, size = 1.9, inherit.aes = FALSE) +
  scale_fill_manual(values = pal_dens, name = "Initial propagule density", drop = FALSE) +
  scale_colour_manual(values = pal_cls, name = "Life strategy", drop = TRUE) +
  facet_wrap(~ panel, nrow = 2, ncol = 3, scales = "free") +
  labs(x = NULL, y = "Relative abundance (%)") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = "grey80", linewidth = 0.25),
        strip.text = element_text(size = 6.6, face = "bold", lineheight = 1.02),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 5.6, angle = 45, hjust = 1, vjust = 1),
        axis.title = element_text(size = 8),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.box.spacing = unit(0.06, "cm"),
        legend.spacing.y = unit(0.02, "cm"),
        legend.title = element_text(size = 6.6),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.24, "cm"),
        legend.margin = margin(t = -2, b = 1),
        panel.spacing.x = unit(0.18, "cm"),
        panel.spacing.y = unit(0.34, "cm"),

        plot.margin = margin(15, 6, 6, 16))

g <- ggplotGrob(p)
pan <- g$layout[grepl("^panel", g$layout$name), c("t", "l"), drop = FALSE]
pan <- pan[order(pan$t, pan$l), , drop = FALSE]
stopifnot(nrow(pan) == 6L)
top_t <- min(pan$t)
g <- gtable_add_rows(g, unit(9, "pt"), pos = top_t - 1)
for (k in seq_len(6L)) {
  tg <- textGrob(paste0("(", letters[k], ")"),
                 x = unit(0, "npc"), hjust = 0,
                 y = unit(2, "pt"), vjust = 0, gp = gpar(fontsize = 7))
  g <- gtable_add_grob(g, tg, t = pan$t[k], l = pan$l[k],
                       clip = "off", name = paste0("panel_tag_", k))
}

pdf(file.path(OUT, "FigS16_composite.pdf"), width = 6.30, height = 6.40)
grid.newpage(); grid.draw(g); dev.off()
png(file.path(OUT, "FigS16_composite_p1.png"), width = 6.30, height = 6.40,
    units = "in", res = 300)
grid.newpage(); grid.draw(g); dev.off()
cat("  -> FigS16_composite.pdf（6.30 × 6.40 in）\n")

leg_src <- ggplot(val, aes(x = taxon, y = ra, fill = density)) +
  geom_boxplot(width = 0.7) +
  geom_point(data = track, aes(x = taxon, y = ypos, colour = strategy_class),
             shape = 15, size = 2.4, inherit.aes = FALSE) +
  scale_fill_manual(values = pal_dens, name = "Initial propagule density", drop = FALSE) +
  scale_colour_manual(values = pal_cls, name = "Life strategy", drop = TRUE) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "center",
        legend.spacing.y = unit(0.05, "cm"),
        legend.title = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.30, "cm"))
gt <- ggplotGrob(leg_src)
idx <- which(gt$layout$name == "guide-box-bottom")
stopifnot(length(idx) == 1L)
pdf(file.path(OUT, "FigS16_legend.pdf"), width = 5.6, height = 1.35)
grid.newpage(); grid.draw(gt$grobs[[idx]]); dev.off()
png(file.path(OUT, "FigS16_legend.png"), width = 5.6, height = 1.35, units = "in", res = 300)
grid.newpage(); grid.draw(gt$grobs[[idx]]); dev.off()
cat("  -> FigS16_legend.pdf\n")

outv <- val
outv$kingdom <- as.character(outv$kingdom)
outv$density <- DENS[match(as.character(outv$density), dens_labs)]
outv$ra <- round(outv$ra, 5)
outv <- outv[, c("kingdom", "ecosystem", "taxon", "strategy_class",
                 "density", "sample", "ra")]
write.csv(outv, file.path(OUT, "FigS16_lefse_boxplot_values.csv"), row.names = FALSE)
cat(sprintf("  -> FigS16_lefse_boxplot_values.csv  %d rows\n", nrow(outv)))

per_tax <- do.call(rbind, lapply(kin_levels, function(kk) {
  d <- res[res$kingdom == kk & res$taxon %in% top_by_kin[[kk]], , drop = FALSE]
  do.call(rbind, lapply(split(d, d$ecosystem), function(x) {
    x <- x[order(-x$LDA), ]
    data.frame(kingdom = kk, ecosystem = x$ecosystem[1], taxon = x$taxon,

               strategy_class = classify_strategy(x$strategy, kk == "16S"),
               strategy_level = x$strategy_level,
               LDA = round(x$LDA, 2), kw_H = round(x$kw_H, 2), kw_df = x$kw_df,
               kw_P = signif(x$kw_P, 3), kw_q = signif(x$P.adj, 3),
               group_enriched = x$group_enriched, strategy_refs = x$strategy_refs,
               stringsAsFactors = FALSE)
  }))
}))
write.csv(per_tax, file.path(STATS, "FigS16_lefse_taxon_stats.csv"), row.names = FALSE)

cap <- c(
  "Fig. S16. Genus-level dose-response to initial propagule density identified by LEfSe in three soil ecosystems.",
  "",
  sprintf("Box plots show the relative abundance (%%) of the top-ranked bacterial/archaeal 16S (top row; %s) and fungal ITS (bottom row; %s) genera along the propagule-density gradient. Columns are the three ecosystems (cropland, forest, grassland). Within each panel, boxes are grouped by genus and coloured by initial propagule density (10%%, 40%%, 70%% and 100%%; n = 24 samples per density and genus, comprising six global-change treatments x four replicates).",
          paste(top_by_kin[["16S"]], collapse = "; "),
          paste(top_by_kin[["ITS"]], collapse = "; ")),
  "Genera were selected as those identified by LEfSe (Segata et al. 2011) as differentially abundant among density levels in at least two of the three ecosystems (Benjamini-Hochberg adjusted q < 0.05), ranked by their summed linear discriminant analysis (LDA) score; the same genera are shown for all three ecosystems of a row so that the dose-response is directly comparable across ecosystems. The coloured square above each genus encodes its documented life-strategy class (see legend), so that each panel reads at two levels at once: how the abundance of a genus shifts along the density gradient (box colours) and whether that genus is a fast-growing r-selected, a slow-growing K-selected, a stress-tolerant or an otherwise unassigned taxon (square colours).",
  "Life-strategy classes were assigned by matching the literature-based annotations compiled in Table S8: r-selected (fast-growing copiotrophs in bacteria; saprotrophs, decomposers and ruderal fungi), K-selected (slow-growing oligotrophs, predators, nitrifiers, methanotrophs and sulfate reducers in bacteria; mycorrhizal symbionts and other biotrophs in fungi), stress-tolerant (spore-forming and desiccation-tolerant bacteria; xerotolerant and stress-tolerant fungi), and other/unassigned. Consistent with the strategy-level summary reported for the ANCOM-BC2 analysis, r-selected genera reach their highest relative abundance at the lowest propagule density and decline towards 100%, whereas K-selected genera increase with density.",
  "Lower initial propagule density corresponds to a larger share of sterilised soil and therefore to a higher resource surplus per propagule, whereas 100% density is competition-dominated with no resource surplus. No significance symbols are used; the Kruskal-Wallis chi-square statistic with its degrees of freedom (df = 3), the exact p value, the Benjamini-Hochberg adjusted q value, the LDA score and the pairwise Wilcoxon rank-sum sub-tests are reported for every taxon in Table S8 (LEfSe version).",
  "",
  sprintf("Generated: %s; LEfSe via microeco %s on rarefied amplicon counts; source data: Table S8 (LEfSe version).",
          format(Sys.Date(), "%Y-%m-%d"), as.character(packageVersion("microeco")))
)
writeLines(cap, file.path(OUT, "figS16_lefse_caption_draft.txt"))

cat("== done ->", file.path(OUT, "FigS16_composite.pdf"), "==\n")
cat("figure written directly by R\n")
