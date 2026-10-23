DATA  <- "data"
OUT   <- file.path("output", "FigS13")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "FigS13")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})

genecopies = read.csv(gzfile(file.path(DATA, "metagenome/gene_copies2.csv.gz")), row.names = 1)
group      = read.csv(file.path(DATA, "metagenome/group.csv"))
KO         = read.csv(file.path(DATA, "metagenome/KO.csv"))

genecopies_sum = aggregate(genecopies$gene_copies,
                           by = list(genecopies$sample, genecopies$KO), FUN = "sum")
colnames(genecopies_sum) = c("sample","ko","gene_copies")
heat_matrix <- spread(genecopies_sum, key = "sample", value = "gene_copies", fill = 0)
rownames(heat_matrix) <- heat_matrix[,1]
input = merge(heat_matrix, KO, by.x = "ko", by.y = "KEGG_KO")
rownames(input) <- input$Gene_name

.hg <- function(grp1m, grp1sd, grp1n, grp2m, grp2sd, grp2n) {
  sp <- sqrt(((grp1n-1)*grp1sd^2 + (grp2n-1)*grp2sd^2) / (grp1n+grp2n-2))
  (grp1m - grp2m) / sp * (1 - 3/(4*(grp1n+grp2n)-9))
}

resultslis <- c()
for (cyc in c("carbon","nitrogen","phosphorus","sulfur")) {
  tempheat <- subset(input, input$cycle == cyc)
  tempheat <- tempheat[order(tempheat$class, tempheat$Gene_name), ]
  heat_data <- tempheat[, c(2:49)]
  for (row in 1:nrow(heat_data)) {
    gene <- as.data.frame(t(heat_data[row, ]))
    gene$name <- rownames(gene)
    genecv <- merge(gene, group, by.x = "name", by.y = "sample")
    for (treat in c("C","N","pH","T","W")) {
      genecv_fin <- genecv[genecv$treat %in% c("CK", treat), ]
      for (sat in c("10%","100%")) {
        f <- genecv_fin[genecv_fin$satru %in% sat, ]
        v_ck <- f[f$treat %in% "CK", 2]; v_tr <- f[f$treat %in% treat, 2]
        stopifnot(length(v_ck) == 4, length(v_tr) == 4)
        g <- .hg(sd(v_ck), mean(v_ck), 4, sd(v_tr), mean(v_tr), 4)
        resultslis <- rbind(resultslis, c(colnames(gene)[1], sat, g,
                                          cyc, tempheat$class[row], treat))
      }
    }
  }
}
colnames(resultslis) <- c("Gene_ID","satru","es","cyc","class","treat")
resultslis <- as.data.frame(resultslis, stringsAsFactors = FALSE)
resultslis$es   <- as.numeric(resultslis$es)
resultslis$type <- paste0(resultslis$treat, resultslis$satru)
write.csv(resultslis, file.path(STATS, "figS13_hedges_g_long.csv"), row.names = FALSE)

gmatrix <- spread(resultslis[, c("Gene_ID","type","es")], key = "type", value = "es", fill = 0)
gmatrix <- na.omit(gmatrix)
input2 = merge(gmatrix, KO, by.x = "Gene_ID", by.y = "Gene_name")
rownames(input2) <- input2$Gene_ID

PAL      <- colorRampPalette(c("#4D6AD0", "white", "#FF9D7E"))(50)
COL_SAT  <- c("10%" = "#FFD460", "100%" = "#2D4059")
COL_TR   <- c(C = "#8DAF89", N = "#F0A0A0", pH = "#FFB347", W = "#FFEA70", D = "#8a9cDa")
.RELBL   <- c(CK = "CK", C = "C", N = "N", pH = "pH", T = "W", W = "D")
STRIP_KO <- function(x) gsub("\\(K[0-9]+\\)$", "", x)

CELL  <- 12
GAPC  <- 4
GAPR  <- 4
BAND  <- 9

GAPS_ROW_ORIG <- list("carbon" = c(3,23), "nitrogen" = c(1,5,12,13,14,18,21),
                      "phosphorus" = 4,       "sulfur" = c(5,6))

.draw_panel <- function(cyc, out_file) {
  tempheat <- subset(input2, input2$cycle == cyc)
  tempheat <- tempheat[order(tempheat$class, tempheat$Gene_ID), ]
  hd  <- tempheat[, c(2:11)]
  cls <- tempheat$class
  ng  <- nrow(hd)

  gap_after_row <- GAPS_ROW_ORIG[[cyc]]
  gap_cls <- which(cls[-1] != cls[-ng])
  if (!identical(as.integer(sort(gap_after_row)), as.integer(gap_cls))) {
    cat(sprintf("[note] %s: class boundary (%s) differs from gaps_row (%s); gaps kept as hard-coded\n",
                cyc, paste(gap_cls, collapse = ","), paste(gap_after_row, collapse = ",")))
  }

  ypos <- -cumsum(c(0, rep(CELL, ng - 1) +
                 ifelse(seq_len(ng - 1) %in% gap_after_row, GAPR, 0))) / CELL
  xpos <- (seq_len(10) - 1 + rep(seq_len(5) - 1, each = 2) * (GAPC/CELL))

  genes <- STRIP_KO(rownames(hd))
  stopifnot(!any(duplicated(genes)))

  df <- data.frame(x = rep(xpos, each = ng), y = rep(ypos, times = 10),
                   v = as.vector(as.matrix(hd)), stringsAsFactors = FALSE)
  df$lab <- sprintf("%.2f", df$v)
  rng <- range(df$v)

  tr_of_col  <- rep(c("C","N","pH","T","W"), each = 2)
  density_of_col <- rep(c("10%","100%"), times = 5)
  band_y1 <- min(ypos) - (2.2 * BAND) / CELL; band_h <- BAND / CELL
  density_df <- data.frame(x = xpos, y = band_y1, fill = I(unname(COL_SAT[density_of_col])))
  tr_x   <- tapply(xpos, tr_of_col, range)
  tr_df  <- data.frame(xmin = sapply(tr_x, min) - 0.47, xmax = sapply(tr_x, max) + 0.47,
                       y = band_y1 - 1.15 * band_h,
                       fill = I(unname(COL_TR[.RELBL[names(tr_x)]])), stringsAsFactors = FALSE)

  x_lab  <- max(xpos) + 0.47 + 0.35
  x_cbar <- x_lab + 2.75
  x_tick <- x_cbar + 0.35
  x_end  <- x_tick + 1.25
  x_min  <- -0.7
  bar_len <- max(6, 0.30 * (-min(ypos)))
  seg <- 60
  bar_df <- data.frame(y = -(seq_len(seg) - 0.5) * bar_len / seg,
                       col = PAL[round(seq(50, 1, length.out = seg))])
  tick_v <- pretty(rng, 6)
  tick_v <- tick_v[tick_v >= min(rng) & tick_v <= max(rng)]
  tick_y <- -bar_len * (tick_v - min(rng)) / diff(range(rng))

  p <- ggplot(df, aes(x, y)) +
    geom_tile(aes(fill = v), colour = "black", linewidth = 0.2, width = 0.97, height = 0.97) +
    geom_text(aes(label = lab), size = 4/2.845, colour = "black") +
    scale_fill_gradientn(colours = PAL, limits = rng, guide = "none") +
    scale_x_continuous(limits = c(x_min, x_end), expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(ypos) - 3.0, 0.7), expand = c(0, 0)) +
    coord_fixed(ratio = 1) +

    geom_text(data = data.frame(x = x_lab, y = ypos, lab = genes),
              aes(x, y, label = lab), inherit.aes = FALSE,
              size = 6/2.845, fontface = "italic", colour = "black", hjust = 0) +

    geom_tile(data = density_df, aes(x, y, fill = fill), inherit.aes = FALSE,
              width = 0.97, height = 0.92, colour = "black", linewidth = 0.15) +

    geom_rect(data = tr_df, aes(xmin = xmin, xmax = xmax, ymin = y - band_h*0.46,
                                ymax = y + band_h*0.46, fill = fill),
              inherit.aes = FALSE, colour = "black", linewidth = 0.15) +

    annotate("text", x = x_lab, y = band_y1, label = "Density",
             hjust = 0, size = 5.5/2.845, colour = "black") +
    annotate("text", x = x_lab, y = band_y1 - 1.15*band_h, label = "Treatments",
             hjust = 0, size = 5.5/2.845, colour = "black") +

    geom_tile(data = bar_df, aes(x = x_cbar, y = y), inherit.aes = FALSE,
              fill = bar_df$col, colour = NA, width = 0.3, height = bar_len/seg * 1.05) +
    annotate("rect", xmin = x_cbar - 0.15, xmax = x_cbar + 0.15,
             ymin = -bar_len, ymax = 0,
             fill = NA, colour = "black", linewidth = 0.2) +
    geom_text(data = data.frame(x = x_tick, y = tick_y, lab = sprintf("%.1f", tick_v)),
              aes(x, y, label = lab), inherit.aes = FALSE,
              size = 6/2.845, hjust = 0, colour = "black") +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

  W <- (x_end - x_min) * CELL; H <- (0.7 - (min(ypos) - 3.0)) * CELL
  ggplot2::ggsave(out_file, p, width = W/72, height = H/72, units = "in", limitsize = FALSE)
  cat(sprintf("[out] %s  %d genes x 10 cols  (%.1f x %.1f pt)\n",
              basename(out_file), ng, W, H))
  invisible(list(df = df, hd = hd, genes = genes, cls = cls))
}

CYC_F <- c("carbon" = "carbon", "nitrogen" = "nitrogen", "phosphorus" = "phosphorus", "sulfur" = "sulfur")
for (cyc in names(CYC_F)) {
  res <- .draw_panel(cyc, file.path(OUT, paste0(cyc, "_effect_size_matrix.pdf")))
  out_df <- data.frame(gene = res$genes, class = res$cls, res$hd, check.names = FALSE)
  write.csv(out_df, file.path(OUT, sprintf("figS13_matrix_%s.csv", CYC_F[[cyc]])), row.names = FALSE)
}

pal_es <- colorRampPalette(c("#4D6AD0", "white", "#FF9D7E"))(200)
seg_h <- 100
barh <- data.frame(x = seq(0, 1, length.out = seg_h),
                   col = pal_es[round(seq(1, 200, length.out = seg_h))])
tick_v <- seq(-1, 1, 0.2)
tick_x <- (tick_v + 1) / 2
cb <- ggplot(barh, aes(x = x, y = 0)) +
  geom_tile(aes(fill = I(col)), width = 1.02/seg_h, height = 1) +
  geom_segment(data = data.frame(x = tick_x),
               aes(x = x, xend = x, y = 0.5, yend = 0.62),
               inherit.aes = FALSE, linewidth = 0.3, colour = "black") +
  geom_text(data = data.frame(x = tick_x, v = tick_v),
            aes(x = x, y = 0.78, label = sprintf("%.1f", v)),
            size = 8/2.845, colour = "black") +
  scale_x_continuous(limits = c(-0.04, 1.04), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.5, 0.95), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
ggsave(file.path(OUT, "FigS13_effectsize_colorbar.pdf"),
       cb, width = 9, height = 1.4, units = "cm")

mk_disc <- function(df, pal, name) {
  ggplot(df, aes(x = .data[[name]], y = y, fill = .data[[name]])) +
    geom_col(width = 0.9) +
    scale_fill_manual(values = pal, name = name) +
    theme_void() + theme(legend.position = "right",
                         legend.text = element_text(size = 8),
                         legend.title = element_text(size = 9, face = "bold"))
}
density_df <- data.frame(Density = factor(c("10%","100%"), levels = c("10%","100%")), y = 1)
tr_df  <- data.frame(Treatments = factor(c("C","N","pH","W","D"),
                                         levels = c("C","N","pH","W","D")), y = 1)
ggsave(file.path(OUT, "FigS13_Density_legend.pdf"),
       cowplot::ggdraw(cowplot::get_legend(mk_disc(density_df, COL_SAT, "Density"))),
       width = 4.2, height = 2.6, units = "cm")
ggsave(file.path(OUT, "FigS13_Treatments_legend.pdf"),
       cowplot::ggdraw(cowplot::get_legend(mk_disc(tr_df, COL_TR, "Treatments"))),
       width = 5.6, height = 3.4, units = "cm")

cat("DONE figS10 hedges-g panels\n")
