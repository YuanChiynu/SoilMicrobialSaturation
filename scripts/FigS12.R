DATA  <- "data"
OUT   <- file.path("output", "FigS12")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "FigS12")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

.std01 <- function(x) {
  x <- as.numeric(x); r <- range(x, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

.boot_cv_ci <- function(v, B = 2000, conf = 0.95, seed = 20260910) {
  v <- v[is.finite(v)]
  if (length(v) < 2) return(c(NA_real_, NA_real_))
  set.seed(seed)
  cvb <- vapply(seq_len(B), function(i) {
    s <- sample(v, replace = TRUE)
    m <- mean(s)
    if (!is.finite(m) || m == 0) NA_real_ else stats::sd(s) / m * 100
  }, numeric(1))
  cvb <- cvb[is.finite(cvb)]
  if (!length(cvb)) return(c(NA_real_, NA_real_))
  stats::quantile(cvb, c((1 - conf) / 2, 1 - (1 - conf) / 2), names = FALSE)
}

{
  genecopies = read.csv(gzfile(file.path(DATA, "metagenome/reads_ra.csv.gz")), row.names = 1)
  group = read.csv(file.path(DATA, "metagenome/group.csv"))
  KO = read.csv(file.path(DATA, "metagenome/KO_annotated.csv"))

  genecopies_sum = aggregate(genecopies$RA,by = list(genecopies$sample,genecopies$KO),FUN = "mean")
  colnames(genecopies_sum) = c("sample","ko","gene_copies")

  heat_matrix <- spread(genecopies_sum,key = "sample",value = "gene_copies",fill = 0)
  rownames(heat_matrix) <- heat_matrix[,1]
  input = merge(heat_matrix,KO,by.x = "ko",by.y ="KEGG_KO" )

  rownames(input) <- make.unique(input$Gene_name)

  samples_all <- as.character(group$sample)
  stopifnot(length(samples_all) == 48,
            all(samples_all %in% colnames(input)))
  dens_vec <- as.character(group$satru)
  trt_vec  <- as.character(group$treat)
  stopifnot(length(unique(dens_vec)) == 2, length(unique(trt_vec)) == 6)

  CYCLES <- c("carbon","nitrogen","phosphorus","sulfur")

  cv_values_rows <- list()
  for (cyc in CYCLES) {
    tempheat <- subset(input,input$cycle==cyc)
    tempheat <- tempheat[order(tempheat$class,tempheat$Gene_name),]
    heat_data <- tempheat[, samples_all]
    samples_10  <- group$sample[group$satru == "10%"]
    samples_100 <- group$sample[group$satru == "100%"]
    for (row in 1:nrow(heat_data)) {
      gene <- rownames(heat_data)[row]
      v10  <- as.numeric(heat_data[row, samples_10])
      v100 <- as.numeric(heat_data[row, samples_100])

      m10 <- mean(v10); s10 <- sd(v10); cv10 <- (s10 / m10) * 100
      m100 <- mean(v100); s100 <- sd(v100); cv100 <- (s100 / m100) * 100

      ci10  <- .boot_cv_ci(v10,  B = 2000)
      ci100 <- .boot_cv_ci(v100, B = 2000)

      p_fm <- NA_real_; d_fm <- NA_real_
      if (requireNamespace("cvequality", quietly = TRUE) &&
          is.finite(s10) && is.finite(s100) && s10 > 0 && s100 > 0) {
        tst <- tryCatch(cvequality::asymptotic_test2(
                          k = 2, n = c(length(v10), length(v100)),
                          s = c(s10, s100), x = c(m10, m100)),
                        error = function(e) NULL)
        if (!is.null(tst)) { d_fm <- tst$D_AD; p_fm <- tst$p_value }
      }

      cv_values_rows[[length(cv_values_rows) + 1]] <- data.frame(
        gene_name = gene, cycle = cyc, class = tempheat$class[row],
        satru = "10%",  n = length(v10),  mean = m10,  sd = s10,
        genecv = cv10,  cv_lo = ci10[1],  cv_hi = ci10[2],
        FeltzMiller_D = d_fm, p_FM = p_fm, stringsAsFactors = FALSE)
      cv_values_rows[[length(cv_values_rows) + 1]] <- data.frame(
        gene_name = gene, cycle = cyc, class = tempheat$class[row],
        satru = "100%", n = length(v100), mean = m100, sd = s100,
        genecv = cv100, cv_lo = ci100[1], cv_hi = ci100[2],
        FeltzMiller_D = d_fm, p_FM = p_fm, stringsAsFactors = FALSE)
    }
  }
  cv_values <- do.call(rbind, cv_values_rows)
  write.csv(cv_values, file.path(STATS, "figS12_gene_cv_values.csv"), row.names = FALSE)

  test_df <- cv_values[cv_values$satru == "10%",
                       c("gene_name","cycle","class","FeltzMiller_D","p_FM")]
  cv10_df <- cv_values[cv_values$satru == "10%", c("gene_name","genecv")]
  colnames(cv10_df)[colnames(cv10_df) == "genecv"] <- "CV_10"
  cv100_df <- cv_values[cv_values$satru == "100%", c("gene_name","genecv")]
  colnames(cv100_df)[colnames(cv100_df) == "genecv"] <- "CV_100"
  test_df <- merge(test_df, cv10_df, by = "gene_name")
  test_df <- merge(test_df, cv100_df, by = "gene_name")
  test_df$df <- 1
  test_df$q_BH_cycle <- NA_real_
  for (cyc in unique(test_df$cycle)) {
    idx <- test_df$cycle == cyc
    test_df$q_BH_cycle[idx] <- p.adjust(test_df$p_FM[idx], method = "BH")
  }
  test_df$q_BH_all <- p.adjust(test_df$p_FM, method = "BH")
  test_df <- test_df[, c("gene_name","cycle","class","df",
                         "CV_10","CV_100","FeltzMiller_D","p_FM",
                         "q_BH_cycle","q_BH_all")]
  write.csv(test_df, file.path(STATS, "figS12_cv_test.csv"), row.names = FALSE)
  cat("[legacy] CV csv written (still referenced by Table S7; no CV panel in the figure):", nrow(cv_values), "rows /\n")
  cat("         ", nrow(test_df), "gene-level test rows\n")

  cell_rows <- list(); samp_rows <- list(); test_rows <- list()
  pt_rows   <- list()
  const_genes <- character(0)
  for (cyc in CYCLES) {
    tempheat <- subset(input, input$cycle == cyc)
    tempheat <- tempheat[order(tempheat$class, tempheat$Gene_name), ]
    M <- as.matrix(tempheat[, samples_all])
    stopifnot(identical(colnames(M), samples_all))

    Z <- t(apply(M, 1, .std01))
    Z <- matrix(Z, nrow = nrow(M), ncol = length(samples_all),
                dimnames = list(rownames(M), samples_all))

    for (gi in seq_len(nrow(Z))) {
      gene <- rownames(Z)[gi]
      v <- as.numeric(Z[gi, ])

      if (!is.finite(diff(range(v))) || diff(range(v)) == 0) {
        const_genes <- c(const_genes, paste0(cyc, "/", gene))
        next
      }
      d_samp <- NULL
      for (gp in sort(unique(dens_vec))) {
        idx <- which(dens_vec == gp)
        ns <- length(idx)
        if (ns < 2) next
        tv <- trt_vec[idx]
        Dm <- as.matrix(stats::dist(v[idx], method = "euclidean"))

        mx <- outer(tv, tv, FUN = "!=")
        diag(mx) <- FALSE
        dd <- Dm[mx]
        nu <- length(dd) / 2

        cell_rows[[length(cell_rows) + 1]] <- data.frame(
          cycle = cyc, gene_name = gene, class = tempheat$class[gi],
          density = gp,
          n = ns, n_comp = length(dd), n_pairs = nu,
          mean_dist = mean(dd), sd_dist = stats::sd(dd),
          stringsAsFactors = FALSE)

        per <- rowSums(Dm * mx) / pmax(rowSums(mx), 1)
        samp_rows[[length(samp_rows) + 1]] <- data.frame(
          cycle = cyc, gene_name = gene, density = gp,
          sample = samples_all[idx], dist = as.numeric(per),
          stringsAsFactors = FALSE)
        d_samp <- rbind(d_samp, data.frame(
          density = gp, treat = tv, sample = samples_all[idx],
          dist = as.numeric(per), stringsAsFactors = FALSE))
      }

      d_samp$density <- factor(d_samp$density, levels = sort(unique(dens_vec)))
      d_samp$treat   <- factor(d_samp$treat,   levels = sort(unique(trt_vec)))
      stopifnot(nrow(d_samp) == 48)
      tst <- tryCatch(suppressWarnings(aov(dist ~ density * treat, data = d_samp)),
                      error = function(e) NULL)
      if (is.null(tst) || !is.finite(stats::sd(d_samp$dist)) ||
          stats::sd(d_samp$dist) == 0) {
        test_rows[[length(test_rows) + 1]] <- data.frame(
          cycle = cyc, gene_name = gene, class = tempheat$class[gi],
          response = "sample mean pairwise Euclidean distance (per-gene 0-1 standardised)",
          term = c("density","treat","density:treat"),
          term_label = c("RS (Propagule density)","T (Treatment)","RS x T (interaction)"),
          df1 = c(1,5,5), df2 = 36, F_val = NA_real_, p_value = NA_real_,
          partial_eta2 = NA_real_, stringsAsFactors = FALSE)

        pt_rows[[length(pt_rows) + 1]] <- data.frame(
          cycle = cyc, gene_name = gene, class = tempheat$class[gi],
          mean_10 = NA_real_, mean_100 = NA_real_,
          n_10 = sum(d_samp$density == "10%"), n_100 = sum(d_samp$density == "100%"),
          diff = NA_real_, t_value = NA_real_, df = NA_real_, p_welch = NA_real_,
          stringsAsFactors = FALSE)
        next
      }

      tt <- tryCatch(suppressWarnings(stats::t.test(dist ~ density, data = d_samp,
                                                    var.equal = FALSE)),
                     error = function(e) NULL)
      if (is.null(tt)) {
        pt_rows[[length(pt_rows) + 1]] <- data.frame(
          cycle = cyc, gene_name = gene, class = tempheat$class[gi],
          mean_10 = NA_real_, mean_100 = NA_real_,
          n_10 = sum(d_samp$density == "10%"), n_100 = sum(d_samp$density == "100%"),
          diff = NA_real_, t_value = NA_real_, df = NA_real_, p_welch = NA_real_,
          stringsAsFactors = FALSE)
      } else {
        .m <- stats::setNames(as.numeric(tt$estimate),
                              c("10%","100%"))
        pt_rows[[length(pt_rows) + 1]] <- data.frame(
          cycle = cyc, gene_name = gene, class = tempheat$class[gi],
          mean_10 = unname(.m["10%"]), mean_100 = unname(.m["100%"]),
          n_10 = sum(d_samp$density == "10%"), n_100 = sum(d_samp$density == "100%"),
          diff = unname(.m["10%"] - .m["100%"]),
          t_value = unname(tt$statistic), df = unname(tt$parameter),
          p_welch = unname(tt$p.value), stringsAsFactors = FALSE)
      }
      s <- summary(tst)[[1]]; ss <- as.data.frame(s); ss$term <- trimws(rownames(s))
      resid <- ss[ss$term == "Residuals", ]; keep <- ss$term != "Residuals"
      .tlab <- ifelse(ss$term[keep] == "density", "RS (Propagule density)",
               ifelse(ss$term[keep] == "treat", "T (Treatment)",
                                               "RS x T (interaction)"))
      test_rows[[length(test_rows) + 1]] <- data.frame(
        cycle = cyc, gene_name = gene, class = tempheat$class[gi],
        response = "sample mean pairwise Euclidean distance (per-gene 0-1 standardised)",
        term = ss$term[keep], term_label = .tlab,
        df1 = ss$Df[keep], df2 = resid$Df[1],
        F_val = ss$`F value`[keep], p_value = ss$`Pr(>F)`[keep],
        partial_eta2 = ss$`Sum Sq`[keep] / (ss$`Sum Sq`[keep] + resid$`Sum Sq`[1]),
        stringsAsFactors = FALSE)
    }
  }
  eu_cell  <- do.call(rbind, cell_rows)
  eu_samp  <- do.call(rbind, samp_rows)
  eu_test  <- do.call(rbind, test_rows)
  eu_pt    <- do.call(rbind, pt_rows)

  eu_pt$q_BH_cycle <- NA_real_
  for (cc in unique(eu_pt$cycle)) {
    ix <- eu_pt$cycle == cc
    eu_pt$q_BH_cycle[ix] <- stats::p.adjust(eu_pt$p_welch[ix], method = "BH")
  }
  eu_pt$q_BH_all <- stats::p.adjust(eu_pt$p_welch, method = "BH")
  eu_pt$signif <- ifelse(is.na(eu_pt$p_welch), "ns",
                         ifelse(eu_pt$p_welch < 0.05, "*", "ns"))
  eu_pt <- eu_pt[, c("cycle","gene_name","class","mean_10","mean_100","n_10","n_100",
                     "diff","t_value","df","p_welch","q_BH_cycle","q_BH_all","signif")]
  if (length(const_genes)) {
    cat("[exclude] genes constant across all samples (zero after standardisation; excluded):",
        paste(const_genes, collapse = ", "), "\n")
  }

  .k12  <- function(d) paste(d$cycle, d$gene_name, d$density, sep = "|")
  .ns12 <- tapply(seq_len(nrow(eu_samp)), .k12(eu_samp), length)
  eu_dens <- do.call(rbind, lapply(split(eu_cell, .k12(eu_cell)), function(cc) {
    k <- .k12(cc)[1]
    data.frame(
      cycle = cc$cycle[1], gene_name = cc$gene_name[1], class = cc$class[1],
      density = cc$density[1],
      mean_dist = cc$mean_dist,
      se_dist   = cc$sd_dist / sqrt(cc$n_pairs),
      sd_dist   = cc$sd_dist,
      n_samp    = unname(.ns12[k]), stringsAsFactors = FALSE)
  }))

  write.csv(eu_cell, file.path(CACHE, "figS12_euclid_cell.csv"),   row.names = FALSE)
  write.csv(eu_samp, file.path(CACHE, "figS12_euclid_sample.csv"), row.names = FALSE)
  write.csv(eu_dens, file.path(STATS, "figS12_euclid_density.csv"),row.names = FALSE)
  write.csv(eu_test, file.path(STATS, "figS12_euclid_test.csv"),   row.names = FALSE)
  write.csv(eu_pt,   file.path(STATS, "figS12_euclid_density_pairtest.csv"),
            row.names = FALSE)

  stopifnot(nrow(eu_pt) == nrow(eu_test) / 3,
            all(eu_pt$n_10 == 24 & eu_pt$n_100 == 24, na.rm = TRUE),
            all(eu_pt$df > 22 & eu_pt$df < 47 | is.na(eu_pt$df)),
            all(eu_pt$p_welch >= 0 & eu_pt$p_welch <= 1 | is.na(eu_pt$p_welch)),
            all(!is.na(eu_pt$q_BH_cycle)), all(!is.na(eu_pt$q_BH_all)))
  cat("[euclid] cell:", nrow(eu_cell), "rows | sample:", nrow(eu_samp),
      "| density:", nrow(eu_dens), "| test:", nrow(eu_test), "rows\n")
  cat("        per-gene 10%-100% Welch t:", nrow(eu_pt), "genes | p<0.05:",
      sum(eu_pt$p_welch < 0.05, na.rm = TRUE), "| q_BH_cycle<0.05:",
      sum(eu_pt$q_BH_cycle < 0.05, na.rm = TRUE), "\n")
  cat("        RS (density) significant genes: ",
      sum(eu_test$term == "density"  & eu_test$p_value < 0.05, na.rm = TRUE), "/",
      sum(eu_test$term == "density"), "\n")

  my_theme <- theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))
  .legend_done_eu <- FALSE
  for (cyc in CYCLES) {
    pdat <- subset(eu_dens, cycle == cyc)

    pdat$gene_name <- sub("\\.[0-9]+$", "", as.character(pdat$gene_name))
    pdat$density <- factor(pdat$density, levels = c("10%","100%"))

    dtop <- max(pdat$mean_dist + pdat$se_dist, na.rm = TRUE)
    if (!is.finite(dtop) || dtop <= 0) dtop <- 1

    p_sub <- ggbarplot(pdat, x="gene_name", y="mean_dist", fill = "density",
                       palette = c("#ffd460","#2d4059"),
                       position = position_dodge(0.9), width = 0.35)+
      xlab(NULL)+
      ylab("Mean pairwise Euclidean distance")+
      labs(fill = "Propagule density(%)")+
      coord_flip()+
      ylim(0,(dtop*1.2))+

      geom_errorbar(aes(ymin = mean_dist - se_dist, ymax = mean_dist + se_dist,
                        group = density), width = 0.35,
                    linewidth = 0.35, position = position_dodge(width = 0.9))+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+

      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.text.y = element_text(size = 6,color="black"))+
      theme(axis.title.y = element_text(size = 8, color="black"))+
      theme(axis.title.x = element_text(size = 8, color="black"))+
      theme(legend.text = element_text(size = 8, color="black"))+
      theme(legend.title = element_text(size = 8, color="black"))+
      theme(legend.position ="none")+
      theme(axis.text.x = element_text(size = 6,color="black",vjust = 1.0, hjust = 1.0))+
      my_theme

    pt_cyc <- subset(eu_pt, cycle == cyc)
    pt_cyc$gene_name <- sub("\\.[0-9]+$", "", as.character(pt_cyc$gene_name))
    pt_cyc$sym <- ifelse(is.na(pt_cyc$p_welch), "",
                  ifelse(pt_cyc$p_welch < 0.001, "***",
                  ifelse(pt_cyc$p_welch < 0.01, "**",
                  ifelse(pt_cyc$p_welch < 0.05, "*", ""))))
    star_df <- pt_cyc[pt_cyc$sym != "", c("gene_name", "sym"), drop = FALSE]
    if (nrow(star_df)) {
      gtop <- tapply(pdat$mean_dist + pdat$se_dist, pdat$gene_name,
                     max, na.rm = TRUE)
      star_df$star_y <- unname(gtop[star_df$gene_name]) + 0.035 * dtop
      stopifnot(all(is.finite(star_df$star_y)),
                all(star_df$star_y < dtop * 1.2),
                !anyDuplicated(star_df$gene_name))
      p_sub <- p_sub + geom_text(data = star_df,
                                 aes(x = gene_name, y = star_y, label = sym),
                                 size = 1.6, color = "black",
                                 inherit.aes = FALSE)
      cat("        [", cyc, "] significant genes: ", nrow(star_df), "/",
          length(unique(pdat$gene_name)), "\n", sep = "")
    }

    .h8 <- c("carbon" = 18.5, "nitrogen" = 9.6, "phosphorus" = 4.6, "sulfur" = 6.6)[cyc]
    ggsave(paste0(cyc, "gene_euclidean_distance.pdf"), p_sub, width = 6, height = .h8,
           path = OUT, units = "cm", useDingbats = FALSE)
    if (!.legend_done_eu) {
      .save_legend(p_sub, "FigS12_euclid_legend", OUT)
      .legend_done_eu <- TRUE
    }
  }
}

message("FigS12 outputs in ", OUT, ":")
print(list.files(path = OUT, pattern = "\\.pdf$", full.names = TRUE))
message("Euclidean-distance table -> ",
        file.path(CACHE, "figS12_euclid_test.csv"), " and 3 more csv files")
message("legacy CV csv (referenced by Table S7) -> ",
        file.path(CACHE, "figS12_gene_cv_values.csv"), " / ",
        file.path(CACHE, "figS12_cv_test.csv"))
