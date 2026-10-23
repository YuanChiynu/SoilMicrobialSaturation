DATA  <- "data"
OUT   <- file.path("output", "Fig1")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "Fig1")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

DENSITY_LEVELS   <- c("10re", "40re", "70re", "100re")
DENSITY_LABELS   <- c("10%", "40%", "70%", "100%")
TREAT_LEVELS <- c("CK", "C", "N", "pH", "T", "W")
TREAT_LABELS <- c("CK", "C", "N", "pH", "W", "D")

.disp_sat <- function(x) {
  x <- as.character(x)
  m <- match(x, DENSITY_LEVELS)
  ifelse(is.na(m), x, DENSITY_LABELS[m])
}
.disp_treat <- function(x) {

  x <- as.character(x)
  m <- match(x, TREAT_LEVELS)
  ifelse(is.na(m), x, TREAT_LABELS[m])
}

read_alpha <- function(kingdom) {
  if (kingdom == "16S") {
    group <- read.csv(file.path(DATA, "amplicon/prokaryote/group.csv"), header = TRUE)
    otu   <- read.csv(file.path(DATA, "amplicon/prokaryote/OTU_rare_16S_all.csv"),
                      header = TRUE, row.names = 1)

    eco_map <- c(C = "Cropland", F = "Forest", G = "Grassland")
    group$ecosystem <- unname(eco_map[group$ecosystem])
  } else if (kingdom == "ITS") {
    group <- read.csv(file.path(DATA, "amplicon/fungi/group.csv"), header = TRUE)
    group <- group[, c(1, 3:5)]
    otu   <- read.csv(file.path(DATA, "amplicon/fungi/OTU_rare_ITS.csv"),
                      header = TRUE, row.names = 1)

    group$ecosystem <- ifelse(group$ecosystem == "Corpland", "Cropland",
                              as.character(group$ecosystem))
    group$ecosystem <- ifelse(group$ecosystem == "Grass", "Grassland",
                              as.character(group$ecosystem))
  } else {
    stop("kingdom must be '16S' or 'ITS'")
  }
  rare <- t(otu)
  richness <- specnumber(rare)
  shannon  <- diversity(rare, index = "shannon")
  alpha <- data.frame(sample = rownames(rare), richness = richness, shannon = shannon,
                      stringsAsFactors = FALSE)
  out <- merge(alpha, group, by = "sample")
  out$ecosystem <- factor(out$ecosystem,
                          levels = c("Cropland", "Forest", "Grassland"))
  out$group     <- factor(out$group,
                          levels = DENSITY_LEVELS, labels = DENSITY_LABELS)
  out$group1    <- factor(out$group1,
                          levels = TREAT_LEVELS, labels = TREAT_LABELS)
  out$kingdom   <- kingdom
  out
}

compute_letters <- function(df, value_col) {
  rows <- list()
  for (ec in unique(as.character(df$ecosystem))) {
    for (gp in unique(as.character(df$group))) {
      sub <- df[df$ecosystem == ec & df$group == gp, ]
      if (nrow(sub) < 2 || length(unique(sub$group1)) < 2) next
      mod <- suppressWarnings(aov(stats::reformulate("group1", value_col),
                                  data = sub))
      lsd <- tryCatch(agricolae::LSD.test(mod, "group1", alpha = 0.05),
                      error = function(e) NULL)
      if (is.null(lsd)) next
      grp_tab <- lsd$groups
      for (trt in rownames(grp_tab)) {
        sub_trt <- sub[sub$group1 == trt, value_col]
        ymax    <- mean(sub_trt) + sd(sub_trt) / sqrt(length(sub_trt))
        rows[[length(rows) + 1]] <- data.frame(
          ecosystem = ec, group = gp, group1 = trt,
          letter = grp_tab[trt, "groups"],
          ymax   = ymax,
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, rows)
  out$group <- factor(out$group, levels = levels(df$group))
  out
}

compute_anova <- function(df, value_col) {
  rows <- list()
  for (ec in unique(as.character(df$ecosystem))) {
    sub <- df[df$ecosystem == ec, ]
    if (length(unique(sub$group)) < 2 || length(unique(sub$group1)) < 2) next
    mod <- suppressWarnings(aov(stats::reformulate("group*group1", value_col),
                                data = sub))
    s <- tryCatch(summary(mod)[[1]], error = function(e) NULL)
    if (is.null(s) || nrow(s) < 4) next
    ss <- as.data.frame(s)

    ss$term <- trimws(rownames(s))
    resid   <- ss[ss$term == "Residuals", ]
    keep    <- ss$term != "Residuals"
    rows[[length(rows) + 1]] <- data.frame(
      ecosystem    = ec,
      term         = ss$term[keep],
      df1          = ss$Df[keep],
      df2          = resid$Df[1],
      F_val        = ss$`F value`[keep],
      p_value      = ss$`Pr(>F)`[keep],
      partial_eta2 = ss$`Sum Sq`[keep] / (ss$`Sum Sq`[keep] + resid$`Sum Sq`[1]),
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

compute_cv_with_test <- function(df, value_col) {
  cv_tab <- list()
  for (ec in levels(df$ecosystem)) {
    for (gp in levels(df$group)) {
      for (trt in levels(df$group1)) {
        v <- df[df$ecosystem == ec & df$group == gp & df$group1 == trt, value_col]
        v <- v[!is.na(v)]
        if (length(v) < 2 || mean(v) == 0) next
        cv_tab[[length(cv_tab) + 1]] <- data.frame(
          ecosystem = ec, group = gp, group1 = trt,
          cv = sd(v) / mean(v) * 100, n = length(v),
          stringsAsFactors = FALSE)
      }
    }
  }
  cv_tab <- do.call(rbind, cv_tab)

  cv_sum <- do.call(rbind, lapply(
    split(cv_tab, paste(cv_tab$ecosystem, cv_tab$group, sep = "|")),
    function(d) data.frame(
      ecosystem = d$ecosystem[1], group = d$group[1],
      mean_cv = mean(d$cv),
      se_cv = if (nrow(d) >= 2) sd(d$cv) / sqrt(nrow(d)) else NA_real_,
      sd_cv = if (nrow(d) >= 2) sd(d$cv) else NA_real_,
      n_trt = nrow(d), stringsAsFactors = FALSE)))
  fm_tab <- list()
  if (requireNamespace("cvequality", quietly = TRUE)) {
    for (ec in levels(df$ecosystem)) {
      groups <- levels(df$group)
      n_gp <- sapply(groups, function(gp)
        sum(cv_tab$ecosystem == ec & cv_tab$group == gp))
      s_gp <- sapply(groups, function(gp) {
        d <- cv_tab[cv_tab$ecosystem == ec & cv_tab$group == gp, ]
        if (nrow(d) >= 2) sd(d$cv) else NA_real_ })
      x_gp <- sapply(groups, function(gp) {
        d <- cv_tab[cv_tab$ecosystem == ec & cv_tab$group == gp, ]
        if (nrow(d) >= 1) mean(d$cv) else NA_real_ })
      if (length(n_gp) >= 2 && all(n_gp >= 2) &&
          !any(is.na(s_gp)) && !any(is.na(x_gp)) && all(s_gp > 0)) {
        tst <- tryCatch(cvequality::asymptotic_test2(
                          k = length(groups), n = n_gp, s = s_gp, x = x_gp),
                        error = function(e) NULL)
        if (!is.null(tst)) {
          fm_tab[[length(fm_tab) + 1]] <- data.frame(
            ecosystem = ec, df = length(groups) - 1,
            D = tst$D_AD, p_FM = tst$p_value,
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (length(fm_tab) == 0) {
    fm_table <- data.frame(ecosystem = character(0), df = integer(0),
                           D = numeric(0), p_FM = numeric(0),
                           stringsAsFactors = FALSE)
  } else {
    fm_table <- do.call(rbind, fm_tab)
  }
  list(cv_table = cv_tab, cv_summary = cv_sum, fm_table = fm_table)
}

.write_cv_legacy <- function(df, kingdom) {
  for (pc in c("richness", "shannon")) {
    res    <- compute_cv_with_test(df, pc)
    cv_tab <- res$cv_table
    cv_sum <- res$cv_summary
    fm_tab <- res$fm_table
    write.csv(cv_tab,
              file.path(STATS, paste0("cv_table_", kingdom, "_", pc, ".csv")),
              row.names = FALSE)
    eco_levels <- levels(df$ecosystem)
    cv_w <- do.call(rbind, lapply(eco_levels, function(ec) {
      d <- cv_sum[cv_sum$ecosystem == ec, ]
      data.frame(ecosystem = ec,
                 CV_10    = d$mean_cv[d$group == "10%"],
                 CV_40    = d$mean_cv[d$group == "40%"],
                 CV_70    = d$mean_cv[d$group == "70%"],
                 CV_100   = d$mean_cv[d$group == "100%"],
                 CV_10_SE  = d$se_cv[d$group == "10%"],
                 CV_40_SE  = d$se_cv[d$group == "40%"],
                 CV_70_SE  = d$se_cv[d$group == "70%"],
                 CV_100_SE = d$se_cv[d$group == "100%"],
                 n_trt     = d$n_trt[d$group == "10%"],
                 stringsAsFactors = FALSE)
    }))
    fm_w <- fm_tab
    colnames(fm_w)[colnames(fm_w) == "D"] <- "FeltzMiller_D"
    combined <- merge(cv_w, fm_w, by = "ecosystem", all.x = TRUE)
    combined$kingdom <- kingdom
    combined$metric  <- pc
    combined <- combined[, c("kingdom", "ecosystem", "metric",
                             "CV_10", "CV_10_SE", "CV_40", "CV_40_SE",
                             "CV_70", "CV_70_SE", "CV_100", "CV_100_SE",
                             "FeltzMiller_D", "df", "p_FM")]
    write.csv(combined,
              file.path(STATS, paste0("cv_bootstrap_", kingdom, "_", pc, ".csv")),
              row.names = FALSE)
  }
}

compute_density_emmeans <- function(df, value_col) {
  if (!requireNamespace("emmeans", quietly = TRUE))
    stop("emmeans package required: install.packages('emmeans')")
  kd <- as.character(df$kingdom[1])
  em_rows <- list(); te_rows <- list()
  for (ec in levels(df$ecosystem)) {
    sub <- df[df$ecosystem == ec, ]
    if (nrow(sub) < 4) next

    sub$group  <- factor(sub$group,  levels = levels(df$group))
    sub$group1 <- factor(sub$group1, levels = levels(df$group1))
    mod <- suppressWarnings(aov(stats::reformulate("group*group1", value_col),
                                data = sub))
    s  <- summary(mod)[[1]]
    ss <- as.data.frame(s)
    ss$term <- trimws(rownames(s))
    ss$p_value <- ss[["Pr(>F)"]]
    ss$partial_eta2 <- ss[["Sum Sq"]] / (ss[["Sum Sq"]] + ss[["Sum Sq"]][nrow(ss)])
    keep <- ss$term != "Residuals"
    te_rows[[length(te_rows) + 1]] <- data.frame(
      kingdom = kd, ecosystem = ec,
      term = ss$term[keep], df1 = ss$Df[keep], df2 = ss$Df[!keep],
      F_val = ss[["F value"]][keep], p_value = ss$p_value[keep],
      partial_eta2 = ss$partial_eta2[keep],
      stringsAsFactors = FALSE)

    em <- tryCatch(
      as.data.frame(summary(emmeans::emmeans(mod, ~ group),
                            level = 0.95, infer = c(TRUE, TRUE))),
      error = function(e) {
        message("[1C] emmeans failed: ", kd, " / ", ec, " - ", conditionMessage(e))
        NULL
      })
    if (is.null(em)) next

    lsd <- tryCatch(agricolae::LSD.test(mod, "group", alpha = 0.05),
                    error = function(e) NULL)
    let <- rep(NA_character_, nrow(em))
    if (!is.null(lsd)) {
      gt  <- lsd$groups
      key <- trimws(rownames(gt))
      let <- unname(as.character(gt[["groups"]])[
               match(as.character(em$group), key)])
    }
    em_rows[[length(em_rows) + 1]] <- data.frame(
      kingdom = kd, ecosystem = ec,
      density = as.character(em$group),
      emmean = em$emmean, se = em$SE, df = em$df,
      lower = em$lower.CL, upper = em$upper.CL,
      letter = let,
      n = as.integer(table(sub$group)[as.character(em$group)]),
      stringsAsFactors = FALSE)
  }
  if (length(em_rows) == 0) return(NULL)
  em <- do.call(rbind, em_rows)
  te <- do.call(rbind, te_rows)
  em$ecosystem  <- factor(em$ecosystem,  levels = c("Cropland", "Forest", "Grassland"))
  te$ecosystem  <- factor(te$ecosystem,  levels = c("Cropland", "Forest", "Grassland"))
  em$density <- factor(em$density, levels = DENSITY_LABELS)
  list(emmeans = em, tests = te)
}

{
  library(permute); library(vegan); library(ggplot2); library(ggpubr)
  library(agricolae); library(stats)

  .legend_done_treat <- FALSE
  .legend_done_satru <- FALSE

  for (kingdom in c("ITS", "16S")) {
    df <- read_alpha(kingdom)
    write.csv(df, file.path(CACHE, paste0("alpha_long_", kingdom, ".csv")),
              row.names = FALSE)

    palette_treat <- c("#dddddd", "#8DaF89", "#F0a0a0", "#FFB347", "#FFEA70", "#8a9cDa")
    palette_density <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")

    for (metric in c("richness", "shannon")) {
      letters_df <- compute_letters(df, metric)
      anova_df   <- compute_anova(df, metric)
      if (!is.null(letters_df))
        write.csv(letters_df,
                  file.path(STATS, paste0("letters_lsd_", kingdom, "_", metric, ".csv")),
                  row.names = FALSE)
      if (!is.null(anova_df))
        write.csv(anova_df,
                  file.path(STATS, paste0("anova_overall_", kingdom, "_", metric, ".csv")),
                  row.names = FALSE)

      for (ec in unique(as.character(df$ecosystem))) {
        sub        <- df[df$ecosystem == ec, ]
        if (nrow(sub) == 0) next
        sub_letter <- if (!is.null(letters_df))
          letters_df[letters_df$ecosystem == ec, ] else NULL
        an <- if (!is.null(anova_df))
          anova_df[anova_df$ecosystem == ec, ] else NULL

        sig_d <- if (!is.null(an) && nrow(an) >= 1) {
          r <- an$p_value[trimws(an$term) == "group"]
          if (length(r) == 0 || is.na(r)) 1 else r[1]
        } else 1
        sig_d <- symnum(sig_d,
                        cutpoints = c(0, 1e-3, 1e-2, 0.05, 1),
                        symbols = c("***", "**", "*", "ns"))
        sig_t <- if (!is.null(an) && nrow(an) >= 1) {
          r <- an$p_value[trimws(an$term) == "group1"]
          if (length(r) == 0 || is.na(r)) 1 else r[1]
        } else 1
        sig_t <- symnum(sig_t,
                        cutpoints = c(0, 1e-3, 1e-2, 0.05, 1),
                        symbols = c("***", "**", "*", "ns"))

        p <- ggbarplot(sub, x = "group", y = metric,
                       add = "mean_se",
                       color = "black",
                       fill = "group1",
                       palette = palette_treat,
                       position = position_dodge(0.8)) +
          geom_point(data = sub, aes(x = group, y = .data[[metric]],
                                     fill = group1, group = group1),
                     shape = 21, colour = "black", stroke = 0.12,
                     size = 0.8, alpha = 0.95,
                     position = position_jitterdodge(jitter.width = 0.09,
                                                     dodge.width = 0.8),
                     inherit.aes = FALSE, show.legend = FALSE) +

          scale_fill_manual(values = setNames(palette_treat, levels(sub$group1))) +

          { if (!is.null(sub_letter) && nrow(sub_letter) > 0)
              geom_text(data = transform(sub_letter,
                                         group  = factor(group,  levels = levels(sub$group)),
                                         group1 = factor(group1, levels = levels(sub$group1))),
                        aes(x = group, y = ymax, label = letter, group = group1),
                        position = position_dodge(width = 0.8),
                        size = 2.1, vjust = -1.1, inherit.aes = FALSE) } +
          ylab(metric) +

          xlab("Propagule density(%)") +
          scale_y_continuous(
            expand = expansion(mult = c(0, 0.3))) +

          labs(title = bquote(RS ~ .(sig_d) ~~ T ~ .(sig_t)),
               fill = "Treatments") +
          theme(plot.title = element_text(size = 7, hjust = 0)) +
          theme(panel.background = element_rect(fill = "transparent", colour = "black")) +
          theme(axis.text.y = element_text(size = 6, color = "black")) +
          theme(axis.title.y = element_text(size = 8, color = "black")) +
          theme(axis.title.x = element_text(size = 8, color = "black")) +

          theme(axis.text.x = element_text(size = 6, color = "black",
                                          vjust = 1.0, hjust = 0.5)) +
          theme(legend.position = "none")

        ggsave(paste0(kingdom, "_", "", ec, "_", metric, "by_treatment", ".pdf"),
               p, width = 7.4, height = 3.95, path = OUT,
               useDingbats = FALSE, units = "cm")
        if (!.legend_done_treat) {
          .save_legend(p, "Fig1_treatment_legend", OUT)
          .legend_done_treat <- TRUE
        }
      }
    }
  }
}

{
  library(permute); library(vegan); library(ggplot2); library(ggpubr)
  library(agricolae)
  if (!exists("palette_density")) palette_density <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")

  .fmt_y <- function(v) sprintf("%.1f", v)
  .brk_y <- function(lim) {
    step <- if (lim[2] > 0.65) 0.2 else 0.1
    seq(0, lim[2], by = step)
  }

  .std01 <- function(x) {
    x <- as.numeric(x); r <- range(x, na.rm = TRUE)
    if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
    (x - r[1]) / diff(r)
  }

  compute_euclid_dist <- function(df, metrics = c("richness", "shannon")) {
    kd <- as.character(df$kingdom[1])
    pair_rows <- list(); samp_rows <- list(); cell_rows <- list(); test_rows <- list()
    letters_rows <- list()
    for (ec in levels(df$ecosystem)) {
      sub <- df[df$ecosystem == ec, ]
      if (nrow(sub) < 6) next

      Z   <- sapply(metrics, function(m) .std01(sub[[m]]))
      Z   <- matrix(Z, nrow = nrow(sub), ncol = length(metrics))
      grp <- as.character(sub$group)
      trt <- as.character(sub$group1)
      sid <- paste0(trt, ".", ave(seq_len(nrow(sub)), grp, trt, FUN = seq_along))

      for (gp in levels(sub$group)) {
        idx <- which(grp == gp)
        if (length(idx) < 2) next
        n_i <- length(idx)
        tt  <- trt[idx]
        Dm  <- as.matrix(stats::dist(Z[idx, , drop = FALSE], method = "euclidean"))

        mx  <- outer(tt, tt, FUN = "!=")
        diag(mx) <- FALSE
        ri  <- row(mx)[mx]; ci <- col(mx)[mx]
        pair_rows[[length(pair_rows) + 1]] <- data.frame(
          kingdom = kd, ecosystem = ec, group = gp,
          sample_i = sid[idx][ri], sample_j = sid[idx][ci],
          treatment_i = tt[ri], treatment_j = tt[ci],
          same_treatment = FALSE,
          dist = as.numeric(Dm[mx]), stringsAsFactors = FALSE)

        per <- rowSums(Dm * mx) / pmax(rowSums(mx), 1)
        samp_rows[[length(samp_rows) + 1]] <- data.frame(
          kingdom = kd, ecosystem = ec, group = gp,
          sample = sid[idx], treatment = tt,
          dist = as.numeric(per), stringsAsFactors = FALSE)
      }

      for (gp in levels(sub$group)) for (tt in levels(sub$group1)) {
        idx <- which(grp == gp & trt == tt)
        if (length(idx) < 2) next
        dd <- as.numeric(stats::dist(Z[idx, , drop = FALSE], method = "euclidean"))
        cell_rows[[length(cell_rows) + 1]] <- data.frame(
          kingdom = kd, ecosystem = ec, group = gp, group1 = tt,
          n = length(idx), n_pairs = length(dd),
          mean_dist = mean(dd), sd_dist = stats::sd(dd), stringsAsFactors = FALSE)
      }
    }
    pairs <- do.call(rbind, pair_rows)
    samp  <- do.call(rbind, samp_rows)
    cell  <- do.call(rbind, cell_rows)

    .key    <- function(d) paste(d$ecosystem, d$group, sep = "|")
    .n_samp <- tapply(seq_len(nrow(samp)), .key(samp), length)
    .n_trt  <- tapply(as.character(samp$treatment), .key(samp), function(x) length(unique(x)))
    dens <- do.call(rbind, lapply(split(pairs, .key(pairs)), function(d) {
      k  <- .key(d)[1]
      nu <- nrow(d) / 2
      data.frame(
        kingdom = d$kingdom[1], ecosystem = d$ecosystem[1], group = d$group[1],
        n = unname(.n_samp[k]),
        n_pairs_ordered = nrow(d), n_pairs_unordered = nu,
        mean_dist = mean(d$dist),
        se_dist   = stats::sd(d$dist) / sqrt(nu),
        sd_dist   = stats::sd(d$dist),
        n_trt     = unname(.n_trt[k]), stringsAsFactors = FALSE)
    }))
    dens <- dens[order(match(dens$ecosystem, levels(df$ecosystem)),
                       match(dens$group, levels(df$group))), ]
    rownames(dens) <- NULL

    for (ec in levels(df$ecosystem)) {
      d <- samp[samp$ecosystem == ec, ]
      if (nrow(d) < 8 || length(unique(d$group)) < 2) next
      d$group <- factor(d$group, levels = levels(df$group))
      m <- suppressWarnings(aov(dist ~ group, data = d))
      s <- summary(m)[[1]]; ss <- as.data.frame(s); ss$term <- trimws(rownames(s))
      resid <- ss[ss$term == "Residuals", ]; keep <- ss$term != "Residuals"

      lsd <- tryCatch(agricolae::LSD.test(m, "group", alpha = 0.05),
                      error = function(e) NULL)
      if (!is.null(lsd)) {
        gt  <- lsd$groups
        key <- trimws(rownames(gt))

        .ms  <- lsd$statistics[["MSerror"]]
        .dfe <- lsd$statistics[["Df"]]
        if (is.null(.dfe)) .dfe <- resid$Df[1]
        .lsd <- lsd$statistics[["LSD"]]

        for (g in levels(d$group)) {
          if (!(g %in% key)) next
          letters_rows[[length(letters_rows) + 1]] <- data.frame(
            kingdom = kd, ecosystem = ec, density = g,
            mean_c_i = unname(gt[g, "dist"]),
            MSE = as.numeric(.ms), df_error = as.numeric(.dfe),
            LSD = as.numeric(.lsd),
            letter = unname(as.character(gt[g, "groups"])),
            stringsAsFactors = FALSE)
        }
      }
      test_rows[[length(test_rows) + 1]] <- data.frame(
        kingdom = kd, ecosystem = ec,
        response = "sample mean pairwise Euclidean distance (ecosystem-wise 0-1 standardised richness+shannon)",
        term = ss$term[keep], df1 = ss$Df[keep], df2 = resid$Df[1],
        F_val = ss$`F value`[keep], p_value = ss$`Pr(>F)`[keep],
        partial_eta2 = ss$`Sum Sq`[keep] / (ss$`Sum Sq`[keep] + resid$`Sum Sq`[1]),
        stringsAsFactors = FALSE)
    }
    list(pairs = pairs, cell = cell, sample = samp, density = dens,
         letters = if (length(letters_rows)) do.call(rbind, letters_rows) else NULL,
         test = if (length(test_rows)) do.call(rbind, test_rows) else NULL)
  }

  for (kingdom in c("ITS", "16S")) {
    df  <- read_alpha(kingdom)
    res <- compute_euclid_dist(df)

    write.csv(res$pairs,   file.path(CACHE, paste0("euclid_pairs_",   kingdom, ".csv")), row.names = FALSE)
    write.csv(res$cell,    file.path(CACHE, paste0("euclid_cell_",    kingdom, ".csv")), row.names = FALSE)
    write.csv(res$sample,  file.path(CACHE, paste0("euclid_sample_",  kingdom, ".csv")), row.names = FALSE)
    write.csv(res$density, file.path(STATS, paste0("euclid_density_", kingdom, ".csv")), row.names = FALSE)
    if (!is.null(res$letters))
      write.csv(res$letters, file.path(STATS, paste0("euclid_letters_", kingdom, ".csv")), row.names = FALSE)
    if (!is.null(res$test))
      write.csv(res$test,  file.path(STATS, paste0("euclid_test_",    kingdom, ".csv")), row.names = FALSE)
    cat("[1B]", kingdom, "| within-ecosystem pairwise Euclidean distance after 0-1 standardisation (richness + shannon)\n")
    cat("      pairs =", nrow(res$pairs), "| sample =", nrow(res$sample),
        "| density =", nrow(res$density),
        "| letters =", if (is.null(res$letters)) 0 else nrow(res$letters),
        "| test =", if (is.null(res$test)) 0 else nrow(res$test), "\n")
    if (!is.null(res$test))
      for (i in seq_len(nrow(res$test)))
        cat("      ", as.character(res$test$ecosystem[i]), "|", res$test$term[i],
            "| F =", round(res$test$F_val[i], 2),
            "| df =", res$test$df1[i], "/", res$test$df2[i],
            "| p =", formatC(res$test$p_value[i], format = "e", digits = 2), "\n")
    if (!is.null(res$letters))
      for (ec_i in levels(df$ecosystem)) {
        .lt <- res$letters[res$letters$ecosystem == ec_i, ]
        if (nrow(.lt) == 0) next
        cat("      LSD", ec_i, ":",
            paste(.lt$density, .lt$letter, sep = "=", collapse = " "), "\n")
      }

    dens <- res$density

    dtop <- max(dens$mean_dist + dens$se_dist, na.rm = TRUE)
    if (!is.finite(dtop)) dtop <- 1

    HAS_LSD <- !is.null(res$letters) && nrow(res$letters) > 0
    LY   <- if (HAS_LSD) 1.46 else 1.25
    RS_Y <- if (HAS_LSD) 1.42 else 1.19

    for (eco in levels(df$ecosystem)) {
      d  <- dens[dens$ecosystem == eco, ]
      d$group  <- factor(.disp_sat(d$group),  levels = DENSITY_LABELS)
      d  <- d[order(d$group), ]

      p_sub <- ggplot() +
        geom_col(data = d, aes(x = group, y = mean_dist, fill = group),
                 width = 0.7, colour = "black") +

        geom_errorbar(data = d, aes(x = group, ymin = mean_dist - se_dist,
                                    ymax = mean_dist + se_dist),
                      width = 0.18, linewidth = 0.4, colour = "black") +
        scale_fill_manual(values = setNames(palette_density, DENSITY_LABELS)) +
        scale_y_continuous(limits = c(0, dtop * LY),
                           breaks = .brk_y, labels = .fmt_y,
                           expand = expansion(mult = c(0, 0.02))) +

        ylab("Mean pairwise distance") +
        xlab("Propagule density(%)") +
        labs(fill = "Propagule density(%)") +
        theme(panel.background = element_rect(fill = "transparent", colour = "black")) +

        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text.y = element_text(size = 6, color = "black")) +
        theme(axis.title.y = element_text(size = 8, color = "black")) +
        theme(axis.title.x = element_text(size = 8, color = "black")) +
        theme(axis.text.x = element_text(size = 6, color = "black",
                                         angle = 45, vjust = 1.0, hjust = 1.0)) +
        theme(legend.position = "none")

      if (!is.null(res$letters)) {
        .lt <- res$letters[res$letters$ecosystem == eco, ]
        if (nrow(.lt) > 0) {
          .lt$density <- factor(.disp_sat(.lt$density), levels = DENSITY_LABELS)
          .lt <- merge(.lt, d[, c("group", "mean_dist", "se_dist")],
                       by.x = "density", by.y = "group", all.x = TRUE)
          p_sub <- p_sub + geom_text(data = .lt,
                                     aes(x = density, y = mean_dist + se_dist,
                                         label = letter),
                                     vjust = -0.45, size = 2.1, colour = "black")
        }
      }

      if (!is.null(res$test)) {
        .st <- res$test[res$test$ecosystem == eco, ]
        if (nrow(.st) >= 1) {
          .p1 <- .st$p_value[1]
          .sym1 <- ifelse(.p1 < 0.001, "***", ifelse(.p1 < 0.01, "**",
                           ifelse(.p1 < 0.05, "*", "ns")))
          p_sub <- p_sub + annotate("text", x = 0.55, y = dtop * RS_Y,
                                    label = paste0("RS", .sym1),
                                    hjust = 0, vjust = 1, size = 2.1, colour = "black")
        }
      }

      ggsave(paste0(kingdom, eco, "euclidean_distance.pdf"), p_sub,
             width = 4.8, height = 4.2, path = OUT, units = "cm",
             useDingbats = FALSE)
      if (!isTRUE(get0(".legend_done_satru"))) {
        .save_legend(p_sub, "Fig1_density_legend", OUT)
        .legend_done_satru <- TRUE
      }
    }
  }

  for (kingdom in c("ITS", "16S")) .write_cv_legacy(read_alpha(kingdom), kingdom)
}

{
  library(ggplot2)

  KINGDOM_TITLE <- c(`16S` = "Prokaryotes (16S)", ITS = "Fungi (ITS)")
  if (!exists("palette_density")) {
    palette_density <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")
  }

  for (kingdom in c("ITS", "16S")) {
    df  <- read_alpha(kingdom)
    res <- compute_density_emmeans(df, "richness")
    if (is.null(res)) { message("[1C] skip ", kingdom, ": no emmeans result"); next }

    write.csv(res$emmeans,
              file.path(STATS, paste0("density_emmeans_", kingdom, "_richness.csv")),
              row.names = FALSE)
    write.csv(res$tests,
              file.path(STATS, paste0("density_model_tests_", kingdom, "_richness.csv")),
              row.names = FALSE)

    write.csv(res$emmeans[, c("kingdom", "ecosystem", "density", "emmean",
                              "se", "df", "lower", "upper", "letter", "n")],
              file.path(STATS, paste0("density_letters_lsd_", kingdom, "_richness.csv")),
              row.names = FALSE)
    cat("[1C]", kingdom, "| marginal-mean rows:", nrow(res$emmeans),
        "| density main effect p:",
        paste(formatC(res$tests$p_value[res$tests$term == "group"],
                      format = "e", digits = 2), collapse = ", "), "\n")
    cat("[1C]", kingdom, "| density LSD letters:",
        paste(res$emmeans$density, res$emmeans$letter, sep = "=", collapse = " "), "\n")

    p <- ggplot() +

      geom_col(data = res$emmeans,
               aes(x = ecosystem, y = emmean, fill = density),
               position = position_dodge(width = 0.8), width = 0.7,
               colour = "black") +

      geom_errorbar(data = res$emmeans,
                    aes(x = ecosystem, ymin = lower, ymax = upper,
                        group = density),
                    position = position_dodge(width = 0.8), width = 0.15,
                    linewidth = 0.4, colour = "black") +

      geom_text(data = res$emmeans,
                aes(x = ecosystem, y = upper, label = letter, group = density),
                position = position_dodge(width = 0.8),
                vjust = -0.45, size = 1.8, colour = "black",
                inherit.aes = FALSE, show.legend = FALSE) +
      scale_fill_manual(values = setNames(palette_density, DENSITY_LABELS)) +
      scale_y_continuous(limits = c(0, NA),
                         expand = expansion(mult = c(0, 0.16))) +
      labs(x = "Ecosystem", y = "Richness",
           title = unname(KINGDOM_TITLE[kingdom]), fill = "Propagule density(%)") +
      ggpubr::theme_pubr() +
      theme(panel.background = element_rect(fill = "transparent", colour = "black"),
            plot.title = element_text(size = 7, hjust = 0, colour = "black"),
            axis.text  = element_text(size = 6, colour = "black"),
            axis.title = element_text(size = 8, colour = "black"),
            legend.position = "none")

    ggsave(paste0("Fig1_density_main_effect_", kingdom, "_richness.pdf"),
           p, width = 7.4, height = 4, units = "cm", path = OUT, useDingbats = FALSE)
  }
}

{
  library(vegan); library(agricolae)

  for (kingdom in c("ITS", "16S")) {
    df <- read_alpha(kingdom)
    for (ecosys in unique(as.character(df$ecosystem))) {
      tempdata <- subset(df, ecosystem == ecosys)
      if (nrow(tempdata) < 2) next
      rich_data <- tempdata[, c("richness", "group", "group1")]

      mod1 <- aov(richness ~ group * group1, data = rich_data)
      cat(paste0(kingdom, "_richness_ecosystem", ecosys,
                 "_group: resource-surplus level, group1: treatment (with S x T interaction)-------\n"))
      print(summary(mod1))

      shnnon_data <- tempdata[, c("shannon", "group", "group1")]
      mod1 <- aov(shannon ~ group * group1, data = shnnon_data)
      cat(paste0(kingdom, "_shannon_ecosystem", ecosys,
                 "_group: resource-surplus level, group1: treatment-------\n"))
      print(summary(mod1))
    }
  }

  PClist <- list.files(path = file.path(DATA, "soil_properties/nitrogen/redone"),
                       pattern = "*.csv")
  reresult <- c()
  for (PCmatrix in PClist) {
    PCdata <- read.csv(file.path(DATA, "soil_properties/nitrogen/redone", PCmatrix),
                       header = TRUE)
    colnames(PCdata)[1] <- "value"
    savename <- substr(PCmatrix, 1, nchar(PCmatrix) - 4)
    for (ecosys in unique(PCdata$ecosystem)) {
      tempdata <- subset(PCdata, ecosystem == ecosys)
      mod1 <- aov(value ~ group1 + group, data = tempdata)
      cat(paste0(savename, "ecosystem", ecosys,
                 "_group1: resource-surplus level, group: treatment-------\n"))
      print(summary(mod1))
    }
  }
}
