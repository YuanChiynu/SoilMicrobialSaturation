DATA  <- "data"
OUT   <- file.path("output", "FigS4")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "FigS4")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

.norm_eco <- function(x) { x <- as.character(x); x[x == "Corplands"] <- "Croplands"; x }

.relabel_treat <- function(x) {
  x <- as.character(x)
  map <- c(CK = "CK", C = "C", N = "N", pH = "pH", T = "W", W = "D")
  out <- map[x]
  ifelse(is.na(out), x, out)
}

.ylab_main <- function(v) switch(v,
  "NH4" = bquote(NH[4]^"+"),
  "NO3" = bquote(NO[3]^"-"),
  "ph"  = bquote(pH),
  bquote(.(v)))

.sym <- function(p) ifelse(is.na(p), "ns",
                    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                    ifelse(p < 0.05,  "*", "ns"))))

.fmtp <- function(p) ifelse(is.na(p), "NA",
                     ifelse(p < 1e-3, formatC(p, format = "e", digits = 2),
                            sprintf("%.3f", p)))

.std01 <- function(x) {
  x <- as.numeric(x); r <- range(x, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

.st_lab <- function(d) {
  m  <- summary(aov(value ~ group1 * group, data = d))[[1]]
  rn <- trimws(rownames(m))
  ps <- m[rn == "group1", "Pr(>F)"]
  pt <- m[rn == "group",  "Pr(>F)"]
  paste0('RS^"', .sym(ps), '"~~T^"', .sym(pt), '"')
}
.st_rows <- function(d, v, e) {
  m  <- summary(aov(value ~ group1 * group, data = d))[[1]]
  rn <- trimws(rownames(m)); ir <- which(rn == "Residuals")
  keep <- setdiff(seq_along(rn), ir)
  data.frame(
    variable = v, ecosystem = e,
    term = c("RS (propagule density)", "T (treatment)", "RS x T"),
    df1  = m[keep, "Df"], df2 = m[ir, "Df"],
    F_value = m[keep, "F value"], p = m[keep, "Pr(>F)"],
    partial_eta2 = m[keep, "Sum Sq"] / (m[keep, "Sum Sq"] + m[ir, "Sum Sq"]),
    symbol = .sym(m[keep, "Pr(>F)"]),
    stringsAsFactors = FALSE)
}

PCLIST <- c(file.path(DATA, "soil_properties/nitrogen/redone", c("NH4.csv", "NO3.csv")),
            file.path(DATA, "soil_properties",          c("DOC.csv", "AP.csv", "ph.csv")))

{
  library(agricolae)

  st_rows       <- list()
  anova_rows    <- list()
  pairwise_rows <- list()
  lsd_mismatch  <- 0L
  reresult      <- NULL

  for (PCmatrix in PCLIST) {
    savename <- tools::file_path_sans_ext(basename(PCmatrix))
    PCdata   <- read.csv(PCmatrix, header = TRUE)
    colnames(PCdata)[1] <- "value"
    PCdata$ecosystem <- .norm_eco(PCdata$ecosystem)

    for (ecosys in c("Croplands", "Forests", "Grasslands")) {
      if (!ecosys %in% unique(PCdata$ecosystem)) next
      tempdata <- subset(PCdata, ecosystem == ecosys)
      tempdata$group1 = factor(tempdata$group1, levels = c("10%","40%","70%","100%"))
      tempdata$group  = factor(.relabel_treat(tempdata$group),
                               levels = c("CK","C","N","pH","W","D"))
      tempdata$ecosystem <- factor(tempdata$ecosystem,
                                   levels = c("Croplands","Forests","Grasslands"))

      st_rows[[length(st_rows) + 1]] <- .st_rows(tempdata, savename, ecosys)

      for (satru in levels(tempdata$group1)) {
        sub_phi <- subset(tempdata, group1 == satru)
        sub_phi$group <- factor(sub_phi$group, levels = c("CK","C","N","pH","W","D"))
        mod1 <- aov(value ~ group, data = sub_phi)
        re   <- LSD.test(mod1, "group", alpha = 0.05)
        re1  <- re$groups
        re1$eco        <- ecosys
        re1$satruation <- satru
        re1$treat      <- rownames(re1)
        re1$PC         <- savename
        reresult <- rbind(reresult, re1)

        sm <- summary(mod1)[[1]]
        rn <- trimws(rownames(sm))
        ig <- which(rn == "group"); ir <- which(rn == "Residuals")
        if (length(ig) == 1 && length(ir) == 1) {
          anova_rows[[length(anova_rows) + 1]] <- data.frame(
            variable = savename, ecosystem = ecosys, density = satru,
            df1 = sm[ig, "Df"], df2 = sm[ir, "Df"],
            F_value = sm[ig, "F value"], p = sm[ig, "Pr(>F)"],
            stringsAsFactors = FALSE)
        }
        pt <- tryCatch(pairwise.t.test(sub_phi$value, sub_phi$group,
                                       p.adjust.method = "none", pool.sd = TRUE),
                       error = function(e) NULL)
        if (!is.null(pt)) {
          pmat <- pt$p.value; gnm <- rownames(re1$groups)
          for (i in seq_len(nrow(pmat))) for (j in seq_len(ncol(pmat))) {
            pv <- pmat[i, j]; if (is.na(pv)) next
            a <- rownames(pmat)[i]; b <- colnames(pmat)[j]
            ma <- mean(sub_phi$value[sub_phi$group == a])
            mb <- mean(sub_phi$value[sub_phi$group == b])
            pairwise_rows[[length(pairwise_rows) + 1]] <- data.frame(
              variable = savename, ecosystem = ecosys, density = satru,
              group1 = a, group2 = b, mean1 = ma, mean2 = mb,
              mean_diff = ma - mb, p = pv, stringsAsFactors = FALSE)

            if (a %in% gnm && b %in% gnm) {
              lsd_sig <- re1$groups[gnm == a, "groups"] != re1$groups[gnm == b, "groups"]
              if (lsd_sig != (pv < 0.05)) lsd_mismatch <- lsd_mismatch + 1L
            }
          }
        }
      }
    }
  }

  write.csv(do.call(rbind, st_rows), file.path(STATS, "figS4_ST_effects.csv"),
            row.names = FALSE)
  write.csv(do.call(rbind, anova_rows), file.path(STATS, "figS4_anova_by_density.csv"),
            row.names = FALSE)
  write.csv(do.call(rbind, pairwise_rows), file.path(STATS, "figS4_pairwise_lsd.csv"),
            row.names = FALSE)
  write.csv(reresult, file = file.path(STATS, "soil_properties_multiple_comparisons.csv"), row.names = FALSE)

  cat(sprintf("[FigS4-A] RS/T/RSxT rows = %d | anova rows = %d | pairwise rows = %d | LSD letter mismatch = %d\n",
              length(st_rows) * 3L, length(anova_rows), length(pairwise_rows), lsd_mismatch))
}

{
  library(ggpubr); library(ggplot2)

  palette_treat <- c("#dddddd","#8DaF89","#F0a0a0","#FFB347","#FFEA70","#8a9cDa")
  .lab <- read.csv(file.path(STATS, "soil_properties_multiple_comparisons.csv"), stringsAsFactors = FALSE)
  .legend_done_treat <- FALSE

  for (PCmatrix in PCLIST) {
    savename <- tools::file_path_sans_ext(basename(PCmatrix))
    PCdata   <- read.csv(PCmatrix, header = TRUE)
    colnames(PCdata)[1] <- "value"
    PCdata$ecosystem <- .norm_eco(PCdata$ecosystem)

    for (ecosys in c("Croplands", "Forests", "Grasslands")) {
      if (!ecosys %in% unique(PCdata$ecosystem)) next
      pre_output <- paste0(savename, "_", ecosys)
      tempdata <- subset(PCdata, ecosystem == ecosys)
      tempdata$group1 = factor(tempdata$group1, levels = c("10%","40%","70%","100%"))
      tempdata$group  = factor(.relabel_treat(tempdata$group),
                               levels = c("CK","C","N","pH","W","D"))

      .ymax0  <- max(tempdata$value, na.rm = TRUE)
      .sub <- .lab[.lab$PC == savename & .lab$eco == ecosys, ]
      if (nrow(.sub)) {
        .sub$satruation <- factor(.sub$satruation, levels = levels(tempdata$group1))
        .sub$treat      <- factor(.sub$treat,      levels = levels(tempdata$group))
        .K <- 85.545; .T1MIN <- 0.055; .MAXPUSH <- 3L
        .FONT   <- 4.0
        .SPACING <- 4.69
        .SEP  <- .FONT * 1.3733 + 2.0
        .NEED <- 7.66 + .FONT * 1.3167 + 1.5
        Rk  <- max(1.16 * .ymax0, 1.30 * max(.sub$value, na.rm = TRUE))
        for (.it in 1:12) {
          .u  <- .K / Rk
          .t1 <- max(.T1MIN * .ymax0, 1.5 / .u)
          .sep <- .SEP / .u
          .ss <- .sub[order(as.integer(.sub$satruation), as.integer(.sub$treat)), ]
          .n  <- nrow(.ss)
          .wd <- 0.53 * .FONT * nchar(as.character(.ss$groups))
          .yy <- numeric(.n); .grp <- integer(.n)
          for (.i in seq_len(.n)) {
            .grp[.i] <- as.integer(.ss$satruation[.i])
            .y <- .ss$value[.i] + .t1
            for (.k in 0:.MAXPUSH) {
              .bad <- FALSE
              for (.j in seq_len(.i - 1L)) {
                if (.i == 1L || .grp[.j] != .grp[.i]) next
                if ((.i - .j) * .SPACING >= (.wd[.i] + .wd[.j]) / 2) next
                if (abs(.yy[.j] - .y) < .sep) { .bad <- TRUE; break }
              }
              if (!.bad || .k == .MAXPUSH) break
              .y <- .y + .sep
            }
            .yy[.i] <- .y
          }
          .ss$y <- .yy
          .sub  <- .ss
          Rk <- max(1.16 * .ymax0, (max(.sub$y) + .NEED / .u) / 0.99)
        }
        .ylim_top <- Rk
        .st_y     <- 0.99 * Rk
      } else {
        .ylim_top <- 1.16 * .ymax0
        .st_y     <- 1.14 * .ymax0
      }

      p_sub <- ggbarplot(tempdata, x="group1", y="value", add = "mean_se",
                         color = "black", fill = "group",
                         palette = palette_treat,
                         position = position_dodge(0.8), width = 0.7) +
        ylab(.ylab_main(savename)) +
        xlab("") +
        ylim(0, .ylim_top) +
        labs(fill = "Treatments") +

        annotate("text", x = 0.55, y = .st_y, label = .st_lab(tempdata),
                 hjust = 0, vjust = 1, size = 2.1, parse = TRUE) +
        theme(panel.background = element_rect(fill = "transparent", colour = "black")) +
        theme(axis.text.y  = element_text(size = 6, color = "black")) +
        theme(axis.title.y = element_text(size = 8, color = "black")) +
        theme(axis.title.x = element_text(size = 8, color = "black")) +
        theme(legend.text  = element_text(size = 8, color = "black")) +
        theme(axis.text.x  = element_text(size = 6, color = "black",
                                          vjust = 1.0, hjust = 1.0)) +
        theme(legend.position = "none")

      p_sub <- p_sub +
        geom_point(data = tempdata, aes(x = group1, y = value, fill = group, group = group),
                   shape = 21, colour = "black", stroke = 0.12,
                   size = 0.8, alpha = 0.95,
                   position = position_jitterdodge(jitter.width = 0.09, dodge.width = 0.8),
                   inherit.aes = FALSE, show.legend = FALSE) +
        scale_fill_manual(values = setNames(palette_treat, levels(tempdata$group)))

      if (nrow(.sub)) {
        p_sub <- p_sub + geom_text(data = .sub,
                     aes(x = satruation, y = y, label = groups, group = treat),
                     position = position_dodge(width = 0.8),
                     size = 1.406, vjust = -0.6, inherit.aes = FALSE)
      }

      ggsave(paste0(pre_output, "_fin.pdf"), p_sub,
             width = 6.55, height = 4.49, path = OUT, units = "cm", useDingbats = FALSE)
      if (!.legend_done_treat) {
        .save_legend(p_sub, "FigS4_treatment_legend", OUT)
        .legend_done_treat <- TRUE
      }
    }
  }
}

{
  library(ggplot2)

  palette_density <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")
  DENSITY_LABELS    <- c("10%", "40%", "70%", "100%")
  TRT_LABELS    <- c("CK", "C", "N", "pH", "W", "D")
  VARS_S4       <- c("DOC", "NH4", "NO3", "ph", "AP")
  .legend_done_satru <- FALSE

  .fmt_y <- function(v) sprintf("%.1f", v)
  .brk_y <- function(lim) {
    step <- if (lim[2] > 0.65) 0.2 else 0.1
    seq(0, lim[2], by = step)
  }

  pair_rows <- list(); samp_rows <- list(); cell_rows <- list()
  for (PCmatrix in PCLIST) {
    vname  <- tools::file_path_sans_ext(basename(PCmatrix))
    PCdata <- read.csv(PCmatrix, header = TRUE)
    colnames(PCdata)[1] <- "value"
    PCdata$ecosystem <- .norm_eco(PCdata$ecosystem)
    for (ecosys in c("Croplands", "Forests", "Grasslands")) {
      if (!ecosys %in% unique(PCdata$ecosystem)) next
      sub <- subset(PCdata, ecosystem == ecosys)
      sub$group1 <- factor(sub$group1, levels = DENSITY_LABELS)
      sub$group  <- factor(.relabel_treat(sub$group), levels = TRT_LABELS)

      sub$z <- .std01(sub$value)

      sub$sid <- paste0(as.character(sub$group), ".",
                        ave(seq_len(nrow(sub)), sub$group1, sub$group, FUN = seq_along))

      for (gp in DENSITY_LABELS) {
        idx <- which(as.character(sub$group1) == gp)
        if (length(idx) < 2) next
        n_i <- length(idx)
        ti  <- as.character(sub$group[idx])
        Dm  <- as.matrix(stats::dist(sub$z[idx], method = "euclidean"))

        mx  <- outer(ti, ti, FUN = "!=")
        diag(mx) <- FALSE
        ri  <- row(mx)[mx]; ci <- col(mx)[mx]
        pair_rows[[length(pair_rows) + 1]] <- data.frame(
          variable = vname, ecosystem = ecosys, density = gp,
          sample_i = sub$sid[idx][ri], sample_j = sub$sid[idx][ci],
          treatment_i = ti[ri], treatment_j = ti[ci],
          same_treatment = FALSE,
          dist = as.numeric(Dm[mx]), stringsAsFactors = FALSE)

        per <- rowSums(Dm * mx) / pmax(rowSums(mx), 1)
        samp_rows[[length(samp_rows) + 1]] <- data.frame(
          variable = vname, ecosystem = ecosys, density = gp,
          sample = sub$sid[idx], treatment = ti,
          dist = as.numeric(per), stringsAsFactors = FALSE)
      }

      for (gp in DENSITY_LABELS) for (tt in TRT_LABELS) {
        idx <- which(as.character(sub$group1) == gp & as.character(sub$group) == tt)
        if (length(idx) < 2) next
        dd <- as.numeric(stats::dist(sub$z[idx], method = "euclidean"))
        cell_rows[[length(cell_rows) + 1]] <- data.frame(
          variable = vname, ecosystem = ecosys, density = gp, treatment = tt,
          n = length(idx), n_pairs = length(dd),
          mean_dist = mean(dd), sd_dist = stats::sd(dd),
          stringsAsFactors = FALSE)
      }
    }
  }
  pairs <- do.call(rbind, pair_rows)
  samp  <- do.call(rbind, samp_rows)
  cell  <- do.call(rbind, cell_rows)
  .ord  <- function(d) d[order(match(d$variable, VARS_S4),
                               match(d$ecosystem, c("Croplands", "Forests", "Grasslands")),
                               match(d$density, DENSITY_LABELS)), ]
  samp <- .ord(samp); rownames(samp) <- NULL

  .key    <- function(d) paste(d$variable, d$ecosystem, d$density, sep = "|")
  .n_samp <- tapply(seq_len(nrow(samp)), .key(samp), length)
  .n_trt  <- tapply(as.character(samp$treatment), .key(samp), function(x) length(unique(x)))
  dens <- do.call(rbind, lapply(split(pairs, .key(pairs)), function(d) {
    k  <- .key(d)[1]
    nu <- nrow(d) / 2
    data.frame(
      variable = d$variable[1], ecosystem = d$ecosystem[1], density = d$density[1],
      n = unname(.n_samp[k]),
      n_pairs_ordered = nrow(d), n_pairs_unordered = nu,
      mean_dist = mean(d$dist),
      se_dist   = stats::sd(d$dist) / sqrt(nu),
      sd_dist   = stats::sd(d$dist),
      n_treat   = unname(.n_trt[k]), stringsAsFactors = FALSE)
  }))
  dens <- .ord(dens); rownames(dens) <- NULL

  test_rows <- list()
  for (v in VARS_S4) for (e in c("Croplands", "Forests", "Grasslands")) {
    d <- samp[samp$variable == v & samp$ecosystem == e, ]
    if (nrow(d) < 8 || length(unique(d$density)) < 2) next
    d$density <- factor(d$density, levels = DENSITY_LABELS)
    m <- suppressWarnings(aov(dist ~ density, data = d))
    s <- summary(m)[[1]]; ss <- as.data.frame(s); ss$term <- trimws(rownames(s))
    resid <- ss[ss$term == "Residuals", ]; keep <- ss$term != "Residuals"
    test_rows[[length(test_rows) + 1]] <- data.frame(
      variable = v, ecosystem = e,
      term = "RS (propagule density)",
      df1 = ss$Df[keep], df2 = resid$Df[1],
      F_value = ss$`F value`[keep], p = ss$`Pr(>F)`[keep],
      partial_eta2 = ss$`Sum Sq`[keep] / (ss$`Sum Sq`[keep] + resid$`Sum Sq`[1]),
      symbol = .sym(ss$`Pr(>F)`[keep]), stringsAsFactors = FALSE)
  }
  eu_test <- do.call(rbind, test_rows)

  mc_rows <- list()
  for (v in VARS_S4) for (e in c("Croplands", "Forests", "Grasslands")) {
    d <- samp[samp$variable == v & samp$ecosystem == e, ]
    if (nrow(d) < 8 || length(unique(d$density)) < 2) next
    d$density <- factor(d$density, levels = DENSITY_LABELS)
    m <- suppressWarnings(aov(dist ~ density, data = d))
    lt <- tryCatch(agricolae::LSD.test(m, "density", alpha = 0.05, group = TRUE),
                   error = function(e) NULL)
    if (is.null(lt) || is.null(lt$groups)) next

    MSE  <- lt$statistics[["MSerror"]]
    df2c <- lt$statistics[["Df"]]
    stopifnot(is.finite(MSE), is.finite(df2c), df2c > 0)
    tcrit <- stats::qt(0.975, df2c)

    grp  <- lt$groups
    lv   <- rownames(grp)
    stopifnot(length(lv) == 4, setequal(lv, DENSITY_LABELS))
    .gmu <- stats::setNames(as.numeric(grp[[1]]), lv)
    .glt <- stats::setNames(as.character(grp[[2]]), lv)
    .gn  <- stats::setNames(as.numeric(lt$means$r), rownames(lt$means))
    stopifnot(all(!is.na(.gmu)), all(nzchar(.glt)), all(.gn > 0),
              all(.gn[DENSITY_LABELS] > 0))
    .prow <- list()
    for (i in seq_len(length(lv) - 1)) for (j in (i + 1):length(lv)) {
      a <- lv[i]; b <- lv[j]
      n1 <- unname(.gn[a]); n2 <- unname(.gn[b])
      sed  <- sqrt(MSE * (1 / n1 + 1 / n2))
      md   <- unname(.gmu[a] - .gmu[b])
      .prow[[length(.prow) + 1]] <- data.frame(
        variable = v, ecosystem = e, level1 = a, level2 = b,
        mean1 = unname(.gmu[a]), mean2 = unname(.gmu[b]),
        n1 = n1, n2 = n2,
        mean_diff = md, se_diff = sed,
        t_value = md / sed, p = 2 * stats::pt(-abs(md / sed), df2c),
        lsd_critical = tcrit * sed,

        letters1 = .glt[a], letters2 = .glt[b],
        stringsAsFactors = FALSE)
    }
    .p <- do.call(rbind, .prow)
    .p$p_adj_BH <- stats::p.adjust(.p$p, method = "BH")
    .p$MSE <- MSE; .p$df2 <- df2c; .p$alpha <- 0.05
    .p$signif <- ifelse(.p$p < 0.05, "*", "ns")
    mc_rows[[length(mc_rows) + 1]] <- .p
  }
  eu_mc <- do.call(rbind, mc_rows)
  eu_mc <- eu_mc[, c("variable", "ecosystem", "level1", "level2", "mean1", "mean2",
                     "n1", "n2", "mean_diff", "se_diff", "t_value", "p", "p_adj_BH",
                     "lsd_critical", "letters1", "letters2", "signif", "MSE", "df2", "alpha")]

  lvl_rows <- list()
  for (v in VARS_S4) for (e in c("Croplands", "Forests", "Grasslands")) {
    d <- samp[samp$variable == v & samp$ecosystem == e, ]
    if (nrow(d) < 8 || length(unique(d$density)) < 2) next
    d$density <- factor(d$density, levels = DENSITY_LABELS)
    m <- suppressWarnings(aov(dist ~ density, data = d))
    lt <- tryCatch(agricolae::LSD.test(m, "density", alpha = 0.05, group = TRUE),
                   error = function(e) NULL)
    if (is.null(lt) || is.null(lt$groups)) next
    .m  <- lt$means
    .g  <- lt$groups
    lvl_rows[[length(lvl_rows) + 1]] <- data.frame(
      variable = v, ecosystem = e,
      density = rownames(.g),
      mean = as.numeric(.g[[1]]),
      letters = as.character(.g[[2]]),
      n = as.numeric(.m$r),
      rank_mean = rank(-as.numeric(.g[[1]])),
      stringsAsFactors = FALSE)
  }
  eu_letters <- do.call(rbind, lvl_rows)
  eu_letters <- eu_letters[order(match(eu_letters$variable, VARS_S4),
                                 match(eu_letters$ecosystem,
                                       c("Croplands", "Forests", "Grasslands")),
                                 match(eu_letters$density, DENSITY_LABELS)), ]
  rownames(eu_letters) <- NULL
  stopifnot(nrow(eu_letters) == 60,
            all(nzchar(eu_letters$letters)),
            all(eu_letters$n > 0))

  .lm <- merge(eu_mc, eu_letters[, c("variable", "ecosystem", "density", "letters")],
               by.x = c("variable", "ecosystem", "level1"),
               by.y = c("variable", "ecosystem", "density"), all.x = TRUE)
  .lm <- merge(.lm, eu_letters[, c("variable", "ecosystem", "density", "letters")],
               by.x = c("variable", "ecosystem", "level2"),
               by.y = c("variable", "ecosystem", "density"),
               all.x = TRUE, suffixes = c("_a", "_b"))
  .shared <- mapply(function(x, y) length(intersect(strsplit(x, "")[[1]],
                                                    strsplit(y, "")[[1]])) > 0,
                    .lm$letters_a, .lm$letters_b)
  stopifnot(all(.shared == (.lm$p >= 0.05)))

  write.csv(pairs,   file.path(CACHE, "figS4_euclid_pairs.csv"),   row.names = FALSE)
  write.csv(samp,    file.path(CACHE, "figS4_euclid_sample.csv"),  row.names = FALSE)
  write.csv(dens,    file.path(STATS, "figS4_euclid_density.csv"), row.names = FALSE)
  write.csv(eu_test, file.path(STATS, "figS4_euclid_test.csv"),    row.names = FALSE)
  write.csv(eu_mc,   file.path(STATS, "figS4_euclid_multiple_comparison.csv"),
            row.names = FALSE)
  write.csv(eu_letters, file.path(STATS, "figS4_euclid_LSD_letters.csv"),
            row.names = FALSE)
  write.csv(cell,    file.path(CACHE, "figS4_euclid_cell.csv"),    row.names = FALSE)

  .chk <- aggregate(list(np = pairs$dist),
                    by = list(variable = pairs$variable, ecosystem = pairs$ecosystem,
                              density = pairs$density), FUN = length)
  .chk <- merge(.chk, dens[, c("variable", "ecosystem", "density",
                               "n_pairs_ordered", "n_pairs_unordered")],
                by = c("variable", "ecosystem", "density"))
  stopifnot(nrow(dens) == 60, nrow(eu_test) == 15, nrow(cell) == 360,
            nrow(samp) == 1440,
            nrow(eu_mc) == 15 * 6,
            nrow(eu_letters) == 15 * 4,
            all(eu_mc$level1 != eu_mc$level2),
            all(eu_mc$n1 >= 20 & eu_mc$n1 <= 28),
            all(eu_mc$n2 >= 20 & eu_mc$n2 <= 28),
            all(is.finite(eu_mc$p)) | all(is.na(eu_mc$p)),
            all(abs(eu_mc$t_value) - abs(eu_mc$mean_diff / eu_mc$se_diff) < 1e-9),
            all(eu_mc$lsd_critical > 0),
            all(eu_mc$alpha == 0.05),
            all(!pairs$same_treatment),
            all(.chk$np == .chk$n_pairs_ordered),
            all(.chk$np == 2 * .chk$n_pairs_unordered),
            all(dens$n_pairs_unordered < dens$n * (dens$n - 1) / 2),
            all(dens$n %in% c(23L, 24L, 25L)),
            all(dens$n_pairs_ordered == 2 * dens$n_pairs_unordered),
            all(eu_test$df1 == 3), all(eu_test$df2 == 92),
            all(dens$n >= 20 & dens$n <= 28))

  cat(sprintf("[FigS4-B] euclid pairs = %d | sample = %d | density = %d | test = %d | LSD mc = %d | aux cell = %d\n",
              nrow(pairs), nrow(samp), nrow(dens), nrow(eu_test), nrow(eu_mc), nrow(cell)))
  cat("[FigS4-B] ordered pairs per group (n):",
      paste(sprintf("%s/%s/%s n=%d pairs=%d", dens$variable, dens$ecosystem, dens$density,
                    dens$n, dens$n_pairs_ordered), collapse = "; "), "\n")
  cat("[FigS4-B] RS F/p:",
      paste(sprintf("%s|%s F=%.2f df=%d/%d p=%s", eu_test$variable, eu_test$ecosystem,
                    eu_test$F_value, eu_test$df1, eu_test$df2,
                    formatC(eu_test$p, format = "e", digits = 2)), collapse = "; "), "\n")

  for (v in VARS_S4) {
    dv <- dens[dens$variable == v, ]
    Rvar <- max(dv$mean_dist + dv$se_dist, na.rm = TRUE)
    if (!is.finite(Rvar) || Rvar <= 0) Rvar <- 1
    y_top <- 1.25 * Rvar
    for (e in c("Croplands", "Forests", "Grasslands")) {
      if (!e %in% unique(dv$ecosystem)) next
      d  <- dv[dv$ecosystem == e, ]
      d$density  <- factor(d$density,  levels = DENSITY_LABELS)
      d  <- d[order(d$density), ]

      rr <- eu_test[eu_test$variable == v & eu_test$ecosystem == e, ]
      st_lab <- if (nrow(rr) >= 1) paste0('RS^"', rr$symbol[1], '"') else NA_character_

      p_sub <- ggplot() +
        geom_col(data = d, aes(x = density, y = mean_dist, fill = density),
                 width = 0.7, colour = "black") +

        geom_errorbar(data = d, aes(x = density, ymin = mean_dist - se_dist,
                                    ymax = mean_dist + se_dist),
                      width = 0.18, linewidth = 0.4, colour = "black") +
        scale_fill_manual(values = setNames(palette_density, DENSITY_LABELS)) +
        scale_y_continuous(limits = c(0, y_top), breaks = .brk_y, labels = .fmt_y) +
        ylab("Euclidean distance") +
        xlab("") +
        labs(fill = "Propagule density(%)") +
        theme(panel.background = element_rect(fill = "transparent", colour = "black")) +

        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text.y  = element_text(size = 6, color = "black")) +
        theme(axis.title.y = element_text(size = 8, color = "black")) +
        theme(axis.title.x = element_text(size = 8, color = "black")) +
        theme(legend.text  = element_text(size = 8, color = "black")) +
        theme(axis.text.x  = element_text(size = 6, color = "black",
                                          angle = 45, hjust = 1)) +
        theme(legend.position = "none")

      if (!is.na(st_lab)) {
        p_sub <- p_sub + annotate("text", x = 0.55, y = 0.99 * y_top, label = st_lab,
                                  hjust = 0, vjust = 1, size = 2.1, parse = TRUE)
      }

      .lt <- eu_letters[eu_letters$variable == v & eu_letters$ecosystem == e, ]
      if (nrow(.lt) >= 1) {
        .lt$density <- factor(.lt$density, levels = DENSITY_LABELS)
        .lt <- .lt[order(.lt$density), ]
        .mi <- match(as.character(.lt$density), as.character(d$density))
        stopifnot(!anyNA(.mi), nrow(.lt) == 4)
        .lt$y <- d$mean_dist[.mi] + d$se_dist[.mi] + 0.02 * y_top
        .lt$y <- pmin(.lt$y, 0.985 * y_top)
        stopifnot(all(is.finite(.lt$y)))
        p_sub <- p_sub + annotate("text", x = as.character(.lt$density), y = .lt$y,
                                  label = .lt$letters, size = 1.4, vjust = 0)
      }

      ggsave(paste0(v, "_", e, "_euclidean_distance.pdf"), p_sub,
             width = 3.15, height = 4.49, path = OUT, units = "cm", useDingbats = FALSE)
      if (!.legend_done_satru) {
        .save_legend(p_sub, "FigS4_density_legend", OUT)
        .legend_done_satru <- TRUE
      }
    }
  }
}
