DATA  <- "data"
OUT   <- file.path("output", "FigS5")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "FigS5")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

.norm_eco <- function(x) {
  x <- as.character(x)
  x[x %in% c("Corpland", "Cropland")] <- "Croplands"
  x[x %in% c("Forest")]               <- "Forests"
  x[x %in% c("Grass", "Grassland")]   <- "Grasslands"
  x
}
.relabel_treat <- function(x) {
  x <- as.character(x)
  map <- c(CK = "CK", C = "C", N = "N", pH = "pH", T = "W", W = "D")
  out <- map[x]
  ifelse(is.na(out), x, out)
}

.eco_levels <- function(x) {
  lv <- c("Croplands", "Forests", "Grasslands")
  lv[lv %in% unique(as.character(x))]
}
TREAT_LEVELS <- c("CK", "C", "N", "pH", "W", "D")

.fmtp <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-3, formatC(p, format = "e", digits = 2), sprintf("%.3f", p)))
}

.fmtp_fig <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 1e-3, "< 0.001", sprintf("= %.3f", p)))
}

.sym_star <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "ns"))))
}

.std01 <- function(x) {
  x <- as.numeric(x); r <- range(x, na.rm = TRUE)
  if (!is.finite(diff(r)) || diff(r) == 0) return(rep(0, length(x)))
  (x - r[1]) / diff(r)
}

YEXP <- 0.34

.cv_fm_test <- function(sub_df, value_col, grp_col,
                        sats = c("10%", "40%", "70%", "100%")) {
  if (!requireNamespace("cvequality", quietly = TRUE)) return(NULL)
  n <- s <- x <- numeric(0)
  for (st in sats) {
    v <- sub_df[sub_df[[grp_col]] == st, value_col]
    v <- v[is.finite(v)]
    n <- c(n, length(v)); s <- c(s, sd(v)); x <- c(x, mean(v))
  }
  if (length(n) < 2 || any(n < 2) || any(!is.finite(s)) || any(s <= 0)) return(NULL)
  tst <- tryCatch(cvequality::asymptotic_test2(k = length(sats), n = n, s = s, x = x),
                  error = function(e) NULL)
  if (is.null(tst)) return(NULL)
  data.frame(df = length(sats) - 1, D = tst$D_AD, p_FM = tst$p_value,
             stringsAsFactors = FALSE)
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

library(ggpubr)
library(agricolae)

MARKERS <- list(
  its = list(file = "ITS_copynum.csv",
             lab  = expression(atop("log10 fungal ITS",
                                    atop("gene copies", "per g dry soil"))),
             cvlab = "Fungal CV (%)"),
  `16s` = list(file = "16S_copynum.csv",
             lab  = expression(atop("log10 16S rRNA",
                                    atop("gene copies", "per g dry soil"))),
             cvlab = "Prokaryotic CV (%)")
)

cv_values_all <- list()
cv_test_all   <- list()
lsd_all       <- list()
st_all        <- list()

for (mk in names(MARKERS)) {
  spec <- MARKERS[[mk]]

  raw <- read.csv(file.path(DATA, "qpcr", spec$file), header = TRUE,
                  stringsAsFactors = FALSE)
  dat <- raw[, c(1, 9, 10, 11, 13)]
  colnames(dat) <- c("sample_id", "ecosystem", "group", "group1", "copynum")
  dat$copynum   <- as.numeric(dat$copynum)
  dat$ecosystem <- .norm_eco(dat$ecosystem)

  dat$group1 <- factor(.relabel_treat(dat$group1), levels = TREAT_LEVELS)

  dat$group  <- factor(gsub("re$", "", dat$group),
                       levels = c("10%", "40%", "70%", "100%"))
  dat$copynum_log <- log10(dat$copynum)

  st_mk <- do.call(rbind, lapply(.eco_levels(dat$ecosystem), function(eco) {
    sub <- subset(dat, ecosystem == eco)
    m   <- aov(copynum_log ~ group * group1, data = sub)
    a   <- summary(m)[[1]]
    rn  <- trimws(rownames(a))
    gv  <- function(term, col) a[match(term, rn), col]
    ss_res <- gv("Residuals", "Sum Sq")
    data.frame(marker = mk, ecosystem = eco,
               df1 = c(gv("group", "Df"), gv("group1", "Df"), gv("group:group1", "Df")),
               df2 = gv("Residuals", "Df"),
               term = c("resource surplus (RS)", "treatment (T)", "RS × T"),
               F_value = c(gv("group", "F value"), gv("group1", "F value"),
                           gv("group:group1", "F value")),
               p = c(gv("group", "Pr(>F)"), gv("group1", "Pr(>F)"),
                     gv("group:group1", "Pr(>F)")),
               partial_eta2 = c(gv("group", "Sum Sq"), gv("group1", "Sum Sq"),
                                gv("group:group1", "Sum Sq")) / ss_res,
               stringsAsFactors = FALSE)
  }))
  st_mk$symbol <- .sym_star(st_mk$p)
  st_all[[length(st_all) + 1]] <- st_mk

  for (ecosys in .eco_levels(dat$ecosystem)) {
    tempdata <- subset(dat, ecosystem == ecosys)
    y_top <- ceiling(max(tempdata$copynum_log, na.rm = TRUE))
    p1 <- ggbarplot(tempdata, x="group", y="copynum_log", add = "mean_se",color = "black",fill = "group1",
                  palette =c("#dddddd","#8DaF89","#F0a0a0","#FFB347","#FFEA70","#8a9cDa"), position = position_dodge(0.8))+
    ylab(spec$lab)+
    xlab("Propagule density(%)")+

    scale_y_continuous(expand = expansion(mult = c(0, YEXP)))+
    labs(fill = "Treatments")+

    theme(panel.background = element_rect(fill = "transparent", colour = "black"))+
    theme(axis.text.y = element_text(size = 6,color="black"))+
    theme(axis.title.y = element_text(size = 8, color="black"))+
    theme(axis.title.x = element_text(size = 8, color="black"))+
    theme(legend.position = "none")+

    theme(legend.text  = element_text(size = 8, color = "black"),
          legend.title = element_text(size = 9, color = "black"))+
    theme(axis.text.x = element_text(size = 6,color="black",vjust = 1.0, hjust = 1.0))

    p1 <- p1 +
      geom_point(data = tempdata, aes(x = group, y = copynum_log, fill = group1, group = group1),
                 shape = 21, colour = "black", stroke = 0.12,
                 size = 0.8, alpha = 0.95,
                 position = position_jitterdodge(jitter.width = 0.09, dodge.width = 0.8),
                 inherit.aes = FALSE, show.legend = FALSE) +

      scale_fill_manual(values = setNames(
        c("#dddddd","#8DaF89","#F0a0a0","#FFB347","#FFEA70","#8a9cDa"),
        levels(tempdata$group1)))

    .lt5 <- file.path(CACHE, "figS5_qpcr_multiple_comparisons.csv")
    if (file.exists(.lt5)) {
      .lab5 <- read.csv(.lt5, stringsAsFactors = FALSE)
      .sub5 <- .lab5[.lab5$marker == mk & .lab5$eco == ecosys, ]
      if (nrow(.sub5)) {
        .sub5$satruation <- factor(.sub5$satruation, levels = levels(tempdata$group))
        .sub5$treat      <- factor(.sub5$treat,      levels = levels(tempdata$group1))
        .sub5$y <- log10(.sub5$copynum) + 0.03 * y_top
        p1 <- p1 + geom_text(data = .sub5,
                  aes(x = satruation, y = y, label = groups, group = treat),
                  position = position_dodge(width = 0.8),
                  size = 2.1, vjust = -0.6, inherit.aes = FALSE)
      }
    }

    .st <- st_mk[st_mk$ecosystem == ecosys &
                 st_mk$term %in% c("resource surplus (RS)", "treatment (T)"), ]
    .st <- .st[match(c("resource surplus (RS)", "treatment (T)"), .st$term), ]
    if (nrow(.st) == 2 && all(!is.na(.st$symbol))) {
      .lab_st <- sprintf("RS%s T%s",
                         ifelse(.st$symbol[1] == "ns", " ns", .st$symbol[1]),
                         ifelse(.st$symbol[2] == "ns", " ns", .st$symbol[2]))
      p1 <- p1 + annotate("text", x = 0.50, y = Inf, label = .lab_st,
                          hjust = 0, vjust = 1.45, size = 2.813, colour = "black")
    }
    p1
    savename <- paste0(mk, ecosys, "_qpcr", ".pdf")

    ggsave(savename, p1, width = 9.6, height = 4.3, path = OUT, useDingbats = FALSE, units = "cm")

    if (ecosys != "Grasslands") {
      p1_mid <- p1 +
        theme(axis.text.x  = element_blank(),
              axis.title.x = element_blank())
      ggsave(paste0(mk, ecosys, "_qpcr_mid.pdf"), p1_mid,
             width = 9.6, height = 3.64, path = OUT, useDingbats = FALSE, units = "cm")
    }
    if (!exists(".figS5_leg_treat")) {
      .save_legend(p1, "FigS5_qpcr_treatment_legend", OUT)
      .figS5_leg_treat <- TRUE
    }
  }

  for (eco in .eco_levels(dat$ecosystem)) {
    data_eco <- subset(dat, ecosystem == eco)
    for (satru in c("10%", "40%", "70%", "100%")) {
      v <- data_eco$copynum[data_eco$group == satru]
      v <- v[is.finite(v)]
      ci <- .boot_cv_ci(v, B = 2000)
      cv_values_all[[length(cv_values_all) + 1]] <- data.frame(
        marker = mk, ecosystem = eco, satru = satru, n = length(v),
        pccv = sd(v) / mean(v) * 100, cv_lo = ci[1], cv_hi = ci[2],
        stringsAsFactors = FALSE)
    }
    tst <- .cv_fm_test(data_eco, "copynum", "group")
    if (!is.null(tst)) {
      cv_of <- function(s) {
        v <- data_eco$copynum[data_eco$group == s]; v <- v[is.finite(v)]
        sd(v) / mean(v) * 100
      }
      cv_test_all[[length(cv_test_all) + 1]] <- data.frame(
        marker = mk, ecosystem = eco, k = 4, df = tst$df,
        CV_10 = cv_of("10%"), CV_40 = cv_of("40%"),
        CV_70 = cv_of("70%"), CV_100 = cv_of("100%"),
        FeltzMiller_D = tst$D, p_FM = tst$p_FM,
        stringsAsFactors = FALSE)
    }
  }

  .fmt_y <- function(v) sprintf("%.1f", v)
  .brk_y <- function(lim) {
    step <- if (lim[2] > 0.65) 0.2 else 0.1
    seq(0, lim[2], by = step)
  }
  eu_pairs <- list(); eu_cell <- list(); eu_samp <- list(); eu_test <- list()
  for (eco in .eco_levels(dat$ecosystem)) {
    sub <- subset(dat, ecosystem == eco)
    sub <- sub[order(sub$group, sub$group1, sub$sample_id), ]
    if (nrow(sub) < 8) next
    Z <- .std01(sub$copynum)
    grp <- as.character(sub$group); trt <- as.character(sub$group1)

    for (gp in levels(dat$group)) {
      idx <- which(grp == gp)
      if (length(idx) < 2) next
      n_i <- length(idx)
      tt  <- trt[idx]
      Dm  <- as.matrix(stats::dist(matrix(Z[idx], ncol = 1), method = "euclidean"))

      mx  <- outer(tt, tt, FUN = "!=")
      diag(mx) <- FALSE
      ri  <- row(mx)[mx]; ci <- col(mx)[mx]
      eu_pairs[[length(eu_pairs) + 1]] <- data.frame(
        marker = mk, ecosystem = eco, group = gp,
        sample_i = sub$sample_id[idx][ri], sample_j = sub$sample_id[idx][ci],
        treatment_i = tt[ri], treatment_j = tt[ci],
        same_treatment = FALSE,
        dist = as.numeric(Dm[mx]), stringsAsFactors = FALSE)

      per <- rowSums(Dm * mx) / pmax(rowSums(mx), 1)
      eu_samp[[length(eu_samp) + 1]] <- data.frame(
        marker = mk, ecosystem = eco, group = gp,
        sample = sub$sample_id[idx], treatment = tt,
        dist = as.numeric(per), stringsAsFactors = FALSE)
    }

    for (gp in levels(dat$group)) for (tt in levels(dat$group1)) {
      idx <- which(grp == gp & trt == tt)
      if (length(idx) < 2) next
      dd <- as.numeric(stats::dist(matrix(Z[idx], ncol = 1), method = "euclidean"))
      eu_cell[[length(eu_cell) + 1]] <- data.frame(
        marker = mk, ecosystem = eco, group = gp, group1 = tt,
        n = length(idx), n_pairs = length(dd),
        mean_dist = mean(dd), sd_dist = stats::sd(dd), stringsAsFactors = FALSE)
    }
  }
  eu_pairs_df <- do.call(rbind, eu_pairs)
  eu_cell_df  <- do.call(rbind, eu_cell)
  eu_samp_df  <- do.call(rbind, eu_samp)
  .k5 <- paste(eu_samp_df$ecosystem, eu_samp_df$group, sep = "|")
  .n5 <- tapply(seq_len(nrow(eu_samp_df)), .k5, length)
  .t5 <- tapply(as.character(eu_samp_df$treatment), .k5, function(x) length(unique(x)))
  eu_dens_df <- do.call(rbind, lapply(
    split(eu_pairs_df, paste(eu_pairs_df$ecosystem, eu_pairs_df$group, sep = "|")),
    function(d) {
      k  <- paste(d$ecosystem[1], d$group[1], sep = "|")
      nu <- nrow(d) / 2
      data.frame(
        marker = mk, ecosystem = d$ecosystem[1], group = d$group[1],
        n = unname(.n5[k]),
        n_pairs_ordered = nrow(d), n_pairs_unordered = nu,
        mean_dist = mean(d$dist),
        se_dist   = stats::sd(d$dist) / sqrt(nu),
        sd_dist   = stats::sd(d$dist),
        n_trt     = unname(.t5[k]), stringsAsFactors = FALSE)
    }))
  eu_dens_df <- eu_dens_df[order(match(eu_dens_df$ecosystem, .eco_levels(dat$ecosystem)),
                                 match(eu_dens_df$group, levels(dat$group))), ]
  rownames(eu_dens_df) <- NULL

  for (eco in .eco_levels(dat$ecosystem)) {
    d <- eu_samp_df[eu_samp_df$ecosystem == eco, ]
    if (nrow(d) < 8 || length(unique(d$group)) < 2) next
    d$group <- factor(d$group, levels = levels(dat$group))
    m <- suppressWarnings(aov(dist ~ group, data = d))
    s <- summary(m)[[1]]; ss <- as.data.frame(s); ss$term <- trimws(rownames(s))
    resid <- ss[ss$term == "Residuals", ]; keep <- ss$term != "Residuals"
    eu_test[[length(eu_test) + 1]] <- data.frame(
      marker = mk, ecosystem = eco,
      response = "sample mean pairwise Euclidean distance (ecosystem-wise 0-1 standardised copynum)",
      term = "RS (propagule density)",
      df1 = ss$Df[keep], df2 = resid$Df[1],
      F_value = ss$`F value`[keep], p = ss$`Pr(>F)`[keep],
      partial_eta2 = ss$`Sum Sq`[keep] / (ss$`Sum Sq`[keep] + resid$`Sum Sq`[1]),
      stringsAsFactors = FALSE)
  }
  eu_test_df <- if (length(eu_test)) do.call(rbind, eu_test) else NULL
  if (!is.null(eu_test_df)) eu_test_df$sym <- .sym_star(eu_test_df$p)

  eu_mc_rows <- list(); eu_lvl_rows <- list()
  for (eco in .eco_levels(dat$ecosystem)) {
    d <- eu_samp_df[eu_samp_df$ecosystem == eco, ]
    if (nrow(d) < 8 || length(unique(d$group)) < 2) next
    d$group <- factor(d$group, levels = levels(dat$group))
    m <- suppressWarnings(aov(dist ~ group, data = d))
    lt <- tryCatch(agricolae::LSD.test(m, "group", alpha = 0.05, group = TRUE),
                   error = function(e) NULL)
    if (is.null(lt) || is.null(lt$groups)) next
    grp <- lt$groups; lv <- rownames(grp)
    MSE  <- lt$statistics[["MSerror"]]; df2c <- lt$statistics[["Df"]]
    stopifnot(is.finite(MSE), is.finite(df2c), df2c > 0)
    tcrit <- stats::qt(0.975, df2c)
    .gmu <- stats::setNames(as.numeric(grp[[1]]), lv)
    .glt <- stats::setNames(as.character(grp[[2]]), lv)
    .gn  <- stats::setNames(as.numeric(lt$means$r), rownames(lt$means))
    stopifnot(length(lv) == 4, setequal(lv, levels(dat$group)),
              all(!is.na(.gmu)), all(nzchar(.glt)), all(.gn > 0))
    .prow <- list()
    for (i in seq_len(length(lv) - 1)) for (j in (i + 1):length(lv)) {
      a <- lv[i]; b <- lv[j]
      n1 <- unname(.gn[a]); n2 <- unname(.gn[b])
      sed  <- sqrt(MSE * (1 / n1 + 1 / n2))
      md   <- unname(.gmu[a] - .gmu[b])
      .pv  <- 2 * stats::pt(-abs(md / sed), df2c)
      .prow[[length(.prow) + 1]] <- data.frame(
        marker = mk, ecosystem = eco, level1 = a, level2 = b,
        mean1 = unname(.gmu[a]), mean2 = unname(.gmu[b]),
        n1 = n1, n2 = n2, mean_diff = md, se_diff = sed,
        t_value = md / sed, p = .pv,
        lsd_critical = tcrit * sed,
        letters1 = .glt[a], letters2 = .glt[b],
        signif = ifelse(.pv < 0.05, "*", "ns"),
        MSE = MSE, df2 = df2c, alpha = 0.05, stringsAsFactors = FALSE)
    }
    .p <- do.call(rbind, .prow)
    .p$p_adj_BH <- stats::p.adjust(.p$p, method = "BH")
    eu_mc_rows[[length(eu_mc_rows) + 1]] <- .p
    eu_lvl_rows[[length(eu_lvl_rows) + 1]] <- data.frame(
      marker = mk, ecosystem = eco,
      density = lv, mean = as.numeric(grp[[1]]),
      letters = as.character(grp[[2]]),

      n = stats::setNames(as.numeric(lt$means$r), rownames(lt$means))[lv],
      stringsAsFactors = FALSE)
  }
  eu_mc_df      <- if (length(eu_mc_rows))  do.call(rbind, eu_mc_rows)  else NULL
  eu_letters_df <- if (length(eu_lvl_rows)) do.call(rbind, eu_lvl_rows) else NULL
  if (!is.null(eu_mc_df)) {
    eu_mc_df <- eu_mc_df[, c("marker","ecosystem","level1","level2","mean1","mean2",
                             "n1","n2","mean_diff","se_diff","t_value","p","p_adj_BH",
                             "lsd_critical","letters1","letters2","signif",
                             "MSE","df2","alpha")]
  }
  if (!is.null(eu_mc_df) && !is.null(eu_letters_df)) {

    .lm <- merge(eu_mc_df, eu_letters_df[, c("marker","ecosystem","density","letters")],
                 by.x = c("marker","ecosystem","level1"),
                 by.y = c("marker","ecosystem","density"), all.x = TRUE)
    .lm <- merge(.lm, eu_letters_df[, c("marker","ecosystem","density","letters")],
                 by.x = c("marker","ecosystem","level2"),
                 by.y = c("marker","ecosystem","density"),
                 all.x = TRUE, suffixes = c("_a", "_b"))
    .shared <- mapply(function(x, y) length(intersect(strsplit(x, "")[[1]],
                                                      strsplit(y, "")[[1]])) > 0,
                      .lm$letters_a, .lm$letters_b)
    stopifnot(nrow(eu_mc_df) == 3 * 6, nrow(eu_letters_df) == 3 * 4,
              all(.shared == (.lm$p >= 0.05)))
  }

  write.csv(eu_pairs_df, file.path(CACHE, paste0("figS5_euclid_pairs_",   mk, ".csv")), row.names = FALSE)
  write.csv(eu_cell_df,  file.path(CACHE, paste0("figS5_euclid_cell_",    mk, ".csv")), row.names = FALSE)
  write.csv(eu_samp_df,  file.path(CACHE, paste0("figS5_euclid_sample_",  mk, ".csv")), row.names = FALSE)
  write.csv(eu_dens_df,  file.path(STATS, paste0("figS5_euclid_density_", mk, ".csv")), row.names = FALSE)
  if (!is.null(eu_test_df))
    write.csv(eu_test_df, file.path(STATS, paste0("figS5_euclid_test_", mk, ".csv")), row.names = FALSE)
  if (!is.null(eu_mc_df))
    write.csv(eu_mc_df, file.path(STATS, paste0("figS5_euclid_multiple_comparison_", mk, ".csv")),
              row.names = FALSE)
  if (!is.null(eu_letters_df))
    write.csv(eu_letters_df, file.path(STATS, paste0("figS5_euclid_LSD_letters_", mk, ".csv")),
              row.names = FALSE)
  .chk5 <- aggregate(list(np = eu_pairs_df$dist),
                     by = list(ecosystem = eu_pairs_df$ecosystem, group = eu_pairs_df$group),
                     FUN = length)
  .chk5 <- merge(.chk5, eu_dens_df[, c("ecosystem", "group",
                                       "n_pairs_ordered", "n_pairs_unordered")],
                 by = c("ecosystem", "group"))
  stopifnot(nrow(eu_cell_df) == 3 * 24,
            nrow(eu_samp_df) == 3 * 96,
            nrow(eu_dens_df) == 3 * 4,
            all(!eu_pairs_df$same_treatment),
            all(.chk5$np == .chk5$n_pairs_ordered),
            all(.chk5$np == 2 * .chk5$n_pairs_unordered),
            all(eu_dens_df$n_pairs_ordered == 480),
            all(eu_dens_df$n >= 20 & eu_dens_df$n <= 28))
  if (!is.null(eu_test_df)) stopifnot(nrow(eu_test_df) == 3,
                                      all(eu_test_df$df1 == 3), all(eu_test_df$df2 == 92))
  cat("[2b]", mk, "| euclid pairs:", nrow(eu_pairs_df), "sample:", nrow(eu_samp_df),
      "density:", nrow(eu_dens_df),
      "test:", if (is.null(eu_test_df)) 0 else nrow(eu_test_df), "\n")

  dtop <- max(eu_dens_df$mean_dist + eu_dens_df$se_dist, na.rm = TRUE)
  if (!is.finite(dtop) || dtop <= 0) dtop <- 1
  for (eco in .eco_levels(dat$ecosystem)) {
    d  <- eu_dens_df[eu_dens_df$ecosystem == eco, ]
    d  <- d[order(d$group), ]
    d$group  <- factor(d$group,  levels = levels(dat$group))
    p_sub <- ggplot() +
      geom_col(data = d, aes(x = group, y = mean_dist, fill = group),
               width = 0.7, colour = "black") +

      geom_errorbar(data = d, aes(x = group, ymin = mean_dist - se_dist,
                                  ymax = mean_dist + se_dist),
                    width = 0.18, linewidth = 0.4, colour = "black") +
      scale_fill_manual(values = setNames(
        c("#ffd460", "#f07b3f", "#ea5455", "#2d4059"), levels(dat$group))) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.24)), breaks = .brk_y, labels = .fmt_y) +
      ylab(expression(atop("Mean pairwise", "Euclidean distance")))+
      xlab("")+
      labs(fill = "Propagule density(%)")+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+

      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      theme(axis.text.y = element_text(size = 6, color = "black"))+
      theme(axis.title.y = element_text(size = 8, color = "black"))+
      theme(axis.title.x = element_text(size = 8, color = "black"))+
      theme(axis.text.x = element_text(size = 6, color = "black", angle = 45, hjust = 1))+
      theme(legend.text = element_text(size = 8, color = "black"))+
      theme(legend.title = element_text(size = 9, color = "black"))+
      theme(legend.position = "none")

    if (!is.null(eu_test_df)) {
      .st <- eu_test_df[eu_test_df$ecosystem == eco &
                        eu_test_df$term == "RS (propagule density)", ]
      if (nrow(.st) == 1 && !is.na(.st$sym[1])) {
        .lab_st <- sprintf("RS%s", ifelse(.st$sym[1] == "ns", " ns", .st$sym[1]))
        p_sub <- p_sub + annotate("text", x = 0.55, y = Inf, label = .lab_st,
                                  hjust = 0, vjust = 1.45, size = 2.109, colour = "black")
      }
    }

    if (!is.null(eu_letters_df)) {
      .ltl <- eu_letters_df[eu_letters_df$ecosystem == eco, ]
      if (nrow(.ltl) >= 1) {
        .ltl <- .ltl[match(levels(dat$group), .ltl$density), ]
        .mi <- match(as.character(.ltl$density), as.character(d$group))
        stopifnot(!anyNA(.mi), nrow(.ltl) == 4)
        .ym <- max(d$mean_dist + d$se_dist, na.rm = TRUE)
        .ltl$y <- d$mean_dist[.mi] + d$se_dist[.mi] + 0.04 * .ym
        stopifnot(all(is.finite(.ltl$y)))
        p_sub <- p_sub + annotate("text", x = as.character(.ltl$density), y = .ltl$y,
                                  label = .ltl$letters, size = 1.4, vjust = 0)
      }
    }
    outputname <- paste0(mk, eco, "euclidean_distance.pdf")
    ggsave(outputname, p_sub, width = 3.0, height = 4.30, path = OUT,
           units = "cm", useDingbats = FALSE)
    if (!exists(".figS5_leg_dens")) {
      .save_legend(p_sub, "FigS5_density_legend", OUT)
      .figS5_leg_dens <- TRUE
    }
  }

  for (eco in .eco_levels(dat$ecosystem)) {
    data_eco <- subset(dat, ecosystem == eco)
    for (satru in levels(data_eco$group)) {
      sub_phi <- subset(data_eco, group == satru)
      if (nrow(sub_phi) < 3 || length(unique(sub_phi$group1)) < 2) next
      mod1 <- aov(copynum ~ group1, data = sub_phi)
      re   <- tryCatch(LSD.test(mod1, "group1", alpha = 0.05), error = function(e) NULL)
      if (is.null(re) || is.null(re$groups)) next
      re1  <- re$groups
      re1$marker     <- mk
      re1$eco        <- eco
      re1$satruation <- satru
      re1$treat      <- rownames(re1)
      lsd_all[[length(lsd_all) + 1]] <- re1
    }
  }
}

cv_values <- do.call(rbind, cv_values_all)
write.csv(cv_values, file.path(STATS, "figS5_cv_values.csv"), row.names = FALSE)

cv_test <- if (length(cv_test_all)) do.call(rbind, cv_test_all) else NULL
if (!is.null(cv_test)) {
  write.csv(cv_test, file.path(STATS, "figS5_cv_test.csv"), row.names = FALSE)
}
cat("[legacy] CV csv:", nrow(cv_values), "rows /",
    if (is.null(cv_test)) 0 else nrow(cv_test), "rows (no CV panel in the figure)\n")

if (length(lsd_all)) {
  lsd_df <- do.call(rbind, lsd_all)
  write.csv(lsd_df, file.path(STATS, "figS5_qpcr_multiple_comparisons.csv"), row.names = FALSE)
}

if (length(st_all)) {
  st_df <- do.call(rbind, st_all)
  st_df$p_disp <- .fmtp(st_df$p)
  write.csv(st_df, file.path(STATS, "figS5_qpcr_ST_main_effects.csv"), row.names = FALSE)
  cat("[ST] figS5_qpcr_ST_main_effects.csv:", nrow(st_df), "rows (6 panels x RS/T/RS x T)\n")
}

cap <- c(
  "Fig. S5. Prokaryotic (16S rRNA) and fungal (ITS) gene copy numbers and the",
  "within-group pairwise Euclidean distance of standardised copy numbers across",
  "propagule density levels and ecosystems.",
  "(a) Gene copy numbers are log10-transformed and expressed as copies per gram",
  "of dry soil (copies g-1 dry soil); error bars denote standard error of the mean",
  "(SE; n = 4). Panels are arranged with columns for prokaryotes (16S rRNA) and",
  "fungi (ITS) and rows for Croplands, Forests and Grasslands. Within each panel,",
  "RS and T indicate the main effects of propagule density (RS, resource surplus)",
  "and treatment (T) from two-way ANOVAs on log10-transformed copy numbers",
  "(density x treatment model, including the interaction term; the density x",
  "treatment interaction is reported in Table S7), and symbols denote",
  "significance (*P < 0.05, **P < 0.01, ***P < 0.001; ns, not significant);",
  "lower-case letters denote significant differences among the six treatments",
  "within each propagule density level (LSD test, P < 0.05). Fungal abundance is",
  "reported as fungal ITS gene copy number, which is subject to among-taxon",
  "variation in ITS copy number per genome and should be interpreted as a",
  "relative rather than absolute measure of fungal biomass. (b) Within each",
  "ecosystem, copy numbers of each marker were min-max standardised to [0, 1],",
  "and the pairwise Euclidean distance was computed among the four replicates",
  "within each propagule density x treatment group. Bars show the mean of the",
  "six treatment-level mean pairwise distances for each density; error bars",
  "denote SE (n = 6) and grey points show treatment-level values. RS and T",
  "indicate the main effects of propagule density and treatment from two-way",
  "ANOVAs on the per-replicate mean distance to the other replicates of the",
  "same group (density x treatment model; the interaction term is reported in",
  "Table S7), with symbols as above. Exact test statistics, degrees of freedom,",
  "P values and partial eta-squared for all panels are reported in Table S7.",
  "Treatments: CK, control; C, cellulose; N, nitrogen; pH, acidification;",
  "W, warming; D, drought."
)
writeLines(cap, file.path(CACHE, "figS5_caption_draft.txt"))

message("FigS5 outputs in ", OUT, ":")
print(list.files(path = OUT, pattern = "\\.pdf$", full.names = TRUE))
