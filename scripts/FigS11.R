DATA  <- "data"
OUT   <- file.path("output", "FigS11")
CACHE <- file.path("output", "FigS11")
STATS <- file.path("output", "FigS11", "stats")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
dir.create(STATS, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

library(vegan)
library(tidyverse)
library(ggpubr)
library(scales)

group_df <- read.csv(file.path(DATA, "amplicon/fungi/group.csv"), header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(group_df) <- group_df$sample

otulist <- list.files(path = file.path(DATA, "amplicon/fungi/by_ecosystem_density"), pattern = "\\.csv$", full.names = FALSE)
otulist <- setdiff(otulist, "group.csv")

group_order <- c("CK", "C", "N", "pH", "T", "W")

.relabel_treat_pm <- function(x) {
  m <- c(CK = "CK", C = "C", N = "N", pH = "pH", T = "W", W = "D")
  out <- m[as.character(x)]; ifelse(is.na(out), as.character(x), out)
}
.DISP6 <- c("CK", "C", "N", "pH", "W", "D")

pairwise_adonis_custom <- function(dist_mat, factors) {
  factors <- factor(factors, levels = group_order)
  groups <- levels(droplevels(factors))

  if (length(groups) < 2) {
    warning("fewer than 2 groups; file skipped.")
    return(NULL)
  }

  pairs <- combn(groups, 2)
  result_list <- list()
  idx <- 1

  for (i in 1:ncol(pairs)) {
    g1 <- pairs[1, i]
    g2 <- pairs[2, i]

    subset_idx <- factors %in% c(g1, g2)
    sub_dist <- as.dist(as.matrix(dist_mat)[subset_idx, subset_idx])
    sub_factor <- factor(factors[subset_idx], levels = c(g1, g2))

    ad <- adonis2(sub_dist ~ sub_factor, permutations = 999)
    p_val <- ad$`Pr(>F)`[1]
    r2_val <- ad$R2[1]

    result_list[[idx]] <- data.frame(Group1 = g1, Group2 = g2, P_val = p_val, R2 = r2_val)
    idx <- idx + 1
    result_list[[idx]] <- data.frame(Group1 = g2, Group2 = g1, P_val = p_val, R2 = r2_val)
    idx <- idx + 1
  }

  out <- do.call(rbind, result_list)
  out$Group1 <- factor(.relabel_treat_pm(out$Group1), levels = .DISP6)
  out$Group2 <- factor(.relabel_treat_pm(out$Group2), levels = .DISP6)
  out
}

pairwise_mean_bray <- function(dist_mat, factors) {
  factors <- factor(factors, levels = group_order)
  groups <- levels(droplevels(factors))

  if (length(groups) < 2) {
    warning("fewer than 2 groups; cannot compute mean Bray-Curtis distance.")
    return(NULL)
  }

  dm <- as.matrix(dist_mat)
  pairs <- combn(groups, 2)
  result_list <- list()
  idx <- 1

  for (i in 1:ncol(pairs)) {
    g1 <- pairs[1, i]
    g2 <- pairs[2, i]

    idx1 <- which(factors == g1)
    idx2 <- which(factors == g2)

    mean_dist <- mean(dm[idx1, idx2], na.rm = TRUE)

    result_list[[idx]] <- data.frame(Group1 = g1, Group2 = g2, Mean_Bray = mean_dist)
    idx <- idx + 1
    result_list[[idx]] <- data.frame(Group1 = g2, Group2 = g1, Mean_Bray = mean_dist)
    idx <- idx + 1
  }

  out <- do.call(rbind, result_list)
  out$Group1 <- factor(.relabel_treat_pm(out$Group1), levels = .DISP6)
  out$Group2 <- factor(.relabel_treat_pm(out$Group2), levels = .DISP6)
  out
}

perm_all <- list()
for (otu_file in otulist) {

  safe_name <- tools::file_path_sans_ext(basename(otu_file))
  safe_name <- gsub("%", "pct", safe_name)
  safe_name <- gsub("[\\\\/:*?\"<>|]", "_", safe_name)
  cat("processing: ", otu_file, "\n")

  otu_raw <- read.csv(file.path(DATA, "amplicon/fungi/by_ecosystem_density", otu_file), row.names = 1, check.names = FALSE)
  otu_t <- t(otu_raw)

  sample_names <- rownames(otu_t)
  group_idx <- match(sample_names, group_df$sample)

  if (any(is.na(group_idx))) {
    missing_samples <- sample_names[is.na(group_idx)]
    warning(paste0("samples absent from the grouping file were dropped: ", paste(missing_samples, collapse = ", ")))
  }

  valid_idx <- !is.na(group_idx)
  otu_t <- otu_t[valid_idx, , drop = FALSE]
  sample_group <- group_df$group1[group_idx[valid_idx]]
  sample_group <- factor(sample_group, levels = group_order)

  if (length(unique(na.omit(sample_group))) < 2) {
    warning(paste0("file ", otu_file, " fewer than 2 groups; skipped."))
    next
  }

  set.seed(123)
  min_depth <- min(rowSums(otu_t))
  otu_rare <- otu_t[, colSums(otu_t) > 0, drop = FALSE]

  bray_dist <- vegdist(otu_rare, method = "bray")

  plot_data <- pairwise_adonis_custom(bray_dist, sample_group)

  if (is.null(plot_data) || nrow(plot_data) == 0) {
    warning(paste0("file ", otu_file, " no plot drawn; skipped."))
    next
  }

  bray_mean_df <- pairwise_mean_bray(bray_dist, sample_group)
  plot_data <- left_join(plot_data, bray_mean_df, by = c("Group1", "Group2"))

  plot_data$Star <- case_when(
    plot_data$P_val < 0.001 ~ "***",
    plot_data$P_val < 0.01  ~ "**",
    plot_data$P_val < 0.05  ~ "*",
    TRUE ~ "ns"
  )

  perm_all[[otu_file]] <- data.frame(
    dataset = safe_name,
    Group1 = as.character(plot_data$Group1),
    Group2 = as.character(plot_data$Group2),
    P_val = plot_data$P_val,
    R2 = plot_data$R2,
    Mean_Bray = plot_data$Mean_Bray,
    Star = plot_data$Star,
    stringsAsFactors = FALSE
  )

  get_index <- function(g) match(g, group_order)
  upper_data <- plot_data %>% filter(get_index(Group1) < get_index(Group2))
  lower_data <- plot_data %>% filter(get_index(Group1) > get_index(Group2))

  p_heatmap <- ggplot(plot_data, aes(x = Group1, y = Group2)) +

    labs(x = "Treatment", y = "Treatment") +

    geom_tile(fill = "gray92", color = NA) +

    geom_hline(
      yintercept = seq_along(group_order),
      color = "white",
      linewidth = 1
    ) +
    geom_vline(
      xintercept = seq_along(group_order),
      color = "white",
      linewidth = 1
    ) +

    geom_point(
      data = lower_data,
      aes(size = Mean_Bray, color = R2),
      shape = 16
    ) +

    geom_text(
      data = upper_data,
      aes(label = Star),
      size = 10,
      vjust = 0.75
    ) +

    scale_color_gradient2(
      low = "#4D6AD0", mid = "white", high = "#FF9D7E",
      midpoint = 0.5, limits = c(0, 1),
      name = "R²", na.value = "gray80"
    ) +

    scale_size_continuous(
      range = c(8, 18),
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "Bray-Curtis"
    ) +

    scale_y_discrete(limits = group_order) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 7),
      plot.subtitle = element_text(hjust = 0.5, size = 7, face = "italic"),
      axis.text.x = element_text(size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "gray92", color = NA),
      plot.background = element_rect(fill = "gray92", color = NA),
      legend.position = "none",
      aspect.ratio = 1
    )

  out_name <- paste0(safe_name, "_PERMANOVA_heatmap.pdf")
  ggsave(out_name, p_heatmap, width = 6.5, height = 5.5, path = OUT)

  if (!exists(".permITs_leg")) {
    .save_legend(p_heatmap, "permanovaITS_legend", OUT)
    .permITs_leg <- TRUE
  }

  cat("written: ", out_name, "\n")
}

if (length(perm_all)) {
  perm_tab <- do.call(rbind, perm_all)
  write.csv(perm_tab, file.path(STATS, "permanovaITs_pairs_stats.csv"), row.names = FALSE)
  cat("statistics table: ", file.path(STATS, "permanovaITs_pairs_stats.csv"), nrow(perm_tab), "rows\n")
}

cat("all done.\n")
