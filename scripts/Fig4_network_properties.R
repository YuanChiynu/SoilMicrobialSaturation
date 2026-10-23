DATA  <- "data"
OUT   <- file.path("output", "Fig4")
STATS <- file.path(OUT, "stats")
if (!dir.exists(STATS)) dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
CACHE <- file.path("output", "Fig4")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "utils_legend.R"))

{
  library(Hmisc)
  library(psych)
  library(igraph)

  otu16s<-read.csv(file.path(DATA, "amplicon/prokaryote/OTU_rare_16S_all.csv"), header=TRUE, row.names=1)
  otu16s <- as.data.frame(t(otu16s))
  otu16s$names <- rownames(otu16s)
  group <- read.csv(file.path(DATA, "amplicon/prokaryote/group.csv"), header=TRUE)
  data_16s <- merge(otu16s,group,by.x ="names" ,by.y = "sample")
  rownames(data_16s) <- data_16s$mental
  data_16s <- data_16s[,c(2:44779)]
  colnames(data_16s) <- paste0(colnames(data_16s),"16s")

  metadata16s <- read.delim(file.path(DATA, "amplicon/prokaryote/metadata.txt"), row.name = 1, check.names = FALSE)
  rownames(metadata16s) <- paste0(rownames(metadata16s),"16s")

  otuits<-read.csv(file.path(DATA, "amplicon/fungi/OTU_rare_ITS.csv"), header=TRUE, row.names=1)
  otuits <- as.data.frame(t(otuits))
  otuits$names <- rownames(otuits)
  group <- read.csv(file.path(DATA, "amplicon/fungi/group.csv"), header=TRUE)
  data_its <- merge(otuits,group,by.x ="names" ,by.y = "sample")
  rownames(data_its) <- data_its$mental
  data_its <- data_its[,c(2:9518)]
  colnames(data_its) <- paste0(colnames(data_its),"its")

  metadataits <- read.delim(file.path(DATA, "amplicon/fungi/metadata.txt"), row.name = 1, check.names = FALSE)
  rownames(metadataits) <- paste0(rownames(metadataits),"its")

  data_16s$names <- rownames(data_16s)
  data_its$names <- rownames(data_its)
  multidomin <- merge(data_16s,data_its,by ="names")
  rownames(multidomin) <- multidomin$names
  multidomin <- multidomin[,-1]

  colnames(metadata16s) <- colnames(metadataits)
  taxo <- rbind(metadata16s,metadataits)

  colcounter <- ncol(multidomin)+1
  multidomin$names <- rownames(multidomin)
multidomin_cna <- merge(multidomin,group,by.x ="names" ,by.y = "mental")
  rownames(multidomin_cna) <- multidomin_cna$names

  resultdeprition <- c()
  for (ecosys in unique(multidomin_cna$ecosystem)) {
    tempdata <- subset(multidomin_cna,ecosystem==ecosys)
    for (satruation in unique(multidomin_cna$group)) {
      otu<-subset(tempdata,group==satruation)
      otu <- otu[,c(2:colcounter)]
      otu <- as.data.frame(t(otu))

      otu1 <- otu
      otu1[otu1>0] <- 1
      oturownum <- ncol(otu1)/2
      otu <- otu[which(rowSums(otu1) >= oturownum),]

      rowsum <- as.data.frame(rowSums(otu))
      sum <- colSums(rowsum)
      otu$RA <- rowsum/sum

      otu <- subset(otu[which(otu$RA>0.0001),], select = -RA)

      otu <- t(otu)

      relabund <- otu / rowSums(otu)
      otu <- scale(otu)

      rcorr_otu <- Hmisc::rcorr(as.matrix(otu), type = 'spearman')

      p <- rcorr_otu$P

      p <- p.adjust(p, method = 'BH')

      r <- rcorr_otu$r
      r[which(r < 0)]

      p[p>=0.05] <- -1
      p[p<0.05 & p>=0] <- 1
      p[p==-1] <- 0

      z <- r * p
      z[which(z < 0)]

      g <- graph.adjacency(z, weighted = TRUE, mode = 'undirected', diag = FALSE)
      g

      g <- delete.vertices(g, names(degree(g)[degree(g) == 0]))
      names(degree(g)[degree(g) == 0])

      E(g)$correlation <- E(g)$weight
      E(g)$weight <- abs(E(g)$weight)

      taxonomy <- taxo[as.character(V(g)$name),]
      V(g)$kingdom <- taxonomy$kingdom
      V(g)$phylum <- taxonomy$phylum
      V(g)$class <- taxonomy$class
      V(g)$order <- taxonomy$order
      V(g)$family <- taxonomy$family
      V(g)$genus <- taxonomy$genus
      V(g)$species <- taxonomy$species

      edge <- data.frame(as_edgelist(g))

      df <- as.data.frame(E(g)$correlation)
      df[df>0] <- 1
      df[df<0] <- -1
      colnames(df) <- c('cor')

      edge_list <- data.frame(
        source = edge[[1]],
        target = edge[[2]],
        weight = E(g)$weight,
        correlation = E(g)$correlation,
        cor = df
      )

      head(edge_list)
      out_filename<-paste0(ecosys,satruation, "_edge_screen0.0001_p0.05.csv")
      write.table(edge_list, file.path(CACHE, out_filename), sep = ',', row.names = FALSE, quote = FALSE)

      deprition <- paste0(ecosys,satruation)
      posedgenum <- nrow(subset(edge_list,cor==1))
      negedgenum <- nrow(edge_list)-posedgenum
      deprition <-cbind(deprition,posedgenum,negedgenum)

      node <- data.frame(
        id = names(V(g)),
        kingdom =V(g)$kingdom,
        phylum =V(g)$phylum,
        class =V(g)$class,
        order =V(g)$order,
        family =V(g)$family,
        genus =V(g)$genus,
        species =V(g)$species
      )
      out_filename<-paste0(ecosys,satruation, "_node_screen0.0001_p0.05.csv")
      write.table(node, file.path(CACHE, out_filename), sep = ',', row.names = FALSE, quote = FALSE)

      nodenum <- nrow(node)
      deprition <-cbind(deprition,nodenum)

      num.edges = length(E(g))
      num.edges

      num.vertices = length(V(g))
      num.vertices

      connectance = edge_density(g,loops=FALSE)
      connectance

      average.degree = mean(igraph::degree(g))
      average.degree

      average.path.length = average.path.length(g)
      average.path.length

      diameter = diameter(g, directed = FALSE, unconnected = TRUE, weights = NULL)
      diameter

      edge.connectivity = edge_connectivity(g)
      edge.connectivity

      clustering.coefficient = transitivity(g)
      clustering.coefficient
      no.clusters = no.clusters(g)
      no.clusters

      centralization.betweenness = centralization.betweenness(g)$centralization
      centralization.betweenness

      centralization.degree = centralization.degree(g)$centralization
      centralization.degree

      r_full <- r
      diag(r_full) <- 0
      r_full[is.na(r_full)] <- 0
      pos.mean <- function(v) { x <- v[v > 0]; if (!length(x)) 0 else mean(x) }
      neg.mean <- function(v) { x <- v[v < 0]; if (!length(x)) 0 else mean(x) }
      conn_pos <- apply(r_full, 2, pos.mean)
      conn_neg <- apply(r_full, 2, neg.mean)

      stopifnot(identical(names(conn_pos), colnames(relabund)))
      pos_cohesion <- as.vector(relabund %*% conn_pos)
      neg_cohesion <- as.vector(relabund %*% conn_neg)
      pos_cohesion.mean <- mean(pos_cohesion)
      pos_cohesion.sd   <- sd(pos_cohesion)
      neg_cohesion.mean <- mean(neg_cohesion)
      neg_cohesion.sd   <- sd(neg_cohesion)
      neg_pos_ratio <- ifelse(posedgenum > 0, negedgenum / posedgenum, NA_real_)
      modularity.value <- tryCatch(
        igraph::modularity(g, membership(igraph::cluster_louvain(g, weights = E(g)$weight))),
        error = function(e) NA_real_)
      degree.assortativity <- tryCatch(
        igraph::assortativity_degree(g, directed = FALSE),
        error = function(e) NA_real_)

      network_property <-c(num.edges=num.edges,num.vertices=num.vertices,connectance=connectance,
                           average.degree=average.degree,average.path.length=average.path.length,
                           diameter=diameter,edge.connectivity=edge.connectivity,clustering.coefficient=clustering.coefficient,
                           no.clusters=no.clusters,centralization.betweenness=centralization.betweenness,
                           centralization.degree=centralization.degree,
                           pos_cohesion.mean=pos_cohesion.mean,pos_cohesion.sd=pos_cohesion.sd,
                           neg_cohesion.mean=neg_cohesion.mean,neg_cohesion.sd=neg_cohesion.sd,
                           neg_pos_ratio=neg_pos_ratio,modularity=modularity.value,
                           degree_assortativity=degree.assortativity)
      network_property <- data.frame(network_property)
      out_filename<-paste0(ecosys,satruation, "_network_property.csv")
      write.csv(network_property, file.path(STATS, out_filename), row.names = T)
      resultdeprition <- rbind(resultdeprition,deprition)
    }

  }

  resultdeprition
  write.csv(resultdeprition, file = file.path(STATS, "edge_node_counts.csv"))

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

  eco_names  <- c("Corpland", "Forest", "Grass")
  eco_labs   <- c("Cropland", "Forest", "Grassland")
  density_names  <- c("10re", "40re", "70re", "100re")
  density_labs   <- c("10%", "40%", "70%", "100%")
  density_cols   <- c("#ffd460", "#f07b3f", "#ea5455", "#2d4059")

  set.seed(20260907)

  auc_rows  <- list()
  auc_sim_rows <- list()
  curve_rows <- list()

  for (ei in seq_along(eco_names)) {
    eco <- eco_names[ei]

    gs <- list()
    for (si in seq_along(density_names)) {
      fn <- file.path(CACHE, paste0(eco, density_names[si], "_edge_screen0.0001_p0.05.csv"))
      ed <- read.csv(fn, header = TRUE)
      gs[[si]] <- graph_from_data_frame(ed[, c("source", "target")], directed = FALSE)
    }

    rand_nc <- lapply(seq_along(gs), function(si) {
      ig <- gs[[si]]
      N <- vcount(ig)
      mat <- replicate(n_sim, {
        ord <- sample(N)
        nc_at_fractions(ig, ord, frac_grid)
      })
      mat
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

  auc_summary$ecosystem <- factor(auc_summary$ecosystem, levels = eco_labs)
  auc_summary$density <- factor(auc_summary$density, levels = density_labs)
  write.csv(auc_summary, file = file.path(STATS, "auc_summary.csv"), row.names = FALSE)
  write.csv(auc_sim_df, file = file.path(CACHE, "robustness_auc_sim.csv"), row.names = FALSE)
  write.csv(kw_stats,   file = file.path(STATS, "robustness_kruskal.csv"), row.names = FALSE)
  write.csv(curves_df,  file = file.path(CACHE, "robustness_curves.csv"), row.names = FALSE)

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
  annot_rand <- data.frame(
    ecosystem = auc_summary$ecosystem, density = auc_summary$density,
    panel = "Random attack",
    lab = sprintf("AUC = %.3f\np50 = %s",
                  auc_summary$random_AUC_mean, fmt_p50(auc_summary$random_p50_mean)))
  annot_targ <- data.frame(
    ecosystem = auc_summary$ecosystem, density = auc_summary$density,
    panel = "Sensitivity: degree-first attack",
    lab = sprintf("AUC = %.3f\np50 = %s",
                  auc_summary$targeted_AUC, fmt_p50(auc_summary$targeted_p50)))
  annot_df <- rbind(annot_rand, annot_targ)

  annot_df$ypos <- 0.99 - 0.075 * (as.numeric(annot_df$density) - 1)
  annot_df$panel <- factor(annot_df$panel, levels = c("Random attack", "Sensitivity: degree-first attack"))

  p_comp <- ggplot(d_comp, aes(x = frac, y = nrm, color = density)) +
    geom_ribbon(data = d_comp[!is.na(d_comp$lo), ],
                aes(ymin = lo, ymax = hi, fill = density), alpha = 0.18, color = NA) +
    geom_line(size = 1.05) +

    geom_line(data = quad_df, aes(x = frac, y = fit, color = density),
              linetype = "dashed", linewidth = 0.55, inherit.aes = FALSE) +
    scale_color_manual(values = density_cols) +
    scale_fill_manual(values = density_cols) +
    expand_limits(y = c(0, max(max(d_comp$hi, na.rm = TRUE) * 1.02, 1.30))) +
    geom_text(data = annot_df, aes(x = 0.55, y = ypos, label = lab, color = density),
              size = 2.4, hjust = 0, vjust = 1, show.legend = FALSE,
              inherit.aes = FALSE) +
    facet_grid(panel ~ ecosystem, scales = "free_y") +
    labs(x = "Proportions of removed nodes", y = "Normalized natural connectivity") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(),

          legend.position = "top", legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          strip.text = element_text(size = 7, face = "bold"))

ggsave(file = file.path(OUT, "Roboutness_composite.pdf"), p_comp, width = 6.3, height = 3.7)

  p_rand <- ggplot(d_rand, aes(x = frac, y = nrm, color = density)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = density), alpha = 0.18, color = NA) +
    geom_line(size = 1.1) +

    geom_line(data = subset(quad_df, panel == "Random attack"),
              aes(x = frac, y = fit, color = density),
              linetype = "dashed", linewidth = 0.55, inherit.aes = FALSE) +
    scale_color_manual(values = density_cols) + scale_fill_manual(values = density_cols) +
    facet_wrap(~ ecosystem, ncol = 3) +
    labs(x = "Proportions of removed nodes", y = "Normalized natural connectivity") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), legend.position = "none",
          strip.text = element_text(size = 7, face = "bold"))

  p_targ <- ggplot(d_targ, aes(x = frac, y = nrm, color = density)) +
    geom_line(size = 1.1) +

    geom_line(data = subset(quad_df, panel == "Sensitivity: degree-first attack"),
              aes(x = frac, y = fit, color = density),
              linetype = "dashed", linewidth = 0.55, inherit.aes = FALSE) +
    scale_color_manual(values = density_cols) +
    facet_wrap(~ ecosystem, ncol = 3) +
    labs(x = "Proportions of removed nodes", y = "Normalized natural connectivity") +
    theme_bw(base_size = 8) +
    theme(panel.grid.minor = element_blank(), legend.position = "none",
          strip.text = element_text(size = 7, face = "bold"))

  leg_p <- ggplot(curves_df[curves_df$attack == "random", ],
                  aes(x = frac, y = nrm, color = density)) +
    geom_line() + scale_color_manual(values = density_cols, name = "Propagule density (%)") +
    theme(legend.position = "right")
  if (exists(".save_legend")) {
    .save_legend(leg_p, "Fig4_Robustness_legend", OUT)
  }
}
