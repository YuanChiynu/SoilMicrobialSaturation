DATA  <- "data"

GRID  <- identical(Sys.getenv("GRID"), "1")
OUT   <- if (GRID) file.path("output", "Fig2", "grid") else file.path("output", "Fig2")
STATS <- file.path("output", "Fig2", "stats")
CACHE <- file.path("output", "Fig2")
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(STATS, showWarnings = FALSE, recursive = TRUE)
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

.grid_theme <- function() theme(
  panel.grid.major = element_line(colour = "grey88", linewidth = 0.25),
  panel.grid.minor = element_line(colour = "grey94", linewidth = 0.18))

source(file.path("scripts", "utils_legend.R"))

DENSITY_LEVELS   <- c("10re", "40re", "70re", "100re")
DENSITY_LABELS   <- c("10%", "40%", "70%", "100%")
TREAT_LEVELS <- c("CK", "C", "N", "pH", "T", "W")
TREAT_LABELS <- c("CK", "C", "N", "pH", "W", "D")

.p_star <- function(p) ifelse(is.na(p), "ns",
                         ifelse(p <= 0.001, "***",
                         ifelse(p <= 0.01,  "**",
                         ifelse(p <= 0.05,  "*", "ns"))))
PLOT_CM    <- 4.5

{
  library(permute)
  library(vegan)
  library(ggplot2)
  library(ggpubr)

.layer_under_errorbar <- function() structure(list(), class = "layer_under_errorbar")
ggplot_add.layer_under_errorbar <- function(object, plot, object_name) {
  n <- length(plot$layers)
  if (n >= 3) plot$layers <- c(plot$layers[1:(n - 2)], plot$layers[n], plot$layers[n - 1])
  plot
}
  library(agricolae)

  .perm_terms <- function(otu_mat, group_df) {
    meta <- group_df[match(rownames(otu_mat), as.character(group_df$sample)), , drop = FALSE]
    if (anyNA(meta$sample)) return(NULL)
    meta$group  <- factor(as.character(meta$group),  levels = DENSITY_LEVELS)
    meta$group1 <- factor(as.character(meta$group1), levels = TREAT_LEVELS)
    set.seed(20260915)
    a   <- vegan::adonis2(vegan::vegdist(otu_mat) ~ group * group1, data = meta,
                          permutations = 9999, by = "terms")
    rn  <- trimws(rownames(a))
    idx <- match(c("group", "group1", "group:group1"), rn)
    data.frame(term = c("density", "treatment", "interaction"),
               p_value = a[["Pr(>F)"]][idx], stringsAsFactors = FALSE)
  }

  .st_label <- function(perm) {
    if (is.null(perm) || !nrow(perm)) return(NA_character_)
    s <- perm$p_value[perm$term == "density"]
    t <- perm$p_value[perm$term == "treatment"]
    if (!length(s) || !length(t)) return(NA_character_)
    sprintf("RS%s\nT%s", .p_star(s[1]), .p_star(t[1]))
  }

  .lsd_tab <- function(d) {
    d$site  <- factor(as.character(d$site),  levels = c("Cropland", "Forest", "Grassland"))
    d$Group <- factor(as.character(d$Group), levels = DENSITY_LABELS)
    out <- do.call(rbind, lapply(levels(d$site), function(s) {
      x <- d[d$site == s & !is.na(d$Group), , drop = FALSE]
      if (nlevels(droplevels(x$Group)) < 2) return(NULL)
      m  <- stats::aov(similarity ~ Group, data = x)
      lt <- agricolae::LSD.test(m, "Group", alpha = 0.05)
      g  <- lt$groups
      mn <- tapply(x$similarity, x$Group, mean, na.rm = TRUE)
      se <- tapply(x$similarity, x$Group,
                   function(v) stats::sd(v, na.rm = TRUE) / sqrt(sum(is.finite(v))))
      data.frame(site = s, Group = rownames(g),
                 mean = as.numeric(mn[rownames(g)]),
                 se   = as.numeric(se[rownames(g)]),
                 letter = trimws(as.character(g[[ncol(g)]])),
                 stringsAsFactors = FALSE)
    }))
    out$site  <- factor(out$site,  levels = c("Cropland", "Forest", "Grassland"))
    out$Group <- factor(out$Group, levels = DENSITY_LABELS)
    out
  }

  .nmds_one <- function(otu_path, group_df) {
    otu <- read.csv(otu_path, header = TRUE, row.names = 1)
    otu <- t(otu)
    sumrow <- rowSums(otu)[1]
    otu <- (otu / sumrow)^0.5
    otu_mat <- as.matrix(otu)
    set.seed(20260915)
    nmds1 <- metaMDS(as.data.frame(otu_mat), autotransform = FALSE)
    perm  <- .perm_terms(otu_mat, group_df)
    sample_site <- data.frame(nmds1$point)[1:2]
    sample_site$names <- rownames(sample_site)
    names(sample_site)[1:2] <- c('NMDS1', 'NMDS2')
    sample_site <- merge(sample_site, group_df, by.x = 'names', by.y = 'sample', all.x = TRUE)
    sample_site$group  <- factor(sample_site$group,  levels = DENSITY_LEVELS,   labels = DENSITY_LABELS)
    sample_site$group1 <- factor(sample_site$group1, levels = TREAT_LEVELS, labels = TREAT_LABELS)
    list(site = sample_site, stress = nmds1$stress, perm = perm)
  }

  my_theme <- theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))
  .legend_done_nmds <- FALSE
  for (kingdom in c("16S", "ITS")) {
    if (kingdom == "16S") {
      group  <- read.csv(file.path(DATA, "amplicon/prokaryote/group.csv"), header = TRUE)
      otudir <- file.path(DATA, "amplicon/prokaryote/by_ecosystem")
    } else {
      group  <- read.csv(file.path(DATA, "amplicon/fungi/group.csv"), header = TRUE)
      otudir <- file.path(DATA, "amplicon/fungi/by_ecosystem")
    }
    otulist  <- list.files(path = otudir, pattern = "*.csv")

    cache_f <- file.path(CACHE, sprintf("nmds_cache_%s.rds", kingdom))
    if (!nzchar(Sys.getenv("RECOMPUTE")) && file.exists(cache_f)) {
      res_list <- readRDS(cache_f)
      cat(sprintf("[Fig2 NMDS] %s reuse cache %s (set RECOMPUTE=1 to recompute)\n", kingdom, basename(cache_f)))
    } else {
      res_list <- list()
      for (i in seq_along(otulist)) {
        name <- otulist[i]
        savename <- substr(name, 1, nchar(name) - 4)
        res_list[[savename]] <- .nmds_one(file.path(otudir, name), group)
      }
      saveRDS(res_list, cache_f)
      cat(sprintf("[Fig2 NMDS] %s recomputed and cached to %s\n", kingdom, basename(cache_f)))
    }

    .ellipse_extent <- function(x, y, level = 0.95) {
      ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
      if (length(x) < 5 || stats::sd(x) == 0 || stats::sd(y) == 0)
        return(list(x = range(x), y = range(y)))
      m  <- cbind(x, y); cm <- colMeans(m); S <- stats::cov(m); n <- nrow(m)
      r  <- sqrt(2 * stats::qf(level, 2, n - 1))
      th <- seq(0, 2 * pi, length.out = 200)
      eig <- eigen(S, symmetric = TRUE)
      E <- cbind(cos(th), sin(th)) %*% diag(sqrt(pmax(eig$values, 0))) %*% t(eig$vectors)
      E <- sweep(E * r, 2, cm, "+")
      list(x = range(E[, 1]), y = range(E[, 2]))
    }
    .env_extent <- function(site) {
      xs <- numeric(0); ys <- numeric(0)
      for (g in levels(site$group)) {
        d <- site[site$group == g, ]
        if (nrow(d) < 5) next
        e <- .ellipse_extent(d$NMDS1, d$NMDS2)
        xs <- c(xs, e$x); ys <- c(ys, e$y)
      }
      if (!length(xs)) { xs <- range(site$NMDS1); ys <- range(site$NMDS2) }
      list(x = range(xs), y = range(ys))
    }
    env <- lapply(res_list, function(r) .env_extent(r$site))
    gx <- range(unlist(lapply(env, `[[`, "x")))
    gy <- range(unlist(lapply(env, `[[`, "y")))
    cx <- mean(gx); cy <- mean(gy)
    span <- max(diff(gx), diff(gy)) * 1.04
    xlims <- c(cx - span / 2, cx + span / 2)
    ylims <- c(cy - span / 2, cy + span / 2)
    cat(sprintf("[Fig2 NMDS] %s unified SQUARE range: x=[%.3f, %.3f] y=[%.3f, %.3f]\n",
                kingdom, xlims[1], xlims[2], ylims[1], ylims[2]))

    for (savename in names(res_list)) {
      sample_site  <- res_list[[savename]]$site
      nmds1.stress <- res_list[[savename]]$stress

      st_txt   <- .st_label(res_list[[savename]]$perm)
      if (is.na(st_txt)) { st_txt <- "RSns\nTns"; cat("[warn] unmatched RS/T: ", kingdom, savename, "\n") }

      p <- ggplot(sample_site, aes(x=NMDS1, y=NMDS2, group = group)) +
        geom_point(aes(color = group, shape = group1), size = 1, alpha = 0.65) +

        scale_shape_manual(values = c(CK = 3, C = 15, N = 16, pH = 17, W = 18, D = 4)) +

        scale_color_manual(values = c("#ffd460","#f07b3f","#ea5455","#2d4059"))+
        stat_ellipse(data=sample_site,
                     geom = "polygon",level=0.95,
                     linetype = 1,size=0.8,
                     aes(fill=group),
                     alpha=0.1,
                     show.legend = T)+
        scale_fill_manual(values = c("#ffd460","#f07b3f","#ea5455","#2d4059"))+

        coord_fixed(ratio = 1, xlim = xlims, ylim = ylims, clip = "off") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(color = 'black', fill = 'transparent')) +

        annotate("text", x = xlims[1] + 0.030 * diff(xlims),
                 y = ylims[2] - 0.020 * diff(ylims),
                 label = st_txt, hjust = 0, vjust = 1, size = 2.1, lineheight = 0.95) +
        annotate("text", x = xlims[2] - 0.030 * diff(xlims),
                 y = ylims[2] - 0.020 * diff(ylims),
                 label = paste('Stress =', round(nmds1.stress, 3)),
                 hjust = 1, vjust = 1, size = 2.1) +

        labs(x = 'NMDS 1', y = 'NMDS 2',
             colour = 'Propagule density(%)', fill = 'Propagule density(%)', shape = 'Treatments') +
        theme(axis.text.x = element_text(size = 6,face = "plain"))+
        theme(axis.text.y = element_text(size = 6,face = "plain"))+
        theme(axis.title.x = element_text(size = 8, face = "plain"))+
        theme(axis.title.y = element_text(size = 8, face = "plain"))+
        theme(panel.border = element_rect(fill=NA,color="black", size=1.3, linetype="solid"))+

        theme(axis.ticks.x=element_line(color="black",size=1.1,lineend = "round"))+
        theme(axis.ticks.y=element_line(color="black",size=1.1,lineend = "round"))+
        theme(legend.position = "none")+
        my_theme

      if (GRID) p <- p + .grid_theme()
      p

      out_filename<-sapply(savename, function(x)paste(x,"NMDS.pdf",sep = ""))

      ggsave(p, file=file.path(OUT, out_filename), width = 5.75, height = 5.37, units = "cm")
      if (!.legend_done_nmds) {
        .save_legend(p, "Fig2_NMDS_legend", OUT)
        .legend_done_nmds <- TRUE
      }
    }
  }

  {
    library(permute)
    library(vegan)
    library(reshape)
    library(ggpubr)
    library(ggplot2)
    library(agricolae)
    library(data.table)

    for (marker in c("16S","ITS")) {
      if (marker == "16S") {
        otu_path  <- file.path(DATA, "amplicon/prokaryote/OTU_rare_16S_all.csv")
        grp1_path <- file.path(DATA, "amplicon/prokaryote/similarity_group1.csv")
        grp2_path <- file.path(DATA, "amplicon/prokaryote/similarity_group2.csv")
      } else {
        otu_path  <- file.path(DATA, "amplicon/fungi/OTU_rare_ITS.csv")
        grp1_path <- file.path(DATA, "amplicon/fungi/similarity_group1.csv")
        grp2_path <- file.path(DATA, "amplicon/fungi/similarity_group2.csv")
      }
      ado1 <- read.csv(otu_path, header = TRUE, row.names = 1)
      otusum<-unname(colSums(ado1)[1])

      ado=(t(ado1)/otusum)^0.5

      prok.dist <- vegdist(ado)

      prok.dist <- as.matrix(prok.dist)
      prok.dist <- as.data.frame(prok.dist)
      prok.dist$X1 <- rownames(prok.dist)

      prok.matr<- reshape::melt(prok.dist,id.vars = "X1", variable.name = "X2", value.name = "score")

      colnames(prok.matr)[2] <-"X2"
      colnames(prok.matr)[3] <-"similarity"
      prok.matr$similarity <- 1-prok.matr$similarity
      reselect <-prok.matr

      group1<-read.csv(grp1_path, header=TRUE)
      group2<-read.csv(grp2_path, header=TRUE)
      simi_dis <- merge(reselect, group1, by.x = 'X1',by.y = 'names', all.x = TRUE)
      simi_dis <- merge(simi_dis, group2, by.x = 'X2',by.y = 'names_2', all.x = TRUE)

      fin_dis_zhu<-subset(simi_dis,ecosystem_2==ecosystem)
      fin_dis_zhu<-subset(fin_dis_zhu,group_2==group)
      fin_dis_zhu<-subset(fin_dis_zhu,group1_2!=group1)
      colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="ecosystem")] <-"site"
      colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="group")] <-"Group"

      write.csv(fin_dis_zhu, file = file.path(CACHE, paste0('fin_dis_zhu_', marker, '.csv')))
    }

    if (file.exists(file.path(CACHE, 'fin_dis_zhu_ITS.csv'))) {
      file.copy(file.path(CACHE, 'fin_dis_zhu_ITS.csv'), file.path(CACHE, 'fin_dis_zhu.csv'), overwrite = TRUE)
    }

    colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="ecosystem")] <-"site"
    colnames(fin_dis_zhu)[which(colnames(fin_dis_zhu)=="group")] <-"Group"

    phi<-read.csv(file=file.path(CACHE, 'fin_dis_zhu_16S.csv'), header = T, row.names = 1)

phi$site = factor( phi$site , levels = c("C","F","G"), labels = c("Cropland","Forest","Grassland"))

    phi$Group = factor(phi$Group, levels = DENSITY_LEVELS, labels = DENSITY_LABELS)
    my_theme <- theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

    .lsd16 <- .lsd_tab(phi)
    .top16  <- max(.lsd16$mean + .lsd16$se, na.rm = TRUE)
    .lsd16$ypos <- .lsd16$mean + .lsd16$se + 0.020 * .top16
    .ymax16 <- .top16 * 1.22

    p_16s <- ggbarplot(phi, x="site", y="similarity", add = "mean_se",color = "black",fill = "Group",

                       palette = c("#ffd460","#f07b3f","#ea5455","#2d4059"), position = position_dodge(0.8))+
      ylab("Similarity")+
      xlab("Ecosystem")+

      labs(fill = "Propagule density(%)")+

      ylim(0, .ymax16)+

      stat_compare_means(aes(group=Group), method = "anova", label = "p.signif", label.y = .top16 * 1.12, size = 2.1)+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+
      theme(axis.text.y = element_text(size = 6,color="black"))+
      theme(axis.title.y = element_text(size = 8, color="black"))+
      theme(axis.title.x = element_text(size = 8, color="black"))+
      theme(legend.text = element_text(size = 8, color="black"))+
      theme(axis.text.x = element_text(size = 6,color="black",vjust = 1.0, hjust = 1.0))+
      theme(legend.position = "none")+
      my_theme

    p_16s <- p_16s +
      geom_text(data = .lsd16, aes(x = site, y = ypos, label = letter, group = Group),
                position = position_dodge(width = 0.8), size = 2.1, vjust = 0)
    p_16s

    ggsave('similarity_16s.pdf', p_16s, width = 5.50, height = 5.40, path = OUT, units = "cm")
    .save_legend(p_16s, "Fig2_similarity_16S_legend", OUT)

    phi<-read.csv(file=file.path(CACHE, 'fin_dis_zhu_ITS.csv'), header = T, row.names = 1)

phi$site = factor( phi$site , levels = c("C","S","T"), labels = c("Cropland","Forest","Grassland"))

    phi$Group = factor(phi$Group, levels = DENSITY_LEVELS, labels = DENSITY_LABELS)
    my_theme <- theme(plot.margin = unit(c(2, 2, 2, 2), "pt"))

    .lsdits <- .lsd_tab(phi)
    .topits  <- max(.lsdits$mean + .lsdits$se, na.rm = TRUE)
    .lsdits$ypos <- .lsdits$mean + .lsdits$se + 0.020 * .topits
    .ymaxits <- .topits * 1.22

    its_p <- ggbarplot(phi, x="site", y="similarity", add = "mean_se",color = "black",fill = "Group",

                       palette = c("#ffd460","#f07b3f","#ea5455","#2d4059"), position = position_dodge(0.8))+
      ylab("Similarity")+
      xlab("Ecosystem")+

      labs(fill = "Propagule density(%)")+
      ylim(0, .ymaxits)+
      stat_compare_means(aes(group=Group), method = "anova", label = "p.signif", label.y = .topits * 1.12, size = 2.1)+
      theme(panel.background = element_rect(fill = "transparent", colour = "black"))+
      theme(axis.text.y = element_text(size = 6,color="black"))+
      theme(axis.title.y = element_text(size = 8, color="black"))+
      theme(axis.title.x = element_text(size = 8, color="black"))+
      theme(legend.text = element_text(size = 8, color="black"))+
      theme(axis.text.x = element_text(size = 6,color="black",vjust = 1.0, hjust = 1.0))+
      theme(legend.position = "none")+
      my_theme
    its_p <- its_p +
      geom_text(data = .lsdits, aes(x = site, y = ypos, label = letter, group = Group),
                position = position_dodge(width = 0.8), size = 2.1, vjust = 0)
    its_p

    ggsave('similarity_its.pdf', its_p, width = 5.50, height = 5.40, path = OUT, units = "cm")
    .save_legend(its_p, "Fig2_similarity_ITS_legend", OUT)

    c("#f6ea8c","#f26d5b","#c03546","#492540","#8DD3C7","#CC6666","#CC6633","#CC6600","#CC33FF","#CC33CC","#3333FF","#3333CC","#333399","#333366","#003399")
  }

}
