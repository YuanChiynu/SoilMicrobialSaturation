OUT   <- file.path("output", "Fig4")
CACHE <- OUT
dir.create(OUT,   showWarnings = FALSE, recursive = TRUE)

suppressMessages({ library(igraph); library(ggplot2) })

RECOMPUTE <- identical(Sys.getenv("RECOMPUTE"), "1")
SEED      <- 20261004

RADIAL_POW <- as.numeric(Sys.getenv("RADIAL_POW", unset = "0.5"))

SHAPE_MEAN_LO <- 0.637; SHAPE_MEAN_HI <- 0.697

phylum_cols <- c(
  "Acidobacteriota"  = "#FDD0B2",
  "Actinobacteriota" = "#C3804C",
  "Ascomycota"       = "#F5C76C",
  "Bacteroidota"     = "#1DB6EA",
  "Basidiomycota"    = "#333362",
  "Chloroflexi"      = "#92C44D",
  "Crenarchaeota"    = "#94C8E2",
  "Firmicutes"       = "#F3CA3B",
  "Gemmatimonadota"  = "#D1A5CC",
  "Mortierellomycota"= "#36445B",
  "Myxococcota"      = "#D5E6D8",
  "Planctomycetota"  = "#CE6F8F",
  "Proteobacteria"   = "#6A55A4",
  "Verrucomicrobiota"= "#D96B43"
)
other_col <- "#B8B8BB"
edge_cols <- c("1" = "#F8B4BF", "-1" = "#5D60AB")
all_phylum_cols <- c(phylum_cols, "Others" = other_col)

PT_OF_SIZE   <- 2.1345
PT_OF_STROKE <- 1.42

BORDER_RATIO <- 0.185
DARKEN       <- 0.88
EDGE_ALPHA   <- 0.22

CURVE_SAG    <- 0.12
NSEG         <- 5
PT_LINE      <- 2.134
LIM          <- 1 / 0.9888

PG_W        <- 453.0
PG_H        <- 358.0
STD_PT      <- 6.0
PANEL_REF_PT <- 267.0
PANEL_PT    <- 101.0
K_PT        <- PANEL_PT / PANEL_REF_PT
NODE_DIA_MIN <- 3.12 * K_PT
NODE_DIA_MAX <- 7.80 * K_PT
EDGE_W_LO    <- 0.15 * K_PT
EDGE_W_HI    <- 0.65 * K_PT

norm_phylum <- function(x) sub("^p__", "", x)

eco_names <- c("Corpland", "Forest", "Grass")
eco_labs  <- c("Cropland", "Forest", "Grassland")
density_names <- c("10re", "40re", "70re", "100re")
density_labs  <- c("10%", "40%", "70%", "100%")

.hexdark <- function(hex, k = DARKEN) {
  m <- grDevices::col2rgb(hex) * k
  m <- pmax(pmin(round(m), 255), 0)
  grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
}

.gephi_layout <- function(g, niter = 500, seed = SEED) {
  set.seed(seed)
  lay <- layout_with_fr(g, niter = niter)
  cx <- mean(lay[, 1]); cy <- mean(lay[, 2])
  dx <- lay[, 1] - cx;   dy <- lay[, 2] - cy
  r  <- sqrt(dx * dx + dy * dy)
  th <- atan2(dy, dx)
  n  <- length(r)

  rn <- max(r) * ((rank(r, ties.method = "average") - 0.5) / n)^(RADIAL_POW)
  cbind(rn * cos(th), rn * sin(th))
}

.disc_norm <- function(p) {

  mid <- (apply(p, 2, min) + apply(p, 2, max)) / 2
  p <- sweep(p, 2, mid)
  p / max(sqrt(rowSums(p^2)))
}

.circle_guard <- function(p, rim = 1.0) {

  rr <- sqrt(rowSums(p^2))
  i  <- which(rr > rim + 1e-12)
  if (length(i)) {
    p[i, ] <- p[i, ] * (rim / rr[i])
    cat(sprintf("  [circle_guard] clipped %d / %d nodes (>rim)\n", length(i), nrow(p)))
  }
  p
}

.bezier <- function(x, y, xe, ye, tt = seq(0, 1, length.out = NSEG)) {
  dx <- xe - x; dy <- ye - y
  ch <- sqrt(dx * dx + dy * dy); ch[ch < 1e-9] <- 1e-9
  h  <- 2 * CURVE_SAG * ch
  cxp <- (x + xe) / 2 + (-dy / ch) * h
  cyp <- (y + ye) / 2 + ( dx / ch) * h
  w0 <- (1 - tt)^2; w1 <- 2 * (1 - tt) * tt; w2 <- tt^2
  bx <- outer(x, w0) + outer(cxp, w1) + outer(xe, w2)
  by <- outer(y, w0) + outer(cyp, w1) + outer(ye, w2)
  list(x = as.vector(t(bx)), y = as.vector(t(by)), np = length(tt))
}

nodes_all <- list(); edges_all <- list(); diag_rows <- list()
t0 <- Sys.time()

for (ei in seq_along(eco_names)) {
  for (si in seq_along(density_names)) {
    tag <- paste0(eco_names[ei], density_names[si])
    nf  <- file.path(CACHE, paste0(tag, "_node_screen0.0001_p0.05.csv"))
    ef  <- file.path(CACHE, paste0(tag, "_edge_screen0.0001_p0.05.csv"))
    if (!file.exists(nf) || !file.exists(ef)) { cat("!! missing", tag, "\n"); next }

    node <- read.csv(nf, header = TRUE, stringsAsFactors = FALSE)
    edge <- read.csv(ef, header = TRUE, stringsAsFactors = FALSE)
    g <- graph_from_data_frame(edge[, c("source", "target")], directed = FALSE,
                               vertices = node$id)

    lrf <- file.path(CACHE, paste0("net_layout_", tag, ".rds"))
    if (!RECOMPUTE && file.exists(lrf)) {
      pos <- readRDS(lrf)
    } else {
      pos <- .circle_guard(.disc_norm(.gephi_layout(g)))
      saveRDS(pos, lrf)
    }
    stopifnot(nrow(pos) == vcount(g))

    deg <- as.numeric(degree(g))

    q <- (deg - min(deg)) / max(max(deg) - min(deg), 1e-9)
    dia  <- NODE_DIA_MIN + q * (NODE_DIA_MAX - NODE_DIA_MIN)
    bord <- BORDER_RATIO * dia
    size_aes   <- (dia - bord) / PT_OF_SIZE
    stroke_aes <- bord / PT_OF_STROKE

    m <- match(V(g)$name, node$id)
    phylum <- norm_phylum(node$phylum[m])
    phylum[is.na(phylum) | phylum == "" | phylum == "unclassified"] <- "Others"
    phylum[!phylum %in% names(all_phylum_cols)] <- "Others"

    fill_hex <- unname(all_phylum_cols[phylum])
    nodes_all[[length(nodes_all) + 1]] <- data.frame(
      x = pos[, 1], y = pos[, 2],
      phylum = phylum, degree = deg,
      fillh = fill_hex, bordh = .hexdark(fill_hex),
      siz = size_aes, strk = stroke_aes,
      ecosystem = eco_labs[ei], density = density_labs[si],
      stringsAsFactors = FALSE)

    si_ <- match(edge$source, V(g)$name); ti_ <- match(edge$target, V(g)$name)
    ex <- pos[si_, 1]; ey <- pos[si_, 2]; exe <- pos[ti_, 1]; eye <- pos[ti_, 2]
    w   <- abs(as.numeric(edge$correlation))
    ws  <- 0.1 + 0.9 * (w - min(w)) / max(max(w) - min(w), 1e-12)
    lw  <- (EDGE_W_LO + (EDGE_W_HI - EDGE_W_LO) * (ws - 0.1) / 0.9) / PT_LINE
    sgn <- ifelse(as.numeric(edge$cor) > 0, "1", "-1")
    ech <- unname(edge_cols[sgn])
    bz  <- .bezier(ex, ey, exe, eye)
    np  <- bz$np; ne <- nrow(edge)
    edges_all[[length(edges_all) + 1]] <- data.frame(
      x = bz$x, y = bz$y, id = rep(seq_len(ne), each = np),
      ech = rep(ech, each = np), lw = rep(lw, each = np),
      ecosystem = eco_labs[ei], density = density_labs[si],
      stringsAsFactors = FALSE)

    r <- sqrt(rowSums(pos^2))
    mean_r <- mean(r) / max(r); q50r <- unname(quantile(r, .5)) / max(r)
    q90r <- unname(quantile(r, .9)) / max(r)
    if (mean_r < SHAPE_MEAN_LO || mean_r > SHAPE_MEAN_HI)
      cat(sprintf("  [WARN shape] %s mean_r=%.3f outside disc range [%.3f, %.3f]\n",
                  tag, mean_r, SHAPE_MEAN_LO, SHAPE_MEAN_HI))
    diag_rows[[length(diag_rows) + 1]] <- data.frame(
      tag = tag, n = vcount(g), m = ecount(g),
      mean_r = mean_r, q50 = q50r, q90 = q90r,
      dia_min = min(dia), dia_max = max(dia),
      bord_min = min(bord), bord_max = max(bord),
      lw_min = min(lw) * PT_LINE, lw_max = max(lw) * PT_LINE)
    cat(sprintf("  %-14s n=%4d m=%7d  r:mean=%.3f q50=%.3f q90=%.3f\n",
                tag, vcount(g), ecount(g), mean_r, q50r, q90r))
  }
}

nodes_df <- do.call(rbind, nodes_all); rm(nodes_all)
edges_df <- do.call(rbind, edges_all); rm(edges_all)
diag_df  <- do.call(rbind, diag_rows)
write.csv(diag_df, file.path(OUT, "net_layout_diagnostics.csv"), row.names = FALSE)

nodes_df$ecosystem  <- factor(nodes_df$ecosystem,  levels = eco_labs)
nodes_df$density <- factor(nodes_df$density, levels = density_labs)
edges_df$ecosystem  <- factor(edges_df$ecosystem,  levels = eco_labs)
edges_df$density <- factor(edges_df$density, levels = density_labs)
cat(sprintf("\ncollected: nodes=%d  edge-points=%d (edges=%d)\n",
            nrow(nodes_df), nrow(edges_df), nrow(edges_df) / NSEG))

p <- ggplot() +
  geom_path(data = edges_df,
            aes(x = x, y = y, group = id, colour = ech, linewidth = lw),
            alpha = EDGE_ALPHA, lineend = "round") +
  geom_point(data = nodes_df,
             aes(x = x, y = y, fill = fillh, colour = bordh, size = siz, stroke = strk),
             shape = 21) +
  scale_colour_identity() + scale_fill_identity() +
  scale_linewidth_identity() + scale_size_identity() +
  facet_grid(ecosystem ~ density, switch = "y") +

  coord_fixed(xlim = c(-LIM, LIM), ylim = c(-LIM, LIM), expand = FALSE) +
  theme_void(base_size = STD_PT) +
  theme(strip.text.x      = element_text(size = STD_PT, margin = margin(b = 1.5)),
        strip.text.y.left = element_text(size = STD_PT, angle = 90, margin = margin(r = 1.5)),
        strip.placement   = "outside",

        panel.spacing.x   = unit(5  * K_PT / 0.378, "pt"),
        panel.spacing.y   = unit(13 * K_PT / 0.378, "pt"),
        plot.margin       = margin(2, 2, 2, 2))

ggsave(file.path(OUT, "Networks_composite.pdf"), p,
       width = PG_W / 72, height = PG_H / 72, units = "in")

LG_W <- 453.0; LG_H <- 40.0
NC <- 6; NR <- 3
taxa <- names(all_phylum_cols)
stopifnot(length(taxa) == 15, length(taxa) <= NC * NR - 2)

sw_pt <- 6.0
ux    <- LG_W / NC
uy    <- LG_H / NR

idx <- seq_along(taxa)
leg <- data.frame(
  cx  = ((idx - 1) %/% NR) + 0.5,
  cy  = NR - ((idx - 1) %% NR) - 0.5,
  lab = taxa, col = unname(all_phylum_cols[taxa]), stringsAsFactors = FALSE)
leg$tcx <- leg$cx - 0.34 + (sw_pt / ux) / 2
leg$lx  <- leg$cx - 0.34 + sw_pt / ux + 3.0 / ux
lnk <- data.frame(cy = c(2.5, 1.5), lab = c("Positive links", "Negative links"),
                  col = c(edge_cols[["1"]], edge_cols[["-1"]]), stringsAsFactors = FALSE)
lnk$tx <- 5.5 - 0.34; lnk$lx <- lnk$tx + sw_pt / ux + 3.0 / ux

p_leg <- ggplot() +
  geom_tile(data = leg, aes(x = tcx, y = cy, fill = col),
            width = sw_pt / ux, height = sw_pt / uy) +
  geom_text(data = leg, aes(x = lx, y = cy, label = lab), hjust = 0, size = STD_PT / .pt) +
  geom_segment(data = lnk, aes(x = tx, xend = tx + sw_pt / ux, y = cy, yend = cy,
                               colour = col),
               linewidth = 0.55, lineend = "butt") +
  geom_text(data = lnk, aes(x = lx, y = cy, label = lab), hjust = 0, size = STD_PT / .pt) +
  scale_fill_identity() + scale_colour_identity() +
  scale_x_continuous(limits = c(0, NC), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, NR), expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0),
                       panel.spacing = unit(0, "pt"))
ggsave(file.path(OUT, "Networks_legend.pdf"), p_leg,
       width = LG_W / 72, height = LG_H / 72, units = "in")

cat(sprintf("\nNetworks composite done  (%.1f min)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
cat(sprintf("  page: %.1f x %.1f pt  |  legend: %.1f x %.1f pt\n", PG_W, PG_H, LG_W, LG_H))
