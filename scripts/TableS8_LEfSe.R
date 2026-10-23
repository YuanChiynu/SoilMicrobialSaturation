pdf(NULL); dev.off()

DATA  <- "data"
OUT   <- file.path("output", "TableS8_LEfSe")
STATS <- file.path("output", "TableS8_LEfSe", "stats")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(STATS, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(microeco)
  library(magrittr)
  library(openxlsx)
})

SMOKE <- nzchar(Sys.getenv("LEFSE_SMOKE", unset = ""))
T0 <- Sys.time()
cat(sprintf("[start] %s  pid=%d  R=%s  microeco=%s  SMOKE=%s\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"), Sys.getpid(),
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(packageVersion("microeco")), SMOKE))

TREAT_LEVELS <- c("CK", "C", "N", "pH", "W", "D")
DENSITY_LEVELS   <- c("S10", "S40", "S70", "S100")
ECO_LEVELS   <- c("Cropland", "Forest", "Grassland")
LEVELS       <- c("phylum", "class", "order", "family", "genus", "species")

cap1 <- function(s) paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))

read_csv_clean <- function(path, row_names = NULL) {
  con <- file(path, open = "r", encoding = "UTF-8")
  on.exit(close(con))
  df <- read.csv(con, header = TRUE, check.names = FALSE,
                 row.names = row_names, stringsAsFactors = FALSE)
  names(df) <- sub("^\ufeff", "", names(df))
  df
}

prep_kingdom <- function(kingdom = c("16S", "ITS")) {
  kingdom <- match.arg(kingdom)
  dirn <- if (kingdom == "16S") "prokaryote" else "fungi"
  ot <- read_csv_clean(file.path(DATA, "data", dirn,
                           if (kingdom == "16S") "OTU_rare_16S_all.csv" else "OTU_rare_ITS.csv"),
                 row_names = 1)
  ot <- as.data.frame(t(ot))
  ot$names <- rownames(ot)
  grp <- read_csv_clean(file.path(DATA, "data", dirn, "group.csv"))
  stopifnot("sample" %in% names(grp))
  m <- merge(ot, grp, by.x = "names", by.y = "sample")
  m$ecosystem <- if (kingdom == "16S") {
    c(C = "Cropland", F = "Forest", G = "Grassland")[m$ecosystem]
  } else {
    c(Corpland = "Cropland", Forest = "Forest", Grass = "Grassland")[m$ecosystem]
  }
  m$density <- factor(paste0("S", sub("%re$", "", m$group)), levels = DENSITY_LEVELS)

  m$treatment <- factor(ifelse(m$group1 == "T", "W",
                        ifelse(m$group1 == "W", "D", m$group1)), levels = TREAT_LEVELS)
  stopifnot(!any(is.na(m$density)), !any(is.na(m$treatment)), !any(is.na(m$ecosystem)))

  otu_cols <- setdiff(colnames(ot), "names")
  otu <- t(as.matrix(m[, otu_cols, drop = FALSE]))
  colnames(otu) <- m$names
  sdat <- m[, c("names", "ecosystem", "density", "treatment")]

  tax_raw <- read.delim(file.path(DATA, "data", dirn, "metadata.txt"),
                        header = TRUE, check.names = FALSE, colClasses = "character",
                        fileEncoding = "UTF-8-BOM")
  names(tax_raw) <- sub("^\ufeff", "", names(tax_raw))
  body <- tax_raw
  pref_cols <- sapply(c(kingdom = "k__", phylum = "p__", taxclass = "c__",
                        order = "o__", family = "f__", genus = "g__", species = "s__"),
                      function(pf) {
                        for (i in seq_len(ncol(body))) {
                          hits <- sum(startsWith(body[[i]][seq_len(min(300, nrow(body)))], pf))
                          if (hits > 100) return(i)
                        }
                        if (pf == "k__") {
                          for (i in seq_len(ncol(body))) {
                            hits <- sum(startsWith(body[[i]][seq_len(min(300, nrow(body)))], "d__"))
                            if (hits > 100) return(i)
                          }
                        }
                        NA_integer_
                      })
  strip <- function(x) ifelse(startsWith(x, "s__") | startsWith(x, "g__") |
                               startsWith(x, "p__") | startsWith(x, "f__") |
                               startsWith(x, "c__") | startsWith(x, "o__") |
                               startsWith(x, "k__") | startsWith(x, "d__"),
                             substr(x, 4, nchar(x)), x)
  clean <- function(x) {
    x <- strip(trimws(x))
    x[!is.na(x) & x %in% c("unclassified", "NA", "")] <- NA
    x
  }
  cl_col <- pref_cols[["taxclass"]]; or_col <- pref_cols[["order"]]
  tax <- data.frame(
    id      = body[[1]],
    kingdom = clean(body[[pref_cols[["kingdom"]]]]),
    phylum  = clean(body[[pref_cols[["phylum"]]]]),
    class   = if (is.na(cl_col)) NA_character_ else clean(body[[cl_col]]),
    order   = if (is.na(or_col)) NA_character_ else clean(body[[or_col]]),
    family  = clean(body[[pref_cols[["family"]]]]),
    genus   = clean(body[[pref_cols[["genus"]]]]),
    stringsAsFactors = FALSE)
  has_species <- !is.na(pref_cols[["species"]])
  if (has_species) tax$species <- clean(body[[pref_cols[["species"]]]])
  rownames(tax) <- tax$id

  common <- intersect(rownames(otu), rownames(tax))
  stopifnot(length(common) == nrow(otu))
  tax <- tax[common, , drop = FALSE]
  list(otu = otu, sdat = sdat, tax = tax, has_species = has_species)
}

genus_label <- function(tax) {
  lab <- tax$genus
  lvl <- rep("genus", nrow(tax))
  need <- is.na(lab)
  lab[need] <- paste0(tax$family[need], "_uncl."); lvl[need] <- "family"
  need <- is.na(lab) | lab == "NA_uncl."
  lab[need] <- paste0(tax$phylum[need], "_uncl."); lvl[need] <- "phylum"
  need <- is.na(lab) | lab == "NA_uncl."
  lab[need] <- "unassigned"; lvl[need] <- "unassigned"
  data.frame(label = lab, lowest_rank = lvl, stringsAsFactors = FALSE)
}

species_label <- function(tax) {
  ok <- !is.na(tax$genus) & !is.na(tax$species) &
    !grepl("uncultur|unidentif|unclassified|sp\\.$", tax$species)
  lab <- rep(NA_character_, nrow(tax))
  lab[ok] <- paste(tax$genus[ok], tax$species[ok])
  data.frame(label = lab, lowest_rank = ifelse(ok, "species", NA_character_),
             stringsAsFactors = FALSE)
}

modal_tax <- function(tax_sub) {
  md <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_character_ else names(sort(table(x), decreasing = TRUE))[1]
  }
  c(phylum = md(tax_sub$phylum), class = md(tax_sub$class),
    order = md(tax_sub$order), family = md(tax_sub$family), genus = md(tax_sub$genus))
}

REFS <- data.frame(
  key = c("F07", "F12", "K16", "H17", "B13", "DV18", "D11", "V19", "H12",
          "PS16", "DS20", "SM15"),
  reference = c(
    "Fierer, N., Bradford, M. A., & Jackson, R. B. (2007). Toward an ecological classification of soil bacteria. Ecology, 88(6), 1354-1364.",
    "Fierer, N., Lauber, C. L., Ramirez, K. S., Zaneveld, J., Bradford, M. A., & Knight, R. (2012). Comparative metagenomic, phylogenetic and physiological analyses of soil microbial communities across nitrogen gradients. ISME Journal, 6(5), 1007-1017.",
    "Kielak, A. M., Barreto, C. C., Kowalchuk, G. A., van Veen, J. A., & Kuramae, E. E. (2016). The ecology of Acidobacteria: connecting their phylogenetic lineage with their ecological characteristics. Frontiers in Microbiology, 7, 738.",
    "Ho, A., Di Lonardo, D. P., & Bodelier, P. L. E. (2017). Revisiting life strategy concepts in environmental microbial communities. ISME Journal, 11(7), 1493-1506.",
    "Barnard, R. L., Osborne, C. A., & Firestone, M. K. (2013). Responses of soil bacterial and fungal communities to extreme desiccation and rewetting. ISME Journal, 7(11), 2043-2052.",
    "de Vries, F. T., Griffiths, R. I., Bailey, M., et al. (2018). Soil bacterial networks are less stable under drought than fungal networks. Nature Communications, 9, 3033.",
    "DeBruyn, J. M., Nixon, L. T., Fawaz, M. N., Johnson, A. M., & Radosevich, M. (2011). Global biogeography and quantitative seasonal dynamics of Gemmatimonadetes in soil. Applied and Environmental Microbiology, 77(17), 6295-6300.",
    "Vetrovsky, T., Kohout, P., Kopecky, M., et al. (2019). A meta-analysis of global fungal distribution reveals climate-driven patterns. Nature Communications, 10, 5142.",
    "Hatzenpichler, R. (2012). Diversity, physiology, and niche differentiation of ammonia-oxidizing archaea. Applied and Environmental Microbiology, 78(21), 7501-7510.",
    "Petters, S., Sollinger, A., Bengtsson, M. M., & Urich, T. (2016). The soil microbial food web revisited: Predatory myxobacteria as keystone taxa? ISME Journal, 11, 2665-2675.",
    "Davidov, Y., & Jurkevitch, E. (2004). Diversity and evolution of Bdellovibrio-and-like organisms (BALOs), reclassification of genus Bacteriovorax. Antonie van Leeuwenhoek, 86, 33-46.",
    "Spain, A. M., Krumholz, L. R., & Elshahed, M. S. (2009). Abundance, diversity and novelty of Planctomycetes in a pristine, acidic peatland. ISME Journal, 3, 129-130."
  ), stringsAsFactors = FALSE)

phylum_strategy <- list(
  Proteobacteria    = c("Copiotrophic (r-selected); fast growth on labile C", "F07;F12;H17"),
  Bacteroidota      = c("Copiotrophic (r-selected)", "F07;F12"),
  Actinobacteriota  = c("Copiotrophic; spore-forming, desiccation-tolerant", "F07;B13;DV18"),
  Firmicutes        = c("r-selected; endospore-forming; desiccation-tolerant", "B13;DV18"),
  Acidobacteriota   = c("Oligotrophic (K-selected); slow-growing", "F07;K16"),
  Gemmatimonadota   = c("Slow-growing; adapted to low soil moisture", "D11"),
  Verrucomicrobiota = c("Oligotrophic / slow-growing (K-selected)", "H17"),
  Chloroflexi       = c("Slow-growing (K-selected)", "H17"),
  Planctomycetota   = c("Slow-growing (K-selected); oligotrophic in pristine soils", "H17;SM15"),
  Crenarchaeota     = c("Ammonia-oxidizing archaea; chemolithoautotrophic", "H12"),
  Myxococcota       = c("Predatory (bacterivorous); oligotrophic, K-selected", "PS16;H17"),
  Bdellovibrionota  = c("Obligate predatory (bacterivorous); slow-growing", "DS20"),
  Desulfobacterota  = c("Anaerobic sulfate-reducing; anoxic microsites", "H17"),
  Nitrospirota      = c("Chemolithoautotrophic nitrite oxidizers; slow-growing", "H17"),
  Cyanobacteria     = c("Photoautotrophic; r-selected opportunists", "H17"),
  Patescibacteria   = c("Ultra-small, reduced-genome symbionts/parasites", "H17"),
  Firmicutes_A      = c("r-selected; endospore-forming", "B13"),
  Methylomirabilota = c("Methanotrophic; anaerobic methane oxidation", "H17"),
  Spirochaetota     = c("Fermentative; copiotrophic in anoxic microsites", "H17"),
  Armatimonadota    = c("Slow-growing (K-selected)", "H17"),
  Dependentiae      = c("Obligate parasitic/symbiotic; reduced genome", "H17"),
  Ascomycota        = c("Stress-tolerant / ruderal; many xerotolerant taxa", "DV18;V19"),
  Basidiomycota     = c("Recalcitrant-C decomposers; mycorrhizal", "V19"),
  Mortierellomycota = c("Ruderal soil saprotrophs", "V19"),
  Mucoromycota      = c("Fast-growing ruderal saprotrophs", "V19"),
  Chytridiomycota   = c("Aquatic/soil zoosporic; labile-C decomposers", "V19"),
  Rozellomycota     = c("Parasitic/symbiotic; poorly characterised", "V19"),
  Glomeromycota     = c("Arbuscular mycorrhizal symbionts (obligate biotrophs)", "V19")
)
genus_overrides <- list(
  "Bacillus|Paenibacillus|Sporosarcina|Lysinibacillus" = c("Endospore-forming, desiccation-tolerant", "B13;DV18"),
  "Streptomyces|Nocardia|Nocardioides|Micromonospora|Marmoricola" = c("Spore-forming actinobacteria; desiccation tolerant", "B13"),
  "Penicillium|Aspergillus|Talaromyces|Paecilomyces|Trichoderma"   = c("Xerotolerant moulds (Eurotiomycetes)", "DV18;V19"),
  "Cladosporium|Alternaria"                            = c("Xerotolerant, drought-associated (Dothideomycetes)", "DV18;V19"),
  "Nitrososphaera|Nitrosocosmicus|Nitrosotalea"        = c("Ammonia-oxidizing archaea; chemolithoautotrophic", "H12"),
  "Bradyrhizobium|Mesorhizobium|Rhizobium|Allorhizobium|Sinorhizobium|Ensifer" = c("Slow-growing N2-fixing root symbionts", "H17"),
  "Pseudomonas|Burkholderia|Ralstonia|Massilia|Cupriavidus|Variovorax|Comamonas|Phenylobacterium" = c("Copiotrophic (r-selected)", "F07;F12"),
  "Bdellovibrio|Haliangium|Phaselicystis|Corallococcus|Myxococcus|Sorangium" = c("Predatory (bacterivorous); slow-growing", "PS16;DS20"),
  "Nitrosopumilus|Nitrosarchaeum|Candidatus Nitrosocaldus" = c("Ammonia-oxidizing archaea", "H12"),
  "Nitrospira|Nitrobacter|Nitrolancea"                 = c("Chemolithoautotrophic nitrite oxidizers", "H17"),
  "Mortierella|Umbelopsis|Mucor|Rhizopus"              = c("Ruderal saprotroph", "V19"),
  "Trichoderma|Fusarium|Acremonium|Clonostachys"       = c("Ruderal soil saprotroph / mycoparasite", "V19"),
  "Cenococcum|Rhizopogon|Suillus|Laccaria|Inocybe|Tomentella|Cortinarius" = c("Ectomycorrhizal symbionts", "V19"),
  "Russula|Lactarius|Amanita|Tricholoma|Sebacina"      = c("Ectomycorrhizal / ericoid symbionts", "V19")
)
assign_strategy <- function(phylum, genus) {
  out <- lapply(seq_along(phylum), function(i) {
    g <- if (is.na(genus[i])) "" else genus[i]
    hit <- NULL; lvl <- NA_character_
    for (pat in names(genus_overrides)) {
      if (nzchar(g) && grepl(pat, g, ignore.case = TRUE)) {
        hit <- genus_overrides[[pat]]; lvl <- "genus"; break
      }
    }
    if (is.null(hit)) {
      ph <- if (is.na(phylum[i])) "" else phylum[i]
      hit <- phylum_strategy[[ph]]
      lvl <- if (is.null(hit)) NA_character_ else "phylum"
    }
    if (is.null(hit)) {
      ph <- if (is.na(phylum[i])) "unclassified" else phylum[i]
      hit <- c(sprintf("No published life-strategy assignment for %s; inferred from taxonomy only", ph), NA_character_)
      lvl <- "unassigned"
    }
    data.frame(strategy = hit[1], refs = hit[2], strategy_level = lvl,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

build_level_table <- function(pk, kingdom, level, eco) {
  ids <- pk$sdat$names[pk$sdat$ecosystem == eco]
  stopifnot(length(ids) == 96)
  otu <- pk$otu[, ids, drop = FALSE]
  tax <- pk$tax
  if (level == "genus") {
    gl <- genus_label(tax)
    lab <- gl$label; lrank <- gl$lowest_rank
  } else if (level == "species") {
    if (!pk$has_species) return(NULL)
    sl <- species_label(tax)
    lab <- sl$label; lrank <- sl$lowest_rank
    ok <- !is.na(lab)
    lab <- lab[ok]; lrank <- lrank[ok]
    otu <- otu[ok, , drop = FALSE]; tax <- tax[ok, , drop = FALSE]
  } else {

    lab <- tax[[level]]
    lrank <- rep(level, nrow(tax))
    ok <- !is.na(lab)
    lab <- lab[ok]; lrank <- lrank[ok]
    otu <- otu[ok, , drop = FALSE]; tax <- tax[ok, , drop = FALSE]
  }
  if (is.null(lab) || all(is.na(lab))) return(NULL)
  keep <- !is.na(lab)
  agg <- as.data.frame(rowsum(as.matrix(otu[keep, , drop = FALSE]), group = lab[keep]))
  stopifnot(nrow(agg) == length(unique(lab[keep])))

  labs_u <- rownames(agg)
  ph <- cls <- ordf <- fam <- gen <- rep(NA_character_, length(labs_u))
  lrk <- rep(NA_character_, length(labs_u))
  for (i in seq_along(labs_u)) {
    sel <- which(lab == labs_u[i])
    mt <- modal_tax(tax[sel, , drop = FALSE])
    ph[i] <- mt[["phylum"]]; cls[i] <- mt[["class"]]; ordf[i] <- mt[["order"]]
    fam[i] <- mt[["family"]]; gen[i] <- mt[["genus"]]
    lrk[i] <- lrank[sel][1]
  }

  prev <- rowMeans(agg > 0)
  keep2 <- prev >= 0.10
  agg <- agg[keep2, , drop = FALSE]
  taxmap <- data.frame(taxon = labs_u[keep2], phylum = ph[keep2], class = cls[keep2],
                       order = ordf[keep2], family = fam[keep2],
                       genus = gen[keep2], lowest_rank = lrk[keep2],
                       prevalence = prev[keep2], stringsAsFactors = FALSE)
  rownames(taxmap) <- taxmap$taxon

  tt <- data.frame(Kingdom = rep(if (grepl("16S", kingdom)) "Bacteria" else "Fungi",
                                nrow(taxmap)),
                   stringsAsFactors = FALSE)
  if (level %in% c("phylum", "class", "order", "family")) {
    chain_names <- c("phylum", "class", "order", "family")
    chain <- c(Phylum = "phylum", Class = "class", Order = "order", Family = "family")
    chain <- chain[1:match(level, chain_names)]
    for (r in seq_along(chain)) {
      col <- chain[[r]]
      tt[[names(chain)[r]]] <- taxmap[[col]]
      bad <- is.na(tt[[names(chain)[r]]])
      if (any(bad)) tt[[names(chain)[r]]][bad] <- taxmap$taxon[bad]
    }
  } else {
    tt$Phylum <- taxmap$phylum
    tt$Family <- taxmap$family
    tt$Genus  <- taxmap$genus
    tt$Genus[is.na(tt$Genus)] <- taxmap$taxon[is.na(tt$Genus)]
    if (level == "species") tt$Species <- taxmap$taxon
  }
  rownames(tt) <- taxmap$taxon
  list(abund = agg, taxmap = taxmap, tax_table = tt,
       tax_level_name = cap1(level))
}

run_one <- function(abund, sdat_eco, tax_table, level_name, group_col, tag) {
  st <- sdat_eco
  st$Group <- st[[group_col]]
  rownames(st) <- st$names
  stopifnot(!anyDuplicated(rownames(st)))
  mt <- microeco::microtable$new(otu_table = abund, sample_table = st,
                                 tax_table = tax_table, auto_tidy = TRUE)

  mt$cal_abund()

  set.seed(20261007)
  td <- microeco::trans_diff$new(dataset = mt, method = "lefse", group = "Group",
                                 alpha = 0.05, p_adjust_method = "fdr",
                                 taxa_level = level_name, remove_unknown = FALSE)
  res_diff <- td$res_diff
  cat(sprintf("  [%s] LEfSe significant taxa %d; LDA %.2f-%.2f\n", tag, nrow(res_diff),
              min(res_diff$LDA), max(res_diff$LDA)))

  ab <- as.data.frame(mt$taxa_abund[[level_name]])
  ab <- ab[, rownames(st), drop = FALSE]
  grp <- st$Group
  taxa <- rownames(ab)

  H <- P <- Df <- rep(NA_real_, length(taxa))
  for (i in seq_along(taxa)) {
    x <- as.numeric(ab[i, ])
    if (length(unique(x)) < 2 || length(unique(as.character(grp))) < 2) next
    kt <- suppressWarnings(stats::kruskal.test(x ~ grp))
    H[i] <- unname(kt$statistic); Df[i] <- unname(kt$parameter); P[i] <- kt$p.value
  }
  keep <- !is.na(H)
  taxa <- taxa[keep]; H <- H[keep]; Df <- Df[keep]; P <- P[keep]
  q <- stats::p.adjust(P, method = "fdr")
  kw <- data.frame(taxon = taxa, kw_H = H, kw_df = Df, kw_P = P, kw_q = q,
                   n_samples = ncol(ab), n_groups = length(unique(as.character(grp))),
                   stringsAsFactors = FALSE)

  if (nrow(res_diff) > 0) {
    mm <- merge(kw[, c("taxon", "kw_P")],
                data.frame(taxon = as.character(res_diff$Taxa),
                           lefse_P = as.numeric(res_diff$P.unadj), stringsAsFactors = FALSE),
                by = "taxon")
    if (nrow(mm) > 0) {
      d <- abs(mm$kw_P - mm$lefse_P)
      cat(sprintf("  [%s] KW vs LEfSe P agreement: n=%d, max|dP|=%.3g, exact (<1e-6) rate=%.1f%%\n",
                  tag, nrow(mm), max(d, na.rm = TRUE), 100 * mean(d < 1e-6, na.rm = TRUE)))
    }
  }

  .wilcox_fast <- function(x, g, g1, g2) {
    a <- x[g == g1]; b <- x[g == g2]
    n1 <- length(a); n2 <- length(b); N <- n1 + n2
    if (n1 == 0 || n2 == 0 || length(unique(c(a, b))) < 2)
      return(c(W = NA_real_, P = NA_real_))
    r <- rank(c(a, b))
    W <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
    nt <- table(r)
    sigma2 <- n1 * n2 / 12 * ((N + 1) - sum(nt^3 - nt) / (N * (N - 1)))
    if (!is.finite(sigma2) || sigma2 <= 0) return(c(W = W, P = NA_real_))
    mu <- n1 * n2 / 2
    corr <- sign(W - mu) * 0.5
    z <- (W - mu - corr) / sqrt(sigma2)
    c(W = W, P = 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE)))
  }
  cand <- kw$taxon[!is.na(kw$kw_q) & kw$kw_q < 0.05]
  gs <- levels(factor(grp))
  combos <- utils::combn(gs, 2, simplify = FALSE)

  if (length(cand) > 0) {
    nchk <- min(5L, length(cand))
    maxd <- 0
    for (i in seq_len(nchk)) for (cp in head(combos, 2L)) {
      x <- as.numeric(ab[cand[i], ])
      f <- .wilcox_fast(x, grp, cp[1], cp[2])
      wt <- suppressWarnings(stats::wilcox.test(x[grp == cp[1]], x[grp == cp[2]],
                                                exact = FALSE))
      maxd <- max(maxd, abs(f[["W"]] - unname(wt$statistic)),
                  abs(f[["P"]] - wt$p.value))
    }
    cat(sprintf("  [%s] pairwise Wilcoxon fast path vs stats::wilcox.test: max diff %.3g\n",
                tag, maxd))
    stopifnot(maxd < 1e-10)
  }
  pw_list <- vector("list", length(cand))
  for (i in seq_along(cand)) {
    x <- as.numeric(ab[cand[i], ])
    pw_list[[i]] <- do.call(rbind, lapply(combos, function(cp) {
      f <- .wilcox_fast(x, grp, cp[1], cp[2])
      data.frame(taxon = cand[i], pair = paste(cp, collapse = " vs "),
                 group1 = cp[1], group2 = cp[2], W = unname(f[["W"]]),
                 P = unname(f[["P"]]), stringsAsFactors = FALSE)
    }))
  }
  pw <- do.call(rbind, pw_list)
  if (is.null(pw)) {
    pw <- data.frame(taxon = character(0), pair = character(0), group1 = character(0),
                     group2 = character(0), W = numeric(0), P = numeric(0))
  }
  pw$q_BH <- stats::p.adjust(pw$P, method = "BH")
  list(res_diff = res_diff, kw = kw, pw = pw, abund = ab, groups = grp)
}

resdiff_all <- list(); kw_all <- list(); pw_all <- list(); abg_all <- list()
val_all <- list(); taxmap_all <- list()
k <- 0
for (kingdom in c("16S", "ITS")) {
  pk <- prep_kingdom(kingdom)
  lvls <- if (kingdom == "16S") setdiff(LEVELS, "species") else LEVELS
  for (level in lvls) {
    for (eco in ECO_LEVELS) {
      bt <- build_level_table(pk, kingdom, level, eco)
      if (is.null(bt)) next
      sdat_eco <- pk$sdat[pk$sdat$ecosystem == eco, , drop = FALSE]
      sdat_eco <- sdat_eco[match(colnames(bt$abund), sdat_eco$names), , drop = FALSE]
      stopifnot(all(sdat_eco$names == colnames(bt$abund)))
      taxmap_all[[length(taxmap_all) + 1]] <-
        cbind(data.frame(kingdom = kingdom, level = level, ecosystem = eco,
                         stringsAsFactors = FALSE), bt$taxmap)
      for (factor in c("density", "treatment")) {
        k <- k + 1
        tag <- sprintf("%s/%s/%s/%s", kingdom, level, eco, factor)
        cat(sprintf("[%d] %s ...\n", k, tag))
        rr <- run_one(bt$abund, sdat_eco, bt$tax_table, bt$tax_level_name, factor, tag)
        grpcol <- if (factor == "density") "density" else "treatment"

        rd <- as.data.frame(rr$res_diff, stringsAsFactors = FALSE)
        rd$kingdom <- kingdom; rd$level <- level; rd$ecosystem <- eco; rd$factor <- factor
        rd$taxon_path <- as.character(rd$Taxa)
        rd$taxon <- sub("^.*\\|", "", rd$taxon_path)
        rd <- merge(rd, rr$kw[, c("taxon", "kw_H", "kw_df", "kw_P", "kw_q")],
                    by.x = "taxon_path", by.y = "taxon", all.x = TRUE)
        rd <- merge(rd, bt$taxmap, by = "taxon", all.x = TRUE)
        st <- assign_strategy(rd$phylum, rd$genus)
        rd$strategy <- st$strategy; rd$strategy_level <- st$strategy_level
        rd$strategy_refs <- st$refs
        rd$group_enriched <- as.character(rd$Group)
        names(rd)[names(rd) == "Comparison"]   <- "comparison"
        names(rd)[names(rd) == "Significance"] <- "significance"

        rd <- rd[, c("kingdom", "level", "ecosystem", "factor", "comparison",
                     "taxon", "taxon_path", "phylum", "class", "order", "family", "genus",
                     "lowest_rank", "prevalence",
                     "group_enriched", "LDA", "kw_H", "kw_df", "kw_P", "kw_q",
                     "P.unadj", "P.adj", "significance", "strategy", "strategy_level",
                     "strategy_refs")]
        resdiff_all[[k]] <- rd

        kwd <- cbind(data.frame(kingdom = kingdom, level = level, ecosystem = eco,
                                factor = factor, stringsAsFactors = FALSE), rr$kw)
        names(kwd)[names(kwd) == "taxon"] <- "taxon_path"
        kwd$taxon <- sub("^.*\\|", "", kwd$taxon_path)
        kwd <- merge(kwd, bt$taxmap[, c("taxon", "phylum", "class", "order", "family",
                                        "genus", "lowest_rank", "prevalence")],
                     by = "taxon", all.x = TRUE)
        kwd <- kwd[, c("kingdom", "level", "ecosystem", "factor", "taxon", "taxon_path",
                       "phylum", "class", "order", "family", "genus", "lowest_rank", "prevalence",
                       "kw_H", "kw_df", "kw_P", "kw_q", "n_samples", "n_groups")]
        kw_all[[k]] <- kwd

        pwd <- cbind(data.frame(kingdom = kingdom, level = level, ecosystem = eco,
                                factor = factor, stringsAsFactors = FALSE), rr$pw)
        names(pwd)[names(pwd) == "taxon"] <- "taxon_path"
        pwd$taxon <- sub("^.*\\|", "", pwd$taxon_path)
        pwd <- pwd[, c("kingdom", "level", "ecosystem", "factor", "taxon", "taxon_path",
                       "pair", "group1", "group2", "W", "P", "q_BH")]
        pw_all[[k]] <- pwd

        abg <- do.call(rbind, lapply(levels(factor(rr$groups)), function(g) {
          sel <- rr$groups == g
          m <- rowMeans(rr$abund[, sel, drop = FALSE]) * 100
          s <- apply(rr$abund[, sel, drop = FALSE], 1, stats::sd) * 100
          n <- sum(sel)
          data.frame(taxon = names(m), group = g, N = n, Mean = as.numeric(m),
                     SD = as.numeric(s), SE = as.numeric(s) / sqrt(n),
                     stringsAsFactors = FALSE)
        }))
        abg <- cbind(data.frame(kingdom = kingdom, level = level, ecosystem = eco,
                                factor = factor, stringsAsFactors = FALSE), abg)
        names(abg)[names(abg) == "taxon"] <- "taxon_path"
        abg$taxon <- sub("^.*\\|", "", abg$taxon_path)
        abg <- merge(abg, bt$taxmap[, c("taxon", "phylum", "class", "order", "family",
                                        "genus", "lowest_rank")],
                     by = "taxon", all.x = TRUE)
        abg <- abg[, c("kingdom", "level", "ecosystem", "factor", "taxon", "taxon_path",
                       "phylum", "class", "order", "family", "genus", "lowest_rank",
                       "group", "N", "Mean", "SD", "SE")]
        abg_all[[k]] <- abg

        keep_sig <- !is.na(rd$P.adj) & rd$P.adj < 0.05
        sig_short <- rd$taxon[keep_sig]; sig_path <- rd$taxon_path[keep_sig]
        if (length(sig_path) > 0) {
          vv <- do.call(rbind, lapply(seq_along(sig_path), function(i) {
            data.frame(taxon = sig_short[i], sample = colnames(rr$abund),
                       group = as.character(rr$groups),
                       ra = as.numeric(rr$abund[sig_path[i], ]) * 100,
                       stringsAsFactors = FALSE)
          }))
          vv <- cbind(data.frame(kingdom = kingdom, level = level, ecosystem = eco,
                                 factor = factor, stringsAsFactors = FALSE), vv)
          val_all[[length(val_all) + 1]] <- vv
        }
        if (SMOKE) break
      }
      if (SMOKE) break
    }
    if (SMOKE) break
  }
  if (SMOKE) break
}

resdiff_all <- do.call(rbind, resdiff_all)
kw_all      <- do.call(rbind, kw_all)
pw_all      <- do.call(rbind, pw_all)
abg_all     <- do.call(rbind, abg_all)
val_all     <- if (length(val_all)) do.call(rbind, val_all) else NULL
taxmap_all  <- do.call(rbind, taxmap_all)

cat(sprintf("\nsummary: resdiff %d rows; KW %d rows; pairwise %d rows; group abundance %d rows; sample-level %s rows\n",
            nrow(resdiff_all), nrow(kw_all), nrow(pw_all), nrow(abg_all),
            if (is.null(val_all)) "0" else nrow(val_all)))

ev <- resdiff_all[!is.na(resdiff_all$P.adj) & resdiff_all$P.adj < 0.05, ]
ev <- do.call(rbind, lapply(split(ev, list(ev$kingdom, ev$level, ev$ecosystem, ev$factor)),
                            function(d) head(d[order(-d$LDA), ], 30)))

chk <- tapply(pw_all$pair, list(pw_all$kingdom, pw_all$level, pw_all$ecosystem, pw_all$factor),
              function(x) length(unique(x)))
cat("group pairs (expect density=6 / treatment=15):\n"); print(chk)

w <- function(df, name, dir = OUT) {
  write.csv(df, file.path(dir, name), row.names = FALSE, fileEncoding = "UTF-8")
  cat(sprintf("  -> %-42s %d rows\n", name, nrow(df)))
}
cat("\nCSV written:\n")
w(resdiff_all, "tableS8_lefse_resdiff_all.csv", STATS)
w(kw_all,      "tableS8_lefse_kw_stats.csv", STATS)
w(pw_all,      "tableS8_lefse_pairwise_wilcox.csv", STATS)
w(abg_all,     "tableS8_lefse_abund_by_group.csv")
if (!is.null(val_all)) w(val_all, "tableS8_lefse_values_long.csv")
w(ev,          "tableS8_lefse_evidence.csv", STATS)

cap <- c(
  "Table S8. Differential taxa from phylum to species level identified by LEfSe in three soil ecosystems under a global-change experiment.",
  "",
  "Linear discriminant analysis effect size (LEfSe; Segata et al. 2011) was run separately for each kingdom (bacterial/archaeal 16S rRNA gene, fungal ITS), each taxonomic level (phylum, class, order, family and genus for both kingdoms; species for ITS) and each ecosystem (Cropland, Forest, Grassland) on the rarefied amplicon counts, using the microeco R package (trans_diff, method = 'lefse').",
  "Two grouping factors were tested independently. (i) Initial propagule density: 10%, 40%, 70% and 100% (4 groups). (ii) Global-change treatment: CK (control), C (carbon), N (nitrogen), pH (acidification), W (warming) and D (drought) (6 groups). LEfSe is a multi-group comparison, so these are omnibus tests of differences among the levels of each factor; density x treatment interactions are not testable with LEfSe.",
  "For every taxon we report the Kruskal-Wallis omnibus statistic (kw_H, chi-square), its degrees of freedom (kw_df = k - 1, i.e. 3 for density and 5 for treatment), the corresponding p value (kw_P) and its Benjamini-Hochberg false-discovery-rate adjusted value (kw_q). For taxa significant in the omnibus test (kw_q < 0.05) we additionally report all pairwise Wilcoxon rank-sum sub-tests (W, P, q_BH) between the levels of the factor, following the LEfSe procedure. The linear discriminant analysis score (LDA, log10) is LEfSe's effect-size measure, and group_enriched is the group with the highest median abundance for that taxon.",
  "A taxon is reported as a differential taxon when its LEfSe-adjusted p value (lefse_q, BH-FDR) is < 0.05. 'strategy' assigns the documented growth/trophic strategy compiled from the cited literature at genus level where well established, otherwise inferred at phylum level ('strategy_level' = genus | phylum | unassigned); full citations are given in the 'references' sheet.",
  "OTUs were aggregated to the named genus where available, otherwise to the nearest reliable rank (family_uncl. / phylum_uncl.); ITS species-level results use named species only. At the phylum, class, order and family levels, OTUs were aggregated directly by the annotated rank and OTUs lacking the corresponding annotation were excluded. Taxa present in fewer than 10% of samples within an ecosystem were removed before analysis, matching the prevalence filter applied to the amplicon data.",
  "",
  sprintf("Generated: %s; microeco %s (R %s); LEfSe after Segata et al. (2011) Genome Biology 12:R60.",
          format(Sys.Date(), "%Y-%m-%d"), as.character(packageVersion("microeco")),
          paste(R.version$major, R.version$minor, sep = "."))
)
writeLines(cap, file.path(OUT, "tableS8_lefse_caption_draft.txt"))

methods <- c(
  "# Table S8 (LEfSe version): methods description",
  "",
  "> Purpose: descriptive paragraph for the newly added method.",
  "> Companion to the LEfSe differential-taxa figure and the life-strategy columns.",
  "> Draft text; manuscript assembly is outside this repository.",
  "",
  "## Methods paragraph",
  "",
  "To identify the specific taxa that respond to propagule density and to global-change treatments,",
  "we applied linear discriminant analysis effect size (LEfSe; Segata et al. 2011, Genome Biology 12:R60)",
  "as implemented in the microeco R package (v1.15.0; trans_diff, method = \"lefse\"). Analyses were run",
  "separately for each kingdom (bacterial/archaeal 16S rRNA gene and fungal ITS), each taxonomic level",
  "(phylum, class, order, family and genus for both kingdoms; species for ITS) and each ecosystem",
  "(cropland, forest, grassland), giving 66 independent runs. Within each run, LEfSe first applies a Kruskal-Wallis rank-sum test across the levels",
  "of the grouping factor (initial propagule density, four levels, or global-change treatment, six levels),",
  "then a pairwise Wilcoxon rank-sum sub-test on the taxa that pass the omnibus test, and finally computes",
  "the linear discriminant analysis (LDA) score on a logarithmic scale as the effect size. We report, for",
  "every taxon, the Kruskal-Wallis chi-square statistic with its degrees of freedom (k - 1) and p value,",
  "the Benjamini-Hochberg adjusted q values, the full set of pairwise Wilcoxon statistics (W, P, q) for",
  "the significant candidates, and the LDA score; a taxon is considered differentially abundant at",
  "q < 0.05. Because LEfSe is a multi-group method it cannot model interactions, so the density x",
  "treatment interaction is not tested here and no interaction term is reported.",
  "OTUs were aggregated to the named genus where available, otherwise to the nearest reliable rank",
  "(family- or phylum-level unclassified groups); ITS species-level results use named species only;",
  "at the phylum, class, order and family levels OTUs were aggregated directly by the annotated rank,",
  "and OTUs lacking the corresponding annotation were excluded.",
  "Taxa present in fewer than 10% of samples within an ecosystem were removed before analysis.",
  "Documented life-strategy and trophic assignments",
  "were compiled from the literature for each differential taxon (Table S8, columns strategy,",
  "strategy_level and strategy_refs; full citations in the references sheet).",
  "",
  "## Figure caption (Fig. S16)",
  "",
  "Fig. S16. Taxon-level responses to initial propagule density identified by LEfSe. Box plots show the",
  "relative abundance (%) of the five bacterial/archaeal 16S (top row) and fungal ITS (bottom row) genera",
  "with the highest LDA scores in each ecosystem (cropland, forest, grassland; columns). Within each panel,",
  "boxes are grouped by genus and coloured by initial propagule density (10%, 40%, 70% and 100%; n = 4",
  "replicates each). All displayed genera were identified as differentially abundant by LEfSe among density",
  "levels after Benjamini-Hochberg correction (q < 0.05); the Kruskal-Wallis chi-square statistic with its",
  "degrees of freedom (df = 3), exact p value, adjusted q value, LDA score and pairwise Wilcoxon sub-tests",
  "are given per taxon in Table S8 (LEfSe version).",
  "",
    sprintf("(generated %s; script: scripts/TableS8_LEfSe.R)", format(Sys.time(), "%Y-%m-%d %H:%M"))
)
writeLines(methods, file.path(OUT, "tableS8_lefse_methods_draft.md"), useBytes = TRUE)

README <- data.frame(
  Item = c("Purpose", "Method", "Grouping factors", "Taxonomic levels", "Omnibus test",
           "Degrees of freedom", "Pairwise sub-test", "Effect size", "Significance",
           "Taxonomic aggregation", "Prevalence filter", "Evidence columns",
           "Sheet guide", "Software"),
  Description = c(
    "Differential taxa from phylum to species level identified with LEfSe.",
    "LEfSe (linear discriminant analysis effect size; Segata et al. 2011) via microeco::trans_diff(method = 'lefse'); run separately for each kingdom x taxonomic level x ecosystem.",
    "Two independent factors: (i) initial propagule density, 4 levels (S10/S40/S70/S100); (ii) global-change treatment, 6 levels (CK/C/N/pH/W/D).",
    "16S: phylum, class, order, family, genus; ITS: phylum, class, order, family, genus, species (named species only).",
    "Kruskal-Wallis rank-sum test across the levels of the factor, per taxon (kw_H = chi-square statistic, kw_P = p value); BH-FDR across taxa within each run gives kw_q.",
    "kw_df = k - 1, i.e. 3 for the density factor (4 groups) and 5 for the treatment factor (6 groups); reported explicitly for every taxon.",
    "Pairwise Wilcoxon rank-sum sub-tests on the taxa with kw_q < 0.05, over all pairs of factor levels (6 pairs for density, 15 pairs for treatment); W, P, and BH-adjusted q_BH are reported.",
    "LDA (log10) score from LEfSe; group_enriched = the factor level with the highest median abundance.",
    "A taxon is called differential when lefse_q (BH-FDR adjusted LEfSe p) < 0.05. LEfSe is a multi-group method, so density x treatment interactions cannot be tested here.",
    "OTUs aggregated to the named genus where available, otherwise nearest reliable rank (family_uncl. / phylum_uncl.); ITS species level uses named species only; at phylum/class/order/family, OTUs aggregated directly by the annotated rank and unannotated OTUs excluded.",
    "Taxa present in < 10% of samples within an ecosystem were removed, matching the prevalence filter applied to the amplicon data (prv_cut = 0.10).",
    "'strategy' = documented growth/trophic strategy from the cited literature at genus level where established, otherwise inferred at phylum level ('strategy_level' = genus | phylum | unassigned); 'strategy_refs' keys refer to the 'references' sheet.",
    "LEfSe_16S_phylum ... LEfSe_16S_genus / LEfSe_ITS_phylum ... LEfSe_ITS_species = differential taxa per kingdom x level, with Kruskal-Wallis statistics and LDA. pairwise_wilcoxon = pairwise sub-tests. abundance_by_group = mean/SD/SE per taxon per group. evidence = top-30-by-LDA differential taxa per run with life-strategy and literature keys. references = cited sources. methods = narrative description of the method.",
    sprintf("R %s; microeco %s; LEfSe after Segata et al. (2011) Genome Biology 12:R60.",
            paste(R.version$major, R.version$minor, sep = "."),
            as.character(packageVersion("microeco")))
  ), stringsAsFactors = FALSE)

sh <- list("README" = README)
for (kd in c("16S", "ITS")) {
  for (lv in LEVELS) {
    if (kd == "16S" && lv == "species") next
    d <- resdiff_all[resdiff_all$kingdom == kd & resdiff_all$level == lv, , drop = FALSE]
    if (nrow(d) > 0) sh[[sprintf("LEfSe_%s_%s", kd, lv)]] <- d
  }
}
sh[["pairwise_wilcoxon"]]  <- pw_all
sh[["abundance_by_group"]] <- abg_all
sh[["evidence"]]           <- ev
sh[["references"]]         <- REFS
sh[["methods"]]            <- data.frame(methods_sheet = methods, stringsAsFactors = FALSE)
wb <- openxlsx::createWorkbook()
for (nm in names(sh)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, sh[[nm]])
  openxlsx::freezePane(wb, nm, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb, file.path(STATS, "TableS8_LEfSe.xlsx"), overwrite = TRUE)
cat("  -> TableS8_LEfSe.xlsx（", length(sh), " sheet）\n", sep = "")

cat(sprintf("\n[done] total %.1f min\n", as.numeric(difftime(Sys.time(), T0, units = "mins"))))
