if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("cowplot is required to export standalone legends")
}
library(cowplot)
library(grid)

.save_legend <- function(p, name, OUT,
                         orient = "bottom", w = NULL, h = NULL, dpi = 300) {
  if (is.null(p)) return(invisible(NULL))
  if (is.null(w)) w <- if (identical(orient, "bottom")) 14 else 4.5
  if (is.null(h)) h <- if (identical(orient, "bottom")) 1.6 else 12

  lp <- p + theme(
    legend.position       = orient,
    legend.justification  = "center",
    legend.box.just       = "center",
    legend.box.spacing    = unit(0.1, "cm"),
    legend.margin         = margin(0.1, 0.1, 0.1, 0.1, "cm")
  )
  leg <- tryCatch(get_legend(lp), error = function(e) NULL)

  if (is.null(leg)) {
    message("[legend] no legend found: ", name)
    return(invisible(NULL))
  }

  leg_w <- as.numeric(grid::convertWidth(grid::grobWidth(leg)  + unit(0.4, "cm"), "cm", valueOnly = TRUE))
  leg_h <- as.numeric(grid::convertHeight(grid::grobHeight(leg) + unit(0.4, "cm"), "cm", valueOnly = TRUE))
  if (!is.finite(leg_w) || leg_w <= 0) leg_w <- w
  if (!is.finite(leg_h) || leg_h <= 0) leg_h <- h

  base_path <- sub("\\.(pdf|png)$", "", name)
  cache_dir <- if (exists("CACHE")) CACHE else OUT
  pdf_path <- file.path(OUT, paste0(base_path, ".pdf"))
  png_path <- file.path(cache_dir, paste0(base_path, ".png"))

  ggsave(pdf_path, ggdraw(leg), width = leg_w, height = leg_h, units = "cm")
  ggsave(png_path, ggdraw(leg), width = leg_w, height = leg_h, units = "cm", dpi = dpi)
  message("[legend] legend written: ", name, " (", round(leg_w, 1), "cm x ", round(leg_h, 1), "cm)")
  invisible(NULL)
}

.save_fig_legend <- function(p, base, OUT, w = 9, h = 6, orient = "bottom",
                             leg_h = 3.2, dpi = 300, ...) {
  ggsave(file.path(OUT, paste0(base, ".pdf")),
         p + theme(legend.position = "none"),
         width = w, height = h, units = "cm", ...)
  .save_legend(p, base, OUT, orient = orient, h = leg_h, dpi = dpi)
  invisible(NULL)
}
