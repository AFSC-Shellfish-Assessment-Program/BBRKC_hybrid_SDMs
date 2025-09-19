# helper function for turning sf geometry into lat/lon columns
# Originally written by Josh London here -> https://github.com/r-spatial/sf/issues/231
sfc_as_cols <- function (x, geometry, names = c("x", "y"))
{
  if (missing(geometry)) {
    geometry <- sf::st_geometry(x)
  }
  else {
    geometry <- rlang::eval_tidy(enquo(geometry), x)
  }
  stopifnot(inherits(x, "sf") && inherits(geometry, "sfc_POINT"))
  ret <- sf::st_coordinates(geometry)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[, !names(x) %in% names]
  ret <- setNames(ret, names)
  dplyr::bind_cols(x, ret)
}

theme_fade <- function(...)
  {
  ggplot2::theme(strip.background = element_blank(), panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(), panel.background = element_blank(),
                 panel.border = element_rect(colour = "black", fill = NA,
                                             size = 0.2), legend.key = element_blank(), axis.title = element_text(size = 10))
  }

fct_to_num <-
  function (x)
  {
    as.numeric(levels(x))[x]
  }
