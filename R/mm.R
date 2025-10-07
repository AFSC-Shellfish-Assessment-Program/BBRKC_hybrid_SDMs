#' Fits diffusion-taxis model to tagging data.
#'
#' @param formula Model formula for the habitat preference function that should be
#' specified without an intercept. This formula is used in the basis expansion step
#' performed by `mgcv` internally and so will accept `mgcv` spline notation. For
#' example, both `~ 0 + s(x, k = 3) + te(y,z, k = 4)` and `~0 + x + y + z` would be valid.
#' Note that no penalty on basis function coefficients is applied in this version of
#' the model. To fit a diffusion-only model, set `formula = NULL` (default) and `move_comps = "diffusion"`.
#' @param data Data frame specifying environmental variables at the time of tag release (pop-up).
#' Used in basis expansion and so should include all variables specified in habitat
#' preference function. Should also include a time variable and environmental variables
#' at all locations should be specified for each of these time steps.
#' @param time Name of time variable in `data` used to group tags released within the same time
#' windows (e.g. `"year"`). Should be specified even if all tags were released during one time window.
#' @param gridded_domain `sfc`-class variable where rows are polygons specifying
#' the study region.
#' @param cellsize Numeric length of one side of the polygons within `gridded_doman`.
#' @param deployment_locs Integer vector specifying the indices of `gridded_domain` where
#' tags were deployed.
#' @param release_locs Integer vector specifying the indices of `gridded_domain` where
#' tags released from animals aka "popped-up".
#' @param tags_per_step Numeric vector specifying the number of tags deployed during each time window specified by `time`.
#' @param move_comps Defaults to `"diffusion"` for fitting a diffusion-only model when `formula = NULL`. To
#' also model taxis, set `move_comps = c("diffusion", "taxis")` while ensuring that `formula != NULL` and `data != NULL`.
#' @param d_scaling Spatial scale discretization method. Accepts `"none"` or `"scale"`.
#'
#' @return A list with the following components:
#' \itemize{
#'   \item `opt` - `nlminb` output.
#'   \item `obj` - TMB model object created by `RTMB::MakeADFun`.
#'   \item `AIC` - AIC for fitted model.
#'   \item `model_matrix` - The design matrix derived from the `mgcv` basis expansion.
#'   \item `A` - Adjacency matrix derived from `gridded_domain`.
#'   \item `formula` - Habitat preference model formula.
#'   \item `preference_model` - Habitat preference model output from basis expansion. Used in prediction.
#'   \item `sd` - Estimates and standard errors for parameters from `RTMB`.
#' }
#'
#' @export
#'

mm <- function(formula = NULL,
               data = NULL,
               gridded_domain,
               deployment_locs,
               release_locs,
               time,
               tags_per_step,
               move_comps = "diffusion",
               cellsize,
               d_scaling = "scale"){

  # create adjacency matrix
  st_rook <- function(m, ...) sf::st_relate(m, m, pattern="F***1****", ... )
  grid_A <- st_rook(gridded_domain, sparse = TRUE)
  A <- methods::as(grid_A, "sparseMatrix") |> methods::as("TsparseMatrix")

  # indices of non-zero elements of adjacency matrix
  At_zz = cbind( attr(A,"i"), attr(A,"j") ) + 1

  # Number of tags
  n_tags <- length(deployment_locs)

  # time steps
  times <- sort(unique(data[[time]]))
  n_t <- length(times)
  n_s <- nrow(A)

  if (!is.null(formula) & identical(move_comps, "diffusion")){
    stop("Habitat preference formula provided for diffusion-only model")
  }

  if (is.null(formula) & all(c("diffusion","taxis") %in% move_comps)){
    stop("No formula provided for habitat preference model while taxis is selected.")
  }

  if (!is.null(formula) & all(c("diffusion","taxis") %in% move_comps) & is.null(data)){
    stop("No gridded environmental data provided to model habitat preference.")
  }

  # spline basis expansion---
  if (!is.null(formula)){
    f_mod <- mgcv::gam(formula = as.formula(paste(c("y", as.character(formula)), collapse=" ")),
    data = data)
    X_sz <- mgcv::predict.gam(f_mod, newdata = data, type = "lpmatrix")
  } else { # for diffusion-only case
    f_mod <- NULL
    X_sz <- NULL
  }

  # Input lists for RTMB
  dat = list(
    "A_ss" = A,                 # adjacency matrix
    "X_sz" = X_sz,              # basis expansion from mgcv
    "At_zz" = At_zz,            # indices of non-zero elements of adjacency matrix (indexed from 1)
    "grid0" = deployment_locs,  # deployment grid cells
    "grid1" = release_locs,     # release grid cells
    "n_tags" = tags_per_step,   # n tags per time step
    "n_s" = n_s,                # n grid cells within year
    "n_t" = n_t,                # n time steps
    "move_comps" = move_comps,  # selected movement components to model
    "delta_d" = cellsize,        # length of grid cell side
    "d_scaling" = d_scaling    # scale discretization
  )

  if (all(move_comps == "diffusion")) {
    par_list = list(
      "ln_D" = 1
    )
  } else {
    par_list = list(
      "ln_D" = 1,
      "beta_k" = rep(0, ncol(X_sz))
    )
  }

  # calculate joint negative log likelihood for RTMB
  f <- function(par_list) {
    RTMB::getAll(dat, par_list, warn = FALSE)
    jnll <- 0

    # ensures the basis expansion is indexed properly (by time variable)
    if (all(c("diffusion","taxis") %in% move_comps)) {
      basis_ind <- seq(1, nrow(X_sz))
      basis_ind_split <- split(basis_ind, sort(basis_ind %% n_t))
    } else if (all(move_comps == "diffusion")) {
      h_s <- NULL
    }

    tag_ind <- c(1L, cumsum(n_tags))
    for (ti in seq_len(n_t)) {

      # calculates habitat preference given some parameters
      if (all(c("diffusion","taxis") %in% move_comps)) {
        h_s <- X_sz[ unlist(basis_ind_split[ti]), ] %*% cbind(beta_k)
      }

      # calculates M_dot
      Mrate_gg <- M_dot(A_ss,
                        At_zz,
                        rate_par = ln_D,
                        pref_g = h_s,
                        move_comps,
                        d_scaling,
                        delta_d)

      # extract deployment and release grid indices
      if (ti == 1){
        dep_loc <- grid0[tag_ind[ti]:tag_ind[ti+1]]
        rel_loc <- grid1[tag_ind[ti]:tag_ind[ti+1]]
      } else {
        dep_loc <- grid0[(tag_ind[ti] + 1) : tag_ind[ti+1]]
        rel_loc <- grid1[(tag_ind[ti] + 1) : tag_ind[ti+1]]
      }

      # initial deployment locations
      V <- matrix(0, nrow = n_s, ncol = length(dep_loc))
      V[cbind(dep_loc, seq_along(dep_loc))] <- 1

      M <- Matrix::expm(t(Mrate_gg)) %*% V
      lik <- M[cbind(rel_loc, seq_along(rel_loc))]
      jnll <- jnll - sum(log(lik))
    }
    jnll
  }
  obj <- RTMB::MakeADFun(f, data = dat, par = par_list, silent = FALSE)
  opt <- nlminb(  start = obj$par,
                  objective = obj$fn,
                  gradient = obj$gr,
                  control = list( trace = 1))
  AIC <- 2 * opt$objective + 2 * length(opt$par)

  list(opt = opt,
       obj = obj,
       AIC = AIC,
       model_matrix = X_sz,
       A = A,
       formula = formula,
       preference_model = f_mod,
       sd = RTMB::sdreport(obj))
}


#' Internal function to generate the instantaneous movement rate matrix (M_dot).
#' Inputs are constructed internally to mm(). Largely copied from the RTMB version
#' of Spatio-temporal models for ecologists by Thorson and Kristensen 2024. Translated to
#' RTMB by Chris Cahill https://github.com/spacetime-ecologist/spacetime-ecologists-RTMB
#'
#' @keywords internal
#' @param A_ss Adjacency matrix of `gridded_domain`.
#' @param At_zz Indices of non-zero elements of adjacency matrix (indexed from 1).
#' @param rate_par Diffusion parameter.
#' @param pref_g Habitat preference vector.
#' @param move_comps `c("diffusion", "taxis")` or `"diffusion"`.
#' @param d_scaling Spatial scale discretization method.
#' @param delta_d `cellsize` parameter.
#' @param colsumA_g Column sums of adjacency matrix.


M_dot <- function(A_ss,
                  At_zz,
                  rate_par,
                  pref_g,
                  move_comps,
                  d_scaling,
                  delta_d) {
  Mrate_gg <- RTMB::AD(Matrix::Matrix(0, nrow(A_ss), ncol(A_ss)))
  ones <- matrix(1, ncol = 1, nrow = nrow(A_ss))

  # no scale discretization applied
  if (d_scaling == "none"){
    if (all(c("diffusion", "taxis") %in% move_comps)) {
      d_pref <- (pref_g[At_zz[,2]] - pref_g[At_zz[,1]])
      Mrate_gg[At_zz] <- Mrate_gg[At_zz] + exp(rate_par + d_pref)
    } else if (all(move_comps == "diffusion")) {
      Mrate_gg[At_zz] <- Mrate_gg[At_zz] + exp(rate_par)
    } else {
      stop("Unsupported move_comps")
    }

    # default, ensures Metzler matrix
  } else if (identical(d_scaling, "scale")) {
    D <- exp(rate_par)
    if (all(c("diffusion","taxis") %in% move_comps)) {

      Mrate_gg[At_zz] <- Mrate_gg[At_zz] +
        D / delta_d^2 * exp((pref_g[At_zz[, 2]] - pref_g[At_zz[, 1]]) / delta_d)

    } else if (all(move_comps == "diffusion")) {
      Mrate_gg[At_zz] <- Mrate_gg[At_zz] + D / (delta_d^2)
    } else stop("Unsupported move_comps")
  } else {
    stop("Unknown d_scaling")
  }
  row_sums <- Mrate_gg %*% ones # rowSums()
  diag(Mrate_gg) <- diag(Mrate_gg) - as.vector(row_sums)
  Mrate_gg
}

