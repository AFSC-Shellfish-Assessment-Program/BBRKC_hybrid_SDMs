# internal function to generate the instantaneous movement rate matrix (M_dot)
make_M <- function(A_ss,      # Adjacency matrix
                   At_zz,     # indices of non-zero elements of adjacency matrix (indexed from 1)
                   rate_par,  # diffusion parameter
                   pref_g,    # habitat preference vector
                   move_comps,# c("diffusion", "taxis) or "diffusion"
                   d_scaling, # Choose scale discretization method from "none", "scale1", or "scale2.
                              # "scale3" is default and ensures Metzler matrix
                              # (see chapter 10 in https://github.com/spacetime-ecologist/Edition_1
                              # and in the RTMB implementation https://github.com/spacetime-ecologist/spacetime-ecologists-RTMB)
                   DeltaD,    # length of grid cell size (in km)
                   colsumA_g  # adjacency matrix column sums
                   ) {
  M <- AD(Matrix(0, nrow(A_ss), ncol(A_ss)))

  # no scale discretization applied
  if (d_scaling == "none"){
    if (all(c("diffusion", "taxis") %in% move_comps)) {
      d_pref <- (pref_g[At_zz[,2]] - pref_g[At_zz[,1]])
      M[At_zz] <- M[At_zz] + exp(rate_par + d_pref)
    } else if (identical(move_comps, "diffusion")) {
      M[At_zz] <- M[At_zz] + exp(rate_par)
    } else {
      stop("Unsupported move_comps")
    }

  } else if (d_scaling == "scale1"){
    D <- exp(rate_par)
    if (all(c("diffusion", "taxis") %in% move_comps)) {

      M[At_zz] <- M[At_zz] + (D / colsumA_g[ At_zz[, 1] ] / DeltaD^2) +
        (pref_g[At_zz[, 2]] - pref_g[At_zz[, 1]]) / DeltaD

    } else if (identical(move_comps, "diffusion")) {

      M[At_zz] <- M[At_zz] + (D / colsumA_g[ At_zz[, 1] ] / DeltaD^2)

    } else {
      stop("Unsupported move_comps")
    }

  # default, ensures Metzler matrix
  } else if (identical(d_scaling, "scale2")) {
    D <- exp(rate_par)
    if (all(c("diffusion","taxis") %in% move_comps)) {

      M[At_zz] <- M[At_zz] +
        D / DeltaD^2 * exp((pref_g[At_zz[, 2]] - pref_g[At_zz[, 1]]) / DeltaD)

    } else if (identical(move_comps, "diffusion")) {
      M[At_zz] <- M[At_zz] + D / (DeltaD^2)
    } else stop("Unsupported move_comps")
  } else {
    stop("Unknown d_scaling")
  }

  row_sums <- M %*% matrix(1, nrow(M), 1)
  diag(M)  <- -row_sums[,1]
  M
}

mm <- function(formula  = ~ 0 + depth + temp,       # formula for habitat preference model
               data = out_interannual,              # environmental data with column for y = 0
               gridded_domain = grid2,              # gridded domain (class sfc)
               deployment_locs = grid0,             # grid indices for tag deployments
               release_locs = grid1,                # grid indices for tag pop-up
               time = "year",                       # time index for groups of deployed tags
               tags_per_step = c(13, 13, 37),       # how many tags per time step?
               move_comps = c("diffusion", "taxis"),# set diffusion and taxis arguments
               rtmb_exp = FALSE,                    # experimental - do matrix exponentiation in RMTB or using Matrix::expm (latter is default)
               cellsize = 25,                       # length in km of grid cell
               d_scaling = "scale2",                # scale discretization method
               apply_cov_scaling = FALSE,           # standardize covariates within years?
               DeltaT = rep(19, 3)){                # within-year tag deployment duration (defaults to 19 weeks)

  # checks
  stopifnot(
    "Time variable not found." = time %in% names(data),
    "deployment_locs must equal release_locs" =
      length(deployment_locs) == length(release_locs),
    "length(tags_per_step) must equal number of time steps" =
      length(unique(data[[time]])) == length(tags_per_step)
  )

  # create adjacency matrix
  st_rook <- function(m, ...) sf::st_relate(m, m, pattern="F***1****", ... )
  grid_A <- st_rook(gridded_domain, sparse = TRUE)
  A <- as(grid_A,"sparseMatrix") |> as("TsparseMatrix")

  # indices of non-zero elements of adjacency matrix
  At_zz = cbind( attr(A,"i"), attr(A,"j") ) + 1
  colsumA_g <- colSums(A)

  # Number of tags
  n_tags <- length(deployment_locs)

  # time steps
  times <- sort(unique(data[[time]]))
  n_t <- length(times)
  n_s <- nrow(A)

  # spline basis expansion---
  if (!is.null(formula)){
    f_mod <- mgcv::gam(formula = as.formula(paste(c("y", as.character(formula)), collapse=" ")),
                       data = data)
    X_sz <- mgcv::predict.gam(f_mod, newdata = data, type = "lpmatrix")
    if (apply_cov_scaling){
      basis_ind_split <- split(seq_len(nrow(X_sz)),
                               rep(seq_len(n_t), each = n_s))
      for (ti in seq_along(basis_ind_split)){
        rows <- basis_ind_split[[ti]]
        X_sz[rows, ] <- scale(X_sz[rows, , drop = FALSE])
      }
    }
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
    "rtmb_exp" = rtmb_exp,      # do matrix exponentiation in rtmb or Matrix
    "DeltaD" = cellsize,        # length of grid cell side
    "d_scaling" = d_scaling,    # scale discretization
    "colsumA_g" = colsumA_g,    # adjacency matrix column sums (needed for d_scaling == "scale1")
    "DeltaT" = DeltaT           # within-year tag deployment duration
  )

  if (identical(move_comps, "diffusion")) {
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
    getAll(dat, par_list, warn = FALSE)
    jnll <- 0

    # ensures the basis expansion is indexed properly (by time variable)
    if (all(c("diffusion","taxis") %in% move_comps)) {
      basis_ind_split <- split(seq_len(nrow(X_sz)),
                               rep(seq_len(n_t), each = n_s))
    } else if (identical(move_comps, "diffusion")) {
      h_s <- NULL
    }

    tag_ind <- c(1L, cumsum(n_tags))
    for (ti in seq_len(n_t)) {

      # calculates habitat preference given some parameters
      if (all(c("diffusion","taxis") %in% move_comps)) {
        h_s <- X_sz[basis_ind_split[[ti]], , drop = FALSE] %*% cbind(beta_k)
      }

      # calculates M_dot
      Mrate_gg <- make_M(A_ss,
                         At_zz,
                         rate_par = ln_D,
                         pref_g = h_s,
                         move_comps,
                         d_scaling,
                         DeltaD,
                         colsumA_g)
      tau <- DeltaT[ti] # within year time step for M = exp(DeltaT * Mrate_gg)

      # extract deployment and release grid indices
      if (ti == 1){
        dep_loc <- grid0[tag_ind[ti]:tag_ind[ti+1]]
        rel_loc <- grid1[tag_ind[ti]:tag_ind[ti+1]]
      } else {
        dep_loc <- grid0[(tag_ind[ti] + 1) : tag_ind[ti+1]]
        rel_loc <- grid1[(tag_ind[ti] + 1) : tag_ind[ti+1]]
      }

      # initial deployment locations
      v <- matrix(0, nrow = n_s, ncol = length(dep_loc))
      v[cbind(dep_loc, seq_along(dep_loc))] <- 1

      # select likelihood calcluation via RTMB::expAv or t(Matrix::expm(Mrate_gg)) %*% v
      if (rtmb_exp){
        M <- expAv(A = tau * Mrate_gg, v = v, transpose = TRUE,
                   Nmax = 10000, tol=.Machine$double.eps,
                   uniformization = TRUE)
        lik <- M[cbind(rel_loc, seq_along(rel_loc))]
        jnll <- jnll - sum(log(lik + 1e-25))
      } else {
        M   <- t(Matrix::expm(tau * Mrate_gg)) %*% v
        lik <- M[cbind(rel_loc, seq_along(rel_loc))]
        jnll <- jnll - sum(log(lik))
      }
    }
    jnll
  }

  obj <- MakeADFun(f, data = dat, par = par_list, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  AIC <- 2 * opt$objective + 2 * length(opt$par)

  list(opt = opt,
       obj = obj,
       AIC = AIC,
       model_matrix = X_sz,
       A = A,
       formula = formula,
       preference_model = f_mod,
       sd = tryCatch(sdreport(obj),
                     error = function(e) e))
}
