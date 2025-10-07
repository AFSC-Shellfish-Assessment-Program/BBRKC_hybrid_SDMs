# Projected movement probs-----
project_distribution <- function(year = 2005,
                                 month = 10,
                                 temp_data = pred_grid,
                                 ests = pref_mod$ests,
                                 f_mod = pref_mod$mod,
                                 adj = A_big,
                                 pred_df = pred,
                                 dep_locs = grid0,
                                 move_comps = c("diffusion", "taxis"),
                                 cellsize = 25){

  # filter environmental data
  env_cov <- temp_data %>%
    filter(year == !!year,
           month == !!month)

  # prediction grid
  pred_grid2 <-
    env_cov %>%
    st_as_sfc()

  # convert to data.frame for prediction
  out <- env_cov %>%
    st_set_geometry(NULL) %>%
    mutate(y = 0,
           year = factor(year))

  ## Extract estimated parameters from preference function
  beta_k <- ests[names(ests) != "ln_D"]
  ln_D <- ests[names(ests) == "ln_D"]

  ## model matrix
  X_gk_n <- predict(f_mod, newdata = out, type = "lpmatrix")
  At_zz = cbind( attr(adj,"i"), attr(adj,"j") ) + 1

  h_s <- X_gk_n %*% beta_k
  Mrate <-  M_dot(A_ss = adj,
                  At_zz,
                  rate_par = ln_D,
                  pref_g = h_s,
                  move_comps,
                  d_scaling = "none",
                  delta_d = cellsize)
  M <- expm::expm(Mrate, do.sparseMsg = F)

  ## do the projection into october
  pred_sel <- pred_df %>%
    filter(year == !!year)

  if (year != 2020){
    proj_out <-
      pred_grid2 %>%
      st_as_sf() %>%
      mutate(  projected  = as.vector(rbind(pred_sel$est) %*% M),
               month,
               year)
  } else {
    proj_out <- NA
  }

  return(list(proj_out = proj_out))
}
