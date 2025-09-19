# Projected movement probs-----
project_distribution <- function(year = 2005,
                                 month = 10,
                                 temp_data = pred_grid,
                                 ests = pref_mod$ests,
                                 f_mod = pref_mod$mod,
                                 adj = A_big,
                                 pred_df = pred,
                                 dep_locs = grid0,
                                 crop = ebs,
                                 diffusion_only = F,
                                 move_comps = c("diffusion", "taxis"),
                                 DeltaT = 19){

  # filter environmental data
  env_cov <- temp_data %>%
    filter(year == !!year,
           month == !!month)

  # prediction grid
  pred_grid2 <-
    env_cov %>%
    st_as_sfc()

  # convert to data.frame for prediction----
  out <- env_cov %>%
    st_set_geometry(NULL) %>%
    mutate(y = 0,
           year = factor(year))

  ## Extract estimated parameters from preference function
  beta_k <- ests[names(ests) != "ln_D"]
  diff_par <- ests[names(ests) == "ln_D"]

  ## model matrix------
  X_gk_n <- predict(f_mod, newdata = out, type = "lpmatrix")
  At_zz = cbind( attr(adj,"i"), attr(adj,"j") ) + 1

  pref_g <- X_gk_n %*% beta_k
  Mrate <- make_M(  A_ss = adj,
                    At_zz,
                    rate_par = diff_par,
                    pref_g,
                    move_comps,
                    d_scaling = "scale2",
                    DeltaD = 25,
                    colsumA_g = colSums(adj))
  M <- expm::expm(DeltaT * Mrate, do.sparseMsg = F)

  ## do the projection into october-----
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
