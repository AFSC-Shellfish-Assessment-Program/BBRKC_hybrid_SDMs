effectDT <- function(model = m6,
                     terms = c("depth","temp"),
                     df = out_interannual){
  pref_mod <- model$preference_model
  cov_names <- names(pref_mod$model)
  term_var <- cov_names[cov_names %in% terms & !cov_names %in% "y"]

  # drop diffusion parameter from estimates and covariance matrix
  ests <- model$sd$par.fixed
  ests <- ests[names(ests) != "ln_D"]
  cov_mat <- model$sd$cov.fixed
  cov_mat <- cov_mat[-which(rownames(cov_mat) == "ln_D"), -which(colnames(cov_mat) == "ln_D")]
  X_gk <- model$model_matrix

  out <- NULL
  for (i in term_var){
    hold_var <- cov_names[!cov_names %in% c(i, "y")]

    term_min <-
      df %>%
      dplyr::select(all_of(i)) %>%
      dplyr::summarise_all(min)  %>%
      pull(.)

    term_max <-
      df %>%
      dplyr::select(all_of(i)) %>%
      dplyr::summarise_all(max)  %>%
      pull(.)

    term_seq <- seq(term_min, term_max, length.out = 100)

    hold_var_mean <-
      df %>%
      dplyr::select(all_of(hold_var)) %>%
      dplyr::summarise_all(mean)

    pred_df <- tibble(term_seq,
                      hold_var_mean)
    names(pred_df)[1] <- i
    Xpred <- predict(pref_mod, newdata = pred_df, type = "lpmatrix")

    marginal_effects <-
      tibble(pref = as.numeric(scale(Xpred %*% ests,
                                     center = T, scale = F)),
             pref_se = as.numeric(sqrt(diag(Xpred %*% cov_mat %*% t(Xpred)))),
             pref_lwr = pref - 1.96 * pref_se,
             pref_upr = pref + 1.96 * pref_se,
             term = i) %>%
      bind_cols(., pred_df %>%
                  dplyr::select(x = 1))

    assign("out", rbind(out, marginal_effects))
  }

  return(out)
}
