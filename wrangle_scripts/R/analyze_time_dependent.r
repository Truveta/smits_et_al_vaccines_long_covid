##########
# time-dependent parts

# prepare data for analysis
make_timedep_data <- function(df, outcome, keep_zeroes = TRUE, only_after = FALSE) {
  df_inter <- 
    df |>
    dplyr::mutate(
      event = dplyr::if_else(is.na({{ outcome }}), 0, 1),
      last_encounter =
        dplyr::if_else(
          last_encounter < time_zero_datetime | is.na(last_encounter),
          time_zero_datetime,
          last_encounter
        ),
      outcome_time = dplyr::coalesce({{ outcome }}, last_encounter),
      outcome_time_days = 
        as.numeric(difftime(outcome_time, time_zero_datetime, units = 'days')),
      outcome_time_weeks =
        as.numeric(difftime(outcome_time, time_zero_datetime, units = 'weeks')),
      event = 
        dplyr::if_else(outcome_time_days > 365, 0, event),
      outcome_time_days = 
        dplyr::if_else(outcome_time_days >= 365, 365, outcome_time_days),
      outcome_time_weeks = 
        dplyr::if_else(outcome_time_weeks >= 52, 52, outcome_time_weeks),
      t_vaccine = 
          as.numeric(difftime(age_vaccine_completed, time_zero_datetime, units = 'weeks')),
      t_vaccine = dplyr::if_else(t_vaccine < 0, 0, t_vaccine),
      t_vaccine = 
          dplyr::if_else(t_vaccine >= outcome_time_weeks, NA_real_, t_vaccine),
      t_boost = 
        as.numeric(
          difftime(age_vaccine_boost, time_zero_datetime, units = 'weeks')
        ),
      t_boost = dplyr::if_else(t_boost < 0, 0, t_boost),
      t_boost = dplyr::if_else(t_boost >= outcome_time_weeks, NA_real_, t_boost),
      t_start = 0,
      id = 1:n()
    )  

  # if only dealing with vaccination events which occur after covid
  # this should *implicitly* also do this with booster doses
  if(only_after) {
    df_inter <- 
      df_inter |>
      dplyr::filter(t_vaccine > 0 | is.na(t_vaccine))
  }

  # find the individuals with zero outcome time. 
  # we don't necessarily want to include them
  # but if we do, we cannot violate rule that S(t = 0) = 1
  # so add super tiny amount of time to the outcome
  df_zeroes <-
    df_inter |> 
    dplyr::filter(outcome_time_weeks == 0) |> 
    dplyr::mutate(
      tstart = 0,
      tstop = outcome_time_weeks + 5e-7,
      long_covid = 0,  # impossible to be diagnosed with long covid if they had no encounters!!!
      vaccinated = dplyr::if_else(is.na(t_vaccine), 0, 1),
      boosted = dplyr::if_else(is.na(t_boost), 0, 1)      
    )

  # given the above is to find the zeroes, this removes them
  # can add back later because no actual time dependence!
  df_inter <- 
    df_inter |> 
    dplyr::filter(outcome_time_weeks > 0)

  # initial tmerge to create "bedrock" with id and 
  # tstart, tstop, and outcome (long covid)
  # tmerge works through multiple calls to tmerge, aggregating the final table in place
  init <- 
    survival::tmerge(
      df_inter, 
      df_inter, 
      id = id, 
      long_covid = event(outcome_time_weeks, event)
    )

  # add in time dependent covariate
  temp <-
    df_inter |> 
    dplyr::select(id, t_vaccine)

  df_inter_dat <- 
    survival::tmerge(
      init, 
      temp, 
      id = id, 
      vaccinated = tdc(t_vaccine)
    )

  # want to add a second set of time dependent covariates around boosters
  temp2 <- 
    df_inter |> 
    dplyr::select(id, t_boost)

  df_inter_dat <- 
    survival::tmerge(
      df_inter_dat, 
      temp2, 
      id = id, 
      boosted = tdc(t_boost)
    )
  

  if(keep_zeroes) {
    out <- bind_rows(df_inter_dat, df_zeroes)
  } else {
    out <- df_inter_dat
  }
    
  out
}

# select key columns for analysis, and do some basic transforms
prep_timedep_data <- function(df) {
  df_prep <-
    df |> 
    dplyr::mutate(
      flu_vaccinated = dplyr::if_else(flu_vaccination == 0, 0, 1),
      inpatient_count_sqrt = sqrt(inpatient_count),
      outpatient_count_sqrt = sqrt(outpatient_count),
      blood_panel_count_sqrt = sqrt(blood_panel_count),
      year_month = as.factor(year_month)
    ) |>
    dplyr::select(
      tstart, 
      tstop, 
      long_covid, 
      vaccinated, 
      boosted,
      age_covid, 
      sex, 
      race, 
      ethnicity, 
      year_month:depression, 
      flu_vaccinated,
      inpatient_count_sqrt,
      outpatient_count_sqrt,
      blood_panel_count_sqrt
    )

  df_prep
}


# km for time dependent data
analyze_km_timedep <- function(data) {
  fit_time_dep <- 
    survfit(
      Surv(time = tstart, time2 = tstop, event = long_covid) ~ vaccinated + boosted, 
      data = data
    )

  g <- 
    ggsurvplot(
      fit_time_dep,
      data = data,
      conf.int = TRUE,
      conf.int.style = "step",  # customize style of confidence intervals
      xlab = 'Time (weeks)',
      ylab = 'Survival probability',      
      subtitle = 'vaccinated and boosted modeled as time dependent covariates',
      ylim = c(0.6, 1)
    )

  pp <- 
    g$plot + 
    #ggplot2::coord_cartesian(ylim = c(0.5, 1)) +
    theme_truveta() +
    scale_colour_manual(
      labels = c('unvaccinated', 'primary vaccination', 'primary vaccination and booster'),
      values = c("#E69F00", "#56B4E9", "#009E73") 
    )
  

  out <- list()
  out$km_fit <- fit_time_dep
  out$km_graph <- pp
  out  
}

# fit cox model and provide meaningful outputs for time dependent analysis
analyze_cox_timedep <- function(data) {
  fit_cox_dep <- 
    coxph(
      Surv(time = tstart, time2 = tstop, event = long_covid) ~ 
        . + 
        ns(age_covid, df = 5) - 
        age_covid,
      data = data
    )
    
  fit_cox_dep_tidy <- cox_table(fit_cox_dep)

  fit_cox_dep_hr <-
    fit_cox_dep_tidy |>
    dplyr::filter(term %in% c('vaccinated', 'boosted')) |>
    dplyr::mutate(term = forcats::fct_rev(term)) |>
    ggplot2::ggplot(ggplot2::aes(x = estimate_exp, y = term)) +
    ggplot2::geom_vline(xintercept = 1) +
    ggplot2::geom_pointrange(
      mapping = ggplot2::aes(xmin = conf_low_exp, xmax = conf_high_exp),
      size = 1.15
    ) +
    theme_truveta() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
    labs(
      x = 'Hazard Ratio', 
      y = 'Covariate',
      subtitle = 'vaccinated and boosted modeled as time dependent covariates',
    )
    
  out <- list()
  out$cox_fit <- fit_cox_dep
  out$cox_table <- fit_cox_dep_tidy
  out$cox_hr <- fit_cox_dep_hr
  out    
}


wrap_emmeans_dep <- function(cox_fit, non_nuisance) {
  rg <- ref_grid(cox_fit, non.nuisance = non_nuisance)
  emm <- emmeans::emmeans(object = rg, specs = non_nuisance, type = 'response')

  pai <- 
    pairs(emm, reverse = TRUE) |>
    broom::tidy(conf.int = TRUE)
    
  pai |>
  dplyr::filter(!(contrast %in% c('0 1 / 0 0', '0 1 / 1 0', '1 1 / 0 1'))) |>
  dplyr::filter(
    !(contrast %in% 
      c(
        'vaccinated0 boosted1 / vaccinated0 boosted0', 
        'vaccinated0 boosted1 / vaccinated1 boosted0', 
        'vaccinated1 boosted1 / vaccinated0 boosted1'
      )
     )
    ) |>
  dplyr::mutate(
    comparison = 
      c(
        'vaccinated vs unvaccinated', 
        'vaccinated and boosted vs unvaccinated', 
        'vaccinated and boosted vs vaccinated'
      ),
    comparison = forcats::fct_inorder(comparison),
    comparison = forcats::fct_rev(comparison)
  )
}


# wrap everything into a single function to decrease conginitive load in "main"
time_dependent_analysis <- function(data, outcome, keep_zeroes = TRUE, only_after = FALSE) {
  time <- 
    make_timedep_data(
      data, 
      {{ outcome }}, 
      keep_zeroes = keep_zeroes, 
      only_after = only_after
    )
  
  time_prep <- prep_timedep_data(time)
  
  time_km <- analyze_km_timedep(time_prep)
  
  time_cox <- analyze_cox_timedep(time_prep)
  
  time_table <- 
    wrap_emmeans_dep(
      time_cox$cox_fit,
      non_nuisance = c('vaccinated', 'boosted')
    )
  
  out <- list()
  out$time <- time
  out$time_prep <- time_prep
  out$time_km <- time_km
  out$time_cox <- time_cox
  out$time_table <- time_table
  
  out
}
