########
# time-independent parts

# prepare a dataframe for the "simple" analysis
make_simple_data <- function(df, outcome, keep_zeroes = TRUE) {
  df_simple <- 
    df |>
    dplyr::filter(
      is.na(age_vaccine_completed_years) |
        age_covid >= age_vaccine_completed_years,
      is.na(age_vaccine_boost_years) |
        age_covid >= age_vaccine_boost_years
    ) |>
    dplyr::mutate(
      event = dplyr::if_else(is.na({{ outcome }}), 0, 1),
      vaccinated_at_covid = dplyr::if_else(is.na(age_vaccine_completed_years), 0, 1),
      boosted_at_covid = dplyr::if_else(is.na(age_vaccine_boost_years), 0, 1),
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
      event = dplyr::if_else(outcome_time_days > 365, 0, event),
      outcome_time_days = 
        dplyr::if_else(outcome_time_days >= 365, 365, outcome_time_days),
      outcome_time_weeks = 
        dplyr::if_else(outcome_time_weeks >= 52, 52, outcome_time_weeks)
    )
  
  if(keep_zeroes) {
    df_simple <- 
      df_simple |>
      dplyr::mutate(
        outcome_time_days = 
          dplyr::if_else(outcome_time_days == 0, 5e-7, outcome_time_days),
        outcome_time_weeks = 
          dplyr::if_else(outcome_time_weeks == 0, 5e-7, outcome_time_weeks)
      )
  } else {
    df_simple <-
      df_simple |>
      dplyr::filter(outcome_time_weeks == 0)
  }

  df_simple
}


# select key columns
prep_simple_data <- function(df) {  
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
      vaccinated_at_covid, 
      boosted_at_covid,
      outcome_time_weeks, 
      event, 
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


# estimate kmcurve and ouput useful values
analyze_km_simple <- function(data) {
  fit_km_simple <- 
    survfit(
      Surv(time = outcome_time_weeks, event = event) ~ vaccinated_at_covid + boosted_at_covid,
      data = data
    )

  g <- 
    ggsurvplot(
      fit_km_simple,
      data = data,
      conf.int = TRUE,
      conf.int.style = "step",  # customize style of confidence intervals
      xlab = 'Time (weeks)',
      ylab = 'Survival probability',
      subtitle = 'vaccination/boosting after COVID infection not considered',
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
  out$km_fit <- fit_km_simple 
  out$km_graph <- pp
  out  
}


# fit cox model and provide useful outputs
analyze_cox_simple <- function(data) {
  fit_cox_simple <- 
    coxph(
      Surv(time = outcome_time_weeks, event = event) ~ 
        . + 
        ns(age_covid, df = 5) - 
        age_covid,
      data = data
    )

  fit_cox_simple_tidy <- cox_table(fit_cox_simple)


  fit_cox_simple_hr <-
    fit_cox_simple_tidy |>
    dplyr::filter(term %in% c('vaccinated_at_covid', 'boosted_at_covid')) |>
    dplyr::mutate(term = forcats::fct_rev(term)) |>
    ggplot2::ggplot(ggplot2::aes(x = estimate_exp, y = term)) +
    ggplot2::geom_vline(xintercept = 1) +
    ggplot2::geom_pointrange(
      mapping = ggplot2::aes(xmin = conf_low_exp, xmax = conf_high_exp)
    ) +  
    theme_truveta() +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 10)) +
    labs(
      x = 'Hazard Ratio', 
      y = 'Covariate',
      title = str_wrap('Esimated hazard ratios for time till Long COVID', 30),
      subtitle = str_wrap('considers only patients either vaccinated before COVID or never vaccinated', 60)
    )

  out <- list()
  out$cox_fit <- fit_cox_simple
  out$cox_table <- fit_cox_simple_tidy
  out$cox_hr <- fit_cox_simple_hr
  out    
}




wrap_emmeans <- function(cox_fit, non_nuisance) {
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
        'vaccinated_at_covid0 boosted_at_covid1 / vaccinated_at_covid0 boosted_at_covid0', 
        'vaccinated_at_covid0 boosted_at_covid1 / vaccinated_at_covid1 boosted_at_covid0', 
        'vaccinated_at_covid1 boosted_at_covid1 / vaccinated_at_covid0 boosted_at_covid1'
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
simple_analysis <- function(data, outcome, keep_zeroes = TRUE) {
  simple <- make_simple_data(data, {{ outcome }}, keep_zeroes = keep_zeroes)
  simple_prep <- prep_simple_data(simple)
  
  simple_km <- analyze_km_simple(simple_prep)
  
  simple_cox <- analyze_cox_simple(simple_prep)
  
  sp_table <- 
    wrap_emmeans(
      simple_cox$cox_fit, 
      non_nuisance = c('vaccinated_at_covid', 'boosted_at_covid')
    )
  
  out <- list()
  out$simple <- simple
  out$simple_prep <- simple_prep
  out$simple_km <- simple_km
  out$simple_cox <- simple_cox
  out$sp_table <- sp_table
  
  out
}
