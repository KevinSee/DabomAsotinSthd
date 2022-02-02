# Author: Kevin See
# Purpose: summarize DABOM results
# Created: 2/1/22
# Last Modified: 2/1/22
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(DABOM)
library(PITcleanr)
library(tidyverse)
library(magrittr)
library(readxl)
library(STADEM)
library(writexl)
library(msm)
library(moments)
library(coda)
library(here)

#-----------------------------------------------------------------
# set year
yr = 2021

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data',
          'site_config.rda'))


# for(yr in 2011:2020) {
  cat(paste("Working on", yr, "\n\n"))

  # load compressed detections and biological data
  load(here('analysis/data/derived_data/PITcleanr',
            paste0('Asotin_Sthd_', yr, '.rda')))

  # load JAGS MCMC results
  load(here("analysis/data/derived_data/model_fits",
            paste0('Asotin_DABOM_Sthd_', yr,'.rda')))


  # estimate final spawning location
  tag_summ = summarizeTagData(filter_obs,
                              bio_data = NULL)

  # look at which branch each tag was assigned to for spawning
  xtabs(~ spawn_node, tag_summ)

  # summarize detection probabilities
  detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                     filter_ch = filter_obs)

  # which sites had detection probabilities fixed at 0% or 100%
  detect_summ %>%
    filter(sd == 0)

  # look at all the other sites
  detect_summ %>%
    filter(sd > 0) %>%
    # arrange(desc(n_tags))
    arrange(desc(sd))

  # compile all movement probabilities, and multiply them appropriately
  trans_df = compileTransProbs_ASOTIC(dabom_mod,
                                      parent_child) %>%
    mutate(origin = recode(origin,
                           "2" = "H",
                           "1" = "W"))

  # summarize transition probabilities
  trans_summ = trans_df %>%
    group_by(origin, param) %>%
    summarise(mean = mean(value),
              median = median(value),
              mode = estMode(value),
              sd = sd(value),
              skew = moments::skewness(value),
              kurtosis = moments::kurtosis(value),
              lowerCI = coda::HPDinterval(coda::as.mcmc(value))[,1],
              upperCI = coda::HPDinterval(coda::as.mcmc(value))[,2],
              .groups = "drop") %>%
    mutate(across(c(mean, median, mode, sd, matches('CI$')),
                  ~ if_else(. < 0, 0, .)))

  trans_summ %>%
    filter(origin == "W")

  #-----------------------------------------------------------------
  # total escapement past Asotin weir, by origin
  org_escape <- tibble(year = yr,
                       origin = c("W", "H"),
                       tot_escp = c(1, 0),
                       tot_escp_se = c(0, 0))

  # translate movement estimates to escapement
  escape_post = trans_df %>%
    left_join(org_escape %>%
                group_by(origin) %>%
                summarise(tot_esc_samp = map2(tot_escp,
                                              tot_escp_se,
                                              .f = function(x, y) {
                                                tibble(tot_escp = rnorm(max(trans_df$iter),
                                                                        mean = x,
                                                                        sd = y)) %>%
                                                  mutate(iter = 1:n())
                                              })) %>%
                unnest(cols = tot_esc_samp)) %>%
    mutate(escp = value * tot_escp)

  escape_summ = escape_post %>%
    group_by(origin, location = param) %>%
    summarise(mean = mean(escp),
              median = median(escp),
              mode = estMode(escp),
              sd = sd(escp),
              skew = moments::skewness(escp),
              kurtosis = moments::kurtosis(escp),
              lowerCI = coda::HPDinterval(coda::as.mcmc(escp))[,1],
              upperCI = coda::HPDinterval(coda::as.mcmc(escp))[,2],
              .groups = 'drop') %>%
    mutate(across(c(mean, median, mode, sd, matches('CI$')),
                  ~ if_else(. < 0, 0, .))) %>%
    mutate(across(c(mean, median, mode, sd, skew, kurtosis, matches('CI$')),
                  round,
                  digits = 2)) %>%
    arrange(desc(origin), location) %>%
    tibble::add_column(species = "Steelhead",
                       spawn_year = yr,
                       .before = 0)

  #-----------------------------------------------------------------
  # write results to an Excel file
  save_list = c(list('All Escapement' = escape_summ %>%
                       select(-skew, -kurtosis) %>%
                       mutate(across(mean:mode,
                                     janitor::round_half_up)) %>%
                       mutate(across(sd:upperCI,
                                     round,
                                     digits = 1)) %>%
                       rename(estimate = mean,
                              se = sd) %>%
                       select(-median, -mode),
                     'Detection' = detect_summ %>%
                       mutate(across(-c(node, n_tags),
                                     round,
                                     digits = 3)) %>%
                       rename(estimate = mean,
                              se = sd) %>%
                       select(-median, -mode),
                     'Tag Summary' = tag_summ))

  writexl::write_xlsx(x = save_list,
                      path = here('outgoing/estimates',
                                  paste0('Asotin_Steelhead_', yr, '_', format(Sys.Date(), '%Y%m%d'), '.xlsx')))

# }


