# Author: Kevin See
# Purpose: one script to clean PIT data and run DABOM for Asotin weir version
# Created: 2/1/22
# Last Modified: 2/2/22
# Notes:

# # install some needed packages
# install.packages(c("tidyverse",
#                    "devtools",
#                    "here",
#                    "magritter",
#                    "readxl",
#                    "writexl",
#                    "janitor",
#                    "rjags",
#                    "msm",
#                    "moments",
#                    "coda"))
#
# remotes::install_github("BiomarkABS/PITcleanr")
# remotes::install_github("BiomarkABS/DABOM")

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(DABOM)
library(tidyverse)
library(lubridate)
library(readxl)
library(writexl)
library(rjags)
library(here)

#-----------------------------------------------------------------
# create configuration file
# build configuration table (requires internet connection)
org_config = buildConfig()

# customize some nodes based on DABOM framework
configuration = org_config %>%
  filter(site_code %in% c("ACM",
                          "ASOTIC",
                          "GEORGC",
                          "ACB",
                          "AFC",
                          "CCA")) %>%
  mutate(site_code = if_else(site_code == "AFC",
                             if_else(str_detect(antenna_group,
                                                "Mainstem") |
                                       str_detect(antenna_group,
                                                  "MAINSTEM"),
                                     "AFCM",
                                     if_else(str_detect(antenna_group,
                                                        "North Fork") |
                                               str_detect(antenna_group,
                                                          "NORTH FORK"),
                                             "AFCN",
                                             "AFCS")),
                             site_code),
         node = if_else(str_detect(site_code, "AFC"),
                        paste0(site_code, "B0"),
                        node)) %>%
  # correct a couple rkm values
  mutate(rkm = if_else(site_code == 'ASOTIC',
                       '522.234.004',
                       rkm),
         rkm_total = if_else(site_code == 'ASOTIC',
                             760,
                             rkm_total),
         rkm = if_else(site_code %in% c("AFCN", "AFCS"),
                       "522.234.026",
                       rkm),
         rkm_total = if_else(site_code %in% c("AFCN", "AFCS"),
                             782,
                             rkm_total)) %>%
  mutate(latitude = if_else(site_code == "ASOTIC",
                            unique(latitude[site_code == "ACM"]),
                            latitude),
         longitude = if_else(site_code == "ASOTIC",
                             unique(longitude[site_code == "ACM"]),
                             longitude)) %>%
  filter(site_code != "ACM")

# build parent-child table
parent_child <- tribble(~ parent, ~ child,
                        "ASOTIC", "ACB",
                        "ASOTIC", "GEORGC",
                        "ACB", "CCA",
                        "ACB", "AFCM",
                        "AFCM", "AFCN",
                        "AFCM", "AFCS") %>%
  left_join(configuration %>%
              select(parent = site_code,
                     parent_rkm = rkm) %>%
              distinct()) %>%
  left_join(configuration %>%
              select(child = site_code,
                     child_rkm = rkm) %>%
              distinct())

#-----------------------------------------------------------------
# read in PTAGIS data for selected tags

# set the spawn year
yr = 2021

# assume the user has already queried this data from PTAGIS
# file path and name
ptagis_file = here("analysis/data/raw_data/PTAGIS",
                   paste0("Asotin_Sthd_", yr, ".csv"))

# recode the PTAGIS observations of double tagged fish so that the tag code matches the TagID (not TagOther)
ptagis_obs = readCTH(ptagis_file)

# any orphaned or disowned tags?
qcTagHistory(ptagis_obs, T)

# recode some of the site codes at AFC
ptagis_obs %<>%
  left_join(configuration %>%
              filter(str_detect(site_code,
                                "^AFC")) %>%
              mutate(site_code = "AFC",
                     new_site_code = str_remove(node, "B0$")) %>%
              select(event_site_code_value = site_code,
                     antenna_group_configuration_value = config_id,
                     antenna_id,
                     new_site_code),
            by = c("event_site_code_value",
                   "antenna_group_configuration_value",
                   "antenna_id")) %>%
  mutate(event_site_code_value = if_else(!is.na(new_site_code),
                                         new_site_code,
                                         event_site_code_value)) %>%
  select(-new_site_code)


#--------------------------------------------------------------
# compress and process those observations with PITcleanr
prepped_ch = PITcleanr::prepWrapper(ptagis_file = ptagis_obs,
                                    configuration = configuration,
                                    parent_child = parent_child %>%
                                      addParentChildNodes(configuration = configuration),
                                    ignore_event_vs_release = T,
                                    add_tag_detects = T,
                                    save_file = F,
                                    file_name = here('outgoing/PITcleanr',
                                                     paste0('Asotin_Sthd_', yr, '_PITcleanr.xlsx')))

# filter out detections we can't keep
# in this example, we're using PITcleanr's default auto_keep_obs
filter_obs = prepped_ch %>%
  mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                 auto_keep_obs,
                                 user_keep_obs)) %>%
  filter(user_keep_obs)

# get fish origin, based on PTAGIS marking data
fish_origin = suppressMessages(read_csv(ptagis_file)) %>%
  select(tag_code = `Tag Code`,
         origin = `Mark Rear Type Name`) %>%
  distinct() %>%
  mutate(origin = str_sub(origin, 1, 1),
         origin = recode(origin,
                         "U" = "W"))
#--------------------------------------------------------------
# file path to the default and initial model
basic_mod_file = here('analysis/model_files',
                      "Asotin_DABOM.txt")

writeDABOM(file_name = basic_mod_file,
           parent_child = parent_child,
           configuration = configuration)

# filepath for specific JAGS model code for species and year
final_mod_file = here('analysis/model_files',
                      paste0("Asotin_DABOM_", yr, ".txt"))

# writes species and year specific jags code
fixNoFishNodes(init_file = basic_mod_file,
               file_name = final_mod_file,
               filter_ch = filter_obs,
               parent_child = parent_child,
               configuration = configuration,
               fish_origin = fish_origin)


# Creates a function to spit out initial values for MCMC chains
init_fnc = setInitialValues(filter_obs,
                            parent_child,
                            configuration)

# Create all the input data for the JAGS model
jags_data = createJAGSinputs(filter_ch = filter_obs,
                             parent_child = parent_child,
                             configuration = configuration,
                             fish_origin = fish_origin)

# Tell JAGS which parameters in the model that it should save.
jags_params = setSavedParams(model_file = final_mod_file,
                             time_varying = F)


# run the model
jags = jags.model(final_mod_file,
                  data = jags_data,
                  inits = init_fnc,
                  # n.chains = 1,
                  # n.adapt = 5)
                  n.chains = 4,
                  n.adapt = 5000)

#--------------------------------------
# test the MCMC outcome and summary functions
dabom_mod = coda.samples(jags,
                         jags_params,
                         # n.iter = 10)
                         n.iter = 5000,
                         thin = 10)

#--------------------------------------------------------------
# summarize the results

# summarize detection probability estimates
detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                   filter_ch = filter_obs,
                                   cred_int_prob = 0.95)

# pull out transistion probabilities
trans_df <- extractTransProbs(dabom_mod,
                              parent_child)

# multiply them together appropriately
trans_df = trans_df %>%
  tidyr::pivot_wider(names_from = "child",
                     values_from = "value") %>%
  # multiply some probabilities together
  rowwise() %>%
  mutate(across(c(ACB_bb, CCA, AFCM),
                ~ . * ACB)) %>%
  mutate(across(c(AFCN, AFCS, AFCM_bb),
                ~ . * AFCM)) %>%
  ungroup() %>%
  mutate(iter = 1:n()) %>%
  tidyr::pivot_longer(cols = -c(CHAIN, ITER,
                                iter, origin),
                      names_to = "param",
                      values_to = "value") %>%
  select(chain = CHAIN,
         iter,
         origin,
         param,
         value) %>%
  mutate(origin = recode(origin,
                         "2" = "H",
                         "1" = "W"))


# total escapement past Asotin weir, by origin
# set equal to 1 if we want to just look at proportions
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

escape_summ %>%
  filter(origin == "W")
