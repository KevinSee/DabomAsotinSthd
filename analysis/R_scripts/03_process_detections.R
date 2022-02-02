# Author: Kevin See
# Purpose: clean PTAGIS data with PITcleanr
# Created: 2/1/22
# Last Modified: 2/1/22
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(lubridate)
library(janitor)
library(readxl)
library(magrittr)
library(here)

#-----------------------------------------------------------------
# load configuration and site_df data
load(here('analysis/data/derived_data/site_config.rda'))

# which spawn year are we dealing with?
yr = 2021

# for(yr in 2011:2020) {

# # load and file biological data
# bio_df = read_rds(here('analysis/data/derived_data/Bio_Data_2011_2021.rds')) %>%
#   filter(year == yr)
#
# # any double-tagged fish?
# dbl_tag = bio_df %>%
#   filter(!is.na(tag_other))


# #-----------------------------------------------------------------
# # start date is June 1 of previous year
# start_date = paste0(yr - 1, '0601')
# # when is the last possible observation date?
# max_obs_date = paste0(yr, "0531")

# get raw observations from PTAGIS
# These come from running a saved query on the list of tags to be used
ptagis_file = here("analysis/data/raw_data/PTAGIS",
                   paste0("Asotin_Sthd_", yr, ".csv"))

# recode the PTAGIS observations of double tagged fish so that the tag code matches the TagID (not TagOther)
ptagis_obs = readCTH(ptagis_file)

# any orphaned or disowned tags?
qcTagHistory(ptagis_obs, T)

# compress and process those observations with PITcleanr
prepped_ch = PITcleanr::prepWrapper(ptagis_file = ptagis_obs,
                                    configuration = configuration,
                                    parent_child = parent_child %>%
                                      addParentChildNodes(configuration = configuration),
                                    # min_obs_date = start_date,
                                    # max_obs_date = max_obs_date,
                                    ignore_event_vs_release = T,
                                    add_tag_detects = T,
                                    save_file = T,
                                    file_name = here('outgoing/PITcleanr', paste0('Asotin_Sthd_', yr, '_PITcleanr.xlsx')))


# save some stuff
save(parent_child, configuration,
     # start_date, bio_df,
     prepped_ch,
     file = here('analysis/data/derived_data/PITcleanr',
                 paste0('Asotin_Sthd_', yr, '.rda')))


# rm(start_date, bio_df, prepped_ch,
#    ptagis_file,
#    ptagis_obs,
#    dbl_tag)
# }

#-------------------------------------------
# NEXT STEPS
#-------------------------------------------
# open that Excel file, and filter on the column user_keep_obs, looking for blanks. Fill in each row with TRUE or FALSE, depending on whether that observation should be kept or not. The column auto_keep_obs provides a suggestion, but the biologist's best expert judgment should be used based on detection dates, detection locations before and after, etc.

load(here('analysis/data/derived_data/PITcleanr',
          paste0('UC_Steelhead_', yr, '.rda')))

wdfw_df = read_excel(here('analysis/data/derived_data/WDFW',
                          paste0('Asotin_Sthd_', yr, '.xlsx')))

filter_obs = wdfw_df %>%
  mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                 auto_keep_obs,
                                 user_keep_obs)) %>%
  filter(user_keep_obs)

# construct all valid paths
all_paths = buildPaths(addParentChildNodes(parent_child,
                                           configuration))

tag_path = summarizeTagData(filter_obs,
                            bio_df %>%
                              group_by(tag_code) %>%
                              slice(1) %>%
                              ungroup()) %>%
  select(tag_code, spawn_node) %>%
  distinct() %>%
  left_join(all_paths,
            by = c('spawn_node' = 'end_loc')) %>%
  rename(tag_path = path)

error_tags = filter_obs %>%
  left_join(tag_path) %>%
  rowwise() %>%
  mutate(node_in_path = str_detect(tag_path, node)) %>%
  ungroup() %>%
  filter(!node_in_path) %>%
  select(tag_code) %>%
  distinct()

nrow(error_tags)
error_tags %>%
  left_join(wdfw_df) %>%
  as.data.frame()

prepped_ch = wdfw_df

save(parent_child, configuration, start_date, bio_df, prepped_ch,
     file = here('analysis/data/derived_data/PITcleanr',
                 paste0('UC_Steelhead_', yr, '.rda')))


wdfw_obs <- read_csv(here("analysis/data/raw_data",
                          "Asotin2021FemalesDABOM.csv")) %>%
  rename(tag_code = `Tag Code`) %>%
  pivot_longer(cols = -tag_code,
               names_to = "node",
               values_to = "seen") %>%
  filter(seen == 1) %>%
  mutate(node = recode(node,
                       "ACBDA" = "ACBB0",
                       "ACBUA" = "ACBA0",
                       "CCADA" = "CCAB0",
                       "CCAMA" = "CCAA0",
                       "CCAUA" = "CCAA0",
                       "AFCDA" = "AFCB0",
                       "AFCSA" = "AFCA0",
                       "AFCNA" = "AFCA0")) %>%
  distinct()

prepped_ch %>%
  # filter(!is.na(user_keep_obs)) %>%
  filter(auto_keep_obs) %>%
  select(tag_code, node) %>%
  distinct() %>%
  mutate(pitclr = T) %>%
  full_join(wdfw_obs %>%
              mutate(wdfw = T)) %>%
  filter(node != "ASOTIC") %>%
  filter(is.na(pitclr) | is.na(wdfw))

#-----------------------------------------------------------------
# tag summaries
#-----------------------------------------------------------------
# use auto_keep_obs for the moment
tag_summ = summarizeTagData(prepped_ch %>%
                              mutate(user_keep_obs = if_else(is.na(user_keep_obs),
                                                             auto_keep_obs,
                                                             user_keep_obs)),
                            bio_df %>%
                              group_by(tag_code) %>%
                              slice(1) %>%
                              ungroup())

# any duplicated tags?
sum(duplicated(tag_summ$tag_code))
tag_summ %>%
  filter(tag_code %in% tag_code[duplicated(tag_code)]) %>%
  as.data.frame()

# where are tags assigned?
janitor::tabyl(tag_summ, spawn_node) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals()

# preliminary estimate of node efficiency
node_eff = prepped_ch %>%
  mutate(user_keep_obs = auto_keep_obs) %>%
  filter(user_keep_obs) %>%
  estNodeEff(node_order = buildNodeOrder(addParentChildNodes(parent_child, configuration)))

node_eff %>%
  filter(tags_at_node > 0,
         eff_est < 1)

#-----------------------------------------------------------------
# examine some of the output
#-----------------------------------------------------------------
# which tags have "strange" capture histories?
prepped_ch %>%
  summarise(n_tags = n_distinct(tag_code),
            n_weird = n_distinct(tag_code[direction == "unknown"]),
            n_fix = n_distinct(tag_code[is.na(user_keep_obs)]),
            prop_weird = n_weird / n_tags,
            prop_fix = n_fix / n_tags)

# look at which branch each tag was assigned to for spawning
brnch_df = buildNodeOrder(addParentChildNodes(parent_child, configuration)) %>%
  separate(col = path,
           into = paste("step", 1:max(.$node_order), sep = "_"),
           remove = F) %>%
  mutate(branch_nm = if_else(node == "PRA",
                             "Start",
                             if_else(grepl('LWE', path) | node %in% c("CLK"),
                                     "Wenatchee",
                                     if_else(grepl("ENL", path),
                                             "Entiat",
                                             if_else(grepl("LMR", path),
                                                     "Methow",
                                                     if_else(grepl("OKL", path) | node %in% c("FST"),
                                                             "Okanogan",
                                                             if_else(step_2 != "RIA" & !is.na(step_2),
                                                                     "Downstream",
                                                                     "Mainstem"))))))) %>%
  select(-starts_with("step"))

tag_summ %<>%
  left_join(brnch_df %>%
              select(spawn_node = node,
                     branch_nm))

# how many tags in each branch?
tag_summ %>%
  janitor::tabyl(branch_nm) %>%
  janitor::adorn_pct_formatting() %>%
  arrange(desc(n))

# age comp in each branch, by sex
tag_summ %>%
  filter(!is.na(final_age)) %>%
  ggplot(aes(x = branch_nm,
             fill = as.ordered(final_age))) +
  geom_bar(position = position_fill()) +
  facet_wrap(~ sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1)) +
  labs(x = "Branch",
       y = "Percent of Tags",
       fill = "Age")


# look at run timing between branches
tag_summ %>%
  ggplot(aes(x = trap_date,
             color = branch_nm,
             fill = branch_nm)) +
  geom_density(alpha = 0.2) +
  theme_bw() +
  scale_color_brewer(palette = 'Set1',
                     name = "Branch") +
  scale_fill_brewer(palette = 'Set1',
                    name = "Branch") +
  labs(x = "Trap Date at Priest Rapids")
