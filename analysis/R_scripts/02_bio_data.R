# Author: Kevin See
# Purpose: create tag lists to feed to PTAGIS query
# Created: 4/1/2020
# Last Modified: 1/18/2022
# Notes:

#-----------------------------------------------------------------
# load needed libraries
# library(PITcleanr)
library(tidyverse)
library(readxl)
library(lubridate)
library(janitor)
library(magrittr)
library(writexl)
library(here)

#-----------------------------------------------------------------
# read in list of tags from data sent by Dan Rawding
max_yr = 2021
tag_list <- read_csv(here("analysis/data/raw_data",
              "Asotin2021FemalesDABOM.csv")) %>%
  select(tag_code = `Tag Code`)

# save to be uploaded to PTAGIS
write_delim(tag_list,
            file = here('analysis/data/raw_data/tag_lists',
                        paste0('Asotin_Sthd_Tags_', max_yr, '.txt')),
            delim = '\n',
            col_names = F)



# # save biological data for later
# write_rds(bio_df,
#           file = here('analysis/data/derived_data',
#                       paste0('Bio_Data_', min_yr, '_', max_yr, '.rds')))


