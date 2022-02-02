# Author: Kevin See
# Purpose: Develop configuration file for DABOM
# Created: 2/1/2022
# Last Modified: 2/1/2022
# Notes:
#
# # install some needed packages
# install.packages(c("tidyverse",
#                    "devtools",
#                    "here",
#                    "sf",
#                    "magritter",
#                    "readxl",
#                    "writexl",
#                    "janitor",
#                    "rjags",
#                    "msm",
#                    "moments",
#                    "coda"))
#
# remotes::install_github("BiomarkABS/STADEM")
# remotes::install_github("BiomarkABS/PITcleanr")
# remotes::install_github("BiomarkABS/DABOM")

#-----------------------------------------------------------------
# load needed libraries
library(PITcleanr)
library(tidyverse)
library(magrittr)
library(sf)
library(here)

#-----------------------------------------------------------------
# set starting point
root_site = "ASOTIC"

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
  mutate(node = if_else(site_code == "AFC",
                        if_else(str_detect(antenna_group,
                                           "Mainstem") |
                                  str_detect(antenna_group,
                                             "MAINSTEM"),
                                "AFCB0",
                                "AFCA0"),
                        node)) %>%
  # correct a couple rkm values
  mutate(rkm = if_else(site_code == 'ASOTIC',
                       '522.234.004',
                       rkm),
         rkm_total = if_else(site_code == 'ASOTIC',
                             760,
                             rkm_total)) %>%
  mutate(latitude = if_else(site_code == "ASOTIC",
                            unique(latitude[site_code == "ACM"]),
                            latitude),
         longitude = if_else(site_code == "ASOTIC",
                            unique(longitude[site_code == "ACM"]),
                            longitude)) %>%
    filter(site_code != "ACM")

# Node network for DABOM

# get spatial object of sites used in model
sites_sf = configuration %>%
  group_by(site_code) %>%
  filter(config_id == max(config_id)) %>%
  ungroup() %>%
  select(site_code,
         site_name,
         site_type = site_type_name,
         type = site_type,
         rkm,
         site_description = site_description,
         latitude, longitude) %>%
  distinct() %>%
  filter(!is.na(latitude)) %>%
  st_as_sf(coords = c("longitude",
                      "latitude"),
           crs = 4326) %>%
  st_transform(crs = 5070)

#-----------------------------------------------------------------
# download the NHDPlus v2 flowlines
# do you want flowlines downstream of root site? Set to TRUE if you have downstream sites
dwn_flw = F
nhd_list = queryFlowlines(sites_sf = sites_sf,
                          root_site_code = root_site,
                          min_strm_order = 2,
                          dwnstrm_sites = dwn_flw,
                          dwn_min_stream_order_diff = 4)

# compile the upstream and downstream flowlines
flowlines = nhd_list$flowlines
if(dwn_flw) {
  flowlines %<>%
    rbind(nhd_list$dwn_flowlines)
}


#-----------------------------------------------------------------
# plot the flowlines and the sites
ggplot() +
  geom_sf(data = flowlines,
          aes(color = as.factor(StreamOrde),
              size = StreamOrde)) +
  scale_color_viridis_d(direction = -1,
                        option = "D",
                        end = 0.8) +
  scale_size_continuous(range = c(0.2, 1.2),
                        guide = 'none') +
  geom_sf(data = nhd_list$basin,
          fill = NA,
          lwd = 2) +
  geom_sf(data = sites_sf,
          size = 4,
          color = "black") +
  geom_sf_label(data = sites_sf,
                aes(label = site_code)) +
  geom_sf_label(data = sites_sf %>%
                  filter(site_code == root_site),
                aes(label = site_code),
                color = "red") +
  theme_bw() +
  theme(axis.title = element_blank()) +
  labs(color = "Stream\nOrder")


#-----------------------------------------------------------------
# build parent child table
parent_child = sites_sf %>%
  buildParentChild(flowlines,
                   rm_na_parent = T,
                   add_rkm = F)

# add RKMs from configuration file (since we had to fix at least one from PTAGIS)
parent_child %<>%
  left_join(configuration %>%
              select(parent = site_code,
                     parent_rkm = rkm) %>%
              distinct(),
            by = "parent") %>%
  left_join(configuration %>%
              select(child = site_code,
                     child_rkm = rkm) %>%
              distinct(),
            by = "child") %>%
  distinct()



# confirm each site has the correct path for a PIT tag
parent_child %>%
  buildPaths()

#-----------------------------------------------------------------
# Save file.
save(configuration,
     sites_sf,
     flowlines,
     parent_child,
     file = here('analysis/data/derived_data/site_config.rda'))


#-----------------------------------------------------------------
# Build network diagram
# simple
pc_graph = plotNodes(parent_child,
                     layout = "tree")

pc_nodes_graph = parent_child %>%
  addParentChildNodes(configuration) %>%
  plotNodes()

# control more settings
node_order = buildNodeOrder(parent_child)

nodes = buildNodeGraph(parent_child) %>%
  as_tibble()

edges = parent_child %>%
  left_join(nodes, by = c('parent' = 'label')) %>%
  rename(from = index) %>%
  left_join(nodes, by = c('child' = 'label')) %>%
  rename(to = index) %>%
  select(from, to)

library(ggraph)
node_graph = tidygraph::tbl_graph(nodes = nodes,
                                  edges = edges)

node_p = node_graph %>%
  ggraph(layout = "tree") +
  # ggraph(layout = "partition") +
  # ggraph(layout = "kk") +
  geom_edge_link(arrow = arrow(length = unit(2, 'mm'),
                               type = "closed"),
                 end_cap = circle(4, 'mm')) +
  geom_node_point(size = 7) +
  theme_graph(base_family = 'Times') +
  theme(legend.position = 'none') +
  scale_color_brewer(palette = "Set1",
                     na.value = "black") +
  geom_node_label(aes(label = label),
                  size = 2,
                  label.padding = unit(0.1, 'lines'),
                  label.size = 0.1)

node_p

# save as pdf
library(here)
ggsave(here("analysis/figures/Asotin_DABOM_sites.pdf"),
       node_p,
       width = 9,
       height = 6)
