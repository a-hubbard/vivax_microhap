# Create relatedness network map

# Load required libraries ----------------------------------------------
library(ggraph)
library(rnaturalearth)
library(rnaturalearthhires)
library(tidygraph)
# These libraries will be referenced without the library name and so 
# should be loaded second
library(magrittr)
library(optparse)                                                                
library(sf)
library(tidyverse)

# Parse arguments ------------------------------------------------------
opts <- list(
  make_option("--rel", help = "CSV file containing relatedness values"), 
  make_option("--site_coords", help = "GeoJSON containing site coordinates"), 
  make_option("--out_base", help = "Basename for output figure")
)
arg <- parse_args(OptionParser(option_list = opts))
# Arguments used for development
# arg <- list(
#   rel = "../../results/relatedness/site_rel.csv", 
#   site_coords = "../../results/unfiltered_sites.geojson"
# )

# Read data ------------------------------------------------------------
site_rel <- read_csv(
  arg$rel, 
  col_types = cols(
    .default = col_character(), 
    mean_r = col_double(), 
    frac_signif = col_double(), 
    frac_high_r = col_double()
  ), 
  progress = FALSE
)

# Create network as tbl_graph ------------------------------------------
# Read in node coordinates
node_coords <- st_read(arg$site_coords) %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
  as_tibble() %>%
  select(-geometry)
nodes <- site_rel %>%
  filter(site_a == site_b) %>%
  select(-site_b) %>%
  rename(site = site_a) %>%
  left_join(node_coords, by = c(site = "Site"))
edges <- site_rel %>%
  filter(site_a != site_b) %>%
  # This won't be a directed network, so this is arbitrary
  rename(from = site_a, to = site_b) %>%
  left_join(node_coords, by = c(from = "Site")) %>%
  rename(from_lon = lon, from_lat = lat) %>%
  select(-n_samp, -Country) %>%
  left_join(node_coords, by = c(to = "Site")) %>%
  rename(to_lon = lon, to_lat = lat) %>%
  select(-n_samp, -Country)
network <- tidygraph::tbl_graph(
  nodes = nodes, 
  edges = edges, 
  directed = FALSE, 
  node_key = "site"
)

# Compute map bounds ---------------------------------------------------
node_points <- network %>%
  tidygraph::activate(nodes) %>%
  as_tibble() %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs(4326)
extent <- st_bbox(node_points)
xrng <- extent[["xmax"]] - extent[["xmin"]]
yrng <- extent[["ymax"]] - extent[["ymin"]]
xbuff <- xrng * 0.1
ybuff <- yrng * 0.1
bnds <- c(
  xmin = extent[["xmin"]] - xbuff, 
  xmax = extent[["xmax"]] + xbuff, 
  ymin = extent[["ymin"]] - ybuff, 
  ymax = extent[["ymax"]] + ybuff
)

# Create and save map --------------------------------------------------
countries <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
fig <- ggraph::ggraph(network) +
  geom_sf(data = countries) +
  annotate(
    "text", 
    x = 108.3, 
    y = 14, 
    label = "Vietnam", 
    size = 5, 
    fontface = "bold"
  ) +
  annotate(
    "text", 
    x = 105.5, 
    y = 13, 
    label = "Cambodia", 
    size = 5, 
    fontface = "bold"
  ) +
  ggraph::geom_edge_arc(
    mapping = aes(
      x = from_lon, 
      y = from_lat, 
      xend = to_lon, 
      yend = to_lat, 
      color = frac_high_r
    ), 
    width = 1.5
  ) +
  ggraph::geom_node_point(
    mapping = aes(x = lon, y = lat, size = n_samp), 
    color = "black", 
    fill = "white", 
    shape = 21
  ) +
  geom_sf_label(
    mapping = aes(label = site), 
    data = node_points, 
    nudge_x = 0.2, 
    nudge_y = -0.2
  ) +
  coord_sf(
    xlim = c(bnds[["xmin"]], bnds[["xmax"]]), 
    ylim = c(bnds[["ymin"]], bnds[["ymax"]])
  ) +
  ggraph::scale_edge_color_distiller(
    name = "Frac. High r", 
    palette = "YlOrBr", 
    direction = 1
  ) +
  labs(size = "No. of Samples") +
  theme(
    legend.position = "bottom", 
    panel.background = element_rect(fill = "lightskyblue1"), 
    panel.border = element_rect(
      color = "black", 
      fill = NA, 
      linewidth = 0.7
    )
  )
w <- 6.5
h <- 5.7
ggsave(
  str_c(arg$out_base, ".pdf"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
ggsave(
  str_c(arg$out_base, ".png"), 
  plot = fig, 
  width = w, 
  height = h, 
  units = "in"
)
