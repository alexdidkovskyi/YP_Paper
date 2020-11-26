# !diagnostics off




#' Initial data processing.
#' Including
#' a) Elaboration of the dataset for classification & regression
#' b) Elaboration of the supplementeray datasets with coordinates and subdomain types'

options(max.print = 100)

library(tidyverse)
library(magrittr)
library(dummies)
library(openxlsx)
library(gstat)
library(sf)
library(sp)


### Scaling functions
resc_ <- function(x, x_min, x_max) {
  (x - x_min) / (x_max - x_min)
}
inv_resc_ <- function(sc_x, x_min, x_max) {
  x_min + sc_x * (x_max - x_min)
}


list2env(readRDS('data/preprocessed/preprocessed_observations.RDS'),
         globalenv())

baseline_weather <-
  readRDS('data/preprocessed/baseline_weather_df.RDS')



#' Test data: locations generation
#' Select 3000 locations from the river domain
pools_shape <-
  rgdal::readOGR(dsn = "data/initial/shp", layer = "MORFOLOGIA_OKSHP_WGS84")

area_pol <- raster::area(pools_shape)
pool_l <- length(pools_shape@polygons)
total_points_num <- 3000
points_num <- round(area_pol / sum(area_pol) * total_points_num, 0)

pool_grid <- list()

for (i in 1:pool_l) {
  pool_bnd <- pools_shape@polygons[[i]]@Polygons[[1]]@coords
  pool_grid <-  c(
    pool_grid,
    list(
      data.frame(x = (runif(
        3000, min(pool_bnd[, 1]), max(pool_bnd[, 1])
      )),
      y = (runif(
        3000, min(pool_bnd[, 2]), max(pool_bnd[, 2])
      ))) %>%
        rowwise() %>%
        mutate(in_out = mgcv::in.out(pool_bnd, c(x, y))) %>%
        mutate(poly_id = i) %>%
        filter(in_out == T) %>%
        .[1:points_num[i],] %>%
        mutate(TipoStart =  as.character(pools_shape@data$TIPO[i]))
    )
  )
}

pool_grid <-
  do.call('rbind', pool_grid) %>% data.frame(., stringsAsFactors = F)
pool_grid %<>% bind_cols(dummies::dummy.data.frame(pool_grid %>% dplyr::select(TipoStart)))

pool_grid %<>% distinct() %>% na.omit()

CRS.n <-
  CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
pools_shape_diff <- spTransform(pools_shape, CRS.n)
pool_grid[, c("x_WGS84", "y_WGS84")] <-
  proj4::project(pool_grid[, c("x", "y")], CRS.n)

#' data processing:'is_stuck' features generation
df_loc <- df_loc %>% group_by(IDREF) %>% arrange(EventoStart) %>%
  mutate(is_stuck_1 = (lag(EventoEnd, default = -1) == EventoStart) &
           (lag(Distance_m, default = -1) == 0)) %>%
  mutate(is_stuck_2 = (lag(EventoEnd, 2, default = -1) + 1 == EventoStart) &
           (lag(Distance_m, default = -1) == 0) &
           (lag(Distance_m, 2, default = -1) == 0)) %>%
  mutate(is_stuck_3 = (lag(EventoEnd, 3, default = -1) + 2 == EventoStart) &
           (lag(Distance_m, default = -1) == 0) &
           (lag(Distance_m, 2, default = -1) == 0) &
           (lag(Distance_m, 3, default = -1) == 0)) %>%
  ungroup()



#' Fix subdomain issues.
df_loc$TipoStart[df_loc$poly_id == 3] <- 'Run_rapid'
df_loc$TipoStart[df_loc$poly_id == 7] <- 'Planebed'
df_loc$TipoStart[df_loc$poly_id == 23] <- 'Cascade'
df_loc$TipoStart[df_loc$poly_id == 39] <- 'Steep_poolzone'
df_loc$TipoStart[df_loc$poly_id == 9] <- 'Banks'


#' Remove backward movement
df_loc <- df_loc %>% filter(X_Start > X_End |
                              Y_Start > Y_End |
                              graph_dist <= 1) %>%
  filter(EventoStart == EventoEnd - 1)




df_models_ <-
  df_ch %>% 
  left_join(df_loc) %>% 
  mutate(Event = EventoEnd) %>% 
  left_join(baseline_weather %>% 
              rename(Event = EventoStart)) %>%
  drop_na()


df_models_ <-
  df_models_ %>%
  bind_cols(dummies::dummy.data.frame(df_models_ %>% 
                                        dplyr::select(TipoStart))) %>%
  mutate(graph_vel = graph_dist / duration,
         id_o = 1:nrow(.))



#' Basic datasets generation
#' Features dataset
df_ <-
  df_models_ %>% dplyr::select(
    -TipoEnd,-TipoStart,-X_Start,-Y_Start,-X_End,-Y_End,-ID,-IDREF,-IDNUM,-EventoStart,-Event,
    -EventoEnd,-duration,-graph_vel,-Distance_m,-starts_with('Tipo'),-poly_id,-graph_dist
  )

#' Target dataset generation
df_pred_ <-
  df_models_ %>% dplyr::select(duration, graph_dist, graph_vel, id_o)


#' Locations dataset generation
df_coordinates_and_type_ <-
  df_models_ %>% dplyr::select(
    TipoEnd,
    TipoStart,
    X_Start,
    Y_Start,
    X_End,
    Y_End,
    EventoStart,
    Event,
    poly_id,
    EventoEnd,
    id_o,
    starts_with('Tipo')
  )

df_coordinates_and_type_[, c('X_Start_WGS84', 'Y_Start_WGS84')] <-
  proj4::project(df_coordinates_and_type_[, c("X_Start", "Y_Start")], CRS.n)

x_min <- df_coordinates_and_type_$X_Start_WGS84 %>% min()
x_max <- df_coordinates_and_type_$X_Start_WGS84 %>% max()

y_min <- df_coordinates_and_type_$Y_Start_WGS84 %>% min()
y_max <- df_coordinates_and_type_$Y_Start_WGS84 %>% max()

df_coordinates_and_type_$X_Start_WGS84_scaled <-
  resc_(df_coordinates_and_type_$X_Start_WGS84, x_min, x_max)
df_coordinates_and_type_$Y_Start_WGS84_scaled <-
  resc_(df_coordinates_and_type_$Y_Start_WGS84, y_min, y_max)

pool_grid$x_WGS84_scaled <- resc_(pool_grid$x_WGS84, x_min, x_max)
pool_grid$y_WGS84_scaled <- resc_(pool_grid$y_WGS84, y_min, y_max)

#' Test data generation. 6 cases: low, medium, high elongation x low/high weight

#(size, elongation, platyness)
#(bhh, shh,bhm, shm, bmm,smm)
id_0_case <- c(1089, 1025, 1159, 1679, 1640 , 1238)

df_test_peb_ <- df_[id_0_case,] %>%
  select(-colnames(baseline_weather)[3:7], -id_o) %>%
  bind_cols(baseline_weather %>%
              filter(EventoStart == 15) %>%
              .[rep(1, 6), 2:7])



#' Calculate covered distance using graph distance between points
pools_shape_polygon_union <-
  sp::polygons(pools_shape) %>%  rgeos::gUnaryUnion(.)

pool_border <-
  pools_shape_polygon_union@polygons[[1]]@Polygons[[1]]@coords %>%
  data.frame() %>%
  round(., 8) %>%
  distinct()

pebbles_loc_st <- df_coordinates_and_type_ %>%
  select('X_Start', 'Y_Start') %>%
  dplyr::distinct() %>% set_colnames(c('x', 'y')) %>%
  distinct() %>%
  anti_join(pool_border)

coords_for_triang <-
  rbind(pool_border, pebbles_loc_st) %>% distinct()

m_border <- matrix(c(1:nrow(pool_border),
                     2:nrow(pool_border), 1),
                   nrow(pool_border),
                   2)


data_plsg <- RTriangle::pslg(P = coords_for_triang, S = m_border)

data_traingl <- RTriangle::triangulate(data_plsg, a = 1e-11)

edge.list <- data_traingl$E
edge.num <- dim(edge.list)[1]

graph <- igraph::graph_from_edgelist(edge.list, directed = FALSE)

edge.weights <- apply(edge.list, 1, function(i) {
  geosphere::distHaversine(data_traingl$P[i[1], 1:2],
                           data_traingl$P[i[2], 1:2])
  
})

triangl_dist <- igraph::distances(
  graph,
  v = igraph::V(graph),
  to = igraph::V(graph),
  mode = c("all"),
  weights = edge.weights,
  algorithm = c("automatic")
)


st_ <- nrow(pool_border) + 1
en_ <- nrow(pool_border) + nrow(pebbles_loc_st)

triangl_dist_btw_points <- triangl_dist[st_:en_, st_:en_]

rm(triangl_dist, st_, en_, edge.weights, edge.num)




data_list <- list(
  dfs = list(
    df_ = df_,
    df_pred_ = df_pred_,
    df_coordinates_and_type_ = df_coordinates_and_type_
  ),
  dfs_test = list(
    df_test_peb_ = df_test_peb_,
    pool_grid = pool_grid,
  ),
  shape = list(pools_shape = pools_shape,
               pools_shape_WGS84 = pools_shape_diff),
  triangl_dist_btw_points = triangl_dist_btw_points
)


#' Save all datasets as .RDS
saveRDS(data_list, file = 'data/preprocessed/all_datasets.RDS')
rm(list = ls())
