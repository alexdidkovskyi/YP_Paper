# !diagnostics off


#' Initial data Preprocessing.
#' Including
#' a) Processing of observations gathered during 2016, 2017 & 2018
#' b) Data cleaning 
#' c) Covered distance recalculation

options(java.parameters = "-Xmx4096m")

library(openxlsx)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(geosphere)
library(sp)
library(rgdal)
library(proj4)
library(raster)
library(igraph)


#' Pebbles' locations preprocessing. 2016-2017
#'
iter <- 0
df_list <- list()
for (i in list.files('data/initial/Rilievo2016')) {
  iter <- iter + 1
  df_list[[i]] <-
    read_csv(paste0('data/initial/Rilievo2016/', i)) %>%
    mutate(iter = iter)
}
filenames_df_16_17 <-
  data.frame(name = list.files('data/initial/Rilievo2016'),
             var = 1:length(iter))
df_locations_16_17 <- do.call('rbind', df_list)


#' Pebbles' locations preprocessing. 2016-2017
df_characteristics_16_17 <-
  read.xlsx('data/initial/Pebbles_16_17.xlsx')
df_characteristics_16_17[,-1] %<>% mutate_if(is.character, as.numeric)



#' Pebbles' locations preprocessing. 2018
df_locations_18 <- read.xlsx('data/initial/MERGED2018_XLS.xlsx')

CRS.n <-
  CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')


#' locations' projection to WGS84
df_locations_18[, 9:10] <-
  data.frame(project(df_locations_18[, 9:10], CRS.n, inverse = TRUE))
df_locations_18[, 11:12] <-
  data.frame(project(df_locations_18[, 11:12], CRS.n, inverse = TRUE))



#' Pebbles' dimensions preprocessing. 2018
df_characteristics_18 <-
  read.xlsx('data/initial/DB 2018.xlsx', sheet = 'Foglio1') %>%
  set_colnames(c('IDREF', 'Weight', 'a_axis', 'b_axis', 'c_axis')) %>%
  #' remove outliers: a_axis > 200, b_axis > 135
  mutate(
    a_axis = ifelse(a_axis > 200, 200, a_axis),
    b_axis = ifelse(b_axis > 135, 135, b_axis)
  ) %>%
  mutate(
    Nominal_diameter = (a_axis * b_axis * c_axis) ^ (1 / 3),
    Elongation = b_axis / a_axis,
    Platyness = c_axis / a_axis,
    Sphericity = (b_axis * c_axis / a_axis ^ 2) ^ (1 / 3)
  )



#' Aggregate and filter pebbles' locations and characteristics. 2018
df_loc <- df_locations_16_17[,-13] %>% bind_rows(df_locations_18)
df_loc[, 9:12] <- df_loc[, 9:12] %>% round(., 8)
df_ch  <-
  df_characteristics_16_17[,-10] %>% bind_rows(df_characteristics_18)

good_IDREF <- intersect(df_loc$IDREF, df_ch$IDREF)

#' Remove observations with positioning errors
df_loc %<>% dplyr::filter(IDREF %in% good_IDREF) %>%
  filter(Y_Start >= 45, X_Start >= 9,
         Y_End >= 45, X_End >= 9)
df_ch %<>% dplyr::filter(IDREF %in% good_IDREF)






#' Rearrange pebbles' dimensions
#' Fix cases where a_axis < b_axis, etc.
id_a_b <- which(df_ch$a_axis < df_ch$b_axis)

temp_ <- df_ch$a_axis[id_a_b]
df_ch$a_axis[id_a_b] <- df_ch$b_axis[id_a_b]
df_ch$b_axis[id_a_b] <- temp_


id_b_c <- which(df_ch$b_axis < df_ch$c_axis)
temp_ <- df_ch$a_axis[id_a_c]
df_ch$a_axis[id_a_c] <- df_ch$c_axis[id_a_c]
df_ch$c_axis[id_a_c] <- temp_



id_a_c <- which(df_ch$a_axis < df_ch$c_axis)
temp_ <- df_ch$b_axis[id_b_c]
df_ch$b_axis[id_b_c] <- df_ch$c_axis[id_b_c]
df_ch$c_axis[id_b_c] <- temp_


id_a_b <- which(df_ch$a_axis < df_ch$b_axis)

temp_ <- df_ch$a_axis[id_a_b]
df_ch$a_axis[id_a_b] <- df_ch$b_axis[id_a_b]
df_ch$b_axis[id_a_b] <- temp_



df_ch %<>% mutate(
  Nominal_diameter = (a_axis * b_axis * c_axis) ^ (1 / 3),
  Elongation = b_axis / a_axis,
  Platyness  = c_axis / a_axis,
  Sphericity = (b_axis * c_axis / a_axis ^ 2) ^ (1 /
                                                   3)
)


#df_loc %>% mutate(E = EventoStart- EventoEnd) %>% dplyr::select(E) %>% pull() %>% table()


#' Add missed events for cases where a pebble didn't move (covered distance is 0) and EventoStart, EventoEnd are not consecutive
df_loc_temp <- df_loc %>% filter(EventoEnd - EventoStart >= 2,
                                 Distance_m == 0)

list_new <- list()
for (j in 1:nrow(df_loc_temp)) {
  temp_ <- df_loc_temp[j, ]
  if (dim(temp_)[1] >= 1) {
    a <- apply(temp_, 1, function(x) {
      .df <-
        data.frame(EventoStart = c(as.numeric(x[4]):(as.numeric(x[5]) - 1))) %>% mutate(EventoEnd = EventoStart + 1)
    })
    
    for (i in nrow(temp_)) {
      temp__ <- temp_[rep(i, nrow(a[[i]])),]
      temp__[, c(4, 5)] <- a[[i]]
      list_new[[paste0(j, '_', i)]] <- temp__
    }
    
  }
  
}

df_loc_ <- df_loc %>% filter(EventoEnd - EventoStart <= 1) %>%
  bind_rows(do.call('rbind', list_new)) %>%
  bind_rows(df_loc %>% filter(EventoEnd - EventoStart >= 2,
                              Distance_m > 0))

df_loc <- df_loc_


df_ch %<>% bind_cols(data.frame((
  df_ch %>% dplyr::select(a_axis, b_axis, c_axis) %>% princomp()
)$scores[, 1:2]) %>%
  set_colnames(c('pebble_PC1', 'pebble_PC2')))




#' Adjust covered distance

#' Provided distance is a length of the straight line between two points.
#' Hence it needs to be adjusted accordingly with the direction of the flow
#' For this aim the domain is triangulated and then the distances between locations are estimated as distances on the graph (between graph nodes)


pools_shape <-
  rgdal::readOGR(dsn = "data/initial/shp", layer = "MORFOLOGIA_OKSHP_WGS84")
pools_shape_polygon_union <-
  polygons(pools_shape) %>%  rgeos::gUnaryUnion(.)

pool_border <-
  pools_shape_polygon_union@polygons[[1]]@Polygons[[1]]@coords %>%
  data.frame() %>%
  round(., 8) %>%
  distinct()

pool_intra <- lapply(pools_shape@polygons,
                     function(x)
                       x@Polygons[[1]]@coords) %>%
  do.call('rbind', .) %>%
  unique()  %>% data.frame() %>%
  set_colnames(c('x', 'y')) %>%
  round(., 8) %>%
  distinct() %>%
  anti_join(., pool_border, by = c('x', 'y')) %>%
  distinct()



#' Triangulation

coords_for_triang <- rbind(pool_border, pool_intra) %>% distinct()

pebbles_loc_st <- df_loc %>%
  dplyr::select(X_Start, Y_Start) %>%
  as.data.frame() %>%
  dplyr::distinct() %>% set_colnames(c('x', 'y')) %>%
  distinct() %>%
  anti_join(coords_for_triang)

pebbles_loc_end <- df_loc %>%
  dplyr::select(X_End, Y_End) %>%
  as.data.frame() %>%
  dplyr::distinct() %>% set_colnames(c('x', 'y')) %>%
  distinct() %>%
  anti_join(coords_for_triang)

coords_for_triang <-
  rbind(coords_for_triang,
        rbind(pebbles_loc_st, pebbles_loc_end) %>% distinct())

m_border <- matrix(c(1:nrow(pool_border),
                     2:nrow(pool_border), 1),
                   nrow(pool_border),
                   2)


data_plsg <- RTriangle::pslg(P = coords_for_triang, S = m_border)

data_traingl <- RTriangle::triangulate(data_plsg, a = 1e-11)

png('triang.png')
plot(data_traingl)
dev.off()


edge.list <- data_traingl$E
edge.num <- dim(edge.list)[1]

graph <- igraph::graph_from_edgelist(edge.list, directed = FALSE)


#' weight of each edge is set to the haversine distance between two associated nodes.
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




start_coord <- df_loc %>% dplyr::select(X_Start, Y_Start)
end_coord <- df_loc %>% dplyr::select(X_End, Y_End)

start_id <- apply(start_coord, 1, function(x) {
  which(x[1] == coords_for_triang[, 1] &
          x[2] == coords_for_triang[, 2])
})

end_id <- apply(end_coord, 1, function(x) {
  which(x[1] == coords_for_triang[, 1] &
          x[2] == coords_for_triang[, 2])
})
options("digits" = 9, scipen = 100)


library(mgcv)
in_out <-
  sapply(1:nrow(coords_for_triang), function(r_)
    mgcv::in.out(as.matrix(pool_border), as.numeric(coords_for_triang[r_, ])))


in_out_id <- setdiff(1:nrow(data_traingl$P), which(in_out == F))
df_loc$graph_dist <-
  sapply(1:nrow(df_loc), function(id)
    triangl_dist[start_id[id], end_id[id]])
id_inf <- which(df_loc$graph_dist %>% is.infinite())
df_loc$graph_dist[id_inf] <- df_loc$Distance_m[id_inf]


pools_shape <-
  rgdal::readOGR(dsn = "data/initial/shp", layer = "MORFOLOGIA_OKSHP_WGS84")


polygon_id <-  sapply(1:nrow(df_loc), function(j) {
  which(sapply(1:45, function(i)
    sp::point.in.polygon(
      df_loc[j, 9],
      df_loc[j, 10],
      pools_shape@polygons[[i]]@Polygons[[1]]@coords[, 1],
      pools_shape@polygons[[i]]@Polygons[[1]]@coords[, 2]
    )) > 0)[1]
}) %>% unlist()

df_loc$poly_id <- polygon_id
rm(list = setdiff(ls(), c('df_loc', 'df_ch')))


#' Flow data processing

df_list <- list()
for (i in 3:30) {
  df_list[[i - 2]] <-
    xlsx::read.xlsx(
      file = 'data/initial/Dati idrometeo.xlsx',
      sheetIndex  = i,
      startRow = 3,
      colClasses = 'character',
      stringsAsFactors = F
    ) %>%
    .[, 2:10] %>%
    mutate(EventoStart = i - 2)
  gc()
}


df_weather <- do.call('rbind', df_list) %>%
  set_colnames(gsub('\\.+', '_', colnames(.))) %>%
  dplyr::select(-h_CP_m_)

rm(i, df_list)


#' Safe processed data as .RDS

saveRDS(
  list(
    df_ch = df_ch,
    df_loc = df_loc,
    df_weather = df_weather
  ),
  file = 'data/preprocessed/preprocessed_observations.RDS',
  compress = F
)

rm(list = ls())
