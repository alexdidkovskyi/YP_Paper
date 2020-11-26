
#' Analysis of red pebbles data.
#' Comparison of XGB and FK on the red pebbles data

options(max.print = 100)

library(tidyverse)
library(magrittr)
library(xgboost)
library(openxlsx)
library(sp)
library(fdagstat)

water_depth_threshold <-  35

source('src/AUX_01_models.R')
source('src/AUX_02_functional_models.R')


list2env(readRDS('data/preprocessed/preprocessed_observations.RDS'),
         globalenv())


baseline_weather <-
  readRDS('data/preprocessed/baseline_weather_df.RDS')



data_list <- readRDS('data/preprocessed/all_datasets.RDS')


list2env(data_list$dfs, globalenv())
list2env(data_list$dfs_test, globalenv())
list2env(data_list$shape, globalenv())


#' XGB training
df_$duration <- df_pred_$duration
df__ <-  df_
df_coordinates_and_type__ <-  df_coordinates_and_type_
df_pred__ <- df_pred_
df_$y_ <- df_pred_$graph_dist > 0

clf_data_ <- df_ %>%
  left_join(df_coordinates_and_type_ %>% select(starts_with('TipoStart'),
                                                id_o,
                                                X_Start,
                                                Y_Start)) %>%
  arrange(id_o) %>%
  select(-TipoStart,-id_o) %>%
  filter(cl_ == 1)

y_ <- clf_data_$y_
clf_data_ <-
  clf_data_ %>% select(X_Start, Y_Start, a_axis, b_axis, mean_h, mean_Q, duration)


xgb_train <-
  xgboost::xgb.train(
    params = list(objective = 'binary:logistic', eta = 0.03),
    nrounds = 210,
    data = xgboost::xgb.DMatrix(data = clf_data_ %>% data.matrix(),
                                label = y_),
    verbose = 0
  )

#' FK Training
triangl_dist_btw_points <- data_list$triangl_dist_btw_points

df_coordinates_and_type_$cl_ <- df_$cl_


a <-
  df_coordinates_and_type_ %>% group_by(X_Start, Y_Start) %>% summarise(is_ = 1 %in% cl_, )

b <-
  df_coordinates_and_type_ %>% dplyr::select(X_Start,
                                             Y_Start,
                                             X_Start_WGS84,
                                             Y_Start_WGS84,
                                             poly_id,
                                             TipoStart) %>% distinct()

c_ <- b %>% left_join (a)

dist_matrix_ <-
  triangl_dist_btw_points[which(c_$is_ == 1), which(c_$is_ == 1)]
locations_ <- c_ %>% filter(is_ == 1) %>% mutate(cl_ = 1)
rm(a, b, c_)

targets_ <- locations_ %>%
  left_join(
    df_pred_ %>%
      left_join(df_coordinates_and_type_) %>%
      left_join(df_) %>%
      dplyr::select(
        X_Start,
        Y_Start,
        graph_dist,
        pebble_PC1,
        cl_,
        duration,
        weather_PC1
      )
  ) %>%
  mutate(y_ = as.numeric(graph_dist > 0))


targets_count_ <-
  targets_ %>% dplyr::select(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
  group_by(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
  summarise(count_ = n()) %>% ungroup() %>%
  mutate(id__ = 1:nrow(.))



crvs_params_ <- list(
  type_ = 'domain',
  r_ = 5,
  range_ = 50,
  nr_obs_ = 12,
  lower_bnd_ = 0.01,
  upper_bnd_ = 0.99,
  bandwidth_ = 20,
  max_points_num = 30,
  range_l_ = -50,
  range_r_ = 50,
  range_step_ = 1,
  train_var = 'pebble_PC1',
  target_var = 'y_'
)

type_vgm_ <- 'Bes'
vName_ <- 'Movement Probability'

df_func_ <- generate_curves(
  locations_ = locations_,
  targets_ =  targets_,
  dist_matrix_ = dist_matrix_,
  crvs_params_ = crvs_params_ ,
  targets_count_ = targets_count_
)

df_func_ <- df_func_$df_func
na_count_ <- apply(df_func_, 2, function(x)
  sum(is.na(x)))

id_0 <- which(na_count_ == 0)
df_func_subset_ <- df_func_[, id_0]
coord_ <-
  locations_[id_0, c('X_Start_WGS84', 'Y_Start_WGS84')] %>% set_colnames(c('x', 'y'))




g_UTrK_fit_ <- functional_kriging_(
  coord_ = coord_,
  func_ = df_func_subset_,
  type_vgm_ = type_vgm_,
  vName_ = vName_,
  Nlags_ = 10,
  LagMax_ = 8,
  ArgStep_ = 1,
  intercept_ = F,
  type_drift_ = 'OLS',
  regCoef_ = 1e1,
  plot_ = F
)



#' Red pebbles data processing

sheet_names <-
  openxlsx::getSheetNames('data/initial/dati sassi rossi.xlsx')
list_of_red_data <- list()
for (sheet_name in sheet_names) {
  list_of_red_data[[sheet_name]] <-
    read.xlsx('data/initial/dati sassi rossi.xlsx', sheet = sheet_name)
}

ch_pca <- df_ch[, 2:4] %>% princomp()


Location_Event <-
  data.frame(
    location = c("R1",  "R2",  "R4", "R5",  "R6",  "R8",  "R9"),
    Event    = c(11, 11, 14, 14, 16, 16, 16),
    x  = c(
      9.426869334,
      9.426570343,
      9.426569269,
      9.426353789,
      9.426863182,
      9.426735219,
      9.426352559
    ),
    y  = c(
      45.874334298,
      45.874033173,
      45.874027936,
      45.873826312,
      45.874321206,
      45.874166715,
      45.873828931
    )
  )


Location_Event_ <-
  Location_Event %>% left_join(baseline_weather %>% rename(Event = EventoStart))




sapply(1:nrow(Location_Event_),
       function(id_) {
         which(sapply(1:45,
                      function(id__)
                        mgcv::in.out(
                          pools_shape@polygons[[id__]]@Polygons[[1]]@coords,
                          c(Location_Event$x[id_], Location_Event$y[id_])
                        )
                      == T))
       }) -> tipo_id

Location_Event_$Tipo <- pools_shape$TIPO[tipo_id]

baseline_weather


CRS.n <-
  CRS('+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
library(raster)
Location_Event[, c("x_WGS84", "y_WGS84")] <-
  proj4::project(Location_Event[, c("x", "y")], CRS.n)
detach("package:raster", unload = TRUE)

transform_red_data <- function(data) {
  data <-
    data %>%
    drop_na() %>%
    mutate_all(as.numeric) %>%
    set_colnames(c('a_axis', 'b_axis')) %>%
    mutate(c_axis = 3 / 4 * b_axis,
           diam_ = sqrt(a_axis * b_axis))
  
  data$pebble_PC1 <-
    t(((t(data[, c(1, 2, 3)]) - ch_pca$center) * ch_pca$loadings[, 1])) %>% rowSums()
  
  return(data)
}



l_ <- list()
i <- 1

for (i in 1:7) {
  pre__$diam_
}

d_ <- c()
for (i in 1:7) {
  pre__ <- transform_red_data(list_of_red_data[[i]][-1, 2:3])
  post__ <-
    transform_red_data(list_of_red_data[[i]][-1, 5:6]) %>% drop_na()
  
  d_ <- c(d_, pre__$diam_)
  
  #print(c(min(pre__$diam_), max(pre__$diam_)))
  #print(c(min(post__$diam_), max(post__$diam_)))
  
}


d_ %>% max()
d_ch <-
  df_ch %>% mutate(diam__ = sqrt(a_axis * b_axis)) %>% pull(diam__)

min(d_ch)
mean(d_ >= min(d_ch) - 20)


#' Results elaboration
for (i in 1:7) {
  n_binom <- 100
  pre__ <- transform_red_data(list_of_red_data[[i]][-1, 2:3])
  post__ <-
    transform_red_data(list_of_red_data[[i]][-1, 5:6]) %>% drop_na()
  
  
  
  df_for_xgb <- cbind(
    Location_Event_[rep(i, length(pre__$pebble_PC1)),] %>%
      select(x, y, duration, mean_h, mean_Q) %>%
      rename(X_Start = x, Y_Start = y),
    pre__ %>% select(a_axis, b_axis)
  ) %>% select(X_Start, Y_Start, a_axis, b_axis, mean_h, mean_Q, duration)
  
  xgb_pred_ <- predict(xgb_train, df_for_xgb %>% data.matrix())
  
  xgb_realisations_ <-
    sapply(xgb_pred_, function(p)
      rbinom(n_binom, 1, p))
  
  xgb_mean_w_ <-
    sapply(1:n_binom, function(row_) {
      if (c(0) %in% unique(xgb_realisations_[row_,])) {
        obs_ <- pre__$diam_[xgb_realisations_[row_, ] == 0]
        
        W_ <-  transport::wasserstein1d(obs_, post__$diam_, p = 2)
        return(W_)
      }
    }) %>% unlist() %>% mean(., na.rm = T)
  
  
  xgb_mean_nm_ <-
    apply(xgb_realisations_, 1, function(x)
      sum(x == 0)) %>% mean()
  
  F_ <-
    predictFstat(
      g_UTrK_fit_,
      .newCoordinates = Location_Event[i, 5:6] %>% set_colnames(c('x', 'y')),
      .type = 'OK',
      .what = vName_
    )$Forecast[, 1]
  F_f_ <- approxfun(x = -50:50,
                    y = F_,
                    method = 'linear')
  
  pre__pc <- pre__$pebble_PC1
  pre__pc[pre__pc < -50] <- -50
  
  pre__pc[pre__pc >= 50] <- 50
  FK_pred_ <- F_f_(pre__pc)
  
  
  FK_realisations_ <-
    sapply(FK_pred_, function(p)
      rbinom(n_binom, 1, p))
  
  FK_mean_w_ <-
    sapply(1:n_binom, function(row_) {
      if (c(0) %in% unique(FK_realisations_[row_,])) {
        obs_ <- pre__$diam_[FK_realisations_[row_, ] == 0]
        
        W_ <-  transport::wasserstein1d(obs_, post__$diam_, p = 2)
        return(W_)
      }
    }) %>% unlist() %>% mean(., na.rm = T)
  
  
  FK_mean_nm_ <-
    apply(FK_realisations_, 1, function(x)
      sum(x == 0)) %>% mean()

  l_[[i]] <- c(
    nrow(pre__),
    nrow(post__),
    nrow(post__) / nrow(pre__),
    xgb_mean_nm_,
    xgb_mean_nm_ / nrow(pre__),
    xgb_mean_w_,
    FK_mean_nm_,
    FK_mean_nm_ / nrow(pre__),
    FK_mean_w_,
    sum(pre__$pebble_PC1 > -50)
  )
}


l_ <-
  l_ %>% do.call(rbind, .) %>% data.frame() %>% set_colnames(
    c(
      'Total',
      'Not moved',
      'Proportion Nm',
      'Xgb Not moved',
      'Proportion Xgb Nm',
      'Xgb wasserstein',
      'FK not moved',
      'Proportion FK Nm',
      'FK wasserstein',
      'NA'
    )
  )

colnames(l_)[10] <- 'NA'
l_[, 3:9] <- round(l_[, 3:9], 3)

l_ %>% select(`FK wasserstein`, `Xgb wasserstein`)
l_$`FK wasserstein` - l_$`Xgb wasserstein`


df_for_vis <-
  cbind(Location_Event, l_) %>% mutate(Event = as.character(Event))
df_for_vis$`FK wasserstein`

ggplot() +
  geom_polygon(
    pools_shape,
    mapping  = aes(x = long, y = lat, group = group),
    colour = 'black',
    fill = NA,
    size = 0.7,
    alpha = 0.5
  ) +
  geom_point(data = df_for_vis,
             aes(
               x = x,
               y = y,
               shape = Event,
               col = `FK wasserstein`
             ),
             size = 5) +
  scale_color_viridis_c('Average \nWasserstein \ndistance',
                        alpha = 0.8,
                        limits = c(4, 42)) + theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16)
  ) +
  xlab('longitude') +
  ylab('latitude') + guides(colour = guide_colourbar(order = 1),
                            shape = guide_legend(order = 2))

ggsave(
  filename = paste0('Wasserstein_FK', '.png'),
  path = 'results/images/red pebbles/',
  w = 8,
  h = 6,
  dpi = 225
)

ggplot() +
  geom_polygon(
    pools_shape,
    mapping  = aes(x = long, y = lat, group = group),
    colour = 'black',
    fill = NA,
    size = 0.7,
    alpha = 0.5
  ) +
  geom_point(data = df_for_vis,
             aes(
               x = x,
               y = y,
               shape = Event,
               col = `Xgb wasserstein`
             ),
             size = 5) +
  scale_color_viridis_c('Average\nWasserstein \ndistance',
                        alpha = 0.8,
                        limits = c(4, 42)) + theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16)
  ) +
  xlab('longitude') +
  ylab('latitude') + guides(colour = guide_colourbar(order = 1),
                            shape = guide_legend(order = 2))


ggsave(
  filename = paste0('Wasserstein_XGB', '.png'),
  path = 'results/images/red pebbles/',
  w = 8,
  h = 6,
  dpi = 225
)


ggplot() +
  geom_polygon(
    pools_shape,
    mapping  = aes(x = long, y = lat, group = group),
    colour = 'black',
    fill = NA,
    size = 0.7,
    alpha = 0.5
  ) +
  geom_point(
    data = df_for_vis,
    aes(
      x = x,
      y = y,
      shape = Event,
      col = `FK wasserstein` - l_$`Xgb wasserstein`
    ),
    size = 5
  ) +
  scale_color_viridis_c('Difference', alpha = 0.8) + theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16)
  ) +
  xlab('longitude') +
  ylab('latitude') + guides(colour = guide_colourbar(order = 1),
                            shape = guide_legend(order = 2))


ggsave(
  filename = paste0('Wasserstein_FK_-_XGB', '.png'),
  path = 'results/images/red pebbles/',
  w = 8,
  h = 6,
  dpi = 225
)

library(ggrepel)

ggplot(Location_Event, aes(x = x, y = y)) +
  theme_bw() + theme(
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16)
  ) +
  xlab('longitude') +
  ylab('latitude') +
  geom_point(
    aes(x = x, y = y, col = location),
    size = 8,
    alpha = 0.8,
    position = position_jitter(width = 0.00001, height = 0.00001)
  ) +
  scale_color_viridis_d() + theme(
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16)
  ) +
  geom_polygon(
    pools_shape,
    mapping  = aes(x = long, y = lat, group = group),
    colour = 'black',
    fill = NA,
    size = 0.7,
    alpha = 0.5
  ) +
  geom_label_repel(
    aes(label = location),
    point.padding = 1.5,
    fill = 'gray',
    size = 8,
    alpha = 1
  )


ggsave(
  filename = paste0('Locations_', '.png'),
  path = 'results/images/red pebbles/',
  w = 8,
  h = 6,
  dpi = 225
)


l_
l__ <- list()
i <- 2
for (i in 1:7) {
  n_binom <- 400
  pre__ <- transform_red_data(list_of_red_data[[i]][-1, 2:3])
  post__ <-
    transform_red_data(list_of_red_data[[i]][-1, 5:6]) %>% drop_na()
  
  df_for_xgb <- cbind(
    Location_Event_[rep(i, length(pre__$pebble_PC1)),] %>%
      select(x, y, duration, mean_h, mean_Q) %>%
      rename(X_Start = x, Y_Start = y),
    pre__ %>% select(a_axis, b_axis)
  ) %>% select(X_Start, Y_Start, a_axis, b_axis, mean_h, mean_Q, duration)
  
  xgb_pred_ <- predict(xgb_train, df_for_xgb %>% data.matrix())
  
  xgb_realisations_ <-
    sapply(xgb_pred_, function(p)
      rbinom(n_binom, 1, p))
  
  realizations_ <- sapply(1:n_binom, function(row_) {
    if (c(0) %in% unique(xgb_realisations_[row_,])) {
      obs_ <- pre__$diam_[xgb_realisations_[row_, ] == 0]
    }
  })
  
  realizations_ <-
    sapply(which(sapply(realizations_, length) >= 4)[1:100], function(id)
      realizations_[[id]])
  df_xgb_vis <-
    lapply(1:100, function(j)
      data.frame(diam = realizations_[[j]], id_ = paste0('realization', j))) %>% do.call('rbind', .)  %>%
    bind_rows(data.frame(diam = post__$diam_,
                         id_ = 'Actual'))
  
  
  
  df_xgb_vis$sim_ <- 0 == grepl('real', df_xgb_vis$id_, 'sim')
  
  ggplot(df_xgb_vis) + stat_ecdf(aes(
    x = diam,
    col = sim_,
    group = id_,
    size = sim_,
    alpha = sim_
  )) +
    theme_bw() +
    scale_color_manual('',
                       values = c('gray55', 'black'),
                       labels = c('', '')) +
    scale_alpha_manual('', values = c(0.3, 1)) +
    scale_size_manual('', values = c(1, 2)) +
    xlab('Diameter') +
    ylab('') +
    theme(
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.position = "none"
    ) +
    xlim(0, 100)
  
  
  ggsave(
    filename = paste0('xgb_ecdf_', Location_Event$location[i], '.png'),
    path = 'results/images/red pebbles/',
    w = 8,
    h = 6,
    dpi = 150
  )
  
  F_ <-
    predictFstat(
      g_UTrK_fit_,
      .newCoordinates = Location_Event[i, 5:6] %>% set_colnames(c('x', 'y')),
      .type = 'OK',
      .what = vName_
    )$Forecast[, 1]
  F_f_ <- approxfun(x = -50:50,
                    y = F_,
                    method = 'linear')
  
  pre__pc <- pre__$pebble_PC1
  pre__pc[pre__pc < -50] <- -50
  
  pre__pc[pre__pc >= 50] <- 50
  FK_pred_ <- F_f_(pre__pc)
  
  
  FK_realisations_ <-
    sapply(FK_pred_, function(p)
      rbinom(n_binom, 1, p))
  
  
  
  
  realizations_ <- sapply(1:n_binom, function(row_) {
    if (c(0) %in% unique(FK_realisations_[row_,])) {
      obs_ <- pre__$diam_[FK_realisations_[row_, ] == 0]
    }
  })
  
  realizations_ <-
    sapply(which(sapply(realizations_, length) >= 4)[1:100], function(id)
      realizations_[[id]])
  df_FK_vis <-
    lapply(1:100, function(j)
      data.frame(diam = realizations_[[j]], id_ = paste0('realization', j))) %>% do.call('rbind', .)  %>%
    bind_rows(data.frame(diam = post__$diam_,
                         id_ = 'Actual'))
  
  df_FK_vis$sim_ <- 0 == grepl('real', df_FK_vis$id_, 'sim')
  
  
  ggplot(df_FK_vis) + stat_ecdf(aes(
    x = diam,
    col = sim_,
    group = id_,
    size = sim_,
    alpha = sim_
  )) +
    theme_bw() +
    scale_color_manual('',
                       values = c('gray55', 'black'),
                       labels = c('', '')) +
    scale_alpha_manual('', values = c(0.3, 1)) +
    scale_size_manual('', values = c(1, 2)) +
    xlab('Diameter') +
    ylab('') +
    theme(
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.position = "none"
    ) +
    xlim(0, 100)
  
  ggsave(
    filename = paste0('FK_ecdf_', Location_Event$location[i], '.png'),
    path = 'results/images/red pebbles/',
    w = 8,
    h = 6,
    dpi = 150
  )
  
}

Location_Event_ %>% bind_cols(l_) %>% select(-x_WGS84,-y_WGS84,-n_o_s,-cl_)

openxlsx::write.xlsx(
  Location_Event_ %>% bind_cols(l_) %>% select(-n_o_s,-cl_,-weather_PC1,-weather_PC2),
  file = 'results/tables/red_pebbles_results.xlsx'
)



