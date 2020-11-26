options(max.print = 100)


#' Yellow Pebbles` Movement Classification
#' Functional perspective: Functional Kriging for probability curves generation

library(tictoc)
library(tidyverse)
library(magrittr)
library(dummies)
library(openxlsx)
library(sp)
library(pROC)
library(tictoc)
library(fdagstat)

source('src/AUX_01_models.R')
source('src/AUX_02_functional_models.R')
### classification error - AUC

water_depth_threshold <- 35

list2env(readRDS('data/preprocessed/preprocessed_observations.RDS'),
         globalenv())
baseline_weather <-
  readRDS('data/preprocessed/baseline_weather_df.RDS')
data_list <- readRDS('data/preprocessed/all_datasets.RDS')


#' Initial Part: Global model & flow line prediction --------------------
#' Initial parameters for curves generation and FK CV.
crvs_params_ <- list(
  type_ = 'domain',
  r_ = 5,
  range_ = 50,
  nr_obs_ = 12,
  lower_bnd_ = 0.01,
  # 0 + delta
  upper_bnd_ = 0.99,
  # 1 - delta
  bandwidth_ = 20,
  max_points_num = 30,
  range_l_ = -50,
  range_r_ = 50,
  range_step_ = 1,
  train_var = 'pebble_PC1',
  target_var = 'y_'
)

krig_params_ <- list(
  type_drift_ = 'OLS',
  intercept_ =  T,
  krig_type_ = 'OK',
  type_vgm_ = 'Bes',
  vName_ = 'Movement probability',
  regCoef_ = 1e1,
  Nlags_ = 10,
  ArgStep = 1,
  LagMax_ = 8,
  features_ = T
)


#' Specific subset of the data is used in the curves generation procedure if curves are based on the pebble_PC1 values

list2env(data_list$dfs, globalenv())
list2env(data_list$dfs_test, globalenv())
list2env(data_list$shape, globalenv())
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
        weather_PC1,
        weather_PC2
      )
  ) %>%
  mutate(y_ = as.numeric(graph_dist > 0))


targets_count_ <-
  targets_ %>% dplyr::select(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
  group_by(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
  summarise(count_ = n()) %>% ungroup() %>%
  mutate(id__ = 1:nrow(.))


set.seed(1234)

targets_$y_ <- as.numeric(targets_$graph_dist > 0)
y_ <- targets_$y_ > 0

FK_res_ <- list()

for (iter in 1:7) {
  folds_ <- caret::createFolds(1:nrow(targets_), 5)
  
  {
    res_ <- list()
    pred_test_ <- NULL
    
    pred_curves <-
      matrix(
        0,
        nrow = nrow(targets_),
        ncol = crvs_params_$range_r_ - crvs_params_$range_l_ + 1
      )
    
    
    tic()
    features_loc_ <-
      dummies::dummy(locations_$TipoStart) %>% as.data.frame()
    colnames(features_loc_) <-
      gsub(')+', '_', colnames(features_loc_))
    
    features_tar_ <-
      dummies::dummy(targets_$TipoStart) %>% as.data.frame()
    colnames(features_tar_) <-
      gsub(')+', '_', colnames(features_tar_))
    
    
    for (i in 1:length(folds_)) {
      pred_ <- list()
      
      test_id_ <- folds_[[i]]
      train_id_ <- setdiff(1:nrow(targets_), test_id_)
      
      targets_train_ <- targets_[train_id_,]
      targets_test_ <- targets_[test_id_,]
      
      
      locations_temp_ <- locations_ %>%
        dplyr::select(X_Start, Y_Start) %>%
        left_join(
          targets_train_ %>% dplyr::select(X_Start, Y_Start, poly_id) %>% distinct(),
          by = c("X_Start", 'Y_Start')
        )
      
      locations_train_id_ <- is.na(locations_temp_$poly_id) == F
      locations_train_ <- locations_[locations_train_id_, ]
      
      dist_matrix_train_ <-
        dist_matrix_[locations_train_id_, locations_train_id_]
      
      
      targets_count_train_ <-
        targets_train_ %>% dplyr::select(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
        group_by(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
        summarise(count_ = n()) %>% ungroup() %>%
        mutate(id__ = 1:nrow(.))
      
      
      
      df_func_train_ <-
        generate_curves(
          locations_ = locations_train_,
          targets_ =  targets_train_,
          dist_matrix_ = dist_matrix_train_,
          crvs_params_ = crvs_params_ ,
          targets_count_ = targets_count_train_
        )
      
      
      
      
      points_id <- df_func_train_$points_id
      df_func_train_ <- df_func_train_$df_func
      df_func_train_ <- log(df_func_train_ / (1 - df_func_train_))
      
      na_count_ <-
        apply(df_func_train_, 2, function(x)
          sum(is.na(x)))
      id_not_0 <- which(na_count_ == 0)
      df_func_train_red_ <- df_func_train_[, id_not_0]
      
      
      coord_train_ <-
        locations_train_[id_not_0, c('X_Start_WGS84', 'Y_Start_WGS84')]
      coord_train_ %<>% set_colnames(c('x', 'y'))
      
      coord_test_ <-
        targets_test_ %>% dplyr::select(X_Start_WGS84, Y_Start_WGS84) %>% set_colnames(c('x', 'y'))
      
      
      if (krig_params_$features_ == T) {
        id_col_ <- setdiff(1:ncol(features_loc_),
                           which(
                             colnames(features_loc_) %in% c("TipoStart_Steep_poolzone", 'TipoStart_Pools')
                           ))
        
        coord_train_ <- cbind(coord_train_,
                              features_loc_[locations_train_id_, ][id_not_0, id_col_])
        
        coord_test_ <- cbind(coord_test_,
                             features_tar_[test_id_, id_col_])
      }
      
      g_UTrK_fit_ <- with(
        krig_params_,
        functional_kriging_(
          coord_train_,
          df_func_train_red_,
          vName_ = vName_,
          type_drift_ = type_drift_,
          intercept_ = intercept_,
          LagMax_ = LagMax_,
          regCoef_ = regCoef_,
          type_vgm_ = type_vgm_,
          Nlags_ = Nlags_
        )
      )
      #print(plotVariogram(g_UTrK_fit_))
      
      
      pred_test_func_ <- with(
        krig_params_,
        predictFstat(
          g_UTrK_fit_,
          .newCoordinates = coord_test_,
          .what = vName_,
          .type = krig_type_,
          ErrorCheck = T,
          algIndependent = F
        )
      )
      pred_curves[test_id_,] <- pred_test_func_$Forecast %>% t()
      
      
      pred_ <- sapply(1:nrow(coord_test_), function(j) {
        res_ <- pred_test_func_$Forecast[, j]
        
        approxfun_ <- approxfun(
          seq(
            crvs_params_$range_l_,
            crvs_params_$range_r_,
            crvs_params_$range_step_
          ),
          1 / (1 + exp(1) ^ -res_)
        )
        
        pc1_ <-
          targets_test_ %>% dplyr::select(crvs_params_$train_var) %>% pull() %>% .[j]
        
        if (pc1_ <  crvs_params_$range_l_)
          pc1_ <-  crvs_params_$range_l_
        if (pc1_ >  crvs_params_$range_r_)
          pc1_ <-  crvs_params_$range_r_
        
        val <- approxfun_(pc1_)
        
        return(val)
        
      })
      
      
      
      
      rm(pred_test_func_)
      gc()
      
      res_[[i]] <- list(pred_ = pred_,
                        pred_test_ = pred_test_)
      
    }
    
    toc()
  }
  
  temp_res_ <- rep(0, nrow(targets_))
  for (i in 1:5) {
    temp_res_[folds_[[i]]] <- res_[[i]]$pred_
  }
  
  F1_thr <- get_F1_thr(y_, temp_res_)
  
  #' Youden's J threshold
  
  
  roc_obj_ <- roc(y_, temp_res_)
  
  thry_ <- coords(roc_obj_, "best", "threshold")
  J_thr <-  thry_$threshold
  
  
  #MLmetrics::Accuracy(temp_res_ >= F1_thr, y_)
  res <- list()
  res$F1_conf <-
    caret::confusionMatrix(as.factor(temp_res_  >= F1_thr), as.factor(y_))
  res$AUC <- MLmetrics::AUC(temp_res_ , y_)
  res$J_conf <-
    caret::confusionMatrix(as.factor(temp_res_  >= J_thr), as.factor(y_))
  res$F1_thr <- F1_thr
  res$J_thr <- J_thr
  
  res$missclassf_by_domain <- targets_ %>%
    rename(domain_type = TipoStart) %>%
    select(domain_type) %>%
    mutate(
      y_true = y_,
      y_pred_F1 = temp_res_ >= F1_thr,
      y_pred_J = temp_res_ >= J_thr
    ) %>%
    group_by(domain_type) %>%
    summarise(
      n_ = n(),
      missclassf_F1 = sum(y_true == y_pred_F1),
      missclassf_J = sum(y_true == y_pred_J)
    ) %>%
    mutate(missclassf_F1 = n_ - missclassf_F1,
           missclassf_J = n_ - missclassf_J) %>%
    mutate(missclassf_F1_prop = missclassf_F1 / n_,
           missclassf_J_prop = missclassf_J / n_)
  
  res$pred <- temp_res_
  res$curves <- pred_curves
  
  FK_res_[[iter]] <- res
  
  
}
###

#' Results visualization
i <- 1
pools_shape_df <-
  lapply(1:45, function(i)
    data.frame(pools_shape@polygons[[i]]@Polygons[[1]]@coords) %>% set_colnames(c('x', 'y')) %>%
      mutate(
        Domain = pools_shape$TIPO[i],
        id_pools = as.numeric(pools_shape$id[i]),
        id_iter = i
      )) %>% do.call('rbind', .)

pred_q <-
  sapply(FK_res_, function(x)
    (x$pred > x$F1_thr) == targets_$y_) %>% rowMeans()

summarise_df <- targets_ %>%
  select(X_Start, Y_Start, pebble_PC1) %>%
  mutate(pred_q = pred_q) %>%
  left_join(
    df_coordinates_and_type_ %>%
      select(TipoStart, X_Start, Y_Start, poly_id) %>% distinct()
  )

summarise_df %>% group_by(poly_id) %>% summarise(mean_pred_q = mean(pred_q))



temp_FK <-
  pools_shape_df %>% left_join(
    summarise_df %>% group_by(poly_id) %>% summarise(mean_pred_q = mean(pred_q)),
    by = c('id_iter' = 'poly_id')
  )


ggplot() +
  geom_polygon(
    data = temp_FK,
    aes(
      x = x,
      y = y,
      group = id_iter,
      fill = mean_pred_q
    ),
    color = 'black',
    size = 0.5
  ) +
  scale_fill_viridis_c('', na.value = 'gray') +
  theme_bw() +
  xlab('longitude') +
  ylab('latitude') +
  theme(text = element_text(size = 14))
ggsave(filename = paste0(
  'results/images/XGB & FK CV results/FK_Average_accuracy_by_polygon.png'
))


a__ <- temp_FK %>% rename(mean_FK = mean_pred_q) %>% left_join(temp_XGB %>% rename(mean_XGB = mean_pred_q))

a__ <- a__ %>% distinct() %>% mutate(temp_dif = mean_XGB - mean_FK)



ggplot() +
  geom_polygon(
    data = a__ %>% distinct(),
    aes(
      x = x,
      y = y,
      group = id_iter,
      fill = temp_dif
    ),
    color = 'black',
    size = 0.5
  ) +
  scale_fill_viridis_c('', na.value = 'gray') +
  theme_bw() +
  xlab('longitude') +
  ylab('latitude') +
  theme(text = element_text(size = 14))

###

FK_cv_res <- targets_ %>%
  mutate(pred_quality = sapply(FK_res_, function(x)
    (x$pred > x$F1_thr) == y_) %>% rowMeans()) %>%
  mutate(pred_quality = factor(
    pred_quality,
    levels = sort(unique(pred_quality)),
    labels = c(paste0(0:7, ' out of 7'))
  ))
FK_cv_res$pred_quality

sapply(FK_res_, function(x)
  (x$pred > x$F1_thr) == y_) %>% rowMeans() %>% mean()

ggplot(FK_cv_res, aes(x  = X_Start, y = Y_Start, col = pred_quality)) + geom_point()  +
  scale_color_viridis_d('Correctness') +
  theme_bw() +
  geom_polygon(
    pools_shape,
    mapping  = aes(x = long, y = lat, group = group),
    colour = 'black',
    fill = NA,
    size = 0.7,
    alpha = 0.5
  ) +
  xlab('longitude') +
  ylab('latitude') +
  theme(text = element_text(size = 14))


ggsave(filename = paste0('results/images/XGB & FK CV results/FK_prediction_correctness_J.png'))



####

#' Results delivery
wb <- createWorkbook()


table_prob_res_ <-
  sapply(c(mean, sd), function(f)
    sapply(c('AUC', 'F1_thr', 'J_thr'),
           function(t_)
             sapply(1:7,
                    function(id)
                      pluck(FK_res_, id, t_)) %>% f)) %>% data.frame() %>%
  set_colnames(c('mean', "sd")) %>%
  round(3)

table_prob_res_

addWorksheet(wb, 'AUC')
openxlsx::writeData(
  wb,
  'AUC',
  table_prob_res_,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)


tables_conf_m <-
  lapply(c(mean, sd), function(f) {
    lapply(c('F1_conf', 'J_conf'),
           function(thr_type) {
             lapply(1:7,
                    function(id)
                      pluck(FK_res_,
                            id,
                            thr_type,
                            'table'))  %>%
               simplify2array() %>%
               apply(., c(1, 2), f) %>% round(., 3)
           })
    
  })
tables_conf_m


addWorksheet(wb, 'conf_m_F1')
openxlsx::writeData(
  wb,
  'conf_m_F1',
  tables_conf_m[[1]][[1]] ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)

addWorksheet(wb, 'conf_m_J')
openxlsx::writeData(
  wb,
  'conf_m_J',
  tables_conf_m[[1]][[2]] ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)

addWorksheet(wb, 'sd_conf')
openxlsx::writeData(
  wb,
  'sd_conf',
  data.frame(sd_F1 = tables_conf_m[[2]][[1]][1, ],
             sd_J = tables_conf_m[[2]][[2]][1, ]) ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)





tables_missclassf_by_subdomain <-
  lapply(c(mean, sd), function(f) {
    lapply(1:7,
           function(id)
             pluck(FK_res_,
                   id,
                   'missclassf_by_domain') %>%
             select(-domain_type) %>%
             as.matrix()) %>%
      simplify2array() %>%
      apply(., c(1, 2), f) %>%
      data.frame() %>%
      round(2) %>%
      mutate(domain_type = FK_res_[[1]]$missclassf_by_domain$domain_type)
  })

tables_missclassf_by_subdomain


addWorksheet(wb, 'missclassf_mean')
openxlsx::writeData(
  wb,
  'missclassf_mean',
  tables_missclassf_by_subdomain[[1]]
  ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)

addWorksheet(wb, 'missclassf_sd')
openxlsx::writeData(
  wb,
  'missclassf_sd',
  tables_missclassf_by_subdomain[[2]] ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)





tables_res <-
  lapply(c('F1_conf', 'J_conf'),
         function(thr_type) {
           sapply(c(mean, sd), function(f) {
             sapply(1:7,
                    function(id)
                      pluck(FK_res_,
                            id,
                            thr_type,
                            'byClass'))   %>%
               apply(., c(1), f) %>% round(., 3)
           })
           
         })
tables_res

addWorksheet(wb, 'res_F1')
openxlsx::writeData(
  wb,
  'res_F1',
  tables_res[[1]] %>% set_colnames(c("mean", 'sd'))
  ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)

addWorksheet(wb, 'res_J')
openxlsx::writeData(
  wb,
  'res_J',
  tables_res[[2]] %>% set_colnames(c("mean", 'sd')) ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)





tables_ACC <-
  sapply(c(mean, sd), function(f) {
    sapply(c('F1_conf', 'J_conf'), function(thr_type) {
      sapply(1:7, function(id)
        pluck(FK_res_, id, thr_type, 'overall'))[1, ] %>%
        f
    }) %>% round(., 3)
  })

tables_ACC
addWorksheet(wb, 'ACC')
openxlsx::writeData(
  wb,
  'ACC',
  tables_ACC %>% set_colnames(c("mean", 'sd')) ,
  startRow = 2,
  startCol = 2,
  rowNames = T,
  colNames = T
)




saveWorkbook(
  wb,
  file = paste0('results/tables/FK_res_5_fold_OLS_10_rep_lag_7', '.xlsx'),
  overwrite = T
)





curves_ <-
  lapply(1:7, function(i)
    FK_res_[[i]]$curves) %>% simplify2array() %>% apply(., c(1, 2), mean)


df_curves_ <-
  data.frame(curves_) %>% set_colnames(-50:50) %>% mutate(subdomain = targets_$TipoStart, id = 1:1594) %>% gather(PC1, Probability,-subdomain, -id) %>%
  mutate(Probability = 1 / (1 + exp(1) ^ (-Probability)))






#' Curves visualization
# 50 Curves


df_func_train_ <- generate_curves(
  locations_ = locations_,
  targets_ =  targets_,
  dist_matrix_ = dist_matrix_,
  crvs_params_ = crvs_params_ ,
  targets_count_ = targets_count_
)




points_id <- df_func_train_$points_id
df_func_train_ <- df_func_train_$df_func
#df_func_train_ <- log(df_func_train_/(1 - df_func_train_))

na_count_ <- apply(df_func_train_, 2, function(x)
  sum(is.na(x)))
id_not_0 <- which(na_count_ == 0)
df_func_train_red_ <- df_func_train_[, id_not_0]





df_curves_ <-
  data.frame(t(df_func_train_)) %>% set_colnames(-50:50) %>% mutate(subdomain = locations_$TipoStart, id = 1:991) %>%
  gather(PC1, Probability,-subdomain, -id)

df_curves_$subdomain[df_curves_$subdomain == 'Bars_sedimentbuildupzones'] <-
  'Bars'
df_curves_$subdomain[df_curves_$subdomain == 'Steep_poolzone'] <-
  'pool-zone'
df_curves_$subdomain[df_curves_$subdomain == 'Run_rapid'] <-
  'Run-rapid'

df_curves_$subdomain %>% unique()

ggplot(df_curves_[df_curves_$id %in% sample(1:991, 50),],
       aes(
         x = as.numeric(PC1),
         y = Probability,
         group = id,
         col = subdomain
       )) + geom_line(size = 1) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  xlab("PC1") +
  ylab('Probability of Movement')#+
#scale_color_viridis_d()
ggsave(filename = paste0('results/images/curves generation/50_estimated_curves.png'))


df_curves_ %>% mutate(PC1 = as.numeric(PC1)) %>% group_by(subdomain, PC1) %>% summarise(Probability = mean(Probability, na.rm = T)) %>%
  ggplot(., aes(x = as.numeric(PC1), y = Probability, col = subdomain)) +
  geom_line(size = 1.1) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  xlab("PC1") +
  ylab('Probability of Movement')#+
# scale_color_viridis_d()
ggsave(filename = paste0('results/images/curves generation/mean_predicted_curves.png'))









coord_train_ <-
  locations_[id_not_0, c('X_Start_WGS84', 'Y_Start_WGS84')]





coord_train_ %<>% set_colnames(c('x', 'y'))

coord_test_ <-
  targets_test_ %>% dplyr::select(X_Start_WGS84, Y_Start_WGS84) %>% set_colnames(c('x', 'y'))


if (krig_params_$features_ == T) {
  id_col_ <- setdiff(1:ncol(features_loc_),
                     which(
                       colnames(features_loc_) %in% c("TipoStart_Steep_poolzone", 'TipoStart_Pools')
                     ))
  
  coord_train_ <- cbind(coord_train_,
                        features_loc_[, ][id_not_0, id_col_])
}
coord_train_ <- coord_train_[, 1:2]

krig_params_$LagMax_ <- 8
krig_params_$Nlags_ <- 10
krig_params_$type_drift_ <- 'OLS'


g_UTrK_fit_ <- with(
  krig_params_,
  functional_kriging_(
    coord_train_,
    df_func_train_red_,
    vName_ = vName_,
    type_drift_ = type_drift_,
    intercept_ = intercept_,
    LagMax_ = LagMax_,
    regCoef_ = regCoef_,
    type_vgm_ = type_vgm_,
    Nlags_ = Nlags_
  )
)




plotVariogram(g_UTrK_fit_)

ggsave(filename = paste0('results/images/FK variograms/VGM_LAGMAX_8_OLS_WO_features.png'))


krig_params_$type_drift_ <- 'GLS'
g_UTrK_fit_ <- with(
  krig_params_,
  functional_kriging_(
    coord_train_,
    df_func_train_red_,
    vName_ = vName_,
    type_drift_ = type_drift_,
    intercept_ = intercept_,
    LagMax_ = LagMax_,
    regCoef_ = regCoef_,
    type_vgm_ = type_vgm_,
    Nlags_ = Nlags_
  )
)

plotVariogram(g_UTrK_fit_)
ggsave(filename = paste0('results/images/FK variograms/VGM_LAGMAX_8_GLS.png'))
