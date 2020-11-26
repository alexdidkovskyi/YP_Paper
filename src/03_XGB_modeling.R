# !diagnostics off

#' Regression and Classification of Yellow Pebbles` Movement
#' Scalar Perspective: GLMNET and XGB


options(max.print = 100)

library(tidyverse)
library(magrittr)
library(dummies)
library(xgboost)
library(openxlsx)
library(sp)
library(pROC)
library(tictoc)


water_depth_threshold <-  35

source('src/AUX_01_models.R')
#source('src/results_delivery.R')
#source('src/models.R')
### classification error - AUC
### regression error - RMSE


list2env(readRDS('data/preprocessed/preprocessed_observations.RDS'),
         globalenv())


baseline_weather <-
  readRDS('data/preprocessed/baseline_weather_df.RDS')



data_list <- readRDS('data/preprocessed/all_datasets.RDS')


list2env(data_list$dfs, globalenv())
list2env(data_list$dfs_test, globalenv())
list2env(data_list$shape, globalenv())



df_$duration <- df_pred_$duration
df__ <-  df_
df_coordinates_and_type__ <-  df_coordinates_and_type_
df_pred__ <- df_pred_



clf_data_ <- df_ %>%
  left_join(
    df_coordinates_and_type_ %>% select(
      starts_with('TipoStart'),
      id_o,
      X_Start,
      Y_Start,
      X_Start_WGS84,
      Y_Start_WGS84
    )
  ) %>%
  arrange(id_o) %>%
  left_join(df_pred_) %>%
  select(-TipoStart,-id_o)
colnames_clf_ <-
  clf_data_ %>% select(starts_with('TipoStart')) %>% colnames() %>% gsub('TipoStart', '', .)

clf_data_$domain_type <-
  colnames_clf_[apply(clf_data_ %>% select(starts_with('TipoStart')) , 1, function(x)
    which(x == 1))]



# case_num <- 1
for (case_num in c(1, 1.5, 2, 3)) {
  ### Case 1 - All events, all features
  
  
  if (case_num == 1) {
    text_ <- 'All_events_all_features'
    
    clf_data_for_cv <-
      clf_data_ %>% select(
        -pebble_PC1,
        -pebble_PC2,-is_stuck_1,
        -is_stuck_2,
        -is_stuck_3,-n_o_s,
        -weather_PC1,
        -weather_PC2,-X_Start_WGS84,
        -Y_Start_WGS84,-graph_dist,
        -graph_vel,
        -domain_type,
        -cl_
      )
    #clf_data_for_cv
    y_ <- clf_data_$graph_dist > 0
    clf_data__ <- clf_data_ %>% select(domain_type)
    
    max_depth <- 7
  }
  ### Case 1.5 - typical events, just PCs + locations + subdomains
  
  if (case_num == 1.5) {
    text_ <- 'typical_events_PCs'
    
    
    clf_data_for_cv <-
      clf_data_ %>% filter(cl_ == 1) %>% select(
        pebble_PC1,
        pebble_PC2,
        weather_PC1,
        weather_PC2,
        starts_with('Tipo'),
        X_Start,
        Y_Start
      )
    #clf_data_for_cv
    y_ <- clf_data_[clf_data_$cl_ == 1,]$graph_dist > 0
    clf_data__ <-
      clf_data_ %>% filter(cl_ == 1) %>% select(domain_type)
    max_depth <- 6
    
  }
  
  ### Case 2 0 All events, PCs + subdomains
  
  if (case_num == 2) {
    text_ <- 'All_events_PCs_without_location'
    
    
    clf_data_for_cv <- clf_data_ %>% select(pebble_PC1,
                                            pebble_PC2,
                                            weather_PC1,
                                            weather_PC2,
                                            starts_with('Tipo'))
    #clf_data_for_cv
    y_ <- clf_data_$graph_dist > 0
    clf_data__ <- clf_data_ %>% select(domain_type)
    max_depth <- 6
  }
  
  
  ### Case 3 - typical events, Pebble PC1, locations, subdomains
  
  
  if (case_num == 3) {
    text_ <- 'typical_events_pebble_PC1'
    
    clf_data_for_cv <-
      clf_data_ %>% filter(cl_ == 1) %>% select(pebble_PC1, starts_with('Tipo'), X_Start, Y_Start)
    #clf_data_for_cv
    y_ <- clf_data_[clf_data_$cl_ == 1,]$graph_dist > 0
    
    clf_data__ <-
      clf_data_ %>% filter(cl_ == 1) %>% select(domain_type)
    
    max_depth <- 5
  }
  
  
  ### Estimation of the max_depth:
  #tic()
  #depth_results_df <-
  #  sapply(3:8, function(depth) {
  #
  #    a <- replicate(n = 7, expr = {
  #      xgb_cv <-
  #        xgb.cv(data = xgb.DMatrix(data = clf_data_for_cv %>% as.matrix(), label = y_),
  #             params = list(eta = 0.03,
  #                           objective = 'binary:logistic',
  #                           max_depth = depth),
  #             metrics = 'logloss',
  #             nrounds = 500,
  #             early_stopping_rounds = 25,
  #             print_every_n = 500,
  #             nfold =  5,
  #             prediction = T)
  #
  #      F1_thr <- get_F1_thr(xgb_cv$pred, act = y_)
  #      acc <- MLmetrics::Accuracy(xgb_cv$pred >= F1_thr, y_)
  #      })
  #    return(c(mean(a), sd(a)))
  #})%>% t() %>% data.frame() %>% set_colnames(.,c('mean','sd')) %>% mutate(max_depth = 3:8)
  #toc()
  
  set.seed(1234)
  folds_list <-
    replicate(caret::createFolds(1:nrow(clf_data_for_cv), k = 5),
              n = 7,
              simplify = F)
  folds <- folds_list[[1]]
  
  
  tic()
  xgb_cv_rep <-
    lapply(folds_list, function(folds) {
      xgb_cv <-
        xgb.cv(
          data = xgb.DMatrix(data = clf_data_for_cv %>% as.matrix(), label = y_),
          params = list(
            eta = 0.03,
            objective = 'binary:logistic',
            max_depth = max_depth
          ),
          metrics = 'logloss',
          nrounds = 500,
          early_stopping_rounds = 25,
          print_every_n = 500,
          folds = folds,
          prediction = T
        )
      
      F1_thr <- get_F1_thr(y_, xgb_cv$pred)
      
      
      
      # Youden's J threshold
      roc_obj_ <- roc(y_, xgb_cv$pred)
      
      thry_ <- coords(roc_obj_, "best", "threshold")
      J_thr <-  thry_$threshold
      
      
      #MLmetrics::Accuracy(xgb_cv$pred >= F1_thr, y_)
      res <- list()
      res$F1_conf <-
        caret::confusionMatrix(as.factor(xgb_cv$pred >= F1_thr), as.factor(y_))
      res$AUC <- MLmetrics::AUC(xgb_cv$pred, y_)
      res$J_conf <-
        caret::confusionMatrix(as.factor(xgb_cv$pred >= J_thr), as.factor(y_))
      res$F1_thr <- F1_thr
      res$J_thr <- J_thr
      
      res$missclassf_by_domain <- clf_data__ %>%
        select(domain_type) %>%
        mutate(
          y_true = y_,
          y_pred_F1 = xgb_cv$pred >= F1_thr,
          y_pred_J = xgb_cv$pred >= J_thr
        ) %>%
        group_by(domain_type) %>%
        summarise(
          n_ = n(),
          missclassf_F1 = sum(y_true == y_pred_F1),
          missclassf_J = sum(y_true == y_pred_J)
        ) %>%
        mutate(
          missclassf_F1 = n_ - missclassf_F1,
          missclassf_J = n_ - missclassf_J
        ) %>%
        mutate(missclassf_F1_prop = missclassf_F1 / n_,
               missclassf_J_prop = missclassf_J / n_)
      
      res$pred <- xgb_cv$pred
      
      return(res)
      
    })
  
  toc()
  
  ###
  #i <- 1
  
  
  if (case_num == 1.5) {
    pools_shape_df <-
      lapply(1:45, function(i)
        data.frame(pools_shape@polygons[[i]]@Polygons[[1]]@coords) %>%
          set_colnames(c('x', 'y')) %>%
          mutate(
            Domain = pools_shape$TIPO[i],
            id_pools = as.numeric(pools_shape$id[i]),
            id_iter = i
          )) %>%
      do.call('rbind', .)
    
    pred_q <-
      sapply(xgb_cv_rep, function(x)
        (x$pred > x$F1_thr) == as.numeric(clf_data_$graph_dist[clf_data_$cl_ == 1] >
                                            0)) %>% rowMeans()
    
    summarise_df <- clf_data_for_cv %>%
      select(X_Start, Y_Start, pebble_PC1) %>%
      mutate(pred_q = pred_q) %>%
      left_join(
        df_coordinates_and_type_ %>%
          select(TipoStart, X_Start, Y_Start, poly_id) %>% distinct()
      )
    
    summarise_df %>% group_by(poly_id) %>% summarise(mean_pred_q = mean(pred_q))
    
    
    
    temp_XGB <-
      pools_shape_df %>% left_join(
        summarise_df %>% group_by(poly_id) %>% summarise(mean_pred_q = mean(pred_q)),
        by = c('id_iter' = 'poly_id')
      )
    
    
    ggplot() +
      geom_polygon(
        data = temp_XGB,
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
      'results/images/XGB & FK CV results/xgb_Average_accuracy_by_polygon.png'
    ))
    ###
    
    
    mean(sapply(xgb_cv_rep, function(x)
      (x$pred > x$F1_thr) == y_) %>% rowMeans())
    sum(xgb_cv_res$pred_quality)
    xgb_cv_res <- clf_data_for_cv %>%
      mutate(pred_quality = sapply(xgb_cv_rep, function(x)
        (x$pred > x$F1_thr) == y_) %>% rowMeans()) %>%
      mutate(pred_quality = factor(
        pred_quality,
        levels = sort(unique(pred_quality)),
        labels = c(paste0(0:7, ' out of 7'))
      ))
    
    
    
    ggplot(xgb_cv_res, aes(x  = X_Start, y = Y_Start, col = pred_quality)) + geom_point()  +
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
    
    
    ggsave(filename = paste0(
      'results/images/XGB & FK CV results/xgb_prediction_correctness_',
      text_,
      '.png'
    ))
  }
  ###
  
  wb <- createWorkbook()
  
  
  table_prob_res_ <-
    sapply(c(mean, sd), function(f)
      sapply(c('AUC', 'F1_thr', 'J_thr'),
             function(t_)
               sapply(1:7,
                      function(id)
                        pluck(xgb_cv_rep, id, t_)) %>% f)) %>% data.frame() %>%
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
                        pluck(xgb_cv_rep,
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
               pluck(xgb_cv_rep,
                     id,
                     'missclassf_by_domain') %>%
               select(-domain_type) %>%
               as.matrix()) %>%
        simplify2array() %>%
        apply(., c(1, 2), f) %>%
        data.frame() %>%
        round(2) %>%
        mutate(domain_type = xgb_cv_rep[[1]]$missclassf_by_domain$domain_type)
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
                        pluck(xgb_cv_rep,
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
          pluck(xgb_cv_rep, id, thr_type, 'overall'))[1, ] %>%
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
  
  
  
  
  saveWorkbook(wb,
               file = paste0('results/tables/xgb_', text_, '.xlsx'),
               overwrite = T)
  
}
####
