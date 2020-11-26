#' Auxiliary functions for scalar classification using GLM and XGBoost
#' Including
#' a) Custom function for Elastic CV classification with alpha parameter selection using CV
#' b) Custom function for XGB cv
#' c) Custom functions for CV results delivery (as images and tables)
#' d) Custom functions for optimal threshold selection (MAX_ACC, F1, Youden's J criteria)



#' Scalar CV (XGB, ElasticNet)
scalar_cv <- function(data_,
                      y_,
                      data_test_ = NULL,
                      folds_ = NULL,
                      params_,
                      error_f_,
                      type_p_ = 'clf',
                      type_model_ = 'glm',
                      nrounds_xgb_ = 500,
                      ...) {
  pred_test_ <-
    foldid_ <-
    coef_glm_ <-
    alpha_opt_ <-
    lambda_ <-  coef_xgb_ <-  best_iteration_ <-  NULL
  
  if (type_model_ == 'glm') {
    if (is.null(folds_) == F) {
      foldid_ <- rep(0, length(y_))
      
      for (i_ in 1:length(folds_)) {
        foldid_[folds_[[i_]]] <- i_
      }
      
    }
    
    m_ <- apply(data_ %>% as.matrix(), 2, mean)
    sd_ <- apply(data_ %>% as.matrix(), 2, sd)
    
    data_train_ <- data_ %>% mutate_all(., scale) %>% as.matrix()
    
    err_ <-     sapply(seq(0, 1, 0.1), function(alpha_) {
      cv_glmnet <-
        cv.glmnet(
          data_train_,
          y_,
          family = params_$family,
          foldid = foldid_,
          keep = T,
          alpha = alpha_,
          nfolds = 7
        )
      id_ <- which(cv_glmnet$lambda == cv_glmnet$lambda.min)
      
      if (type_p_ == 'clf') {
        pred_ <- e1071::sigmoid(cv_glmnet$fit.preval[, id_])
        f1_thr_ <- get_F1_thr(y_clf_, pred_)
        output_ <- MLmetrics::Accuracy(pred_ > f1_thr_, y_clf_)
      }
      
      if (type_p_ == 'reg') {
        pred_ <- cv_glmnet$fit.preval[, id_]
        output_ <- MLmetrics::RMSE(pred_, y_)
      }
      
      return(output_)
    })
    
    
    alpha_opt_ <- (which.max(err_) - 1) * 0.1
    if (type_p_ == 'reg') {
      alpha_opt_ <- (which.min(err_) - 1) * 0.1
    }
    
    cv_glmnet <-
      cv.glmnet(
        data_train_,
        y_,
        family = params_$family,
        alpha = alpha_opt_,
        nfolds = 7,
        keep = T
      )
    id_ <- which(cv_glmnet$lambda == cv_glmnet$lambda.min)
    
    pred_ <- cv_glmnet$fit.preval[, id_]
    
    if (type_p_ == 'clf') {
      pred_ <- e1071::sigmoid(pred_)
    }
    
    
    if (is.null(data_test_) == F) {
      data_test_sc_ <- scale(data_test_, m_, sd_)
      pred_test_ <-
        predict(cv_glmnet, data_test_sc_, s = cv_glmnet$lambda.min)
    }
    
    coef_glm_ <-
      coef(cv_glmnet, s = cv_glmnet$lambda.min) %>% as.matrix() %>% data.frame()
    lambda_ <- cv_glmnet$lambda.min
  }
  
  if (type_model_ == 'xgb') {
    data_train_ <- xgb.DMatrix(data = data_ %>% as.matrix(), label = y_)
    
    cv_xgb <-
      xgb.cv(
        data = data_train_,
        nrounds = nrounds_xgb_,
        early_stopping_rounds = 25,
        params = params_$params,
        prediction = T,
        folds = folds_,
        verbose = 0,
        eval_metric = params_$eval_metric
      )
    pred_ <- cv_xgb$pred
    best_iteration_ <- cv_xgb$best_iteration
    
    model_xgb <-
      xgb.train(
        data = data_train_,
        label = y_,
        nrounds = best_iteration_,
        params = params_$params,
        verbose = 0,
        eval_metric = params_$eval_metric
      )
    coef_xgb_ <- xgb.importance(model = model_xgb)
    
    if (is.null(data_test_) == F) {
      pred_test_ <- predict(model_xgb, data_test_ %>% as.matrix())
    }
    
  }
  
  
  if (type_p_ == 'reg') {
    pred_[pred_ < 0] <- 0
    
    if (is.null(data_test_) == F) {
      pred_test_[pred_test_ < 0] <-  0
    }
  }
  
  
  
  return(
    list(
      error_ = error_f_(pred_, y_),
      pred_ = pred_,
      pred_test_ = pred_test_ %>% as.numeric(),
      coef_glm_ = coef_glm_,
      alpha_ = alpha_opt_,
      lambda_ = lambda_,
      coef_xgb_ = coef_xgb_,
      best_iteration_ = best_iteration_
    )
  )
  
}


#' CV results visualization
res_ggplot <- function(data, pools_shape_ = pools_shape) {
  ggplot(data, aes(x = x, y = y, col = pred__)) + theme_bw() +
    theme(
      axis.text = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 16)
    ) +
    geom_point(size = 2, alpha = 0.8)  +
    geom_polygon(
      pools_shape_,
      mapping  = aes(x = long, y = lat, group = group),
      colour = 'black',
      fill = NA,
      size = 0.7,
      alpha = 0.5
    ) +
    xlab('longitude') +
    ylab('latitude')
  
}


#' Get classification threshold according to F1-criteria
get_F1_thr <- function(act,
                       pred,
                       step = 0.01) {
  thr_vec <- sapply(seq(0.01, 0.99, step), function(thr) {
    pred__ <- as.numeric(pred >= thr)
    if (length(unique(pred__)) == 1)  {
      return(0)
    } else {
      return(MLmetrics::F1_Score(as.numeric(act), pred__))
    }
  })
  
  thr <- which.max(thr_vec) / 100
  return(thr)
  
}


get_ACC_thr <- function(act,
                        pred,
                        step = 0.01) {
  thr_vec <- sapply(seq(0.01, 0.99, step), function(thr) {
    pred__ <- as.numeric(pred >= thr)
    if (length(unique(pred__)) == 1)  {
      return(0)
    } else {
      return(MLmetrics::Accuracy(as.numeric(act), pred__))
    }
  })
  
  thr <- which.max(thr_vec) / 100
  return(thr)
}

get_tpnfpn <- function(act,
                       pred) {
  a <- case_when(
    act == 1 & pred == 1 ~ 'TP',
    act == 0 & pred == 1 ~ 'FP',
    act == 1 & pred == 0 ~ 'FN',
    act == 0 & pred == 0 ~ 'TN'
  )
  return(a)
}



#' Generate Test Data
generate_test_data <- function(df_test_peb_,
                               pool_grid,
                               data_colnames) {
  df__ <- merge(df_test_peb_, pool_grid) %>%
    rename(X_Start = x,
           Y_Start = y)
  
  if (setdiff(data_colnames, colnames(df__)) %>% length() > 0) {
    stop("Not enough variables in the test data")
  } else {
    df__  %<>%
      dplyr::select(data_colnames)
    return(df__)
  }
}





#' Results delivery. Visualization
visualize_results <- function(res_,
                              y_,
                              data_test_ = NULL,
                              df_coordinates_and_type_,
                              type_p_,
                              type_model_,
                              type_output = 'Probability',
                              folder_name_,
                              w_ = 7,
                              h_ = 4) {
  if (!dir.exists(paste0('results/results_', type_p_, '/', folder_name_)))
    dir.create(paste0('results/results_', type_p_, '/', folder_name_))
  
  ### Part 1. Test Prediction
  
  
  
  if (type_p_ %in% c('clf', 'func_clf')) {
    thr_ <- list()
    thr_[['F1']] <- get_F1_thr(y_, res_$pred_)
    
    roc_obj_ <- roc(y_, res_$pred_)
    
    thry_ <- coords(roc_obj_, "best", "threshold")
    thr_[['y']] <-  thry_$threshold
  }
  
  if (is.null(res_$pred_test_) == F &
      is.null(data_test_) == F) {
    repeat_num_ <-  nrow(data_test_) / 6
    
    test_data_for_pred_ <- data_test_ %>%
      mutate(
        size = rep(factor(
          c('big', 'small', 'big', 'small', 'big', 'small'),
          levels = c('small', 'big')
        ), repeat_num_),
        elongation = rep(factor(
          c(
            'high-high',
            'high-high',
            'high-medium',
            'high-medium',
            'medium-medium',
            'medium-medium'
          ),
          levels = c('medium-medium', 'high-medium', 'high-high')
        ), repeat_num_)
      ) %>%
      dplyr::select(X_Start, Y_Start, elongation, size) %>%
      rename(x = X_Start,
             y = Y_Start) %>%
      mutate(pred__ = res_$pred_test_)
    
    if (type_p_ %in% c('clf',  'func_clf')) {
      p <-
        res_ggplot(test_data_for_pred_) + facet_wrap(size ~ elongation) +
        scale_color_viridis_c('CV\nMovement\nProbability') + theme(axis.text = element_text(size = 8))
      
      ggsave(
        filename = paste0(
          'results/results_',
          type_p_,
          '/',
          folder_name_,
          "/test_CV_Probability_",
          type_model_,
          '.png'
        ),
        plot = p,
        width = w_,
        height = h_
      )
      
      for (thr_type_ in c('F1', 'y')) {
        thr__ <- thr_[[thr_type_]]
        
        test_data_for_pred_$pred__ <-
          factor(
            as.numeric(res_$pred_test_ >= thr__),
            levels = c(0, 1),
            labels = c('No', 'Yes')
          )
        
        p <-
          res_ggplot(test_data_for_pred_) + facet_wrap(size ~ elongation) +
          scale_color_viridis_d('Will it move?', guide = guide_legend(reverse = TRUE)) +
          theme(axis.text = element_text(size = 8))
        
        ggsave(
          filename = paste0(
            'results/results_',
            type_p_,
            '/',
            folder_name_,
            "/test_outcome_",
            thr_type_,
            '_',
            type_model_,
            '.png'
          ),
          plot = p,
          width = w_,
          height = h_
        )
        
      }
      
    }
    
    
    
    
    if (type_p_ %in% c('reg', 'func_reg')) {
      p <-
        res_ggplot(test_data_for_pred_) + facet_wrap(size ~ elongation) +
        scale_color_viridis_c(paste0('CV\n', type_output)) + theme(axis.text = element_text(size = 8))
      
      ggsave(
        filename = paste0(
          'results/results_',
          type_p_,
          '/',
          folder_name_,
          "/test_CV_",
          type_output,
          '_',
          type_model_,
          '.png'
        ),
        plot = p,
        width = w_,
        height = h_
      )
    }
  }
  
  
  ### 2. CV error: clf - FP/TP/FN/TN & T vs F, reg - abs
  
  
  df_res_and_coord <-
    df_coordinates_and_type_ %>%
    rename(x = X_Start, y = Y_Start) %>%
    mutate(act = y_,
           pred_ = res_$pred_) %>%
    dplyr::select(x, y, act, pred_)
  
  if (type_p_ %in% c('clf', 'func_clf')) {
    for (thr_type_ in c('F1', 'y')) {
      thr__ <- thr_[[thr_type_]]
      
      
      df_res_and_coord$pred___ <-  df_res_and_coord$pred_ >= thr__
      df_res_and_coord$pred__  <-
        get_tpnfpn(df_res_and_coord$act, df_res_and_coord$pred___)
      df_res_and_coord$pred__  <-
        factor(df_res_and_coord$pred__, levels = c('TN', 'FN', 'FP', 'TP'))
      
      p <-
        res_ggplot(df_res_and_coord) + scale_color_viridis_d('',
                                                             na.translate = FALSE,
                                                             guide = guide_legend(reverse = TRUE))
      
      ggsave(
        filename = paste0(
          'results/results_',
          type_p_,
          '/',
          folder_name_,
          "/Error_4cl_",
          thr_type_,
          '_',
          type_model_,
          '.png'
        ),
        plot = p,
        width = w_,
        height = h_
      )
      
      
      
      
      
      df_res_and_coord %<>% mutate(pred__ = factor(as.character(
        ifelse(act == pred___, "True", "False")
      ), levels = c("True", 'False')))
      
      p <-
        res_ggplot(df_res_and_coord) + scale_color_manual(
          'Correct\nPredictions',
          values = c('green', 'red'),
          na.translate = FALSE
        )
      
      ggsave(
        filename = paste0(
          'results/results_',
          type_p_,
          '/',
          folder_name_,
          "/Error_2cl_",
          thr_type_,
          '_',
          type_model_,
          '.png'
        ),
        plot = p,
        width = w_,
        height = h_
      )
      
      
      
    }
    
  } else {
    p <- res_ggplot(df_res_and_coord %>%
                      mutate(pred__ = pred_)) + scale_color_viridis_c('Prediction')
    
    ggsave(
      filename = paste0(
        'results/results_',
        type_p_,
        '/',
        folder_name_,
        "/Prediction_",
        type_model_,
        '.png'
      ),
      plot = p,
      width = w_,
      height = h_
    )
    
    
    
    df_res_and_coord %<>% mutate(pred__ = act - pred_)
    
    p <-
      res_ggplot(df_res_and_coord) + scale_color_viridis_c('Error')
    
    ggsave(
      filename = paste0(
        'results/results_',
        type_p_,
        '/',
        folder_name_,
        "/Error_",
        type_model_,
        '.png'
      ),
      plot = p,
      width = w_,
      height = h_
    )
    
    
  }
  
  
}



generate_test_data <- function(df_test_peb_,
                               pool_grid,
                               data_colnames) {
  df__ <- merge(df_test_peb_, pool_grid) %>%
    rename(X_Start = x,
           Y_Start = y)
  
  if (setdiff(data_colnames, colnames(df__)) %>% length() > 0) {
    stop("Not enough variables in the test data")
  } else {
    df__  %<>%
      dplyr::select(data_colnames)
    return(df__)
  }
}

#' Results delivery. Reporting
write_results <- function(res_,
                          y_,
                          type_p_,
                          type_model_,
                          folder_name_) {
  if (type_p_ %in% c('clf', 'func_clf')) {
    thr_ <- list()
    thr_[['F1']] <- get_F1_thr(y_, res_$pred_)
    
    roc_obj_ <- roc(y_, res_$pred_)
    
    thry_ <- coords(roc_obj_, "best", "threshold")
    thr_[['y']] <-  thry_$threshold
    
    thr_[['ACC']] <- get_ACC_thr(y_, res_$pred_)
    
    pred__ <-
      factor(
        as.numeric(res_$pred_ > thr_[['y']]),
        levels = c(1, 0),
        labels = c('Moved', 'Not Moved')
      )
    pred__F1_ <-
      factor(
        as.numeric(res_$pred_ > thr_[['F1']]),
        levels = c(1, 0),
        labels = c('Moved', 'Not Moved')
      )
    act__  <-
      factor(
        as.numeric(y_),
        levels = c(1, 0),
        labels = c('Moved', 'Not Moved')
      )
    conf_m <- caret::confusionMatrix(pred__, act__)
    conf_m_F1 <- caret::confusionMatrix(pred__F1_, act__)
    
    
    Acc <-
      sapply(thr_, function(x)
        MLmetrics::Accuracy(res_$pred_ > x, y_)) %>% data.frame()
    
    
    res_list_ <- list(
      auc = MLmetrics::AUC(res_$pred_, y_),
      accuracy = Acc,
      thresholds = data.frame(thr_),
      conf_m_F1 = conf_m_F1$table,
      conf_m_F1_overall = conf_m_F1$overall,
      conf_m_F1_byClass = conf_m_F1$byClass,
      conf_m = conf_m$table,
      conf_m_overall = conf_m$overall,
      conf_m_byClass = conf_m$byClass
    )
  }
  if (type_p_ == 'reg' |
      type_p_ == 'func_reg') {
    res_list_ <- list(RMSE = MLmetrics::RMSE(res_$pred_, y_),
                      MAE = MLmetrics::MAE(res_$pred_, y_))
    
  }
  
  if (type_model_ == 'glm') {
    res_list_[["coef_"]] <- res_$coef_glm_
    res_list_[['lambda_']] <- res_$lambda_
    res_list_[['alpha_']] <- res_$alpha_
  }
  
  if (type_model_ == 'xgb') {
    res_list_[["coef_"]] <- res_$coef_xgb_
    res_list_[['best_iteration']] <- res_$best_iteration_
  }
  
  wb <- createWorkbook()
  
  
  for (iter in 1:length(res_list_)) {
    temp_vec <- res_list_[[iter]]
    names(temp_vec) <- as.character(names(temp_vec))
    if (sum(is.na(names(temp_vec))) > 0) {
      names(temp_vec)[is.na(names(temp_vec))] <- 'NA_'
    }
    
    addWorksheet(wb, names(res_list_)[iter])
    if (class(temp_vec) == 'numeric')
      temp_vec <- data.frame(value = temp_vec)
    writeData(
      wb,
      names(res_list_)[iter],
      startCol = 3,
      startRow = 3,
      x = temp_vec,
      colNames = T,
      rowNames = T
    )
    
  }
  
  dir <- paste0("results/results_", type_p_, "/", folder_name_)
  
  if (!dir.exists(dir))
    dir.create(dir)
  saveWorkbook(
    wb,
    file = paste0(dir, "/cv_results_", type_model_, ".xlsx"),
    overwrite = TRUE
  )
  
}
