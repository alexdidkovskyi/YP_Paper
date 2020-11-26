
#' Auxiliary functions for Functional Kringing
#' Including
#' a) Custom function for probability curves generation
#' b) Custom function for Functional Kriging
#' c) Custom function for FK CV
#' d) Custom function for CV of FK with features


#' Generate probability curves
#'
#' @param locations_ A dataframe with all locations
#' @param targets_ A dataframe with information about all movements
#' @param dist_matrix_ A square matrix of distances between every two points
#' @param targets_count_ A list of counts of observations associated with each location
#' @param crvs_params_ A list of parameters for curves generation
#' 
#' @return A list that consists of matrix of curves and list of points ids used for generation of each curve
generate_curves <-
  function(locations_,
           targets_,
           dist_matrix_,
           targets_count_,
           crvs_params_) {
    col_n <- intersect(colnames(locations_), colnames(targets_))
    data_func <- list()
    points_id <- rep(list(NA), nrow(locations_))
    
    for (j in 1:nrow(locations_)) {
      id_ <- which(dist_matrix_[j, ] <= crvs_params_$r_)
      
      
      if (length(id_) >  crvs_params_$max_points_num) {
        dist_ <-
          head(sort(dist_matrix_[j, id_]), crvs_params_$max_points_num)[crvs_params_$max_points_num]
        id_ <- which(dist_matrix_[j, ] <= dist_)
      }
      
      nr_ <-
        nrow(locations_[j, ] %>% left_join(targets_, by = c("X_Start", "Y_Start")))
      
      if (nr_ > 0) {
        p_l_ <- p_l <- locations_[id_, ] %>% left_join(targets_,  by = col_n)
        
        if (crvs_params_$type_  == 'poly') {
          p_l_ <- p_l %>% filter(poly_id == locations_$poly_id[j])
        }
        
        if (crvs_params_$type_  == 'domain') {
          p_l_ <- p_l %>% filter(TipoStart ==  locations_$TipoStart[j])
        }
        
        if (nrow(p_l_) >= crvs_params_$nr_obs_) {
          points_id[[j]] <-
            p_l_ %>% dplyr::select(X_Start, Y_Start) %>% left_join(targets_count_, by = c("X_Start", 'Y_Start')) %>% dplyr::select(id__) %>% pull() %>% unique()
          
          trn_var_  <-
            p_l_ %>% dplyr::select(crvs_params_$train_var) %>% pull()
          trgt_var_ <-
            p_l_ %>% dplyr::select(crvs_params_$target_var) %>% pull()
          
          ks_ <- with(
            crvs_params_,
            ksmooth(
              trn_var_,
              trgt_var_,
              kernel = 'normal',
              bandwidth = bandwidth_,
              range.x = c(range_l_, range_r_)
            )
          )
          
          f_appr <- approxfun(ks_$x, ks_$y)
          pred_ks <-
            with(crvs_params_, f_appr(seq(range_l_, range_r_, range_step_)))
          
          pred_ks[pred_ks < crvs_params_$lower_bnd_] <-
            crvs_params_$lower_bnd_
          pred_ks[pred_ks > crvs_params_$upper_bnd_] <-
            crvs_params_$upper_bnd_
          data_func[[j]] <- pred_ks
        } else {
          data_func[[j]] <-
            with(crvs_params_, rep(NA, (-range_l_ + range_r_) %/% range_step_ + 1))
        }
      } else {
        data_func[[j]] <-
          with(crvs_params_, rep(NA, (-range_l_ + range_r_) %/% range_step_ + 1))
      }
    }
    
    df_func <-
      data_func %>% do.call('rbind', .) %>% t() %>% data.frame()
    
    return(list(df_func = df_func, points_id = points_id))
    
  }


#' Funtional Kriging
#'
#' @param coord_ A data frame with coordinates of curves
#' @param func_ A matrix of curves
#' @param type_vgm_ Variogram type. one of vgm()
#' @param vName_ Variable Name
#' @param plot_ Variogram plot - T/F 
#' @param type_drift_ Drift type
#' @param ... FK parameters
#'
#' @return Returns the fstat
functional_kriging_ <- function(coord_,
                                func_,
                                type_vgm_ = 'Bes',
                                vName_,
                                Nlags_ = 40,
                                LagMax_ = 12,
                                ArgStep_ = 1,
                                regCoef_ = 1e1,
                                plot_ = F,
                                type_drift_ = 'OLS',
                                intercept_ = T) {
  if (type_drift_ != 'GLS') {
    g_UTrK_s <-
      fdagstat::fstat(
        NULL,
        vName = vName_,
        Coordinates = coord_,
        Functions = func_,
        scalar = F
      )
    g_UTrK_s <-
      fdagstat::estimateDrift(
        "~.",
        g_UTrK_s,
        .type = type_drift_,
        Intercept = intercept_ ,
        regCoef = regCoef_
      )
    g_UTrK_s  <-
      fdagstat::fvariogram(
        "~x+y",
        g_UTrK_s ,
        Nlags = Nlags_,
        LagMax = LagMax_,
        ArgStep = ArgStep_,
        useResidual = TRUE,
        comments = T
      )
    
    g_UTrK_s <-
      fdagstat::fitVariograms(g_UTrK_s, model = vgm(type_vgm_))
    if (g_UTrK_s$model$omni[[1]]$range[2]  < 0)
      (g_UTrK_s$model$omni[[1]]$range[2] <-
         abs(g_UTrK_s$model$omni[[1]]$range[2]))
    
    g_UTrK_s <- addCovariance(g_UTrK_s)
  } else {
    g_UTrK_s <-
      fdagstat::fstat(
        NULL,
        vName = vName_,
        Coordinates = coord_,
        Functions = func_,
        scalar = F
      )
    
    g_UTrK_s <-
      fdagstat::fvariogram(
        "~x+y",
        g_UTrK_s ,
        Nlags = Nlags_,
        LagMax = LagMax_,
        ArgStep = ArgStep_,
        useResidual = TRUE,
        comments = T
      )
    
    g_UTrK_s <-
      fdagstat::fitVariograms(g_UTrK_s, model = vgm(type_vgm_))
    if (g_UTrK_s$model$omni[[1]]$range[2]  < 0)
      (g_UTrK_s$model$omni[[1]]$range[2] <-
         abs(g_UTrK_s$model$omni[[1]]$range[2]))
    
    g_UTrK_s <- addCovariance(g_UTrK_s)
    g_UTrK_s <-
      fdagstat::estimateDrift(
        "~.",
        g_UTrK_s,
        .type = type_drift_,
        Intercept = intercept_ ,
        regCoef = regCoef_
      )
    g_UTrK_s  <-
      fdagstat::fvariogram(
        "~x+y",
        g_UTrK_s ,
        Nlags = Nlags_,
        LagMax = LagMax_,
        ArgStep = ArgStep_,
        useResidual = TRUE,
        comments = T
      )
    
    g_UTrK_s <-
      fdagstat::fitVariograms(g_UTrK_s, model = vgm(type_vgm_))
    if (g_UTrK_s$model$omni[[1]]$range[2]  < 0)
      (g_UTrK_s$model$omni[[1]]$range[2] <-
         abs(g_UTrK_s$model$omni[[1]]$range[2]))
    
  }
  return(g_UTrK_s)
}



#' Kriging CV
#'
#' @param locations_ A dataframe with all locations
#' @param targets_ A dataframe with information about all movements
#' @param dist_matrix_ A square matrix of distances between every two points
#' @param crvs_params_ A list of parameters for curves generation
#' @param krig_params_ A list of FK parameters
#' @param folds_ CV folds
#' @param test_data_  Test locations, if provided
#'
#' @return list of CV predictions and Test predictions
kriging_cv_ <- function(targets_,
                        locations_,
                        dist_matrix_,
                        crvs_params_,
                        krig_params_,
                        folds_,
                        test_data_) {
  res_ <- list()
  pred_test_ <- NULL
  
  if (is.null(test_data_) == F) {
    coord_test_data_ <- test_data_[, 1:2] %>% set_colnames(c('x', 'y'))
    pc1_test_data_ <-
      test_data_ %>% dplyr::select(crvs_params_$train_var) %>% pull()
  }
  
  for (i in 1:length(folds_)) {
    test_id_ <- folds_[[i]]
    train_id_ <- setdiff(1:nrow(targets_), test_id_)
    
    targets_train_ <- targets_[train_id_,]
    targets_test_ <- targets_[test_id_,]
    
    
    locations_temp_ <- locations_ %>%
      dplyr::select(X_Start, Y_Start) %>%
      left_join(targets_train_ %>% dplyr::select(X_Start, Y_Start, poly_id) %>% distinct())
    
    locations_train_id_ <- is.na(locations_temp_$poly_id) == F
    locations_train_ <- locations_[locations_train_id_, ]
    
    dist_matrix_train_ <-
      dist_matrix_[locations_train_id_, locations_train_id_]
    
    
    targets_count_train_ <-
      targets_train_ %>% dplyr::select(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
      group_by(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
      summarise(count_ = n()) %>% ungroup() %>%
      mutate(id__ = 1:nrow(.))
    
    
    
    df_func_train_ <- generate_curves(
      locations_ = locations_train_,
      targets_ =  targets_train_,
      dist_matrix_ = dist_matrix_train_,
      crvs_params_ = crvs_params_ ,
      targets_count_ = targets_count_train_
    )
    
    
    df_func_train_ <- df_func_train_$df_func
    
    
    
    na_count_ <- apply(df_func_train_, 2, function(x)
      sum(is.na(x)))
    
    id_0 <- which(na_count_ == 0)
    
    
    coord_train_ <-
      locations_train_[id_0, c('X_Start_WGS84', 'Y_Start_WGS84')]
    coord_train_ %<>% set_colnames(c('x', 'y'))
    df_func_train_red_ <- df_func_train_[, id_0]
    
    coord_test_ <-
      targets_test_ %>% dplyr::select(X_Start_WGS84, Y_Start_WGS84) %>% set_colnames(c('x', 'y'))
    
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
        type_vgm_ = type_vgm_
      )
    )
    
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
    
    
    pred_ <- sapply(1:nrow(coord_test_), function(j) {
      res_ <- pred_test_func_$Forecast[, j]
      
      res_[res_ < crvs_params_$lower_bnd_] <-
        crvs_params_$lower_bnd_
      res_[res_ > crvs_params_$upper_bnd_] <-
        crvs_params_$upper_bnd_
      approxfun_ <- approxfun(
        seq(
          crvs_params_$range_l_,
          crvs_params_$range_r_,
          crvs_params_$range_step_
        ),
        res_
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
    
    if (is.null(test_data_) == F) {
      pred_test_data_func_  <-
        predictFstat(
          g_UTrK_fit_,
          .newCoordinates = coord_test_data_,
          .what = crvs_params_$vName_,
          .type = "UK",
          ErrorCheck = T,
          algIndependent = T
        )
      
      ncol(pred_test_data_func_$Forecast)
      gc()
      
      pred_test_ <- sapply(1:nrow(coord_test_data_),
                           function(j) {
                             res_ <- pred_test_data_func_$Forecast[, j]
                             
                             
                             res_[res_ < crvs_params_$lower_bnd_] <-
                               crvs_params_$lower_bnd_
                             res_[res_ > crvs_params_$upper_bnd_] <-
                               crvs_params_$upper_bnd_
                             approxfun_ <-
                               approxfun(
                                 seq(
                                   crvs_params_$range_l_,
                                   crvs_params_$range_r_,
                                   crvs_params_$range_step_
                                 ),
                                 res_
                               )
                             
                             pc1_ <-
                               targets_test_ %>% dplyr::select(crvs_params_$train_var) %>% pull() %>% .[j]
                             
                             if (pc1_ < range_l_)
                               pc1_ <- range_l_
                             if (pc1_ > range_r_)
                               pc1_ <- range_r_
                             
                             val <- approxfun_(pc1_)
                             
                             return(val)
                             
                           })
      
      rm(pred_test_data_func_)
      gc()
      
    }
    
    res_[[i]] <- list(pred_ = pred_,
                      pred_test_ = pred_test_)
    
  }
  
  return(res_)
}






### -------------------------------------------------------------


#' Kriging CV with features
#'
#' @param locations_ A dataframe with all locations
#' @param targets_ A dataframe with information about all movements
#' @param dist_matrix_ A square matrix of distances between every two points
#' @param crvs_params_ A list of parameters for curves generation
#' @param krig_params_ A list of FK parameters
#' @param folds_ CV folds
#' @param test_data_  Test locations, if provided
#' @param exploration_mode If True all possible combinations of drift type, features, and prediction type are tested
#'
#' @return list of CV predictions and Test predictions
kriging_cv_with_features <- function(targets_,
                                     locations_,
                                     dist_matrix_,
                                     crvs_params_,
                                     krig_params_,
                                     folds_,
                                     thr_ = 0.3,
                                     features_ = F,
                                     error_ = F,
                                     exploration_mode = F) {
  res_ <- list()
  pred_test_ <- NULL
  
  
  
  features_loc_ <-
    dummies::dummy(locations_$TipoStart) %>% as.data.frame()
  colnames(features_loc_) <- gsub(')+', '_', colnames(features_loc_))
  
  features_tar_ <-
    dummies::dummy(targets_$TipoStart) %>% as.data.frame()
  colnames(features_tar_) <- gsub(')+', '_', colnames(features_tar_))
  
  
  for (i in 1:length(folds_)) {
    pred_ <- list()
    
    test_id_ <- folds_[[i]]
    train_id_ <- setdiff(1:nrow(targets_), test_id_)
    
    targets_train_ <- targets_[train_id_,]
    targets_test_ <- targets_[test_id_,]
    
    
    locations_temp_ <- locations_ %>%
      dplyr::select(X_Start, Y_Start) %>%
      left_join(targets_train_ %>% dplyr::select(X_Start, Y_Start, poly_id) %>% distinct())
    
    locations_train_id_ <- is.na(locations_temp_$poly_id) == F
    locations_train_ <- locations_[locations_train_id_, ]
    
    dist_matrix_train_ <-
      dist_matrix_[locations_train_id_, locations_train_id_]
    
    
    targets_count_train_ <-
      targets_train_ %>% dplyr::select(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
      group_by(X_Start, Y_Start, X_Start_WGS84, Y_Start_WGS84) %>%
      summarise(count_ = n()) %>% ungroup() %>%
      mutate(id__ = 1:nrow(.))
    
    
    
    df_func_train_ <- generate_curves(
      locations_ = locations_train_,
      targets_ =  targets_train_,
      dist_matrix_ = dist_matrix_train_,
      crvs_params_ = crvs_params_ ,
      targets_count_ = targets_count_train_
    )
    
    
    points_id <- df_func_train_$points_id
    df_func_train_ <- df_func_train_$df_func
    na_count_ <- apply(df_func_train_, 2, function(x)
      sum(is.na(x)))
    id_not_0 <- which(na_count_ == 0)
    
    
    if (exploration_mode == T) {
      for (error__ in c(T, F)) {
        df_func_train_red_ <- df_func_train_[, id_not_0]
        
        if (error__ == T) {
          vec_of_groups <-  locations_train_$TipoStart
          
          
          
          exch_id_ <-
            get_exchangeable_ids(
              vec_of_groups = vec_of_groups,
              points_id = points_id,
              targets_count_ = targets_count_train_,
              thr = thr_
            )
          
          exch_id_ <- intersect(exch_id_, id_not_0)
          exch_domains_ <-
            unique(locations_train_$TipoStart[exch_id_])
          
          
          means_by_domain_ <- list()
          
          for (domain_ in unique(locations_train_$TipoStart[id_not_0])) {
            domain_id_ <- which(locations_train_$TipoStart == domain_)
            
            if (domain_ %in% exch_domains_) {
              domain_id_ <- intersect(domain_id_, exch_id_)
            } else {
              domain_id_ <- intersect(domain_id_, id_not_0)
            }
            
            
            if (length(domain_id_) == 1) {
              means_by_domain_[[domain_]] <- df_func_train_[, domain_id_]
              
            } else {
              means_by_domain_[[domain_]] <-
                rowMeans(df_func_train_[, domain_id_])
              
            }
            
          }
          
          for (domain_ in unique(locations_train_$TipoStart[id_not_0])) {
            domain_mean_ <- means_by_domain_[[domain_]]
            
            id__ <-
              which(locations_train_$TipoStart[id_not_0] == domain_)
            
            df_func_train_red_[, id__] <-
              df_func_train_red_[, id__] - domain_mean_[[1]]
            
          }
        }
        
        
        
        
        for (features__ in c(T, F)) {
          coord_train_ <-
            locations_train_[id_not_0, c('X_Start_WGS84', 'Y_Start_WGS84')] %>% set_colnames(c('x', 'y'))
          
          coord_test_ <-
            targets_test_ %>% dplyr::select(X_Start_WGS84, Y_Start_WGS84) %>% set_colnames(c('x', 'y'))
          
          
          if (features__ == T) {
            id_col_ <- setdiff(1:ncol(features_loc_),
                               which(
                                 colnames(features_loc_) %in% c("TipoStart_Steep_poolzone", 'TipoStart_Pools')
                               ))
            
            id_col_ <- c(1:4, 6)
            coord_train_ <- cbind(coord_train_,
                                  features_loc_[locations_train_id_, ][id_not_0, id_col_])
            coord_test_ <- cbind(coord_test_,
                                 features_tar_[test_id_, id_col_])
          }
          
          for (type_drift__ in c('OLS', 'GLS')) {
            g_UTrK_fit_ <- with(
              krig_params_,
              functional_kriging_(
                coord_train_,
                df_func_train_red_,
                vName_ = vName_,
                type_drift_ = type_drift__,
                intercept_ = intercept_,
                LagMax_ = 8,
                regCoef_ = regCoef_,
                type_vgm_ = type_vgm_
              )
            )
            
            
            for (krig_type__ in c('OK', 'UK')) {
              pred_test_func_ <- with(
                krig_params_,
                predictFstat(
                  g_UTrK_fit_,
                  .newCoordinates = coord_test_,
                  .what = vName_,
                  .type = krig_type__,
                  ErrorCheck = T,
                  algIndependent = F
                )
              )
              
              
              pred__ <- sapply(1:nrow(coord_test_), function(j) {
                if (error__ == T) {
                  if (targets_test_$TipoStart[j] %in% names(means_by_domain_)) {
                    res_ <-
                      pred_test_func_$Forecast[, j] + means_by_domain_[[targets_test_$TipoStart[j]]]
                  } else {
                    res_ <-
                      pred_test_func_$Forecast[, j] + means_by_domain_ %>% simplify2array() %>% rowMeans()
                  }
                  
                } else {
                  res_ <- pred_test_func_$Forecast[, j]
                }
                
                res_[res_ < crvs_params_$lower_bnd_] <-
                  crvs_params_$lower_bnd_
                res_[res_ > crvs_params_$upper_bnd_] <-
                  crvs_params_$upper_bnd_
                approxfun_ <- approxfun(
                  seq(
                    crvs_params_$range_l_,
                    crvs_params_$range_r_,
                    crvs_params_$range_step_
                  ),
                  res_
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
              
              
              pred_[[paste(error__,
                           features__,
                           type_drift__,
                           krig_type__,
                           sep = '_')]] <-  pred__
            }
          }
          
          
          
          
        }
      }
    } else {
      if (error_ == T) {
        vec_of_groups <-  locations_train_$TipoStart
        
        
        
        exch_id_ <-
          get_exchangeable_ids(
            vec_of_groups = vec_of_groups,
            points_id = points_id,
            targets_count_ = targets_count_train_,
            thr = thr_
          )
        
        exch_id_ <- intersect(exch_id_, id_not_0)
        exch_domains_ <-
          unique(locations_train_$TipoStart[exch_id_])
        
        
        means_by_domain_ <- list()
        
        for (domain_ in unique(locations_train_$TipoStart[id_not_0])) {
          domain_id_ <- which(locations_train_$TipoStart == domain_)
          
          if (domain_ %in% exch_domains_) {
            domain_id_ <- intersect(domain_id_, exch_id_)
          } else {
            domain_id_ <- intersect(domain_id_, id_not_0)
          }
          means_by_domain_[[domain_]] <-
            rowMeans(df_func_train_[, domain_id_])
          
        }
        
        for (domain_ in unique(locations_train_$TipoStart[id_not_0])) {
          domain_mean_ <- means_by_domain_[[domain_]]
          
          id__ <-
            which(locations_train_$TipoStart[id_not_0] == domain_)
          
          df_func_train_red_[, id__] <-
            df_func_train_red_[, id__] - domain_mean_[[1]]
          
        }
      }
      
      coord_train_ <-
        locations_train_[id_not_0, c('X_Start_WGS84', 'Y_Start_WGS84')]
      coord_train_ %<>% set_colnames(c('x', 'y'))
      
      coord_test_ <-
        targets_test_ %>% dplyr::select(X_Start_WGS84, Y_Start_WGS84) %>% set_colnames(c('x', 'y'))
      
      
      if (features_ == T) {
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
          LagMax_ = 8,
          regCoef_ = regCoef_,
          type_vgm_ = type_vgm_
        )
      )
      
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
      
      
      pred_ <- sapply(1:nrow(coord_test_), function(j) {
        if (error_ == T) {
          if (targets_test_$TipoStart[j] %in% names(means_by_domain_)) {
            res_ <-
              pred_test_func_$Forecast[, j] + means_by_domain_[[targets_test_$TipoStart[j]]]
          } else {
            res_ <-
              pred_test_func_$Forecast[, j] + means_by_domain_ %>% simplify2array() %>% rowMeans()
          }
          
        } else {
          res_ <- pred_test_func_$Forecast[, j]
        }
        
        res_[res_ < crvs_params_$lower_bnd_] <-
          crvs_params_$lower_bnd_
        res_[res_ > crvs_params_$upper_bnd_] <-
          crvs_params_$upper_bnd_
        approxfun_ <- approxfun(
          seq(
            crvs_params_$range_l_,
            crvs_params_$range_r_,
            crvs_params_$range_step_
          ),
          res_
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
      
      
    }
    
    rm(pred_test_func_)
    gc()
    
    res_[[i]] <- list(pred_ = pred_,
                      pred_test_ = pred_test_)
    
  }
  
  return(res_)
}
