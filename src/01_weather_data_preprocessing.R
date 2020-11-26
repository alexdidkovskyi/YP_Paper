# !diagnostics off


#' Initial weather processing.
#' Including
#' a) Events' flow characteristics recalculation
#' b) Clustering of the events

library(tidyverse)
library(magrittr)
library(pracma)
library(wavelets)
library(plotly)
list2env(readRDS('data/preprocessed/preprocessed_observations.RDS'),
         globalenv())


#'  h_CP_cm_ (water depth) equal to 35 is the threshold that indicates beginning of an event

water_depth_threshold <- 35

df_weather_ <-
  df_weather %>%
  group_by(EventoStart) %>%
  mutate(id = seq_along(Data),
         is_thr = h_CP_cm_ >= water_depth_threshold) %>%
  filter(is_thr == T) %>%
  mutate(id_delta = id - lag(id, default = NA)) %>%
  mutate(id_delta_diff = id_delta > 1) %>%
  mutate(number_of_subevents = sum(id_delta_diff, na.rm = T) + 1)


#' Average Intencity + Duration

#' Average water depth, average duration, number of subevents
baseline_weather <-
  df_weather_ %>%
  ungroup() %>%
  group_by(EventoStart) %>%
  summarise(
    duration = n(),
    mean_h = mean(h_CP_cm_),
    n_o_s = as.numeric(mean(number_of_subevents) > 1) + 1,
    mean_Q = mean(Q_CP_m3_s_)
  ) %>%
  mutate(id = 1:nrow(.))




#' Clustering of mobilizing events weather 

#' Scale to [0,1]
#' Cluster using hierarchical clustering with ward. linkage. Number of clusters = 3
library(scales)
hclt <-
  hclust(d = dist(
    baseline_weather %>% select(duration, mean_h, mean_Q) %>% mutate_all(., rescale)
  ),
  method = 'ward.D2')
cut_ <- cutree(hclt, 3)


png(
  filename =  paste0(
    'results/images/Mobilizing events/clustering_thr_',
    water_depth_threshold,
    '.png'
  )
)



hclust(d = dist(
  baseline_weather %>% select(duration, mean_h, mean_Q) %>% mutate_all(., rescale)
),
method = 'ward.D2') %>% plot(xlab = '')
dev.off()

png(
  filename =  paste0(
    'results/images/Mobilizing events/pairplot',
    water_depth_threshold,
    '.png'
  )
)
baseline_weather %>% select(duration, mean_h, mean_Q) %>% pairs(
  ,
  pch = 16,
  col = cutree(hclt, 3),
  cex = 2 ,
  cex.labels =  1.6,
  cex.axis = 1.4
)
dev.off()

pca_weather <-
  princomp(baseline_weather %>% select(duration, mean_h, mean_Q) %>% mutate_all(., rescale)) %>% summary()



#' Dimensionality reduction: First and Second PCs
baseline_weather <-
  baseline_weather %>% bind_cols(data.frame(pca_weather$scores[, 1:2]) %>% set_colnames(paste0('weather_PC', 1:2))) %>%
  mutate(cl_ = cut_)

saveRDS(
  baseline_weather %>% select(-id),
  'data/preprocessed/baseline_weather_df.RDS',
  compress = F
)


gc()



###
