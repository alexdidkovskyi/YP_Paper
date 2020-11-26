
#' Additional images
options(max.print = 100)

library(tictoc)
library(tidyverse)
library(magrittr)
library(dummies)
library(openxlsx)
library(sp)
library(pROC)
library(tictoc)
library(fdagstat)

### classification error - AUC

water_depth_threshold <- 35

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


baseline_weather %>% filter(cl_ == 1) %>% View()


p1 <-
  ggplot(baseline_weather, aes(x = duration, y = mean_Q, col = as.factor(cl_))) +
  geom_point(size = 1.5) + xlab('Duration') + ylab('Mean water flow') +
  scale_color_manual("", values = c('black', 'red', 'green')) + theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "none")



p2 <-
  ggplot(baseline_weather, aes(x = duration, y = mean_h, col = as.factor(cl_))) +
  geom_point(size = 1.5) + xlab('Duration') + ylab('Mean water depth') +
  scale_color_manual("", values = c('black', 'red', 'green')) + theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "none")


p3 <-
  ggplot(baseline_weather, aes(x = mean_h, y = mean_Q, col = as.factor(cl_))) +
  geom_point(size = 1.5) + xlab('Mean water depth') + ylab('Mean water flow') +
  scale_color_manual("", values = c('black', 'red', 'green')) + theme_bw() +
  theme(text = element_text(size = 12),
        legend.position = "none")
p4 <-
  egg::ggarrange(
    nrow = 2,
    ncol = 2,
    plots = list(p1, p2, p3),
    byrow = F
  )
ggsave(p4, filename = 'results/images/Mobilizing events/events_clustering.png')


####




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



#ggsave( filename = 'results/additional_pics/curves_generation_1.png')



locations_ %>%
  select(X_Start, Y_Start) %>%
  mutate(dist = dist_matrix_[500,]) %>%
  left_join(targets_) %>%
  filter(dist <= 5) %>%
  arrange(dist) %>% .[1:30, ] -> observations_for_plot


observations_for_plot$t_ <- 1


targets_ %>% left_join(observations_for_plot)
ggplot(
  targets_ %>%
    left_join(observations_for_plot) %>%
    mutate(y__ = case_when(is.na(t_) == T ~ -1,
                           is.na(t_) == F ~ y_)),
  aes(
    x = X_Start,
    y = Y_Start,
    col = as.character(y__),
    alpha = as.character(y__)
  )
) +
  geom_point(size = 1) +
  scale_color_manual(
    '',
    values = c('gray20', 'purple', 'blue'),
    labels = c(-1, 0, 1),
    guide = guide_legend(reverse = TRUE)
  ) +
  theme_bw() +
  scale_alpha_manual('',
                     values = c(0.8, 1, 1)) +
  geom_polygon(
    pools_shape,
    mapping  = aes(x = long,
                   y = lat,
                   group = group),
    colour = 'black',
    fill = NA,
    size = 0.7,
    alpha = 0.5
  ) +
  xlab('longitude') +
  ylab('latitude') +
  theme(text = element_text(size = 14),
        legend.position = 'none') +
  geom_point(
    aes(x = locations_$X_Start[500], y  = locations_$Y_Start[500]),
    col = "blue",
    size = 1
  ) +
  geom_point(
    aes(x = locations_$X_Start[500], y  = locations_$Y_Start[500]),
    col = "blue",
    size = 7,
    color = 'black',
    shape  = 1,
    stroke = 2
  )

ggsave(filename = 'results/images/curves generation/curves_generation_1.png')

curves_for_plot <-
  ksmooth(
    observations_for_plot$pebble_PC1,
    observations_for_plot$graph_dist > 0,
    kernel = 'normal',
    bandwidth = 20,
    range.x = c(-50, 50)
  ) %>% data.frame()

ggplot(observations_for_plot, aes(x = pebble_PC1, y = as.numeric(graph_dist > 0))) +
  geom_jitter(size = 2,
              width = 0.2,
              height = 0) + theme_bw() +
  xlab('PC1') +
  ylab('') +
  geom_line(data = curves_for_plot,
            aes(x = x, y = y),
            col = 'black',
            size = 2)

ggsave(filename = 'results/images/curves generation/curves_generation_2.png')
