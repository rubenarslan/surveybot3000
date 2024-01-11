library(tidyverse)

number_of_items <- 175
number_of_scales <- 40
combinations_items <- choose(number_of_items, 2)
combinations_scales <- choose(number_of_scales, 2)

ses <- data.frame(
  r = rep(seq(-1, 1, by = 0.01), each = 14),
  N = rep(c(50, 60, number_of_scales, combinations_items, combinations_scales, 70, 80, 100, 200, 300, 500, 1000, 2000, 4000), times = 201)) %>%
  mutate(
    se = (1 - r^2)/sqrt(N - 2)
    # https://www.jstor.org/stable/2277400?seq=1#page_scan_tab_contents
  )

ses %>%
  ggplot(aes(r, se, color = N)) +
  geom_point() +
  scale_y_sqrt(breaks = c(0.005, 0.002, seq(0, 1, by = 0.01)))

ses %>%
  filter(N %in% c(number_of_scales, combinations_items, combinations_scales)) %>%
  filter(round(r,2) %in% c(0.73, 0.80, 0.89)) %>%
  ggplot(aes(r, se, color = N)) +
  geom_point() +
  scale_y_sqrt(breaks = c(0.005, 0.002, seq(0, 1, by = 0.01)))

ses %>% filter(N == 300, between(r, -0.7, 0.7)) %>%
  summarise(mean(se))
ses %>% filter(N == 60, between(r, -0.5, 0.5)) %>%
  summarise(mean(se))
# SE lm pearson on cosinesimilarity 0.1154

holdout <- arrow::read_feather("ignore.data-holdout-set-item-similarity-20230710-164559")

items <- holdout %>% select(ItemStemIdA) %>% distinct() %>% sample_n(175) %>% pull(ItemStemIdA)

subset <- holdout %>% filter(ItemStemIdA %in% items, ItemStemIdB %in% items)

r <- broom::tidy(cor.test(subset$Pearson, subset$CosineSimilarity))
(r$conf.high - r$conf.low)/2



## having the LLM estimates is like having a sample of X for the items/rels/scale scores
testset <- arrow::read_feather("ignore.data-test-set-item-similarity-20230710-164559")
# testset <- arrow::read_feather("ignore.item-similarity-20230710-164559.scores.data.training.test.feather")

# SE items residual 0.1103
summary(lm(Pearson ~CosineSimilarity, testset))
ggplot(testset, aes(CosineSimilarity, Pearson)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm")

ses <- data.frame(
  r = rep(seq(-1, 1, by = 0.01), each = 14),
  N = rep(c(50, 60, number_of_scales, combinations_items, combinations_scales, 70, 80, 100, 200, 300, 500, 1000, 2000, 4000), times = 201)) %>%
  mutate(
    se = (1 - r^2)/sqrt(N - 2)
    # https://www.jstor.org/stable/2277400?seq=1#page_scan_tab_contents
  )

manifest_scores_ft <- feather::read_feather("holdout_manifest_scores_ft.feather")
ggplot(manifest_scores_ft, aes(machine_cor, human_cor)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm")
latent_scores2 <- feather::read_feather("holdout_scale_scores_ft.feather")
ggplot(latent_scores2, aes(machine_cor, human_latent_cor)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm")
summary(lm(human_latent_cor ~ machine_cor, latent_scores2))


ggplot(latent_scores2, aes(machine_cor, human_cor)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm")


ses %>% filter(between(r, -.6, .6), round(se, 2) == 0.11) %>% summarise(max(N), mean(N))
ses %>% filter(between(r, -.9, .9), round(se, 2) == 0.17) %>% summarise(max(N), mean(N))
ses %>% filter(between(r, 0, 1), round(se, 2) == 0.05) %>% summarise(max(N), mean(N))
ses %>% filter(between(r, 0.5, 1), round(se, 2) == 0.05) %>% summarise(max(N), mean(N))

