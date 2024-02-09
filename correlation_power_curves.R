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

sim_results <- tibble()
library(lavaan)

for(i in 1:50) {
  items <- holdout %>% select(ItemStemIdA) %>% distinct() %>% sample_n(175) %>% pull(ItemStemIdA)

  subset <- holdout %>% filter(ItemStemIdA %in% items, ItemStemIdB %in% items)

  N <- 400
  subset <- subset %>% mutate(se = (1 - Pearson^2)/sqrt(N - 2))
  se2 <- mean(subset$se^2)

  r <- broom::tidy(cor.test(subset$Pearson, subset$CosineSimilarity))
  (r$conf.high - r$conf.low)/2

  model <- paste0('
    # Latent variables
    PearsonLatent =~ 1*Pearson

    # Fixing error variances based on known standard errors
    Pearson ~~ ',se2,'*Pearson

    # Relationship between latent variables
    PearsonLatent ~~ CosineSimilarity
  ')

  fit <- sem(model, data = subset)

  sim_results <- bind_rows(sim_results,
    standardizedsolution(fit) %>% filter(lhs == "PearsonLatent", rhs ==  "CosineSimilarity")
  )
}
sim_results %>% summarise(mean(est.std), sqrt(mean(se^2)))

## having the LLM estimates is like having a sample of X for the items/rels/scale scores
testset <- arrow::read_feather("ignore.data-test-set-item-similarity-20230710-164559")
# testset <- arrow::read_feather("ignore.item-similarity-20230710-164559.scores.data.training.test.feather")

# SE items residual 0.1103
summary(lm(Pearson ~CosineSimilarity, testset))

## Realistic SEs for a given N given the variability in empirical item correlations
r <- testset$Pearson
N <- rep(c(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), times = length(r))
se = (1 - r^2)/sqrt(N - 2)
tibble(N, se) %>% group_by(N) %>%
  summarise(sqrt(mean(se^2)))

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

r <- latent_scores2$human_latent_cor
N <- rep(c(20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), times = length(r))
se = (1 - r^2)/sqrt(N - 2)
tibble(N, se) %>% group_by(N) %>%
  summarise(sqrt(mean(se^2)))

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



### Precision sym synthetic nomological nets
library(tidyverse)
holdout <- arrow::read_feather("ignore.data-holdout-set-item-similarity-20230710-164559")
llm_holdout_meta <- arrow::read_feather("ignore.llmdata-holdout-set-item-similarity-20230710-164559.feather")
holdout_real <- holdout %>%
  select(ItemStemIdA, ItemStemIdB, Pearson) %>%
  left_join(llm_holdout_meta %>% select(ItemStemIdA = ItemStemId, VariableA = Variable, InstrumentA = instrument, ScaleA = scale_0, SubscaleA = scale_1)) %>%
  left_join(llm_holdout_meta %>% select(ItemStemIdB = ItemStemId, VariableB = Variable, InstrumentB = instrument, ScaleB = scale_0, SubscaleB = scale_1))

ho <- holdout_real %>%
  group_by(InstrumentA, ScaleA, SubscaleA) %>%
  mutate(ScaleA_item_number = n_distinct(ItemStemIdA)) %>%
  group_by(InstrumentB, ScaleB, SubscaleB) %>%
  mutate(ScaleB_item_number = n_distinct(ItemStemIdB))

subs <- ho %>% filter(ScaleA_item_number == 5, ScaleB_item_number == 5, SubscaleA == "Motor Impulsivity", ScaleB == "Satisfaction with Life")

scales <-
  bind_rows(
    subs %>%
      ungroup() %>%
      select(Variable = VariableA, Instrument = InstrumentA, Scale = ScaleA, Subscale = SubscaleA) %>%
      distinct(),
    subs %>%
      ungroup() %>%
      select(Variable = VariableB, Instrument = InstrumentB, Scale = ScaleB, Subscale = SubscaleB) %>%
      distinct() ) %>%
  mutate(
    scale_0 = coalesce(str_c(Scale, Subscale), Scale)) %>%
  # drop_na(Instrument, scale_0) %>%
  mutate(scale = str_replace_all(paste(Instrument, scale_0), "[^a-zA-Z_0-9]", "_")) %>%
  group_by(scale) %>%
  summarise(
    items = list(Variable),
    lvn = paste(first(scale), " =~ ", paste(Variable, collapse = " + "))) %>%
  drop_na()
