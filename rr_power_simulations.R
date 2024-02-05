library(tidyverse)

number_of_items <- 175
number_of_scales <- 40
combinations_items <- choose(number_of_items, 2)
combinations_scales <- choose(number_of_scales, 2)

## Precision simulation for inter-item correlations
holdout <- arrow::read_feather("ignore.data-holdout-set-item-similarity-20230710-164559")
llm_holdout_meta <- arrow::read_feather("ignore.llmdata-holdout-set-item-similarity-20230710-164559.feather")
holdout_llm <- holdout %>%
  select(ItemStemIdA, ItemStemIdB, Pearson, CosineSimilarity) %>%
  left_join(llm_holdout_meta %>% select(ItemStemIdA = ItemStemId, VariableA = Variable, InstrumentA = instrument, ScaleA = scale_0, SubscaleA = scale_1)) %>%
  left_join(llm_holdout_meta %>% select(ItemStemIdB = ItemStemId, VariableB = Variable, InstrumentB = instrument, ScaleB = scale_0, SubscaleB = scale_1))

sim_results <- tibble()
library(lavaan)

for(i in 1:500) {
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



## precision analysis for reliability
cors_llm <- holdout_llm %>%
  select(x = VariableA, y = VariableB, r = CosineSimilarity) %>%
  as.data.frame() |>
  igraph::graph_from_data_frame(directed = FALSE) |>
  igraph::as_adjacency_matrix(attr = "r", sparse = FALSE)
diag(cors_llm) <- 1

cors_real <- holdout_llm %>%
  select(x = VariableA, y = VariableB, r = Pearson) %>%
  as.data.frame() |>
  igraph::graph_from_data_frame(directed = FALSE) |>
  igraph::as_adjacency_matrix(attr = "r", sparse = FALSE)
diag(cors_real) <- 1

subscales <- llm_holdout_meta %>%
  select(instrument, scale_0, scale_1, Variable, ItemStemId) %>%
  distinct() %>%
  mutate(instrument = coalesce(str_c(str_trim(instrument), "_"), ""),
         scale_0 = coalesce(str_c(str_trim(scale_0), "_"), ""),
         scale_1 = coalesce(str_trim(scale_1), ""),
         scale = str_replace_all(paste0(instrument, scale_0, scale_1), "[^a-zA-Z_0-9]", "_")
  )
n_distinct(subscales$scale)

scales <- llm_holdout_meta %>%
  select(instrument, scale_0, scale_1, Variable, ItemStemId) %>%
  distinct() %>%
  mutate(instrument = coalesce(str_c(str_trim(instrument), "_"), ""),
         scale_0 = coalesce(str_c(str_trim(scale_0), "_"), ""),
         scale_1 = "",
         scale = str_replace_all(paste0(instrument, scale_0, scale_1), "[^a-zA-Z_0-9]", "_")
  )

scales <- bind_rows(scales, subscales) %>% distinct() %>% filter(scale != "")
n_distinct(scales$scale)

scales <- scales %>%
  group_by(scale) %>%
  summarise(
    items = list(Variable),
    number_of_items = n_distinct(Variable),
    lvn = paste(first(scale), " =~ ", paste(Variable, collapse = " + "))) %>%
  drop_na()

random_scales <- list()
for(i in 1:200) {
  n_items <- rpois(1, 10)
  n_items <- if_else(n_items < 3, 3, n_items)
  random_scales[[i]] <- llm_holdout_meta %>%
    sample_n(n_items) %>%
    mutate(scale = paste0("random", i)) %>%
    group_by(scale) %>%
    summarise(
      items = list(Variable),
      number_of_items = n_distinct(Variable),
      lvn = paste(first(scale), " =~ ", paste(Variable, collapse = " + "))) %>%
    drop_na()
}

random_scales <- bind_rows(random_scales)
scales <- bind_rows(scales, random_scales)
n_distinct(scales$scale)

find_reverse_items_by_first_item <- function(rs) {
  # negatively correlated with first item
  items <- rs[-1, 1]
  reverse_keyed_items <- names(items)[which(items < 0)]
  reverse_keyed_items
}

reverse_items <- function(rs, reverse_keyed_items) {
  # Reverse the correlations for the reverse-keyed items
  for (item in reverse_keyed_items) {
    # Get the index of the reverse-keyed item
    item_index <- which(rownames(rs) == item)

    # Reverse the correlations
    rs[item_index, ] <- rs[item_index, ] * -1
    rs[, item_index] <- rs[, item_index] * -1

    # Since the diagonal is the correlation of the item with itself, set it back to 1
    rs[item_index, item_index] <- 1
  }
  rs
}

scales <- scales %>% filter(number_of_items >= 3)

scales <- scales %>%
  rowwise() %>%
  mutate(r_real = list(cors_real[items, items]),
         r_llm = list(cors_llm[items, items])) %>%
  mutate(reverse_items = list(find_reverse_items_by_first_item(r_real)),
         r_real_rev = list(reverse_items(r_real, reverse_items)),
         r_llm_rev = list(reverse_items(r_llm, reverse_items))) %>%
  mutate(
    rel_real = list(psych::alpha(r_real_rev, keys = F, n.obs = 400)$feldt)) %>%
  mutate(
    rel_llm = list(psych::alpha(r_llm_rev, keys = F, n.obs = 400)$feldt)) %>%
  mutate(rel_real_alpha = rel_real$alpha$raw_alpha,
         rel_real_upper = rel_real$upper.ci$raw_alpha,
         rel_real_lower = rel_real$lower.ci$raw_alpha,
         rel_real_se = (rel_real_upper - rel_real_lower)/2) %>%
  mutate(rel_llm_alpha = rel_llm$alpha$raw_alpha,
       rel_llm_upper = rel_llm$upper.ci$raw_alpha,
       rel_llm_lower = rel_llm$lower.ci$raw_alpha,
       rel_llm_se = (rel_llm_upper - rel_llm_lower)/2)

qplot(scales$rel_real_se)
qplot(scales$rel_real_alpha, scales$rel_real_se)
qplot(scales$number_of_items, scales$rel_real_se)


scales <- scales %>% mutate(
  alpha_se = mean(diff(unlist(psychometric::alpha.CI(rel_real_alpha, k = number_of_items, N = 400, level = 0.95))))
)

realistic_scales <- scales %>% filter(rel_real_alpha >= 0.3) %>% ungroup()

sim_results <- tibble()
for(i in 1:50) {
  picked_scales <- realistic_scales %>% filter(!str_detect(scale, "random")) %>% sample_n(40)
  subset <-
    bind_rows(picked_scales,
              realistic_scales %>% filter(str_detect(scale, "random")) %>% sample_n(100)
  )

  se2 <- mean(subset$alpha_se^2)

  r <- broom::tidy(cor.test(subset$rel_real_alpha, subset$rel_llm_alpha))
  (r$conf.high - r$conf.low)/2

  model <- paste0('
    # Latent variables
    PearsonLatent =~ 1*rel_real_alpha

    # Fixing error variances based on known standard errors
    rel_real_alpha ~~ ',se2,'*rel_real_alpha

    # Relationship between latent variables
    PearsonLatent ~~ rel_llm_alpha
  ')

  fit <- sem(model, data = subset)

  sim_results <- bind_rows(sim_results,
                           standardizedsolution(fit) %>% filter(lhs == "PearsonLatent", rhs ==  "rel_llm_alpha")
  )
}
sim_results %>% filter(is.finite(se)) %>%
  summarise(mean(est.std), sqrt(mean(se^2, na.rm = T)),
            max(se))


### TODO

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
