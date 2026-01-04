library(reticulate)
setwd("~/research/items_llm/data/")
install_python()
# Import the pandas library
virtualenv_create("my-environment", version = "3.11.13")
use_virtualenv("my-environment")
reticulate::py_install("pandas", "my-environment")
pd <- import("pandas")

# Read the pickle file using the pandas function
pickle_data <- pd$read_pickle("tmp.pkl")


library(tidyverse)
items <- pickle_data %>% filter(is_subscale) %>% select(meta_doi, construct_name, item_text, meta_item_count, embeddings, keying)

items <- items %>% rowwise() %>% mutate(item_count = length(meta_item_count),
                                        mixed_keying = any(keying == "negative") & any(keying == "positive"))

unnested_items <- items %>%
  filter(mixed_keying) %>%
  # Ungroup to avoid rowwise complications and process columns as a whole
  ungroup() %>%

  # Step 1: Force every element in the key list-columns to be a character vector.
  # This is the most critical step to solve the error.
  mutate(
    # For each list element in 'keying', apply this function:
    # If it's NULL or has no length, return an NA of the correct type.
    # Otherwise, ensure it's a character vector.
    keying = map(keying, ~ if (is.null(.x) || length(.x) == 0) NA_character_ else as.character(.x)),
  ) %>%
  unnest(cols = c(construct_name, meta_item_count, item_text, keying, embeddings)) %>%
  select(meta_doi, construct_name, keying, item_text, everything())

test_items <- unnested_items %>% group_by(meta_doi, construct_name) %>% filter(n() > 6)

job_stress <- test_items %>% filter(construct_name == "Job Stress")
df <- job_stress %>% ungroup()

# Pull the matrix out, name the 768 columns, and bind back
emb_mat <- df[["embeddings"]]                    # 10 x 768 matrix
stopifnot(is.matrix(emb_mat), nrow(emb_mat) == nrow(df))

colnames(emb_mat) <- paste0("emb_", seq_len(ncol(emb_mat)))

df_wide <- bind_cols(
  df %>% select(-embeddings),
  as_tibble(emb_mat)
)

mu_pos <- df_wide %>% filter(keying == "positive") %>%
  summarise(across(starts_with("emb_"), mean))

mu_neg <- df_wide %>% filter(keying == "negative") %>%
  summarise(across(starts_with("emb_"), mean))

cor(unlist(mu_pos), unlist(mu_neg))

m <- (mu_pos + mu_neg) / 2
d_vec <- mu_pos - mu_neg
d <- d_vec / sqrt(sum(d_vec^2))

proj_score <- function(x) sum((x - m) * d)

df_scored <- df_wide %>%
  rowwise() %>%
  mutate(
    s_raw = proj_score(c_across(starts_with("emb_"))),
    s = if_else(keying == "negative", -s_raw, s_raw)
  ) %>%
  ungroup()
