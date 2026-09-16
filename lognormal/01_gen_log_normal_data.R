# Simulate datasets

library(readr)
library(tidyverse)

sim_lognormal <- function(n, mu, sigma2, eps = 0, delta =0){
  u <- runif(n)
  is_out <- u < eps
  tibble(
    index = seq_len(n),
    is_out = is_out,
    y = rlnorm(n, meanlog = mu + delta * is_out, sdlog = sqrt(sigma2))
  )
}

mu <- 0
sigma2 <- 0.25
n <- 1000
n_sim <- 50

#outlier weight
pct_w_list <- list(
  "01" = 0.01,
  "05" = 0.05,
  "10" = 0.10
)

delta_seq <- c(-4, -3, -2, -1, 1, 2, 3, 4)

# Simulate data
all_data <- map_dfr(names(pct_w_list), function(out_pct_label){
  
  w <- pct_w_list[[out_pct_label]]

  map_dfr(delta_seq, function(delta_id){

    map_dfr(1:n_sim, function(sim_id){

      sim_lognormal(
        n = n,
        mu = mu,
        sigma2 = sigma2,
        eps = w,
        delta = delta_id
      ) |>
        mutate(
          pct = out_pct_label,
          delta = delta_id,
          sim = sim_id,
          .before = 1
        )

    })
  })
})

write_rds(all_data, "outputs/02_all_data_sims.rds", compress = "gz")











