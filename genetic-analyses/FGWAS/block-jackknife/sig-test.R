trait_1_dir <- "./F-GWAS/jackknife/average_obs"
trait_2_dir <- "./F-GWAS/jackknife/fis1"

n_blocks <- 200

block_file <- function(dir, n) {
  file.path(dir, paste0("block_", n, "_corrs.txt"))
}

get_estimate <- function(file) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE)
  idx <- match("r_direct_avg_NTC", df$correlation)
  as.numeric(df$est[idx])
}

trait_1_blocks <- vapply(seq_len(n_blocks), function(n) block_file(trait_1_dir, n), character(1))
trait_2_blocks <- vapply(seq_len(n_blocks), function(n) block_file(trait_2_dir, n), character(1))

trait_1_est <- vapply(trait_1_blocks, get_estimate, numeric(1))
trait_2_est <- vapply(trait_2_blocks, get_estimate, numeric(1))

diff <- trait_1_est - trait_2_est
mean_diff <- mean(diff)

se_diff <- sqrt((n_blocks - 1) / n_blocks * sum((diff - mean_diff)^2))
z <- mean_diff / se_diff
p <- 2 * pnorm(abs(z), lower.tail = FALSE)

data.frame(
  diff_estimate = mean_diff,
  se_diff = se_diff,
  z = z,
  p = p
)