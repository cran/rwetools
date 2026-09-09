# tests/testthat/helper.R
# Shared test data loaded once for all test files

# Load sample data
sample_csv <- system.file("extdata", "sample_data.csv", package = "rwetools")
test_df  <- read.csv(sample_csv, stringsAsFactors = FALSE)
small_df <- test_df  # full dataset (2870 rows); subsetting caused sparse cell counts

# Variable types based on sample_data.csv schema
cont_vars   <- c("cont1", "cont2", "cont3", "cont4", "cont5", "cont6", "cont7")
binary_vars <- c("binary1", "binary2", "binary3", "binary4", "binary5",
                 "binary6", "binary7", "binary8", "binary9")
class_vars  <- c("cat1", "cat2", "cat3", "cat4")

# Deterministic, structurally valid 1:1 matched cohort for effect tests.
# Rows alternate reference/exposed within each pair.
make_test_matched_pairs <- function(n_pairs = 300L, seed = 1L) {
  set.seed(seed)
  ref_rows <- sample(which(small_df$exposure == 0), n_pairs)
  exp_rows <- sample(which(small_df$exposure == 1), n_pairs)
  stacked <- rbind(small_df[ref_rows, , drop = FALSE],
                   small_df[exp_rows, , drop = FALSE])
  order_idx <- as.vector(rbind(seq_len(n_pairs), n_pairs + seq_len(n_pairs)))
  out <- stacked[order_idx, , drop = FALSE]
  rownames(out) <- NULL
  out$pair_id <- rep(seq_len(n_pairs), each = 2L)
  out$.match_weights <- 1
  out
}
