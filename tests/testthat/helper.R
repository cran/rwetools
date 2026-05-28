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
