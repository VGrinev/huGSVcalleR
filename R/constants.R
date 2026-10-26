# Supported classical statistical tests for SNV calling.
# These tests operate on a frequency table (fr_data) and return a vector
# of p-values or probabilities.
SUPPORTED_CLASSICAL_TESTS <- c("binomial", "counting", "entropy",
                               "fisher", "poisson")

# Supported variant callers for SNV calling.
# Includes classical tests, BAM-based callers (fisherSmyth, gatk) and
# machine learning models (ml).
SUPPORTED_CALLERS <- c(SUPPORTED_CLASSICAL_TESTS,
                       "fisherSmyth", "gatk", "ml")
