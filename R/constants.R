# Supported classical statistical tests for SNV calling.
SUPPORTED_CLASSICAL_TESTS <- c("binomial", "counting", "entropy",
                               "fisher", "poisson")

# All supported criteria for callSNVs().
SUPPORTED_CALLERS <- c(SUPPORTED_CLASSICAL_TESTS, "fisherRSubread", "gatk", "ml")

# Criteria that require fr_data (frequency table).
SUPPORTED_TABLE_CALLERS <- c(SUPPORTED_CLASSICAL_TESTS, "ml")

# Criteria that require bamFile.
SUPPORTED_BAM_CALLERS <- c("fisherRSubread", "gatk")
