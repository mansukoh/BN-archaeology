# ============================================================================
# normalize()
# Normalize continuous variables with the Yeo-Johnson transformation
# (bestNormalize package) and report Shapiro-Wilk p-values before and after.
#
# Arguments
#   dat       : data frame
#   cols      : column indices (or names) of the continuous variables
#   n_shapiro : max sample size for shapiro.test() (subsampled if larger)
#
# Returns a list with
#   dat_trans : data frame with transformed columns (other columns unchanged)
#   summary   : variable, lambda, Shapiro p-value before / after
#   models    : fitted yeojohnson objects (one per variable)
# ============================================================================

normalize <- function(dat, cols, n_shapiro = 5000) {

  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    install.packages("bestNormalize")
  }
  library(bestNormalize)

  # 1) Fit a Yeo-Johnson transformation to each continuous variable
  bn_models <- lapply(dat[, cols], yeojohnson)

  # 2) Apply the fitted transformations
  dat_t <- dat
  dat_t[, cols] <- as.data.frame(mapply(function(m, x) predict(m, newdata = x),
                                        bn_models, dat[, cols], SIMPLIFY = FALSE))

  # 3) Shapiro-Wilk normality test p-values (before and after transformation)
  shap <- function(x, n = n_shapiro){
    x <- x[is.finite(x)]
    if (length(x) > n) x <- sample(x, n)
    if (length(x) < 3) return(NA_real_)
    shapiro.test(x)$p.value
  }
  p_before <- sapply(dat[, cols], shap)
  p_after  <- sapply(dat_t[, cols], shap)

  # 4) Estimated Yeo-Johnson lambda for each variable
  lambdas <- sapply(bn_models, `[[`, "lambda")

  summary_tbl <- data.frame(
    var = names(dat)[cols],
    lambda = as.numeric(lambdas),
    p_before = as.numeric(p_before),
    p_after  = as.numeric(p_after),
    row.names = NULL
  )

  list(
    dat_trans = dat_t,
    summary   = summary_tbl,
    models    = bn_models
  )
}
