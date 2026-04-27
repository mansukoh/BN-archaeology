
# Normalize continuous variables using "bestNormalize" package


normalize <- function(dat, cols, n_shapiro = 5000) {
  
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    install.packages("bestNormalize")
  }
  library(bestNormalize)
  
  
  bn_models <- lapply(dat[, cols], yeojohnson)
  
  dat_t <- dat
  dat_t[, cols] <- as.data.frame(mapply(function(m, x) predict(m, newdata = x),
                                        bn_models, dat[, cols], SIMPLIFY = FALSE))
  
  # 3) Shapiro p-value
  shap <- function(x, n = n_shapiro){
    x <- x[is.finite(x)]
    if (length(x) > n) x <- sample(x, n)
    if (length(x) < 3) return(NA_real_)
    shapiro.test(x)$p.value
  }
  p_before <- sapply(dat[, cols], shap)
  p_after  <- sapply(dat_t[, cols], shap)
  
  # 4) lambda
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


