library(data.table)
library(RollingWindow)
library(DescTools)
library(TTR)
# library(janitor)
# library(PerformanceAnalytics)
# library(erf)
# library(foreach)
# library(parallel)
# library(doParallel)


# PRICE DATA --------------------------------------------------------------
# Import QC daily data
if (interactive()) {
  prices = fread("/home/sn/lean/data/stocks_daily.csv")
} else {
  prices = fread("stocks_daily.csv")
}
setnames(prices, gsub(" ", "_", c(tolower(colnames(prices)))))

# Remove duplicates
prices = unique(prices, by = c("symbol", "date"))

# Remove duplicates - there are same for different symbols (eg. phun and phun.1)
dups = prices[, .(symbol , n = .N),
              by = .(date, open, high, low, close, volume, adj_close,
                     symbol_first = substr(symbol, 1, 1))]
dups = dups[n > 1]
dups[, symbol_short := gsub("\\.\\d$", "", symbol)]
symbols_remove = dups[, .(symbol, n = .N),
                      by = .(date, open, high, low, close, volume, adj_close,
                             symbol_short)]
symbols_remove[n >= 2, unique(symbol)]
symbols_remove = symbols_remove[n >= 2, unique(symbol)]
symbols_remove = symbols_remove[grepl("\\.", symbols_remove)]
prices = prices[symbol %notin% symbols_remove]

# Adjust all columns
prices[, adj_rate := adj_close / close]
prices[, let(
  open = open*adj_rate,
  high = high*adj_rate,
  low = low*adj_rate
)]
setnames(prices, "close", "close_raw")
setnames(prices, "adj_close", "close")
prices[, let(adj_rate = NULL)]
setcolorder(prices, c("symbol", "date", "open", "high", "low", "close", "volume"))

# Remove observations where open, high, low, close columns are below 1e-008
# This step is opional, we need it if we will use finfeatures package
prices = prices[open > 1e-008 & high > 1e-008 & low > 1e-008 & close > 1e-008]

# Sort
setorder(prices, symbol, date)

# Calculate returns
prices[, returns := close / shift(close, 1) - 1]

# Remove missing values
prices = na.omit(prices)

# Set SPY returns as market returns
spy_ret = na.omit(prices[symbol == "spy", .(date, market_ret = returns)])
prices = spy_ret[prices, on = "date"]

# Minimal observations per symbol is 253 days
remove_symbols = prices[, .(symbol, n = .N), by = symbol][n < 253, symbol]
prices = prices[symbol %notin% remove_symbols]

# Free memory
gc()


# PREDICTORS --------------------------------------------------------------
# Rolling beta
setorder(prices, symbol, date)
prices = prices[, beta := RollingBeta(market_ret, returns, 252, na_method = "ignore"),
                by = symbol]

# Highest beta by date - 5% of symbols by highest beta
prices[, beta_rank := frank(abs(beta), ties.method = "dense", na.last = "keep"), by = date]
prices[, beta_rank_pct := beta_rank / max(beta_rank, na.rm = TRUE), by = date]
prices[, beta_rank_largest_99 := 0, by = date]
prices[beta_rank_pct > 0.99, beta_rank_largest_99 := 1, by = date]
prices[, beta_rank_largest_95 := 0, by = date]
prices[beta_rank_pct > 0.95, beta_rank_largest_95 := 1, by = date]
prices[, beta_rank_largest_90 := 0, by = date]
prices[beta_rank_pct > 0.90, beta_rank_largest_90 := 1, by = date]
setorder(prices, symbol, date)

# Momentum predictors
months_size = c(3, 6, 9, 12)
mom_vars = paste0("momentum_", months_size)
f_ = function(x, n) {
  shift(x, 21) / shift(x, n * 21) - 1
}
prices[, (mom_vars) := lapply(months_size, function(x) f_(close, x)), by = symbol]

# Momentum ensambles
weights_ = c(12, 6, 3, 1) / sum(c(12, 6, 3, 1))
prices[, momentum_average := momentum_3 * weights_[1] +
         momentum_6 * weights_[2] +
         momentum_9 * weights_[3] +
         momentum_12 * weights_[4]]

# Dolar volume z-score
dv_cols = paste0("dollar_volume_zscore_winsorized", months_size)
f_ = function(x, y, z) RollingZscore(as.numeric(x * y), z, na_method = "ignore")
prices[, (dv_cols) := lapply(months_size, function(s)  as.vector(f_(close_raw, volume, s * 21))),
       by = symbol]
prices[, (dv_cols) := lapply(.SD, function(x) Winsorize(x, val = c(-5, 5))),
       .SDcols = dv_cols]

# Tecnical indicators
prices[, rsi := RSI(close, n = 14), by = symbol]

# Remove columns we don't need
prices[, c("beta_rank", "beta_rank_pct", "open", "high", "low", "volume") := NULL]

# Save local temporarily
if (!dir.exists("data")) dir.create("data")
fwrite(prices, "data/prices.csv")

# Add file to padobran

