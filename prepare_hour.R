library(data.table)
library(RollingWindow)
library(DescTools)
library(TTR)
library(arrow)
library(dplyr)
library(lubridate)


# SET UP ------------------------------------------------------------------
# Paths
QCDATA = "/home/sn/lean/data"


# SP500 SYMBOLS -----------------------------------------------------------
# Import SP500 symbols from lean
spy_files = list.files(file.path(QCDATA, "equity", "usa", "universes", "etf", "spy"), 
                       full.names = TRUE)
spy_const = lapply(spy_files, fread)
names(spy_const) = gsub("\\.csv", "", basename(spy_files))
spy_const = rbindlist(spy_const, idcol = "id")
spy_const[!is.na(V6)]

# Extract symbols and take unique
spy_symbols = spy_const[, unique(V1)]
spy_symbols = tolower(spy_symbols)
spy_symbols = c(spy_symbols, "spy")


# PRICE DATA --------------------------------------------------------------
# Import hour data for sp500 symbols
prices = open_dataset(file.path(QCDATA, "stocks_hour.csv"), format = "csv") |>
  dplyr::filter(Symbol %in% spy_symbols) |>
  dplyr::rename(symbol = Symbol, date = Date, open = Open, high = High, low = Low, 
                close = Close, volume = Volume, adj_close = `Adj Close`) |>
  dplyr::mutate(date = force_tz(date, "America/New_York")) |>
  collect()
setDT(prices)

# Remove duplicates
prices = unique(prices, by = c("symbol", "date"))

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

# Remove missing values
prices = na.omit(prices)

# Keep observations with at least 2000 hours of observations. We need this because otherwise
# the R studio will crashes
id_n = prices[, .N, by = symbol]
remove_symbols = id_n[N < 2000, unique(symbol)]
prices = prices[symbol %notin% remove_symbols]

# Sort
setorder(prices, symbol, date)

# Calculate returns
prices[, returns := close / shift(close) - 1, by = .(symbol)]

# Remove NA values
prices = na.omit(prices)

# Set SPY returns as market returns
spy_ret = na.omit(prices[symbol == "spy", .(date, market_ret = returns)])
prices = spy_ret[prices, on = "date"]

# Remove columns we don't need
prices = prices[, let(open = NULL, high = NULL, low = NULL)]

# free memory
gc()


# PREDICTORS --------------------------------------------------------------
# Rolling beta
setorder(prices, symbol, date)
prices[, beta := RollingBeta(market_ret, returns, 1000, na_method = "ignore"),
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
months_size = c(15, 22, 22*3, 22*6, 22*12)
mom_vars = paste0("momentum_", months_size)
f_ = function(x, n) {
  shift(x, 5*7) / shift(x, n * 7) - 1
}
prices[, (mom_vars) := lapply(months_size, function(x) f_(close, x)), by = symbol]

# Momentum ensambles
weights_ = c(24, 12, 4, 2, 1) / sum(c(24, 12, 4, 2, 1))
prices[, momentum_average := momentum_15 * weights_[1] +
         momentum_22 * weights_[2] +
         momentum_66 * weights_[3] +
         momentum_132 * weights_[4] +
         momentum_264 * weights_[5]]

# Dollar volume z-score
dv_cols = paste0("dollar_volume_zscore_winsorized", months_size)
f_ = function(x, y, z) RollingZscore(as.numeric(x * y), z, na_method = "ignore")
prices[, (dv_cols) := lapply(months_size, function(s)  as.vector(f_(close_raw, volume, s * 7))),
       by = symbol]
prices[, (dv_cols) := lapply(.SD, function(x) Winsorize(x, val = c(-5, 5))),
       .SDcols = dv_cols]

# Tecnical indicators
prices[, rsi_14 := RSI(close, n = 14*7), by = symbol]
prices[, rsi_22 := RSI(close, n = 22*7), by = symbol]
prices[, rsi_66 := RSI(close, n = 66*7), by = symbol]

# Remove columns we don't need
prices[, c("beta_rank", "beta_rank_pct", "open", "high", "low", "volume") := NULL]

# Create target variable
setorder(prices, symbol, date)
prices[, target := shift(close, 1, type = "lead") / close - 1, by = symbol]
prices = na.omit(prices)

# Create ID column
prices[, id := data.table::rleid(symbol)]

# Save sh file that we qsub on srce
sh_file = sprintf("
#!/bin/bash

#PBS -N ERFHOUR
#PBS -l ncpus=8
#PBS -l mem=10GB
#PBS -J 1-%d
#PBS -o logs_hour
#PBS -j oe

cd ${PBS_O_WORKDIR}
apptainer run image.sif estimate.R
", prices[, length(unique(id))])
sh_file_name = "padobran_hour.sh"
file.create(sh_file_name)
writeLines(sh_file, sh_file_name)

# Save every id (symbol) separately
dirname_ = "data_hour"
if (!dir.exists(dirname_)) dir.create(dirname_)
for (i in prices[, unique(id)]) {
  prices_ = prices[id == i]
  fwrite(prices_, file.path(dirname_, paste0("prices_", i, ".csv")))
}

# Add file to padobran
# scp -r /home/sn/projects_r/alpha_erf/data_hour padobran:/home/jmaric/alpha_erf/data_hour/
