# import data from padobran
# scp -r padobran:/home/jmaric/alpha_erf/results_hour/ /home/sn/data/strategies/alpha_erf/ 

library(data.table)
library(erf)
library(AzureStor)
library(arrow)
library(dplyr)
library(PerformanceAnalytics)
library(ggplot2)
library(TTR)
library(lubridate)


# # CHECKS ------------------------------------------------------------------
# # Check have many symbols are in total
# files = list.files("data_hour", full.names = TRUE, pattern = "csv")
# symbols = lapply(files, fread, select = "symbol") 
# symbols = rbindlist(symbols)
# symbols[, length(unique(symbol))]

# DATA --------------------------------------------------------------------
# Set up
PATH = "/home/sn/data/strategies/alpha_erf/results_hour"
blobendpoint = storage_endpoint(Sys.getenv("BLOB-ENDPOINT-SNP"),
                                key=Sys.getenv("BLOB-KEY-SNP"))
cont = storage_container(blobendpoint, "qc-backtest")
qlcal::calendars
qlcal::setCalendar("UnitedStates/NYSE")

# Import results
files = list.files(PATH, full.names = TRUE)
erf_predictions = lapply(files, fread)
erf_predictions = erf_predictions[lengths(erf_predictions) > 0]

# Combine all results
erf_predictions = rbindlist(erf_predictions)

# Check timezones
attr(erf_predictions$date, "tzone")
Sys.timezone()
erf_predictions[, date := with_tz(date, tz = "America/New_York")]

# Summary
erf_predictions[, .(
  min_date = min(date),
  max_date = max(date),
  number_observations = .N
)]

# Hour prices
QCDATA = "/home/sn/lean/data"
prices = open_dataset(file.path(QCDATA, "stocks_hour.csv"), format = "csv") |>
  dplyr::filter(Symbol %in% erf_predictions[, unique(symbol)]) |>
  dplyr::rename(symbol = Symbol, date = Date, open = Open, high = High, low = Low, 
                close = Close, volume = Volume, adj_close = `Adj Close`) |>
  dplyr::mutate(date = force_tz(date, "America/New_York")) |>
  collect()
setDT(prices)
attr(prices$date, "tz")
prices = unique(prices, by = c("symbol", "date"))
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
prices = prices[open > 1e-008 & high > 1e-008 & low > 1e-008 & close > 1e-008]
setorder(prices, symbol, date)
prices[, returns := close / shift(close, 1) - 1]
prices[, returns_intra := close / open - 1]
prices = na.omit(prices)
spy_ret = na.omit(prices[symbol == "spy", .(date, market_ret = returns)])
prices = spy_ret[prices, on = "date"]
remove_symbols = prices[, .(symbol, n = .N), by = symbol][n < 253, symbol]
prices = prices[symbol %notin% remove_symbols]
gc()

# # Check datetimes
# erf_predictions[, unique(symbol)]
# s = "cars"
# x = erf_predictions[symbol == s, .(symbol, date, targetr)]
# y = prices[symbol == s, .(symbol, date, returns)]
# i = c("2024-01-01", "2024-02-01")
# tail(x[date %between% i], 20); tail(y[date %between% i], 20);
# 
# # Move by one hour !?
# erf_predictions[, date_erf := date]
# prices[, date_prices := date]
# test = prices[erf_predictions, on = c("symbol", "date"), roll = Inf]
# test[, .(symbol, date, date_erf, date_prices)]
# tail(test[symbol == "ctx", .(symbol, date, date_erf, date_prices)], 15)


# PREPARE -----------------------------------------------------------------
# Convert all predictors to numeric
erf_predictions[, names(.SD) := lapply(.SD, as.numeric), .SDcols = -c("symbol", "date")]

# Define quantile ratios
erf_predictions[, qr_995 := upper_quantile_0_995 / lower_quantile_0_995]
erf_predictions[, qr_99 := upper_quantile_0_99 / lower_quantile_0_99]
erf_predictions[, qr_98 := upper_quantile_0_98 / lower_quantile_0_98]
erf_predictions[, qr_95 := upper_quantile_0_95 / lower_quantile_0_95]

# Standardize
setorder(erf_predictions, symbol, date)
erf_predictions[, stand_q95 := roll::roll_scale(qr_95, nrow(.SD), min_obs = 44*7), by = symbol]
erf_predictions[, stand_q99 := roll::roll_scale(qr_99, nrow(.SD), min_obs = 44*7), by = symbol]
erf_predictions[, stand_q98 := roll::roll_scale(qr_98, nrow(.SD), min_obs = 44*7), by = symbol]
erf_predictions[, stand_q995 := roll::roll_scale(qr_995, nrow(.SD), min_obs = 44*7), by = symbol]

# # Merge with prices
# dt = prices[, .(date, symbol, returns)][erf_predictions, on = c("symbol", "date")]

# Generate signals
generate_signal = function(dt, threshold, n = 0, q = c("99", "95", "all")) {
  dt_ = copy(dt)
  dt_[, signal := 0]
  if (q == "995") {
    dt_[stand_q995 > threshold, signal := 1]
  } else if (q == "99") {
    dt_[stand_q99 > threshold, signal := 1]
  } else if (q == "98") {
    dt_[stand_q98 > threshold, signal := 1]
  } else if (q == "95") {
    dt_[stand_q95 > threshold, signal := 1]
  } else {
    dt_[stand_q995 > threshold & stand_q99 > threshold & stand_q98 > threshold & stand_q95 > threshold, signal := 1]
  }
  dt_[, signal := shift(signal, n), by = symbol]
  dt_[, strategy_return := targetr * signal]
  return(na.omit(dt_))
}

# INDIVIDUAL ASSETS -------------------------------------------------------
# Sample
symbol_ = erf_predictions[, sample(unique(symbol), 1)]
# symbol_ = "msft"
dt_ = erf_predictions[symbol == symbol_]

# Plot quantiles and returns
cols = c("date", "upper_quantile_0_99", "lower_quantile_0_99",
         "upper_quantile_0_98", "lower_quantile_0_98",
         "upper_quantile_0_995", "lower_quantile_0_995",
         "upper_quantile_0_95", "lower_quantile_0_95", "targetr")
ggplot(dt_[, ..cols], aes(x = date)) +
  geom_line(aes(y = upper_quantile_0_99, color = "upper_quantile_0_99")) +
  geom_line(aes(y = -lower_quantile_0_99, color = "lower_quantile_0_99")) +
  geom_line(aes(y = upper_quantile_0_95, color = "upper_quantile_0_95")) +
  geom_line(aes(y = -lower_quantile_0_95, color = "lower_quantile_0_95")) +
  geom_line(aes(y = upper_quantile_0_98, color = "upper_quantile_0_98")) +
  geom_line(aes(y = -lower_quantile_0_98, color = "lower_quantile_0_98")) +
  geom_line(aes(y = upper_quantile_0_995, color = "upper_quantile_0_995")) +
  geom_line(aes(y = -lower_quantile_0_995, color = "lower_quantile_0_995")) +
  geom_line(aes(y = targetr, color = "targetr")) +
  scale_color_manual(values = c("upper_quantile_0_99" = "red",
                                "lower_quantile_0_99" = "red",
                                "upper_quantile_0_95" = "blue",
                                "lower_quantile_0_95" = "blue",
                                "upper_quantile_0_98" = "green",
                                "lower_quantile_0_98" = "green",
                                "upper_quantile_0_995" = "purple",
                                "lower_quantile_0_995" = "purple",
                                "targetr" = "black")) +
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal()
n_ = 7*5
ggplot(dt_[, ..cols], aes(x = date)) +
  geom_line(aes(y = EMA(upper_quantile_0_99, n_), color = "upper_quantile_0_99")) +
  geom_line(aes(y = EMA(-lower_quantile_0_99, n_), color = "lower_quantile_0_99")) +
  geom_line(aes(y = EMA(upper_quantile_0_95, n_), color = "upper_quantile_0_95")) +
  geom_line(aes(y = EMA(-lower_quantile_0_95, n_), color = "lower_quantile_0_95")) +
  geom_line(aes(y = EMA(upper_quantile_0_98, n_), color = "upper_quantile_0_98")) +
  geom_line(aes(y = EMA(-lower_quantile_0_98, n_), color = "lower_quantile_0_98")) +
  geom_line(aes(y = EMA(upper_quantile_0_995, n_), color = "upper_quantile_0_995")) +
  geom_line(aes(y = EMA(-lower_quantile_0_995, n_), color = "lower_quantile_0_995")) +
  geom_line(aes(y = EMA(targetr, n_), color = "targetr")) +
  scale_color_manual(values = c("upper_quantile_0_99" = "red",
                                "lower_quantile_0_99" = "red",
                                "upper_quantile_0_95" = "blue",
                                "lower_quantile_0_95" = "blue",
                                "upper_quantile_0_98" = "green",
                                "lower_quantile_0_98" = "green",
                                "upper_quantile_0_995" = "purple",
                                "lower_quantile_0_995" = "purple",
                                "targetr" = "black")) +
  # scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal()

# Scatterplot between signal and returns
ggplot(dt_, aes(x = stand_q95, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_, aes(x = stand_q98, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_, aes(x = stand_q99, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_, aes(x = stand_q995, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_[stand_q95 > 0], aes(x = stand_q95, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_[stand_q98 > 0], aes(x = stand_q98, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_[stand_q99 > 0], aes(x = stand_q99, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_[stand_q995 > 0], aes(x = stand_q995, y = targetr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_minimal()
ggplot(dt_[stand_q95 > 0], aes(x = stand_q95, y = targetr)) +
  geom_point() +
  geom_smooth(method = "gam") +
  theme_minimal()


# Plot cumulative returns
# symbol_ = "invh"
symbol_ = erf_predictions[, sample(unique(symbol), 1)]
back_ = erf_predictions[symbol == symbol_]
back_ = generate_signal(back_, 0, n = 0, q = "995")
back_xts = as.xts.data.table(back_[, .(date, strategy_return, targetr)])
Return.cumulative(back_xts)
SharpeRatio(back_xts)
maxDrawdown(back_xts)
charts.PerformanceSummary(back_xts, main = back_[1, symbol])
charts.PerformanceSummary(back_xts["2010/"], main = back_[1, symbol])

# Add to QC to test the strategy\
qc_data = back_[, .(date, stand_q95, stand_q98, stand_q99, stand_q995)]
qc_data[, date := paste0(as.character(date), " 16:00:00")]
storage_write_csv(qc_data, cont, paste0("erf_", back_[1, symbol], ".csv"))


# INSAMPLE STOCK OPTIMIZATION ---------------------------------------------
# # Parameters
# params = expand.grid(s = erf_predictions[, unique(symbol)],
#                      q = c("995", "99", "98", "95", "all"),
#                      t = seq(0.5, 2, 0.1),
#                      stringsAsFactors = FALSE)
# 
# # Calculate portfolio performance for all parameters
# results = list()
# for (i in 1:nrow(params)) {
#   dt_ = erf_predictions[symbol == params$s[i]]
#   dt_ = generate_signal(copy(dt_), params$t[i], 0, params$q[i])
#   back_xts = as.xts.data.table(dt_[, .(date, strategy_return, targetr)])
#   x = table.AnnualizedReturns(back_xts)
#   x = as.data.table(x, keep.rownames = "var")
#   x = data.table::transpose(x[, 1:2], make.names = "var")
#   x = clean_names(x)
#   results[[i]] = cbind(s = params$s[i], q = params$q[i], t = params$t[i], x)
# }
# results_df = rbindlist(results)
# setnames(results_df, c("symbol", "q", "threshold", "cagr", "std", "SR"))
# 
# # list best
# setorder(results_df, -SR)
# head(na.omit(results_df), 10)
# 
# # list best across symbols
# setorder(results_df, symbol, -SR)
# na.omit(results_df)[, head(.SD), by = symbol]


# PORTFOLIO ---------------------------------------------------------------
# Portfolio returns
portfolio = erf_predictions[, .(symbol, date, stand_q95, stand_q995, stand_q99, targetr)]
portfolio[, signal := 0]
portfolio[shift(stand_q95) > 0.5, signal := 1]
# portfolio[, signal := shift(signal, 1), by = symbol]
portfolio[, weights := signal / nrow(.SD[signal == 1]), by = date]
setorder(portfolio, date)
portfolio_ret = portfolio[, .(returns = sum(targetr * weights, na.rm = TRUE)), by = date]
portfolio_ret = as.xts.data.table(portfolio_ret)
# table.AnnualizedReturns(portfolio_ret)
maxDrawdown(portfolio_ret)
charts.PerformanceSummary(portfolio_ret)
charts.PerformanceSummary(portfolio_ret["2020-01/"])

# Add all to QC
qc_data = erf_predictions[stand_q95 > 2, .(symbol, date, stand_q995, stand_q95)]
qc_data[, date := as.Date(vapply(as.Date(date), qlcal::advanceDate, FUN.VALUE = double(1L), days = 1))]
qc_data[, date := paste0(as.character(date), " 15:59:00")]
qc_data = qc_data[, .(
  symbol = paste0(symbol, collapse = "|"),
  qr_99 = paste0(stand_q995, collapse = "|"),
  qr_95 = paste0(stand_q95, collapse = "|")),
  by = date]
setorder(qc_data, date)
storage_write_csv(qc_data, cont, "erf.csv")
colnames(qc_data)

# PORTFOLIO CROSS SECTION -------------------------------------------------
library(portsort)
# Downsample to monthly frequency
dtm = copy(prices)
dtm[, year_month_id := lubridate::ceiling_date(date, unit = "day")]
dtm = dtm[, .(
  date = tail(date, 1),
  open = head(open, 1),
  high = max(high),
  low = min(low),
  close = tail(close, 1),
  close_raw = tail(close_raw, 1),
  volume_mean = mean(volume, na.rm = TRUE),
  volume = sum(volume, na.rm = TRUE)
), by = c("symbol", "year_month_id")]
dtm = erf_predictions[, .(symbol, date, stand_q95, stand_q98, stand_q99, stand_q995, targetr)][
  dtm, on = c("symbol", "date")]
dtm = na.omit(dtm)
dtm = dtm[close_raw > 1]

# create forward returns
setorder(dtm, symbol, date)
dtm[, ret_forward := shift(close, -1, type = "shift") / close - 1, by = symbol]
dtm[, .(symbol, year_month_id, close, ret_forward)]
dtm = na.omit(dtm)

# Portfolio sort
predictors = c("stand_q95", "stand_q98", "stand_q99", "stand_q995")
results = list()
for (i in seq_along(predictors)) {
  print(i)
  predictor = predictors[i]
  print(predictor)
  
  # sample data we need
  sample = dtm[, .SD, .SDcols = c("symbol", "year_month_id", "ret_forward", predictor)]
  
  # check for duplicates
  dup_index = which(sample[, duplicated(.SD[, .(symbol, year_month_id)])])
  sample = unique(sample, by = c("symbol", "year_month_id"))
  
  # remove NA values for target
  sample = na.omit(sample)
  
  # predictors matrix
  Fa = dcast(sample, year_month_id ~ symbol, value.var = predictor)
  setorder(Fa, year_month_id)
  Fa = as.xts.data.table(Fa)
  cols_with_na = apply(Fa, MARGIN = 2, FUN = function(x) sum(is.na(x)) > as.integer(nrow(Fa) * 0.6))
  Fa = Fa[, !cols_with_na]
  
  # remove all NA values
  rows_with_na <- apply(Fa, MARGIN = 1, FUN = function(x) sum(is.na(x)) > as.integer(ncol(Fa) * 0.99))
  Fa = Fa[!rows_with_na, ]
  
  # Forward returns
  R_forward = as.xts.data.table(dcast(sample, year_month_id ~ symbol, value.var = "ret_forward"))
  R_forward = R_forward[, !cols_with_na]
  R_forward = R_forward[!rows_with_na, ]
  
  # uncondtitional sorting
  dimA = 0:20/20
  psort_out = tryCatch({
    unconditional.sort(
      Fa,
      Fb = NULL,
      Fc = NULL,
      R_forward,
      dimA = dimA,
      dimB = NULL,
      dimC = NULL,
      type = 7
    )
    
  }, error = function(e) NULL)
  
  portsort_dt = table.AnnualizedReturns(psort_out$returns)
  results[[i]] = setDT(portsort_dt, keep.rownames = TRUE)
}

# Inspect results
names(results) = predictors
results_dt = rbindlist(results, idcol = "predictor")
results_dt[grep("Sharpe", rn)]
results_dt[grep("Retu", rn)]

# Create quantile portfolio for every week
portfolio_sort_week = erf_predictions[, .(symbol, date, stand_q95, targetr)]
portfolio_sort_week = na.omit(portfolio_sort_week)
setorder(portfolio_sort_week, date, -stand_q95)
portfolio_sort_week = portfolio_sort_week[, head(.SD, 50), by = date]
portfolio_sort_week = portfolio_sort_week[date < as.Date("2024-10-01")]
portfolio_sort_week[, weight := 1 / nrow(.SD), by = date]
portfolio_returns = portfolio_sort_week[, sum(targetr * weight, na.rm = TRUE), by = date]
charts.PerformanceSummary(as.xts.data.table(portfolio_returns))
charts.PerformanceSummary(as.xts.data.table(portfolio_returns)["2024"])

# Save to Quantconnect
qc_data = portfolio_sort_week[stand_q95 > 0]
qc_data[, date := paste0(as.character(date), " 16:00:00")]
qc_data = qc_data[, .(
  symbol = paste0(symbol, collapse = "|"),
  qr_95 = paste0(stand_q95, collapse = "|")),
  by = date]
setorder(qc_data, date)
storage_write_csv(qc_data, cont, "erf_sort.csv")

# SYSTEMIC RISK -----------------------------------------------------------
# SPY
spy = open_dataset(file.path(QCDATA, "stocks_hour.csv"), format = "csv") |>
  dplyr::filter(Symbol == "spy") |>
  dplyr::select(Date, Open, Close, `Adj Close`) |>
  dplyr::rename(date = Date, open = Open, close_raw = Close, close = `Adj Close`) |>
  dplyr::mutate(open = (close / close_raw) * open) |>
  dplyr::select(date, open, close) |>
  dplyr::mutate(date = force_tz(date, "America/New_York")) |>
  collect()
spy = as.data.table(spy)
spy[, returns := close / shift(close) - 1]
spy[, returns_intra := close / open - 1]
spy = na.omit(spy)
plot(as.xts.data.table(spy[, .(date, close)]))
plot(as.xts.data.table(spy[, .(date, returns)]))
charts.PerformanceSummary(as.xts.data.table(spy[, .(date, returns)]))
charts.PerformanceSummary(as.xts.data.table(spy[, .(date, returns_intra)]))

# Simple aggregation
sr = erf_predictions[, .(symbol, date, stand_q95, stand_q995, targetr)]
setorder(sr, symbol, date)
sr = na.omit(sr)
sr = sr[, .(
  mean = mean(stand_q95, na.rm = TRUE),
  median = median(stand_q95, na.rm = TRUE),
  sd = sd(stand_q95, na.rm = TRUE)
), by = date]
plot(as.xts.data.table(sr),
     main = "Systemic Risk",
     ylab = "Standardized Quantile 95",
     col = c("blue", "red", "green"),
     lwd = 2,
     legend.loc = "topleft")
sr_ema = copy(sr)
sr_ema = sr[, names(.SD) := lapply(.SD, TTR::EMA, n = 7*2), .SDcols = 2:ncol(sr)]
plot(as.xts.data.table(sr_ema),
     main = "Systemic Risk",
     ylab = "Standardized Quantile 95",
     col = c("blue", "red", "green"),
     lwd = 2,
     legend.loc = "topleft")

# Backtest
back = sr_ema[spy[, .(date, returns, returns_intra)], on = "date"]
setorder(back, date)
# predictors = back[, colnames(.SD), .SDcols = -c("date", "returns")]
# back[, (predictors) := lapply(.SD, shift, n = 1), .SDcols = predictors]
back[, strategy := fifelse(shift(mean) > 0, 0, returns)]
# back[, strategy := fifelse(shift(sd) > 1, 0, returns)]
back[is.na(strategy), strategy := 0]
back_xts = as.xts.data.table(na.omit(back[, .(date, strategy, returns)]))
# table.AnnualizedReturns(back_xts)
maxDrawdown(back_xts)
SharpeRatio(back_xts)
charts.PerformanceSummary(back_xts)
charts.PerformanceSummary(back_xts["2020/"])

# Add to Quantconnect
qc_data = sr[spy[, .(date, returns)], on = "date"]
setorder(qc_data, date)
qc_data = qc_data[, .(date = as.character(date), mean)]
qc_data = na.omit(qc_data)
storage_write_csv(qc_data, cont, paste0("erf_risk_hour.csv"))
