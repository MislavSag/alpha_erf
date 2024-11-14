suppressMessages(library(data.table)) 
suppressMessages(library(erf))
suppressMessages(library(janitor))


# setup
SAVEPATH = file.path("results_hour")
if (!dir.exists(SAVEPATH)) dir.create(SAVEPATH, recursive = TRUE)

# Get index
i = as.integer(Sys.getenv('PBS_ARRAY_INDEX'))
# i = 1

# Import data
print("Import data")
file_ = paste0("data_hour/prices_", i, ".csv")
dt = fread(file_, drop = c("close_raw", "close", "returns"))

# Select some date hour interval
print("Define dates")
dates = fread("dates.csv", sep = ",")
tz_ = attr(dates[, date], "tz")
times_1 = data.table::as.ITime(dates[, date])
dates[, date := as.POSIXct(format(date, tz = tz_, usetz = TRUE), tz = "America/New_York")]
times_2 = data.table::as.ITime(dates[, date])
all(times_1 == times_2) # have to be true
setorder(dates, date)
dates = dates[, date]

# Loop over date and symbols and train erf model extract predictions
print("Variables and params")
quantile_levels = c(0.0005, 0.01, 0.02, 0.05, 0.95, 0.98, 0.99, 0.995)
setorder(dt, symbol, date)
predictors = dt[, colnames(.SD), 
                .SDcols = -c("id", "symbol", "date", "target")]

# Check if already estimated
s = dt[, data.table::first(symbol)]
file_name = file.path(SAVEPATH, paste0(s, ".csv"))
if (file.exists(file_name)) next

print(dt)
print(s)
print(dates)

# estimation
l = list()
for (i in 1:length(dates)) { # 1:length(dates
  d = dates[i]
  # d = dates[10000]
  
  # Train data
  dtd = dt[date < d]
  if (nrow(dtd) < 7 * 252) return(NULL)
  if (d - dt[, max(date)] > 1) return(NULL)
  
  # Test data
  test_data = dt[date == d]
  if (nrow(test_data) == 0) return(NULL)
  
  # Fit model for upper
  train_data_upper = dtd[target > 0]
  erf_model_upper = erf(
    X = as.matrix(train_data_upper[, .SD, .SDcols = predictors]),
    Y = train_data_upper[, target],
    min.node.size = 5,
    lambda = 0.001,
    intermediate_quantile = 0.8
  )
  
  # Fit model for lower
  train_data_lower = dtd[target < 0]
  erf_model_lower = erf(
    X = as.matrix(train_data_lower[, .SD, .SDcols = predictors]),
    Y = -train_data_lower[, target],
    min.node.size = 5,
    lambda = 0.001,
    intermediate_quantile = 0.8
  )
  
  # Predict
  erf_predictions_upper = predict(erf_model_upper, as.matrix(test_data[, .SD, .SDcols = predictors]), quantiles = quantile_levels)
  erf_predictions_lower = predict(erf_model_lower, as.matrix(test_data[, .SD, .SDcols = predictors]), quantiles = quantile_levels)
  erf_predictions_upper = clean_names(as.data.frame(erf_predictions_upper))
  colnames(erf_predictions_upper) = paste0("upper_", colnames(erf_predictions_upper))
  erf_predictions_lower = clean_names(as.data.frame(erf_predictions_lower))
  colnames(erf_predictions_lower) = paste0("lower_", colnames(erf_predictions_lower))
  l[[i]] = cbind(
    symbol = s,
    date = d,
    erf_predictions_upper,
    erf_predictions_lower,
    targetr = test_data[, target]
  )
}
  
# Clean and save
x = rbindlist(l)
fwrite(x, file_name)
