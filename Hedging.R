### DEFAULT STUFF STARTS ###

values = function(x){
  x <- as.matrix(x)
  res <- matrix(NA, nrow = ncol(x), ncol = 4)
  
  for (i in 1:ncol(x)){
    res[i,] <- matrix(c(mean(x[,i]), sd(x[,i]), quantile(x[,i], 0.025), quantile(x[,i], 0.975)), nrow=1)
  }
  
  colnames(res) = c("Mean", "Std.Dev.", "5%", "95%")
  res
}

require(zoo)
require(boot)
require(PerformanceAnalytics)
require(rugarch)
require(FinTS)
require(rmgarch)
require(timeSeries) 
require(pastecs)
require(tseries)
require(Hmisc)
require(gogarch)
require(parallel)
require(fpp)
require(texreg)
require(broom)    # not installed
require(lubridate)
require(ggplot2)
require(scales)
require(gridExtra)
require(ggthemes)
require(ggfortify)    # not installed
require(bbplot)       # not installed
require(tidyverse)    # not installed
require(onewaytests)
library(abind)

data_df = read.csv(file = "/Users/szonyid/Documents/MacBook/School/Corvinus/TDK/TDK/R_Code/Data.csv")

data_df$Date <- as.Date(data_df$Date)
data_df <- data_df[!duplicated(data_df$Date), ]
data_zoo <- read.zoo(data_df, order.by = data_df$Date)

# Take as individual time series
sp500 = data_zoo[, "S.P.500", drop = FALSE]
bond = data_zoo[, "BOND", drop = FALSE]
oil = data_zoo[, "OIL", drop = FALSE]
gold = data_zoo[, "GOLD", drop = FALSE]
eurusd = data_zoo[, "EUR.USD", drop = FALSE]

# Calculate log-returns for M-GARCH analysis
sp500_ret = 100 * CalculateReturns(sp500, method="log")
bond_ret = 100 * CalculateReturns(bond, method="log")
oil_ret = 100 * CalculateReturns(oil, method="log")
gold_ret = 100 * CalculateReturns(gold, method="log")
eurusd_ret = 100 * CalculateReturns(eurusd, method="log")

# Remove first NA observation
sp500_ret = sp500_ret[-1,]
bond_ret = bond_ret[-1,]
oil_ret = oil_ret[-1,]
gold_ret = gold_ret[-1,]
eurusd_ret = eurusd_ret[-1,]

# Check if there's any remaining NA's
summary(oil_ret)

# Linear interpolation
oil_ret <- na.approx(oil_ret)

# Combine variables
y_ret = merge(sp500_ret, bond_ret, oil_ret, gold_ret, eurusd_ret)

# Model specification
dcc_garch_spec = ugarchspec(mean.model = list(armaOrder = c(1,0)), 
                            variance.model = list(garchOrder = c(1,1), 
                                                  model = "sGARCH"), 
                            distribution.model = "std")

dcc_garch = dccspec(uspec = multispec(replicate(5, dcc_garch_spec)), 
                    dccOrder = c(1,1), 
                    distribution = "mvt")

# Set number of model refits
r.e = 20

# DCC and GO model forecasts
cl = makePSOCKcluster(4)
dcc_garch_roll = dccroll(dcc_garch, y_ret, n.ahead = 1,
                         n.start = 4452,
                         refit.every= r.e, refit.window = "moving", cluster = cl)
stopCluster(cl) 

go_garch_spec = gogarchspec(mean.model = list( model='AR',lag=1), variance.model = list(model = "sGARCH", garchOrder = c(1, 1),
                                                                                        variance.targeting = FALSE), distribution = 'maniga', ica = 'fastica')

cl = makePSOCKcluster(4)
go_garch_roll = gogarchroll(go_garch_spec, y_ret, n.ahead = 1,
                            n.start = 4452,
                            refit.every = r.e, refit.window = "moving", rseed = 1, cluster = cl)
stopCluster(cl) 




data <- y_ret

colnames(data) <- c("sp500", "bond", "oil", "gold", "eurusd")

H = rcov(go_garch_roll)

data_total_length <- length(data) / ncol(data)

length_split <- data_total_length - length(H) / 25

splitted_data_out <- data[-c(1:length_split),]

k = ncol(data) 
NAMES = matrix(NA, ncol = k, nrow = k) 

for (i in 1:k){
  NAMES[i,] = paste0(colnames(data), "/", colnames(data)[i]) 
}

DM = array(NA, c(k, k, nrow(splitted_data_out)))

for (i in 1:k){
  for (j in 1:k){
    DM[i,j,] = H[i,j,] / H[j,j,] 
  }
}

HR = array(1, c(k, 4, k)) 

for (i in 1:k){
  for (j in 1:k){
    HR[i,,j] = values(DM[i,j,])
  }
}

colnames(HR) = c("Mean", "Std.Dev.", "5%", "95%")
dimnames(HR)[[1]] = colnames(data)

for (i in 1:k){
  if (i == 1){
    HRatio = HR[,,i]
  } else {
    HRatio = rbind(HRatio, HR[,,i])
  }
}

rownames(HRatio) = c(t(NAMES))
pval = HE = matrix(NA, ncol = k, nrow = k)
colnames(HE) = rownames(HE) = colnames(pval) = rownames(pval) = colnames(data)

portfolio_return = array(NA, c(k, k, nrow(splitted_data_out)))

for (i in 1:k) {
  for (j in 1:k) {
    portfolio_return[i,j,] = splitted_data_out[,i] - DM[i,j,] * splitted_data_out[,j]
    
    HE[i,j] = 1 - var(portfolio_return[i,j,]) / var(splitted_data_out[,i])
    
    pval[i,j] = var.test(x = portfolio_return[i,j,], y = splitted_data_out[,i], ratio = 1)$p.value
  }
}

# Convert the zoo object to a data.frame
splitted_data_out_df <- data.frame(index(splitted_data_out), coredata(splitted_data_out))
colnames(splitted_data_out_df) <- c("date", colnames(splitted_data_out))

# Use aperm() to rearrange dimensions of portfolio_return
portfolio_return_2D <- aperm(portfolio_return, c(3, 1, 2))
# Reshape to 2D matrix
portfolio_return_2D <- matrix(portfolio_return_2D, nrow = nrow(splitted_data_out), ncol = k)
colnames(portfolio_return_2D) <-  c("sp500", "bond", "oil", "gold", "eurusd")

hedged_data <- data.frame(index(portfolio_return_2D), coredata(portfolio_return_2D))
hedged_data <- hedged_data[,-1]
unhedged_data <- splitted_data_out_df
unhedged_data <- unhedged_data[,-1]

# create empty data frame to store p-values
results <- data.frame(p.value = numeric(), stringsAsFactors = FALSE)

# create list of tickers for the hedged portfolio
tickers <- c("bond", "oil", "gold", "eurusd")

# get unhedged S&P 500 returns
unhedged_returns <- unhedged_data[,1]

# loop through tickers
for (ticker in tickers) {
  
  # extract returns for hedged portfolio
  hedged_returns <- hedged_data[[ticker]]
  
  # melt the data and add a column indicating the hedging status
  melted_hedged <- data.frame(returns = hedged_returns, group = "Hedged")
  melted_unhedged <- data.frame(returns = unhedged_returns, group = "Unhedged")
  
  # bind the data frames
  bond_returns <- rbind(melted_hedged, melted_unhedged)
  
  # run the Brown-Forsythe test using bf.test and extract the p-value
  bf_test <- bf.test(returns ~ group, data = bond_returns)
  p_value <- bf_test$p.value
  
  # append p-value to data frame
  results <- rbind(results, data.frame(p.value = p_value, stringsAsFactors = FALSE))
}

table <- cbind(HRatio, c(t(HE)), c(t(pval)))
table <- table[-which(table[,1] == 1),]
table <- table[grep("sp500/", rownames(table)),]
table <- cbind(table, results)
colnames(table) = c("Mean", "Std.Dev.", "5%", "95%", "HE", "F_p-value", "BF_p-value")

round(table,3)








# Create a data frame with the summary statistics and hedging effectiveness
summary_df <- data.frame(
  Asset_Pair = rep(c("S&P 500 / Bond", "S&P 500 / Oil", "S&P 500 / Gold", "S&P 500 / EUR/USD"), each = 2),
  Model = rep(c("DCC-GARCH", "GO-GARCH"), times = 4),
  Hedge_Ratio_Mean = c(-0.818, -1.269, 0.121, 0.116, 0.003, -0.044, 0.301, 0.488),
  Hedge_Ratio_SD = c(1.055, 0.998, 0.068, 0.078, 0.218, 0.139, 0.585, 0.316),
  Hedging_Effectiveness = c(7, -1.5, 1.7, 1.8, 1.8, -0.9, 10.4, 2.5),
  Ftest_pvalue = c(0.224, 0.801, 0.776, 0.758, 0.757, 0.883, 0.068, 0.674),
  BFtest_pvalue = c(0.443, 0.441, 0.794, 0.708, 0.918, 0.912, 0.340, 0.342)
)


ggplot(summary_df, aes(x = Asset_Pair, y = Hedging_Effectiveness, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label=round(Hedging_Effectiveness,1)), position=position_stack(vjust=0.5), size=3.5, color="white", fontface="bold") +
  ggtitle("Hedging Effectiveness by Asset Pair and M-GARCH Model") +
  xlab("") + 
  ylab("Hedging Effectiveness (%)") +
  scale_fill_manual("", 
                    breaks = c("DCC-GARCH", "GO-GARCH"),
                    values = c("seagreen", "firebrick1")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.key.size = unit(3, 'cm'),
        legend.key.width= unit(1, 'cm'),
        legend.text = element_text(size = 9, face = "bold"))



library()

# Generate sample data
data <- tibble(
  Asset = c(rep("Bond", 4), rep("Oil", 4), rep("Gold", 4), rep("EUR/USD", 4)),
  Model = rep(c("Original", "Refit = 50", "Refit = 100", "Normal Dist."), 4),
  HE = c(7, 7, 7, 7, 1.7, 1.8, 1.9, 0.8, 1.8, 1.8, 1.8, 2, 10.4, 10.4, 10.4, 10.8)
)

# Define the order of the columns
cols <- c("Bond", "Oil", "Gold", "EUR/USD")

# Create the grouped bar chart
ggplot(data, aes(x = Asset, y = HE, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(aes(label = HE), position = position_dodge(width = 0.9), vjust = -0.5, fontface = "bold") +
  scale_fill_manual("", 
                    breaks = c("Original", "Refit = 50", "Refit = 100", "Normal Dist."),
                    values = c("seagreen", "firebrick1", "deepskyblue", "darkorange")) +
  labs(title = "Hedging Effectiveness by Asset Pair and Simulation Type for DCC-GARCH",
       x = "",
       y = "Hedging Effectiveness",
       fill = "Simulation Type") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.key.size = unit(1.5, 'lines'),
        legend.key.width = unit(0.5, 'lines'),
        legend.text = element_text(size = 9, face = "bold")) +
  scale_x_discrete(limits = cols) +
  scale_y_continuous(limits = c(-2, 12), expand = c(0, 0))





