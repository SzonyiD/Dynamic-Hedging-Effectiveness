#####################################
### DYNAMIC HEDGING EFFECTIVENESS ### 
### SZŐNYI DÁVID ###
### 2022.01.28. ###
#####################################

#################
### FUNCTIONS ###
#################

SummaryStatistics = function(data){
  moments = matrix(NA, ncol=ncol(data), nrow=14)
  colnames(moments)=colnames(data)
  rownames(moments)=c("Mean","Variance","Skewness","","Kurtosis","","JB","","ERS","","Q(20)","","Q2(20)","")
  for (i in 1:ncol(data)){
    moments[1,i] = mean(data[,i])
    moments[2,i] = var(data[,i])
    skew = moments::agostino.test(data[,i])
    moments[3,i] = skew$statistic[1]
    moments[4,i] = skew$p.value
    kurt = moments::anscombe.test(data[,i])
    moments[5,i] = kurt$statistic[1]-3
    moments[6,i] = kurt$p.value
    jb = moments::jarque.test(as.numeric(data[,i]))
    moments[7,i] = jb$statistic
    moments[8,i] = jb$p.value
    ers = urca::ur.ers(data[,i],type="DF-GLS",model="constant")
    moments[9,i] = ers@teststat
    moments[10,i]= ers@testreg$coefficients[1,4]
    bt = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20)
    moments[11,i] = bt$statistic
    moments[12,i] = bt$p.value
    bt2 = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20, sqrd.res=T)
    moments[13,i] = bt2$statistic
    moments[14,i] = bt2$p.value
  }
  
  cc=c(4,6,8,10,12,14)
  moments = round(moments,3)
  moments1 = moments
  for (j in 1:k){
    for (i in 1:length(cc)){
      i = cc[i]
      if (moments[i,j]<=0.01) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"***",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.05) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"**",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.10) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"*",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else {
        moments1[(i-1),j] = format(round(moments[(i-1),j],3),nsmall=3)
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      }
    }
  }
  
  for (j in 1:k){
    i = 9
    if (moments[i,j]<=-2.57) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"***",sep="")
    } else if (moments[i,j]<=-1.96) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"**",sep="")
    } else if (moments[i,j]<=-1.62) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"*",sep="")
    } else {
      moments1[i,j] = format(round(moments[i,j],3),nsmall=3)
    }
  }
  moments1
}

#####################
### DATA ANALYSIS ###
#####################

setwd("/Users/szonyid/Documents/MacBook/School/Corvinus/TDK/TDK/R_Code/Data.csv")

library(PerformanceAnalytics)
library(rugarch)
library(FinTS)
library(rmgarch)
library(timeSeries) 
library(pastecs)
library(tseries)
library(Hmisc)
library(gogarch)
library(parallel)
library(fpp)
library(texreg)
library(broom)
library(lubridate)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggthemes)
library(ggfortify)
library(bbplot)
library(tidyverse)


# Four digits for numbers
options(digits=4)

# Import data from CSV
data_df = read.csv(file = "/Users/szonyid/Documents/MacBook/School/Corvinus/TDK/TDK/R_Code/Data.csv")

# Figure 1. Time series plots
names = list("US 10-year T-Note", "EUR/USD", "Gold", "WTI Crude Oil", "S&P 500")

Figure_1 <- data_df %>%
  gather(variable, value, S.P.500:EUR.USD) %>%
  group_nest(variable) %>%
  mutate(plot = map2(data, names,
                     ~ { ggplot(.x, aes(x = as.Date(Date), y = value)) +
                         geom_line(color = "dodgerblue4", size = 0.5) +
                         ggtitle(.y) +
                         xlab("") + 
                         ylab("") +
                         scale_x_date(
                           breaks = seq(as.Date("2001-01-01"), 
                                        as.Date("2023-04-01"), 
                                        by = "2 years"),
                           labels = date_format("%b %y")) +
                         scale_colour_canva() + 
                         theme_bw() +
                         theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                               panel.grid.minor = element_blank(),
                               panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
                               axis.text = element_text(colour = "black", face = "bold"))
                        }
                    )
        ) %>%
  walk2(.x = .$variable,
        .y = .$plot,
        .f = ~ ggsave(filename = paste0("I_", .x, ".png"),
                      path = "/Users/szonyid/Documents/MacBook/School/Corvinus/TDK/TDK/Figures",
                      plot = .y)
       )
                  

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
y_ret_sq = as.data.frame(y_ret ^ 2)
data_sq = rownames_to_column(y_ret_sq, "Date")

# Figure 2. Squared daily log-returns

Figure_2 <- data_sq %>%
  gather(variable, value, sp500_ret:eurusd_ret) %>%
  group_nest(variable) %>%
  mutate(plot = map2(data, names,
                     ~ { ggplot(.x, aes(x = as.Date(Date), y = value)) +
                         geom_line(color = "dodgerblue4", size = 0.5) +
                         ggtitle(.y) +
                         xlab("") + 
                         ylab("") +
                         scale_x_date(
                           breaks = seq(as.Date("2001-01-01"), 
                                        as.Date("2023-04-01"), 
                                        by = "2 years"),
                           labels = date_format("%b %y")) +
                         scale_colour_canva() + 
                         theme_bw() +
                         theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                               panel.grid.minor = element_blank(),
                               panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
                               axis.text = element_text(colour = "black", face = "bold"))
                        }
                     )
        ) %>%
  walk2(.x = .$variable,
        .y = .$plot,
        .f = ~ ggsave(filename = paste0("II_", .x, ".png"),
                      path = "/Users/szonyid/Documents/MacBook/School/Corvinus/TDK/TDK/Figures",
                      plot = .y)
        )

# Specify the number of columns for functions
k = ncol(y_ret)

# Table 1. Descriptive statistics
path_to_tables = '/Users/szonyid/Documents/MacBook/School/Corvinus/TDK/TDK/Tables'

table_1 <- data.frame(SummaryStatistics(y_ret))
write.csv(table_1, file.path(path_to_tables, 'Table_1.csv'), row.names = TRUE)

# Table 2. Pearson correlations between daily log-returns 
table_2 <- format(round(data.frame(cor(y_ret)), 3), nsmall = 2)
write.csv(table_2, file.path(path_to_tables, 'Table_2.csv'), row.names = TRUE)
rcorr(y_ret) # Check significance

# Splitting the data at 80% for in-sample model specification

length(y_ret) / 5 * 0.8
y_ret_in_sample = y_ret[1:4451,]

# Table 3. Finding the best model

#################
### DCC MODEL ###
#################

# DCC specification - GARCH(1,1) for conditional correlations (MVT!)

dcc_garch_spec = ugarchspec(mean.model = list(armaOrder = c(1,0)), 
                          variance.model = list(garchOrder = c(1,1), 
                                                model = "sGARCH"), 
                          distribution.model = "std")

dcc_garch = dccspec(uspec = multispec(replicate(5, dcc_garch_spec)), 
                            dccOrder = c(1,1), 
                            distribution = "mvt")

cl = makePSOCKcluster(4)
dcc_garch_fit = dccfit(dcc_garch, data = y_ret_in_sample, cluster = cl)
stopCluster(cl)

# Table 4. DCC-GARCH parameter estimates
table_4 <- dcc_garch_fit # Copying manually as text, then CSV then LaTeX

# Extracting correlation series
dcc_rho12 = rcor(dcc_garch_fit)[1,2,]
dcc_rho13 = rcor(dcc_garch_fit)[1,3,]
dcc_rho14 = rcor(dcc_garch_fit)[1,4,]
dcc_rho15 = rcor(dcc_garch_fit)[1,5,]

dcc_correlations = merge(as.zoo(dcc_rho12), 
                         as.zoo(dcc_rho13),
                         as.zoo(dcc_rho14),
                         as.zoo(dcc_rho15))

index(dcc_correlations) = index(y_ret)
dcc_correlations_df = rownames_to_column(as.data.frame(dcc_correlations), "Date")
colnames(dcc_correlations_df) <- c("Date", "dcc_bond", "dcc_oil", "dcc_gold", "dcc_eurusd")

#########################
### ROLLING DCC-MODEL ###
#########################

# Set number of model refits
r.e = 20

cl = makePSOCKcluster(4)

dcc_garch_roll = dccroll(dcc_garch, y_ret, n.ahead = 1,
                         n.start = 4452,
                         refit.every= r.e, refit.window = "moving", cluster = cl)

stopCluster(cl) 

# Dynamic conditional correlations
corr_dcc_12 = rcor(dcc_garch_roll)[1,2,]
corr_dcc_13 = rcor(dcc_garch_roll)[1,3,]
corr_dcc_14 = rcor(dcc_garch_roll)[1,4,]
corr_dcc_15 = rcor(dcc_garch_roll)[1,5,]
 
# Hedge ratios
hedge_dcc_12 = rcov(dcc_garch_roll)[1,2,] / rcov(dcc_garch_roll)[2,2,]
hedge_dcc_13 = rcov(dcc_garch_roll)[1,3,] / rcov(dcc_garch_roll)[3,3,]
hedge_dcc_14 = rcov(dcc_garch_roll)[1,4,] / rcov(dcc_garch_roll)[4,4,]
hedge_dcc_15 = rcov(dcc_garch_roll)[1,5,] / rcov(dcc_garch_roll)[5,5,]

pw_dcc_12 = (rcov(dcc_garch_roll)[2,2,] - rcov(dcc_garch_roll)[1,2,]) / (rcov(dcc_garch_roll)[1,1,] - 2 * rcov(dcc_garch_roll)[1,2,] + rcov(dcc_garch_roll)[2,2,])
pw_dcc_13 = (rcov(dcc_garch_roll)[3,3,] - rcov(dcc_garch_roll)[1,3,]) / (rcov(dcc_garch_roll)[1,1,] - 2 * rcov(dcc_garch_roll)[1,3,] + rcov(dcc_garch_roll)[3,3,])  
pw_dcc_14 = (rcov(dcc_garch_roll)[4,4,] - rcov(dcc_garch_roll)[1,4,]) / (rcov(dcc_garch_roll)[1,1,] - 2 * rcov(dcc_garch_roll)[1,4,] + rcov(dcc_garch_roll)[4,4,])
pw_dcc_15 = (rcov(dcc_garch_roll)[5,5,] - rcov(dcc_garch_roll)[1,5,]) / (rcov(dcc_garch_roll)[1,1,] - 2 * rcov(dcc_garch_roll)[1,5,] + rcov(dcc_garch_roll)[5,5,])

pw_dcc_12 = ifelse(pw_dcc_12 > 1, 1, pw_dcc_12)
pw_dcc_13 = ifelse(pw_dcc_13 > 1, 1, pw_dcc_13)
pw_dcc_14 = ifelse(pw_dcc_14 > 1, 1, pw_dcc_14)
pw_dcc_15 = ifelse(pw_dcc_15 > 1, 1, pw_dcc_15)

pw_dcc_12 = ifelse(pw_dcc_12 < 0, 0, pw_dcc_12)
pw_dcc_13 = ifelse(pw_dcc_13 < 0, 0, pw_dcc_13)
pw_dcc_14 = ifelse(pw_dcc_14 < 0, 0, pw_dcc_14)
pw_dcc_15 = ifelse(pw_dcc_15 < 0, 0, pw_dcc_15)

dcc_dynamic_correlations = merge(as.zoo(corr_dcc_12), 
                          as.zoo(corr_dcc_13),
                          as.zoo(corr_dcc_14),
                          as.zoo(corr_dcc_15))

dcc_optimal_hedge = merge(as.zoo(hedge_dcc_12), 
                         as.zoo(hedge_dcc_13),
                         as.zoo(hedge_dcc_14),
                         as.zoo(hedge_dcc_15))

dcc_portfolio_weights = merge(as.zoo(pw_dcc_12), 
                              as.zoo(pw_dcc_13),
                              as.zoo(pw_dcc_14),
                              as.zoo(pw_dcc_15))

full_data_total_length <- length(y_ret) / ncol(y_ret)

length_split <- full_data_total_length - 1112 # out-of-sample observations

index(dcc_dynamic_correlations) = index(y_ret)[-c(1:length_split)]
dcc_dynamic_correlations_df = rownames_to_column(as.data.frame(dcc_dynamic_correlations), "Date")
colnames(dcc_dynamic_correlations_df) <- c("Date", "dcc_corr_bond", "dcc_corr_oil", "dcc_corr_gold", "dcc_corr_eurusd")

index(dcc_optimal_hedge) = index(y_ret)[-c(1:length_split)]
dcc_optimal_hedge_df = rownames_to_column(as.data.frame(dcc_optimal_hedge), "Date")
colnames(dcc_optimal_hedge_df) <- c("Date", "dcc_hedge_bond", "dcc_hedge_oil", "dcc_hedge_gold", "dcc_hedge_eurusd")

index(dcc_portfolio_weights) = index(y_ret)[-c(1:length_split)]
dcc_portfolio_weights_df = rownames_to_column(as.data.frame(dcc_portfolio_weights), "Date")
colnames(dcc_portfolio_weights_df) <- c("Date", "dcc_pw_bond", "dcc_pw_oil", "dcc_pw_gold", "dcc_pw_eurusd")



##################
### ADCC MODEL ###
##################

# aDCC specification - GARCH(1,1) for conditional correlations (MVT!)

adcc_garch_spec = ugarchspec(mean.model = list(armaOrder = c(1,0)), 
                          variance.model = list(garchOrder = c(1,1), 
                                                model = "gjrGARCH"), 
                          distribution.model = "std")

adcc_garch = dccspec(uspec = multispec( replicate(5, adcc_garch_spec)), 
                            dccOrder = c(1,1), 
                            distribution = "mvt", model="aDCC")

cl = makePSOCKcluster(10)
adcc_garch_fit = dccfit(adcc_garch, data = y_ret, cluster = cl)
stopCluster(cl)

# Table 3. DCC-GARCH parameter estimates
table_3_2 <- adcc_garch_fit # Copying manually as text, then csv then LaTeX

# Extracting correlation series

adcc_rho12 = rcor(adcc_garch_fit)[1,2,]
adcc_rho13 = rcor(adcc_garch_fit)[1,3,]
adcc_rho14 = rcor(adcc_garch_fit)[1,4,]
adcc_rho15 = rcor(adcc_garch_fit)[1,5,]

par(mfrow=c(2,2))
plot(timeSeries(adcc_rho12), main="aDCC between S&P 500 / Bond", ylab="") 
plot(timeSeries(adcc_rho13), main="aDCC between S&P 500 / Oil", ylab="") 
plot(timeSeries(adcc_rho14), main="aDCC between S&P 500 / Gold", ylab="") 
plot(timeSeries(adcc_rho15), main="aDCC between S&P 500 / EUR/USD", ylab="") 
par(mfrow=c(1,1))

##########################
### ROLLING aDCC-MODEL ###
##########################

cl = makePSOCKcluster(10)

adcc_garch_roll = dccroll(adcc_garch, y_ret, n.ahead=1,
                            forecast.length = f.l,
                            refit.every = r.e, refit.window = "moving", cluster = cl)

stopCluster(cl) 

# Hedge ratios
hedge_adcc_12 = rcov(adcc_garch_roll)[1,2,] / rcov(adcc_garch_roll)[2,2,]
hedge_adcc_13 = rcov(adcc_garch_roll)[1,3,] / rcov(adcc_garch_roll)[3,3,]
hedge_adcc_14 = rcov(adcc_garch_roll)[1,4,] / rcov(adcc_garch_roll)[4,4,]
hedge_adcc_15 = rcov(adcc_garch_roll)[1,5,] / rcov(adcc_garch_roll)[5,5,]

# Optimal hedge ratio plots (aDCC-GARCH)
par(mfrow=c(2,2))
plot(timeSeries(hedge_adcc_12), main="Optimal hedge ratio: S&P 500 / Bond", ylab="")  
plot(timeSeries(hedge_adcc_13), main="Optimal hedge ratio: S&P 500 / Oil", ylab="")  
plot(timeSeries(hedge_adcc_14), main="Optimal hedge ratio: S&P 500 / Gold", ylab="")  
plot(timeSeries(hedge_adcc_15), main="Optimal hedge ratio: S&P 500 / EUR/USD", ylab="")  
par(mfrow=c(1,1))

# Hedge performance 
adcc_port_12 = y_ret.l[,1] - hedge_adcc_12 * y_ret.l[,2]
adcc_he_12 = ( var(y_ret.l[,1]) - var(adcc_port_12) ) / var(y_ret.l[,1])

adcc_port_13 = y_ret.l[,1] - hedge_adcc_13 * y_ret.l[,3]
adcc_he_13 = ( var(y_ret.l[,1]) - var(adcc_port_13) ) / var(y_ret.l[,1])

adcc_port_14 = y_ret.l[,1] - hedge_adcc_14 * y_ret.l[,4]
adcc_he_14 = ( var(y_ret.l[,1]) - var(adcc_port_14) ) / var(y_ret.l[,1])

adcc_port_15 = y_ret.l[,1] - hedge_adcc_15 * y_ret.l[,5]
adcc_he_15 = ( var(y_ret.l[,1]) - var(adcc_port_15) ) / var(y_ret.l[,1])

# Matrices
mem_bond[2,1] = mean(hedge_adcc_12)
mem_bond[2,2] = min(hedge_adcc_12)
mem_bond[2,3] = max(hedge_adcc_12)
mem_bond[2,4] = adcc_he_12

mem_oil[2,1] = mean(hedge_adcc_13)
mem_oil[2,2] = min(hedge_adcc_13)
mem_oil[2,3] = max(hedge_adcc_13)
mem_oil[2,4] = adcc_he_13

mem_gold[2,1] = mean(hedge_adcc_14)
mem_gold[2,2] = min(hedge_adcc_14)
mem_gold[2,3] = max(hedge_adcc_14)
mem_gold[2,4] = adcc_he_14

mem_eurusd[2,1] = mean(hedge_adcc_15)
mem_eurusd[2,2] = min(hedge_adcc_15)
mem_eurusd[2,3] = max(hedge_adcc_15)
mem_eurusd[2,4] = adcc_he_15

######################
### GO-GARCH MODEL ###
######################



# GO specification - GARCH(1,1) for conditional correlations (MANIG!)


go_garch_spec = gogarchspec(mean.model = list( model='AR',lag=0), variance.model = list(model = "sGARCH", garchOrder = c(1, 1),
          variance.targeting = FALSE), distribution = 'manig', ica = 'fastica')

cl = makePSOCKcluster(4)
go_garch_fit = gogarchfit(go_garch_spec, y_ret_in_sample, gfun = 'tanh', rseed = 77, cluster = cl)
stopCluster(cl)

go_rho12 = rcor(go_garch_fit)[1,2,]
go_rho13 = rcor(go_garch_fit)[1,3,]
go_rho14 = rcor(go_garch_fit)[1,4,]
go_rho15 = rcor(go_garch_fit)[1,5,]

go_correlations = merge(as.zoo(go_rho12), 
                         as.zoo(go_rho13),
                         as.zoo(go_rho14),
                         as.zoo(go_rho15))

index(go_correlations) = index(y_ret_in_sample)
go_correlations_df = rownames_to_column(as.data.frame(go_correlations), "Date")
colnames(go_correlations_df) <- c("Date", "go_bond", "go_oil", "go_gold", "go_eurusd")

# Table 5. GO-GARCH parameter estimates
table_5_1 <- go_garch_fit
table_5_2 <- coef(go_garch_fit) 

# Figure 4_1. News impact correlation surfaces (DCC-GARCH)
ni_dcc = nisurface(dcc_garch_fit, type = "cor", pair = c(1, 2),  plot = TRUE)

# Figure 5_1. News impact correlation surfaces (GO-GARCH)
ni_go = nisurface(go_garch_fit, type = "cor", pair = c(1, 2), plot = TRUE)

##############################
### ROLLING GO-GARCH MODEL ###
##############################

cl = makePSOCKcluster(4)

go_garch_roll = gogarchroll(go_garch_spec, y_ret, n.ahead = 1,
                            n.start = 4452,
                            refit.every = r.e, refit.window = "moving", rseed = 77, cluster = cl)

stopCluster(cl) 
 
# Dynamic conditional correlations
corr_go_12 = rcor(go_garch_roll)[1,2,]
corr_go_13 = rcor(go_garch_roll)[1,3,]
corr_go_14 = rcor(go_garch_roll)[1,4,]
corr_go_15 = rcor(go_garch_roll)[1,5,]

# Hedge ratios  
hedge_go_12 = rcov(go_garch_roll)[1,2,] / rcov(go_garch_roll)[2,2,]
hedge_go_13 = rcov(go_garch_roll)[1,3,] / rcov(go_garch_roll)[3,3,]
hedge_go_14 = rcov(go_garch_roll)[1,4,] / rcov(go_garch_roll)[4,4,]
hedge_go_15 = rcov(go_garch_roll)[1,5,] / rcov(go_garch_roll)[5,5,]

pw_go_12 = (rcov(go_garch_roll)[2,2,] - rcov(go_garch_roll)[1,2,]) / (rcov(go_garch_roll)[1,1,] - 2 * rcov(go_garch_roll)[1,2,] + rcov(go_garch_roll)[2,2,])
pw_go_13 = (rcov(go_garch_roll)[3,3,] - rcov(go_garch_roll)[1,3,]) / (rcov(go_garch_roll)[1,1,] - 2 * rcov(go_garch_roll)[1,3,] + rcov(go_garch_roll)[3,3,])  
pw_go_14 = (rcov(go_garch_roll)[4,4,] - rcov(go_garch_roll)[1,4,]) / (rcov(go_garch_roll)[1,1,] - 2 * rcov(go_garch_roll)[1,4,] + rcov(go_garch_roll)[4,4,])
pw_go_15 = (rcov(go_garch_roll)[5,5,] - rcov(go_garch_roll)[1,5,]) / (rcov(go_garch_roll)[1,1,] - 2 * rcov(go_garch_roll)[1,5,] + rcov(go_garch_roll)[5,5,])

pw_go_12 = ifelse(pw_go_12 > 1, 1, pw_go_12)
pw_go_13 = ifelse(pw_go_13 > 1, 1, pw_go_13)
pw_go_14 = ifelse(pw_go_14 > 1, 1, pw_go_14)
pw_go_15 = ifelse(pw_go_15 > 1, 1, pw_go_15)

pw_go_12 = ifelse(pw_go_12 < 0, 0, pw_go_12)
pw_go_13 = ifelse(pw_go_13 < 0, 0, pw_go_13)
pw_go_14 = ifelse(pw_go_14 < 0, 0, pw_go_14)
pw_go_15 = ifelse(pw_go_15 < 0, 0, pw_go_15)

go_dynamic_correlations = merge(as.zoo(corr_go_12), 
                                 as.zoo(corr_go_13),
                                 as.zoo(corr_go_14),
                                 as.zoo(corr_go_15))

go_optimal_hedge = merge(as.zoo(hedge_go_12), 
                          as.zoo(hedge_go_13),
                          as.zoo(hedge_go_14),
                          as.zoo(hedge_go_15))

go_portfolio_weights = merge(as.zoo(pw_go_12), 
                              as.zoo(pw_go_13),
                              as.zoo(pw_go_14),
                              as.zoo(pw_go_15))

index(go_dynamic_correlations) = index(y_ret)[-c(1:length_split)]
go_dynamic_correlations_df = rownames_to_column(as.data.frame(go_dynamic_correlations), "Date")
colnames(go_dynamic_correlations_df) <- c("Date", "go_corr_bond", "go_corr_oil", "go_corr_gold", "go_corr_eurusd")

index(go_optimal_hedge) = index(y_ret)[-c(1:length_split)]
go_optimal_hedge_df = rownames_to_column(as.data.frame(go_optimal_hedge), "Date")
colnames(go_optimal_hedge_df) <- c("Date", "go_hedge_bond", "go_hedge_oil", "go_hedge_gold", "go_hedge_eurusd")

index(go_portfolio_weights) = index(y_ret)[-c(1:length_split)]
go_portfolio_weights_df = rownames_to_column(as.data.frame(go_portfolio_weights), "Date")
colnames(go_portfolio_weights_df) <- c("Date", "go_pw_bond", "go_pw_oil", "go_pw_gold", "go_pw_eurusd")



####################
### PLOT FACTORY ###
####################



# Plot dynamic conditional correlations manually (DCC - GO-GARCH)

correlations_df = merge(dcc_dynamic_correlations_df, go_dynamic_correlations_df)

ggplot(correlations_df, 
       aes(x = as.Date(Date))) +
  geom_line(aes(y = dcc_corr_bond, color = "DCC-GARCH"), size = 0.6) +
  geom_line(aes(y = go_corr_bond, color = "GO-GARCH"), size = 0.6) +
  ggtitle("S&P 500 / Bond") +
  xlab("") + 
  ylab("") +
  scale_x_date(
    breaks = seq(as.Date("2018-01-01"), 
                 as.Date("2023-04-01"), 
                 by = "6 months"),
    labels = date_format("%b %y")) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.2)) +
  scale_colour_manual("", 
                      breaks = c("DCC-GARCH", "GO-GARCH"),
                      values = c("seagreen", "firebrick1")) +
  #scale_colour_canva() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.key.size = unit(3, 'cm'),
        legend.key.width= unit(1, 'cm'),
        legend.text = element_text(size = 9, face = "bold"))


# Plot dynamic hedging ratios manually (DCC - GO-GARCH)


hedging_df = merge(dcc_optimal_hedge_df, go_optimal_hedge_df)

ggplot(hedging_df, 
       aes(x = as.Date(Date))) +
  geom_line(aes(y = dcc_hedge_bond, color = "DCC-GARCH"), size = 0.6) +
  geom_line(aes(y = go_hedge_bond, color = "GO-GARCH"), size = 0.6) +
  ggtitle("S&P 500 / Bond") +
  xlab("") + 
  ylab("") +
  scale_x_date(
    breaks = seq(as.Date("2018-01-01"), 
                 as.Date("2023-04-01"), 
                 by = "6 months"),
    labels = date_format("%b %y")) +
  #scale_y_continuous(breaks = seq(-3, 5, by = 1)) +
  scale_colour_manual("", 
                      breaks = c("DCC-GARCH", "GO-GARCH"),
                      values = c("seagreen", "firebrick1")) +
  # scale_colour_canva() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.key.size = unit(3, 'cm'),
        legend.key.width= unit(1, 'cm'),
        legend.text = element_text(size = 9, face = "bold"))


# Plot dynamic portfolio weights manually (DCC - GO-GARCH)


pw_df = merge(dcc_portfolio_weights_df, go_portfolio_weights_df)

ggplot(pw_df, 
       aes(x = as.Date(Date))) +
  geom_line(aes(y = dcc_pw_eurusd, color = "DCC-GARCH"), size = 0.6) +
  geom_line(aes(y = go_pw_eurusd, color = "GO-GARCH"), size = 0.6) +
  ggtitle("S&P 500 / EUR/USD") +
  xlab("") + 
  ylab("") +
  scale_x_date(
    breaks = seq(as.Date("2018-01-01"), 
                 as.Date("2022-01-01"), 
                 by = "6 months"),
    labels = date_format("%b %y")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_colour_manual("", 
                      breaks = c("DCC-GARCH", "GO-GARCH"),
                      values = c("seagreen", "firebrick1")) +
  # scale_colour_canva() + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey85", linetype = "dashed"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.key.size = unit(3, 'cm'),
        legend.key.width= unit(1, 'cm'),
        legend.text = element_text(size = 9, face = "bold"))


# Table 6. Correlations between dynamic conditional correlations


corr_model = matrix(1:4, ncol=4, byrow=TRUE )

k = 1
for (j in 2:5){
corr_model[1,k] = (cor(rcor(dcc_garch_roll)[1,j,], rcor(go_garch_roll)[1,j,]) )
k = k+1
}

table_6 <- corr_model
write.csv(table_6, file.path(path_to_tables, 'Table_6.csv'), row.names = TRUE)


# Table 7. Correlations between dynamic hedge ratios


corr_hedge = matrix(1:4, ncol=4, byrow=TRUE )

corrhedge = cbind(hedge_dcc_12,hedge_go_12,
                  hedge_dcc_13,hedge_go_13,
                  hedge_dcc_14,hedge_go_14,
                  hedge_dcc_15,hedge_go_15)
                 
k = 1
for (j in seq(1,7,2)){
  corr_hedge[1,k] = cor(corrhedge[,j], corrhedge[,j + 1])
  k = k + 1
}

table_7 <- corr_hedge
write.csv(table_7, file.path(path_to_tables, 'Table_7.csv'), row.names = TRUE)


# Table 8. Correlations between dynamic portfolio weights


corr_pw = matrix(1:4, ncol=4, byrow=TRUE )

corrpw = cbind(pw_dcc_12,pw_go_12,
               pw_dcc_13,pw_go_13,
               pw_dcc_14,pw_go_14,
               pw_dcc_15,pw_go_15)

k = 1
for (j in seq(1,7,2)){
  corr_pw[1,k] = cor(corrpw[,j], corrpw[,j + 1])
  k = k + 1
}

table_8 <- corr_pw
write.csv(table_8, file.path(path_to_tables, 'Table_8.csv'), row.names = TRUE)