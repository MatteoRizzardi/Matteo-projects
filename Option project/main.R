######################
#     Point 1        #
######################

#Require readxl package

if (require("readxl") == FALSE) {
  install.packages("readxl", repos = "http://cran.rstudio.com/")
}

library(readxl)

#Call the functions script
source("functions.R")

#Read the raw dataset
raw_data <- read_excel("group_2.xlsx")

#Interest rate adjustment
adj_data <- daily_int(raw_data)

#Clean the dataset using log-moneyness criterion
moneyness_data <- ds_clean_logmoneyness(adj_data, -0.2, 0.2)

#Clean the dataset using Merton's constraints criterion
noarb_data <- Merton(moneyness_data)

#Graphical approach on monotonicity and convexity with respect to strike K
plot(noarb_data$strike, noarb_data$option_price,
  ylab = "Call price", xlab = "Strike price",
  main = "Monotonicity and convexity of C with respect to K")
#the call prices look monotonically non-increasing in K, and they are also a
#convex function of K.

#Final number and percentage of eliminated calls
nec <- abs(nrow(raw_data)) - nrow(noarb_data)
perc <- (abs(nrow(raw_data)) - nrow(noarb_data)) / nrow(raw_data)
eliminated_calls <- list(
  "Total number of eliminated calls" = nec,
  "Percentage of calls eliminated" = perc,
  "Number of remaining calls" = nrow(noarb_data),
  "Percentage of remaining calls" = (1 - perc))
eliminated_calls




#######################
#       Point 2       #
#######################


#############
# Point 2.1 #
#############

#Implied volatility with Newton Raphson

noarb_data$implied_volatility <- rep(0, nrow(noarb_data))

for (i in 1:nrow(noarb_data)) {
  noarb_data$implied_volatility[i] <- NR_implied_vol(
    noarb_data$underlying_price[i],
    noarb_data$strike[i],
    noarb_data$daily_interest[i],
    noarb_data$ttm[i],
    noarb_data$option_price[i],
    0.5
  )[1]
  # initial guess: sigma = 0.5.
}

#Volatility smile Newton-Raphson

plot(noarb_data$strike, noarb_data$implied_volatility,
  xlab = "Strike", ylab = "NR Implied volatility",
  main = "NR Volatility Smile"
)
abline(v = noarb_data$underlying_price[1])
legend("topright",
  legend = "Underlying price", col = "black", lwd = 2,
  cex = 0.5
) #adding the underlying price on the plot with abline

#Implied volatility with optim

noarb_data$optim_implied_volatility <- rep(0, nrow(noarb_data))

noarb_data$optim_implied_volatility <- implied_vol_optim(
  noarb_data$underlying_price,
  noarb_data$strike,
  noarb_data$daily_interest,
  noarb_data$ttm,
  par = 0.5,
  noarb_data$option_price
)
#Error 52 is an error from the "L-BFGS-B" method: it tells us that convergence
#was not achieved within the number of specified iterations.
#However, by comparing the implied volatilities obtained with the Newton-Raphson
#method, we consider these convergence errors negligible (Check "Optim_NR_MSE").

#Volatility smile with optim

plot(noarb_data$strike, noarb_data$optim_implied_volatility,
  xlab = "Strike", ylab = "Optim Implied Volatility",
  main = "Optim Volatility Smile"
)
abline(v = noarb_data$underlying_price[1]) # adding the underlying price on the
legend("topright",
  legend = "Underlying price", col = "black", lwd = 2,
  cex = 0.5
)

#MSE between the implied volatilities computed with the two methods.

Optim_NR_MSE <- MSE(
  noarb_data$implied_volatility,
  noarb_data$optim_implied_volatility
)

Optim_NR_MSE 
#the values of implied volatilities obtained with both
#methods are almost the same. 


#############
# Point 2.2 #
#############

#Initialization

atm_option <- atm_opt(noarb_data)
atm_sigma <- as.numeric(atm_option$implied_volatility)
h <- 0.0001

S_vector <- seq(600, 1000, 1) #vector of the varying underlying prices

#Approximated gamma

approx_gammaC <- approx_gamma_central(S_vector,
  atm_option$strike,
  atm_option$daily_interest,
  atm_option$ttm,
  atm_sigma,
  h = h
)

#Exact gamma

atm_d1 <- (log(S_vector / atm_option$strike) +
  (atm_option$daily_interest + 0.5 * atm_sigma^2) *
    atm_option$ttm) / (atm_sigma * sqrt(atm_option$ttm))

exact_gamma <- dnorm(atm_d1) / (S_vector * atm_sigma
  * sqrt(atm_option$ttm))

#Plotting the behaviour of the gammas, the approximation looks good.

plot(S_vector, exact_gamma,
  xlab = "Underlying price", ylab = "Gamma",
  col = "red", lwd = 3, main = "Gamma behaviour"
)
lines(S_vector, approx_gammaC, col = "blue", type = "p", lwd = 1)
legend("topright",
  legend = c("Exact Gamma", "Approx Gamma"), col = c("red", "blue"),
  cex = 0.5, pch = 16
)




#######################
#       Point 3       #
#######################

calibrated_sigma <- sigma_calibration(
  noarb_data$underlying_price,
  noarb_data$strike,
  noarb_data$daily_interest,
  noarb_data$ttm,
  0.10, #initial guess
  noarb_data$option_price
)

calibrated_sigma

K1 <- noarb_data$underlying_price[1] * 0.975
K2 <- noarb_data$underlying_price[1] * 1.025
ttm_strat <- 60

Bear_strat_price <- Bear_price(
  noarb_data$underlying_price[1],
  K1, K2, noarb_data$daily_interest[1],
  ttm_strat, calibrated_sigma
)

Bear_strat_price 

#We sell the call with the lower strike (K1), and we buy the call with the
#higher strike (K2). In a no-arbitrage context, the call we sell is more
#expensive than the call we buy: thus, we receive money at time t0, for 
#a payoff at capital T that is always lower than or equal to zero.




#######################
#       Point 4       #
#######################

mu <- 0.001 #We arbitrarily choose a small average 5-day rate of return of the
# asset.
ttm_eu <- 30
t0 <- 0
N <- 6
nsim <- 500
M <- 1000

set.seed(123)

#We always take the median of M simulations for a more robust result.

MC_price_interval_with_Mreps <- M_simulations_MC_prices(
  noarb_data$underlying_price[1],
  noarb_data$daily_interest[1],
  mu, calibrated_sigma, ttm_eu, t0, N, 
  nsim, M)

MC_price_interval_with_Mreps

#Analyzing the MonteCarlo price behavior as the number of simulation changes

set.seed(123)

nsim_vec <- c(10, 50, 100, 150, 200, 400, 800, 1000, 2500, 5000)

#Adding M repetitions:
#plotting the MonteCarlo price behavior as the number of simulation changes

mat <- M_simulations_MC_behaviour_plot(
  noarb_data$underlying_price[1],
  noarb_data$daily_interest[1],
  mu, calibrated_sigma, ttm_eu, t0, N,
  nsim_vec, M) # taking the median of M = 1000 LB, MC_price, and UB processes.

mat
