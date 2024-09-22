######################
#     Point 1        #
######################

#Function to transform interest rate from yearly to daily basis
daily_int <- function(dataset){
  
  dataset$daily_interest <- (1 + dataset$interest_rate)^(1/252)-1
  
  return(dataset)
}

#Function to obtain a bound on the moneyness
ds_clean_logmoneyness <- function(dataset, low_lim, up_lim){
  
  dataset$log_moneyness <- rep(0, nrow(dataset))
  
  for(i in 1:nrow(dataset)){
    
    dataset$log_moneyness[i]<-log(dataset$underlying_price[i]/dataset$strike[i])
  }
  
  new_dataset <- subset(dataset,dataset$log_moneyness >= low_lim &
                         dataset$log_moneyness <= up_lim)
  
  return(new_dataset)
}

#Function to control and apply Merton's constraints
Merton <- function(dataset){
  
  cond_call <- dataset$underlying_price >= dataset$option_price &
    dataset$option_price >= pmax(dataset$underlying_price - 
  dataset$strike * exp(-dataset$daily_interest*dataset$ttm),0)
  
  merton_data <- dataset[cond_call,] #keeping rows that satisfy the above 
  #condition
  nec <- abs(dim(dataset[1]) - dim(merton_data[1])) #number of eliminated calls
  cat("Number of eliminated calls with Merton:", nec[1], fill = TRUE)
  
  return(merton_data)
}




#######################
#       Point 2       #
#######################

#Black and Scholes function: application of the Black and Scholes' formula.
#It apperars in the NR_implied_vol function for the first time.

option_priceBS <- function(S0, K, r, ttm, sigma) {
  
  d1 <- (log(S0/K) + (r+0.5*sigma^2)*ttm)/(sigma*sqrt(ttm))
  d2 <- d1 - sigma*sqrt(ttm)
  
  call_price <- S0*pnorm(d1) - K*exp(-r*ttm)*pnorm(d2)
  
  return(call_price)
}


#Vega function: application of the vega formula.
#It appears in the NR_implied_vol function for the first time.

vega <- function(S0, K, r, ttm, sigma) {
  
  d1 <- (log(S0/K) + (r+0.5*sigma^2)*ttm)/(sigma*sqrt(ttm))
  
  vega_value <- S0 * dnorm(d1) * sqrt(ttm)
  return(vega_value)
  
}


#Newton-Raphson method for implied volatility: this function replicates the
#Newton-Raphson method, applying the formula x_i = x_(i-1) - f_old/df_old.
#f_old is calculated with option_priceBS function, while df_old with vega 
#function.

NR_implied_vol <- function(S0, K, r, ttm, mkt_price, x){
  
  f_old <- option_priceBS(S0, K, r, ttm, x) - mkt_price
  cond <- 10^-5
  iter <- 0
  
  while (abs(f_old) >= cond) {
    
    f_old <- option_priceBS(S0, K, r, ttm, x) - mkt_price
    df_old <- vega(S0=S0, K=K, r=r, ttm=ttm, sigma=x)
    
    if (is.nan(f_old) || is.nan(df_old)) {
      return(NA)  #return NA if NaN or near-zero vega
    }
    x <- x - f_old/df_old
    iter <- iter+1
  }
  return(list("implied_vol" = x, "number of iterations" = iter))
}


#Optim for implied volatility: we need two functions, BS_mkt_MSE calculates
#the Mean Squared Error between the BS price and the option market price. In the
#second function the former has been called.

BS_mkt_MSE <- function(S0, K, r, ttm, par, mkt_price){
  
  res <- (option_priceBS(S0 = S0, K = K, r = r, ttm = ttm, sigma = par)
          - mkt_price)^2 #MSE between B&S price and mkt price.
  
  return(mean(res))
}

implied_vol_optim <- function(S0, K, r, ttm, par, mkt_price){
  
  par_vec <- rep(0, length(mkt_price))#initialize vector of implied volatilities
  
  for(i in c(1:length(mkt_price))){
    
    dummy <- optim(par = par, fn = BS_mkt_MSE, lower = 0.01, upper = 1,
                   method ="L-BFGS-B", S0 = S0[i], K = K[i],
                   ttm = ttm[i], r = r[i], mkt_price = mkt_price[i])
    
    #Numerical computation brings us to choose "lower = 0.01" instead of 0.001.
    #We think this is the case because, in the latter case, we allow sigma (par)
    #to be very small, and many optimal solutions would fall in such value
    #(sigma = 0.001), which would be unreasonable.
    #By restricting the possible values of sigma between 0.01 and 1, we find
    #solution coinciding with the ones found with Newton-Raphson method.
    
    #Printing the convergence status of the function for each iteration.
    if(dummy$convergence==0){
      cat("\n", "okay") #succesful computation of optim
    }else{
      cat("\n objective function ", c(i, dummy$convergence)) #convergence error
    }
    
    par_vec[i] <- dummy$par
  }
  return(par_vec)
}

#Function to check the MSE between Optim and Newton-Raphson methods

MSE <- function(dataset1, dataset2){
  
  mse <- mean((as.numeric(dataset1)- as.numeric(dataset2))^2)
  
  return(mse)
}


#Nearest at the money option function: the "most at the money"
#option's data are required in the 2.2 point. This function takes them.

atm_opt <- function(dataset) {
  
  min_SE <- NA
  
  for (i in 1:nrow(dataset)){
    
    SE <- (dataset$underlying_price[i] - dataset$strike[i])^2
    
    if(is.na(min_SE)) {
      
      min_SE <- SE
    }else if (SE < min_SE){
      
      min_SE <- SE
      atm_opt <- dataset[i,]
    }
  }
  return(atm_opt)
}


#Function for approximated gamma with (2nd order) central difference method
#gamma = lim h -> 0: (f(x+h) - 2f(x) + (f(x-h)))/h^2

approx_gamma_central <- function(S0, K, r, ttm, sigma, h) {
  
  gamma <- (option_priceBS(S0+h, K, r, ttm, sigma) - 
            2*(option_priceBS(S0, K, r, ttm, sigma)) +
            (option_priceBS(S0-h, K, r, ttm, sigma)))/(h^2)
  return(gamma)
}




#######################
#       Point 3       #
#######################

#Calibration of sigma function: this function returns the volatility that 
#minimizes the MSE between the B&S price and the Market price of the options.

sigma_calibration <- function(S0, K, r, ttm, sigma_par, opt_price_mkt) {
    
    sigma_ok <- optim(par = sigma_par, fn = BS_mkt_MSE, lower = 0.01, upper = 1,
                      
                      method ="L-BFGS-B", S0 = S0, K = K,
                      
                      ttm = ttm, r = r, mkt_price = opt_price_mkt)
    
    return(sigma_ok$par)
    
}


#Bear strategy pricing function:
#In a no-arbitrage context, for call options, if K1 is lower than K2, by selling
#call(K1) and buying call(K2), the price of the strategy will be negative,
#meaning that we receive money to hold the portfolio.

Bear_price <- function(S0, K1, K2, r, ttm, sigma){
  
  call_1 <- option_priceBS(S0, K1, r, ttm, sigma)
  
  call_2 <- option_priceBS(S0, K2, r, ttm, sigma)
  
  price <- call_2 - call_1 #buying call_2 and selling call_1
  
  #call_1 is more expensive than call_2, because it has a lower strike, thus
  #we will receive money at time t0 to hold the portfolio.
  
  return(price)
}



#######################
#       Point 4       #
#######################


#Function to estimate the price of a Floating Arithmetic Asian Call Option 
#using Monte Carlo simulation with Euler discretization of a Geometric Brownian 
#Motion (GBM) to simulate the underlying assets' prices.
#Important note: 
#N: number of time steps.
#nsim: number of simulations for each option price.

MC_price_GBM <- function(S0, r, mu, sigma, ttm, t0, N, nsim) {
  
  und_price <- matrix(NA, nrow = nsim, ncol = N+1) #initializing the matrix
  #containing the evolution of the underlying asset.
  
  und_price[,1] <- S0 #at time t0, the price is given.
  
  delta_t <- (ttm-t0)/N #time-step is 5 days.
  
  #Euler discretization of a GBM.
  for(i in c(2:(N+1))){
    
    und_price[,i] <- und_price[,i-1] + mu*und_price[,i-1] * delta_t +
      sigma * sqrt(delta_t) * und_price[,i-1] * rnorm(nsim)
    
  }
  
  ar_mean <- rowMeans(und_price) #taking the mean of the simulated prices
  
  sub <- ar_mean - und_price[,N+1] #computing the difference between
  #the means and the final value ST.
  
  final_payoffs <- pmax(sub, 0) #taking the maximum between 0 and the diff.
  
  sim_prices <- exp(-r*(ttm-t0)) * final_payoffs #discounted payoffs.
  
  mc_price <- mean(sim_prices) #MonteCarlo: approximating the Expected value
  #under Q with the function "mean".
  
  LB <- mc_price + qnorm(0.025)*sqrt(var(sim_prices)/nsim)
  UB <- mc_price + qnorm(0.975)*sqrt(var(sim_prices)/nsim)
  
  res <- c(LB, mc_price, UB)
  names(res) <- c("Lower Bound", "Montecarlo Price", "Upper Bound")
  
  return(res)
  
}


#Function to perform multiple simulations of Monte Carlo pricing for
#our option. This function improves the robustness of the analysis:
#we compute the Monte Carlo prices M times and then we take the median result.
#Important note:
#M: Number of simulations to perform.

M_simulations_MC_prices <- function(S0, r, mu, sigma, ttm, t0, N, nsim, M){
  
  matr_LB <- matrix(NA, nrow = 1, ncol = M)
  matr_MC <- matrix(NA, nrow = 1, ncol = M)
  matr_UB <- matrix(NA, nrow = 1, ncol = M) #initialize matrices to store
  #lower bounds, Monte Carlo prices, and upper bounds for each simulation.
  
  #here we perform M simulations
  for(i in 1:M){
    
    Option_price_behavior <- as.matrix(MC_price_GBM(S0, r, mu, sigma, ttm, 
                                                    t0, N, nsim))
    #call MC_price_GBM function to obtain option price behavior for 
    #each simulation.
    
    matr_LB[,i] <- Option_price_behavior[1,]
    matr_MC[,i] <- Option_price_behavior[2,]
    matr_UB[,i] <- Option_price_behavior[3,] #store lower bounds, 
    #Monte Carlo prices, and upper bounds for each simulation.
  }
  
  
  median_mat <- matrix(NA, nrow = 1, ncol = 3)
  median_mat[1] <- quantile(matr_LB, probs = 0.5)
  median_mat[2] <- quantile(matr_MC, probs = 0.5)
  median_mat[3] <- quantile(matr_UB, probs = 0.5) #calculate the median of lower
  #bounds, Monte Carlo prices, and upper bounds across all simulations.
  
  rownames(median_mat) <- paste("nsim: ", nsim)
  colnames(median_mat) <- c("Lower Bound", "Montecarlo Price", "Upper Bound")  
  return(median_mat)
}
  

#Function to analyze Monte Carlo behavior for different numbers of simulations
#of underlying prices.
#Important note:
#nsim_vec: Vector containing different numbers of simulations (S_t).

MC_behavior <- function(S0, r, mu, sigma, ttm, t0, N, nsim_vec){
  
  Mc_price <- matrix(NA, nrow = length(nsim_vec), ncol = 3)
  #initialize a matrix to store Monte Carlo prices for 
  #different numbers of simulations.
  
  #loop through each number of simulations.
  for (i in 1:(length(nsim_vec))){
    
    Mc_price[i,] <- MC_price_GBM(S0, r, mu, sigma, ttm, t0, N, nsim_vec[i])
    #call MC_price_GBM function to obtain Monte Carlo prices for 
    #the current number of simulations.

  }
  
  #return the matrix containing Monte Carlo prices for different 
  #numbers of simulations.

  return(Mc_price)
  
}


#Function to perform multiple simulations of Monte Carlo pricing for our
#option and plot the results.

M_simulations_MC_behaviour_plot <- function(S0, r, mu, sigma, ttm, t0, N, 
                                            nsim_vec, M){
  
  matr_LB <- matrix(NA, nrow = length(nsim_vec), ncol = M)
  matr_MC <- matrix(NA, nrow = length(nsim_vec), ncol = M)
  matr_UB <- matrix(NA, nrow = length(nsim_vec), ncol = M)#initialize matrices 
  #to store lower bounds, Monte Carlo prices, and upper bounds for 
  #each simulation.
  
  for(i in 1:M){
    
    Option_price_behavior <- MC_behavior(S0, r, mu, sigma, ttm, t0, N, nsim_vec)
    #call MC_behavior function to obtain option price behavior for 
    #each simulation.
    
    matr_LB[,i] <- Option_price_behavior[,1]
    matr_MC[,i] <- Option_price_behavior[,2]
    matr_UB[,i] <- Option_price_behavior[,3]#store lower bounds, Monte Carlo 
    #prices, and upper bounds for each simulation.
  }
  
  median_mat <- matrix(NA, nrow = length(nsim_vec), ncol = 3)#initialization 
  
  for (i in 1:nrow(median_mat)){
    
    median_mat[i,1] <- quantile(matr_LB[i,], probs = 0.5)
    median_mat[i,2] <- quantile(matr_MC[i,], probs = 0.5)
    median_mat[i,3] <- quantile(matr_UB[i,], probs = 0.5) #calculate the 
    #median of lower bounds, Monte Carlo prices, and upper bounds across 
    #all simulations.
    
  }
  
  #plot the results with average MC price for each nsim, the confidence 
  #intervals, the outliers and the boxplots containing all the MC prices.
  plot(nsim_vec, median_mat[,2], ylim = c(-1,30),
       xlab = "Number of Simulations", ylab = "MonteCarlo Price", col = "red",
       lwd = 2, main = "Behavior of MC Price as simulations increase")
  lines(nsim_vec, median_mat[,1], col = "blue")
  lines(nsim_vec, median_mat[,3], col = "blue")
  boxplot(t(matr_MC), add = TRUE, at = nsim_vec, boxwex = 15, col = "orange",
          axes = F)
  legend("topright", legend = c("Median Monte Carlo prices",
                                "Confidence Intervals",
                                "Monte Carlo prices"),
         col = c("red", "blue", "orange"), cex = 0.5, pch = 16)
  
  
  #at the end we return the matrix containing median values of lower bounds, 
  #Monte Carlo prices, and upper bounds.
  rownames(median_mat) <- paste("nsim: ", nsim_vec)
  colnames(median_mat) <- c("Lower Bound", "Montecarlo Price", "Upper Bound")
  return(median_mat)
}
