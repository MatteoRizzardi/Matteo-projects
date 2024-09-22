In the Option Project you can see a work I have created during my master's degree in Quantitative Finance. In there, you can find two R files: one is for the functions, one is the main.
The excel file is the original dataset, and the pdf file is a comment to the work. Please download the entire foulder and run only the main file, everything will work.

Here is what I had to do:

I was given a dataset of 58 call options and I had to:

- Clean it up using log-moneyness, Merton's, monotonicity and convexity constraints.
- Construct the volatility smile using both the Newton-Raphson algorithm (constructed manually by us) and the R function "Optim".
- Using the volatility nearest at-the-money level, compare the behaviour of the approximated Gamma (using finite difference method) with the exact Gamma for varying values of the underlying asset S0.
- Calibrate the volatility for B&S using MSE, and price a Bear Strategy with it.
- Using the abovementioned calibrated sigma, compute the Monte Carlo price for a Floating Arithmentic Call Asian Option (payoff in the paper). We used the Euler simulation scheme for a Geometric Brownian Motion to simulate the behavior of the underlying asset S0.
- Analyse the behavior of the Monte Carlo price as the number of simulations increases.

This work was a success: we wrote about 700 lines of code, which was of great help for our Rskills.
