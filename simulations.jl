#######################################Intializations##############################################

W_t  = zeros(Simulations, Steps); # Brownian motion
N_t  = zeros(Simulations, Steps); # Poisson process
X_t  = zeros(Simulations, Steps); # Compound Poisson process
S_t_1  = zeros(Simulations, Steps); # Price of asset 1 based on equation (4.5.3)
S_t_2  = zeros(Simulations, Steps); # Price of asset 2 based on equation (4.5.3)
C_t  = zeros(Simulations, Steps); # Price of the contract based on equation (4.5.3)

##Optimal hedging with one hedging asset. 
one_phi_t = zeros(Simulations, Steps); # Stores phi if asset 1 was the only available hedging asset. Equation (4.3.21)
two_phi_t = zeros(Simulations, Steps); # Stores phi if asset 2 was the only available hedging asset. Equation (4.3.21)
one_Risk_1_asset  = zeros(Simulations, Steps); # Stores the change in the value of the portfolio 
#if asset 1 was the only available hedging asset. Equation (4.3.10)
two_Risk_1_asset  = zeros(Simulations, Steps); # Stores the change in the value of the portfolio 
#if asset 2 was the only available hedging asset. Equation (4.3.10)
one_ExpectedRisk_1_asset = zeros(Steps); # Average change in the value of the portfolio 
#if asset 1 was the only available hedging asset
two_ExpectedRisk_1_asset = zeros(Steps); # Average change in the value of the portfolio 
#if asset 2 was the only available hedging asset

phi_t_1 = zeros(Simulations, Steps); # Stores phi_1 if both assets 1 and 2 were available as hedging assets. Equation (4.5.1)
phi_t_2 = zeros(Simulations, Steps); # Stores phi_2 if both assets 1 and 2 were available as hedging assets. Equation (4.5.1)
Risk_2_assets  = zeros(Simulations, Steps); # Stores the change in the value of the portfolio 
#if both assets 1 and 2 were available as hedging assets. Equation (4.4.6)
ExpectedRisk_2_assets = zeros(Steps); # Average change in the value of the portfolio 
#if both assets 1 and 2 were available as hedging assets


mean_S1 = zeros(Steps); # Average price of asset 1
mean_S2 = zeros(Steps); # Average price of asset 2
mean_C = zeros(Steps); # Average price of contract

current_C = C_0 .+ zeros(Simulations); # Initializing the value of the contract
current_N = zeros(Simulations); # Initializing the Poisson Process
current_S_1 = S_0_1 .+ zeros(Simulations); # Initializing the value of the risky asset 1
current_S_2 = S_0_2 .+ zeros(Simulations); # Initializing the value of the risky asset 2
current_X  = zeros(Simulations); # Initializing the value of the poisson process
current_W  = zeros(Simulations); # Initializing the value of the brownian motion
