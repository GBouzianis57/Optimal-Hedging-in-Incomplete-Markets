module calc

using Distributions

function computations(T, Steps, Simulations, sigma_1, sigma_2, sigma_H, beta_1, beta_2, beta_H, S_0_1, S_0_2, C_0, m, jump, p)
    
#######################################Initialization##############################################

W_t  = zeros(Simulations, Steps); # Brownian motion
N_t  = zeros(Simulations, Steps); # Poisson process
X_t  = zeros(Simulations, Steps); # Compound Poisson process
S_t_1  = zeros(Simulations, Steps); # Price of asset 1 based on equation (4.5.3)
S_t_2  = zeros(Simulations, Steps); # Price of asset 2 based on equation (4.5.3)
C_t  = zeros(Simulations, Steps); # Price of the contract based on equation (4.5.3)

##Optimal hedging with one hedging asset. 
one_phi_t = zeros(Simulations, Steps); # Stores phi if asset 1 was the only available hedging asset. Equation (4.3.21)
two_phi_t = zeros(Simulations, Steps); # Stores phi if asset 2 was the only avaislable hedging asset. Equation (4.3.21)
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

#######################################Calculations##############################################
    
#Note: Here sigma is used for jump volatilities and beta is used to brownian volatilities. In the paper itself, in the "Simulations" section, sigma is used for the brownian volatilities and beta for the jump volatilities. The formulas in the code are adjusted based to the notation that is described in this note. Here i=1 and j=2.
    
#Note: Here we simulate the processes based on the selected number of scenarios. If the formula that is used for any variable is unclear,
#look back in the "Initialization" section above, and see the comments next to the variable. It gives the equation numbers in parallel to the paper. 
    
dt = T/Steps # Instantaneous change in time based on the selected values of T and Steps. 

for Step_Id = 1 : Steps
    
        one_phi_t[:, Step_Id] .=  (current_C./current_S_1) .* ( beta_1*beta_H + (m* ( p*(exp(sigma_1*jump[1])-1)*(exp(sigma_H*jump[1])-1) + (1-p)*(exp(sigma_1*jump[2])-1)*(exp(sigma_H*jump[2])-1) ) ) ) / ( beta_1^2 + m*( p*(exp(sigma_1*jump[1])-1)^2 +  (1-p)*(exp(sigma_1*jump[2])-1)^2  ))
    
        two_phi_t[:, Step_Id] .=  (current_C./current_S_2) .* ( beta_2*beta_H + (m* ( p*(exp(sigma_2*jump[1])-1)*(exp(sigma_H*jump[1])-1) + (1-p)*(exp(sigma_2*jump[2])-1)*(exp(sigma_H*jump[2])-1) ) ) ) / ( beta_2^2 + m*( p*(exp(sigma_2*jump[1])-1)^2 +  (1-p)*(exp(sigma_2*jump[2])-1)^2  ))
    
        
        E = beta_2^2 + m*( p*(exp(sigma_2*jump[1])-1)^2 +  (1-p)*(exp(sigma_2*jump[2])-1)^2 )
    
        F = beta_1^2 + m*( p*(exp(sigma_1*jump[1])-1)^2 +  (1-p)*(exp(sigma_1*jump[2])-1)^2 )
    
        G = (beta_1*beta_2 + m*( p*(exp(sigma_1*jump[1])-1)*(exp(sigma_2*jump[1])-1) +  (1-p)*(exp(sigma_1*jump[2])-1)*(exp(sigma_2*jump[2])-1) ) )^2
        
        D = (E * F) - (G)
        
        A1 = beta_1*beta_2 + m*( p*(exp(sigma_1*jump[1])-1)*(exp(sigma_2*jump[1])-1) +  (1-p)*(exp(sigma_1*jump[2])-1)*(exp(sigma_2*jump[2])-1) )
        
        B1 = beta_H*beta_2 + m*( p*(exp(sigma_H*jump[1])-1)*(exp(sigma_2*jump[1])-1) +  (1-p)*(exp(sigma_H*jump[2])-1)*(exp(sigma_2*jump[2])-1) )
    
        C1 = beta_H*beta_1 + m*( p*(exp(sigma_H*jump[1])-1)*(exp(sigma_1*jump[1])-1) +  (1-p)*(exp(sigma_H*jump[2])-1)*(exp(sigma_1*jump[2])-1) )
    
        D1 = beta_2^2 + m*( p*(exp(sigma_2*jump[1])-1)^2 +  (1-p)*(exp(sigma_2*jump[2])-1)^2 )

        A2 = beta_1*beta_2 + m*( p*(exp(sigma_1*jump[1])-1)*(exp(sigma_2*jump[1])-1) +  (1-p)*(exp(sigma_1*jump[2])-1)*(exp(sigma_2*jump[2])-1) )
        
        B2 = beta_H*beta_1 + m*( p*(exp(sigma_H*jump[1])-1)*(exp(sigma_1*jump[1])-1) +  (1-p)*(exp(sigma_H*jump[2])-1)*(exp(sigma_1*jump[2])-1) )
    
        C2 = beta_H*beta_2 + m*( p*(exp(sigma_H*jump[1])-1)*(exp(sigma_2*jump[1])-1) +  (1-p)*(exp(sigma_H*jump[2])-1)*(exp(sigma_2*jump[2])-1) )
    
        D2 = beta_1^2 + m*( p*(exp(sigma_1*jump[1])-1)^2 +  (1-p)*(exp(sigma_1*jump[2])-1)^2 )
        
        phi_t_1[:, Step_Id] .= (current_C ./ current_S_1) .* ( ( (C1)*(D1) - (A1)*(B1) ) / D )
    
        phi_t_2[:, Step_Id] .= (current_C ./ current_S_2) .* ( ( (C2)*(D2) - (A2)*(B2) ) / D )
    
        t = (Step_Id) * dt;
    
        W_t[:, Step_Id] = current_W + sqrt(dt)*randn(Simulations)
              
        N_t[:, Step_Id] = current_N .+ (rand(Simulations) .< (1 - exp(-m * dt)));
    
        Y_t =  rand(Binomial(1, p), Simulations) 
    
        Y_t[Y_t .== 1] .= jump[1] # Jumps of type (1) with probability p
    
        Y_t[Y_t .== 0] .= jump[2] # Jumps of type (2) with probability 1-p
    
        X_t[N_t[:, Step_Id] .!= current_N, Step_Id] .= current_X[ N_t[:, Step_Id] .!= current_N ] .+  Y_t[N_t[:, Step_Id] .!= current_N]
    
        X_t[N_t[:, Step_Id] .== current_N, Step_Id] .= current_X[ N_t[:, Step_Id] .== current_N ]
    
        S_t_1[:, Step_Id] .= S_0_1 .* exp.(beta_1.* W_t[:, Step_Id] .- (1/2)*(beta_1^2)*t .+ sigma_1.* X_t[:, Step_Id] .- m*t*( p*(exp(jump[1]*sigma_1)) + (1-p)*(exp(jump[2]*sigma_1)) - 1))
    
        S_t_2[:, Step_Id] .= S_0_2 .* exp.(beta_2.* W_t[:, Step_Id] .- (1/2)*(beta_2^2)*t .+ sigma_2.* X_t[:, Step_Id] .- m*t*( p*(exp(jump[1]*sigma_2)) + (1-p)*(exp(jump[2]*sigma_2)) - 1))
    
        C_t[:, Step_Id] .=  C_0 .* exp.(beta_H.* W_t[:, Step_Id] .- (1/2)*(beta_H^2)*t .+ sigma_H.* X_t[:, Step_Id] .- m*t*( p*(exp(jump[1]*sigma_H)) + (1-p)*(exp(jump[2]*sigma_H)) - 1))
        
        one_Risk_1_asset[:, Step_Id] .= ( (beta_H .* (W_t[:, Step_Id] .- current_W) .+ (exp.(sigma_H.*(X_t[:, Step_Id] .- current_X)) .- 1)) .* current_C) .- ( (beta_1 .* (W_t[:, Step_Id] .- current_W)  +  (exp.(sigma_1.*(X_t[:, Step_Id] .- current_X)) .- 1) ) .* one_phi_t[:, Step_Id] .* current_S_1) 
    
        two_Risk_1_asset[:, Step_Id] .= ( (beta_H .* (W_t[:, Step_Id] .- current_W) .+ (exp.(sigma_H.*(X_t[:, Step_Id] .- current_X)) .- 1)) .* current_C) .- ( (beta_2 .* (W_t[:, Step_Id] .- current_W)  +  (exp.(sigma_2.*(X_t[:, Step_Id] .- current_X)) .- 1) ) .* two_phi_t[:, Step_Id] .* current_S_2) 

        Risk_2_assets[:, Step_Id] .= ( (beta_H .* (W_t[:, Step_Id] .- current_W) .+ (exp.(sigma_H.*(X_t[:, Step_Id] .- current_X)) .- 1)) .* current_C) .- ( (beta_2 .* (W_t[:, Step_Id] .- current_W)  +  (exp.(sigma_2.*(X_t[:, Step_Id] .- current_X)) .- 1) ) .* phi_t_2[:, Step_Id] .* current_S_2) .- ( (beta_1 .* (W_t[:, Step_Id] .- current_W)  +  (exp.(sigma_1.*(X_t[:, Step_Id] .- current_X)) .- 1) ) .* phi_t_1[:, Step_Id] .* current_S_1)   
       
        one_ExpectedRisk_1_asset[Step_Id] = mean(one_Risk_1_asset[:, Step_Id])
    
        two_ExpectedRisk_1_asset[Step_Id] = mean(two_Risk_1_asset[:, Step_Id])
    
        ExpectedRisk_2_assets[Step_Id] = mean(Risk_2_assets[:, Step_Id])
        
        current_W = W_t[:, Step_Id]
            
        current_N = N_t[:, Step_Id];
    
        current_X = X_t[:, Step_Id];
    
        current_S_1 = S_t_1[:, Step_Id];
    
        current_S_2 = S_t_2[:, Step_Id];
    
        current_C = C_t[:, Step_Id];
    
        mean_S1[Step_Id] = mean(S_t_1[:, Step_Id])
    
        mean_S2[Step_Id] = mean(S_t_2[:, Step_Id])
    
        mean_C[Step_Id] = mean(C_t[:, Step_Id])
   
end #for
    
return(W_t, N_t, X_t, S_t_1, S_t_2, C_t, one_phi_t, two_phi_t, one_Risk_1_asset, two_Risk_1_asset, one_ExpectedRisk_1_asset, two_ExpectedRisk_1_asset, phi_t_1, phi_t_2, Risk_2_assets, ExpectedRisk_2_assets)
end #function

end #module