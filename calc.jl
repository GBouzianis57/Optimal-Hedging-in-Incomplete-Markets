dt = T/Steps
    
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
   
    end