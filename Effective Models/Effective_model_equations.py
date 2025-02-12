## File contains each ODE to be solved step by step to get the effective model 

# Import all packages here 

import numpy as np 
import scipy.integrate as scint 
import pandas as pd 


def run_model(param, degtag, tspan):
    
    # Load the parameters as a dictionary
    p = param.valuesdict()

    degtag = degtag # This parameter will dictate if we have degradation tag or not 
    C0 = p['C_0'] # Initial cell population
    A_0 = p['A_0'] # Initial amino acid concentration
    N_plasmid = 100  # Plasmid copy number for a high copy plasmid
 

    def model_odes(x, t, param, degtag = 1):

        # Load the parameters as a dictionary
        p = param.valuesdict()

        ## Unpack and define state variables
        M = x[0]; P = x[1]; Pm = x[2]; Pc = x[3] ; A = x[4]; C = x[5]

        ## Unpack the parameters 
        k_tx = p['k_tx'] # Transcription rate, nM/min
        d_m = p['d_m'] # mRNA degradation rate constant, 1/min
        k_gr = p['k_gr'] # logistic growth rateconstant, 1/min
        C_max = p['C_max'] # max cell count
        k_tl = p['k_tl'] # Protein/min-mRNA, translation elongation rate constant
        d_p = p['d_p'] # 1/min, rate constant of active degradation without deg-tag
        k_fold = p['k_fold'] # 1/min, rate constant of folding of polypeptide
        b_tag = p['b_tag'] # Basal deg-tag degradation rate 
        b_fold = p['b_fold'] # Basal folding rate 
        n = 5.5 # Fixed parameter for RMF delta 
        d_dil = k_gr # 1/min, rate constant of dilution
        k_rep = p['k_rep'] # Amino acid replenishment rate, nM/min
        k_lag = p['k_lag'] # Peptide chain degradation rate constant, 1/min
        n_gamma_resources = p['n_gamma_resources'] # Exponent of gamma
        n_gamma_rate = p['n_gamma_rate'] # Exponent of gamma
        n_gamma_deg = p['n_gamma_deg'] # Exponent of gamma

        # deg-tag term dependent on parameter degtag
        d_tag = p['d_tag'] * degtag # 1//min, rate constant of degradation due to deg-tag


        ## Define the RMF variables 
        f = C/C_max
        P_total = P + Pm 
        y = f * (1 - f) 

        # Define RMFs
        alpha = 1 - f
        delta = f**n/(1 + f**n)
        y_resources = np.power(y, n_gamma_resources)
        y_rate = np.power(y, n_gamma_rate)
        y_deg = np.power(y, n_gamma_deg)
        
        # Transcription rate 
        k_tx = k_tx * y_resources 
        beta_m = N_plasmid * k_tx 

        # mRNA degradation rate 
        d_m = d_m * alpha
        
        # Protein folding rate
        k_fold = k_fold * (y_rate  + b_fold)

        # Protein degradation rate
        d_p = d_p * delta
        d_dil = d_dil * alpha
        d_tag = d_tag * (y_deg + b_tag) 
       
        # Amino acid replenishment rate 
        k_rep = k_rep * y_rate 

        # Translation rate 
        k_tl = k_tl * y_resources 
        
        ## Define ODEs here: 
        dMdt = beta_m - (d_m + d_dil) * M
        dPdt = k_tl * M * A - (d_p + d_dil + d_tag + k_fold) * P
        dPmdt = k_fold * P - (d_p + d_dil + d_tag) * Pm
        dPcdt = (d_p + d_tag) * P_total - (k_lag + d_dil) * Pc
        dAdt = k_lag * Pc + k_rep - k_tl * M * A - d_dil * A 
        dCdt = k_gr * C * (1 - f)
        
        
        return [dMdt, dPdt, dPmdt, dPcdt, dAdt, dCdt]
    
    X0 = [0, 0, 0, 0, A_0, C0]
    
    sol = scint.odeint(model_odes, X0, tspan, args = (param, degtag))

    return sol