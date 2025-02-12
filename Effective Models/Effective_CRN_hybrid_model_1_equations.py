## File contains  each ODE to be solved step by step to get the effective model 

# Here we show how to write the equations for the hybrid-CRN model - 1
# We have divided translation into 3 CRNs: (1) Amino acid charging (2) Translation inititaion 
# and (3) Translation elongation
# 

# Import all the packages here 

import numpy as np 
import scipy.integrate as scint 
import pandas as pd 


def run_model(param, degtag, tspan):

    p = param.valuesdict()
    degtag = degtag # This parameter will dictate if we have degradation tag or not 
    C_0 = p['C_0']
    A_0 = p['A_0']
    N_plasmid = 100  # Plasmid copy number for a high copy plasmid
    

    def model_odes(x, t, param, degtag = 1):

        p = param.valuesdict()
        
        ## Unpack and define state variables
        M = x[0]; C_tic = x[1]; P = x[2]; Pm = x[3]; Pc = x[4] ; A = x[5]; C = x[6]; Caa = x[7]

        ## Unpack the parameters 
        k_tx = p['k_tx'] # Transcription rate, nM/min
        d_m = p['d_m'] # mRNA degradation rate constant, 1/min
        k_gr = p['k_gr'] # logistic growth rate constant, 1/min
        C_max = p['C_max'] # max cell counts
        k_tli_b = p['k_tli_b'] # Translation initiation complexbinding rate constant, 1/min.nM^2
        k_tli_u = p['k_tli_u'] # Translation initiation complex unbinding rate constant, 1/min
        k_tl = p['k_tl'] # Protein/min-mRNA, translation elongation rate constant
        d_p = p['d_p'] # 1/min, rate constant of active degradation without deg-tag
        k_fold = p['k_fold'] # 1/min, rate constant of folding of polypeptide
        b_fold = p['b_fold'] # Basal folding rate 
        b_tag = p['b_tag'] # Basal deg-tag degradation rate 
        n = 5.5
        d_dil = k_gr # 1/min, rate of dilution
        k_rep = p['k_rep'] # Amino acid replenishment rate, nM/min
        k_lag = p['k_lag'] # Peptide chain degradation rate constant, 1/min
        k_tlaa_b = p['k_tlaa_b'] # Amino acid charging forward rate constant, 1/min.nM
        k_tlaa_u = p['k_tlaa_u'] # Amino acid charging reverese rate constant, 1/min
        n_gamma_resources = p['n_gamma_resources'] # Exponent of gamma
        n_gamma_rate = p['n_gamma_rate'] # Exponent of gamma
        tRNA_max = p['tRNA_max'] # Max cap on tRNA resource
        Ribo_max = p['Ribo_max']  # Max cap on Ribosome resource
    
        # deg-tag term dependent on blank step 
        d_tag = p['d_tag'] * degtag # 1//min, rate constant of degradation due to deg-tag

        # Define the RMF variables 
        f = C/C_max
        Pt = P + Pm
        y_square = f * (1 - f)
        y_resources = np.power(y_square, n_gamma_resources)
        y_rate = np.power(y_square, n_gamma_rate)
        
        # Conservation equations for resources
        tRNA_total = tRNA_max * y_resources
        tRNA_free = tRNA_total - Caa

        Ribo_total = Ribo_max * y_resources
        Ribo_free = Ribo_total - C_tic
        
        # Transcription rate 
        k_tx = k_tx * y_resources
        beta_m = N_plasmid * k_tx

        # mRNA degradation rate 
        d_m = d_m * (1 - f)
        
        # Protein folding rate
        k_fold = k_fold * (y_rate + b_fold)

        # Protein degradation rate
        d_p = d_p * f**n / (1 + f**n)
        d_dil = d_dil * (1 - f)
        d_tag = d_tag * (y_rate + b_tag) 
       
        # Amino acid replenishment rate 
        k_rep = k_rep * y_rate

        # Amino acid charging
        r_aa_charging = k_tlaa_b * A * tRNA_free - k_tlaa_u * Caa 

        # Translation inititaiton and elongation
        r_translation_inititation = k_tli_b * Ribo_free * M - k_tli_u * C_tic 
        r_translation_elongation = k_tl * Caa * C_tic 
        
        ## Define ODEs here: 
        dMdt = beta_m - (d_m + d_dil) * M - r_translation_inititation + r_translation_elongation
        dC_ticdt = r_translation_inititation - r_translation_elongation
        dPdt = r_translation_elongation - (d_p + d_dil + d_tag + k_fold) * P
        dPmdt = k_fold * P - (d_p + d_dil + d_tag) * Pm
        dPcdt = (d_p + d_tag) * Pt - (k_lag + d_dil) * Pc
        dAdt = k_lag * Pc + k_rep - r_aa_charging - d_dil * A
        dCdt = k_gr * C * (1 - f)
        dCaadt = r_aa_charging - r_translation_elongation
        
        
        return [dMdt, dC_ticdt, dPdt, dPmdt, dPcdt, dAdt, dCdt, dCaadt]
    
    
    X0 = [0, 0, 0, 0, 0, A_0, C_0, 0]

    sol = scint.odeint(model_odes, X0, tspan, args = (param, degtag))

    return sol