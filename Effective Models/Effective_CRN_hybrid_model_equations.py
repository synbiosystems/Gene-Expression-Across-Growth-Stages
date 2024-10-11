## File contains  each ODE to be solved step by step to get the effective model 

# Here we show how to write the equations for the hybrid-CRN model 
# We have divided translation into 3 CRNs: (1) Amino acid charging (2) Translation inititaion 
# and (3) Translation elongation
# 

# Import all the packages here 

import numpy as np 
import scipy.integrate as scint 
import pandas as pd 


def all_steps(param, degtag, tspan, N_steps = 1000):
    
    p = param.valuesdict()
    degtag = degtag # This parameter will dictate if we have degradation tag or not 
    C0 = p['C_0']
    N_plasmid = 100  # Plasmid copy number for a high copy plasmid
    

    def step_blank(x, t, param, degtag = 1):

        # This step is to run the model without the deg-tag and the RMF on protein folding to get an estimate of 
        # the maximum protein concentration possible 

        p = param.valuesdict()
        
        M = x[0]; R = x[1]; C_tic = x[2]; P = x[3]; Pm = x[4]; Pc = x[5] ; A = x[6]; C = x[7]; Caa = x[8]

        k_tx = p['k_tx'] # Transcription rate, nM/min
    
        d_m = p['d_m'] # mRNA degradation rate, 1/min

        k_gr = p['k_gr'] # logistic growth rate, 1/min
        C_max = p['C_max'] # counts
        
        k_tli_b = p['k_tli_b'] # Translation initiation complex formation rate, 1/min.nM^2
        k_tli_u = p['k_tli_u'] # Translation initiation complex unbinding rate, 1/min
        k_tl = p['k_tl'] # Protein/min-mRNA, translation elongation rate
        d_p = p['d_p'] # 1/min, rate of active degradation without deg-tag
        d_tag = p['d_tag'] * degtag * 0 # 1//min, rate of degradation due to deg-tag
        k_fold = p['k_fold'] # 1/min, rate of folding of polypeptide
        b_tag = p['b_tag'] # Basal deg-tag degradation rate 
        b_tl = p['b_tl'] # Basal translation rate 
        Kp = p['Kp'] # Binding coefficient for protease, nM
        R_min = p['R_min'] # Minimum cap on species R, nM
        k_R = p['k_R'] # Synthesis/degradation rate of R, 1/min


        d_dil = k_gr # 1/min, rate of dilution
        
        
        k_rep = p['k_rep'] # Amino acid replenishment rate, nM/min
        k_lag = p['k_lag'] # Peptide chain degradation rate, 1/min
        k_tlaa_b = p['k_tlaa_b'] # Amino acid charging forward rate, 1/min.nM
        k_tlaa_u = p['k_tlaa_u'] # Amino acid charging reverese rate, 1/min
    
        

        # Defining the RMF variables 

        f = C/C_max
       
        Pt = P + Pm

        y = f * (1 - f) 
        
        # Transcription rate 
        k_tx = k_tx * y 
        beta_m = N_plasmid * k_tx

        # mRNA degradation rate 
        d_m = d_m * (1 - f)
        
        # Protein folding rate, no RMF in this step 
        k_fold = k_fold #* (y * (1 - p)  + b_fold)

        # Protein degradation rate
        d_p = d_p * f
        d_dil = d_dil * (1 - f)
        d_tag = d_tag * (y + b_tag) 
       
        # Amino acid replenishment rate 
        k_rep = k_rep * y 

        # Species concentration rate 
        k_R = k_R * (R + Caa - R_min) * (1 - 2 * f)

        # Amino acid charging
        aa_charging = 1 * k_tlaa_b * A * R - k_tlaa_u * Caa * 1

        # Translation inititaiton and elongation
        transation_inititation = y * k_tli_b * Caa * M - k_tli_u * C_tic * 1
        translation_elongation = k_tl * C_tic * (1 - f + b_tl)
        

        dMdt = beta_m - (d_m + d_dil) * M - transation_inititation +  translation_elongation
        dRdt = - aa_charging +  translation_elongation + k_R
        dC_ticdt = transation_inititation - translation_elongation
        dPdt = translation_elongation - (d_p + d_dil + d_tag + k_fold) * P
        dPmdt = k_fold * P - (d_p + d_dil + d_tag) * Pm
        dPcdt = (d_p + d_tag) * Pt - (k_lag + d_dil) * Pc
        dAdt = k_lag * Pc + k_rep - aa_charging - d_dil * A
        dCdt = k_gr * C * (1 - f)
        dCaadt = aa_charging - transation_inititation
        
        
        return [dMdt, dRdt, dC_ticdt, dPdt, dPmdt, dPcdt, dAdt, dCdt, dCaadt]
    
    AA_0 = p['A_0']
    R_0 = p['R_0']
    X0 = [0, R_0, 0, 0, 0, 0, AA_0, C0, 0]

    sol_blank = scint.odeint(step_blank, X0, tspan, args = (param, degtag))

    # Estimating P_max
    P_max_est = np.max(sol_blank[:,3] + sol_blank[:,4])
   

    def step_final(x, t, param, degtag = 1):

        p = param.valuesdict()
        
        M = x[0]; R = x[1]; C_tic = x[2]; P = x[3]; Pm = x[4]; Pc = x[5] ; A = x[6]; C = x[7]; Caa = x[8]

        k_tx = p['k_tx'] # Transcription rate, nM/min
    
        d_m = p['d_m'] # mRNA degradation rate, 1/min

        k_gr = p['k_gr'] # logistic growth rate, 1/min
        C_max = p['C_max'] # counts
        
        k_tli_b = p['k_tli_b'] # Translation initiation complex formation rate, 1/min.nM^2
        k_tli_u = p['k_tli_u'] # Translation initiation complex unbinding rate, 1/min
        k_tl = p['k_tl'] # Protein/min-mRNA, translation elongation rate
        d_p = p['d_p'] # 1/min, rate of active degradation without deg-tag
        d_tag = p['d_tag'] * degtag * 1 # 1//min, rate of degradation due to deg-tag
        k_fold = p['k_fold'] # 1/min, rate of folding of polypeptide
        b_fold = p['b_fold'] # Basal folding rate 
        b_tag = p['b_tag'] # Basal deg-tag degradation rate 
        b_tl = p['b_tl'] # Basal translation rate 
        Kp = p['Kp'] # Binding coefficient for protease, nM
        R_min = p['R_min'] # Minimum cap on species R, nM
        k_R = p['k_R'] # Synthesis/degradation rate of R, 1/min


        d_dil = k_gr # 1/min, rate of dilution
        
        
        k_rep = p['k_rep'] # Amino acid replenishment rate, nM/min
        k_lag = p['k_lag'] # Peptide chain degradation rate, 1/min
        k_tlaa_b = p['k_tlaa_b'] # Amino acid charging forward rate, 1/min.nM
        k_tlaa_u = p['k_tlaa_u'] # Amino acid charging reverese rate, 1/min
    
        

        # Defining the RMF variables 

        f = C/C_max
       
        Pt = P + Pm
        p_ratio = Pt/P_max_est

        y = f * (1 - f) 
        
        # Transcription rate 
        k_tx = k_tx * y 
        beta_m = N_plasmid * k_tx

        # mRNA degradation rate 
        d_m = d_m * (1 - f)
        
        # Protein folding rate
        k_fold = k_fold * (y * (1 -  p_ratio)  + b_fold)

        # Protein degradation rate
        d_p = d_p * f
        d_dil = d_dil * (1 - f)
        d_tag = d_tag * (y + b_tag) 
       
        # Amino acid replenishment rate 
        k_rep = k_rep * y 

        # Species concentration rate 
        k_R = k_R * (R + Caa - R_min) * (1 - 2 * f)

        # Amino acid charging
        aa_charging = 1 * k_tlaa_b * A * R - k_tlaa_u * Caa * 1

        # Translation inititaiton and elongation
        transation_inititation = y * k_tli_b * Caa * M - k_tli_u * C_tic * 1
        translation_elongation = k_tl * C_tic * (1 - f + b_tl)
        

        dMdt = beta_m - (d_m + d_dil) * M - transation_inititation +  translation_elongation
        dRdt = - aa_charging +  translation_elongation + k_R
        dC_ticdt = transation_inititation - translation_elongation
        dPdt = translation_elongation - (d_p + d_dil + d_tag + k_fold) * P
        dPmdt = k_fold * P - (d_p + d_dil + d_tag) * Pm
        dPcdt = (d_p + d_tag) * Pt - (k_lag + d_dil) * Pc
        dAdt = k_lag * Pc + k_rep - aa_charging - d_dil * A
        dCdt = k_gr * C * (1 - f)
        dCaadt = aa_charging - transation_inititation
        
        
        
        return [dMdt, dRdt, dC_ticdt, dPdt, dPmdt, dPcdt, dAdt, dCdt, dCaadt]

    AA_0 = p['A_0']
    R_0 = p['R_0']
    X0 = [0, R_0, 0, 0, 0, 0, AA_0, C0, 0]

    sol_final = scint.odeint(step_final, X0, tspan, args = (param, degtag))

    return sol_final, P_max_est