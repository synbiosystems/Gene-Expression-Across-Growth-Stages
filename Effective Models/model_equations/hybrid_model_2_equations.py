## File contains  each ODE to be solved step by step to get the effective model 

# Here we show how to write the equations for the hybrid-CRN model - 2
# We have divided translation into 3 CRNs: (1) Amino acid charging (2) Translation inititaion 
# and (3) Translation elongation
# In this hybrid model we also split deg-tag to account for the binding and unbinding of proteain 
# folded and unfolded, with Protease

# Import all the packages here 

import numpy as np 
import scipy.integrate as scint  


def run_model(param, degtag, tspan):

    p = param.valuesdict()
    degtag = degtag # This parameter will dictate if we have degradation tag or not 
    C0 = p['C_0']
    AA_0 = p['A_0']
    N_plasmid = 100  # Plasmid copy number for a high copy plasmid
    

    def model_odes(x, t, param, degtag = 1):
        
        # load parameters into dictionary
        p = param.valuesdict()
        
        ## Unpack and define state variables
        M = x[0]; C_tic = x[1]; P = x[2]; Pm = x[3]; Pc = x[4] ; 
        A = x[5]; C = x[6]; Caa = x[7]; C_deg_unfolded = x[8]; C_deg_folded = x[9]

        ## Unpack the parameters 
        k_tx = p['k_tx'] # Transcription rate, nM/min
        d_m = p['d_m'] # mRNA degradation rate constant, 1/min
        k_gr = p['k_gr'] # logistic growth rate constant, 1/min
        C_max = p['C_max'] # max cell counts
        k_tli_b = p['k_tli_b'] # Translation initiation complex formation rate constant, 1/min.nM^2
        k_tli_u = p['k_tli_u'] # Translation initiation complex unbinding rate constant, 1/min
        k_tl = p['k_tl'] # Protein/min-mRNA, translation elongation rate constant
        d_p = p['d_p'] # 1/min, rate constant of active degradation without deg-tag
        k_fold = p['k_fold'] # 1/min, rate constant of folding of polypeptide
        b_fold = p['b_fold'] # Basal folding rate  
        n = 5.5 # Fixed parameter for RMF delta
        d_dil = k_gr # 1/min, rate constant of dilution
        k_syn = p['k_syn'] # Amino acid synthesis rate, nM/min
        k_lag = p['k_lag'] # Peptide chain degradation rate constant, 1/min
        k_tlaa_b = p['k_tlaa_b'] # Amino acid charging forward rate constant, 1/min.nM
        k_tlaa_u = p['k_tlaa_u'] # Amino acid charging reverese rate constant, 1/min
        k_protease_b_unfolded = p['k_protease_b_unfolded'] # sfYFP-Protease binding rate constant, 1/min.nM
        k_protease_b_folded = p['k_protease_b_folded'] # sfYFP-Protease binding rate constant, 1/min.nM
        k_protease_u = p['k_protease_u'] # sfYFP-Protease unbinding rate constant, 1/min.nM
        k_c_deg = p['k_c_deg'] # Protease-sfYFP complex degradatoin rate constant, 1/min
        n_gamma_resources = p['n_gamma_resources'] # Exponent of gamma for tx, tl resources
        n_gamma_rate = p['n_gamma_rate'] # Exponent of gamma for rates
        n_gamma_deg = p['n_gamma_deg'] # Exponent of gamma for protease
        tRNA_max = p['tRNA_max'] # Max cap on tRNA resource
        Ribo_max = p['Ribo_max'] # Max cap on Ribo resource
        Protease_max = p['Protease_max'] # Max cap on Protease resource
            

        # Define the RMF variables 
        f = C/C_max
        Pt = P + Pm
        gamma_base = f * (1 - f)
        gamma_resources = np.power(gamma_base, n_gamma_resources)
        gamma_deg = np.power(gamma_base, n_gamma_deg)
        gamma_rate = np.power(gamma_base, n_gamma_rate)

        # Setting up conservation laws
        tRNA_total = tRNA_max * gamma_resources
        tRNA_free = tRNA_total - Caa
        
        Ribo_total = Ribo_max * gamma_resources
        Ribo_free = Ribo_total - C_tic
        
        Protease_total = Protease_max * gamma_deg
        Protease_free = Protease_total - C_deg_folded - C_deg_unfolded
        
        # Defining deg-tag dependence
        Protease_total = Protease_total * degtag
        k_protease_b_unfolded = p['k_protease_b_unfolded'] * degtag 
        k_protease_b_folded = p['k_protease_b_folded'] * degtag 
        k_protease_u = p['k_protease_u'] * degtag 
        k_c_deg = p['k_c_deg'] * degtag 

        # Transcription rate 
        k_tx = k_tx * gamma_resources
        beta_m = N_plasmid * k_tx

        # mRNA degradation rate 
        d_m = d_m * (1 - f)
        
        k_fold = k_fold * (gamma_rate + b_fold)

        # Protein degradation rate
        d_p = d_p * f**n / (1 + f**n)
        d_dil = d_dil * (1 - f)
       
        # Amino acid synthesis rate 
        k_syn = k_syn * gamma_rate


        # Amino acid charging
        aa_charging = k_tlaa_b * A * tRNA_free - k_tlaa_u * Caa 

        # Translation inititaiton and elongation
        r_transation_inititation = k_tli_b * Ribo_free * M - k_tli_u * C_tic 
        r_translation_elongation = k_tl * C_tic * Caa

        # Protease mediated degradation rate:
        protease_binding_rate_unfolded = k_protease_b_unfolded * Protease_free * P - k_protease_u * C_deg_unfolded
        protease_binding_rate_folded = k_protease_b_folded * Protease_free * Pm - k_protease_u * C_deg_folded
        
        ## Define ODEs here: 

        dMdt = beta_m - (d_m + d_dil) * M - r_transation_inititation + r_translation_elongation
        dC_ticdt = r_transation_inititation - r_translation_elongation
        dPdt = r_translation_elongation - (d_p + d_dil + k_fold) * P - protease_binding_rate_unfolded
        dPmdt = k_fold * P - (d_p + d_dil) * Pm - protease_binding_rate_folded
        dPcdt = (C_deg_unfolded + C_deg_folded) * k_c_deg + d_p * Pt - (k_lag + d_dil) * Pc
        dAdt = k_lag * Pc + k_syn - aa_charging - d_dil * A
        dCdt = k_gr * C * (1 - f)
        dCaadt = aa_charging - r_translation_elongation
        dC_deg_unfoldeddt = protease_binding_rate_unfolded - C_deg_unfolded * k_c_deg
        dC_deg_foldeddt = protease_binding_rate_folded - C_deg_folded * k_c_deg
        
        
        return [dMdt, dC_ticdt, dPdt, dPmdt, dPcdt, dAdt, dCdt, dCaadt, dC_deg_unfoldeddt, dC_deg_foldeddt]
    
    
    X0 = [0, 0, 0, 0, 0, AA_0, C0, 0, 0, 0]

    sol = scint.odeint(model_odes, X0, tspan, args = (param, degtag))

    return sol