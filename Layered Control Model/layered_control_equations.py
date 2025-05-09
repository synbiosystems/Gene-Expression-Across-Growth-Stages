# This file contains the minimal Layered Control Model as it is in the paper 


# Import here 

import numpy as np 
import scipy.integrate as scint 

# Creating a sample array just in case the user does not pass it
t_final = 960 ; N_steps = 96 
timepoints = np.linspace(0, t_final, N_steps)

def run_model(x0, param, loop = 'open', tspan = timepoints):

    N_plasmid_high = 100
    N_plasmid_low = 10

    def model_ODEs(t, x, p, loop):

        # unpack the states: 
        M_cin = x[0]; P_cin = x[1]; C_cin = x[2]; M_G = x[3]; R = x[4]; Rm = x[5]; P_lac = x[6]; 
        Pm_lac = x[7]; P_fp = x[8]; Pm_fp = x[9]; C = x[10]; TI_cin = x[11]; TI_lac = x[12]; TI_fp = x[13]

        # Unpack the parameters
        p = param.valuesdict()

        # CinR cassette: 
        beta_A = p['beta_A'] # Transcription rate of first cassette
        K_x = p['K_x'] # activation coefficient of the chemical inducer x
        x_i = p['x_i'] # the chemical inducer that activates P_rhl/lac promoter
        k_tli_b_cin = p['k_tli_b_cin'] # Translational resource-binding rate constant for CinR
        k_tl_cin = p['k_tl_cin'] # Translation elongation rate of CinR
        K_r = p['K_r'] # Complex formation rate of CinR-Cin complex
        K_lac = p['K_lac'] # Repression constant for P_rhl/lac promoter

        # Common parameters
        k_tli_u = p['k_tli_u'] # Translational resource-unbinding rate constant
        d_m = p['d_m'] # mRNA degradation rate
        d_p = p['d_p'] # Protein degradation rate

        # FP cassette 
        beta_B = p['beta_B'] # Transcription rate of second cassette
        K_cin = p['K_cin'] # Activation constant for P_cin promoter
        l0 = p['l0_Pcin'] ; # leak coefficient for P_cin promoter
        k_fold_r = p['k_fold_r'] # Folding rate of sRNA
        d_r = p['d_r'] # degradation rate of sRNA
        K_R = p['K_R'] # Repression constant for attenuator
        k_tli_b_lac = p['k_tli_b_lac'] # Translational resource-binding rate constant for LacI
        k_tl_lac = p['k_tl_lac'] # Translation elongation rate of LacI
        k_fold_lac = p['k_fold_lac'] # Folding rate of LacI
        k_tli_b_fp = p['k_tli_b_fp'] # Translational resource-binding rate constant for sfYFP
        k_tl_fp = p['k_tl_fp'] # Translation elongation rate of sfYFP
        k_fold_fp = p['k_fold_fp'] # Folding rate of sfYFP

        # RMF related parameters
        b_fold = p['b_fold']
        K_P_total = p['K_P_total']
        K_R_total = p['K_R_total'] # This parameter is redundant. Not used in the model at all.
        Rsrc_max = p['Rsrc_max']
        n_gamma_resources_tx = p['n_gamma_resources_tx'] # Exponent of gamma
        n_gamma_resources_tl = p['n_gamma_resources_tl'] # Exponent of gamma
        n_gamma_rate = p['n_gamma_rate'] # Exponent of gamma
        
        # Load the hill coefficients in case they are passed in lmfit object:
        n_Cin_actn = p.get('n_Cin_actn', 2)
        n_Rhl_actn = p.get('n_Rhl_actn', 1)
        n_Hill_delta = p.get('n_Hill_delta', 8.5)
        n_Hill_psi = p.get('n_Hill_psi', 3)
        n_Hill_trans = p.get('n_Hill_trans', 1)
        n_Hill_cis = p.get('n_Hill_cis', 1)

        # Growth related parameters 
        C_max = p['C_max']
        k_gr = p['k_gr']

        # Define the feedback functions
        f_Hill_trans = K_lac**n_Hill_trans/(K_lac**n_Hill_trans + Pm_lac**n_Hill_trans) # Hill function for trans repression
        f_Hill_cis = K_R**n_Hill_cis/(K_R**n_Hill_cis + Rm**n_Hill_cis) # Hill function for cis repression

        # Define conditions where specific feedback should be applied
        if loop == 'open': 
            f_cis = 1
            f_trans = 1
        
        elif loop == 'cis': 
            f_cis = f_Hill_cis 
            f_trans = 1

        elif loop == 'trans':
            f_cis = 1
            f_trans = f_Hill_trans

        elif loop == 'layered':
            f_cis = f_Hill_cis
            f_trans = f_Hill_trans

        ## apply RMFs: 

        f = C/C_max
        gamma_base = f * (1 - f) 

        gamma_resources_tx = np.power(gamma_base, n_gamma_resources_tx)
        gamma_resources_tl = np.power(gamma_base, n_gamma_resources_tl)
        gamma_rate = np.power(gamma_base, n_gamma_rate)
        delta = f**n_Hill_delta/(1 + f**n_Hill_delta)
        alpha = 1 - f

        P_total = P_cin + C_cin + P_lac + Pm_lac + P_fp + Pm_fp

        psi_P = K_P_total**n_Hill_psi/(K_P_total**n_Hill_psi + P_total**n_Hill_psi)

        # Rsrc Conservation equation
        Rsrc_total = Rsrc_max * gamma_resources_tl
        Rsrc_free = Rsrc_total - (TI_cin + TI_lac + TI_fp)

        # dilution 
        d_dil = k_gr * alpha

        # transcription
        beta_A = beta_A * gamma_resources_tx * N_plasmid_low 
        beta_B = beta_B * gamma_resources_tx * N_plasmid_high  

        # mRNA degradation 
        d_m = d_m * alpha + d_dil 
        d_r = d_r * alpha + d_dil 

        # protein folding
        k_fold_fp = k_fold_fp * (gamma_rate * (psi_P) + b_fold)
        k_fold_lac = k_fold_lac * (gamma_rate * (psi_P) + b_fold)
        k_fold_r = k_fold_r 

        # protein degradation as function of growth
        d_p = d_p * delta * (1 - psi_P)

        # Translation rates:
        # CinR Cassette:
        r_translation_initiation_cin = k_tli_b_cin * Rsrc_free * M_cin - k_tli_u * TI_cin
        r_translation_elongation_cin = k_tl_cin * TI_cin

        # sfYFP cassette: 
        r_translation_initiation_lac = k_tli_b_lac * Rsrc_free * M_G - k_tli_u * TI_lac
        r_translation_elongation_lac = k_tl_lac * TI_lac
        r_translation_initiation_fp = k_tli_b_fp * Rsrc_free * M_G - k_tli_u * TI_fp
        r_translation_elongation_fp = k_tl_fp * TI_fp

        odes = []

        dMcindt = beta_A * (x_i**n_Rhl_actn / (K_x + x_i**n_Rhl_actn)) * f_trans - d_m * M_cin - r_translation_initiation_cin + r_translation_elongation_cin
        odes.append(dMcindt)

        dPcindt = r_translation_elongation_cin - (d_p + d_dil + K_r) * P_cin
        odes.append(dPcindt)

        dCcindt = K_r * P_cin - (d_p + d_dil) * C_cin
        odes.append(dCcindt)

        dMgdt = beta_B * ((C_cin**n_Cin_actn/(K_cin**n_Cin_actn + C_cin**n_Cin_actn)) + l0) * f_cis - d_m * M_G - (r_translation_initiation_lac + r_translation_initiation_fp) + (r_translation_elongation_lac + r_translation_elongation_fp)
        odes.append(dMgdt)

        dRdt = beta_B * ((C_cin**n_Cin_actn/(K_cin**n_Cin_actn + C_cin**n_Cin_actn)) + l0) * f_cis - (d_r + k_fold_r) * R
        odes.append(dRdt)

        dRmdt = k_fold_r * R - d_r * Rm
        odes.append(dRmdt)

        dPlacdt = r_translation_elongation_lac - (d_p + d_dil + k_fold_lac) * P_lac
        odes.append(dPlacdt)

        dPmlacdt = k_fold_lac * P_lac - (d_p + d_dil) * Pm_lac
        odes.append(dPmlacdt)

        dPfpdt = r_translation_elongation_fp - (d_p + d_dil + k_fold_fp) * P_fp
        odes.append(dPfpdt)

        dPmfpdt = k_fold_fp * P_fp - (d_p + d_dil) * Pm_fp
        odes.append(dPmfpdt)

        dCdt = k_gr * C * (1 - f)
        odes.append(dCdt)

        dTI_cindt = r_translation_initiation_cin - r_translation_elongation_cin 
        odes.append(dTI_cindt)

        dTI_lacdt = r_translation_initiation_lac - r_translation_elongation_lac 
        odes.append(dTI_lacdt)

        dTI_yfpdt = r_translation_initiation_fp - r_translation_elongation_fp 
        odes.append(dTI_yfpdt)

        return np.array(odes)
    
    # Unpack the initial conditions
    p = param.valuesdict()
    
    x0[2] = p['CinR_0']
    x0[7] = p['LacI_0']
    x0[9] = p['FP_0']
    x0[0] = p['Mcin_0']
    x0[3] = p['Mg_0']
    x0[10] = p['C_0']

    if loop == 'open':
        x0[0] = p['Mcin_0'] * 6 # More mRNA for Cin in open compared to trans since trans has repression
        x0[2] = p['CinR_0'] * 2 # More Cin in open compared to trans since trans has repression
        x0[3] = p['Mg_0'] * 4 # More YFP initial condition in open compared to cis and layered and more than trans
    
    if loop == 'cis': 
        x0[0] = p['Mcin_0'] * 6 # More mRNA for Cin in cis compared to trans since trans has repression
        x0[2] = p['CinR_0'] * 2 # More Cin in cis compared to trans since trans has repression

    
    if loop == 'trans': 
        x0[3] = p['Mg_0'] * 4 # More YFP initial condition in trans compared to cis and layered but less than open
    
    if loop == 'layered':
        x0[0] = p['Mcin_0'] * 3 # More mRNA for Cin in layered compared to trans since less leaky LacI made
        x0[2] = p['CinR_0'] * 3 # More Cin in layered compared to trans since less leaky LacI made
    
    
    
    ### NEW: First always run open loop internally
    ## This step is used only while fitting parameters. Once you get a good fit, you can skip this step
    ## and just run the model with the parameters you got from fitting.
    ## This is done to get the new K_lac and K_R values based on the maximum of Pm_lac and Rm
    ## and C_cin over time.
    sol_open = scint.solve_ivp(
        model_ODEs, [tspan[0], tspan[-1]], np.copy(x0),
        args=(param, 'open'), method='LSODA', t_eval=tspan
    )
    sol_open_y = sol_open.y.T

    # Extract max state
    Pm_lac_max = np.max(sol_open_y[:, 7]) # Maximum of Pm_lac over time
    Rm_max = np.max(sol_open_y[:, 5])     # Maximum of Rm over time
    C_cin_max = np.max(sol_open_y[:, 2])     # Maximum of C_cin over time

    # Step 2: Set new K_lac and K_R properly in lmfit Parameters
    scaling_factor_lac = p['r_Lac']  # (adjust if needed)
    scaling_factor_R = p['r_R']  # (adjust if needed)
    scaling_factor_cin = p['r_Cin']  # (adjust if needed)
    new_K_lac = (Pm_lac_max) * scaling_factor_lac
    new_K_R = (Rm_max) * scaling_factor_R
    new_K_cin = (C_cin_max) * scaling_factor_cin

    ## Uncomment when you want to update the new K_lac, K_R and K_cin values
    # param['K_lac'].set(value=new_K_lac)
    # param['K_R'].set(value=new_K_R)
    # param['K_cin'].set(value=new_K_cin)


    sol_final_run = scint.solve_ivp(model_ODEs, [tspan[0], tspan[-1]], x0, args = (param, loop),
    method = 'LSODA', t_eval = tspan)

    # Extract results in odeint format
    sol = sol_final_run.y.T  
    
    P_total = sol[:,1] + sol[:,2] + sol[:,6] + sol[:,7] + sol[:,8] + sol[:,9]

    return sol, P_total
