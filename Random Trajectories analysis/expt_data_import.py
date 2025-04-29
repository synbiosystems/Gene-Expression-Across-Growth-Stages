## File to import all the GEAGS expt data so that same code need not be called all the time we need the expt data

import pandas as pd
import numpy as np 

def Get_OD_Data():
    # Importing growth data
    OD_data = pd.read_csv("expt_growth_data.csv")
    time = OD_data["time(min)"].to_numpy()
    OD_blank = 0.078*1
    A1 = OD_data["A1"].to_numpy() - OD_blank
    A2 = OD_data["A2"].to_numpy() - OD_blank
    A3 = OD_data["A3"].to_numpy() - OD_blank
    B1 = OD_data["B1"].to_numpy() - OD_blank
    B2 = OD_data["B2"].to_numpy() - OD_blank
    B3 = OD_data["B3"].to_numpy() - OD_blank

    C_OD = 1e9
    C1 = A1 * C_OD
    C2 = A2 * C_OD
    C3 = A3 * C_OD
    C4 = B1 * C_OD
    C5 = B2 * C_OD
    C6 = B3 * C_OD
    C = [C1, C2, C3, C4, C5, C6]

    C1_max = np.max(C1)
    C2_max = np.max(C2)
    C3_max = np.max(C3)
    C4_max = np.max(C4)
    C5_max = np.max(C5)
    C6_max = np.max(C6)
    C_max = [C1_max, C2_max, C3_max, C4_max, C5_max, C6_max]

    C1_0 = np.min(C1)
    C2_0 = np.min(C2)
    C3_0 = np.min(C3)
    C4_0 = np.min(C4)
    C5_0 = np.min(C5)
    C6_0 = np.min(C6)
    C_0 = [C1_0, C2_0, C3_0, C4_0, C5_0, C6_0]

    k_gr = [0.01692117, 0.01527194, 0.01514989, 0.01656997, 0.01770599,
        0.01571999]
    k_gr1 = k_gr[0]
    k_gr2 = k_gr[1]
    k_gr3 = k_gr[2]
    k_gr4 = k_gr[3]
    k_gr5 = k_gr[4]
    k_gr6 = k_gr[5]

    C_max_avg = np.mean(C_max[:3])
    C_0_avg = np.mean(C_0[:3])
    k_gr_avg = np.mean(k_gr[:3])

    return [C, C_max, C_0, k_gr, C_max_avg, C_0_avg, k_gr_avg]

def Get_FLOD_Data():
    # Importing experimental data

    geags_data = pd.read_csv("FL_by_OD_expt_data.csv")
    time = geags_data["Time (min)"].to_numpy()
    A1 = geags_data["A1"].to_numpy()
    A2 = geags_data["A2"].to_numpy()
    A3 = geags_data["A3"].to_numpy()
    B1 = geags_data["B1"].to_numpy()
    B2 = geags_data["B2"].to_numpy()
    B3 = geags_data["B3"].to_numpy()

    A12 = A1 - A1[0]
    A1_non_leaky = A12[np.argwhere(A12 >= 0)]
    t12 = time[len(A1_non_leaky) - 1]
    time12 = np.linspace(0,t12,len(A1_non_leaky))

    A22 = A2 - A2[0]
    A2_non_leaky = A22[np.argwhere(A22 >= 0)]
    t22 = time[len(A2_non_leaky) - 1]
    time22 = np.linspace(0,t22,len(A2_non_leaky))

    A32 = A3 - A3[0]
    A3_non_leaky = A32[np.argwhere(A32 >= 0)]
    t32 = time[len(A3_non_leaky) - 1]
    time32 = np.linspace(0,t32,len(A3_non_leaky))

    A_non_leaky = [A1_non_leaky, A2_non_leaky, A3_non_leaky]
    time_A = [time12, time22, time32]

    B12 = B1 - B1[0]
    B1_non_leaky = B12[np.argwhere(B12 >= 0)]
    tB12 = time[len(B1_non_leaky) - 1]
    timeB12 = np.linspace(0,tB12,len(B1_non_leaky))

    B22 = B2 - B2[0]
    B2_non_leaky = B22[np.argwhere(B22 >= 0)]
    tB22 = time[len(B2_non_leaky) - 1]
    timeB22 = np.linspace(0,tB22,len(B2_non_leaky))

    B32 = B3 - B3[0]
    B3_non_leaky = B32[np.argwhere(B32 >= 0)]
    tB32 = time[len(B3_non_leaky) - 1]
    timeB32 = np.linspace(0,tB32,len(B3_non_leaky))

    B_non_leaky = [B1_non_leaky, B2_non_leaky, B3_non_leaky]
    time_B = [timeB12, timeB22, timeB32]

    A_avg = (A1_non_leaky[:88] + A2_non_leaky[:88] + A3_non_leaky[:88])/3
    B_avg = (B1_non_leaky[:88] + B2_non_leaky[:88] + B3_non_leaky[:88])/3
    avg_fold_change = np.mean((np.max(B1_non_leaky)/np.max(A1_non_leaky),np.max(B2_non_leaky)/np.max(A2_non_leaky),np.max(B3_non_leaky)/np.max(A3_non_leaky)))
    
    return [A_non_leaky, time_A, B_non_leaky, time_B, avg_fold_change]