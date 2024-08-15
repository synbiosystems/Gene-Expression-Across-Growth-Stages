## Function used to set the model parameters and initial conditions

# Import here 

from bioscrape.types import Model
from bioscrape.simulator import py_simulate_model

def run_CRN_model(param, model_degtag, model_no_tag, timepoints): 
    
    params = param.to_dict()

    Ribo_min = params['Ribo_min']

    tRNA_min = params['tRNA_min']
    
    RNAP_min = params['RNAP_min']

    Et_min = params['Et_min']

    protease_min = params['protease_min']
    
    n = 5.5

    model_degtag.set_species({'dna_X': 100,
    'protein_RNAP_machinery': params['RNAP_IC'],
    'protein_sigma_machinery': params['Sigma_IC'],
    'protein_NT_units': params['NT_IC'],
    'protein_Ribo_machinery': params['Ribo_IC'],
    'protein_AA_units0tl': params['AA_IC'],
    'protein_tRNA_machinery': params['tRNA_IC'],
    'protein_RNAase_machinery': params['RNAase_IC'],
    'protein_Et_machinery': params['Et_IC'],
    'protein_protease_machinery': params['protease_IC'],
    'protein_unfolded_protein_degtag': 0,
    'rna_T_X':0,
    'protein_X':0,
    'cell_count_count': 6666666.67})

    model_no_tag.set_species({'dna_X': 100,
    'protein_RNAP_machinery': params['RNAP_IC'],
    'protein_sigma_machinery': params['Sigma_IC'],
    'protein_NT_units': params['NT_IC'],
    'protein_Ribo_machinery': params['Ribo_IC'],
    'protein_AA_units0tl': params['AA_IC'],
    'protein_tRNA_machinery': params['tRNA_IC'],
    'protein_RNAase_machinery': params['RNAase_IC'],
    'protein_Et_machinery': params['Et_IC'],
    'protein_protease_machinery': params['protease_IC'],
    'protein_unfolded_protein_degtag': 0,
    'rna_T_X':0,
    'protein_X':0,
    'cell_count_count': 6666666.67})

    model_degtag.set_params({'k__logistic_cell_growth': params['k_gr'], 
    'c_max__logistic_cell_growth': params['C_max'],
    'k_rnap__bacterial_transcription':params['k_RNAP'],
    'rnap_min__bacterial_transcription': RNAP_min,
    'k_ribo__bacterial_translation':params['k_Ribo'],
    'ribo_min__bacterial_translation': Ribo_min,
    'k_tRNA__bacterial_translation':params['k_tRNA'],
    'tRNA_min__bacterial_translation': tRNA_min,
    'k_Et__bacterial_translation':params['k_Et'],
    'Et_min__bacterial_translation': Et_min,
    'k_protease__bacterial_translation':params['k_protease'],
    'protease_min__bacterial_translation': protease_min,

    'k_tx_1b__bacterial_transcription':params['k_tx_1b'],
    'k_tx_1u__bacterial_transcription':params['k_tx_1u'],
    'k_tx_2u__bacterial_transcription':params['k_tx_2u'],
    'k_tx_2b__bacterial_transcription':params['k_tx_2b'],
    'k_tx_3__bacterial_transcription': params['k_tx_3'],
    'k_tx_4b__':params['k_tx_4b'],
    'k_tx_4u__':params['k_tx_4u'],
    'k_tx_5__': params['k_tx_5'],

    'k_tl_1b__bacterial_translation':params['k_tl_1b'],
    'k_tl_1u__bacterial_translation':params['k_tl_1u'],
    'k_tl_2__bacterial_translation':params['k_tl_2'],
    'k_tl_3__bacterial_translation':params['k_tl_3'],
    'k_tl_4__bacterial_translation':params['k_tl_4'],
    'k_tl_5__bacterial_translation':params['k_tl_5'],
    'k_tl_6b__bacterial_translation':params['k_tl_6b'],
    'k_tl_6u__bacterial_translation':params['k_tl_6u'],
    'k_tl_7__bacterial_translation':params['k_tl_7'],
    'b_tl_7__bacterial_translation':params['b_tl_7'],
    'k_tl_8__bacterial_translation':params['k_tl_8'],
    'P_max__bacterial_translation':params['P_max'],
    'b_tl_8__bacterial_translation':params['b_tl_8'],
    'k_tl_9__': params['k_tl_9'],
    'k_tl_10__': params['k_tl_10'],
    'k_tl_11__': params['k_tl_11'],
    'k_tl_12__bacterial_translation': params['k_tl_12'],
    'k_tl_13b__':params['k_tl_13b'],
    'k_tl_13u__':params['k_tl_13u'],
    'n__': n, 
    'n__mrna_degradation': n})

    model_no_tag.set_params({'k__logistic_cell_growth': params['k_gr'], 
    'c_max__logistic_cell_growth': params['C_max'],
    'k_rnap__bacterial_transcription':params['k_RNAP'],
    'rnap_min__bacterial_transcription': RNAP_min,
    'k_ribo__bacterial_translation':params['k_Ribo'],
    'ribo_min__bacterial_translation': Ribo_min,
    'k_tRNA__bacterial_translation':params['k_tRNA'],
    'tRNA_min__bacterial_translation': tRNA_min,
    'k_Et__bacterial_translation':params['k_Et'],
    'Et_min__bacterial_translation': Et_min,

    'k_tx_1b__bacterial_transcription':params['k_tx_1b'],
    'k_tx_1u__bacterial_transcription':params['k_tx_1u'],
    'k_tx_2u__bacterial_transcription':params['k_tx_2u'],
    'k_tx_2b__bacterial_transcription':params['k_tx_2b'],
    'k_tx_3__bacterial_transcription': params['k_tx_3'],
    'k_tx_4b__':params['k_tx_4b'],
    'k_tx_4u__':params['k_tx_4u'],
    'k_tx_5__': params['k_tx_5'],

    'k_tl_1b__bacterial_translation':params['k_tl_1b'],
    'k_tl_1u__bacterial_translation':params['k_tl_1u'],
    'k_tl_2__bacterial_translation':params['k_tl_2'],
    'k_tl_3__bacterial_translation':params['k_tl_3'],
    'k_tl_4__bacterial_translation':params['k_tl_4'],
    'k_tl_5__bacterial_translation':params['k_tl_5'],
    'k_tl_6b__bacterial_translation':params['k_tl_6b'],
    'k_tl_6u__bacterial_translation':params['k_tl_6u'],
    'k_tl_7__bacterial_translation':params['k_tl_7'],
    'b_tl_7__bacterial_translation':params['b_tl_7'],
    'k_tl_8__bacterial_translation':params['k_tl_8'],
    'P_max__bacterial_translation':params['P_max'],
    'b_tl_8__bacterial_translation':params['b_tl_8'],
    #'k_tl_9__': params['k_tl_9'],
    'k_tl_10__': params['k_tl_10'],
    'k_tl_11__': params['k_tl_11'],
    'k_tl_12__bacterial_translation': params['k_tl_12'],
    'n__': n, 
    'n__mrna_degradation': n})

    sol_deg = py_simulate_model(timepoints = timepoints, Model = model_degtag)
    sol_no_deg = py_simulate_model(timepoints = timepoints, Model = model_no_tag)

    return sol_deg, sol_no_deg