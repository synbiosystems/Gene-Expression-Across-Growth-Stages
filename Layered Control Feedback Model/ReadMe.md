## This folder contains all the file required to simulate a Layered Feedback Control circuit applying the GEAGS framework

Here, we simulated the dynamics of 4 constructs with different feedback strategies as described in the publication: https://doi.org/10.1038/s41467-022-33058-6 ('Layered feedback control overcomes performance trade-off in synthetic biomolecular networks') <br> <br>

File names starting with paper_ contain the experimental data obtained from the paper. However, the file paper_controller_profile.txt contains the simulation data from a model that was described in the paper. <br> <br>

LCM_model_equations.py contains the ODEs and the ODE solver required to run the Layered Control Model (LCM). LCM_implementation_main.ipynb demonstrates how the model was used to explain the dynamics of the 4 constructs. LCM)param_file_042025.csv contains the model parameters
