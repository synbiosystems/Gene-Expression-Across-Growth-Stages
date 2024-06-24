# GEAGS-Files
All the files required to run the GEAGS CRN model 

GEAGS_propensities_4.py: 
This files contains all the propensities that are used in the CRN. GEAGS_mechanisms_4 imports them from here to define the mechanisms 

GEAGS_mechanisms_4.py: 
This file contains all the mechanisms created for the CRN. Imports the propensities and exports the mechanism to the main file

GEAGS_model_biocrnpyler_4.ipynb: 
This jupyter notebook is the main file where the model is defined. It imports the mechanisms. This file contains the Mixture definition. The Mixture defines all the species used in the CRN and the mechanisms. Then the file saves the models as SBML files which can be then imported into other files to be modified or to play around with

GEAGS_model_comparison_biocrnpyler_clean.ipynb: 
This file imports the SBML file and runs the model using bioscrape. Made a  separate file for doing this so that the main model defining file is not affected. All the hand fitting exercises are done using files like this where we use LMFIT to define the parametres and their bounds and if they have to be used for fitting using LMFIT in any other file. 

GEAGS_growth_data.csv: 
this file contains the OD data for the experiment 

GEAGS_data_2.csv: 
this file contains the FL/OD data that has been processed 

GEAGS_parameters.txt: 
this file contains all the parameters used in the CRN model 

