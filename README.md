# GEAGS-Files
All the files required to run the GEAGS CRN model 

BioCRNpyler_Propensities.py: 
This files contains all the propensities that are used in the CRN. BioCRNpyler_Mechanisms.py imports them from here to define the mechanisms 

BioCRNpyler_Mechanisms.py: 
This file contains all the mechanisms created for the CRN. Imports the propensities and exports the mechanism to the main file

BioCRNpyler_Model_Building.ipynb: 
This jupyter notebook is the main file where the model is defined. It imports the mechanisms. This file contains the Mixture definition. The Mixture defines all the species used in the CRN and the mechanisms. Then the file saves the models as SBML files which can be then imported into other files to be modified or to played around with. The file also shows how to import the SBML file and run the model using bioscrape in order to check the model outputs as and when you build or modify the model. 

BioCRNpyler_Model_Comparison.ipynb: 
This file imports the SBML file and runs the model using bioscrape. Made a  separate file for doing this so that the main model defining file is not affected. All the hand fitting exercises are done using files like this where we use LMFIT to define the parametres and their bounds and if they have to be used for fitting using LMFIT in any other file. 

GEAGS_growth_data.csv: 
this file contains the OD data for the experiment 

GEAGS_data_2.csv: 
this file contains the FL/OD data that has been processed 

GEAGS_parameters.txt: 
this file contains all the parameters used in the CRN model with the respective units and their descriptions. This file is used to build the model. The parameter values can be modified using bioscrapes functions and the same is demonstrated in the BioCRNpyler_Model_Comparison.ipynb file. 

GEAGS_param_file_070824.csv: 
This file has all the parameters, including the Initial conditions, as an array, which can be easily loaded into BioCRNpyler_Model_Comparison.ipynb to modify the existing parameters. 
