# GEAGS-Files
All the files required to run the GEAGS CRN model and replicate different results published in the paper. 

## Installation: 
To run the files in the repository, user can clone this repository via git. <br>
Once cloned, we recommend installing the dependencies in a new Virtual Environment as some of the package versions may not be the latest. <br> 
After setting up the environment, user can install the dependencies using the command: pip install -r requirements.txt

## File Information:

-- biocrnpyler_files: The file contains 2 .py files that are used to define the BioCRNpyler propensities and mechanisms. geags_parameters.txt contains all the parameters used in the CRN model, including the respective units and their descriptions. This file is used to build the model. The parameter values can be modified using bioscrapes functions, which are demonstrated in the Run_GEAGS_model.ipynb file.  <br> <br>
-- Build_GEAGS_model.ipynb : This jupyter notebook is the main file where the model is defined. It imports the mechanisms. This file contains the Mixture definition. The Mixture defines all the species used in the CRN and the mechanisms. Then the file saves the models as SBML files, which can be then imported into other files to be modified or played around with. The file also shows how to import the SBML file and run the model using bioscrape to check the model outputs as and when you build or modify the model. <br> <br>
-- Run_GEAGS_model.ipynb: This file imports the SBML file and runs the model using bioscrape. I made a separate file for doing this so that the main model-defining file is not affected. All the hand fitting exercises are done using files like this, where we use LMFIT to define the parameters and their bounds and if they have to be used for fitting using LMFIT in any other file. <br> <br>
-- experiment_data: Contains the data files used to plot the experimental data as well as the utility .py file used to import the data from the csv files.  <br> <br>
-- model_param_file_030525.csv:  This file has all the parameters, including the Initial conditions, as an array, which can be easily loaded into Run_GEAGS_model.ipynb to modify the existing parameters.
