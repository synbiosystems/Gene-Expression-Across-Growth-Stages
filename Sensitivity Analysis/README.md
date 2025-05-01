## This folder contains the files used for sensitivity analysis 

The files Sensitivity_analysis_degtag.ipynb and Sensitivity_analysis_no_degtag.ipynb are two jupyter notebooks used for generating the sensitivity matrices (SSMs) using the Bioscrape toolbox  <br> 
Created two separate files so that they can be run in parallel, as the estimation of the sensitivity matrix can be computationally expensive depending on the simulation step size. Running each file takes about 90 minutes on a Mac Mini M2 Pro <br>
These notebooks export the SSMs as a numpy datafile. The numpy datafile used in the paper figure is large and GitHub doesn't support the upload of files of that size, so please refer to the data provided with the pubilcation if needed.  <br>
The file Sensitivity_matrix_data_export_utility.ipynb shows how to extract data from the datafile and export it for plotting purposes   <br>
The simulation data files are used to normalize the parameters modified by RMFs in the SSM matrix <br>
For more information on the sensitivity analysis and the toolbox, you can read about the tool at:   <br> https://github.com/biocircuits/bioscrape/blob/master/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb
