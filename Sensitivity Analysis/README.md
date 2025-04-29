## This folder contains the files used for sensitivity analysis 

The files Sensitivity_analysis_degtag.ipynb and Sensitivity_analysis_no_degtag.ipynb are two jupyter notebooks used for generating the sensitivity matrices (SSMs) using the Bioscrape toolbox  <br> 
Created two separate files so that they can be run in parallel, as the estimation of the sensitivity matrix can be computationally expensive depending on the simulation step size. <br>
These notebooks export the SSMs as a numpy datafile  <br>
The file Sensitivity_matrix_analysis.ipynb shows how to extract data from the datafile and export it for plotting purposes   <br>
The simulation data files are used to normalize the SSM matrix <br>
For more information on the sensitivity analysis and the toolbox, you can read about the tool at :   <br> https://github.com/biocircuits/bioscrape/blob/master/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb
