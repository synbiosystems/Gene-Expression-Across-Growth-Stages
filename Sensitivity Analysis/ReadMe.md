### File Information: 

-- __Sensitivity_analysis_degtag.ipynb__ and __Sensitivity_analysis_no_degtag.ipynb__ : Jupyter notebooks used for generating the sensitivity matrices (SSMs) using the Bioscrape toolbox  <br> 
Created two separate files so that they can be run in parallel, as the estimation of the sensitivity matrix can be computationally expensive depending on the simulation step size. Running each file takes about 90 minutes on a Mac Mini M2 Pro <br>
These notebooks export the SSMs as a numpy datafile. The numpy datafile used in the paper figure is large and GitHub doesn't support the upload of files of that size, so please refer to the data provided with the pubilcation if needed.  <br> <br>
-- __Sensitivity_matrix_data_export_utility.ipynb__ : Notebook to extract data from the datafiles and export it for plotting purposes  <br> <br>
-- __simulation_data__ : Contains the simulation data files are used to normalize the parameters modified by RMFs in the SSM matrix. Also, this is the folder where the SSM datafile is exported to. <br> <br>
#### __IMPORTANT__ : 
The implementation of sensitivity analysis for GEAGS model has a modification in the algorithm and the modified sensitivity analysis can be found in : https://github.com/hariKRN2000/bioscrape <br> 
Make sure the user has the above version of bioscrape installed. The requirements.txt file already contains a git implementation of installation of bioscrape from the above repository. <br> <br>
For more information on the sensitivity analysis and the toolbox, you can read about the tool at:   <br> https://github.com/biocircuits/bioscrape/blob/master/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb
