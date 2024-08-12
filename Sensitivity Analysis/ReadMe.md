## This folder contains the files used for sensitivity analysis 

The files Sensitivity_analysis_degtag.ipynb and Sensitivity_analysis_no_degtag.ipynb are two jupyter notebooks used for generating the sensitivity matrices (SSMs) using the Bioscrape toolbox  <br> 
Used two separate notebooks to generate the two sensitivity matrices in parallel because running the sensitivity analysis function can be computationally expensive, so ran the two functions in parallel in two notebooks. <br>
These notebooks export the SSMs as a numpy datafile  <br>
The file Sensitivity_matrix_analysis.ipynb shows how to extract data from the datafile and export it for plotting purposes   <br>
For more information on the sensitivity analysis and the toolbox, you can read about the tool at :   <br> https://github.com/biocircuits/bioscrape/blob/master/examples/Sensitivity%20Analysis%20using%20Bioscrape.ipynb
