## Changelog

### 1.1.0
* FEAT: Integrate the AlphaPeptDeep package into AlphaViz to predict the peptide properties, such as the MS2 intensity prediction for DDA data and retention time/ion mobility prediction for DIA data.
* FEAT: Update the Quality Control tab to make it lighter than the previous version. Variables for plots (m/z, im, rt, length, etc. distributions) can now be selected.
* FEAT: Reorganize the Settings panel and add more customization options.
* DOCS: Extend Jupyter notebooks for manual visual inspection of timsTOF DDA and DIA data analyzed by MaxQuant and DIA-NN, respectively.


### 1.0.5

* FEAT: Enable the possibility to read DIA data analysed by the DIA-NN software analysis tool.
* FEAT: Extend the "Quality Control" tab with a summary statistics table and additional plots.
* FEAT: Add a "Targeted Mode" tab to work with raw data only.
* DOCS: Update the AlphaViz tutorial with all functionality implemented.
* FEAT: Create detailed notebooks demonstrating visualisation of DDA and DIA data analysed with MaxQuant and DIA-NN respectively.
* FIX: Correct the error when reading the MQ output files (missing columns).

### 0.0.1

* FEAT: Initial creation of alphaviz.
* FEAT: Read DDA Bruker data processed with MaxQuant software.
