# SWTH_translocation
Code for analysis of Swainson's thrush translocation experiment

### 01_translocation_dataFormatting.Rmd

* formatting survival and phenotypic data from detections
* basic data checks and visualization
* making maps of detection data (Fig 3A, Fig 4)

### 02_translocation_analysis_orientation.Rmd

* analysis of migratory orientations
* data checks and visualization
* plotting analysis outcome (Fig. 3C)

### 03_translocation_analysis_survival.R

* CJS model implementation
* uses output from 01_translocation_dataFormatting.Rmd

### 04_translocation_analysis_survival.Rmd

* visualizing output of CJS models (Fig. 2)
* formatting CJS parameter and estimate tables