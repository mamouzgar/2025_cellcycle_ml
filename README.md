# 2025_cellcycle_ml
Example code with test data for running ML algorithms of "A deep single cell mass cytometry approach to capture canonical and noncanonical cell cycle states". 

System requirements:
tidyverse 2.0.0
pROC 1.18.5
nnet 7.3.19
caret 7.0.1
patchwork 1.3.0
rstudioapi 0.17.1
magrittr 2.0.3

Installation:
 The script will detect and install missing packages. Install time for packages will vary with user but should be <1-3 minutes (performed on a MacOS 2.4 GHz 8-Core Intel Core i9).

Demo:
Expected final output are the plots used in the figure with general trends similar to the primary manuscript. Minor differences are expected due to downsampling.

Instructions for use:
Simply download the data and script into the same directory and run the script. Code will run directly on the data provided in this github rep. Can be easily reproduced using other data by just changing the feature variables, class label (cell_line in this case), and input data.


