# 2025_cellcycle_ml
Example code with test data for running ML algorithms of "A deep single cell mass cytometry approach to capture canonical and noncanonical cell cycle states". 


MELD is a python package with installation instruction available at https://github.com/KrishnaswamyLab/MELD
Reticulate and conda can be used to run MELD in R. Information on installing reticulate can be found at https://rstudio.github.io/reticulate/ and conda at https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html 
After correct installations, Reticulate in R can be used to load the newly created conda environment and run python packages like MELD directly through R. 


After users have correctly installed the software, we have also provided an example starter function in the Github page to help users run MELD in R:


```
meld_using_reticulate_example_starter_function = function(df_features_only, labels_vector){
     meld_likelihood_score <- reticulate::import("meld")
     np <- reticulate::import("numpy")
     pd <- reticulate::import("pandas")
     pd_df = pd$pandas$DataFrame(df_features_only)
     sample_densities = meld_likelihood_score$meld$MELD()$fit_transform(X = pd_df, 
     sample_labels = np$array(labels_vector))
     sample_likelihoods = meld_likelihood_score$meld$utils$normalize_densities(sample_densities)
     return(sample_likelihoods)
)
```

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
Simply download the data and script into the same directory and run the script. Code will run directly on the data provided in this github repo. Can be easily reproduced using other data by just changing the feature variables, class label (cell_line in this case), and input data.



