## Updated July 2025. 2025_feb_25
## author: amouzgar@stanford.edu
## Nature Communications review: please do not distribute

list.of.packages <- c("tidyverse", "pROC","nnet","caret",'patchwork',"rstudioapi",'dynutils','pheatmap')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("variancePartition")

#############################################################################
## script summary
#############################################################################
#' This script provides example code to run: 
#' (1) Multinomial Log-linear Models and multi-class ROC for cell line predictions using the minimal, core, and larger CC panels
#' (2) Detect feature-level differences between cell lines using anova to calculate variance explained by cell lines
#' (3) Example code to partition variance using main effects and a two-way interaction
#' (4) Code to compute an aberrancy score using nearest-neighbor analysis
#' (5) Dimensionality reduction (umap)
#' (6) Raw code to run single-cell perturbation scores (requires Reticulate)

########################################################################################################
## load packages and variables
########################################################################################################
set.seed(0)  # For reproducibility
library(tidyverse)
library(caret)

## cell line varaibles
cell_line_levels =  c('NALM6','U937','HEL','293T','JURKAT')


## cell cycle molecules, core panel
features_cellcycle = c('Ki67','pRbS780', 'pH3_s10','CDT1', 'IdU', 'Geminin', 'PLK1', 'DNA',  'CyclinB1', 'PCNA', 'SLBP')%>% make.names()
## cell cycle molecules, exmaple additional molecules
features_other = c( 'EZH2','HH3','CyclinE','FoxM1','CDC20','RUNX1','CTCF','WGA','MLL1')
## combined core and additional molecules
features_cellcycle_full = c(features_cellcycle, features_other)

## examples for minimal panel, core panel, and an example larger panel
features_cellcycle_minimal = features_cellcycle %>% .[!. %in% c('IdU','DNA')]
features_cellcycle_core = features_cellcycle
features_cellcycle_curated = features_cellcycle_full

## read example data
filepath=dirname(rstudioapi::getSourceEditorContext()$path)
cell_line_df =read.csv(paste0(filepath, '/2025_feb-example_data.csv'),sep=',')

####################################################
## function to extract ROC curve infomration 
####################################################
##' extract_multiclass_roc_info:  builds multiple ROC curve to compute the multi-class AUC as defined by Hand and Till
##' This function will extract relevant ROC information to construct ROC curves.
extract_multiclass_roc_info = function(data, predictions_obj){
  lapply(1:ncol(predictions_obj), function(predict_res_probs_NUMBER){
    ## identify cell line of interest
    predict_res_probs_NUMBER_CELL_LINE = colnames(predictions_obj)[predict_res_probs_NUMBER]
    message(predict_res_probs_NUMBER)
    
    ## calculate multiclass ROC
    roc_data_complete <- pROC::multiclass.roc(data %>% .$cell_line, predictions_obj[,predict_res_probs_NUMBER])
    
    ## returns the coordinates of the ROC curve between 0 and 100% for each cell line comparison
    roc_coords_complete = lapply(roc_data_complete$rocs, function(roc_data_1){
      roc_coords= pROC::coords(roc_data_1, seq(0,1,0.01)) %>%
        mutate(group1 = roc_data_1$levels[1],
               group2 = roc_data_1$levels[2])
    }) %>% bind_rows()
    roc_coords_complete$TASK = predict_res_probs_NUMBER_CELL_LINE
    return(roc_coords_complete)
  }) %>% bind_rows()%>%
    ## add specificity and comparison metadata
    mutate(specificity_min1 = 1-specificity,
           comparison = paste(group1, group2, sep = '-vs-'))
}

########################################################################################################
## train models and extract ROC information for cell line comparisons
########################################################################################################

data = cell_line_df %>% dplyr::select(all_of(features_cellcycle_full)) ## feature-only matrix

## perform train-test split using caret
trainIndex <- caret::createDataPartition(cell_line_df$cell_line, p = .5, list = FALSE, times = 1)
trainSet <- data[trainIndex, ]
testSet <- data[-trainIndex, ]

# Scale the training data
preProc <- preProcess(trainSet, method = c("center", "scale"))
trainSet_scaled <- predict(preProc, trainSet)
trainSet_scaled <- cbind(trainSet_scaled, cell_line =  cell_line_df[trainIndex, ]$cell_line)

# Scale the test data using the same parameters as the training data
testSet_scaled <- predict(preProc, testSet)
testSet_scaled <- cbind(testSet_scaled, cell_line =  cell_line_df[-trainIndex, ]$cell_line)

# View the scaled datasets
head(trainSet_scaled)
head(testSet_scaled)

###########################################################################################################################
## Fit Multinomial Log-linear Models to compare cell lines using the different feature sets defined in variables section ##
###########################################################################################################################

#' (1) use FULL feature set
multinom_res = nnet::multinom(formula =  as.formula(paste0('cell_line~0+',paste(features_cellcycle_curated, collapse ='+'))), data = trainSet_scaled)
predict_res_probs = predict(multinom_res, testSet_scaled,type = 'probs')
multiclass_roc_full = extract_multiclass_roc_info(data = testSet_scaled,predictions_obj = predict_res_probs)
multiclass_roc_full$analysis = 'complete' ## add some metadata

## predict on test set
predict_res = predict(multinom_res, testSet_scaled)
conf_mat_results_complete = caret::confusionMatrix(table(predict_res, testSet_scaled %>% .$cell_line))

#' (2) use core feature set
multinom_res = nnet::multinom(formula =  as.formula(paste0('cell_line~0+',paste(features_cellcycle_core, collapse ='+'))), data = trainSet_scaled)
predict_res_probs = predict(multinom_res, testSet_scaled,type = 'probs')
multiclass_roc_core = extract_multiclass_roc_info(data = testSet_scaled,predictions_obj = predict_res_probs)
multiclass_roc_core$analysis = 'core' ## add some metadata

## predict on test set
predict_res = predict(multinom_res, testSet_scaled)
conf_mat_results_core = caret::confusionMatrix(table(predict_res, testSet_scaled %>% .$cell_line))

#' (3) use minimal feature set
multinom_res = nnet::multinom(formula =  as.formula(paste0('cell_line~0+',paste(features_cellcycle_minimal, collapse ='+'))), data = trainSet_scaled)
predict_res_probs = predict(multinom_res, testSet_scaled,type = 'probs')
multiclass_roc_minimal = extract_multiclass_roc_info(data = testSet_scaled,predictions_obj = predict_res_probs)
multiclass_roc_minimal$analysis = 'core, no DNA' ## add some metadata

## predict on test set
predict_res = predict(multinom_res, testSet_scaled)
conf_mat_results_minimal = caret::confusionMatrix(table(predict_res, testSet_scaled %>% .$cell_line))


###################################################
## combine Multinomial Log-linear Models results ##
###################################################
multiclass_roc_relevant_task = multiclass_roc_full%>%
  bind_rows(multiclass_roc_core,multiclass_roc_minimal) %>%
  dplyr::filter(TASK == group1 | TASK == group2) 

## plot ROC results
p1 = ggplot(multiclass_roc_relevant_task, aes(x=specificity_min1, y = sensitivity, color = comparison ))+
  theme_bw() +
  theme(text = element_text(size = 6,family = 'sans')) +
  geom_abline(1,intercept = 0 , size = 0.5, color  ='darkgray',  linetype = 'dashed')+
  geom_line(aes(group =paste0(TASK, comparison)), linewidth = 0.5) +
  xlab('1-specificity') +
  ylab('sensitivity') +
  facet_wrap(~analysis)
p1

################################################
## sensitivity, specificity, accuracy summary ##
################################################
## compute sensitivity, specificity, and balanced accuracy by class for the different panel subsets
conf_mat_results = conf_mat_results_complete$byClass %>%
  data.frame(check.rows = F, check.names = F) %>%
  mutate(analysis = 'complete',
         cell_line = rownames(.)) %>%
  bind_rows(conf_mat_results_core$byClass  %>%
              data.frame(check.rows = F, check.names = F)%>%
              mutate(analysis = 'core',
                     cell_line = rownames(.))) %>%
  bind_rows(conf_mat_results_minimal$byClass  %>%
              data.frame(check.rows = F, check.names = F)%>%
              mutate(analysis = 'minimal',
                     cell_line = rownames(.))) %>%
  mutate(cell_line = gsub('Class: ','',cell_line))


## data wrangling for each metric
conf_mat_results_plot =conf_mat_results%>%
  dplyr::select(Sensitivity, Specificity, cell_line, analysis, `Balanced Accuracy`) %>%
  gather(key = metric, value =value, -cell_line,-analysis ) %>%
  mutate(metric =factor(metric, levels = c('Sensitivity','Specificity','Balanced Accuracy')),
         cell_line = factor(cell_line, levels =cell_line_levels))

## plot the sensitivity, specificity, and balanced accuracy by class results
p2 =  ggplot(conf_mat_results_plot, aes(x = analysis, y = value))+
  theme_bw() +
  theme(text = element_text(size = 6,family = 'sans')) +
  geom_boxplot(size = 0.25) +
  geom_point(aes(color =cell_line), size = 0.25) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_brewer(palette = 'Dark2') +
  facet_wrap(~metric) 
p2

#########################################################################################
## Detect feature-level differences between cell lines using anova:
#########################################################################################
## anova analysis: fit linear model using stats::lm and then
## compute analysis of variance (anova) on the fitted lm object for every feature 
#########################################################################################
## get anova table information for each feature, which is input into calc_var_explained to calculate the proportion of variance by cell line for each feature
my_anovas=lapply(features_cellcycle_full,function(feat){
  # glm_res = glm(formula =  as.formula(paste0(feat, '~0+cell_line')), data = glm_input)
  glm_res = lm(formula =  as.formula(paste0(feat, '~0+ cell_line')), data = testSet_scaled) ## construct linear model
  anova_res = anova(glm_res)  ## anova calculation
  anova_res$feature = feat
  return(anova_res)
})



## function to calculate the proportion of variance by cell line for each feature
calc_var_explained = function(data, anova_table, feature){
  # Get the total sum of squares
  total_ss <- sum((data[[feature]] - mean(data[[feature]]))^2)
  # Extract the Sum of Squares (SS) for the categorical variable
  ss_category <- anova_table$`Sum Sq`[1]
  # Calculate the proportion of variance explained by the categorical variable
  variance_explained <- ss_category / total_ss
  return(variance_explained*100)
}

## begin calculating variance explained for each feature using the calc_var_explained function
my_var_explained= lapply(my_anovas, function(anova_table){
  calc_var_explained(testSet_scaled, anova_table, unique(anova_table$feature)) %>% data.frame(var_expl=.)%>% mutate(feature = unique(anova_table$feature))
}) %>% bind_rows()

## plot the results
p3 = ggplot(my_var_explained, aes(x= var_expl, y = reorder(feature, var_expl))) +
  theme_bw()+
  theme(text = element_text(size = 6,family = 'sans')) +
  geom_col()

## plot results, expected to be a little different from original manuscript due to much smaller data sample size but general trends are the same
patchwork::wrap_plots(p1,p2,p3)

############################################################################################
## (3) Example code to peform variance partition analysis ##
############################################################################################

varpart_input = cell_line_df  %>% na.omit()%>% group_by(cell_line)%>% sample_n(500)%>% ungroup() 
varpart_input_mat = varpart_input %>% dplyr::select(all_of(features_cellcycle_full), cell.id)%>% column_to_rownames('cell.id')%>% t()  ## sc matrix
varpart_input_md = varpart_input %>% dplyr::select(cell.id, gate, cell_line) %>% column_to_rownames('cell.id') ## sc metadata

form <- ~  (1 | cell_line)  + (1 | gate) + (1 | cell_line:gate)## random effect models for cell line and gate and their interaction as a two-way interaction example

varPart <- variancePartition::fitExtractVarPartModel(varpart_input_mat, form, varpart_input_md)
vp <- variancePartition::sortCols(varPart)
variancePartition::plotPercentBars(vp)
variancePartition::plotVarPart(vp)

############################################################################################
## (4) Code to compute an unsupervised aberrancy estimate using nearest-neighbor analysis ##
############################################################################################
## We provide example code to compute an unsupervised aberrancy score. 
## The threshold for assigning aberrancy is up to the user.
##' full_data: single-cell dataframe
##' LM_col: column name (string) containing your identifier for your group of interest (eg, "treatment_group") as the landmarks. 
##' LM_group: column name (string) containing your landmark group (eg, "NoTx") 
##' features: vector of features to use in your analysis c('feature1','feature2','feature3'). The columns must be in the full_data  dataframe that you provide. 
##' k_ave: number of neighbors to use
##' return_knn_graph: Default is FALSE. Returns the knn graph.
##' return_nearest_cellid: Default is FALSE. Returns the nearest cell id.
##' 
## load knn aberrancy function
source(paste0(filepath, '/knn_aberrancy_fxn.R'))

## as a toy example, we use JURKAT cells as landmarks for all cell lines.
knn_res = compute_knn_aberrancy(full_data = cell_line_df,LM_col = 'cell_line',LM_group = c('JURKAT'),k_ave = 50,features = features_cellcycle,return_knn_graph = F,return_nearest_cellid = F ) 
table(cell_line_df$cell.id == names(knn_res))## output order should match the original order of the dataframe
cell_line_df$aberrancy_score = knn_res ## insert the score as a new column

ggplot(cell_line_df, aes(x=aberrancy_score, color = cell_line))+
  theme_bw()+
  geom_density()

############################################################################################
## (5) Dimensional reduction ##
############################################################################################
umap_input = cell_line_df %>% ungroup()%>% mutate_at(all_of(features_cellcycle), scale) %>%  group_by(cell_line) %>%sample_n(1500) %>% ungroup()
umap_res=uwot::umap(umap_input %>%  dplyr::select(all_of(features_cellcycle)) ,spread = 6,n_neighbors = 10, verbose =T)
umap_df = data.frame(umap_res) %>% ## extract coordinates
  bind_cols(umap_input)

ggplot(umap_df, aes(x=X1, y=X2, color = cell_line))+
  geom_point()


############################################################################################
## (6) Raw code to run single-cell perturbation scores (requires Reticulate) ##
############################################################################################
## This code will calculate single-cell perturbation scores using Meld. 
## Please see the original package details for how to install it.
reticulate::use_condaenv(condaenv = 'meld') ## note: must setup your environment 

## compute meld on cell line identities as an example
## we normalize and downsample the data
meld_input = cell_line_df %>% ungroup()%>% mutate_at(all_of(features_cellcycle), dynutils::scale_minmax) %>%  group_by(cell_line) %>%sample_n(1000) %>% ungroup()

## this function uses reticulate to calculate meld via python. 
## you must configure your environment yourself for this to work.
## please see the original manuscript for installation. Burkhardt et al. 2021
calculate_meld = function (df_features_only, labels) {
  meld_likelihood_score <- reticulate::import("meld")
  np <- reticulate::import("numpy")
  pd <- reticulate::import("pandas")
  pd_df = pd$pandas$DataFrame(df_features_only)
  sample_densities = meld_likelihood_score$meld$MELD()$fit_transform(X = pd_df, 
                                                                     sample_labels = np$array(labels))
  sample_likelihoods = meld_likelihood_score$meld$utils$normalize_densities(sample_densities)
  return(sample_likelihoods)
}

meld_score = calculate_meld(meld_input %>% dplyr::select(all_of(features_cellcycle)),labels =  meld_input$cell_line )

## bind the meld scores to your data
meld_output =meld_input %>%
  bind_cols(meld_score)

## calculate and plot correlations between cell line scores
pheatmap::pheatmap(cor(meld_score))








