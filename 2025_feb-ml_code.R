## 2025_feb_25
## author: amouzgar@stanford.edu
## Nature Communications review: please do not distribute

list.of.packages <- c("tidyverse", "pROC","nnet","caret",'patchwork',"rstudioapi")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


########################################################################################################
## load packages and variables
########################################################################################################
set.seed(0)  # For reproducibility
library(tidyverse)
library(caret)

cell_line_levels =  c('NALM6','U937','HEL','293T','JURKAT')

features_cellcycle = c('Ki67','pRbS780', 'pH3_s10','CDT1', 'IdU', 'Geminin', 'PLK1', 'DNA',  'CyclinB1', 'PCNA', 'SLBP')%>% make.names()
features_other = c( 'EZH2','HH3','CyclinE','FoxM1','CDC20','RUNX1','CTCF','WGA','MLL1')
features_cellcycle_full = c(features_cellcycle, features_other)

features_cellcycle_curated = features_cellcycle_full
features_cellcycle_core = features_cellcycle
features_cellcycle_minimal = features_cellcycle %>% .[!. %in% c('IdU','DNA')]

## read example data
filepath=dirname(rstudioapi::getSourceEditorContext()$path)
cell_line_df =read.csv(paste0(output_ml, '/2025_feb-example_data.csv'),sep=',')

####################################################
## function to extract ROC curve infomration 
####################################################
extract_multiclass_roc_info = function(data, predictions_obj){
  lapply(1:ncol(predictions_obj), function(predict_res_probs_NUMBER){
    predict_res_probs_NUMBER_CELL_LINE = colnames(predictions_obj)[predict_res_probs_NUMBER]
    message(predict_res_probs_NUMBER)
    roc_data_complete <- pROC::multiclass.roc(data %>% .$cell_line, predictions_obj[,predict_res_probs_NUMBER])
    
    roc_coords_complete = lapply(roc_data_complete$rocs, function(roc_data_1){
      roc_coords= pROC::coords(roc_data_1, seq(0,1,0.01)) %>%
        mutate(group1 = roc_data_1$levels[1],
               group2 = roc_data_1$levels[2])
    }) %>% bind_rows()
    roc_coords_complete$TASK = predict_res_probs_NUMBER_CELL_LINE
    return(roc_coords_complete)
  }) %>% bind_rows()%>%
    mutate(specificity_min1 = 1-specificity,
           comparison = paste(group1, group2, sep = '-vs-'))
}

####################################################
## extract ROC information 
####################################################

data = cell_line_df %>%
  dplyr::select(all_of(features_cellcycle_full))
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


trainSet_scaled =  CELL_LINE_TRAJECTORY_INFERENCE_df%>%
  ungroup() %>%
  mutate_at(features_cellcycle_curated, scale) %>%
  dplyr::filter(cell.id %in% train_ids)
testSet_scaled = CELL_LINE_TRAJECTORY_INFERENCE_df%>%
  ungroup() %>%
  mutate_at(features_cellcycle_curated, scale) %>%
  dplyr::filter(!cell.id %in% train_ids)
## FULL feature set
multinom_res = nnet::multinom(formula =  as.formula(paste0('cell_line~0+',paste(features_cellcycle_curated, collapse ='+'))), data = trainSet_scaled)
predict_res_probs = predict(multinom_res, testSet_scaled,type = 'probs')
multiclass_roc_full = extract_multiclass_roc_info(data = testSet_scaled,predictions_obj = predict_res_probs)
multiclass_roc_full$analysis = 'complete'
predict_res = predict(multinom_res, testSet_scaled)
conf_mat_results_complete = caret::confusionMatrix(table(predict_res, testSet_scaled %>% .$cell_line))

## core feature set
multinom_res = nnet::multinom(formula =  as.formula(paste0('cell_line~0+',paste(features_cellcycle_core, collapse ='+'))), data = trainSet_scaled)
predict_res_probs = predict(multinom_res, testSet_scaled,type = 'probs')
multiclass_roc_core = extract_multiclass_roc_info(data = testSet_scaled,predictions_obj = predict_res_probs)
multiclass_roc_core$analysis = 'core'
predict_res = predict(multinom_res, testSet_scaled)
conf_mat_results_core = caret::confusionMatrix(table(predict_res, testSet_scaled %>% .$cell_line))

## minimal feature set
multinom_res = nnet::multinom(formula =  as.formula(paste0('cell_line~0+',paste(features_cellcycle_minimal, collapse ='+'))), data = trainSet_scaled)
predict_res_probs = predict(multinom_res, testSet_scaled,type = 'probs')
multiclass_roc_minimal = extract_multiclass_roc_info(data = testSet_scaled,predictions_obj = predict_res_probs)
multiclass_roc_minimal$analysis = 'core, no DNA'
predict_res = predict(multinom_res, testSet_scaled)
conf_mat_results_minimal = caret::confusionMatrix(table(predict_res, testSet_scaled %>% .$cell_line))



multiclass_roc_relevant_task = multiclass_roc_full%>%
  bind_rows(multiclass_roc_core,multiclass_roc_minimal) %>%
  dplyr::filter(TASK == group1 | TASK == group2) 


p1 = ggplot(multiclass_roc_relevant_task, aes(x=specificity_min1, y = sensitivity, color = comparison ))+
  theme_bw() +
  theme(text = element_text(size = 6,family = 'sans')) +
  geom_abline(1,intercept = 0 , size = 0.5, color  ='darkgray',  linetype = 'dashed')+
  geom_line(aes(group =paste0(TASK, comparison)), linewidth = 0.5) +
  xlab('1-specificity') +
  ylab('sensitivity') +
  facet_wrap(~analysis)
p1

#############################################
## sensitivity, specificity, accuracy summary
#############################################
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

conf_mat_results_plot =conf_mat_results%>%
  dplyr::select(Sensitivity, Specificity, cell_line, analysis, `Balanced Accuracy`) %>%
  gather(key = metric, value =value, -cell_line,-analysis ) %>%
  mutate(metric =factor(metric, levels = c('Sensitivity','Specificity','Balanced Accuracy')),
         cell_line = factor(cell_line, levels =cell_line_levels))

p2 =  ggplot(conf_mat_results_plot, aes(x = analysis, y = value))+
  theme_bw() +
  theme(text = element_text(size = 6,family = 'sans')) +
  geom_boxplot(size = 0.25) +
  geom_point(aes(color =cell_line), size = 0.25) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_brewer(palette = 'Dark2') +
  facet_wrap(~metric) 
p2
###################################################################
##  anova analysis
###################################################################

my_anovas=lapply(features_cellcycle_full,function(feat){
  # glm_res = glm(formula =  as.formula(paste0(feat, '~0+cell_line')), data = glm_input)
  glm_res = lm(formula =  as.formula(paste0(feat, '~0+ cell_line')), data = testSet_scaled)
  anova_res = anova(glm_res)
  anova_res$feature = feat
  return(anova_res)
})

calc_var_explained = function(data, anova_table, feature){
  # Get the total sum of squares
  total_ss <- sum((data[[feature]] - mean(data[[feature]]))^2)
  # Extract the Sum of Squares (SS) for the categorical variable
  ss_category <- anova_table$`Sum Sq`[1]
  # Calculate the proportion of variance explained by the categorical variable
  variance_explained <- ss_category / total_ss
  return(variance_explained*100)
}

my_var_explained= lapply(my_anovas, function(anova_table){
  calc_var_explained(testSet_scaled, anova_table, unique(anova_table$feature)) %>% data.frame(var_expl=.)%>% mutate(feature = unique(anova_table$feature))
}) %>% bind_rows()
p3 = ggplot(my_var_explained, aes(x= var_expl, y = reorder(feature, var_expl))) +
  theme_bw()+
  theme(text = element_text(size = 6,family = 'sans')) +
  geom_col()
## plot results, expected to be a little different from original manuscript due to much smaller data sample size but general trends are the same
patchwork::wrap_plots(p1,p2,p3)



