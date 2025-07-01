## We provide example code to compute an unsupervised aberrancy score. 

## RANN::nn2 has weird behavior with the query search with the indexes of the query dataset (full) and landmark dataset.
## to fix this, the index of the query (fill) and landmark datasets must match the first rows. All remaining indexes won't matter
wrangle_data_for_knn_index_resolution = function(full_data,LM_col=NULL,LM_group=NULL, LM_data=NULL){
  if(!'cell.id' %in% colnames(full_data) ){
    stop('input data must have a cell.id column')
  }
  if( !is.null(LM_data)){
    if(!'cell.id' %in% colnames(LM_data) ){
      stop('landmark data must have a cell.id column')
    }
    QUERY_DF = LM_data %>%
      bind_rows(full_data %>% dplyr::filter(!cell.id %in%  LM_data$cell.id))
  } else {
    
    if(!LM_col %in% colnames(full_data) ){
      stop(paste0('input data does not contain the provided column for the "LM_col" arguement: ',LM_col))
    }
    
    if(!LM_group %in% unique(full_data[[LM_col]]) ){
      stop(paste0(LM_col, ' does not contain ',LM_group))
    }
    # LM_col = 'condition'
    # LM_group = 'WT'
    LM_df = full_data %>% .[.[[LM_col]]==LM_group,]
    PROJECT_DF = full_data %>% .[!.[[LM_col]]==LM_group,]
    QUERY_DF = bind_rows(LM_df,PROJECT_DF)
  }
  
  return(QUERY_DF)
}
## The threshold for aberrancy is up to the user.
##' full_data: single-cell dataframe
##' LM_col: column name (string) containing your identifier for your group of interest (eg, "treatment_group") as the landmarks. 
##' LM_group: column name (string) containing your landmark group (eg, "NoTx") 
##' features: vector of features to use in your analysis c('feature1','feature2','feature3'). The columns must be in the full_data  dataframe that you provide. 
##' k_ave: number of neighbors to use
##' return_knn_graph: Default is FALSE. Returns the knn graph.
##' return_nearest_cellid: Default is FALSE. Returns the nearest cell id.


compute_knn_aberrancy= function(full_data, LM_col=NULL, LM_group, features, k_ave =2, return_knn_graph=FALSE, return_nearest_cellid=FALSE){
  if (is.null(LM_col)){
    message('performing search  on full dataset...')
    df_landmark_features_only=full_data%>% dplyr::select(all_of(features))
    full_data_features_only=full_data%>% dplyr::select(all_of(features))
    restructured_data=full_data
  } else {
    message('restructuring data for knn search...')
    restructured_data = wrangle_data_for_knn_index_resolution(full_data, LM_col, LM_group)
    df_landmark_features_only = restructured_data %>% .[.[[LM_col]]==LM_group,] %>% dplyr::select(all_of(features))
    full_data_features_only = restructured_data %>% dplyr::select(all_of(features))
  }
  
  
  message('performing knn search...')
  if (is.null(LM_col)){
    knn_res_full = RANN::nn2(data = df_landmark_features_only, searchtype = 'standard',eps=0,treetyp = 'kd',
                             # radius = 100, ## for exact search
                             k = k_ave+1)
  } else {
    knn_res_full = RANN::nn2(data = df_landmark_features_only, query = full_data_features_only,searchtype = 'standard',eps=0,treetyp = 'kd',
                             # radius = 100, ## for exact search
                             k = k_ave+1)
  }
  
  if(return_knn_graph ==T){
    return(knn_res_full)
  }
  
  if (return_nearest_cellid==T){
    k_ave =1
    score_res_project=knn_res_full$nn.idx[,1:k_ave] %>% data.frame(check.rows = F, check.names = T)
    score_res_project$cell.id = restructured_data$cell.id
    df_landmark_features_only$cell.id = restructured_data %>% .[.[[LM_col]]==LM_group,] %>%.$cell.id
    return( list(df_landmark_features_only=df_landmark_features_only, score_res_project=score_res_project) )
  }
  
  message('extracting closest cell and returning scores with original indexing...')
  if (k_ave==1  ){ 
    message('returning distance of nearest landmark neighbor')
    knn_res_full$nn.dists[,1] = ifelse(knn_res_full$nn.idx[,1] == 1:nrow(knn_res_full$nn.idx), knn_res_full$nn.dists[,2], knn_res_full$nn.dists[,1])
    score_res_project= knn_res_full$nn.dists[,1] %>% unlist()
  } else if (k_ave > 1 ){
    knn_res_full$nn.dists[,1] = ifelse(knn_res_full$nn.idx[,1] == 1:nrow(knn_res_full$nn.idx), knn_res_full$nn.dists[,(k_ave+1)], knn_res_full$nn.dists[,1])
    # return(knn_res_full$nn.dists)
    score_res_project= rowMeans(knn_res_full$nn.dists[,1:k_ave])
  } else {
    score_res_project = knn_res_full$nn.dists[,1] %>% unlist()  #%>% data.frame(dist_from_nearest_cell = .) 
  }
  names(score_res_project) = restructured_data$cell.id
  # table(score_res_project$cell.id[match(full_data$cell.id, score_res_project$cell.id)] ==full_data$cell.id)
  score_res_project=score_res_project[match(full_data$cell.id, names(score_res_project))]
  return(score_res_project)
}
