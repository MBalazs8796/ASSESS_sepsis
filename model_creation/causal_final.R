library(readr)
library(grf)
library(tidyverse)
library(caret)
library(rpart)
library(rpart.plot)
library(doParallel)
library(foreach)
library(purrr)
library(yaml)

calc_ATE_k <- function(Y, W, X, k = 5, n_cores = parallel::detectCores() - 1){
  
  library(grf)
  library(caret)
  library(dplyr)
  library(doParallel)
  library(foreach)
  
  set.seed(42)
  
  folds <- createFolds(Y, k = k, list = TRUE, returnTrain = TRUE)
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  results <- foreach(i = 1:k,
                     .packages = c("grf","dplyr"),
                     .combine = "c") %dopar% {
                       
                       tryCatch({
                         
                         train_idx <- folds[[i]]
                         test_idx <- setdiff(seq_along(Y), train_idx)
                         
                         X_train <- X[train_idx,,drop=FALSE]
                         X_test  <- X[test_idx,,drop=FALSE]
                         
                         Y_train <- Y[train_idx]
                         Y_test  <- Y[test_idx]
                         
                         W_train <- W[train_idx]
                         W_test  <- W[test_idx]
                         
                         if(length(unique(W_test)) < 2){
                           warning(paste("Fold",i,"skipped: no treatment variation"))
                           return(list(list(
                             fold=i,ATE=NA,top_vars=NA,split_summary=NULL,
                             case=NULL,aggr=NULL
                           )))
                         }
                         
                         Y.forest <- regression_forest(X_train,Y_train)
                         Y.hat <- predict(Y.forest,X_test)$predictions
                         
                         W.forest <- regression_forest(X_train,W_train)
                         W.hat <- predict(W.forest,X_test)$predictions
                         
                         cf <- causal_forest(
                           X_test, Y_test, W_test,
                           Y.hat = Y.hat,
                           W.hat = W.hat,
                           num.trees = 2000,
                           tune.parameters = "all",
                           honesty = FALSE,
                           seed = 123
                         )
                         
                         ate_fold <- suppressWarnings(
                           average_treatment_effect(cf, target.sample="treated")[["estimate"]]
                         )
                         
                         imp <- variable_importance(cf)
                         var_names <- colnames(X_test)
                         
                         imp_df <- data.frame(var=var_names,imp=imp)
                         
                         top_vars <- imp_df %>%
                           arrange(desc(imp)) %>%
                           slice_head(n=10) %>%
                           pull(var) %>%
                           paste(collapse=", ")
                         
                         tau_hat <- predict(cf)$predictions
                         tau_hat[!is.finite(tau_hat)] <- NA
                         
                         case_df <- data.frame(tau=tau_hat,X_test)
                         
                         if(sum(!is.na(case_df$tau)) < 10){
                           warning(paste("Fold",i,"too few valid tau predictions"))
                           return(list(list(
                             fold=i,
                             ATE=ate_fold,
                             top_vars=top_vars,
                             split_summary=NULL,
                             case=NULL,
                             aggr=NULL
                           )))
                         }
                         
                         breaks <- quantile(case_df$tau,
                                            probs=seq(0,1,0.2),
                                            na.rm=TRUE)
                         
                         breaks <- unique(breaks)
                         
                         if(length(breaks) < 2){
                           case_df$group <- NA
                         } else{
                           case_df$group <- cut(case_df$tau,
                                                breaks=breaks,
                                                include.lowest=TRUE)
                         }
                         
                         agg_res <- aggregate(tau~group,case_df,mean)
                         
                         num_cols <- names(case_df)[sapply(case_df,is.numeric)]
                         
                         summary_df <- case_df %>%
                           group_by(group) %>%
                           summarise(
                             group_size=n(),
                             across(all_of(num_cols),~mean(.x,na.rm=TRUE)),
                             .groups="drop"
                           )
                         
                         split_summary <- tryCatch({
                           
                           trees <- min(2000,cf$num.trees)
                           
                           splits <- do.call(rbind,lapply(1:trees,function(tid){
                             
                             tr <- get_tree(cf,tid)
                             
                             df <- data.frame(var=character(),
                                              split_value=double())
                             
                             for(node in tr[[3]]){
                               if(!node$is_leaf){
                                 df <- rbind(df,data.frame(
                                   var=colnames(X)[node$split_variable],
                                   split_value=node$split_value
                                 ))
                               }
                             }
                             df
                           }))
                           
                           if(is.null(splits) || nrow(splits)==0) return(NULL)
                           
                           splits %>%
                             group_by(var) %>%
                             summarise(
                               mean_split=mean(split_value),
                               sd_split=sd(split_value),
                               n_splits=n(),
                               .groups="drop"
                             )
                           
                         },error=function(e) NULL)
                         
                         list(list(
                           fold=i,
                           ATE=ate_fold,
                           top_vars=top_vars,
                           split_summary=split_summary,
                           case=summary_df,
                           aggr=agg_res
                         ))
                         
                       }, error=function(e){
                         
                         warning(paste("Fold",i,"failed:",e$message))
                         
                         list(list(
                           fold=i,ATE=NA,top_vars=NA,
                           split_summary=NULL,case=NULL,aggr=NULL
                         ))
                         
                       })
                     }
  
  stopCluster(cl)
  
  ate_vals <- sapply(results, function(x) {
    val <- x$ATE
    if (is.null(val) || !is.numeric(val) || length(val) == 0) {
      return(NA_real_)
    }
    as.numeric(val[1])
  })
  
  mean_ate <- mean(ate_vals, na.rm = TRUE)
  
  sd_ate <- if(sum(!is.na(ate_vals)) > 1) {
    sd(ate_vals, na.rm = TRUE)
  } else {
    NA_real_
  }
  
  return(list(
    mean_ATE = mean_ate,
    sd_ATE = sd_ate,
    fold_results = results
  ))
}

calc_ATE_k_old <- function(Y, W, X, k = 5, n_cores = parallel::detectCores() - 1){   
  set.seed(42)
  folds <- createFolds(Y, k = k, list = TRUE, returnTrain = TRUE)
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  results <- foreach(i = 1:k, .combine = 'c', .packages = c("grf", "caret", "dplyr", "purrr")) %dopar% {
    tryCatch({
      train_idx <- folds[[i]]
      test_idx <- setdiff(1:length(Y), train_idx)
      
      X_train <- X[train_idx, , drop = FALSE]
      X_test  <- X[test_idx, , drop = FALSE]
      Y_train <- Y[train_idx]
      Y_test  <- Y[test_idx]
      W_train <- W[train_idx]
      W_test  <- W[test_idx]
      
      if (length(unique(W_test)) < 2) {
        warning(paste("Skipping fold", i, "- only one treatment level present"))
        return(data.frame(fold = i, ATE = NA, top_vars = NA, split_summary="2"))
      }
      
      Y.forest = regression_forest(X_train, Y_train)
      Y.hat = predict(Y.forest, X_test)$predictions
      W.forest = regression_forest(X_train, W_train)
      W.hat = predict(W.forest, X_test)$predictions
      
      cf <- tryCatch({
        causal_forest(
          X_test, Y_test, W_test,
          Y.hat = Y.hat, W.hat = W.hat,
          num.trees = 2000,
          tune.parameters = "all",
          honesty = FALSE,
          seed = 123
        )
      }, error = function(e) {
        warning(paste("Fold", i, "failed to fit:", e$message))
        return(NULL)
      })
      
      if (is.null(cf)) {
        return(data.frame(fold = i, ATE = NA, top_vars = NA))
      }
      
      tree_splits <- do.call(rbind, lapply(1:2000, function(tree_id) {
        t <- get_tree(cf, tree_id)
        df <- data.frame(var = character(), split_value = double())
        for (variable in t[3]) {
          for (b in variable){
            if(b$is_leaf == FALSE){
              df <- rbind(df, data.frame(var = names(X)[b$split_variable], split_value = b$split_value))
            }
          }
          return (df)
        }}))
      if (is.null(tree_splits) || nrow(tree_splits) == 0) {
        split_summary_fold <- data.frame(var = NA, mean_split = NA, sd_split = NA, n_splits = 0)
      } else {
        split_summary_fold <- tree_splits %>%
          group_by(var) %>%
          summarise(
            mean_split = mean(split_value, na.rm = TRUE),
            sd_split = sd(split_value, na.rm = TRUE),
            n_splits = n(),
            .groups = "drop"
          )
      }
      ate_fold <- tryCatch({
        average_treatment_effect(cf, target.sample = "treated")[["estimate"]]
      }, error = function(e) {
        warning(paste("ATE failed in fold", i, ":", e$message))
        return(NA)
      })
      
      varimp <- tryCatch({
        imp <- variable_importance(cf)
        var_names <- colnames(X_test)
        imp_df <- data.frame(var = var_names, imp = imp)
        top_names <- imp_df %>%
          arrange(desc(imp)) %>%
          slice_head(n = 10) %>%
          pull(var)
        paste(top_names, collapse = ", ")
      }, error = function(e) {
        warning(paste("Var importance failed in fold", i))
        return(NA)
      })
      tau_hat <- predict(cf)$predictions
      case_df <- data.frame(tau = tau_hat, X_test)
      case_df_ordered <- case_df[order(case_df$tau), ]
      tau_hat <- predict(cf)$predictions
      
      tau_hat[!is.finite(tau_hat)] <- NA
      
      case_df <- data.frame(tau = tau_hat, X_test)
      
      if (sum(!is.na(case_df$tau)) < 10) {
        warning(paste("Fold", i, "too many NA in tau"))
        return(list(
          fold = i,
          ATE = ate_fold,
          top_vars = varimp,
          split_summary = split_summary_fold,
          case = NA,
          aggr = NA
        ))
      }
      
      breaks <- quantile(case_df$tau, probs = seq(0,1,0.2), na.rm = TRUE)
      
      breaks <- unique(breaks)
      
      if (length(breaks) < 2) {
        warning(paste("Fold", i, "not enough variation in tau"))
        case_df$group <- NA
      } else {
        case_df$group <- cut(case_df$tau, breaks = breaks, include.lowest = TRUE)
      }
      agg_res <- aggregate(tau ~ group, case_df, mean)

      cols <- colnames(case_df %>% 
                         select(where(is.numeric)) %>% 
                         select(where(~ any(!is.na(.x)))))
      summary_df <- case_df %>%
        group_by(group) %>%
        summarise(
          group_size = n(),
          across(all_of(cols), ~ mean(.x, na.rm = TRUE))
          )
      
      list(fold = i, ATE = ate_fold, top_vars = varimp,  split_summary = split_summary_fold, case = summary_df, aggr = agg_res, stringsAsFactors = FALSE)
    }, error = function(e) {
      warning(paste("Fold", i, "crashed:", e$message))
      return(list(fold=i, ATE=NA, top_vars=NA, split_summary=NULL, case=NULL, aggr=NULL))
    })
    }
  
  stopCluster(cl)
  
  mean_ate <- mean(results$ATE, na.rm = TRUE)
  sd_ate <- sd(results$ATE, na.rm = TRUE)
  
  return(list(
    mean_ATE = mean_ate,
    sd_ATE = sd_ate,
    fold_results = results
  ))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    success <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!success) {
      stop(sprintf("Failed to create directory: %s", path))
    }
  }
  invisible(path)
}

calc_exclusive_ATE_k <- function(selection_criteria, exl_selectors, first_input, negative_df){
  inner_input <- first_input %>%
    filter(if_all(all_of(exl_selectors), ~ . == TRUE),    
           if_all(all_of(setdiff(selection_criteria, exl_selectors)), ~ . == FALSE))  
  inner_input <- bind_rows(negative_df, inner_input)
  Y <- inner_input$vasopressor_free_days 
  W <- apply(inner_input[,exl_selectors], 1, all)
  return(
    calc_ATE_k(Y, W, inner_input[, !colnames(inner_input) %in% selection_criteria] %>%
                 select(-vasopressor_free_days), k = 5) 
  )
}

calc_ATE <- function(Y, W, X){
  Y.forest = regression_forest(X, Y)
  Y.hat = predict(Y.forest)$predictions
  W.forest = regression_forest(X, W)
  W.hat = predict(W.forest)$predictions
  cf.raw = causal_forest(X, Y, W,
                         Y.hat = Y.hat, W.hat=W.hat, num.trees = 4000,
                         tune.parameters = "all")
  varimp = variable_importance(cf.raw)
  selected.idx = which(varimp >= mean (varimp))
  cf = causal_forest( X[ , selected.idx ], Y, W,
                      Y.hat = Y.hat, W.hat = W.hat, num.trees = 4000,
                      tune.parameters = "all")
  tau_hat <- predict(cf)$predictions
  case_df <- data.frame(tau = tau_hat, X)
  case_df_ordered <- case_df[order(case_df$tau), ]
  case_df$group <- cut(case_df$tau, breaks = quantile(case_df$tau, probs = seq(0,1,0.2)), include.lowest = TRUE)
  agg_res <- aggregate(tau ~ group, case_df, mean)
  
  cols <- colnames(case_df %>% 
                     select(where(is.numeric)) %>% 
                     select(where(~ any(!is.na(.x)))))
  summary_df <- case_df %>%
    group_by(group) %>%
    summarise(
      group_size = n(),
      across(all_of(cols), ~ mean(.x, na.rm = TRUE))
    )
  return(list(average_treatment_effect( cf ), X[ , selected.idx ], cf, summary_df))
}

importance_ranker <- function(Y, W, X){
  Y.forest = regression_forest(X, Y)
  Y.hat = predict(Y.forest)$predictions
  W.forest = regression_forest(X, W)
  W.hat = predict(W.forest)$predictions
  cf.raw = causal_forest(X, Y, W,
                         Y.hat = Y.hat, W.hat=W.hat, num.trees = 4000)
  varimp = variable_importance(cf.raw)
  n_df <- data.frame(
    column_numbers = colnames(X),
    importance = varimp 
  )
  return(n_df[order(-n_df$importance),])
}


calc_exl_importance <- function(selection_criteria, exl_selectors, first_input, negative_df){
  inner_input <- first_input %>%
    filter(if_all(all_of(exl_selectors), ~ . == TRUE),  
           if_all(all_of(setdiff(selection_criteria, exl_selectors)), ~ . == FALSE)) 
  inner_input <- bind_rows(negative_df, inner_input)
  Y <- inner_input$vasopressor_free_days
  W <- apply(inner_input[,exl_selectors], 1, all)
  return(importance_ranker(Y, W, inner_input[, !colnames(inner_input) %in% selection_criteria] %>% select(-vasopressor_free_days)))
}


calc_exclusive_ATE <- function(selection_criteria, exl_selectors, first_input, negative_df){
  inner_input <- first_input %>%
    filter(if_all(all_of(exl_selectors), ~ . == TRUE),    
           if_all(all_of(setdiff(selection_criteria, exl_selectors)), ~ . == FALSE)) 
  inner_input <- bind_rows(negative_df, inner_input)
  Y <- inner_input$vasopressor_free_days
  W <- apply(inner_input[,exl_selectors], 1, all)
  return(calc_ATE(Y, W, inner_input[, !colnames(inner_input) %in% selection_criteria] %>% select(-vasopressor_free_days)))
}

is_boolean <- function(x){
  return (typeof(x) == "logical")
}

selection_criteria <- c(
  "vasopressor_prod",
  "fluid_admin_prod",
  "lactate_prod",
  "antibiotics_prod",
  "blood_check_prod"
)

config <- read_yaml("config.yaml")

first_input <- read_csv(config$final_dataset_path)

first_input <- first_input %>% select(-hadm_id, -subject_id, -apsiii, -baseline, -early_death_flag, -time_zero, -group_label)


first_input$GCS <- first_input$`GCS - Eye Opening` + first_input$`GCS - Verbal Response` + first_input$`GCS - Motor Response`

first_input <- first_input %>% select(-`GCS - Eye Opening`, -`GCS - Verbal Response`, -`GCS - Motor Response`)

first_input$vasopressor_prod <- first_input$vasopressor_init_prod & first_input$vasopressor_recom_prod

first_input <- first_input %>% select(-vasopressor_init_prod, -vasopressor_recom_prod, -vasopressor_recom )

first_input



c_target <- "vasopressor_free_days"
first_input <- first_input %>% select(-icu_free_days, -icu_los,                                          
                                      -hosp_free_days, -hospital_expire_flag, -icu_expire_flag, -hosp_los,
                                      -one_year_mort, -one_month_mort, -readmission_during_same_stay, -readmission_within_90)

first_input <- first_input %>% mutate_if(is_boolean, as.numeric)

sapply(first_input, typeof)

for(i in colnames(first_input)){
  if (i %in% selection_criteria){
    first_input[,i] = lapply(first_input[,i], as.logical)
  }
}


n_df = first_input %>%
  filter(if_all(all_of(selection_criteria), ~ . == FALSE))

exl_combinations <- list(
  c("antibiotics_prod"),
  c("fluid_admin_prod", "lactate_prod", "antibiotics_prod", "vasopressor_prod"),
  c("fluid_admin_prod", "lactate_prod", "antibiotics_prod", "blood_check_prod"),
  c("fluid_admin_prod", "lactate_prod", "antibiotics_prod"),
  c("lactate_prod", "antibiotics_prod"),
  c("fluid_admin_prod", "lactate_prod"),
  c("fluid_admin_prod", "antibiotics_prod"),
  c("fluid_admin_prod"),
  c("lactate_prod")
  
)

first_input %>%
  select(where(~ !is.numeric(.))) %>%
  names()

run_pipeline_old <- function(target_df, name) {
message(name)
for (exl_selectors in exl_combinations) {
  message("Running for: ", paste(exl_selectors, collapse = ", "))
 
  c_list <- calc_exclusive_ATE(selection_criteria, exl_selectors, target_df, n_df)
  importance <- calc_exl_importance(selection_criteria, exl_selectors, target_df, n_df)

  filename_stem <- paste0(substr(exl_selectors, 1, 1), collapse = "")
  base_path <- paste0(config$importance_path, name, '/', c_target, "/")
  ensure_dir(base_path)

  write.csv(
    c_list[[4]],
    file = paste0(base_path, filename_stem, "_groups",".csv"),
    row.names = FALSE
  )
  
  
  write.csv(
    importance,
    file = paste0(base_path, filename_stem, ".csv"),
    row.names = FALSE
  )
  

  res <- data.frame(
    estimate = c_list[[1]][["estimate"]],
    err = c_list[[1]][["std.err"]],
    ci = c_list[[1]][["std.err"]] * 1.96,
    imp_var = paste(colnames(c_list[[2]]), collapse = ",")
  )
  
  write.csv(
    res,
    file = paste0(base_path, filename_stem, "_estimate_res.csv"),
    row.names = FALSE
  )
  
  c_list <- calc_exclusive_ATE_k(selection_criteria, exl_selectors, target_df, n_df)
  
  res_k <- data.frame(
    mean_ATE = c_list$mean_ATE,                 
    sd_ATE = c_list$sd_ATE,                    
    imp_var = NA,                              
    fold_estimates = paste(c_list$k_fold_results, collapse = ", ")
  )
  
  split_df <-  data.frame(
    vars = c_list$fold_results$split_summary.var,
    mean_value = c_list$fold_results$split_summary.mean_split,
    sd_value = c_list$fold_results$split_summary.sd_split,
    density = c_list$fold_results$split_summary.n_splits
  )
  for (i in 1:5) {
    write.csv(
      c_list$fold_results[ i, "case"][[1]],
      file = paste0(base_path, filename_stem, "_fold_case", i, ".csv"),
      row.names = FALSE
    )
    write.csv(
      c_list$fold_results[ i, "aggr"][[1]],
      file = paste0(base_path, filename_stem, "_fold_aggr", i, ".csv"),
      row.names = FALSE
    )
    
    write.csv(
      c_list$fold_results[ i, "ATE"],
      file = paste0(base_path, filename_stem, "_fold_ATE", i, ".csv"),
      row.names = FALSE
    )
    
    write.csv(
      c_list$fold_results[ i, "top_vars"],
      file = paste0(base_path, filename_stem, "_fold_top_vars", i, ".csv"),
      row.names = FALSE
    )
  }
}}

run_pipeline <- function(target_df, name) {
  message(name)
  for (exl_selectors in exl_combinations) {
    message("Running for: ", paste(exl_selectors, collapse = ", "))

    c_list <- calc_exclusive_ATE_k(selection_criteria, exl_selectors, target_df, n_df)
    filename_stem <- paste0(substr(exl_selectors, 1, 1), collapse = "")
    base_path <- paste0(config$importance_path, name, '/', c_target, "/")
    ensure_dir(base_path)
    
    fold_ates <- sapply(c_list$fold_results, function(x) x$ATE)
    
    res_k <- data.frame(
      mean_ATE = c_list$mean_ATE,
      sd_ATE = c_list$sd_ATE,
      imp_var = NA,
      fold_estimates = paste(fold_ates, collapse = ", ")
    )
    
    write.csv(
      res_k,
      file = paste0(base_path, filename_stem, "_kfold_estimate_res.csv"),
      row.names = FALSE
    )
    
    
    split_list <- lapply(c_list$fold_results, function(x) x$split_summary)
    split_df <- do.call(rbind, split_list)
    
    if(!is.null(split_df) && nrow(split_df) > 0){
      
      split_df_out <- data.frame(
        vars = split_df$var,
        mean_value = split_df$mean_split,
        sd_value = split_df$sd_split,
        density = split_df$n_splits
      )
      
      write.csv(
        split_df_out,
        file = paste0(base_path, filename_stem, "_split_summary.csv"),
        row.names = FALSE
      )
    }
    
    
    for(i in seq_along(c_list$fold_results)){
      
      fold <- c_list$fold_results[[i]]
      
      if(!is.null(fold$case)){
        write.csv(
          fold$case,
          file = paste0(base_path, filename_stem, "_fold_case", i, ".csv"),
          row.names = FALSE
        )
      }
      
      if(!is.null(fold$aggr)){
        write.csv(
          fold$aggr,
          file = paste0(base_path, filename_stem, "_fold_aggr", i, ".csv"),
          row.names = FALSE
        )
      }
      
      write.csv(
        data.frame(ATE = fold$ATE),
        file = paste0(base_path, filename_stem, "_fold_ATE", i, ".csv"),
        row.names = FALSE
      )
      
      write.csv(
        data.frame(top_vars = fold$top_vars),
        file = paste0(base_path, filename_stem, "_fold_top_vars", i, ".csv"),
        row.names = FALSE
      )
      
    }
  }}

lactate_df = first_input %>% select(-chf_history)

run_pipeline(
  lactate_df[!is.na(lactate_df$first_lactate_mes) & lactate_df$first_lactate_mes >= 4, ] %>% select(-first_lactate_mes),
  "high_lactate")
run_pipeline(
  lactate_df[!is.na(lactate_df$first_lactate_mes) & lactate_df$first_lactate_mes < 4, ] %>% select(-first_lactate_mes),
  "low_lactate")

first_input = first_input %>% select(-first_lactate_mes)

run_pipeline(
  first_input[first_input$chf_history == TRUE, ] %>% select(-chf_history),
  "chf_history")
run_pipeline(
  first_input[first_input$chf_history == FALSE, ] %>% select(-chf_history),
  "no_chf_history")

first_input = first_input %>% select(-chf_history)
med_sofa <- median(first_input$sofa_24hours, na.rm = TRUE)


run_pipeline(first_input[!is.na(first_input$sofa_24hours) & first_input$sofa_24hours >= med_sofa, ], "high_sofa")
run_pipeline(first_input[!is.na(first_input$sofa_24hours) & first_input$sofa_24hours < med_sofa, ], "low_sofa")

message("Done!")
