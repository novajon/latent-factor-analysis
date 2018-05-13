## Descriptive Analytics Project: Josefine, Dimitri, Jonathan
# install.packages("VIM")
library("VIM")
rm(list=ls())

## @knitr part1
##########################################################
#                   Model Initialisation                 #
##########################################################
###TODO: Include checks: data must be of class data.frame, 
### ncol of mapping table must be 2, mapping table must be of kind matrix
my.pls.model <- function(data,inner_model_mappings,outer_model_mappings){
  
  if(!is.data.frame(outer_model_mappings) || !is.data.frame(inner_model_mappings) || !all(names(outer_model_mappings)==c("Source","Target")) || !all(names(inner_model_mappings)==c("Source", "Target"))){
    stop("Inner or Outer Model Mappings is not a data frame with source and target")
  }
  
  if(!is.matrix(data)){
    stop("Data is not a matrix")
  }
  
  # reflective:latent variables (LV) --> manifest variables (MV)
  # formative: manifest variables (MV) --> latent variables (LV)
  LV <- unique(unlist(c(inner_model_mappings)))
  MV <- setdiff(unique(unlist(c(t(outer_model_mappings)))),LV)
  
  if (!all(MV %in% colnames(data))){
    stop("Data does not contain all manifest variables!")
  }
    
  #### inner model
  inner_model <- matrix(c(0), nrow=length(LV), ncol=length(LV), dimnames = list(LV, LV))
  inner_model <- fillModel (inner_model, inner_model_mappings, LV)
  #### outer model
  outer_model <- matrix(c(0), nrow=length(MV), ncol=length(LV), dimnames = list(MV,LV))
  blocks <- getBlocks(outer_model_mappings, inner_model)
  outer_weights <- fillModel(outer_model, outer_model_mappings, LV)
  
  result <- list()
  result$MV <- MV
  result$LV <- LV
  result$inner_model <- inner_model
  result$blocks <- blocks
  result$outer_model <- outer_weights
  return(result)
}

fillModel <- function (model, connections, LV) {
  apply(connections,1,function(x){
      if (x["Source"] %in% LV){
        model[x["Target"], x["Source"]]<<-1
      } else {
        model[x["Source"], x["Target"]]<<-1
      }
    })
  return (model)
}

getBlocks<- function (outer_model_mappings, inner_model){
    blocks = list()
    LV <- colnames(inner_model)
    for (lv in LV) {
      blocks[[lv]]$predecessors <- colnames(inner_model)[inner_model[lv,]==1]
      blocks[[lv]]$successors <- rownames(inner_model)[inner_model[,lv]==1]
      if (any(outer_model_mappings[,"Source"]==lv) && any(outer_model_mappings[,"Target"]==lv)) {
        # catch case where both reflective and formative
        stop("Latent variable cannot be reflective and formative")
      }
      if (any(outer_model_mappings[,"Source"]==lv)) {
        blocks[[lv]]$MV <- outer_model_mappings[,"Target"][(which (outer_model_mappings[,"Source"]==lv))]   
        blocks[[lv]]$mode <- "A"
      }
      if (any(outer_model_mappings[,"Target"]==lv)){ 
        blocks[[lv]]$MV <- outer_model_mappings[,"Source"][(which (outer_model_mappings[,"Target"]==lv))]   
        blocks[[lv]]$mode <- "B"
      }
    }
    return (blocks)
}

##########################################################
#                     Imputation                         #
##########################################################

imputeMean <- function (data, dimension=1) {
  for (i in 1:dim(data)[dimension]) {
    # data_arr either row or column (dimension 1 or 2)
    if (dimension==1) data_arr <- data [i, ]
    if (dimension==2) data_arr <- data [, i]
    missing_index <- which(!is.na(data_arr))
    arr_mean <- round(sum(data_arr[missing_index])/(dim(data)[3-dimension]-sum(is.na(data_arr))), 0)
    if (dimension==1) data [i,] [which(is.na(data_arr))] <- arr_mean
    if (dimension==2) data [,i] [which(is.na(data_arr))] <- arr_mean
  }
  return (data)
}

imputeSpecialMean <- function (data) {
  data_imputed <- as.matrix(data)
  total_avg_rating = mean(data_imputed, na.rm=TRUE)
  for (j in 1:dim(data)[1]){
    for (i in 1:dim(data)[2]){
      if (is.na(data[j,i])){
        data_imputed[j,i] <- round(mean(data[,i], na.rm=TRUE) + (mean(data[j,], na.rm=TRUE) - total_avg_rating))
      }
    }
  }
}

##### k-nearest neighbour
imputeNearestNeighbour <- function (data, neighboursize = 10) {
  return(kNN(data, imp_var=F, k=neighboursize))
}
##### delete rows with missing values
imputeDeletion <- function (data) {
  return (na.omit(data))  
}

##########################################################
#                   Data Treatment                       #
##########################################################
dataTreatment <- function(data, impute_method=imputeDeletion){
  
  data_imputed <- impute_method(data)
  # data_imputed <- data_imputed[!apply(data_imputed,1,var)<=0.05,]
  data_normalised <- scale(data_imputed)
  
  return(data_normalised)
}

##########################################################
#                Inner Model Approximation               #
##########################################################

centroidWeighting <- function (inner_model, blocks, Y_approx) {
  # blocks not needed
  correlation <- cor(Y_approx)
  symmetric_model <- inner_model + t(inner_model)
  return(sign(correlation) *  symmetric_model)
}

factorialWeighting <- function (inner_model, blocks, Y_approx) {
  # blocks not needed
  correlation <- cor(Y_approx)
  symmetric_model <- inner_model + t(inner_model)
  return(correlation *  symmetric_model)
}

pathWeighting <- function (inner_model, blocks, Y_approx) {
  correlation <- cor(Y_approx)
  E_approx <- correlation * inner_model
  for(lv in colnames(inner_model)) {
    if (length(blocks[[lv]]$predecessors)>0) {
      linear_model <- lm(data = as.data.frame(Y_approx), formula = as.formula(paste(lv, " ~ ", paste(blocks[[lv]]$predecessors, collapse = " + "))))
      for (regressor in 2:length(linear_model$coefficients)) {
        E_approx[names(linear_model$coefficients)[regressor], lv] <- linear_model$coefficients[regressor]
      }
    }
  }
  return (E_approx)
}

##########################################################
#                Model Calculation Functions             #
##########################################################

calculateCoefficients <- function (y_value, inner_model) {
  # inner model weights
  pathCoeff <- matrix(c(0), nrow = length(colnames(inner_model)), ncol=length(colnames(inner_model)), dimnames = list(colnames(inner_model), colnames(inner_model)))
  for (regressand in colnames(inner_model)) {
    regressors <- colnames(inner_model)[inner_model[regressand, ]==1]
    if (length(regressors)>0) {
      lv_string <- paste(regressors, collapse="+")
      regression <- lm(data=as.data.frame(y_value), formula = as.formula(paste(regressand, " ~ ", lv_string)))
      sapply(regressors, function(x) {
        pathCoeff[x, regressand] <<- regression$coefficients[x]
      })
    }
  }
  return (pathCoeff)
}

calculateTotalEffects <- function(path_coefficients) {
  total_effects <- path_coefficients
  for (i in 2:length(path_coefficients)){
    matrix_factor <- path_coefficients
    for (j in 2:i){ matrix_factor <- matrix_factor %*% path_coefficients}
    total_effects <- total_effects + matrix_factor
  }
  return(total_effects)
}

##########################################################
#           Partial Least Squares Algorithm              #
##########################################################

my.pls <- function(data_normalised, model, threshold = 0.00001, inner_approximation_method=centroidWeighting){
  # 1) Initiate outer weights (set to 1)
  outer_model <- model$outer_model
  inner_model <- model$inner_model
  outer_model_weights <- outer_model

  if(!is.matrix(outer_model) || !is.matrix(outer_model)){
    stop("Inner or Outer Model is not a matrix")
  }
  
  result_model <- list()
  continue = T
  # i <- 1
  while (continue) {
    outer_model_weights_old <- outer_model_weights
    # print (i)
    # i <- i+1
    # 2) Outer Approximation: Approxiate Y (LV's)
    Y_approx_denormal <- data_normalised %*% outer_model_weights
    # normalize
    Y_approx <- scale(Y_approx_denormal)
    # 3) Approximate  inner model connections
    E_approx <- inner_approximation_method(inner_model, model$blocks, Y_approx)
    
    # 4) Inner approximation: Approxiate Z
    Z_approx_denormal <- Y_approx %*% E_approx
    # normalize
    Z_approx <- scale(Z_approx_denormal)
    # 5) Outer weight recalculation (w)
    # outer_model_weights <- cov(data_normalised, Z_approx) * outer_model
    ## divide into blocks, calculate & combine again
    for (lv in model$LV) {
      ### divide
      model$blocks[[lv]][["Z_approx"]] <- Z_approx[,lv]
      model$blocks[[lv]][["data"]] <- data_normalised[,as.vector(model$blocks[[lv]]$MV)]
      ### calculate
      if (model$blocks[[lv]]$mode=="A") {
        #### MODE A: COR(X, Y)
        model$blocks[[lv]][["outer_weights"]] <- cor(model$blocks[[lv]]$data, model$blocks[[lv]]$Z_approx)
      } else if (model$blocks[[lv]]$mode=="B") {
        #### MODE B: VAR(X)^-1 * COR(Y,X)
        var_inv <- solve(var (model$blocks[[lv]]$data))
        corr <- cor(model$blocks[[lv]]$Z_approx, model$blocks[[lv]]$data)
        model$blocks[[lv]][["outer_weights"]] <- var_inv %*% t(corr)
      }
      ### combine blocks
      outer_model_weights[rownames(model$blocks[[lv]]$outer_weights), lv] <- model$blocks[[lv]]$outer_weights
    }
    
    lv_scores <- data_normalised %*% outer_model_weights
    lv_scaled <- scale(lv_scores)
    sd <- attr(lv_scaled, "scaled:scale")
    outer_model_weights <- t(t(outer_model_weights) * (1/sd))
    
    # 6) Check for convergence
    continue <- !all(abs(outer_model_weights_old - outer_model_weights)<threshold)
  }
  #Final Y Matrix
  Y_final <- scale(data_normalised %*% outer_model_weights)
  
  # Estimation of Path coefficients & Total Effects
  path_coefficients <- calculateCoefficients(Y_final, model$inner_model)
  total_effects <- calculateTotalEffects(path_coefficients)
  
  # Computation of loadings
  crossloadings <- cor(data_normalised,Y_final);
  loadings <- crossloadings * outer_model
  crossloadings[crossloadings<0.6] <- 0
  
  result_model$weights <- outer_model_weights
  result_model$path_coefficients <- path_coefficients
  result_model$total_effects <- total_effects
  result_model$crossloadings <- crossloadings
  result_model$loadings <- loadings
  result_model$scores <- Y_approx
  result_model$blocks <- model$blocks
  return(result_model)
  
}

##########################################################
#                   Assessment Measures                  #
##########################################################
assessmentMeasures <- function (data, initial_model, result) {
  measures <- list()
  ## general tests
  measures$structural_model_tests$r_square <- rSquared(initial_model, result)
  measures$structural_model_tests$redundancy <- redundancy(initial_model, result)
  measures$structural_model_tests$goodness_of_fit <- gof(initial_model, result)
  ## reflective tests
  measures$reflective <- assessOuterReflective (data, initial_model, result)
  ## formative tests
  measures$formative <- assessOuterFormative(data, initial_model)
  return (measures)
}

assessOuterReflective <- function(data, initial_model, result){
  ### Examining covariance matrices by blocks
  measures <- list()
  measures$eigenvalue_analysis <- eigenvalueAnalysis(data, initial_model)
  measures$cronbachs_alpha <- cronbachsAlpha(data, initial_model)
  measures$cross_loadings_validation <- crossLoadingsValidation(result, initial_model)
  measures$dillon_goldstein <- dillonGoldstein(initial_model, result)
  measures$communality <- communality(initial_model, result)
  measures$insignificant_manifest_variables_based_on_communality <- outerLoadingsRelevanceTest(initial_model, result)
  return(measures)
}

assessOuterFormative <- function(data, initial_model){
  measures <- list()
  measures$collin_measure <- multiCollinearityPCAAnalysis(data, initial_model)
  return(measures)
}


########

# Cronbach Alpha
cronbachsAlpha <- function(data, initial_model){
  alphas <-list()
  for(lv in initial_model$LV){
    if(initial_model$blocks[[lv]]$mode == "A"){
      K = length(initial_model$blocks[[lv]]$MV)
      total = apply(data[,as.vector(initial_model$blocks[[lv]]$MV)],1,sum)
      alphas[[lv]] = (K/(K-1)) * (1 - sum(apply(data[,as.vector(initial_model$blocks[[lv]]$MV)], 2, var))/var(total))
    }
  }
  return(alphas)
}

# Dillon Goldstein Rho
dillonGoldstein <-function(initial_model, result){
  loadings = result$loadings
  rhos <- list()
  for(lv in initial_model$LV){
    rhos[[lv]] = sum(loadings[as.vector(initial_model$blocks[[lv]]$MV),lv])^2 / (sum(loadings[as.vector(initial_model$blocks[[lv]]$MV),lv])^2 + sum(1 - loadings[as.vector(initial_model$blocks[[lv]]$MV),lv]^2))
  }
  return(rhos)
}

# Communality Indexes
#TODO check block length > 1 --> comm indexes only make sense if block size > 1
communality <- function(initial_model, result){ 
  output <- list()
  output$communalities = result$loadings^2
  output$average_variance_extracted = matrix(nrow = 2, ncol=length(initial_model$LV), dimnames = list(c("# Observations", "Average Variance Extracted"), as.vector(initial_model$LV)))
  for(lv in initial_model$LV){
    output$average_variance_extracted[2,lv] = mean(output$communalities[as.vector(initial_model$blocks[[lv]]$MV),lv])
    output$average_variance_extracted[1,lv] = length(initial_model$blocks[[lv]]$MV)
  }
  
  num_indic = 0
  for(lv in initial_model$LV){
    num_indic = num_indic + length(initial_model$blocks[[lv]]$MV)
  }
  output$mean_communality = sum(output$communalities)/num_indic
  
  return(output)
  
}

# R Squared Indexes
rSquared <- function(initial_model, result){
  output <- list()
  output$r_squared <- list()
  n = 0
  tot=0
  
  for (lv in rownames(initial_model$inner_model)) {
    regressors <- initial_model$blocks[[lv]]$predecessors
    if(length(regressors) > 0){
      n = n + 1
      vector = result$scores[,lv]
      data = result$scores[,regressors]
      predicted = data %*% as.matrix(result$path_coefficients[regressors,lv])
      m = mean(result$scores[,lv])
      ss_tot = sum((vector - m)^2)
      ss_reg = sum((predicted - m)^2)
      ss_err = sum((vector - predicted)^2)
      output$r_squared[[lv]] <- 1 - ss_err/ss_tot
      tot = tot + output$r_squared[[lv]]
    }
  }
  
  output$mean_r_squared = tot/n
  return(output)
}

# Redundancy Indexes
redundancy <- function(initial_model, result){
  output <- list()
  output$redundancy <- list()
  tot_redundancy = 0
  n = 0
  
  for(lv in initial_model$LV){
    if(length(initial_model$blocks[[lv]]$predecessors) > 0){
      n = n + 1
      output$redundancy[[lv]] = rSquared(initial_model, result)$r_squared[[lv]] * communality(initial_model, result)$average_variance_extracted[2,lv]
      tot_redundancy = tot_redundancy + output$redundancy[[lv]] 
    }
  }
  
  if(n>0){
    output$tot_redundancy = tot_redundancy/n
  }
  return(output)
}

# Goodness of Fit
gof <- function(initial_model, result){
  output <- list()
  output$gof = sqrt(communality(initial_model, result)$mean_communality * rSquared(initial_model, result)$mean_r_squared)
  
  return(output)
}

# indicator significance based on communality and ave
# returns insignificant manifest variables
outerLoadingsRelevanceTest <- function(initial_model, result){
  output <- list()
  for(lv in initial_model$LV){
    output[[lv]] <- list()
    output[[lv]]$low_communalites <- initial_model$blocks[[lv]]$MV[which(result$loadings[as.vector(initial_model$blocks[[lv]]$MV),lv] < 0.4)]
    output[[lv]]$ave_candidates <- initial_model$blocks[[lv]]$MV[which(result$loadings[as.vector(initial_model$blocks[[lv]]$MV),lv] < 0.7)]
    temp <- vector()
    for(mv in output[[lv]]$ave_candidates){
      if(recalculate_ave(initial_model, result, mv, lv) < 0.5){
        temp = c(temp, T)
      } else {temp = c(temp, F)}
    }
    output[[lv]]$ave_filtered = output[[lv]]$ave_candidates[temp]
  }
  return(output)
}

# eigenvalue analysis to check that first eigenvalue is "much larger" than 1
  eigenvalueAnalysis <-function(data, initial_model){
  highest_eigens <- list()
  for (lv in initial_model$LV) {
    # TODO: correct Mode selected??
    if(initial_model$blocks[[lv]]$mode == "A"){
      highest_eigens[[lv]] <- eigen(cov(data[,as.vector(initial_model$blocks[[lv]]$MV)]))$values[1]
    }
  }
  return(highest_eigens)
}

# check that no indivator has higher loadings on another LV than on the LV it intends to measure
crossLoadingsValidation <- function(result, initial_model){
  checks <- list()
  for(lv in initial_model$LV){
    # TODO: correct Mode selected??
    if(initial_model$blocks[[lv]]$mode == "A"){
      for(mv in initial_model$blocks[[lv]]$MV){
        checks[[mv]] = all(result$crossloadings[mv, lv] >= result$crossloadings[mv, ])
      }
    }
  }
  return(checks)
}

# PCA analysis for formative blocks
#TODO: compare results of Dimitri(internet) against teacher
multiCollinearityPCAAnalysis <- function(data, initial_model){
  measures <- list()
  for(lv in initial_model$LV){
    # TODO: correct Mode selected??
    if(initial_model$blocks[[lv]]$mode == "B"){
      eigens = eigen(cov(data[,as.vector(initial_model$blocks[[lv]]$MV)]))$values
      #measures[[lv]] <- sqrt(head(eigens, 1)/tail(eigens, 1))
      #measures[[lv]] <- eigens[1]/eigens[2]
      measures[[lv]] <- head(eigens, 1)/tail(eigens, 1)
    }
  }
  return(measures)
}

# Help Function
recalculate_ave <- function(initial_model, result, mv, lv){
  output <- list()
  squared_loadings = result$loadings[,lv]^2
  return(mean(squared_loadings[as.vector(initial_model$blocks[[lv]]$MV[which(as.vector(initial_model$blocks[[lv]]$MV)!=mv)])]))
}

##########################################################
#             Bootstrapping & Param Testing              #
##########################################################

### BOOTSTRAPPING & HYPOTHESIS TESTING
bootstrap <- function (data, sample_number=500, sample_size=nrow(data)) {
  # n° samples
  B = sample_number
  # smaple size
  N = sample_size
  
  bootstrap_model <- list()
  for (i in 1:B) {
    sample <- round(runif(N, 1, dim(data)[1]),0)
    sample_data <- data[sample,]
    bootstrap_model[[i]] <- my.pls(sample_data, initial_model, inner_approximation_method = pathWeighting, threshold = 0.000000001)
  }
  return(bootstrap_model)
}

# Hypothesis: beta-values = 0
testSignificance <- function (model_list, significance_level = 0.05) {
  #test parameter for significance level...
  
  attribute = "path_coefficients"
  num_elements <- length(model_list)
  avg_model <- calculateAverage(model_list, attribute)
  std_dev_model <- calculateStdDev(model_list, attribute)
  # TODO: num_elements or sample_size?
  sig_score <- avg_model/std_dev_model
  sig_levels <- list()
  significance <- matrix (c(F), ncol=ncol(sig_score), nrow=nrow(sig_score), dimnames = list(rownames(sig_score), colnames(sig_score)))
  for (ccol in 1:ncol(avg_model)) {
    param_number <- sum(avg_model[,ccol]!=0)
    sig_level <- qt(c(significance_level/2, 1-significance_level/2), df=num_elements-param_number)[2]
    sig_levels[[colnames(avg_model)[ccol]]] <- sig_level
    ############### compare model to boostrap ###################
    significance[,ccol] <- abs(sig_score[,ccol]) > sig_level
  }
  
  results <- list()
  results$significance <- significance
  results$sig_score <- sig_score
  results$sig_level <- sig_levels
  return (results)
}


plotBootstrap <- function (bootstrap_model, attribute="path_coefficients") {
  avg <- calculateAverage(bootstrap_model, attribute)
  std <- calculateStdDev(bootstrap_model, attribute)
  connections <- sum(avg!=0)
  # Density plot
  par(mfrow=c(3, 4))
  for (i in 1:nrow(avg)) {
    for (j in 1:nrow(avg)) {
      values <-  c()
      for (model in 1:length(bootstrap_model)) {
        values [model] <- bootstrap_model[[model]][[attribute]][i, j]
      }
      if (any(values!= 0.0)) {
        d <- density(values)
        plot(d, type="n", main=paste(colnames(avg)[i],"->",colnames(avg)[j]))
        polygon(d, col="gray", border="black")
        abline(v = avg[i,j], col = "blue", lwd = 2)
        abline(v = avg[i,j] + std[i,j], col = "blue", lwd = 1)
        abline(v = avg[i,j] - std[i,j], col = "blue", lwd = 1)
      }
    }
  }

  return ()
}

## calculate average of multiple models from list
calculateAverage <- function(model_list, attribute) {
  num_elements <- length(model_list)
  sum_model <- model_list[[1]][[attribute]]
  for (i in 2:length(model_list)) {
    sum_model <- model_list[[i]][[attribute]] + sum_model
  }
  avg_model <- sum_model/num_elements
  return(avg_model)
}

## calculate std dev of multiple models from list
calculateStdDev <- function (model_list, attribute) {
  return(sqrt(calculateVariance(model_list, attribute)))
}

## calculate variance of multiple models from list
calculateVariance <- function (model_list, attribute) {
  num_elements <- length(model_list)
  avg_model <- calculateAverage(model_list, attribute)
  sum_square_diff <- (model_list[[1]][[attribute]]-avg_model)^2
  for (i in 2: length(model_list)) {
    sum_square_diff <- sum_square_diff + (model_list[[i]][[attribute]]-avg_model)^2
  }
  # TODO: Should we divide by num_elements-1 because it's a sample variance?
  # sample variance
  var_model <- sum_square_diff/(num_elements-1)
  return (var_model)
}

##########################################################
#                      Consistent PLS                    #
##########################################################

consistentPLS <- function(data, initial_model, inner_approximation_method = pathWeighting, threshold = 0.000000001){
  result <- my.pls(data, initial_model, inner_approximation_method = pathWeighting, threshold = 0.000000001)
  output <- list()
  output$rhos = vector(length=length(initial_model$LV))
  names(output$rhos) = initial_model$LV
  
  # 1) determine (inconsistent) latent variable correlations
  r = cor(result$scores)
  
  for(lv in initial_model$LV){
    # w seems correct
    w = result$weights[as.vector(initial_model$blocks[[lv]]$MV), lv]
    # s seems correct
    S = cov(data)[as.vector(initial_model$blocks[[lv]]$MV), as.vector(initial_model$blocks[[lv]]$MV)]
    
    #rhos seem correct
    output$rhos[[lv]] = (t(w)%*%w)^2*((t(w) %*% (S - diag(diag(S))) %*% w)/(t(w) %*% (w%*%t(w) - diag(diag(w%*%t(w)))) %*% w))
    # set rho = 1 if mode == B
    if (model$blocks[[lv]]$mode != "A") output$rhos[[lv]] <- 1
  }
  output$consist_cor = r/sqrt(output$rhos%*%t(output$rhos))
  diag(output$consist_cor) <- diag(length(lv))
  
  #consistent loadings
  
  output$consist_loadings = matrix(c(0), nrow = length(initial_model$MV), ncol = length(initial_model$LV), dimnames = list(initial_model$MV, initial_model$LV))
  
  for(lv in initial_model$LV){
    w = result$weights[as.vector(initial_model$blocks[[lv]]$MV), lv]
    rho = output$rhos[[lv]]
    consist_loadings = w * as.numeric(sqrt(rho)/t(w)%*%w)
    names(consist_loadings) = initial_model$blocks[[lv]]$MV
    for(mv in initial_model$blocks[[lv]]$MV){
      output$consist_loadings[mv, lv] = consist_loadings[mv]
    }
  }
  
  #consistent path coefficients
  output$consist_path_coeffs = matrix(c(0), nrow = length(initial_model$LV), ncol = length(initial_model$LV), dimnames = list(initial_model$LV, initial_model$LV))
  for(lv in initial_model$LV){
    predecessors <- initial_model$blocks[[lv]]$predecessors
    if(length(predecessors) > 0) {
      
      r_pred = output$consist_cor[as.vector(predecessors), as.vector(predecessors)]
      r_pred_inv = solve(r_pred)
      r_lv_pred = output$consist_cor[as.vector(predecessors), lv]
      
      consist_path_coeffs = as.vector(r_pred_inv %*% r_lv_pred)
      names(consist_path_coeffs) = predecessors
      
      for(lv2 in predecessors){
        output$consist_path_coeffs[lv2, lv] = consist_path_coeffs[lv2]
      }
    }
  }
  
  #Consistent R squares
  n = 0
  tot=0
  
  for (lv in rownames(initial_model$inner_model)) {
    regressors <- initial_model$blocks[[lv]]$predecessors
    if(length(regressors) > 0){
      n = n + 1
      vector = result$scores[,lv]
      data_reg = result$scores[,regressors]
      predicted = data_reg %*% as.matrix(output$consist_path_coeffs[regressors,lv])
      m = mean(vector)
      ss_tot = sum((vector - m)^2)
      ss_err = sum((vector - predicted)^2)
      output$consistent_r_squared[[lv]] <- 1 - ss_err/ss_tot
      tot = tot + output$consistent_r_squared[[lv]]
    }
  }
  
  output$consistent_mean_r_squared = tot/n
  
  return(output)
}

## @knitr part2
##########################################################
#                       MAIN                             #
##########################################################
  # User parameters:
    wd_directory <- "C:\\Users\\jonat\\OneDrive for Business\\Dokumente\\Master-Studium\\Advanced Analytics, NOVA IMS\\Kurse MAA2016\\2_Descriptive Analytics\\Project\\workspace"
    file_name <- "bank.csv"
    inner_model_file <- "inner_model_mappings.txt"
    outer_model_file <- "outer_model_mappings.txt"     #formative: "outer_model_mappings_formative.txt"
    convergence_threshold <- 0.000000001
    # 1. imputeDeletion 2. imputeMean 3. imputeNearestNeighbour 4. imputeSpecialMean
    imputation_method <- imputeDeletion
    # 1. centroidWeighting 2. factorialWeighting 3. pathWeighting
    inner_approx_method <- centroidWeighting
    bootstrap_sample_number <- 500
    # significance level for t-test
    significance_level <- 0.05
  #
  setwd(wd_directory)
    
  raw_data <- read.csv(file_name, header = T, sep=",", row.names = 1)
  inner_model_mappings <- read.table(inner_model_file, header = T, sep=",")
  outer_model_mappings <- read.table(outer_model_file, header = T, sep=",")
  
  # imputation method to be passed as parameter (standard: imputeDeletion method)
  data <- dataTreatment(raw_data, impute_method = imputation_method)
  initial_model <- my.pls.model(data, inner_model_mappings = inner_model_mappings, outer_model_mappings = outer_model_mappings)
  
  # inner_approx standard: centroid, weighting thresshold: 0.00001
  model <- my.pls(data, initial_model, inner_approximation_method = inner_approx_method, threshold = convergence_threshold)
  
  #Bootstrapping Part1
  bootstrap_sample_size <- dim(data)[1]
  bootstrap_model <- bootstrap (data, sample_number = bootstrap_sample_number)

## @knitr part3
  #Bootstrapping Part2
  bootstrap_plot <- plotBootstrap(bootstrap_model)
  significance <- testSignificance(bootstrap_model, significance_level)
  
  #Assessment Measures
  assessment <- assessmentMeasures(data, initial_model, model)
  
  #Build consistent model
  consistent_model <- consistentPLS(data, initial_model)

  