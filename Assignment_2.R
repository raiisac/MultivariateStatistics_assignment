update.packages("rlang")
library(xfun)
library(class)#for the knn.cv function
library(tidyverse)
library(rprojroot)
library(kableExtra)
library(MASS)#lda function
library(biotools)
library(ISLR2)
library(patchwork)
library(randomForest)
library(viridis)
library(GPArotation)
library(psych) #EFA
library(maptools)
library(lavaan)
library(HDclassif)
library(mclust)
library(fdm2id)

# this package does not work for the newest version of R -->
# only load if the local data files are not there
if (!file.exists(find_root_file("data/rect_constr.csv",
     criterion = has_file("MultivariateStatistics_assignment.Rproj")))) {
  library(smacof)
}

knitr::opts_chunk$set(echo = TRUE)

#'
#'
#' All R scripts and the data can be found on [this GitHub repository](https://github.com/raiisac/MultivariateStatistics_assignment).
#'
#' # Task 1
## ----load, include=FALSE, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------------------------
#load the data and functions
data(College)

#test the difference between centroids
Wilks_manova <- function(data, group = "Private"){
  lm.out <- lm(cbind(
    sapply(X = colnames(data)[colnames(data) != group],
           FUN = function(x) get(x)))~
       get((group)),
   data = data)
  s <- summary(manova(lm.out), test = c("Wilks"))
  return(s$stats[1, "Pr(>F)"])#return the pvalue
}
#test the differences in var-cov matrix between groups
boxM_test <- function(data, group = "Private"){
  out <- boxM(data[,colnames(data) != group], data[,group])
  return(out$p.value)#return the pvalue
}

#execute lda with LOOCV and save the hit rate, sensitivity, and specificity
lda_analysis <- function(data, group = "Private", datatype){
  K <- length(unique(levels(data[,group])))
  positivelevel <- levels(data[,group])[K]# by default, it is assumed the last level is the positive level
  negativelevel <-  levels(data[,group])[-K]
  lda.out <- lda(
    formula(paste0(group, "~", paste(colnames(data)[colnames(data) != group],
                                     collapse = "+"))),
    data = data,
    prior = rep(1/K, K),
    CV = TRUE)
  tab <- table(data[,group], lda.out$class)
  out <- data.frame(hitrate = sum(diag(tab))/sum(tab),
                    sensitivity = tab[positivelevel, positivelevel] /
                      sum(tab[positivelevel,]),
                    specificity = tab[negativelevel, negativelevel] /
                      sum(tab[negativelevel,]),
                    data = datatype,
                    method = "LDA"
  )
  return(out)
}
#execute qda with LOOCV and save the hit rate, sensitivity, and specificity
qda_analysis <- function(data, group = "Private", datatype){
  K <- length(unique(levels(data[,group])))
  positivelevel <- levels(data[,group])[K]# by default, it is assumed the last level is the positive level
  negativelevel <-  levels(data[,group])[-K]
  qda.out <- qda(
    formula(paste0(group, "~", paste(colnames(data)[colnames(data) != group],
                                     collapse = "+"))),
    data = data,
    prior = rep(1/K, K),
    CV = TRUE)
  tab <- table(data[,group], qda.out$class)
  out <- data.frame(hitrate = sum(diag(tab))/sum(tab),
                    sensitivity = tab[positivelevel, positivelevel] /
                      sum(tab[positivelevel,]),
                    specificity = tab[negativelevel, negativelevel] /
                      sum(tab[negativelevel,]),
                    data = datatype,
                    method = "QDA"
  )
  return(out)
}

knn_analysis <-  function(data, group = "Private", datatype, testk){
  #In this function, we test a knn approach with several different values of K and we choose the one with the highest hit rate. If there are ties for the best hit rate, we choose the lowest K
  #testk is a vector for all values of K that should be tried.
  K <- length(unique(levels(data[,group])))
  positivelevel <- levels(data[,group])[K]# by default, it is assumed the last level is the positive level
  negativelevel <-  levels(data[,group])[-K]
  output <- data.frame(hitrate = NA,
                    sensitivity = NA,
                    specificity = NA,
                    data = datatype,
                    method = "KNN",
                    K = testk)
  for (x in testk) {
    knn.out <- knn.cv(train = data[,colnames(data) != group],
                      cl = data[, group],
                      k = x,
                      prob = TRUE
                      )
    tab <- table(data[,group], unname(knn.out))
    out <- data.frame(hitrate = sum(diag(tab))/sum(tab),
                    sensitivity = tab[positivelevel, positivelevel] /
                      sum(tab[positivelevel,]),
                    specificity = tab[negativelevel, negativelevel] /
                      sum(tab[negativelevel,]),
                    data = datatype,
                    method = "KNN",
                    K = x
    )
  output[x, ] <- out
  }
  graph <- output %>%
    pivot_longer(cols = c("hitrate","sensitivity", "specificity"),
                 names_to = "performance.measure", values_to = "value") %>%
    ggplot() +
    geom_line(aes(color = `performance.measure`, x = K, y = value)) +
    ylab("Performance") + theme_bw() + ylim(0.6, 1)
  selected <- which(output$hitrate == max(output$hitrate))[1]
  return(list(graph,output[selected,]))
}

bagging_analysis <- function(data, group = "Private", datatype, ntree){
  K <- length(unique(levels(data[,group])))
  positivelevel <- levels(data[,group])[K]# by default, it is assumed the last level is the positive level
  negativelevel <-  levels(data[,group])[-K]
  nvar <- ncol(data)-1
  bagging.out <- randomForest(
    formula(paste0(group, "~", paste(colnames(data)[colnames(data) != group],
                                     collapse = "+"))),
    data = data,
    mtry = nvar,
    ntree = ntree,
    importance = TRUE)
  tab <- table(data[,group], bagging.out$predicted) #the predicted values of the input data based on out-of-bag samples.
  out <- data.frame(hitrate = sum(diag(tab))/sum(tab),
                    sensitivity = tab[positivelevel, positivelevel] /
                      sum(tab[positivelevel,]),
                    specificity = tab[negativelevel, negativelevel] /
                      sum(tab[negativelevel,]),
                    data = datatype,
                    method = "bagging"
    )
  graph <- bagging.out$importance %>%
    as.data.frame %>%
    dplyr::select(MeanDecreaseGini) %>%
    mutate(Var = rownames( bagging.out$importance )) %>%
    ggplot() + geom_point(aes(x=MeanDecreaseGini, y=Var)) +
  scale_y_discrete(limits = names(sort(bagging.out$importance[,4]))) +
    theme_bw() + ylab("")
  return(list(graph, out))
}

RF_analysis <- function(data, group = "Private", datatype, ntree, mtry){
  K <- length(unique(levels(data[,group])))
  nvar <- ncol(data)-1
  positivelevel <- levels(data[,group])[K]# by default, it is assumed the last level is the positive level
  negativelevel <-  levels(data[,group])[-K]
  if (mtry != "optimize") {
    bagging.out <- randomForest(
      formula(paste0(group, "~", paste(colnames(data)[colnames(data) != group],
                                     collapse = "+"))),
      data = data,
      mtry = mtry,
      ntree = ntree,
      importance = TRUE)
    tab <- table(data[,group], bagging.out$predicted) #the predicted values of the input data based on out-of-bag samples.
    out <- data.frame(hitrate = sum(diag(tab))/sum(tab),
                    sensitivity = tab[positivelevel, positivelevel] /
                      sum(tab[positivelevel,]),
                    specificity = tab[negativelevel, negativelevel] /
                      sum(tab[negativelevel,]),
                    data = datatype,
                    method = "randomForest"
      )
    graph <- bagging.out$importance %>%
      as.data.frame %>%
      dplyr::select(MeanDecreaseGini) %>%
      mutate(Var = rownames( bagging.out$importance )) %>%
      ggplot() + geom_point(aes(x=MeanDecreaseGini, y=Var)) +
    scale_y_discrete(limits = names(sort(bagging.out$importance[,4]))) +
      theme_bw() + ylab("")
    return(list(graph, out))

  }else{
    output <- data.frame(hitrate = NA,
                    sensitivity = NA,
                    specificity = NA,
                    data = datatype,
                    method = "randomForest",
                    mtry = 1:nvar)
    for (mtry in 1:nvar) {
    bagging.out <- randomForest(
      formula(paste0(group, "~", paste(colnames(data)[colnames(data) != group],
                                     collapse = "+"))),
      data = data,
      mtry = mtry,
      ntree = ntree,
      importance = TRUE)
    tab <- table(data[,group], bagging.out$predicted) #the predicted values of the input data based on out-of-bag samples.
    out <- data.frame(hitrate = sum(diag(tab))/sum(tab),
                    sensitivity = tab[positivelevel, positivelevel] /
                      sum(tab[positivelevel,]),
                    specificity = tab[negativelevel, negativelevel] /
                      sum(tab[negativelevel,]),
                    data = datatype,
                    method = "randomForest",
                    mtry = mtry
      )
  output[mtry, ] <- out
  }
  graphopt <- output %>%
    pivot_longer(cols = c("hitrate","sensitivity", "specificity"),
                 names_to = "performance.measure", values_to = "value") %>%
    ggplot() +
    geom_line(aes(color = `performance.measure`, x = mtry, y = value)) +
    ylab("Performance") + theme_bw() + ylim(0.6, 1)
  selected <- which(output$hitrate == max(output$hitrate))[1]
  bagging.out <- randomForest(
      formula(paste0(group, "~", paste(colnames(data)[colnames(data) != group],
                                     collapse = "+"))),
      data = data,
      mtry = output[selected,'mtry'],
      ntree = ntree,
      importance = TRUE)
  graph <- bagging.out$importance %>%
      as.data.frame %>%
      dplyr::select(MeanDecreaseGini) %>%
      mutate(Var = rownames( bagging.out$importance )) %>%
      ggplot() + geom_point(aes(x=MeanDecreaseGini, y=Var)) +
    scale_y_discrete(limits = names(sort(bagging.out$importance[,4]))) +
      theme_bw() + ylab("")
  return(list(graph, output[selected, ], graphopt))
  }
}

#'
#' The data consists of `r ncol(College)` variables of which the categorical variable *Private* denotes whether the school is private (*Private == Yes*) or not (*Private == No*). There are `r nrow(College)` schools in the data; `r sum(College$Private == "Yes")` (`r round(sum(College$Private == "Yes")/nrow(College)*100,2)`%) private schools and `r sum(College$Private == "No")` (`r round(sum(College$Private == "No")/nrow(College)*100,2)`%) non-private schools.
#'
#' A detailed outline of all functions that are used for this task can be found in the appendix and on [this GitHub repository](https://github.com/raiisac/MultivariateStatistics_assignment). Before we start to build models to distinguish private and non-private schools, Wilk's lambda test is performed where the null hypothesis states that the centroid (mean) in both groups is the same. Both in the log-transformed and the original (centered or standardized) data, this null hypothesis is rejected. This means that the centroids of private and non-private schools differ significantly.
#'
#' The test of Box to test for equality of within-group covariance matrices shows that $H_0$ of equal covariance matrices across groups is not supported by the data. However, we still report on the results from LDA (which assumes equal covariances) since the alternative, QDA, does not always perform better because the lower bias of the QDA classifier may not outweigh its higher model complexity.
#'
#' Table \@ref(tab:compareresults) shows all results for the different performance measures (hit rate, sensitivity and specificity) and the different datasets. Overall, sensitivity (proportion of positive observations that are being correctly classified) tends to be higher than specificity (proportion of negative observations that are being correctly classified) meaning it's generally easier to detect private schools than non-private schools.
#'
#' As expected, it does not matter whether the data is standardized or centered for the QDA and LDA models and log-transforming skewed variables yields a higher hit rate in LDA and QDA. The hit rate is higher in LDA than in QDA which suggests that the true decision boundary is likely close to linear. The extra complexity of QDA is not worthwhile.
#'
#' It is curious that KNN works so badly on the centered data with log-transformed skewed variables.
#' The highest overall hit rate, for all datasets, is achieved in the random forest models. Bagging achieves approximately the same performance as LDA although it sacrifices specificity for higher sensitivity.
#'
#'
## ----classification, include=TRUE, echo = TRUE, message=FALSE, warning=FALSE, cache = TRUE-----------------------------------------------------------------
College_logtranfskewed <- College %>%
  mutate(Apps = log(Apps),
         Accept = log(Accept),
         Enroll = log(Enroll),
         Top10perc = log(Top10perc),
         F.Undergrad = log(F.Undergrad),
         P.Undergrad = log(P.Undergrad),
         Books = log(Books),
         Personal = log(Personal),
         Expend = log(Expend))
College_centered <- College %>%
  dplyr::select(-Private) %>%
  scale(center = TRUE, scale = FALSE) %>%
  as.data.frame() %>%
  mutate(Private = College$Private)
College_logtranfskewed_centered <- College_logtranfskewed %>%
  dplyr::select(-Private) %>%
  scale(center = TRUE, scale = FALSE) %>%
  as.data.frame() %>%
  mutate(Private = College$Private)
College_std <- College %>%
  dplyr::select(-Private) %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  mutate(Private = College$Private)
College_logtranfskewed_std <- College_logtranfskewed %>%
  dplyr::select(-Private) %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame() %>%
  mutate(Private = College$Private)
testlist <- list(center = College_centered,
                 center_logtransfskewed = College_logtranfskewed_centered,
                 std = College_std,
                 std_logtransfskewed = College_logtranfskewed_std
                 )
#test the difference between centroids
wilks_results <- testlist %>% purrr::map(function(x) Wilks_manova(x, "Private"))
#In each dataset, H0:mu_{yes} = mu_{no} is rejected -->
# centroids of Private and non-private schools differ significantly
#test the difference between centroids
boxm_results <- testlist %>% purrr::map(function(x) boxM_test(x, "Private"))
#LDA
lda_results <- map2(testlist, names(testlist),
                    function(x,y){lda_analysis(x,"Private", y)}) %>%
  do.call("rbind", .)
#QDA
qda_results <- map2(testlist, names(testlist),
                    function(x,y){qda_analysis(x,"Private", y)}) %>%
  do.call("rbind", .)
#KNN
knn_output <- map2(testlist, names(testlist),
                    function(x,y){knn_analysis(x,"Private", y, 1:100)})
knn_results <- lapply(knn_output, function(element){
  element[[2]] # The first element contains graphs
})
knn_results <- knn_results %>%
  do.call("rbind", .)
#bagging
bagging_output <- map2(testlist, names(testlist),
                    function(x,y){bagging_analysis(x, "Private", y,
                                              ntree = 2000)})
bagging_results <- lapply(bagging_output, function(element){
  element[[2]] # The first element contains graphs
})
bagging_results <- bagging_results %>%
  do.call("rbind", .)
#random forests
rf_output <- map2(testlist, names(testlist),
                    function(x,y){RF_analysis(x, "Private", y,
                                              ntree = 2000, mtry = "optimize")})
rf_results <- lapply(rf_output, function(element){
  element[[2]] # The first element contains graphs
})
rf_results <- rf_results %>%
  do.call("rbind", .)

#'
## ----compareresults, echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
as.data.frame(lda_results) %>%
  bind_rows(as.data.frame(qda_results)) %>%
  bind_rows(as.data.frame(knn_results)) %>%
  bind_rows(as.data.frame(bagging_results)) %>%
  bind_rows(as.data.frame(rf_results)) %>%
  dplyr::select(-K, -mtry) %>%
  pivot_wider(names_from = data,
              values_from = c(hitrate, sensitivity, specificity)) -> A
A %>% knitr::kable(
  caption = 'A performance comparison of the different methods',
  digits = 3,
  col.names = c("Method",rep(c("Centered", "Centered Logtransf",
                   "Standardized", "Standardized Logtransf"), 3)),
  booktabs = TRUE,
  valign = 't'
) %>%
  row_spec(0, angle = 90) %>%
  add_header_above(c(" ", "Hit rate" = 4, "Sensitivity" = 4, "Specificity" = 4)) %>%
  column_spec(2, color = "white",
              background = spec_color(unlist(A[,2]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(3, color = "white",
              background = spec_color(unlist(A[,3]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(4, color = "white",
              background = spec_color(unlist(A[,4]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(5, color = "white",
              background = spec_color(unlist(A[,5]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(6, color = "white",
              background = spec_color(unlist(A[,6]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(7, color = "white",
              background = spec_color(unlist(A[,7]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(8, color = "white",
              background = spec_color(unlist(A[,8]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(9, color = "white",
              background = spec_color(unlist(A[,9]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(10, color = "white",
              background = spec_color(unlist(A[,10]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(11, color = "white",
              background = spec_color(unlist(A[,11]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(12, color = "white",
              background = spec_color(unlist(A[,12]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13])))) %>%
  column_spec(13, color = "white",
              background = spec_color(unlist(A[,13]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(A[,2:13]),
                                                     max(A[,2:13]))))

#'
#'
#'
#' For the K-nearest neighbour (KNN) approach, the function knn_cv from the *class* package is used to execute KNN with LOOCV and test values of K between 1 and 100. The model with the highest hit rate is chosen as the final model in Table \@ref(tab:compareresults). If there are multiple models with the same hit rate, the smallest K is chosen, Figure \@ref(fig:KNNgraphs) shows the results in each of the datasets. The chosen K is `r knn_results[1,'K']` for the centered data, `r knn_results[2,'K']` for the centered data where the skewed variables are log-transformed, `r knn_results[3,'K']` for the standardized data, `r knn_results[4,'K']` for the standardized data where the skewed variables are log-transformed.
#'
## ----KNNgraphs, echo=FALSE, include=TRUE, message=FALSE, warning=FALSE, fig.width = 8, fig.height=5, fig.cap="Evolution of hit rate, sensitivity and specificity for different K in the KNN approach. A shows the centered data, B the centered data with log-transformed skewed variables, C shows the standardized data, and D the standardized data with log-transformed skewed variables."----
knn_graphs <- lapply(knn_output, function(element){
  element[[1]] # The first element contains graphs
  })
knn_graphs[[1]] + knn_graphs[[2]] + knn_graphs[[3]] + knn_graphs[[4]] +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")

#'
#' For the bagging and  Random forest (RF) approach, the function randomForest from the *randomForest* package is used. 2000 bootstrapped data and trees are drawn in all models.  Figure \@ref(fig:bagginggraphs) and \@ref(fig:rfgraphs) show the importance of different variables in the final models. *F.Undergrad* and *outstate* are the most important one in all bagging and RF models. The difference in MeanDecreaseGini (the mean decrease in the Gini index) between the most important variables (*F.Undergrad* and *Outstate*) and the other variables is much larger in the bagging model than in the RF models. *Enroll*, which is not at all important in the bagging models, is consistently the third most important variable in the RF models.
#'
#' For the RF results, we try different values for the mtry parameter which is the number of variables randomly sampled as candidates at each split. The model with the highest hit rate is chosen as the final model in Table \@ref(tab:compareresults). If there are multiple models with the same hit rate, the smallest mtry is chosen, Figure \@ref(fig:rfgraphsopt) shows the results in each of the datasets. The chosen value for mtry is `r rf_results[1,'mtry']` for the centered data (Figure \@ref(fig:rfgraphsopt)A), `r rf_results[2,'mtry']` for the centered data where the skewed variables are log-transformed (Figure \@ref(fig:rfgraphsopt)B), `r rf_results[3,'mtry']` for the standardized data (Figure \@ref(fig:rfgraphsopt)C), `r rf_results[4,'mtry']` for the standardized data where the skewed variables are log-transformed (Figure \@ref(fig:rfgraphsopt)D).
#'
## ----bagginggraphs, echo=FALSE, include=TRUE, message=FALSE, warning=FALSE, fig.width = 8, fig.height=5, fig.cap="Variance importance plots for the bagging models. A shows the centered data, B the centered data with log-transformed skewed variables, C shows the standardized data, and D the standardized data with log-transformed skewed variables."----
bagging_graphs <- lapply(bagging_output, function(element){
  element[[1]] # The first element contains graphs
  })
bagging_graphs[[1]] + bagging_graphs[[2]] + bagging_graphs[[3]] + bagging_graphs[[4]] +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")

#'
## ----rfgraphs, echo=FALSE, include=TRUE, message=FALSE, warning=FALSE, fig.width = 8, fig.height=5, fig.cap="Variance importance plots for the random forest models. A shows the centered data, B the centered data with log-transformed skewed variables, C shows the standardized data, and D the standardized data with log-transformed skewed variables."----
rf_graphs <- lapply(rf_output, function(element){
  element[[1]] # The first element contains graphs
  })
rf_graphs[[1]] + rf_graphs[[2]] + rf_graphs[[3]] + rf_graphs[[4]] +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")

#'
## ----rfgraphsopt, echo=FALSE, include=TRUE, message=FALSE, warning=FALSE, fig.width = 8, fig.height=5, fig.cap="Optimization of the mtry parameter for the random forest models. A shows the centered data, B the centered data with log-transformed skewed variables, C shows the standardized data, and D the standardized data with log-transformed skewed variables."----
rf_graphs <- lapply(rf_output, function(element){
  element[[3]] # The first element contains graphs
  })
rf_graphs[[1]] + rf_graphs[[2]] + rf_graphs[[3]] + rf_graphs[[4]] +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")

#'
#' \pagebreak
#'
#'
#' # Task 2
#'
## ----load2, include=TRUE, message=FALSE, warning=FALSE, cache = TRUE---------------------------------------------------------------------------------------
load(find_root_file("data/fashion.Rdata",
     criterion = has_file("MultivariateStatistics_assignment.Rproj")))
#centered variables
zshopping <- scale(train.data, center = TRUE, scale = FALSE)
##### HDclassif ######
set.seed(20)
#common dimension hddc() models "AkjBkQkD" and "AkjBQkD" with
# com_dim=2,3,4,5, or 6 fitted on the training data.
names <- list()
ari_values <- list()
bic <- list()
for (model in c("AkjBkQkD", "AkjBQkD")) {
  for (dim in 2:6)
  {
    hddc.out <- hddc(zshopping, K = 3, model = model, com_dim = dim,
                     d_select = "BIC", init = "kmeans")
    ari <- adjustedRandIndex(hddc.out$class, train.target)
    model_name = paste(model,"_",toString(dim), sep = "")
    names <- append(names, model_name)
    ari_values <- append(ari_values, ari)
    bic <- append(bic, hddc.out$BIC)
  }
}

###### mclust model ######
set.seed(201)
#Mclust() models "VVE","VEV","EVV","VVV" fitted on the first 2,3,4,5, or 6
#unstandardized principal components extracted from the training data.
prcomp.out <- prcomp(zshopping)
comp <- zshopping %*% prcomp.out$rotation
names2 <- list()
ari_values2 <- list()
bic2 <- list()
for (model in c("VVE","VEV","EVV","VVV")) {
  for (c in 2:6)
  {
    mclust.out <- Mclust(comp[,1:c], G = 3, modelNames = model)
    ari <- adjustedRandIndex(mclust.out$class, train.target)
    model_name = paste(model, "_", toString(c), sep = "")
    names2 <- append( names2, model_name)
    ari_values2 <- append(ari_values2, ari)
    bic2 <- append(bic2,mclust.out$bic)
  }
}
###### the best model (highest ARI value) selected is AkjBQkD_4 ######
hddc.out <- hddc(zshopping, K = 3, model = "AkjBQkD", com_dim = 4, d_select = "BIC",
               init = "kmeans")
ari <- adjustedRandIndex(hddc.out$class, train.target)
tab <- table(hddc.out$class, train.target)
#switch class labels to have maximum correspondance
mapping <- mapClass(hddc.out$class, train.target)
#compute principal components
loading <- prcomp.out$rotation %*% diag(prcomp.out$sd)
Z <- predict(prcomp.out, zshopping)

#'
## ----showhddc, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------
#compare models
compare_hddc <- do.call(rbind, Map(data.frame, A = names, B = ari_values, C = bic))
compare_hddc %>% kable(caption = 'A comparison of HDDC models',
  digits = 3,
  col.names = c("Method_dimensions","ARI", "BIC"),
  booktabs = TRUE) %>%
  column_spec(2, color = "white",
              background = spec_color(unlist(compare_hddc[,2]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(compare_hddc[,2]),
                                                     max(compare_hddc[,2])))) %>%
  column_spec(3, color = "white",
              background = spec_color(unlist(compare_hddc[,3]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(compare_hddc[,3]),
                                                     max(compare_hddc[,3]))))

#'
## ----showmclust, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------------------------------------------
compare_mclust <- do.call(rbind, Map(data.frame, A = names2, B = ari_values2, C = bic2))
compare_mclust %>% kable(caption = 'A comparison of mclust models',
  digits = 3,
  col.names = c("Method_dimensions","ARI", "BIC"),
  booktabs = TRUE) %>%
  column_spec(2, color = "white",
              background = spec_color(unlist(compare_mclust[,2]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(compare_mclust[,2]),
                                                     max(compare_mclust[,2])))) %>%
  column_spec(3, color = "white",
              background = spec_color(unlist(compare_mclust[,3]), option = "viridis",
                                      direction = -1,
                                      scale_from = c(min(compare_mclust[,3]),
                                                     max(compare_mclust[,3]))))

#'
#'
## ----bestplots, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "A comparison of the best models with the true, actual labels."----------
par(mfrow = c(1,2))
#plot hddc solution in the space of the first two principal components
plot(Z, xlab = "PC1", ylab = "PC2", main = "Best model")
points(Z[hddc.out$class == 1,], col = "red")
points(Z[hddc.out$class == 3,], col = "green")
points(Z[hddc.out$class == 2,], col = "blue")

#plot actual label in the space of the first two principal components
plot(Z, xlab = "PC1", ylab = "PC2", main = "Actual label")
points(Z[train.target == 7,], col = "red")
points(Z[train.target == 0,], col = "green")
points(Z[train.target == 1,], col = "blue")

#'
#' To evaluate the extent to which the model-based unsupervised clustering can recover true class labels, we use the Adjusted Rand Index (ARI) based on the true and predicted class labels. Model with a higher ARI is then considered as a better model.
#'
#' When fitting the common dimension hddc models on the training data, with a 3-cluster  model specified as "AkjBkQkD" or "AkjBQkD", dimensions from 2 to 6 (Table \@ref(tab:showhddc)), we can see that the best performing model is "AkjBQkD" with 4 dimensions (ARI = `r round(max(compare_hddc$B),3)`). It has a BIC of `r round(compare_hddc[compare_hddc$B == max(compare_hddc$B),"C"],3)`, which is however not the lowest value among all models. The model "AkjBQkD" with 2 dimensions has a BIC of `r round(compare_hddc[compare_hddc$C == min(compare_hddc$C),"C"],3)` and an ARI of `r round(compare_hddc[compare_hddc$C == min(compare_hddc$C),"B"],3)`.
#'
#' In addition, when fitting the Mclust models "VVE","VEV","EVV","VVV" on the first 2 to 6 unstandardized principal components (Table \@ref(tab:showmclust)), we can see that the best performing model is "VVE" fitted on 6 principal components (ARI = `r round(max(compare_hddc$B),3)`). Comparatively, the best hddc model outperforms the Mclust model in terms of ARI. Also, the BIC values of Mclust are all much higher than those of hdcc models.
#'
#' When visualizing the predicted class labels for the best clustering model in the space of the first two principal components, we can see that the space is divided into 3 parts colored each in red, green and blue (Figure \@ref(fig:bestplots)). There is little overlap among the clusters, providing a visual illustration that the clustering worked well. On the other hand, when visualizing the actual class labels in the space of the first two principal components, we see some clear overlaps between two clusters.
#'
#' \clearpage
#'
#' # Task 3
## ----load3, include=FALSE, message=FALSE, warning=FALSE, echo=FALSE----------------------------------------------------------------------------------------
if (file.exists(find_root_file("data/rect_constr.csv",
     criterion = has_file("MultivariateStatistics_assignment.Rproj")))) {
  rect_constr <- read.csv(
    find_root_file("data/rect_constr.csv",
                   criterion =
                     has_file("MultivariateStatistics_assignment.Rproj")))
} else {
  data(rect_constr)
}
if (file.exists(find_root_file("data/rectangles.csv",
     criterion = has_file("MultivariateStatistics_assignment.Rproj")))) {
  rectangles <- read.csv(
    find_root_file("data/rectangles.csv",
                   criterion =
                     has_file("MultivariateStatistics_assignment.Rproj")))
} else {
  data(rectangles)
}
library(smacof)
data(rectangles)
data(rect_constr)

#'
## ----smacofs, results='hide', echo=TRUE, message=FALSE, warning=FALSE--------------------------------------------------------------------------------------
#ratio
m1 <- smacofSym(delta = rectangles, ndim = 2, type = "ratio", init = "torgerson")
#interval
m2 <- smacofSym(delta = rectangles, ndim = 2, type = "interval", init = "torgerson")
#ordinal
m3 <- smacofSym(delta = rectangles, ndim = 2, type = "ordinal", init = "torgerson")
#spline
m4 <- smacofSym(delta = rectangles, ndim = 2, type = "mspline",
                spline.degree = 4, spline.intKnots = 4, init = "torgerson")
#stress-1 values
round(c(m1$stress, m2$stress, m3$stress, m4$stress), 3)

#'
## ----gofsm, results='hide', echo=TRUE, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------
#check goodness of fit using stress-norms and permutation test
#stress norm
set.seed(1)
rstress.m1 <- randomstress(n = 16, ndim = 2, nrep = 500, type = "ratio")
rstress.m2 <- randomstress(n = 16, ndim = 2, nrep = 500, type = "interval")
rstress.m3 <- randomstress(n = 16, ndim = 2, nrep = 500, type = "ordinal")
rstress.m4 <- randomstress(n = 16, ndim = 2, nrep = 500, type = "mspline")
#distribution of stress for random data
distr_stress1 <- mean(rstress.m1) - 2 * sd(rstress.m1)
distr_stress2 <- mean(rstress.m2) - 2 * sd(rstress.m2)
distr_stress3 <- mean(rstress.m3) - 2 * sd(rstress.m3)
distr_stress4 <- mean(rstress.m4) - 2 * sd(rstress.m4)
#permutation test
set.seed(1)
perm.m1 <- permtest(m1, nrep = 500)
perm.m2 <- permtest(m2, nrep = 500)
perm.m3 <- permtest(m3, nrep = 500)
perm.m4 <- permtest(m4, nrep = 500)
#stability of solution using jackknife
jack.car <- jackmds(m3)

#'
## ----MDSresults, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE------------------------------------------------------------------------------------
data.frame(
  method = c("ratio", "interval", "ordinal", "mspline"),
  stress1 = c(m1$stress, m2$stress, m3$stress, m4$stress),
  stressnorm = c(distr_stress1, distr_stress2, distr_stress3, distr_stress4),
  perm = c(mean(perm.m1$stressvec), mean(perm.m2$stressvec),
           mean(perm.m3$stressvec), mean(perm.m4$stressvec))
) %>%
  kable(booktabs = TRUE,
        caption = "Stress-1, mean-2*standard deviation for stress norm test, and mean for the permutation test.",
        col.names = c("measurement levels", "stress-1", "stress norm","permutation test"),
        digits = 3)

#'
## ----shephard, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "Residual plots and Shepard diagrams for ordinal and spline MDS."---------
par(mfrow = c(2,2))
plot(m3,plot.type = "resplot", main = "residual plot ordinal MDS")
plot(m3,plot.type = "Shepard", main = "Shepard diagram ordinal MDS")
plot(m4,plot.type = "resplot", main = "residual plot spline MDS")
plot(m4,plot.type = "Shepard", main = "Shepard diagram spline MDS")

#'
## ----ordinalMDS, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "Configuration ordinal MDS."--------------------------------------------
#configuration ordinal MDS
plot(m3, plot.type = "conf", main = "configuration ordinal MDS", pch = '',
     label.conf = list(label = TRUE, pos = 3, col = 1, cex = 1.5))

#'
## ----distrstress, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "Distribution of stress."----------------------------------------------
#plot distribution stress
par(mfrow = c(1,2))
hist(rstress.m3, main = "stress random data")
hist(perm.m3$stressvec, main = "stress permuted data")

#'
#'
## ----include=TRUE, echo=TRUE, message=FALSE, warning=FALSE-------------------------------------------------------------------------------------------------
#compute MDS biplot: run multivariate linear regression of
#external variables on configuration
zrect_constr <- scale(rect_constr, center = TRUE)
#biplot MDS: regression of standardized external variables on configuration dimensions
bimorse <- biplotmds(m3, zrect_constr, scale = TRUE)
#coefficients
bimorse_coef <- round(coef(bimorse), 3)
#print R-square of regressions
bimorse_r2 <- round(bimorse$R2vec, 3)
#print correlations of dimensions and external variables
bimorse_corr <- round(cor(m3$conf, zrect_constr), 3)

#'
## ----stability, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "Jackknife plot and biplot.", fig.height = 8, fig.width=5----------------
par(mfrow = c(2,1))
plot(jack.car, xlim = c(-1.2, 1.2), ylim = c(-1, 1))
# project external variables in the MDS solution
plot(bimorse, main = "Biplot Vector Representation", vecscale = 0.8,
  xlim = c(-1.2, 1.2), ylim = c(-1, 2.5), vec.conf = list(col = "brown"),
  pch = '', cex = 0.5)

#'
#'
## ----table3, include=TRUE, echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------------------------------------------------
tab <- bimorse_coef %>%
  rbind(bimorse_corr) %>%
  rbind(bimorse_r2)
rownames(tab)[5] <- "Rsquared"
tab %>%
  kable(booktabs = TRUE,
        caption = "Ordinal MDS dimensions and their link with the external variables (width and height).") %>%
  group_rows(1,2, group_label = "Coefficients") %>%
  group_rows(3,4, group_label = "Correlations")

#'
#' The results show that stress-1 is lowest using ordinal MDS (`r round(m3$stress,3)`, Table \@ref(tab:MDSresults)) or spline MDS (`r round(m4$stress,3)`). A residual plot and shephard plot for these is shown in Figure \@ref(fig:shephard). Since stress-1 is between 0.05 and 0.1 for ordinal MDS, it is a fair fit. The configuration can be seen in Figure \@ref(fig:ordinalMDS) where we see that rectangles 10 and 14 are perceived to be very similar (as are the pairs 11 & 15 and 12 & 16). Two options to better assess goodness of fit, are stress norm and permutation test (results for all measurement levels in Table \@ref(tab:MDSresults)). For ordinal MDS, we can see that the stress-1 obtained on the real dataset (`r round(m3$stress,3)`) is much lower than `r round(distr_stress3,3)` (the mean $-2*$ standard deviation of the stress random data). Similarly, it is much lower than the mean of the distribution from the stress permuted data (`r round(mean(perm.m3$stressvec),3)`, Figure \@ref(fig:distrstress)). We can thus conclude that ordinal MDS has a satisfactory fit.
#'
#' From the jackknife plot for the chosen, ordinal MDS (Figure \@ref(fig:stability)), we can see that the solution is quite stable, because the positions of the points do not change much. The stability measure reported by jackmds is `r round(jack.car$stab, 3)`, also confirming that the solution is very stable.
#'
#' When projecting the external variables height and width in the MDS solution (Table \@ref(tab:table3)), we can see from the coefficients and the plot (Figure \@ref(fig:stability)) that, if D1 increases one unit, the predicted average height of rectangles will decrease 1.623 standard deviations, and if D2 increases one unit, the predicted average height will increase 0.334 standard deviations. On the other hand,  if D2 increases one unit, the predicted average width of rectangles will increase 2.778 standard deviations, and if D1 increases one unit, the predicted average width will increase 0.096 standard deviations. We also see that D1 and D2 account for 98.7% of the variation in the height around its mean, and 94.5% of the variation in the width around its mean. The correlations of D1 and D2 with the external variables show that the first dimension has a strong negative correlation with the height (cor = -0.987), and that the second dimension has rather strong positive correlation with the width (0.970).
#'

