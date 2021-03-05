#' @export
aracne_comhub <- function(comhub_object, network_cutoff = 100000){
  mim <- minet::build.mim(dataset = transpose_expression_data(comhub_object),
                          estimator = "spearman")

  net <- minet::aracne(mim)

  dimnames(net) <- list(rownames(comhub_object$python_object$expression_data),
                        rownames(comhub_object$python_object$expression_data))

 comhub_object$results$aracne <- melt_result(net, comhub_object, network_cutoff)

 return(comhub_object)

}
#' @export
clr_comhub <- function(comhub_object, network_cutoff = 100000){
  mim <- minet::build.mim(dataset = transpose_expression_data(comhub_object),
                          estimator = "spearman")
  net <- minet::clr(mim)

  dimnames(net) <- list(rownames(comhub_object$python_object$expression_data),
                        rownames(comhub_object$python_object$expression_data))

  comhub_object$results$clr <- melt_result(net, comhub_object, network_cutoff)

  return(comhub_object)

}
#' @export
tigress_comhub <- function(comhub_object, network_cutoff = 100000){
  edgepred <- tigress::tigress(expdata = transpose_expression_data(comhub_object),
                               tflist = unlist(comhub_object$python_object$transcription_factors))[[5]]

  comhub_object$results$tigress <- melt_result(edgepred, comhub_object, network_cutoff)

  return(comhub_object)
}

#' @export
genie_comhub <- function(comhub_object,  network_cutoff = 100000){
  vim <- genie_main(t(comhub_object$python_object$expression_data),
                    rownames(comhub_object$python_object$expression_data),
                    comhub_object$python_object$transcription_factors[,1])

  genie_result <- genie_link_list(vim,
                                 rownames(comhub_object$python_object$expression_data),
                                 comhub_object$python_object$transcription_factors[,1])

  comhub_object$results$genie3 <- build_genie_result(genie_result, network_cutoff)

  return(comhub_object)
}

#' @export
pcc_comhub <- function(comhub_object, network_cutoff = 100000){
  comhub_object$results$pcc <- comhub_object$python_object$pcc()

  return(comhub_object)
}

#' @export
elasticnet_bootstrap_comhub <- function(comhub_object, bootstrap = 100, network_cutoff = 100000){
  comhub_object$results$elasticnet_bootstrap <- comhub_object$python_object$elasticnet_bootstrap(bootstrap = as.integer(bootstrap))

  return(comhub_object)
}

#' @export
pairwise_correlation <- function(comhub_object){
  comhub_object$edge_cutoff <- comhub_object$python_object$pairwise_correlation(network_files = unname(comhub_object$results),
                                                   names = names(comhub_object$results), plot = F)

  return(comhub_object)
}
#' @export
get_tf_outdegree <- function(comhub_object){
  comhub_object$tf_outdegree <- comhub_object$python_object$get_tf_outdegree(network_files = unname(comhub_object$results),
                                                                                names = names(comhub_object$results), edge_cutoff = comhub_object$edge_cutoff)

  return(comhub_object)
}

#' @export
community <- function(comhub_object){
  comhub_object$community <- apply(comhub_object$tf_outdegree, 1, mean)

  return(comhub_object)
}
