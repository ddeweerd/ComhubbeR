#' @import tidyverse
#' @import magrittr
#' @import tigress
#' @import minet
#' @import dplyr
#' @import reshape2
#' @import reticulate

transpose_expression_data <- function(comhub_object){
  comhub_object$python_object$expression_data %>%
    t(.) %>%
    as.data.frame(.)
}


melt_result <- function(net, comhub_object, network_cutoff){
  reshape2::melt(net) %>%
    dplyr::filter(., value > 0) %>%
    dplyr::filter(Var1 %in% comhub_object$python_object$transcription_factors[,1]) %>%
    dplyr::arrange(., desc(value)) %>%
    magrittr::set_colnames(., c('TF', 'target', 'confidence')) %>%
    rapply(., as.character, classes = "factor", how = "replace") %>%
    mutate_at(c(1,2), as.character) %>%
    dplyr::slice(., 1:network_cutoff)
}

build_genie_result <- function(genie_result, network_cutoff){
  gene_names <- genie_result[[2]]
  results <- genie_result[[1]]
  col1 <- unlist_genie(1, results, gene_names)
  col2 <- unlist_genie(2, results, gene_names)

  data.frame(col1, col2, unlist(sapply(results, function(x)x[3])), stringsAsFactors = F) %>%
    set_colnames(c('TF', 'target', 'confidence')) %>%
    dplyr::slice(., 1:network_cutoff)
}

unlist_genie <- function(col, results, gene_names){
  indici <- sapply(results, function(x)x[col]) %>%
    unlist(.) + 1

  gene_names[indici]
}

#' @export
consensus_network <- function(comhub_object){
  comhub_object$consensus_network <- Reduce(rbind, comhub_object$results) %>%
    dplyr::select(-confidence) %>%
    dplyr::count(TF, target) %>%
    dplyr::arrange(., dplyr::desc(n))

  return(comhub_object)
}
#' @export
combine_comhub_objects <- function(comhub_objects){
  c1 <- comhub_objects[[1]]
  c1$results <-  lapply(comhub_objects, function(x)x$results) %>%
    unlist(., recursive = F)

  return (c1)
}

#' @export
retrieve_community <- function(comhub_object){
  revive_python_object(comhub_object) %>%
   pairwise_correlation() %>%
   get_tf_outdegree() %>%
   community() %>%
   consensus_network()
}
#' @export
revive_python_object <- function(comhub_object){
  comhub_object$python_object <- comhub_main(comhub_object$network_name,
                                             comhub_object$expression_data,
                                             comhub_object$transcription_factors)

  return (comhub_object)
}
