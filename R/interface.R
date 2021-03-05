#' @export
initialize_comhub <- function(network_name, expression_data, transcription_factors){

  transcription_factors <- transcription_factors[transcription_factors %in% rownames(expression_data)]

  transcription_factors <- as.data.frame(transcription_factors)

  c_core <- comhub_main(network_name = network_name,
              expression_data = as.data.frame(expression_data),
              transcription_factors = transcription_factors)

  construct_comhub_object(c_core = c_core,
                          expression_data <- as.data.frame(expression_data),
                          network_name = network_name,
                          transcription_factors = transcription_factors)

}
#' @export
MODifieR_input_to_comhub <- function(MODifieR_input,
                                     network_name,
                                     transcription_factors = NULL){
  if (is.null(transcription_factors)){
    transcription_factors <- ComhubbeR::transcription_factors
  }

  initialize_comhub(network_name = network_name,
                    expression_data = MODifieR_input$annotated_exprs_matrix,
                    transcription_factors = transcription_factors)
}

#' @export
summary.comhub <- function(c){
  cat("Network name is", c$network_name, "\n")
}


construct_comhub_object <- function(c_core, expression_data, network_name, transcription_factors){
  c <- list()
  c$python_object <- c_core

  c$expression_data <- expression_data

  c$network_name <- network_name

  c$transcription_factors <- transcription_factors

  class(c) <- "comhub"

  return(c)
}
