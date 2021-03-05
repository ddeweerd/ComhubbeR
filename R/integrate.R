#' @export
retrieve_tfs_and_genes  <- function(MODifieR_module, MODifieR_input, comhub_object, ppi_network, ORcutoff = 2){

  module_genes <- MODifieR_module$module_genes

  if (class(MODifieR_module)[3] %in% c("DiffCoEx", "WGCNA", "MODA")){
    background_genes <- get_input_background(MODifieR_input)
  }else{
    background_genes <- get_module_background(ppi_network)
  }

  network = comhub_object$consensus_network %>%
    .[1:comhub_object$edge_cutoff,]

  or_pval_tf <- try(get_tf_OR(module_genes = module_genes,
                              network = network,
                              nbgenes = length(background_genes)))

  significant_tfs <- try(get_significant_tfs(or_pval_tf = or_pval_tf,
                                             comhub_hub_scores = comhub_object$community,
                                             tfs = unique(network$TF),
                                             ORcutoff = ORcutoff))
  if (class(significant_tfs) == "try-error"){
    return (significant_tfs)
  }

  tfs_and_genes <- try(append_genes_to_tfs(significant_tfs = significant_tfs,
                                           network = network,
                                           module_genes = module_genes))
  if(length(tfs_and_genes) == 0){
    stop("No significant hubs detected", call. = F)
  }
  tfs_and_genes
}

get_tf_OR <- function(module_genes, network, nbgenes) {
  tfs <- unique(network$TF)
  pval <- c()
  or <- c()
  for (tf in tfs) {

    TFtargetgenes <- network[(network$TF == tf),]$target
    TFmodule_genes <- intersect(TFtargetgenes, module_genes)

    #Contingency table
    inTFinModule <- length(TFmodule_genes)
    notTFnotModule <- nbgenes - length(union(module_genes, TFtargetgenes))
    inTFnotModule <- length(TFtargetgenes) - length(TFmodule_genes)
    notTFinModule <- length(module_genes) - length(TFmodule_genes)

    dat <- matrix(c(inTFinModule, inTFnotModule, notTFinModule, notTFnotModule), ncol=2, nrow=2)
    colnames(dat) <- c('inModule', 'notModule')
    rownames(dat) <- c('inTF', 'notTF')
    dat <- as.table(dat)

    # Fisher's exact test
    test <- fisher.test(dat)
    pval <- c(pval, test$p.value)
    or <- c(or, test[["estimate"]][["odds ratio"]])
  }

  or_pval_tf <- cbind(or, pval)
  rownames(or_pval_tf) <- tfs
  colnames(or_pval_tf) <- c('OR', 'pvalue')
  return(or_pval_tf)
}

get_significant_tfs <- function(or_pval_tf, comhub_hub_scores, tfs, ORcutoff) {
  or_pval_tf <- na.omit(cbind(or_pval_tf, comhub_hub_scores[tfs]))
  colnames(or_pval_tf) <- c('OR', 'pvalue', 'comhub hub score')
  or_pval_tf <- or_pval_tf[order(or_pval_tf[,'comhub hub score']),]

  hubscore <- or_pval_tf[,'comhub hub score']
  #Make into parameter?
  hub_cutoff <- .9
  hubcutoff <- min(hubscore[cumsum(hubscore) >= max(cumsum(hubscore)) * hub_cutoff]) #Hubs standing for 10% of interactions.

  hubs <- or_pval_tf[hubscore >= hubcutoff,]
  if (is.null(dim(hubs)) || nrow(hubs) == 0){
    stop("No significant hubs detected", call. = F)
  }
  hubs <- hubs[hubs[ ,'pvalue'] <0.05, ]
  if (is.null(dim(hubs)) || nrow(hubs) == 0){
    stop("No significant hubs detected", call. = F)
  }
  significanthubs <- hubs[hubs[,'OR'] >= ORcutoff, ]

  return(significanthubs)
}

append_genes_to_tfs <- function(significant_tfs, network, module_genes) {
  hubtargets <- list()
  allhubtargets <- c()
  for (hub in rownames(significant_tfs)) {
    hubtargetgenes <- network[(network$TF == hub), ]$target
    hubmodule_genes <- intersect(hubtargetgenes, module_genes)
    hubtargets[[hub]] <- hubmodule_genes

    allhubtargets <- append(allhubtargets, hubmodule_genes)
  }
  hubtargets[['all']] <- unique(allhubtargets)

  hubtargets
}

sort_comhub_genes <- function(module_name, module_genes, con){
  input_name <- as.character(MODifieRDB::MODifieR_module_from_db(module_name, con = con)$settings$MODifieR_input)
  input_data <- MODifieRDB::MODifieR_input_from_db(input_name, con = con)$diff_genes
  subset_genes <- input_data[(input_data$gene %in% module_genes), ]

  genes <- subset_genes$pvalue
  names(genes) <- subset_genes$gene
  sort(genes, decreasing = T)
}

get_input_background <- function(MODifieR_input){
  rownames(MODifieR_input$annotated_exprs_matrix)
}

get_module_background <- function(ppi_network){
  unique(unlist(ppi_network[,1:2]))

}



