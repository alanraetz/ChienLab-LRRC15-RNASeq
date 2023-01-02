#Maps genes Ids to symbols
map_symbols <- function(geneIDs) {
  geneIDs <- gsub("\\..*", "", geneIDs)
  symbols <- mapIds(
    org.Hs.eg.db,
    keys=geneIDs,
    column="SYMBOL",
    keytype="ENSEMBL",
    multiVals="first")
  return(symbols)
}

#Map non-empty symbols to gene Ids
map_symbols_adv <- function(geneIDs) {
  temp <- gsub("\\..*", "", geneIDs)
  symbols <- mapIds(
    org.Hs.eg.db,
    keys=temp,
    column="SYMBOL",
    keytype="ENSEMBL",
    multiVals="first")
  
  for(i in 1:length(geneIDs)){
    if(is.na(symbols[i])) {
      symbols[i] <- geneIDs[i]
    } 
  }
  
  return(symbols)
}

#Map Ensembl IDs
map_ids <- function(dds, symbols) {
  ids <- mapIds(
    org.Hs.eg.db,
    keys=symbols,
    column="ENSEMBL",
    keytype="SYMBOL",
    multiVals="first")
  
  pos <- gsub("\\..*", "", rownames(dds))
  ids <- which(pos %in% ids)
  return(rownames(dds[ids,]))
}

#Selects the n lowest p-adj values 
get_genes <- function(res, n){
  topSigGenes <- head(order(res$padj), n)
  topSigGenes <- rownames(res[topSigGenes,])
  return(topSigGenes)
}

#Replace empty values
replace_NA <- function(intgroup, val){
  for(i in 1:length(intgroup)){
    if(is.na(intgroup[i])) {intgroup[i] <- val}
  }
  return(intgroup)
}

