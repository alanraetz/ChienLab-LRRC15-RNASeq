#Dennis Minn

library(EnhancedVolcano)
library(ggrepel)
#Creates Volcano with LABELLED gene symbols 
#n <- max number of genes displayed
#res <- results
mk_volcano <- function(res, n, p, selectLab, fc, xlim, ylim, unit, title, xlabel, ylabel){
  
  if(missing(fc)){fc <- 1}
  if(missing(unit)){unit <- ""; print("UNITS ARE MISSING")}
  if(missing(title)){title <- ""}
  
  if(!missing(n) && !missing(selectLab)){
    return("INVALID ARGUMENT, please provide EITHER n to select the lowest p-adj genes OR selectLab to select specific geneIDs")
  }
  
  #add row names as a column (dplyr removes row.names)
  res$name <- map_symbols_adv(rownames(res))
  
  if(!missing(n)){
    #selects labels and makes corresponding subtitle
    labels <- head(res[order(res$padj),],n)
    labels <- dplyr::filter(as.data.frame(labels), padj <= p)
  }

  if(!missing(selectLab)){
    labels <- which(rownames(res) %in% selectLabs)
    labels <- df[labels,]
  }
  
  if(nrow(labels) == n){
    subtitle <- paste0("Volcano Plot with lowest ", n," p-adj genes labelled \n p-adj <= ",p," & |log2(foldchange)| >= ",fc)
  } else {
    subtitle <- paste0("Volcano Plot with significant p-adj genes labelled \n p-adj <= ",p," & |log2(foldchange)| >= ",fc)
  }
  
  if(missing(xlim)){
    xmin <- res$log2FoldChange[which.min(res$log2FoldChange)]
    xmax <- res$log2FoldChange[which.max(res$log2FoldChange)]
    xlim <- c(xmin, xmax)
  }
  
  if(missing(ylim)){
    ymax <- res$padj[which.min(res$padj)]
    ylim <- c(0, -log10(ymax))
  }
  #creates x axis label
  x_lab <- paste("log2(foldchange) ", unit)
  
  if(!missing(xlabel)){
      x <- xlabel
  }
  if(!missing(ylabel)){
      y <- ylabel
  }

  EnhancedVolcano(
    res,
    lab = res$name,                         #labels all genes
    selectLab = labels$name,                #specifies labels w/ desired genes
    xlab = (xlabel),                         #x-axis label
    ylab = ("-log10(p-adj)"),               #y-axis label
    x = 'log2FoldChange',                   #x position
    y = 'padj',                             #y position
    title = title,                          #sets title
    titleLabSize = 10,                      #title font size
    subtitle = subtitle,                    #sets description
    subtitleLabSize = 8,                    #subtitle font size
    #xlim = xlim,                            #boundaries of x-axis (log2FoldChange)
    ylim = ylim,                            #boundaries of y-axis (p-adj)
    pCutoff = p,                            #sets p-adj threshold
    FCcutoff = fc,                          #sets lfc threshold 
    legendVisible = FALSE,                  #hides legend (should be described in subtitle)
    boxedlabels = TRUE,                     #boxes labels
    drawConnectors = TRUE,                  #draws arrows
    widthConnectors = .2,                   #connector line width 
    lengthConnectors = unit(0.005, 'npc'),  #eliminates arrow
    transcriptLabSize = 2,                  #font size
    col = c('grey','grey','grey','red3'))   #highlights genes above both threshold
}
