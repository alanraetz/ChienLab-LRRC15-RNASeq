#Dennis Minn

library(ggplot2)
library(ggrepel)
#Creates MAplot with labels on the "n" lowest p-adj genes 
#res <- results table
#n <- max number of genes displayed
#selectLab <- list of geneIDs to be displayed 
#p <- p-adj threshold (labelled red) 
#
#Optional Arguments
####################
#unit <- Comparison made (PeerGroup/ReferenceGroup)
#title <- title of graph
#subtitle <- subtitle of graph
#ylim <- limits height of the y-axis
#by <- specifies the "ticks" of your y-axis
mk_MAplot <- function(res, n, p, selectLabs, unit, title, subtitle, ylim, by){
  
    if(missing(unit)){unit <- ""; print("UNITS ARE MISSING")}
    if(missing(title)){title <- ""}
    if(!missing(n) && !missing(selectLabs)){
      return("INVALID ARGUMENT, please provide EITHER n to select the lowest p-adj genes OR selectLab to select specific geneIDs")
    }
      
    #creates symbols and significant columns
    res$symbol <- map_symbols_adv(rownames(res))
    df <- dplyr::mutate(as.data.frame(res), significant = padj <= p)
    
    #sets NA values to FALSE and 0 
    df$significant <- replace_NA(df$significant,  FALSE)
    df$log2FoldChange <- replace_NA(df$log2FoldChange, 0)
    
    #Creates dataframe for desired genes
    if(!missing(selectLabs)){
      labels <- which(rownames(res) %in% selectLabs)
      labels <- df[labels,]
    } 
    
    if(!missing(n)){
      labels <- head(df[order(res$padj),],n)
      labels <- dplyr::filter(as.data.frame(labels), significant)
    }
   
    #Creates subtitle
    if(missing(subtitle) ){
      if(!missing(n)){
        n <- nrow(labels)
        subtitle <- paste0("MAPlot with lowest ", n," p-adj genes labelled \n p-adj <= ", p)
      } else {
        subtitle <- paste0("p-adj <= ", p)
      }
    }
    
    y_lab <- paste("log2(foldchange) ", unit)
    
    #plots graph
    ggplot(df, aes(x=log2(baseMean), y=log2FoldChange, colour=significant)) +
      
      ggtitle(title) +
      
      #sets graphs size and spacing    
      {if(!missing(ylim)){
        if(!missing(by)){
          scale_y_continuous(breaks = seq(ylim[1],ylim[2],by), limits = c(ylim[1],ylim[2]))
        } else {
          scale_y_continuous(breaks = seq(ylim[1],ylim[2],1), limits = c(ylim[1],ylim[2]))
        } 
      }} +
            
      #circles and labels significant genes
      labs(y = y_lab, subtitle = subtitle, font="Times New Roman") +
      geom_point(size = .5, cex = .35) +
      geom_point(data=labels, aes(x=log2(baseMean), y=log2FoldChange),size = .65, cex = .35, shape=1, color="black") +
      ggrepel::geom_label_repel(data=labels, aes(label=symbol), color="black", size=2, segment.size=.1) +
      
      #color = significant_column, grey = FALSE, red = TRUE
      scale_color_manual(values = c("grey", "red"))  +
      
      #removes legend
      theme(legend.position = "none") 
    

}