#Dennis Minn

#Requires 

#goal: Subsets countdata based on specific significance OR desired genes  
#param rld: RangedSummarizedExperiment object / normalized count data
#           created by DESeqDataSet object and vst() / (preferably) rlog()
#           can use unnormalized data (not recommended) 
#           DESeq2 creates ddsTxi object
#
#(p and res)/gene_list: parameter for selecting genes (choose 1 option)
#param p: decimal representing the p-adj threshhold for significant genes
#param res: results table created from any DESeqDataSet object and DESeq2's results()
#           (optional parameter iff matrix is created using a genelist)
#
#param gene_list: list/column vector of desired genes   
#
#return: input matrix for heatmap
mk_matrix <- function(res, rld, p, gene_list){
  
  #Creates matrix based on input gene_list 
  if(!missing(gene_list)){
    labels <- which(rownames(rld) %in% gene_list)   #selects genes with matching name
    mat <- assay(rld)[labels, ]                     #obtains data
    mat <- mat - rowMeans(mat)                      #better highlight expression of genes
    rownames(mat) <- map_symbols_adv(rownames(mat)) #maps symbols 
    return(mat)
  }
  
  #Creates matrix based on p-adj values in results table 
  if(!missing(res) && !missing(p)){
    
    #reorders results table
    temp <- res[order(res$padj),] 
    n <- 1 
    
    #Counts number of genes below p-adj threshhold
    while(temp[n,5] <= p & n <= 100){
      n <- n+1
    }
    topSigGenes <- head(order(res$padj), n-1)       #Selects genes under p-adj threshhold, 
                                                    #only lowest 100 p-adj genes if the amount of genes exceeds 100
    mat <- assay(rld)[topSigGenes, ]                #creates gene matrix based on position
    mat <- assay(rld);
    mat <- mat - rowMeans(mat)                      #better highlight expression of genes
    rownames(mat) <- map_symbols_adv(rownames(mat)) #maps symbols
    return(mat)
  } 
  print("invalid parameters")
}

#goal: Creates heatmap
#param mat: matrix created by mk_mat() 
#param df: dataframe object that holds characteristic of samples (look at pheatmap API for specifics and how to initialize)
#param sample_names: list/column vector of sample names
#param title <- string that represents the title of heatmap
# 
#recommend including the following in title:
#specific characteristic compared 
#number of genes in heatmap
#(if used in input matrix) p-adj threshold 
mk_pheatmap <- function(mat, df, sample_names, title){
  
  #creates empty string for title to prevent errors
  if(missing(title)){title <- ""}
  
  pheatmap(
    mat,                        #Matrix of count data 
    annotation_color=df,        #Creates legend
    main = title,               #Creates title
    labels_col = sample_names,  #Sample Names
    cellheight = 4,             #Height of displayed data
    cellwidth = 8,              #Width of displayed data
    fontsize = 4)               #Font size of genes
}
#DISCLAIMER:
#Please adjust pdf according graph size
#100 genes will not fit default height pdf 
