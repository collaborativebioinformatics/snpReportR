#Brandon Blobner 
# Example: gene.names<-c("SCNN1A", "SCNN1B", "SCNN1G", "SCNN1D")

query_get_and_filter <- function(gene.names){
  print(gene.names)
  
  for(i in 1:length(gene.names)){
    # Disease query
    assign(paste(gene.names[i],"Disease", sep="."),query.Disease(gene.names[i]))
    
    # Associated diseases for report
    gene.disease<-get(paste(gene.names[i],"Disease", sep="."))
    assign(paste(gene.names[i],"Associated.Diseases", sep="."), 
           gene.disease[-which(duplicated(gene.disease$Gene.alleles.diseases.name)),which(colnames(gene.disease)=="Gene.alleles.diseases.name")])
    
    # Expression query
    assign(paste(gene.names[i],"Expression", sep="."),query.Expression(gene.names[i]))
    
    # Highly expressed tissues for report
    gene.expression<-get(paste(gene.names[i],"Expression", sep="."))
    assign(paste(gene.names[i],"Top.Expression", sep="."),
           gene.expression[head(order(gene.expression$Gene.rnaSeqResults.expressionScore, decreasing = T), n=5),
                           which(colnames(gene.expression)=="Gene.rnaSeqResults.tissue")])
    
    # SNPS query
    assign(paste(gene.names[i],"SNPS", sep="."),query.Snps(gene.names[i]))
    
    # Publication query
    assign(paste(gene.names[i],"Publication", sep="."),query.Publication(gene.names[i]))
    
    # Recent publications
    gene.publication<-get(paste(gene.names[i],"Publication", sep = "."))
    assign(paste(gene.names[i], "Recent.Publications", sep = "."),
           gene.publication[1:5,which(colnames(gene.publication)%in%c("Gene.publications.firstAuthor",
                                                                      "Gene.publications.year",
                                                                      "Gene.publications.pubMedId"))])
  }
  
  results_names <- ls(pattern = paste(gene.names, collapse = "|"))
  results <- lapply(results_names,  function(x) get(x))
  names(results) <- results_names
  return(results)
  
  # Associated diseases
  # assign(gene.disease, get(paste(gene.names[i], "Disease", sep=".")))
  # assign(paste(gene.names[i],"Associated.Diseases", sep="."), +
  # gene.disease[-which(duplicated(gene.disease$Gene.alleles.diseases.name)),which(colnames(gene.disease)=="Gene.alleles.diseases.name")])
}