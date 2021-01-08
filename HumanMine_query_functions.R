#Brandon Blobner 
# Library InterMineR, the HumanMine API
library(InterMineR)


query.Disease<-function(gene.symbol) {
# Set databse to HumanMine
im <- initInterMine(mine=listMines()["HumanMine"])

# Set template to query Gene_Loation
queryDisease = getTemplateQuery(
  im = im, 
  name = "Gene_Alleles_Disease2"
)

# Set query constraints
queryDisease$where[[1]][["value"]] <- paste(gene.symbol)

# Query gene paths
runQuery(im, queryDisease)
}

## Query expression

query.Expression<-function(gene.symbol) {
im <- initInterMine(mine=listMines()["HumanMine"])

# Set template to query Gene_Loation
queryExpression = getTemplateQuery(
  im = im, 
  name = "Gene_ExpressionProteinAtlas"
)

# Set query constraints
queryExpression$where[[2]][["value"]] <- paste(gene.symbol)

# Query gene paths
runQuery(im, queryExpression)
}


## Set template to query Gene snps
query.Snps<-function(gene.symbol) {
  
im <- initInterMine(mine=listMines()["HumanMine"])

querySnps = getTemplateQuery(
  im = im, 
  name = "Gene_SigSNP"
)

# Set query constraints
querySnps$where[[1]][["value"]] <- paste(gene.symbol)

# Query gene paths
runQuery(im, querySnps)
}


query.Publication<-function(gene.symbol) {
im <- initInterMine(mine=listMines()["HumanMine"])

# Set template to query Gene_Loation
queryPublication = getTemplateQuery(
  im = im, 
  name = "Gene_To_Publications"
)

# Set query constraints
queryPublication$where[[1]][["value"]] <- paste(gene.symbol)

# Query gene paths
runQuery(im, queryPublication)

}

