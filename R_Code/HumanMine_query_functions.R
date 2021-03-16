#Brandon Blobner 
# Library InterMineR, the HumanMine API
library(InterMineR)


#' Query HumanMine for associated diseases
#'
#' @param gene.symbol character(1) string, specifying gene symbol
#'
#' @return Dataframe of associated diseases
#' @export
#'
#' @examples query.Disease("SCNN1A")
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

#' Query HumanMine for tissue expression data
#'
#' @param gene.symbol character(1) string, specifying gene symbol
#'
#' @return Dataframe of tissue expression data
#' @export
#'
#' @examples query.expression("SCNN1A")
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
#' Query HumanMine for snps in genes of interest
#'
#' @param gene.symbol character(1) string, specifying gene symbol
#'
#' @return Dataframe of snp data for genes of interest
#' @export
#'
#' @examples query.Snps("SCNN1A")
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


#' Query HumanMine for publications relevant to genes of interest
#'
#' @param gene.symbol character(1) string, specifying gene symbol
#'
#' @return Dataframe of publications and relevant data
#' @export
#'
#' @examples query.Publications("SCNN1A")
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

