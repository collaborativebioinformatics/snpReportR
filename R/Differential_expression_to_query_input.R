# Differential Expression to query input

# Read in differntial expression results
diff.ex<-read.table(here::here("data-raw/edgeR_test_results.tsv"), header = T)

# Filter significant results
diff.ex<-diff.ex[which(diff.ex$PValue<0.05),]

# Sort significant results by log fold change
diff.ex<-diff.ex[order(diff.ex$logFC, decreasing = T),]

# Extract top (greatest logFC) five upregulated and downregulated genes
top.diff.ex.up<-head(diff.ex, n=5)
top.diff.ex.down<-tail(diff.ex, n=5)
top.diff.ex<-rbind(top.diff.ex.up, top.diff.ex.down)

# Extract gene ID
gene.names<-top.diff.ex[,1]
