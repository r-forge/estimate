###
### filterCommonGenes.R
###


##-----------------------------------------------------------------------------
filterCommonGenes <- function(input.f, output.f, id=c("GeneSymbol","EntrezID")) {

    ## Check arguments
    stopifnot((is.character(input.f) && length(input.f) == 1 && nzchar(input.f)) ||
              (inherits(input.f, "connection") && isOpen(input.f, "r")))
    stopifnot((is.character(output.f) && length(output.f) == 1 && nzchar(output.f)))
    id <- match.arg(id)   
     
    ## Load estimate.RData
    data(common_genes)
   
    ## Read input data
    input.data <- read.table(input.f, header=TRUE, row.names=1, sep="\t", 
                             quote="", stringsAsFactors=FALSE)
     
    if(id == "GeneSymbol"){
        input.f1 <- merge(common_genes, input.data, by.x="GeneSymbol", by.y="row.names")
    } else {
        input.f1 <- merge(common_genes, input.data, by.x="EntrezID", by.y="row.names")
    }
    rownames(input.f1) <- input.f1$GeneSymbol
    input.f1 <- input.f1[,-1:-ncol(common_genes)]
    x.mismatched.genes <- nrow(common_genes) - nrow(input.f1)
    print(paste("This dataset includes ", nrow(input.f1), "genes.", sep=""))
    print(paste(x.mismatched.genes, "genes were mismatched.", sep=""))
    outputGCT(input.f1, output.f)
 }
