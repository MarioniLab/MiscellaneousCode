library(matrixStats)
library(edgeR)
library(scater)
library(scran)
library(biomaRt)
register(SerialParam())

# GINI calculation function

giniCoeff <- function(x, na.rm = FALSE){
  if (!is.numeric(x)){
    warning("'x' is not numeric; returning NA")
    return(NA)
  }
  if (!na.rm && any(na.ind = is.na(x)))
    stop("'x' contain NAs")
  if (na.rm)
    x = x[!na.ind]
  if (all(x == 0)) {
    gini <- 0
  } else {
  N <- length(x)
  gini <- sum((2*seq_len(N) - N  -1)*sort(x))/(N-1)/sum(x)   # Problem with the function: For vector c(-25, 25) it would return NaN  --> Pseudocount?
  }
  return(gini)
}

# Function for technical dispersion

giniNull <- function(x, design=NULL, niterations=10000, npts=30){
  if (is.null(design)){
    design <- matrix(1, nrow=ncol(x), ncol=1)
  }
  sce_spike <- x[isSpike(x), ]
  conv_spike <- convertTo(sce_spike, type = "edgeR", get.spikes=TRUE)
  disp_spike <- estimateDisp(conv_spike, design=design)
  means_spike <- log(rowMeans(counts(sce_spike))+1)
  trendDisp <- approxfun(means_spike, disp_spike$trended.dispersion, rule=2)

  popsize <- ncol(x)
  all.means <- log(rowMeans(counts(x)+1))
  index.means <- seq(min(all.means), max(all.means), length.out=npts)
  all.gini <- as.list(numeric(length(index.means)))
  for (index in seq_len(npts)){
    disp = trendDisp(index.means[index])
    m <- matrix(rnbinom(niterations*popsize, mu=exp(index.means[index])-1, size=1/disp), ncol=popsize)  
    current_gini <- as.list(apply(m, 1, function(x) giniCoeff(x)))
    all.gini[[index]] <- current_gini
  }
return(list(all.gini, index.means))    #return a list with the first element all.gini and the second element the means of each index
}

# Calculate Gini Indizes with p-vals and FDR
giniCal <- function(x, nulloutput, na.rm=FALSE){ 
  if (class(x) == "SCESet") {
    norm_count <- cpm.default(counts(x), lib.size=sizeFactors(x)*1e6)
  } else {stop(paste0("unable to find an inherited method for function 'giniCal' for signature '", class(x), "'"))
  }
  gini_indexes <- apply(norm_count, 1, FUN=giniCoeff, na.rm=na.rm)
  p_value_sce <- numeric(nrow(x))
  all.means <- log(rowMeans(counts(x)+1))

  for (gene in seq_len(nrow(x))){
    gene.mean <- all.means[gene]
    gene_gini <- gini_indexes[gene]
  
    # Set higher and lower indixes
    index_p_values <- numeric(2)
    
    lower_index <- max(which(nulloutput[[2]] <= gene.mean))
    null_ginis <- nulloutput[[1]][[lower_index]]
    index_p_values[1] <- (sum(gene_gini <= null_ginis)+1) / (length(null_ginis)+1)
    
    higher_index <- min(which(nulloutput[[2]] >= gene.mean))
    null_ginis <- nulloutput[[1]][[higher_index]]
    index_p_values[2] <- (sum(gene_gini <= null_ginis)+1) / (length(null_ginis)+1)

    current.means <- c(nulloutput[[2]][lower_index], nulloutput[[2]][higher_index])

    
    if (current.means[1] == current.means[2]) {
     p_value_sce[gene] <- index_p_values[1]
    } else {
     p_value_sce[gene] <- approx(current.means, index_p_values, gene.mean)$y
    }
  } 
return(data.frame(mean = all.means, GINI = gini_indexes, p.value = p_value_sce, 
                  FDR = p.adjust(p_value_sce, method = "fdr"), row.names = rownames(fData(x))))
}

