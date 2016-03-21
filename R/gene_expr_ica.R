#' Custom ICA function for analyzing gene expression data.
#'
#' Performing ICA on a dataset and create a list object with results.
#'
#' @param phenotype.mx Phenotype matrix with diemnsions g x N
#' @param info.df Dataframe that holds sample covariates (ex. population, gender, age, etc...)
#' @param check.covars Column names of info.df which hold the covariates
#' that should be used for association testing with IC coefficients.
#' @param k.est Number of components to be estimated or method to estimate it.
#' @param scale.pheno Logical value specifying the scaling of row of the phenotype.mx.
#' @param h.clust.cutoff is the cutoff value used in hierarchical clustering. Default is set to 0.3.
#' @param n.runs Number of runs for estimating k. Default value is set to 5.
#' @param max.iter Maximum iterations for estimating k for each run. Default value is set to 10.
#' @param n.cores Number of cores to be used for estimating k. Default is set to 1.
#' @param cor.threshold Threshold for significant correlation calling. Default is set to 0.05.
#' @param similarity.measure How to measure the similarity between ICs.
#' If set to "peaks" only gene weights that are greater than 1 sd are used to calculate similarity.
#' @return List with the following entries.
#' @keywords keywords
#'
#' @import mclust
#' @export
gene_expr_ica <- function(phenotype.mx = NULL, 
                          info.df = NULL, 
                          check.covars = NULL,
                          k.est = NULL, 
                          scale.pheno = FALSE, 
                          h.clust.cutoff = 0.3,
                          n.runs = 5, max.iter = 10, 
                          n.cores = NULL, cor.threshold = 0.05,
                          similarity.measure = "peaks", var.cutoff = 99){

    if(is.null(phenotype.mx)){
        stop("Error: Phenotype matrix is missing \n")

    }

    if(is.null(colnames(phenotype.mx))){
      message("<phenotype.mx> is missing column names, set to default. \n")
      colnames(phenotype.mx) <- paste("sample",c(1:ncol(phenotype.mx)), sep = "_")
    }

    if(is.null(rownames(phenotype.mx))){
      message("<Phenotype.mx> is missing row names, set to default. \n")
      rownames(phenotype.mx) <- paste("feature",c(1:nrow(phenotype.mx)), sep = "_")
    }

    if(is.null(n.cores)){
        n.cores = 1
    }

    pheno.nrow <- nrow(phenotype.mx)
    pheno.ncol <- ncol(phenotype.mx)
    message("Original dimensions of <phenotype.mx> = ", pheno.nrow , " x ", pheno.ncol, "\n")

    if (pheno.nrow < pheno.ncol){
      message("[Caution] Number of samples exceeding number of measured features,
              please check rows and columns of <phenotype.mx> \n")
      message("- If you are from the future and have more samples than measured features,
              disregard the above message and please proceed. \n")
    }

    if(!is.null(info.df)){

      message("Checking dimensions and column/row names \n")

      if(nrow(info.df) == ncol(phenotype.mx)){
        message("[Good] Number of samples in <info.df> and <phenotype.mx> match \n")
      } else {
        stop("Error: Sample numbers in <info.df> and <phenotype.mx> don't match. Stopping script")
      }

      matching.names <- sum(rownames(info.df) %in% colnames(phenotype.mx))

      if(matching.names == ncol(phenotype.mx)){
        message("[Good] All samples in <phenotype.mx> are accounted for in <info.df> \n")
      } else {
        stop("Missing sample Information: Check rownames of <info.df> and column names of <phenotype.mx>")
      }
    }

    if(!is.null(info.df) & is.null(check.covars)){
      message("- Sample info detected but missing input for <check.covars> \n")
      message("- Using column names of <info.df> as <check.covars> \n")
      check.covars <- colnames(info.df)
    }

    message("------ Pre-processing Data ------- \n")
    # removing 0 variance genes and scaling and centering the phenotype matirx
    phenotype.mx <- pre_process_data(phenotype.mx, scale.pheno = scale.pheno)


    if(is.null(k.est)){
      message("Missing input for <k.est>, using the number of principal components explaining 99% of total variance \n")
      pca.pheno <- prcomp(t(phenotype.mx))
      percent <- (cumsum(pca.pheno$sdev^2) /sum(pca.pheno$sdev^2)) * 100
      k.est <- which(percent > var.cutoff)[1]
      message(k.est," components needed to explain more than ",var.cutoff,"% of the variance \n")

      if(k.est == 1){
        stop("1 component explains more than",var.cutoff,"% of the variance,
        check your data or set <k.est> to a number bigger than 1 \n")
      }


    }

    ica.list <- list()
    ica.result <- list()


    message("------ Running ICA ------- \n")
    message("* This may take some time depending on the size of your dataset \n")
    if(n.runs >1){
      message("- Generating ",n.runs," replicates using ",n.cores," core(s): ",k.est," Components are estimated in each run\n")

      ica.list <- parallel::mclapply(1:n.runs,function(x) fastICA_gene_expr(phenotype.mx, k.est,
                                                               fun = "logcosh",                            # function that should be used to estimate ICs, default is logcosh
                                                               alpha = 1, scale.pheno = FALSE,                  # row.norm is set to false since the phenotype.mx is scaled separately
                                                               maxit=500, tol = 0.0001, verbose = FALSE), mc.cores = n.cores)

      for(i in 1:length(ica.list)){

        ica.list[[i]]$peak.mx <- apply(ica.list[[i]]$S, 2, function(x) 1*(abs(x) > 2*sd(x)))

      }

      # After running ICA sevaral times
      #combined.A <- do.call(rbind, lapply(ica.list, function(x) x$A))
      combined.S <- do.call(cbind, lapply(ica.list, function(x) x$S)) # combine all components into a single matrix
      peak.matrix <- do.call(cbind, lapply(ica.list, function(x) x$peak.mx)) # combine all peak position matrices as well

      if(similarity.measure == "peaks"){
        if(0 %in% apply(peak.matrix, 2, sum)){
          cor.mx <- cor(combined.S)
        } else{
          peak.component <- peak.matrix * combined.S             # do element-wise multiplication to only save values for peaks
          cor.mx <- cor(peak.component) # calculate correlation between components (only with their peak values)
        }
      } else {
        cor.mx <- cor(combined.S)
      }


      # clustering of the components
      dissimilarity <- 1 - abs(cor.mx) # create dissimilarity matrix
      cor.dist <- as.dist(dissimilarity) # convert into distance matrix format for clustering
      h.clust <- hclust(cor.dist)          # run hierarchical clustering
      groups <- cutree(h.clust, h=h.clust.cutoff)   # cut tree at height 0.3 (absolute correlation > 0.7)


      group.table <- table(groups)     # count member components for each group
      multi.component.group <- which(group.table >= (n.runs * 0.9)) # get groups with more than 2 members

      k.update <- length(multi.component.group)
      message("- ",k.update," Replicating Components Estimated \n")
      if(k.update < 1){
        stop('None of the ICs replicated. You could try to the increase sample size or <n.runs>. \n')
      }
      Avg.S <- matrix(0,nrow = dim(combined.S)[1],ncol = k.update)
      #Avg.A <- matrix(0,nrow = k.update, ncol = dim(combined.A)[2])
      # for each group calculate the average component
      for(i in 1:length(multi.component.group)){
          group.members <- which(groups %in% multi.component.group[i]) # get component indexes for groups with multiple components
          sub.group.S <- combined.S[,group.members]                # subset matrix for those components
          #sub.group.A <- combined.A[group.members,]
          group.cor.mx <- cor(sub.group.S)                       # calculate correlation between them
          match.sign <- ifelse(group.cor.mx[,1] < 0, -1,1 )  # in order to average the components signs need to be matched (positive vs negative)
          avg.component <- (sub.group.S %*% match.sign) / length(match.sign) # calculate the mean component
          #avg.mixing <- (match.sign %*% sub.group.A ) / length(match.sign)

          Avg.S[,i] <- avg.component
          #Avg.A[i,] <- avg.mixing
      }

      hclust.list <- list()
      hclust.list$plot <- h.clust
      hclust.list$cutoff <- h.clust.cutoff
      hclust.list$dist.mx <- dissimilarity
      ica.result$hclust <- hclust.list
      ica.result$S <- Avg.S
      ica.result$A <- solve(t(Avg.S) %*% Avg.S) %*% t(Avg.S) %*% phenotype.mx

      rm(Avg.S)

    } else if (n.runs ==1){
      message("Running ICA once with ", k.est," components to be estimated")

      ica.result <- fastICA_gene_expr(phenotype.mx, k.est,
                                      fun = "logcosh",                            # function that should be used to estimate ICs, default is logcosh
                                      alpha = 1, scale.pheno = FALSE,                  # row.norm is set to false since the phenotype.mx is scaled separately
                                      maxit=500, tol = 0.0001, verbose = FALSE)
      k.update <- k.est
    }

    rownames(ica.result$S) <- rownames(phenotype.mx)                   # Setting appropriate names for signals and mixing matrix
    colnames(ica.result$S) <- paste("IC",c(1:dim(ica.result$S)[2]),sep="")
    colnames(ica.result$A) <- colnames(phenotype.mx)
    rownames(ica.result$A) <- paste("IC",c(1:dim(ica.result$A)[1]),sep="")

    # Attaching the sample info dataframe to the ica list
    if(!is.null(info.df)){
      ica.result$info.df <- info.df[colnames(phenotype.mx),]
    }
    message("------ Post-ICA Processing ------- \n")
    message("- Labeling peak genes in each IC \n")
    ica.result$peaks <- apply(ica.result$S, 2, peak_detection)
    # peaks are defined as gene contributions that are larger than 2 standard deviations
    ica.result$peak.mx <- apply(ica.result$S, 2, function(x) 1*(abs(x) > 2*sd(x)))

    message("- Calculating variance explained by each IC \n")
    # get the total variance by the sums of squares of the scaled phenotype.mx
    total.var <- sum(phenotype.mx^2)

    # applying IC component-wise variance calculations
    var.IC <- sapply(1:dim(ica.result$A)[1],
                     function (x) IC_variance_calc(ica.result$S[,x], ica.result$A[x,]))

    # % variance explained by each IC
    percent.var <- (var.IC / total.var) * 100

    message("- Sanity Check : Total % of variance explained by ",k.update," ICs = ", sum(percent.var), "\n")

    message("- Creating index based on Variance explained \n")
    ica.result$order <- order(percent.var,decreasing = T) # ordering the ICs based on the amount of variance they explain
    ica.result$percent.var <- percent.var

    # Checking correlation between IC coefficients and measured covariates
    if(!is.null(check.covars) & !is.null(info.df)){
        message("- Checking associations between ICs and covariates \n")
        # Anova analysis for covariates vs ICA weights (A matrix)
        ica.result$cov.pval.mx <- component_association_test(ica.result$A,info.df,check.covars)
        # Multiple Hypothesis correction
        corrected.threshold <- cor.threshold / (dim(ica.result$cov.pval.mx)[1] * dim(ica.result$cov.pval.mx)[2])
        corr.idx <- which(ica.result$cov.pval.mx < corrected.threshold, arr.ind = T)
    } else{
        corr.idx <- NULL
    }

    if(length(corr.idx) != 0 ){
        correlated.ic <- unique(corr.idx[,1])

        covariate.corr.df <- data.frame(matrix(nrow = length(correlated.ic),ncol = 3))
        colnames(covariate.corr.df) <- c("IC","Covariate.idx","Covariate.Name")
        for( c in 1:length(correlated.ic)){
            ic.index <- correlated.ic[c]
            covar.index <- which.min(ica.result$cov.pval.mx[ic.index,])
            covariate.corr.df[c,"IC"] <- ic.index
            covariate.corr.df[c,"Covariate.idx"] <- covar.index
            covariate.corr.df[c,"Covariate.Name"] <- names(covar.index)
        }
        ica.result$cov.corr.idx <- covariate.corr.df

        rm(covariate.corr.df, covar.index, ic.index, c)
    } else {
        ica.result$cov.corr.idx <- NULL     # in case there are no apparent covariates
    }



    if(is.null(ica.result$cov.corr.idx)){
        sig <- rep(0,k.update)
        correlated.ic <- NULL
    } else{
        sig <- rep(0,k.update)
        #correlated.ic <- unique(ica.result$cov.corr.idx$IC)
        sig[correlated.ic] <- 1
    }

    # Turn off warnings to prevent unneccessary panic
    options(warn=-1)
    mclust.result <- apply(ica.result$A, 1, function(x) mclust::Mclust(x))
    # Turn warnings back on to prevent disasters
    options(warn=0)

    if(!is.null(info.df)){
      ica.result$ica.stat.df <- data.frame("N.peaks"=sapply(ica.result$peaks, function(x) length(x)), # Number of peaks for each IC
                                           "n.clust"= sapply(mclust.result, function(x) x$G),      # Number of predicted clusters
                                           "percent.var" = percent.var,                                     # Percent variance explained
                                           "corr.ic" = factor(sig), "idx" = c(1:k.update))             # if correlated with covariate = 1 , 0 otherwise

    } else {
      ica.result$ica.stat.df <- data.frame("N.peaks"=sapply(ica.result$peaks, function(x) length(x)), # Number of peaks for each IC
                                           "n.clust"= sapply(mclust.result, function(x) x$G),      # Number of predicted clusters
                                           "percent.var" = percent.var,                                     # Percent variance explained
                                           "idx" = c(1:k.update))             # if correlated with covariate = 1 , 0 otherwise

    }


    # which IC has more than 1 predicted clusters?
#    multi.clust <- which(ica.result$ica.stat.df$n.clust > 1)

    # get union between multi clusters and correlated ICs
#    hf.vec <- sort(union(multi.clust,correlated.ic))

#    hf.vec.names <-paste("IC",hf.vec,sep="")
    # name the ICs as IC#
#    n.hf <- length(hf.vec.names)
#    message("- ",n.hf," out of ",k.update," ICs marked as confounding factors = \n")
#    cat(hf.vec.names,"\n\n")

    # Creating a matrix indicating which genes are influenced by a given IC
    ica.result$ica.confeti.mx <- matrix(0,nrow = nrow(ica.result$S), ncol = ncol(ica.result$S))
    # rownames = gene names
    rownames(ica.result$ica.confeti.mx) <- rownames(ica.result$S)
    # column names = ICs
    colnames(ica.result$ica.confeti.mx) <- colnames(ica.result$S)

    #message("Creating CONFETI matrix for regression \n")
    for ( i in 1:ncol(ica.result$S)){
#        k <- hf.vec[i]
#        ic.name <- hf.vec.names[i]
        peak.temp <- names(ica.result$peaks[[i]])
        ica.result$ica.confeti.mx[peak.temp,i] <- 1
    }
    message("------ Process completed without any interuptions ------- :) \n")
    attr(ica.result, 'method') <- "ica"
    return(ica.result)
}
