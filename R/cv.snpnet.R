#' Cross-validation for snpnet
#' 
#' Does k-fold cross-validation for snpnet, and returns a value for lambda
#' 
#' @importFrom foreach '%dopar%'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#'
#' @export

cv.snpnet <- function(genotype.pfile, phenotype.file, phenotype, family = NULL, 
                      covariates = NULL, weights = NULL, alpha = 1, nlambda = 100, 
                      lambda.min.ratio = NULL, full.lams = NULL, foldids = NULL,
                      nfolds = 5, p.factor = NULL, ncores = 1, parallel = FALSE,
                      status.col = NULL, mem = NULL, configs = NULL) {
  
  `%dopar%` <- if (parallel) foreach::`%dopar%` else foreach::`%do%`
  
  if (!parallel | ncores == 1) {
    cores_per_fold <- ncores
  } else{
    cores_per_fold <- max(1, floor(ncores/nfolds))
    options(cores = cores_per_fold)
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    
    on.exit(parallel::stopCluster(cl))
  }
  
  ### Set up cross-validation folds in phenotype.file
  
  ids <- readIDsFromPsam(paste0(genotype.pfile, '.psam'))
  configs <- setupConfigs(configs, genotype.pfile, phenotype.file, phenotype, 
                          covariates, alpha, nlambda, split.col = NULL, 
                          p.factor, status.col, mem)
  
  cv.phenotype.file <- paste0(configs[['cv.full.prefix']], '.tsv')
  
  phe <- readPheMaster(phenotype.file, ids, 
                       family, covariates, phenotype, status.col, 
                       split.col = NULL, configs)
  
  if (is.null(foldids)){
    folds <- cut(seq(1, nrow(phe)), breaks=nfolds, labels=FALSE)
    foldids <- sample(folds, length(folds))
  }
  
  phe %>%
    dplyr::mutate(fold = foldids, tmp = 'val') %>%
    tidyr::pivot_wider(names_prefix = 'fold', names_from = fold, 
                       values_from = tmp, values_fill = list(tmp = 'train')) %>%
    data.table::fwrite(cv.phenotype.file, sep='\t')
  
  
  ### Get lambda path (hack for now...)
  full.lams.file <- file.path(configs$results.dir, "full.lams.txt")
  if (file.exists(full.lams.file)){
    full.lams <- scan(full.lams.file)
  } else{
    full.lams <- snpnet(genotype.pfile, phenotype.file, phenotype, family, 
                        covariates, weights, alpha, nlambda, lambda.min.ratio,
                        full.lams, split.col = NULL, p.factor, status.col, 
                        mem, configs, lambda_only = TRUE)
    
    write(full.lams, file = full.lams.file, ncolumns = 1)
  }

  cv_configs <- configs
  cv_configs$nCores <- cores_per_fold
  if (!is.null(mem)) {
    mem_cv <- mem*cores_per_fold
  } else {
    mem_cv <- NULL
  }
  cv_configs$mem <- mem_cv
  
  cvout = foreach::foreach(i = unique(foldids), .combine = 'rbind',
                           .packages = c("glmnet", "glmnetPlus", "snpnet")) %dopar% 
    {
      cv_configs$results.dir <- paste0(configs$results.dir, "/fold", i)
      outfile <- file.path(cv_configs$results.dir, "metric.val.txt")
      
      if (file.exists(outfile)) {
        out <- scan(outfile)
      } else {
        cv_configs$gcount.full.prefix <- NULL
        out <- snpnet(genotype.pfile, cv.phenotype.file, phenotype, family, 
                      covariates, weights, alpha, nlambda, lambda.min.ratio,
                      full.lams, split.col = paste0("fold", i), p.factor, status.col, 
                      mem_cv, cv_configs)$metric.val
        
        write(out, file = outfile, ncolumns = 1)
      }

      out
    }
  
  cvm <- apply(cvout, 2, mean)
  cvsd <- apply(cvout, 2, sd)
  lambda.min <- full.lams[which.max(cvm)]
  lambda_na <- apply(cvout, 2, function(x) !all(is.na(x)))
  
  ### Run snpnet on full training data set
  fit.lams = full.lams[lambda_na]
  
  if (!is.null(mem)) mem <- mem*ncores
  configs$nCores <- ncores
  configs$mem <- mem
  
  snpnet.object = snpnet(genotype.pfile, phenotype.file, phenotype, family, 
                         covariates, weights, alpha, nlambda, lambda.min.ratio,
                         fit.lams, split.col = NULL, p.factor, status.col, 
                         mem, configs)
  
  if (is.na(cvm[which.max(cvm) + 1]))
    warning("Cross-validation may have stopped early. Consider increasing stopping.lag.")
  
  if(! configs[['save']]) cleanUpCVFiles(configs)
  
  list(snpnet.object = snpnet.object, 
       cvm = cvm, 
       cvsd = cvsd, 
       cvout = cvout,
       lambda.min = lambda.min,
       full.lams = full.lams,
       fit.lams = fit.lams)
  
}