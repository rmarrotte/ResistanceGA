# Run Mixed effects models, recover parameter estimates
#' Run maximum likelihood population effects mixed effects model (MLPE)
#'
#' Runs MLPE as detailed by Clarke et al. (2002). This function will run the model and return lmer object
#'
#' @param resistance Path to pairwise resistance distance matrix (resistances.out) from CS results. Alternatively, provide the pairwise resistances created from optimizing with `gdistance` (result of Run_gdistance).
#' @param pairwise.genetic Lower half of pairwise genetic distance matrix
#' @param REML Logical. If TRUE, mixed effects model will be fit using restricted maximum likelihood. Default = FALSE
#' @param ID The to_from ID list for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param ZZ The sparse matrix object for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param scale Specify whether the pairwise distance values be scaled and centered (Default = TRUE)
#' @return A lmer object from the fitted model
#' @details An AIC value will only be returned if \code{REML = FALSE}

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage MLPE.lmm(resistance, 
#'                 pairwise.genetic, 
#'                 REML = FALSE, 
#'                 scale = TRUE)
#' @references Clarke, R. T., P. Rothery, and A. F. Raybould. 2002. Confidence limits for regression relationships between distance matrices: Estimating gene flow with distance. Journal of Agricultural, Biological, and Environmental Statistics 7:361-372.

MLPE.lmm <-
  function(results_df, # assume resistance is a 4 column df: Pop1, Pop2, response, resistance
           REML = FALSE,
           ZZ,
           scale = T) {
    
    # Scale?      
    if(scale) {
      results_df$resistance <- scale(results_df$resistance, center = TRUE, scale = TRUE)
    }       
    
    # If pops are not BOTH found in pop1 and pop2, ignore the ZZ mat
    if(any(!is.na(match(results_df$pop1,results_df$pop2)))){
      # Fit model
      mod <- lFormula(response ~ resistance + (1 | pop1),
              data = results_df,
              REML = REML)
      mod$reTrms$Zt <- ZZ
    }else{
      mod <- lFormula(response ~ resistance + (1 | pop1) + (1 | pop2),
               data = results_df,
               REML = REML)
    }
    dfun <- do.call(mkLmerDevfun, mod)
    opt <- optimizeLmer(dfun)
    MOD <- (mkMerMod(environment(dfun), opt, mod$reTrms, fr = mod$fr))
    return(MOD)
  }

MLPE.lmm_coef <-
  function(res_list, # List of all distances
           genetic.dist, # df: pop1, pop2, response
           out.dir = NULL,
           ZZ) {
    
    COEF.Table <- array()
    for (i in 1:length(res_list)) {
      dat <- res_list[[i]]
      MOD <- MLPE.lmm(results_df=dat,ZZ=ZZ)
      Mod.Summary <- summary(MOD)
      COEF <- Mod.Summary$coefficients
      row.names(COEF) <- c("Intercept", resist.names[i])
      COEF.Table <- rbind(COEF.Table, COEF)
    }
        
    if (is.null(out.dir)) {
      COEF.Table <- (COEF.Table[-1, ])
    } else {
      COEF.Table <- COEF.Table[-1, ]
      write.table(
        COEF.Table,
        file = paste0(out.dir, "MLPE_coeff_Table.csv"),
        sep = ",",
        row.names = T,
        col.names = NA)
      return(COEF.Table)
    }
 }
