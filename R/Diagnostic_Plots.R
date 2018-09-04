##################################
#' Create diagnostic plots
#'
#' This function will generate mixed effect model diagnostic plots following optimization
#'
#' @param resistance.mat Path to CIRCUITSCAPE "_resistances.out" file, or costDistance object created from running gdistance
#' @param genetic.dist Vector of pairwise genetic distances (lower half of pairwise matrix). Can be input as CS.inputs$response
#' @param XLAB Label for x-axis (Defaults to "Estimated resistance")
#' @param YLAB Label for y-axis (Defaults to "Genetic distance")
#' @param plot.dir Directory to output TIFF of diagnostic plots
#' @param type Specify whether the optimized surface is "continuous" or "categorical"
#' @param name The name to be attached to the output file. Must be specified when getting diagnostic plots for gdistance models
#' @param ID The to_from ID list for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @param ZZ The sparse matrix object for the MLPE model. The function will automatically create this object, but it can be specified directly from the output of CS.prep or gdist.prep (Default = NULL)
#' @return A multipanel panel .tif including histogram of residuals and qqplot of fitted mixed effects model

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage Diagnostic.Plots(resistance.mat, genetic.dist, XLAB,YLAB, plot.dir, type, name, ID, ZZ)

Diagnostic.Plots <-
  function(results_df,
           XLAB = "Estimated resistance",
           YLAB = "Genetic distance",
           plot.dir,
           type = "categorical",
           name = NULL,
           ZZ) {
    
    if (is.null(name)) {
        NAME <- gsub(pattern = "*_resistances.out", "", x = (basename(resistance.mat)))
    }else{
      NAME <- name
    }
    
    # scale and keep
    results_df$cs.unscale <- results_df$resistance
    results_df <- scale(results_df$resistance, center = TRUE, scale = TRUE)     
     
    # Assign value to layer
    LAYER <- assign("LAYER", value = results_df$resistance)
      
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
    
    #######
    # Make diagnostic plots
    if (type == "categorical") {
      tiff(
        filename = paste0(plot.dir, NAME, "_DiagnosticPlots.tif"),
        width = 279,
        height = 215,
        units = "mm",
        compression = c("lzw"),
        bg = "white",
        res = 300
      )
      par(
        mfrow = c(2, 1),
        oma = c(0, 4, 0, 0) + 0.1,
        mar = c(4, 4, 1, 1) + 0.1
      )
      hist(residuals(Mod), xlab = "Residuals", main = "")
      qqnorm(resid(Mod), main = "")
      qqline(resid(Mod))
      dev.off()
      par(mfrow = c(1, 1))
    } else {
      tiff(
        filename = paste0(plot.dir, NAME, "_DiagnosticPlots.tif"),
        width = 279,
        height = 215,
        units = "mm",
        compression = c("lzw"),
        bg = "white",
        res = 300
      )
      par(
        mfrow = c(2, 2),
        oma = c(0, 4, 0, 0) + 0.1,
        mar = c(4, 4, 1, 1) + 0.1
      )
      plot(results_df$response ~ results_df$cs.unscale,
           xlab = XLAB,
           ylab = YLAB)
      abline(lm(results_df$response ~ results_df$cs.unscale))
      plot(residuals(Mod) ~ results_df$cs.unscale,
           xlab = XLAB,
           ylab = "Residuals")
      abline(lm(residuals(Mod) ~ results_df$cs.unscale))
      hist(residuals(Mod), xlab = "Residuals", main = "")
      qqnorm(resid(Mod), main = "")
      qqline(resid(Mod))
      dev.off()
      par(mfrow = c(1, 1))
    }
  }
