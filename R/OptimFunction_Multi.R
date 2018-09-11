#' Optimize multiple resistance surfaces simultaneously
#'
#' Create composite resistance surface by simultaneously optimizing multiple categoricla and continuous surfaces. This optimization function is designed to be called from GA
#'
#' @param PARM Parameters to transform conintuous surface or resistance values of categorical surface. Should be a vector with parameters specified in the order of resistance surfaces. These values are selected during optimization if called within GA function.
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @param Min.Max Define whether the optimization function should minimized ('min') or maximized ('max')
#' @param quiet Logical, if TRUE, AICc and iteration time will not be printed to the screen at the completion of each iteration. Default = FALSE
#' @return AIC value from mixed effect model
#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Resistance.Opt_multi <- function(PARM,
                                 CS.inputs = NULL,
                                 gdist.inputs = NULL,
                                 GA.inputs,
                                 Min.Max,
                                 quiet = FALSE) {
  # Set global vars
  if(!is.null(gdist.inputs)){
    results_df <- gdist.inputs$response
    n <- gdist.inputs$n.Pops
    ZZ <- gdist.inputs$ZZ
  }else if(!is.null(CS.inputs)){
    results_df <- CS.inputs$response
    n <- CS.inputs$n.Pops
    ZZ <- CS.inputs$ZZ
  }else{stop("Need a proper method")}
  
  #print(PARM)
  t1 <- proc.time()[3]
  
  method <- GA.inputs$method
  EXPORT.dir <- GA.inputs$Write.dir
  ######
  #   r <- GA.inputs$Resistance.stack
  File.name = "resist_surface"

  r <- Combine_Surfaces(
    PARM = PARM,
    GA.inputs = GA.inputs,
    out = GA.inputs$Write.dir,
    File.name = File.name,
    rescale = FALSE
  )
  
  if(cellStats(r, "mean") == 0) { # Skip iteration
    obj.func.opt <- -99999
  } else { # Continue with iteration
    if(!is.null(CS.inputs)){
      cd <- try(Run_CS(CS.inputs, r = r ), TRUE)
    }else if(!is.null(gdist.inputs)){
      cd <- try(Run_gdistance(gdist.inputs, r), TRUE)
    }else{
      stop("Specify a proper distance method")
    }
    
    if(isTRUE(class(cd) == 'try-error') || isTRUE(exists('obj.func.opt'))) {
      obj.func.opt <- -99999
    }else{
      # Run mixed effect model on each Circuitscape effective resistance
      mod_mlpe <- suppressWarnings(MLPE.lmm(results_df = cd,
                                            form = GA.inputs$form,
                                            ZZ = ZZ,
                                            REML = FALSE,
                                            scale = T))
      if (method == "AIC") {            
        obj.func <- suppressWarnings(AIC(mod_mlpe))
        obj.func.opt <- obj.func * -1
      } else if (method == "R2") {
        obj.func <- suppressWarnings(r.squaredGLMM(mod_mlpe))
        obj.func.opt <- obj.func[[1]]
      } else {
        obj.func <- suppressWarnings(logLik(mod_mlpe))
        obj.func.opt <- obj.func[[1]]
      }
    }
  }
  
  rt <- proc.time()[3] - t1
  if (quiet == FALSE) {
    cat(paste0("\t", "Iteration took ", round(rt, digits = 2), " seconds"), "\n")
    cat(paste0("\t", method, " = ", round(obj.func.opt, 4)), "\n", "\n")
  }
  
  return(obj.func.opt)
  # OPTIM.DIRECTION(Min.Max)*(obj.func) # Function to be minimized/maximized
  
}