#' Simultaneous optimization of multiple resistance surfaces
#'
#' Optimize multiple resistance surfsaces simultaneously using genetic algorithms
#'
#' @param CS.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function. Defined if optimizing using CIRCUITSCAPE
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{gdist.prep}} function. Defined if optimizing using gdistance
#' @param GA.inputs Object created from running \code{\link[ResistanceGA]{GA.prep}} function
#' @return This function optimizes multiple resistance surfaces, returning a Genetic Algorithm (GA) object with summary information. Diagnostic plots of model fit are output to the "Results/Plots" folder that is automatically generated within the folder containing the optimized ASCII files. A text summary of the optimization settings and results is printed to the results folder.
#' @usage MS_optim(CS.inputs, gdist.inputs, GA.inputs)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
MS_optim <- function(CS.inputs = NULL,
                     gdist.inputs = NULL,
                     GA.inputs) {
  if (!is.null(GA.inputs$scale)) {
    stop(
      "This function should NOT be used if you intend to apply kernel smoothing to your resistance surfaces"
    )
  }
  k.value <- GA.inputs$k.value
  
  
  if (GA.inputs$parallel != FALSE) {
    warning(
      "\n CIRCUITSCAPE cannot be optimized in parallel. \n Ignoring parallel arguement. \n If you want to optimize in parallel, use least cost paths or commute time with gdistance.",
      immediate. = TRUE
    )
  }
  
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
  
  t1 <- proc.time()[3]
  
  multi.GA_nG <- ga(
    type = "real-valued",
    fitness = Resistance.Opt_multi,
    population = GA.inputs$population,
    selection = GA.inputs$selection,
    mutation = GA.inputs$mutation,
    pcrossover = GA.inputs$pcrossover,
    crossover = GA.inputs$crossover,
    pmutation = GA.inputs$pmutation,
    Min.Max = GA.inputs$Min.Max,
    GA.inputs = GA.inputs,
    CS.inputs = CS.inputs,
    gdist.inputs = gdist.inputs,
    lower = GA.inputs$ga.min,
    upper = GA.inputs$ga.max,
    popSize = GA.inputs$pop.size,
    maxiter = GA.inputs$maxiter,
    run = GA.inputs$run,
    #parallel = GA.inputs$parallel,
    keepBest = GA.inputs$keepBest,
    seed = GA.inputs$seed,
    suggestions = GA.inputs$SUGGESTS,
    quiet = GA.inputs$quiet
  )
  rt <- proc.time()[3] - t1
  
  if(dim(multi.GA_nG@solution)[1] > 1) {
    multi.GA_nG@solution <- t(as.matrix(multi.GA_nG@solution[1,]))
  }
  
  Opt.parm <- GA.opt <- multi.GA_nG@solution
  for (i in 1:GA.inputs$n.layers) {
    if (GA.inputs$surface.type[i] == "cat") {
      ga.p <-
        GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
      parm <- ga.p / min(ga.p)
      Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                     1])] <- parm
      
    } else {
      parm <-
        GA.opt[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i + 1])]
      Opt.parm[(GA.inputs$parm.index[i] + 1):(GA.inputs$parm.index[i +
                                                                     1])] <- parm
    }
  }
  multi.GA_nG@solution <- Opt.parm
  
  
  RAST <-
    Combine_Surfaces(
      PARM = multi.GA_nG@solution,
      GA.inputs = GA.inputs,
      rescale = TRUE,
      p.contribution = TRUE
    )
  
  p.cont <- RAST$percent.contribution
  RAST <- RAST$combined.surface
  
  NAME <- paste(GA.inputs$parm.type$name, collapse = ".")
  names(RAST) <- NAME
  
  # Run CS or gdistance
  if(is.null(CS.inputs)){
    results_df <- Run_gdistance(gdist.inputs, RAST)
  }else{
    results_df <- Run_CS(CS.inputs, 
                         r = RAST, 
                         EXPORT.dir = GA.inputs$Results.dir)
  }
  
  ifelse(length(unique(RAST)) > 15,
         type <- "continuous",
         type <- "categorical")
  
  Diagnostic.Plots(results_df,
                   plot.dir = GA.inputs$Plots.dir,
                   type = type,
                   ZZ = ZZ)
  
  # Run model
  MLPE_mod <- MLPE.lmm(results_df = results_df,
                       form = GA.inputs$form,
                       REML = F,
                       ZZ = ZZ)
  fit.stats <- suppressWarnings(r.squaredGLMM(MLPE_mod))        
  aic <- AIC(MLPE_mod)        
  LL <- logLik(MLPE_mod)
  
  if (k.value == 1) {
    k <- 2
  } else if (k.value == 2) {
    k <- length(Opt.parm) + 1
  } else if (k.value == 3) {
    k <- length(Opt.parm) + length(GA.inputs$layer.names) + 1
  } else {
    k <- length(GA.inputs$layer.names) + 1
  }
  
  AICc <- (-2 * LL) + (2 * k) + ((2 * k) * (k + 1)) / (n - k - 1)
  
  
  # Get parameter estimates
  MLPE.results <- MLPE.lmm_coef(
    res_list = list(results_df),
    form = GA.inputs$form,
    out.dir = GA.inputs$Results.dir,
    ZZ = ZZ
  )
  
  Result.txt(
    GA.results = multi.GA_nG,
    GA.inputs = GA.inputs,
    method = ifelse(!is.null(CS.inputs),"CIRCUITSCAPE","gDistance"),
    Run.Time = rt,
    fit.stats = fit.stats,
    optim = GA.inputs$method,
    k = k,
    aic = aic,
    AICc = AICc,
    LL = LL[[1]]
  )
  
  write.table(p.cont, file = paste0(GA.inputs$Results.dir, "Percent_Contribution.csv"), sep = ",",
              row.names = F,
              col.names = T)
  
  # save(multi.GA_nG, 
  #      file = paste0(GA.inputs$Results.dir, NAME, ".rda"))
  
  saveRDS(multi.GA_nG, 
          file = paste0(GA.inputs$Results.dir, NAME, ".rds"))
  
  #file.remove(list.files(GA.inputs$Write.dir, full.names = TRUE))
  
  #unlink(GA.inputs$Write.dir, recursive = T, force = T)
  
  k.df <- data.frame(surface = NAME, k = k)
  
  cd.list <- list(results_df)
  names(cd.list) <- NAME
  
  AICc.tab <- data.frame(surface = NAME,
                         obj = multi.GA_nG@fitnessValue,
                         k = k,
                         AIC = aic,
                         AICc = AICc,
                         R2m = fit.stats[[1]],
                         R2c = fit.stats[[2]],
                         LL = LL)
  
  colnames(AICc.tab) <-
    c(
      "Surface",
      paste0("obj.func_", GA.inputs$method),
      "k",
      "AIC",
      "AICc",
      "R2m",
      "R2c",
      "LL"
    )
  
  out <- list(GA.summary = multi.GA_nG,
              MLPE.model = MLPE_mod,
              AICc.tab = AICc.tab,
              cd = cd.list,
              percent.contribution = p.cont,
              k = k.df)
  
  return(out)
}