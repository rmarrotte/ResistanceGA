#' Prepare and bundle input CIRCUITSCAPE model parameters to run in Julia
#'
#' This function will prepare objects needed for running optimization functions
#'
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix). Not necessary if only executing Julia run.
#' @param CS_Point.File Provide a \code{\link[sp]{SpatialPoints}} object containing sample locations. Alternatively, specify the path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param covariates Data frame of additional covariates that you want included in the MLPE model during opitmization.
#' @param JULIA_HOME Path to the folder containing the Julia binary (See Details)
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide the path to the properly formatted "pairs_to_include.txt" file here. Currently only "include" method is supported.
#' @param parallel (Logical; Default = FALSE) Do you want to run CIRCUITSCAPE in parallel?
#' @param cores If `parallel = TRUE`, how many cores should be used for parallel processing?
#' @param cholmod (Logical; Default = TRUE). Should the cholmod solver be used? See details.
#' @param precision (Logical; Default = FALSE). Should experimental single precision method be used? See details.
#' @param run_test (Logical; Default = TRUE). Should a test of Julia Circuitscape be conducted? (This can take several seconds to complete)
#' @param write.files (Default = NULL). If a directory is specified, then the .ini and .asc files used in the CS.jl run will be exported.
#' @param write.criteria Criteria for writing .ini and .asc files. If a time in seconds is not specified, then all files will be written if a \code{write.files} directory is specified.
#' @param silent (Default = TRUE) No updates or logging of CIRCUITSCAPE will occur. May be useful to set to FALSE to debug. 

#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage jl.prep(n.Pops, 
#' response, 
#' CS_Point.File, 
#' covariates = NULL,
#' JULIA_HOME,
#' Neighbor.Connect, 
#' pairs_to_include, 
#' platform, 
#' parallel, 
#' cores,
#' cholmod,
#' precision, 
#' run_test,
#' write.files = NULL,
#' write.criteria = NULL,
#' silent = TRUE)
#' @details 
#' This function requires that Julia is properly installed on your system. Upon first running of this function, the Circuitscape.jl library will be downloaded and tested. (see https://github.com/Circuitscape/Circuitscape.jl for more details). This may take some time.
#' 
#' Using \code{cholmod} (see https://github.com/Circuitscape/Circuitscape.jl)
#' "The cholesky decomposition is a direct solver method, unlike the algebraic multigrid method used by default in both the old and the new version. The advantage with this new direct method is that it can be much faster than the iterative solution, within a particular problem size. 
#' Word of caution: The cholesky decomposition is not practical to use beyond a certain problem size because of phenomenon called fill-in, which results in loss of sparsity and large memory consumption."
#' The cholmod solver can only be used when \code{precision} `= FALSE` (double precision).
#' 
#' If \code{precision} is TRUE, then the EXPERIMENTAL single precision method will be used. Single precision usually uses less memory, but is likely to reduce accuracy. NOTE: Preliminary testing of single precision mode in a Windows pc resulted in extremely slow runs.
#' 
#' \code{JULIA_HOME} is where the Julia binary files are stored. Usually in a `bin` directory within the Julia install directory.
#' 
jl.prep <- function(response_df, #dataframe: pop1, pop2, response"
                    sites_sp,
                    Neighbor.Connect = 8,
                    cholmod = TRUE,
                    precision = FALSE,
                    run_test = TRUE,
                    write.files = NULL,
                    write.criteria = NULL,
                    silent = TRUE) {
  
  
  # Determine if CIRCUITSCAPE package is installed
  if(julia_installed_package("Circuitscape") == 'nothing') {
    julia_install_package('Circuitscape')
    julia_call('Pkg.test', "Circuitscape")
  } 
  
  # Format inputs -----------------------------------------------------------
  if(isTRUE(precision)) {
    precision <- 'single'
  } 
  
  if(isTRUE(cholmod) && (precision == 'single')) {
    stop(cat(paste0('\n', 
                    "jl.prep ERROR:", '\n',
                    "CHOLMOD solver only works when using double precision. Set either `cholmod = FALSE` OR `precision = FALSE` to proceed", 
                    '\n', '\n' ))
    )
  }
  
  if(isTRUE(cholmod)) {
    solver <- 'cholmod'
  } else {
    solver <- NULL
  }
  
  # Point file
  sites_df <- data.frame("ID"=1:length(sites_sp),coordinates(sites_sp))
  write.table(sites_df,"sites.txt", row.names = F, col.names = F)
  
  # response data
  if (!is.null(response_df)) {
    TEST.response <- ncol(response_df) == 3
    if (TEST.response == FALSE) {
      stop("The object 'response' is not in the form of a 3 column dataframe: pop1, pop2, GD")
    }
  }
  
  # Sort response df
  response_df <- response_df[order(response_df$pop1,response_df$pop2),]
  row.names(response_df) <- NULL 
  
  # Make factors
  response_df$pop1 <- factor(response_df$pop1,levels=sort(unique(c(response_df$pop1,response_df$pop2)))) # Necessary for ZZ
  response_df$pop2 <- factor(response_df$pop2,levels=sort(unique(c(response_df$pop1,response_df$pop2))))
  
  # Take pairs and write them in the CS format to the work dir
  pairs_to_include <- response_df[,c("pop1","pop2")]
  colnames(pairs_to_include) <- c("mode","include")
  write.table(pairs_to_include,"pairs_to_include.txt",quote = F,row.names = F) 
  pairs_to_include <- "pairs_to_include.txt"  
  
  # Make ZZ Mat
  suppressWarnings(ZZ <- ZZ.mat(response_df[,c("pop1","pop2")]))
  
  list(
    ZZ = ZZ,
    n.Pops = length(levels(response_df$pop1)),
    response_df = response_df,
    CS_Point.File = "sites.txt",
    Neighbor.Connect = Neighbor.Connect,
    solver = solver,
    precision = precision,
    write.files = write.files,
    write.criteria = write.criteria,
    silent = silent
  )
}