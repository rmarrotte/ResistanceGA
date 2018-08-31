#' Prepare and bundle input CIRCUITSCAPE model parameters
#'
#' This function will prepare objects needed for running optimization functions
#'
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param CS_Point.File The path to the Circuitscape formatted point file. See Circuitscape documentation for help.
#' @param CS.program The path to the CIRCUITSCAPE executable file (cs_run.exe) on a Windows PC. If using a Linux or Mac system, provide the full path to the "csrun.py" file. See details below.
#' @param Neighbor.Connect Select 4 or 8 to designate the connection scheme to use in CIRCUITSCAPE (Default = 8)
#' @param pairs_to_include Default is NULL. If you wish to use the advanced CIRCUITSCAPE setting mode to include or exclude certain pairs of sample locations, provide the path to the properly formatted "pairs_to_include.txt" file here. Currently only "include" method is supported.
#' @param platform What computing platform are you using ("pc", "other"). This code has only been tested on Windows PC!!!
#' @param parallel (Logical) If using Linux / Ubuntu, do you want to run CIRCUITSCAPE in parallel?
#' @param cores If using Linux / Ubuntu and `parallel = TRUE`, how many cores should be used for parallel processing?

#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
#' @usage CS.prep(n.Pops, 
#' response, 
#' CS_Point.File, 
#' CS.program, 
#' Neighbor.Connect, 
#' pairs_to_include, 
#' platform, 
#' parallel, 
#' cores)
#' @details \code{CS.program} Example of path to CIRCUITSCAPE executible on Windows:
#'
#' '"C:/Program Files/Circuitscape/cs_run.exe"'
#'
#' ***NOTE: Double quotation used***
#' This is the current default for \code{CS.program}, but the directory may need to be changed depending upon your installation of CIRCUITSCAPE
#'
#' Linux
#' To call CIRCUITSCAPE from R on with Linux, first change file permissions from the command line terminal (shortcut: control + alt+ t) ` sudo chomod 755 /usr/local/bin/csrun.py `
#' Then specify \code{CS.program} as `csrun.py`
#'
#' Only with Linux, \code{parallel} can be set to \code{TRUE}, and the number of cores to run in parallel can be specified with \code{cores}.
#'
#' The Linux and Mac versions are in development. Please let me know if you encounter errors.

CS.prep <- function(n.Pops,
                    response = NULL,
                    sites_sp,
                    CS.program = '"C:/Program Files/Circuitscape/cs_run.exe"',
                    Neighbor.Connect = 8,
                    pairs_to_include = NULL,
                    platform = 'pc',
                    parallel = FALSE,
                    cores = NULL) {
  CS.exe_Test <- gsub("\"", "", CS.program)
  
  # CS path
  if (platform == 'pc') {
    if (!file.exists(gsub("\"", "", CS.program))) {
      stop("The specified path to 'cs_run.exe' is incorrect")
    }
  }
  
  # Point file
  sites_df <- data.frame("ID"=1:length(sites_sp),coordinates(sites_sp))
  write.table(sites_df,"sites.txt", row.names = F, col.names = F)
  
  # response data
  if (!is.null(response)) {
    TEST.response <- (is.vector(response) || ncol(response) == 1)
    if (TEST.response == FALSE) {
      stop("The object 'response' is not in the form of a single column vector")
    }
  }
  
  if (platform == 'pc') {
    parallel = FALSE
    cores = NULL
  }
  
  # Take pairs and write them in the CS format to the work dir
  if (!is.null(pairs_to_include)) {
    colnames(pairs_to_include) <- c("mode","include")
    pairs_to_include <- pairs_to_include[order(pairs_to_include$mode,pairs_to_include$include),]
    write.table(pairs_to_include,"pairs_to_include.txt",quote = F,row.names = F)
    colnames(pairs_to_include) <- c("pop1","pop2")
    ID <- pairs_to_include
    ID$pop1 <- factor(ID$pop1,levels=sort(unique(c(ID$pop1,ID$pop2)))) # Necessary for ZZ
    ID$pop2 <- factor(ID$pop2,levels=sort(unique(c(ID$pop1,ID$pop2))))
    pairs_to_include <- "pairs_to_include.txt"
    
  } # close pairs to include statement
  
  # Make to-from population list if include is null
  if (!exists(x = "ID")) {
    ID <- To.From.ID(n.Pops)
  }
  
  # Make ZZ Mat
  suppressWarnings(ZZ <- ZZ.mat(ID))
  
  # Make input list
  list(
    ID = ID,
    ZZ = ZZ,
    response = response,
    CS_Point.File = "sites.txt",
    CS.program = CS.program,
    Neighbor.Connect = Neighbor.Connect,
    n.Pops = n.Pops,
    platform = platform,
    pairs_to_include = pairs_to_include,
    parallel = parallel,
    cores = cores
  )
}
