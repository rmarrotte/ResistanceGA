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

CS.prep <- function(response_df, #dataframe: pop1, pop2, response"
                    sites_sp,
                    CS.program = '"C:/Program Files/Circuitscape/cs_run.exe"',
                    Neighbor.Connect = 8,
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
  if (!is.null(response_df)) {
    TEST.response <- ncol(response_df) == 3
    if (TEST.response == FALSE) {
      stop("The object 'response' is not in the form of a 3 column dataframe: pop1, pop2, GD")
    }
  }
  
  if (platform == 'pc') {
    parallel = FALSE
    cores = NULL
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
  suppressWarnings(ZZ <- ZZ.mat(response_df[,"pop1","pop2"]))
  
  # Make input list
  list(
    ZZ = ZZ,
    n.Pops = length(levels(response_df$pop1)),
    response_df = response_df,
    CS_Point.File = "sites.txt",
    CS.program = CS.program,
    Neighbor.Connect = Neighbor.Connect,
    platform = platform,
    pairs_to_include = pairs_to_include,
    parallel = parallel,
    cores = cores
  )
}
