#' Prepare data for optimization using \code{gdistance}
#'
#' Creates a necessary input for optimizing resistance surfaces based on pairwise cost distances, implemented using the \code{gdistance} library
#'
#' @param n.Pops The number of populations that are being assessed
#' @param response Vector of pairwise genetic distances (lower half of pairwise matrix).
#' @param samples Either provide the path to a .txt file containing the xy coordinates, or provide a matrix with x values in column 1 and y values in column 2. Alternatively, you can provide a \code{\link[sp]{SpatialPoints}} object
#' @param transitionFunction The function to calculate the gdistance TransitionLayer object. See \code{\link[gdistance]{transition}}. Default = function(x) 1/mean(x)
#' @param directions Directions in which cells are connected (4, 8, 16, or other). Default = 8
#' @param longlat Logical. If true, a \code{\link[gdistance]{geoCorrection}} will be applied to the transition  matrix. Default = FALSE
#' @param method Specify whether pairwise distance should be calulated using the \code{\link[gdistance]{costDistance}} or \code{\link[gdistance]{commuteDistance}} (Default) functions. \code{\link[gdistance]{costDistance}} calculates least cost path distance, \code{\link[gdistance]{commuteDistance}} is equivalent (i.e. nearly perfectly correlated with) resistance distance calculated by CIRCUITSCAPE.
#' @return An R object that is a required input into optimization functions

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>

#' @usage gdist.prep(n.Pops, 
#'                   response = NULL,
#'                   samples,
#'                   transitionFunction = function(x)  1 / mean(x),
#'                   directions = 8,
#'                   longlat = FALSE,
#'                   method = 'commuteDistance')

gdist.prep <-
  function(response_df, # dataframe: pop1, pop2, GD"
           sites_sp,
           transitionFunction = function(x) 1 / mean(x),
           directions = 8,
           longlat = FALSE,
           method = 'commuteDistance') {
    
    
    if (method != 'commuteDistance') {
      method <- 'costDistance'
    }
    
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
    
    # Make ZZ Mat
    suppressWarnings(ZZ <- ZZ.mat(response_df[,"pop1","pop2"]))
    
    ret <- list(response_df = response_df,
                sites_sp = sites_sp,
                transitionFunction = transitionFunction,
                directions = directions,
                ZZ = ZZ,
                longlat = longlat,
                method = method)
   
  }
