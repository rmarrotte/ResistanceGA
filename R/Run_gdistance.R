#' Get cost distance using gdistance
#' Execute gdistance
#'
#' @param gdist.inputs Object created from running \code{\link[ResistanceGA]{CS.prep}} function
#' @param r Accepts two types of inputs. Provide either the path to the resistance surface file (.asc) or specify an R RasterLayer object
#' @param scl scale the correction values (default is TRUE). No scaling should be done if the user wants to obtain absolute distance values as output. See \code{\link[gdistance]{geoCorrection}} for details
#' @return A costDistance matrix object from gdistance
#' @usage Run_gdistance(gdist.inputs, r, scl)

#' @export
#' @author Bill Peterman <Bill.Peterman@@gmail.com>
Run_gdistance <- function(gdist.inputs, 
                          r, 
                          scl = TRUE) {
  if (class(r)[1] != 'RasterLayer') {
    r <- raster(r)
  }
  
  # Make transition matrix
  tr <- transition(
    x = r,
    transitionFunction = gdist.inputs$transitionFunction,
    directions = gdist.inputs$directions
  )
  
  # Correct it
  if(gdist.inputs$method == 'costDistance'){
    trCorrect <- geoCorrection(tr, "c", scl = scl)
  }else{
    trCorrect <- geoCorrection(tr, "r", scl = scl)
  }
  
  # All pairs ?
  if(any(is.na(match(gdist.inputs$response$pop1,
                      gdist.inputs$response$pop2)))){
    #print("Running gdistance without all pairs!")
    # Gdistance internal function
    my_gdistance <- function(pair_IDs,method,trCorrect,response_df,sites_sp,scl){
      require(gdistance)
      pair_IDs <- as.numeric(response_df[pair_IDs,1:2])
      pairs_sp <- sites_sp[pair_IDs,]
      # Cost dist
      if (method == 'costDistance') {
        costDistance(trCorrect, pairs_sp)
        # Commute dist
      }else if (method == 'commuteDistance') {
        commuteDistance(trCorrect, pairs_sp) / 1000
      }else{
        stop("Method must be 'costDistance' or 'commuteDistance'")
      }
    }
    # Run each pair sequentially or in parallel
    if(gdist.inputs$ncores == 1){
      ret <- lapply(1:dim(response_df)[1],
                    FUN = my_gdistance, 
                    method = gdist.inputs$method,
                    trCorrect = trCorrect,
                    response_df = gdist.inputs$response_df,
                    sites_sp = gdist.inputs$sites_sp,
                    scl = scl)
    }else{
      ret <- sfLapply(1:dim(response_df)[1],
                      fun = my_gdistance, 
                      method = gdist.inputs$method,
                      trCorrect = trCorrect,
                      response_df = gdist.inputs$response_df,
                      sites_sp = gdist.inputs$sites_sp,
                      scl = scl)
    }
  }else{
    # Run all pairs
    # Cost dist
    if (gdist.inputs$method == 'costDistance') {
      ret <- costDistance(trCorrect, sites_sp)
      # Commute dist
    }else if (gdist.inputs$method == 'commuteDistance') {
      ret <- commuteDistance(trCorrect, sites_sp) / 1000
    }else{
      stop("Method must be 'costDistance' or 'commuteDistance'")
    }
    ret <- as.numeric(as.dist(ret))
  }
  
  # Merge the results to the response df
  results_df <- gdist.inputs$response_df
  results_df <- data.frame(results_df,"resistance"=unlist(ret))
  
  # Check for -1 and NA
  if(any(results_df$resistance == -1)){
    results_df <- results_df[results_df$resistance != -1,]
    print("Warning! -1 found in output and removed")
  }  
  if(any(is.na(results_df$resistance == -1))){
    results_df[is.na(results_df$resistance)] <- 0
    print("Warning! NA found in output and replaced with 0")
  }
  
  # Return results
  return(results_df)
}
