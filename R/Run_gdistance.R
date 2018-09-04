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
  
   
  tr <- transition(
    x = r,
    transitionFunction = gdist.inputs$transitionFunction,
    directions = gdist.inputs$directions
  )
  
  ret <- c()
    for(i in 1:dim(gdist.inputs$response_df)[1]){
      pairID <- as.numeric(gdist.inputs$response_df[i,1:2])
      pairCoords <- coordinates(gdist.inputs$sites_sp[pairID,])
      # Commute dist
      if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'costDistance') {
        trC <- geoCorrection(tr, "c", scl = scl)
        ret <- c(ret,costDistance(trC, pairCoords))
      }    
      if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'commuteDistance') {
        trR <- geoCorrection(tr, "r", scl = scl)
        ret <- c(ret,commuteDistance(trR, pairCoords) / 1000)
      }
    }
    rm(i)
  
  # Merge the results to the response df
  results_df <- gdist.inputs$response_df
  results_df <- data.frame(results_df,"resistance"=ret)
  
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
