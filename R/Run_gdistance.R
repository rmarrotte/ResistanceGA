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
  
  # All pairs
  if(is.null(gdist.inputs$pairs_to_include)){
    if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'costDistance') {
      trC <- geoCorrection(tr, "c", scl = scl)
      ret <- costDistance(trC, gdist.inputs$samples)
    }
    
    if (gdist.inputs$longlat == TRUE | gdist.inputs$directions >= 8 & gdist.inputs$method == 'commuteDistance') {
      trR <- geoCorrection(tr, "r", scl = scl)
      ret <- commuteDistance(trR, gdist.inputs$samples) / 1000
    }    
  }else{ # Specific pairs
    # For each pair
    ret <- c()
    for(i in 1:dim(gdist.inputs$pairs_to_include)[1]){
      pairID <- as.numeric(gdist.inputs$pairs_to_include[i,])
      pairCoords <- gdist.inputs$samples[pairID,]
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
  }
  return(ret)
}
