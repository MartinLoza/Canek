
##Rad2Deg##
#Convert radians to degrees
# INPUT : Angle in radians
#
# OUTPUT : Angle in degree
Rad2Deg <- function(Angle_rad = NULL){
  return(Angle_rad*180/(pi))
}

## Cosine Normalization
CosNorm <- function(data){
  sizeFactor <- sqrt(colSums(data^2))
  normData <- t(t(data)/sizeFactor)
  return(normData)
}

#' RowScale
#'
#' @param matrix Numeric matrix to scale.
#' @param center Center the matrix.
#' @param scale Scale the matrix
#' @param filterRowZeros Filter the rows which contain only zeros.
#'
#' @return Scaled by rows matrix.
RowScale <- function(matrix, center = TRUE, scale = FALSE, filterRowZeros = FALSE){

  if((filterRowZeros == FALSE) & scale){
    if(length(WhichRowZeros(matrix = matrix)$names) != 0)
      stop("Rows with only zeros were found. Scaling is not possible.")
  }

  if(filterRowZeros){
    matrix <- FilterRowZeros(matrix = matrix)
  }

  if(center){
    rm <- matrixStats::rowMeans2(x = matrix, na.rm = TRUE)
    matrix <- matrix - rm
  }

  if(scale){
    rsd <- matrixStats::rowSds(matrix, na.rm = TRUE)
    matrix <- matrix/rsd
  }

  return(matrix)
}
