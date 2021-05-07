
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
