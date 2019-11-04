# Define a function for getting which mz values are 
# within an independent list of mz values
which_within <- function(x, y, mzd = 0.01){
  lapply(x, function(z){
    which(y >= z - mzd & y <= z + mzd)
  }
  )
}