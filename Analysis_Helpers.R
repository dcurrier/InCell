# # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ## # # # # # ##  ### # ##
f_cvmts <- function(x, y) {
  # (1) Sort x and y
  xS <- sort(x)
  yS <- sort(y)
  # (2) Get Ranks within x and y
  i <- rank(xS)
  j <- rank(yS)
  xDF <- data.frame(val = xS, ID = "x")
  yDF <- data.frame(val = yS, ID = "y")
  # (3) Join x and y
  DF	<- rbind(xDF, yDF)
  # (4) Get Ranks within Combined Data
  r <- rank(DF$val)[which(DF$ID == "x")]
  s <- rank(DF$val)[which(DF$ID == "y")]
  # (5) Calculate N and M
  N <- length(x)
  M <- length(y)
  # (6) Calculate U
  U <- N*sum((r -i)^2) + M*sum((s - j)^2)
  # (7) Calculate T
  T <- U/N/M/(N+M) - (4*M*N-1)/(6*(M+N))
  # (8) Scale by Maximum Theoretical Value
  Tmax <- (2*M*N + 1)/(6*(M+N))
  Tscaled <- T / Tmax
  # (9) Correct Sign
  Tscaled <- Tscaled * sign(median(y) - median(x))
  return(Tscaled)
}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #

## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
f_ks <- function(x, y) {
  ksTwo <- as.vector(ks.test(x, y, exact=F)$statistic)
  ksGreater <- as.vector(ks.test(x, y, alternative="gr", exact=F)$statistic)
  if (ksTwo == ksGreater) {
    ks <- ksTwo
  }else{
    ks <- -ksTwo
  }
  return(ks)
}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #

## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #
is.even <- function(x){
  if( !(typeof(x) %in% c("integer", "double")) ) stop("argument 'x' must be an integer")
  if( typeof(x) == "double" && floor(x) != x ) stop("argument 'x' must be an integer")
  x %% 2 == 0
}
## # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # ### # # # # # ##  # ## # # # # # #