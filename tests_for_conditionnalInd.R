


# tests for continuous data

# Kernel conditional independence test
library(CondIndTests)
kci.test <- function(x, y, S, suffStat){
  
  dm = suffStat$dm
  dm = data.frame(dm)
  
  Y = dm[,x]
  E = dm[,y]
  X = dm[,S] # conditional var
  
  if (length(X) != 0){
    pv = KCI(Y, E, X)$pvalue # super long ... 
  }else  pv = KCI(Y, E, rep(0,dim(dm)[1]))$pvalue
  
  
  
  return(pv)

}

  


Jonckheere.Terpstra.test <- function(x, y, S, suffStat){
  # suffStat : list(dm)
  require(bnlearn)
  dm = suffStat$dm
  dm =  data.frame(dm)
  for (v in 1:dim(dm)[2]){
    dm[, v] = as.ordered(dm[, v])
  }
  
  # cat(paste("mon message : S =  ", S, " !"))
  # cat(paste("S est null ? :  ", is.null(S), " !"))
  # cat(paste("str(S) ? ", str(S) , " !"))
  
  
  if (length(S) == 0){
    p = bnlearn::ci.test(x = dm[,x] , y = dm[,y],  test = "jt" )$p.value
    # cat(paste("et p =", p, " !"))  
    } else {  p = bnlearn::ci.test(x = dm[,x] , y = dm[,y], z =  dm[,S],  test = "jt" )$p.value}
  
  p = as.numeric(p)
  return(p)
} 


Mutual.Information.test <- function(x, y, S, suffStat){
  # suffStat : list(dm)
  dm = suffStat$dm
  dm =  data.frame(dm)
  for (v in 1:dim(dm)[2]){
    dm[, v] = as.factor(dm[, v])
  }
  if (length(S) == 0){
    p = ci.test(x = dm[,x] , y = dm[,y], test = "mi" )$p.value
  } else   p = ci.test(x = dm[,x] , y = dm[,y], z = dm[,S], test = "mi" )$p.value
  
  return(p)
} 



Mutual.Information.test <- function(x, y, S, suffStat){
  # suffStat : list(dm)
  dm = suffStat$dm
  dm =  data.frame(dm)
  for (v in 1:dim(dm)[2]){
    dm[, v] = as.factor(dm[, v])
  }
  if (length(S) == 0){
    p = ci.test(x = dm[,x] , y = dm[,y], test = "mi" )$p.value
  } else   p = ci.test(x = dm[,x] , y = dm[,y], z = dm[,S], test = "mi" )$p.value
  
  return(p)
}




