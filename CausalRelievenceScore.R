
library(stringr)

compute.score <- function(u, alpha, ns){

  res = list()
  temp = which(unlist(lapply(u, is.character)))
  u = u[temp]
  u = u[str_detect(unlist(u), "x=")]
  
  
  
  
  if (length(u) == 0){
    res$ScorePertienceCausale = NA
    res$ScorePertienceCausale2 = NA
    res$sumAsupp = NA
    res$sumArest = NA
    res$ns = ns
    res$pvals = NA
    res$normalizedpvals = NA
    res$ne = 0 
    return(res)}
    ucut = unlist(str_split(u[which(unlist(lapply(u, is.character)))], "\\n"))
 
  ucut2 =  ucut[!str_detect(ucut, "Order")]
 
  
  
  ucut2b = ucut2[ucut2 != ""]

  
  
  extractTest <- function(v){
    x = as.numeric(str_match(v, "x=\\s*(.*?)\\s*y")[2])
    y = as.numeric(str_match(v, "y=\\s*(.*?)\\s*S")[2])
    S  = as.numeric(unlist(str_split(str_match(v, "S=\\s*(.*?)\\s*:")[2], " ")))
    pval = as.numeric(str_match(v, "pval =\\s*(.*?)\\s* ")[2])
    return(list(x = x, y = y, S = S, pval = pval))}

  listTests = lapply(ucut2b,extractTest)
  
  if (length(listTests) == 1){
    res$ScorePertienceCausale = NA
    res$ScorePertienceCausale2 = NA
    res$sumAsupp = NA
    res$sumArest = NA
    res$ns = ns
    res$pvals = NA
    res$normalizedpvals = NA
    res$ne = NA 

    message("impossible computation")
    return(res)}
   test_x = unlist(lapply(listTests, function(e) e$x))
  test_y = unlist(lapply(listTests, function(e) e$y))
  test_S = lapply(listTests, function(e) e$S)
  test_pval = unlist(lapply(listTests, function(e) e$pval))
  
  nvar = length(table(c(test_x, test_y)))
  EdgeSpec = matrix(NA, nrow = nvar, ncol = nvar)

  
  for (i in 1:nvar){
    for (j in 1:nvar){
      indEdgesij = (test_x == i)&(test_y == j)
      pvaluesEdgesij = test_pval[indEdgesij]
      
      indEdgesji = (test_x == j)&(test_y == i)
      pvaluesEdgesji = test_pval[indEdgesji]
      
      EdgeSpec[i,j] = max(c(pvaluesEdgesij, pvaluesEdgesji), na.rm = T)
    } 
  }
  
  pvalTOComputeScore = EdgeSpec
  pvalTOComputeScore[upper.tri(EdgeSpec, diag = TRUE)] <- NA
  
  pvalTOComputeScore
  
  ls = length(listTests)
  

  
  pvalTOComputeScore < alpha
  
  aa = which(pvalTOComputeScore >= alpha, arr.ind = T) 
  ab = which(pvalTOComputeScore < alpha, arr.ind = T)
  
  
  S = 0
  
  sumaa = 0
  sumab = 0 
  sums2 = 0
  
  normalizedpvals = pvalTOComputeScore
  
  if (length(aa) != 0){
  for (arrow in 1:dim(aa)[1]){
    temp = aa[arrow,]
    sumaa = sumaa + (pvalTOComputeScore[temp[1], temp[2]] - alpha)/(1 - alpha)
    sums2 = sums2 + pvalTOComputeScore[temp[1], temp[2]]
    
    normalizedpvals[temp[1], temp[2]] = (pvalTOComputeScore[temp[1], temp[2]] - alpha)/(1 - alpha)
  }
    sumaa = sumaa/dim(aa)[1] 
  } else (sumaa = 0)
  
  
  
  if (length(ab) != 0){
    for (arrow in 1:dim(ab)[1]){
      temp = ab[arrow,]
      sumab = sumab + (alpha - pvalTOComputeScore[temp[1], temp[2]])/(alpha)
      sums2 = sums2 + pvalTOComputeScore[temp[1], temp[2]]
      normalizedpvals[temp[1], temp[2]] =  (alpha - pvalTOComputeScore[temp[1], temp[2]])/(alpha)
        
    }
    sumab = sumab/dim(ab)[1]
  } else (sumab = 0)
  
  
  ns = nvar 

  ScorePertienceCausale = (sumaa + sumab)
  res$ScorePertienceCausale = ScorePertienceCausale
  res$ScorePertienceCausale2 = 1 - sums2/(ns*(ns-1)/2)
  
  res$sumAsupp = sumaa
  res$sumArest = sumab
  res$ns = ns
  res$pvals = pvalTOComputeScore
  res$normalizedpvals = normalizedpvals
  res$ne = dim(ab)[1] 
  return(res)
}
