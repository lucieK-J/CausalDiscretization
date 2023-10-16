

library(corrplot)
require(arules)
require(stringr)
library(igraph)
library(RColorBrewer)
library(tidyr) # pour extract_numeric
library(SCCI)
library(bnlearn)
library(sna)
library(caret)
library(entropy)






plot2 <- function(x, ...){
  if (is.null(x)){
    plot.new()
    text(x = 0.5, y = 0.5, labels = "X")
  }
  else plot(x, ...)
}

round2 <- function(x, ...){
  if (is.character(x)){
    return(x)
  } else return(round(x, ...))
}

barDiagAndDiscret <- function(data, mapping, output = output...) {
  print(mapping) 
  X = eval_data_col(data, mapping$x)
  print(X)

  
  discr = mdlp(cbind(X, output))
  cuts = discr$cutp
  print(cuts)
  p <- ggplot(data, mapping) + geom_histogram(alpha=0.5)

  if (cuts != "All"){
    p <- p +
      geom_vline(xintercept = cuts[[1]])
  }
  p
}



reformate <- function(ADJ, sortiePC){
  # this function standardizes ADJ and PC output objects
  # ADJ contains all possible variables but PC output may have deletions 
  # it's a matter of completing with empty columns. 
  dat1=as(sortiePC, "amat")
  # ref = paste0("V", colnames(ADJ))
  ref = colnames(ADJ)
  d = length(ref)
  namessortiePC = colnames(dat1)
  missingVarInd = which(!ref %in% namessortiePC)
  dat2 = matrix(data = NA, nrow = d, ncol = d)
  for (i in 1:d){
    for (j in 1:d){
      if ((i %in% missingVarInd)|(j %in% missingVarInd)){
        dat2[i,j] = 0 
      } else dat2[i,j] = dat1[namessortiePC == paste0("V", i), namessortiePC == paste0("V", j)]
      
    }
  }
  
  return(dat2)
}


reformate2 <- function(ADJ, sortiePC, autreMat){
  # this function standardizes ADJ and outputPC objects
  # vesion 2 also converts matrices such as pvals 
  # ADJ contains all possible variables but PC output can have deletions
  # here, it's a matter of completing with empty columns.
  
  dat1=as(sortiePC, "amat")
  
  ref = colnames(ADJ)
  d = length(ref)
  namessortiePC = colnames(dat1)
  missingVarInd = which(!ref %in% namessortiePC)
  dat2 = matrix(data = NA, nrow = d, ncol = d)
  autreMat2 =  matrix(data = NA, nrow = d, ncol = d)
  for (i in 1:d){
    for (j in 1:d){
      if ((i %in% missingVarInd)|(j %in% missingVarInd)){
        dat2[i,j] = 0 
        autreMat2[i,j] = NA 
      } else {dat2[i,j] = dat1[namessortiePC == paste0("V", i), namessortiePC == paste0("V", j)]
      autreMat2[i,j] =  autreMat[namessortiePC == paste0("V", i), namessortiePC == paste0("V", j)] }
      
    }
  }
  return(list(Greformat = dat2, Mreformat = autreMat2))
}


evaluation <-function(ADJ, sortiePC){
  
  
  d = dim(ADJ)[1]
  
  if (is.numeric(sortiePC)){
    dat2 = sortiePC
  } else dat2 = reformate(ADJ, sortiePC)
  
  dat2 = t(dat2)
  shd_pcalg = pcalg::shd(as(ADJ, "graphNEL"), as(dat2, "graphNEL"))
  
  hdist_sna = hdist(ADJ, dat2= dat2, g1=NULL, g2=NULL, normalize=FALSE, diag=FALSE, mode="digraph")
  
  trulyPresent = sum(ADJ)
  
  trueAbsent = d*d-d-trulyPresent
  
  trulyPresentANDfound = 0
  trulyAbsentANDErased = 0
  
  for (i in 1:d){
    for (j in 1:d){
      if (i == j){next}
      true1 = ADJ[i,j] == 1
      pred1 = dat2[i,j] == 1
      if ( true1&pred1 ){
        trulyPresentANDfound = trulyPresentANDfound + 1
      }
      
      if ( (!true1)&(!pred1) ){
        trulyAbsentANDErased = trulyAbsentANDErased + 1
      }
      
    }}
  
  truePos = trulyPresentANDfound / trulyPresent
  
  trueNeg = trulyAbsentANDErased / trueAbsent
  

  
  dat2Sqlt =  ADJSqlt = matrix(0, dim(ADJ)[1], dim(ADJ)[2])
  for (i in 1:dim(ADJSqlt)[1]){
    for (j in 1:dim(ADJSqlt)[2]){
      if (i>j){next}
      ADJSqlt[i,j] = max(ADJ[i,j], ADJ[j,i])
      dat2Sqlt[i,j] = max(dat2[i,j], dat2[j,i])
    }
  }
  
  trulyPresentSqlt = sum(ADJSqlt)
  
  trueAbsentSqlt = d*(d-1)/2-trulyPresentSqlt
  
  trulyPresentANDfoundSqlt = 0
  trulyAbsentANDErasedSqlt = 0
  
  for (i in 1:d){
    for (j in 1:d){
      if (i == j){next}
      if (i>j){next}
      true1 = ADJSqlt[i,j] == 1
      pred1 = dat2Sqlt[i,j] == 1
      
      if ( true1&pred1 ){
        trulyPresentANDfoundSqlt = trulyPresentANDfoundSqlt + 1
      }
      
      if ( (!true1)&(!pred1) ){
        
        trulyAbsentANDErasedSqlt = trulyAbsentANDErasedSqlt + 1
      }
      
    }}
  
  truePosSqlt = trulyPresentANDfoundSqlt / trulyPresentSqlt
  
  trueNegSqlt = trulyAbsentANDErasedSqlt / trueAbsentSqlt
  
  return(list(shd_pcalg = shd_pcalg, hdist_sna = hdist_sna, truePos = truePos, trueNeg = trueNeg, 
              truePosSqlt = truePosSqlt, trueNegSqlt = trueNegSqlt))
}

boundaryFrom2interval <- function(cats){
  temp = str_split(cats, ",")
  temp[[1]][2]
  temp2 = temp[[2]][1]
  boundary = readr::parse_number(temp2)
  return(boundary)
}


mod.pc <-function(...){
  arguments <- list(...)
  nlevInput = arguments$suffStat$nlev
  labelsInput = arguments$labels
  
  dm = arguments$suffStat$dm
  if( is.null(dm)){return(do.call(pc, arguments))}
  # do.call(pc, arguments)
  
  d = dim(dm)[2]
  indToErase = rep()
  for( v in 1:d ){
    if (length(table(dm[,v])) == 1){
      indToErase = c(indToErase, v)
    }
  }
  
  if( length(indToErase) == 0 ){
    return( do.call(pc, arguments) )
  }
  
  dm2 = dm[,-indToErase]
  labels2 = arguments$labels[-indToErase]
  nlev2 = arguments$suffStat$nlev[-indToErase]
  arguments$suffStat$dm = dm2
  arguments$labels = labels2
  arguments$suffStat$nlev = nlev2
  
  
  if (length(arguments$suffStat$nlev)<2){
    ADJ0 = matrix(data = 0, length(nlevInput), length(nlevInput))
    colnames(ADJ0) =  rownames(ADJ0) = labelsInput
    return(ADJ0 )
  }
  
  return( do.call(pc, arguments) )
}

