

source("CausalRelievenceScore.R")
require(arules)
require(stringr)
options(encoding = "latin1") 
source(file = "STUCCO.R")
library(tidyr)


random.discretize <- function(x, seed = seed, m){
  if (!is.null(seed)){
    set.seed(seed)
  }
  varRange = range(x)
  Boundaries = c(0,1, runif(m-2, range(x)[1], range(x)[2]))
  res = base::cut(x, breaks = Boundaries)
  return(res)
}


causal.discretization <- function(data, m, Vtodiscr, esp, alpha, onlyInit = F, criteria = "spc", test, truth, seed = NULL, randomInit = T){
  
  datacont = data
  
   if (!is.null(seed)){
    set.seed(seed)
  }

  a = proc.time()
  V = dim(data)[2]
  N =  dim(data)[1]
  namesVars = colnames(data)
  Vtodiscr = c()
  nvar = V
  for (v in 1:V){
      if ((class(data[,v]) != "character") ){
      if (length(table(data[,v])) == 1){stop('No variable should be constant')}
      
      if (randomInit){
        data[,v] = random.discretize(data[,v], seed = seed, m)
      } else data[,v] = arules::discretize(data[, v], method = "frequency", breaks = m)
      
      Vtodiscr = c(Vtodiscr, v)
    } else { data[,v] = factor(data[,v]) }
  }
  sum(is.na(as.matrix(data))) 
  
  data_initale = data
  if (onlyInit){
    return(data_initale)
  }

  
  newdata = data

  boundaries = vector(mode = "list", length(Vtodiscr))
  names(boundaries) = Vtodiscr
  
  continue = T
  k = 1
  SCOREfusion = rep()
  SCORE = rep()
  TP =  rep()
  TN = rep()
  itnewboundVar = rep()
  itnewboundVal = rep()
  itnewfusionVar = rep()
  itnewfusionVal = rep()
  
  while( continue ){
  
    intervalCandidats = vector(mode = "list", length(Vtodiscr))
    names(intervalCandidats) = names(Vtodiscr)

    for (v in Vtodiscr){

      VcurentSupports = table(newdata[,v])
      VallsupportCandidat = rep(NA, (length(VcurentSupports)-1))
      VallnamesCandidat = rep(NA, (length(VcurentSupports)-1))
      if (length(VcurentSupports) == 1){ next}
      for (boundary in 1:(length(VcurentSupports)-1)){
        temp = VcurentSupports[c(boundary, boundary+1)]
        temp2 = paste0(names(temp),collapse = "")
        temp3 = str_split(temp2, ",")
        temp4 = temp3[[1]][c(1,3)]
        valboundary = boundaryFrom2interval(names(temp))
        
        if (valboundary %in% boundaries[[which(v == Vtodiscr)]]){
          VallsupportCandidat[[boundary]] = NA
          VallnamesCandidat[[boundary]] = "must not be merged"
          next()
        }
        
        nomCandidat = paste0(temp4, collapse = ", ")
        supportCandidat = sum(temp)
        VallsupportCandidat[[boundary]] = supportCandidat
        VallnamesCandidat[[boundary]] = nomCandidat
      }
      intervalCandidats[[v]] = list()
      intervalCandidats[[v]]$support = VallsupportCandidat
      intervalCandidats[[v]]$names = VallnamesCandidat
      
    }
    
    if (sum(unlist(lapply(intervalCandidats, is.null))) == length(intervalCandidats)){
      continue = F
      next
    }
    
    minimSupportV = rep(NA, length(Vtodiscr))
    for (v in Vtodiscr){
      minimSupportV[v == Vtodiscr] = min(intervalCandidats[[v]]$support, na.rm = T)
    }
    
    if (sum(is.finite(minimSupportV)) == 0){
      continue = F 
      break()
    }
    
    varCand = Vtodiscr[which(minimSupportV == min(minimSupportV, na.rm = T))][1]
    varCand
    
    
    candidateInterV = which.min(intervalCandidats[[varCand]]$support)
    

    candidateInter = list(v = varCand, name = intervalCandidats[[varCand]]$names[candidateInterV],
                          replacementofIntervals = c(candidateInterV, candidateInterV+1), support = intervalCandidats[[varCand]]$support[candidateInterV])
    
    
    
    candidateInter
    var = candidateInter$v
    table(newdata[,var])
    names(newdata)
    names(newdata[,-var])
    cats = levels(newdata[,var])[candidateInter$replacementofIntervals]
    groups = list(list(varindex = var, cat = cats[1]), list(varindex = var, cat = cats[2]))
    dataforstucco = data
    dataforstucco[,var] = newdata[,var]
    dm = newdata
    for (l in 1:nvar){
      dm[,l] = as.integer(newdata[,l])-1
    }
    dm = as.matrix(dm)
    
    
    
    res_mattodm = mattodm(m = data_initale)
    data = res_mattodm$dm
    nlev = res_mattodm$nlev
    suffStat = list(dm = dm, nlev  = nlev, adaptDF = F)
    
    
    pc.fit <-  mod.pc(suffStat = suffStat,
                      indepTest = get(test), 
                      alpha = alpha, labels = colnames(dm), verbose = F)
    
    
    invisible(u <- evaluate::evaluate({
      function(){pc.fit <-  mod.pc(suffStat = suffStat,
                                   indepTest = get(test),
                                   alpha = alpha, labels = colnames(dm), verbose = TRUE) }
    }, new_device = F) )
    
    
    score = compute.score(u = u, alpha = alpha, ns = d)
    
    
    if (criteria == "spc"){
      score = score$ScorePertienceCausale2
    }
    
    if (criteria == "sumArest"){
      score = score$sumArest
    }
    score
    
    
    
    SCORE = c(SCORE, score)
  
    
    ev = evaluation(truth, pc.fit)
    TP = c(TP, ev$truePosSqlt)
    TN = c(TN, ev$trueNegSqlt)

    newdatafusion = newdata
    
    sum(is.na(as.matrix(newdata)))
    sum(is.na(as.matrix(newdatafusion)))
 
    newdatafusion[,var] = factor(newdatafusion[,var])
    temp = levels(factor(newdatafusion[,var]))
    
    length(temp)
    length(table(newdatafusion[,var]))
    
    if( !is.na(sum(which(temp %in% cats)[1])) ){
      temp = temp[-which(temp %in% cats)[1]]
      temp[which(temp %in% cats)] = candidateInter$name
    }
    
    length(temp)
    length(table(newdatafusion[,var]))
    sum(is.na(as.matrix(newdatafusion)))
    
    sum(is.na(as.matrix(newdata)))
    sum(is.na(as.matrix(newdatafusion)))
    
    
    newdatafusion[,var] = as.character(newdatafusion[,var])
    
    newdatafusion[,var][which(newdatafusion[,var] %in% cats)] = candidateInter$name
    
    sum(is.na(as.matrix(newdatafusion)))

    dmfusion = newdatafusion

    isFusionPossible  = T 

    
    for (l in 1:nvar){
      dmfusion[,l] = as.numeric(factor(newdatafusion[,l]))-1
      
      if (length(table(dmfusion[,l]))<2){ isFusionPossible = F}
    }
    
    if (isFusionPossible){
      
      suffStat$dm = dmfusion

      ufusion <- evaluate::evaluate({
        function(){pc.fit <- pc.fit <-  mod.pc(suffStat = suffStat,
                                               indepTest = get(test),
                                               alpha = alpha, labels = colnames(dm), verbose = TRUE) }
      }, new_device = F) 
      
      
      scoreFusion = compute.score(u = ufusion, alpha, ns = d)
      if (criteria == "spc"){
        scoreFusion = scoreFusion$ScorePertienceCausale2
      }
      
      if (criteria == "sumArest"){
        scoreFusion = scoreFusion$sumArest
      }
      
      
      SCOREfusion = c(SCOREfusion, scoreFusion)
      
      
      fusion = scoreFusion - score > esp
      
    } else {
      scoreFusion = NA 
      fusion = FALSE
      SCOREfusion = c(SCOREfusion, scoreFusion)
    }
    
      bornes = c(readr::parse_number(str_split(candidateInter$name, ',')[[1]][1]), readr::parse_number(str_split(candidateInter$name, ',')[[1]][2]))
    
    
    if (fusion){
      temp = levels(newdata[,var])
      temp = temp[-which(temp %in% cats)[1]]
      temp[which(temp %in% cats)] = candidateInter$name
      
      
      newdata[,var] = as.character(newdata[,var])
      newdata[,var][which(newdata[,var] %in% cats)] = candidateInter$name
      newdata[,var] = factor(newdata[,var], levels = temp)
      
      itnewboundVar = c(itnewboundVar, NA)
      itnewboundVal = c(itnewboundVal, NA)
      itnewfusionVar = c(itnewfusionVar, var)
      itnewfusionVal = c(itnewfusionVal, candidateInter$name)
      
    } else {
    newBound = boundaryFrom2interval(cats)
    boundaries[[which(var == Vtodiscr)]] = c(boundaries[[which(var == Vtodiscr)]], newBound) # ; message("one boundary found")
    
    itnewboundVar = c(itnewboundVar, var)
    itnewboundVal = c(itnewboundVal, newBound)
    itnewfusionVar = c(itnewfusionVar, NA)
    itnewfusionVal = c(itnewfusionVal, NA)
    }
    
    k = k + 1
    
  }
  
  b = proc.time()
  
  return(list(boundaries = boundaries, newdata = newdata, Vtodiscr = Vtodiscr, time = a-b, 
              SCOREfusion = SCOREfusion, SCORE = SCORE,
              itnewboundVar = itnewboundVar, itnewboundVal = itnewboundVal,
              data_initale = data_initale, 
              itnewfusionVar = itnewfusionVar, 
              itnewfusionVal = itnewfusionVal,
              TP = TP, TN = TN))
}



causal.discretization.with.several.inits <- function(nInits, data, m, Vtodiscr, esp, alpha, onlyInit = F, criteria, test, truth, seed = NULL, randomInit = F, display = T, SIM){
  
  
  
  Res_Causal_for_each_init = list()
  for (initIndex in 1:nInits){
    print(paste0("---------------- init number ", initIndex, ' over ' , nInits))
    truth = SIM$ADJ
    Res_Causal = causal.discretization(data, m, Vtodiscr, esp, alpha, onlyInit = F, criteria = "spc", test, truth, seed = initIndex, randomInit = T)
    Res_Causal_for_each_init[[initIndex]] = Res_Causal

    
    
    
    
  }
  
  Res_Causal_for_each_init
  finalScoresInits = rep(NA, nInits)
  for (initIndex in 1:nInits){
    scores = Res_Causal_for_each_init[[initIndex]]$SCORE
    finalScoresInits[initIndex] = scores[length(scores)]
  }
  
  
  ibestDisc = which.max(finalScoresInits)
  nbit = length(Res_Causal_for_each_init[[ibestDisc]]$SCORE)
  
  evOptimised = list(truePosSqlt = Res_Causal_for_each_init[[ibestDisc]]$TP[nbit], trueNegSqlt = Res_Causal_for_each_init[[ibestDisc]]$TN[nbit])
  
  scoreOptimised =  Res_Causal_for_each_init[[ibestDisc]]$SCORE[nbit]
  

  SIM$DATAdiscr
  SIM$dm
  SIM$NameVars
  SIM$nlev
  suffStat = list(dm = SIM$dm, nlev = SIM$nlev, adaptDF = F)
  
  pc.fit.optimal <-  mod.pc(suffStat = suffStat,
                            indepTest = get(test), 
                            alpha = alpha, labels = SIM$NameVars, verbose = FALSE)

  
  invisible(u.optimal <- evaluate::evaluate({
    function(){pc.fit <-  mod.pc(suffStat = suffStat,
                                 indepTest = get(test), 
                                 alpha = alpha, labels = SIM$NameVars, verbose = TRUE) }
  }, new_device = F)) 
  
  scoreOptimal = compute.score(u = u.optimal, alpha = alpha, ns = d)$ScorePertienceCausale2
  
  evOptimal = evaluation(truth, pc.fit.optimal)
  
  
  if (display){
    nbrow = sqrt(nInits)
    if (floor(nbrow) == nbrow){nbrcol = nbrow} else nbrow = floor(nbrow) +1 ; nbrcol = nbrow
    par(mfrow = c(1,1))
    par(mar=c(5, 4, 4, 18), xpd=F)
    maxIt = lapply(Res_Causal_for_each_init,function(x){length(x$TP)})
    maxIt = max(unlist(maxIt))
    if (test == "pSCCI"){ylim = c(0.989,0.992)}
    else {ylim = c(0,1)}
    plot(x = NA, xlim= c(1,maxIt), ylim = ylim, ylab = " ", xlab = "iteration", main = paste(structure, 'm = ',m, ' ', test))
    
    for (initIndex in 1:nInits){
      Res_Causal = Res_Causal_for_each_init[[initIndex]]
      
      if (initIndex == ibestDisc){
        finalIt = length(Res_Causal$SCORE)
        points(finalIt, Res_Causal$SCORE[finalIt], pch = 1, col = 2, cex = 2)
      }
      
      TP = Res_Causal$TP
      TN = Res_Causal$TN
      
      
      tempGriser1 = which(TP == 1 & TN == 1 )
      tempGriser2 = split(tempGriser1, cumsum(c(1, diff(tempGriser1) != 1)))
      tempGriser3 = tempGriser2
      for (inter in 1:length(tempGriser2)){
        tempGriser3[[inter]] = range( tempGriser2[[inter]])
      }
      
      
      points(Res_Causal$SCORE, type = 'o', lwd = 1.5, col ="grey")
      
      
      
      Res_Causal$SCOREfusion
      
      
      Res_Causal$boundaries
      Res_Causal$itnewboundVar
      Res_Causal$itnewboundVal
      fusionits = which(!is.na(Res_Causal$itnewfusionVal))
      
      textVar = Res_Causal$itnewfusionVar[fusionits]
      textVal = Res_Causal$itnewfusionVal[fusionits]
      labels = paste0("var ", textVar, '\n',textVal)

    }
    

    
    points(Res_Causal_for_each_init[[ibestDisc]]$SCORE, type = 'o', lwd = 1.5)
    Res_Causal_for_each_init[[ibestDisc]]$TP[finalIt]
    Res_Causal_for_each_init[[ibestDisc]]$TN[finalIt]
    
    # adding info about the target optimal discretisation 
    
    abline(h = scoreOptimal, col = 1, lty = 2)
    
    par(xpd=TRUE)
    evOptimal
    
    
    
    legend('topright', inset=c(-2, 0), legend = c("crs trajectories", "selected crs trajectories", "final optimised discretization crs"), lty = c(1,1,NA), col = c(1,'grey',2), pch = c(1,1,1))
    text(x = maxIt + maxIt*1, y = 0.5, labels = paste0("optimized discretization TP rate:", Res_Causal_for_each_init[[ibestDisc]]$TP[finalIt], 
                                                       "\n optimized discretization TN rate:", Res_Causal_for_each_init[[ibestDisc]]$TN[finalIt],
                                                       "\n perfect discretization TP rate:", evOptimal$truePosSqlt, 
                                                       "\n perfect discretization TN rate:", evOptimal$trueNegSqlt))
    
  }
  
  
  
  
  RES = list(RESBestIt = Res_Causal_for_each_init[[ibestDisc]], scoreOptimised = scoreOptimised, evOptimised = evOptimised )
  return(RES)
  
}



