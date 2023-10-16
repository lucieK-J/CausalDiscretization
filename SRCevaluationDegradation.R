

alpha = 0.01
n = 1000

nbaddBound = 100


structure = "vstructure"; seed = 38; n = 100; nbaddBound = 60 # paper



param = c(structure,n, seed, alpha)


source(file = "CausalDiscretization.R")
source(file = "SimulationGenerator.R")
source(file = "usefulFunctions.R")
source(file = "tests_for_conditionnalInd.R")
require(arules)
require(stringr)
library(igraph)
library(RColorBrewer)
library(tidyr) # for extract_numeric
library(SCCI)
library(bnlearn)
library(sna)
library(caret)
library(entropy)





SIM = simulation(structure, n, seed)




experiences <- function(SIM, alpha){
  
  
  mapply(assign, names(SIM), SIM, MoreArgs = list(envir = globalenv()) )
  EXP = vector(mode = "list", 5)
  
  
  EXP$SIM = SIM
  
  ADJ = SIM$ADJ
  d = dim(SIM$DATAcont)[2]
  
  for (dataType in c("cont", "optimal", "pert", "pertS")){
   
    EXP[[dataType]] = list()
    
    if (dataType == "cont"){
      data = SIM$DATAcont
      nlev = NULL
      labels = SIM$NameVars
    }
    
    
    if (dataType == "optimal"){
      data = SIM$dm
      nlev =SIM$nlev
      labels = SIM$NameVars
    }
    
    
    if (dataType == "init"){
      Vtodiscr = NameVars
      esp = 0.001
      data = data.frame(DATAcont)
      colnames(data) = NameVars
      data_initale = causal.discretization(data, floor(n/100), Vtodiscr, esp, alpha, onlyInit = T)
      # save(file = "data_initale",data_initale)
      res_mattodm = mattodm(m = data_initale)
      data = res_mattodm$dm
      nlev = res_mattodm$nlev
      labels = SIM$NameVars
    }
    
    if (dataType == "pert"){
      DATA = SIM$DATA
      SEP = SIM$SEP
      
      ind = which(ADJ != 0, arr.ind = T)[1]
      
      SEPP = SEP
      
      if (!is.na(ind)){
        SEPP[[ind]] = SEP[[ind]] [-2]
        

        DATAdiscrP = data.frame(DATAcont)
        
        
        for (v in 1:dim(DATAcont)[2]){
          datacont =  unlist(DATA[[v]])
          
          DATAcont[,v] = datacont
          DATAdiscrP[,v] = cut(datacont, SEPP[[v]])
        }
        
        
        dmP = DATAdiscrP
        
        
        nlevP = rep(NA, dim(DATAcont)[2])
        for (v in 1:dim(DATAcont)[2]){
          dmP[,v] = as.numeric(dmP[,v]) - 1
          nlevP[v] = length(table(dmP[,v]))
        }
        
        data = dmP
        nlev = nlevP
      }}
    
    if (dataType == "pertS"){
      DATA = SIM$DATA
      SEP = SIM$SEP

      ind = which(ADJ != 0, arr.ind = T)[1]
      
      SEPPS = vector(mode = "list", d)
      names(SEPPS) = paste0("V",1:d)
      for (v in 1:length(SEPPS)){
        SEPPS[[v]] = c(0,1)
      }
      
      nbsep = rbinom(n = 1,size = d, prob = neT/d) 
      
      for (i in 1:nbsep){
        vars = base::sample(size = 2,1:d)
        for (j in vars){
          SEPPS[[j]] = c(SEPPS[[j]] , runif(1))
        }
      }
      
      for (v in 1:length(SEPPS)){
        SEPPS[[v]] = SEPPS[[v]][order(SEPPS[[v]] )]
      }
      

      DATAdiscrPS = data.frame(DATAcont)
      
      
      for (v in 1:dim(DATAcont)[2]){
        datacont =  unlist(DATA[[v]])
        
        DATAcont[,v] = datacont
        DATAdiscrPS[,v] = cut(datacont, SEPPS[[v]])
      }
      
      
      DATAdiscrPS
      
      dmPS = DATAdiscrPS
      
      
      nlevPS = rep(NA, dim(DATAcont)[2])
      for (v in 1:dim(DATAcont)[2]){
        dmPS[,v] = factor(dmPS[,v])
        dmPS[,v] = as.numeric(dmPS[,v]) - 1
        nlevPS[v] = length(table(dmPS[,v]))
      }
      
      data = dmPS
      nlev =nlevPS
    }
      
    tests = c("gaussCItest", "disCItest", "pSCCI", "Mutual.Information.test")
    for (test in tests){
      
      
      if(test == "gaussCItest"){
        suffStat = list(C = cor(data), n = n)
      }
      
      
      
      if (test %in% c("disCItest", "Jonckheere.Terpstra.test", "Mutual.Information.test","pSCCI")){
        suffStat = list(dm = data, nlev  = nlev, adaptDF = F)
        if (dataType == "cont"){
          next
        }
      }
      
      
      EXP[[dataType]][[test]] = list()
      
      pc.fit <- mod.pc(suffStat = suffStat,
                       indepTest = get(test),
                       alpha = alpha, labels = labels, verbose = TRUE)

      
      EXP[[dataType]][[test]]$plot = pc.fit
      
      ev = evaluation(ADJ, pc.fit)
      EXP[[dataType]][[test]]$ev = ev
      
      
      upc.fit <- evaluate::evaluate({
        function(){pc.fit <- mod.pc(suffStat = suffStat,
                                    indepTest = get(test), ## indep.test: partial correlations
                                    alpha = alpha, labels = labels, verbose = TRUE) }
      }, new_device = F)
      
      
      EXP[[dataType]][[test]]$spc = compute.score(u = upc.fit, alpha = alpha, ns = d)
    }
    
  }
  
  
  
  
  return(EXP)
  
}

EXP = experiences(SIM, alpha)

degradation <- function(EXP, nbaddBound  = 10, alpha, seed.degr = 1, param){
  
  set.seed(seed.degr)
  d = dim(EXP$SIM$ADJ)[1]
  ADJ = EXP$SIM$ADJ
  DATAcont = EXP$SIM$DATAcont

  
  temp = unlist(EXP$SIM$SEP)
  nbbound = length(temp[!temp %in% c(0,1)])
  
  
  
  RANGEG = vector(mode = "list", length = nbaddBound + nbbound + 1)
  
  names(RANGEG) = c(paste0("D+",nbaddBound:1,"b"),"D*", paste0("D-", 1:nbbound,"b") )
  
  RANGEG[["D*"]] = list()
  RANGEG[["D*"]]$SEP = EXP$SIM$SEP
  RANGEG[["D*"]]$DATAdiscr = EXP$SIM$DATAdiscr
  
  SEPminus = EXP$SIM$SEP
  candVar = 1:d
  for ( j in 1:d){
    if (length(SEP[[j]])==2){
      candVar = candVar[-(which(candVar) == j)]
    }
  }
  
  for (i in 1:nbbound){
    
    varAl = sample(candVar,1)
    
    if( length(candVar) == 1){varAl = candVar }
    
    interAl = sample(2:(length(SEPminus[[varAl]])-1),1)
    
    
    if(length(SEPminus[[varAl]]) == 3){
      interAl = 2
    } 
    
    SEPminus[[varAl]] = SEPminus[[varAl]][-interAl]
    
    RANGEG[[nbaddBound +1+i]] = list()
    
    RANGEG[[nbaddBound +1+i]]$SEP = SEPminus
    
    discret <- function(DATAcont, SEP = SEPminus){
      d = dim(DATAcont)[2]   
      DATAdiscr = data.frame(DATAcont)
      
      for (v in 1:d){
        DATAdiscr[,v] = cut(DATAcont[,v], SEP[[v]])
      }
      return(DATAdiscr)
    }
    
    RANGEG[[nbaddBound+1+i]]$DATAdiscr = discret(DATAcont, SEP = SEPminus)
    if(sum(is.na(RANGEG[[nbaddBound+1+i]]$DATAdiscr)) > 0){
      stop()
    }
    
    
    varAlnbLeft = length(SEPminus[[varAl]])
    if ( varAlnbLeft == 2){
      candVar = candVar[-which(candVar == varAl)]
    }
    
    
  }
  
  
  
  SEPplus = EXP$SIM$SEP
  
  
  for (i in 1:nbaddBound){
    
    varAl = sample(1:d,1)
    
    newBound = runif(1)
    SEPplus[[varAl]] = c(newBound, SEPplus[[varAl]])[order( c(newBound, SEPplus[[varAl]]))]
    
    RANGEG[[nbaddBound + 1 - i]] = list()
    
    RANGEG[[nbaddBound + 1 - i]]$SEP = SEPplus
    
    
    RANGEG[[nbaddBound + 1 - i]]$DATAdiscr = discret(DATAcont, SEP = SEPplus)
    
  }
  
  
  for (i in 1:length(RANGEG)){
    RANGEG[[i]]$DATAdiscr
    
    tests = c("disCItest", "pSCCI", "Mutual.Information.test")
    
    for (test in tests){
      RANGEG[[i]][[test]] = list()
      
      labels = EXP$SIM$NameVars
      nlevP = rep(NA, dim(RANGEG[[i]]$DATAdiscr)[2])
      
      dm = RANGEG[[i]]$DATAdiscr
      for (v in 1:dim(RANGEG[[i]]$DATAdiscr)[2]){
        dm[,v] = as.numeric(dm[,v]) - 1
        nlevP[v] = length(table(dm[,v]))
      } 
      
      suffStat = list(dm = dm, nlev  = nlev, adaptDF = F)
      
      pc.fit <- mod.pc(suffStat = suffStat,
                       indepTest = get(test),
                       alpha = alpha, labels = labels, verbose = TRUE)
      
      # plot(pc.fit, main = " ")
      RANGEG[[i]][[test]]$plot = pc.fit
      
      ev = evaluation(ADJ, pc.fit)
      
      RANGEG[[i]][[test]]$ev = ev
      if(is.null((ev))){
        RANGEG[[i]][[test]]$ev = list(NA)
      }
      
      upc.fit <- evaluate::evaluate({
        function(){pc.fit <- mod.pc(suffStat = suffStat,
                                    indepTest = get(test),
                                    alpha = alpha, labels = labels, verbose = TRUE) }
      }, new_device = F)
      
      RANGEG[[i]][[test]]$spc = compute.score(u = upc.fit, alpha = alpha, ns = d)
      
    }
    
    
  }
  
  for (test in tests){
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    spc_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$ScorePertienceCausale2} ))
    sumAsupp_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$sumAsupp} ))
    sumArest_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$sumArest} ))
    
    truePos_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$ev$truePosSqlt} ))
    trueNeg_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$ev$trueNegSqlt} ))
    
    tempGriser1 = which(truePos_RANGEG == 1 & trueNeg_RANGEG == 1 )
    tempGriser2 = split(tempGriser1, cumsum(c(1, diff(tempGriser1) != 1)))
    tempGriser3 = tempGriser2
    for (inter in 1:length(tempGriser2)){
      tempGriser3[[inter]] = range( tempGriser2[[inter]])
    }
    
    
    postscript(file = paste0("results/degradation/degradationSPC",paste0(param, collapse = "_"), "_", test,".eps"), 
               pointsize = 18, width = 20, height = 5)
    
    
    plot(1:length(RANGEG), spc_RANGEG, xaxt = "n", yaxt = "n", type = 'o', ylim = c(0,1), xlab =" ", ylab = " ") 
   
    axis(1, at =1:length(RANGEG),
         labels = names(RANGEG),las = 3,cex.axis=0.5)
    axis(2, at =c(0,0.25,0.5,0.75,1),cex.axis=0.5)
    
       lines(1:length(RANGEG), spc_RANGEG, type = 'o')
    
    abline(v = nbaddBound + 1 )

    lines(1:length(RANGEG),truePos_RANGEG, col = 4, lwd = 1.5)
    lines(1:length(RANGEG),trueNeg_RANGEG, col = 2, lwd = 1.5, lty = 2)
    
    legend('topleft', legend = c("crs", "TP", "TN"), col = c(1,4,2), lty = c(1,1,2), pch = c(1,NA,NA), lwd = c(1,1.5, 1.5), bg="white", cex=0.5)
    
    dev.off()
    
    
    
    Graphes_interessants = 1:length(RANGEG)
    #for (g in Graphes_interessants){
    #  plot(RANGEG[[g]][[test]]$plot, main = names(RANGEG)[g]) 
    #}
    
    
    
  }
  
  
  
}

degradation(EXP, nbaddBound  = nbaddBound, alpha, seed.degr = 1, param)



degradation.etude.edge <- function(EXP, nbaddBound  = 10, alpha, seed.degr = 1, test = "Mutual.Information.test", param){
  
  set.seed(seed.degr)
  d = dim(EXP$SIM$ADJ)[1]
  ADJ = EXP$SIM$ADJ
  DATAcont = EXP$SIM$DATAcont
 
  temp = unlist(EXP$SIM$SEP)
  nbbound = length(temp[!temp %in% c(0,1)])
  
  
  
  RANGEG = vector(mode = "list", length = nbaddBound + nbbound + 1)
  
  names(RANGEG) = c(paste0("D+",nbaddBound:1,"b"),"D*", paste0("D-", 1:nbbound,"b") )
  
  RANGEG[["D*"]] = list()
  RANGEG[["D*"]]$SEP = EXP$SIM$SEP
  RANGEG[["D*"]]$DATAdiscr = EXP$SIM$DATAdiscr
  
  SEPminus = EXP$SIM$SEP
  candVar = 1:d
  for ( j in 1:d){
    if (length(SEP[[j]])==2){
      candVar = candVar[-(which(candVar) == j)]
    }
  }
  
  for (i in 1:nbbound){
    
    varAl = sample(candVar,1)
    
    if( length(candVar) == 1){varAl = candVar }
    
    interAl = sample(2:(length(SEPminus[[varAl]])-1),1)
  
    
    if(length(SEPminus[[varAl]]) == 3){
      interAl = 2
    } 
    
    SEPminus[[varAl]] = SEPminus[[varAl]][-interAl]
    
    RANGEG[[nbaddBound +1+i]] = list()
    
    RANGEG[[nbaddBound +1+i]]$SEP = SEPminus
    
    discret <- function(DATAcont, SEP = SEPminus){
      d = dim(DATAcont)[2]   
      DATAdiscr = data.frame(DATAcont)
      
      for (v in 1:d){
        DATAdiscr[,v] = cut(DATAcont[,v], SEP[[v]])
      }
      return(DATAdiscr)
    }
    
    RANGEG[[nbaddBound+1+i]]$DATAdiscr = discret(DATAcont, SEP = SEPminus)
    if(sum(is.na(RANGEG[[nbaddBound+1+i]]$DATAdiscr)) > 0){
      stop()
    }
    
    
    varAlnbLeft = length(SEPminus[[varAl]])
    if ( varAlnbLeft == 2){
      candVar = candVar[-which(candVar == varAl)]
    }
    
    
  }
  
  
  SEPplus = EXP$SIM$SEP
  
  
  for (i in 1:nbaddBound){
    
    varAl = sample(1:d,1)
    
    newBound = runif(1)
    SEPplus[[varAl]] = c(newBound, SEPplus[[varAl]])[order( c(newBound, SEPplus[[varAl]]))]
    
    RANGEG[[nbaddBound + 1 - i]] = list()
    
    RANGEG[[nbaddBound + 1 - i]]$SEP = SEPplus
    
    
    RANGEG[[nbaddBound + 1 - i]]$DATAdiscr = discret(DATAcont, SEP = SEPplus)
    
  }

  
  for (i in 1:length(RANGEG)){
    RANGEG[[i]]$DATAdiscr

    RANGEG[[i]][[test]] = list()
    
    labels = EXP$SIM$NameVars
    nlevP = rep(NA, dim(RANGEG[[i]]$DATAdiscr)[2])
    
    dm = RANGEG[[i]]$DATAdiscr
    for (v in 1:dim(RANGEG[[i]]$DATAdiscr)[2]){
      dm[,v] = as.numeric(dm[,v]) - 1
      nlevP[v] = length(table(dm[,v]))
    } 
    
    suffStat = list(dm = dm, nlev  = nlev, adaptDF = F)
    
    pc.fit <- mod.pc(suffStat = suffStat,
                     indepTest = get(test), 
                     alpha = alpha, labels = labels, verbose = TRUE)
    
    RANGEG[[i]][[test]]$plot = pc.fit
    
    ev = evaluation(ADJ, pc.fit)
    
    RANGEG[[i]][[test]]$ev = ev
    
    
    upc.fit <- evaluate::evaluate({
      function(){pc.fit <- mod.pc(suffStat = suffStat,
                                  indepTest = get(test), 
                                  alpha = alpha, labels = labels, verbose = TRUE) }
    }, new_device = F) 
    
    RANGEG[[i]][[test]]$spc = compute.score(u = upc.fit, alpha = alpha, ns = d)
    
    
    
  } 
  
  
  spc_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$ScorePertienceCausale} ))
  sumAsupp_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$sumAsupp} ))
  sumArest_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$sumArest} ))
  
  truePos_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$ev$truePosSqlt} ))
  trueNeg_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$ev$trueNegSqlt} ))
  
  tempGriser1 = which(truePos_RANGEG == 1 & trueNeg_RANGEG == 1 )
  tempGriser2 = split(tempGriser1, cumsum(c(1, diff(tempGriser1) != 1)))
  tempGriser3 = tempGriser2
  
  for (inter in 1:length(tempGriser2)){
    tempGriser3[[inter]] = range( tempGriser2[[inter]])
  }
  
  plot.area <- function(tempGriser3, ylim = c(-1,3)){
    for (inter in 1:length(tempGriser3)){
      polygon(x =  rep(tempGriser3[[inter]],each = 2),                           # X-Coordinates of polygon
              y = c(ylim[1], ylim[2], ylim[2], ylim[1]),                             # Y-Coordinates of polygon
              col = "seashell", border = NA)
      if (tempGriser3[[inter]][1] == tempGriser3[[inter]][2]){
        abline(v = tempGriser3[[inter]][1], col= "seashell", lwd = 2 )
      }
    }
  }
  

  
  possible_edges = which(upper.tri(ADJ), arr.ind = T)
  EDGES = vector(mode = 'list', dim(possible_edges)[1])
  names(EDGES) = paste0("V",possible_edges[,1],"-V",possible_edges[,2])
  
  for (edge in 1:dim(possible_edges)[1]){
    ind_edge = possible_edges[edge,]
    EDGES[[edge]] = list()
    EDGES[[edge]]$pval = rep(NA, length(RANGEG))
    for (i in 1:length(RANGEG)){
      pvals = RANGEG[[i]][[test]]$spc$pvals
      
      if ( typeof(pvals) == "double" ){
        pvals2 = reformate2(ADJ = ADJ, 
                            sortiePC = RANGEG[[i]][[test]]$plot,
                            autreMat = pvals)$Mreformat
        
        
        EDGES[[edge]]$pval[i] = pvals2[ind_edge[2], ind_edge[1]]
      }
      
    }
    
  }
  
 
  
  for (edge in 1:dim(possible_edges)[1]){
    postscript(file = paste0("results/degradation/degradationEdegeWiseSPCedge",edge,paste0(param, collapse = "_"), "_", test,".eps"), 
               pointsize = 18, width = 15, height = 7)
    pval_RANGEG =  EDGES[[edge]]$pval
    statut = "truly absent"
    if ( ADJ[possible_edges[edge,1],possible_edges[edge,2] ] == 1){
      statut = "truly present"
    }
    
    if ( ADJ[possible_edges[edge,2] , possible_edges[edge,1]] == 1){
      statut = "truly present"
    }
    
    
    plot(1:length(RANGEG), pval_RANGEG, xaxt = "n" , yaxt = "n", type = 'o', ylim = c(0,1), 
         main = " ", xlab =" ", ylab = " ")
        
    lines(1:length(RANGEG), pval_RANGEG, xaxt = "n", type = 'o')
    axis(1, at = 1:length(RANGEG),
         labels = names(RANGEG), las = 3, cex.axis=1.5)
    axis(2, at =c(0,0.25,0.5,0.75,1), cex.axis=1.5)
    abline(v = nbaddBound + 1 )
    abline(h = alpha)
    par(xpd = T )
    text(length(RANGEG) +5 , alpha, expression(alpha), cex = 1.6)
    dev.off()
  }

  
  
}


degradation.etude.edge(EXP, nbaddBound  = nbaddBound, alpha, seed.degr = 1, test = "Mutual.Information.test", param)

