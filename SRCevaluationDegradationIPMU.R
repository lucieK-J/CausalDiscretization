
# Ici on realise une etude de la degradation du score pour IPMU

alpha = 0.01
n = 100
seed.degr = 74
nbaddBound = 60
N = 100
indexExample = 43

# N = indexExample
seed.degr = 1
structure = "mediator"; seed = 6; n=500 ; indexExample = 1# , 5# pas mal
# seed = 14 # à voir
# seed = 5
test = "Mutual.Information.test"
# test = "pSCCI"

# structure = "vstructure"; seed = 38; n = 100; nbaddBound = 60 # paper AISTAT
# structure = "vstructure"; seed = 2; n = 100; nbaddBound = 60; indexExample = 43 ; seed.degr = 74# bien pour IPMU mais pourquoi la courbe remonte ??? 
# verifier que je ne calcule plus l'index pour des variable contante.. 
# structure = "vstructure"; seed = 2; n = 100; nbaddBound = 60


# structure = "diamond";  seed = 2#  seed = 1

param = c(structure,n, seed, alpha)


source(file = "CausalDiscretization.R")
source(file = "SimulationIPMUGenerator.R")
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





SIM = simulationIPMU(structure, n, seed = seed, N = N)


SIM

SIM$RULES


discret <- function(DATAcont, SEP = SEPminus){
  d = dim(DATAcont[[1]])[2]   
  N = length(DATAcont)
  DATAdiscr = vector(mode = "list", length = length(DATAcont))
  names(DATAdiscr) = paste0("DATAdiscr", 1:length(DATAcont))
  
  for (dist in 1:N){
    DATAdiscr[[dist]] = data.frame(DATAcont[[dist]])
    
    for (v in 1:d){
      DATAdiscr[[dist]][,v] = cut(DATAcont[[dist]][,v], SEP[[v]])
    }
  }
  return(DATAdiscr)
}



degradationIPMU <- function(SIM, nbaddBound  = nbaddBound, alpha, seed.degr = seed.degr, param, nosave){
  # on se passe de la fonction experiences. 
  # on suppose quel que chose comme SIM = simulationIPMU(structure, n, seed, N = 100) 
  N = length(SIM$dm)
  set.seed(seed.degr)
  d = dim(SIM$ADJ)[1]
  ADJ = SIM$ADJ
  DATAcont = SIM$DATAcont

  SEP = SIM$SEP
  temp = unlist(SEP)
  nbbound = length(temp[!temp %in% c(0,1)])
  
  
  
  RANGEG = vector(mode = "list", length = nbaddBound + nbbound + 1)
  
  names(RANGEG) = c(paste0("D+",nbaddBound:1,"b"),"D*", paste0("D-", 1:nbbound,"b") )
  
  RANGEG[["D*"]] = list()
  RANGEG[["D*"]]$SEP = SIM$SEP
  RANGEG[["D*"]]$DATAdiscr = SIM$DATAdiscr
  
  SEPminus = SIM$SEP
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
    

    
    RANGEG[[nbaddBound+1+i]]$DATAdiscr = discret(DATAcont, SEP = SEPminus)
    for (dist in 1:length(DATAcont)){
    if(sum(is.na(RANGEG[[nbaddBound+1+i]]$DATAdiscr[[dist]])) > 0){
      stop()
    }
    }
    
    
    varAlnbLeft = length(SEPminus[[varAl]])
    if ( varAlnbLeft == 2){
      candVar = candVar[-which(candVar == varAl)]
    }
    
    
  }
  
  
  
  SEPplus = SIM$SEP
  
  
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
    
    for (dist in 1:N){
      RANGEG[[i]][[paste("dist",dist)]] = list()
      
    for (test in tests){
      RANGEG[[i]][[paste("dist",dist)]][[test]] = list()
      
      
      labels = SIM$NameVars
      nlevP = rep(NA, dim(RANGEG[[i]]$DATAdiscr[[dist]])[2])
      
      dm = RANGEG[[i]]$DATAdiscr[[dist]]
      for (v in 1:dim(RANGEG[[i]]$DATAdiscr[[dist]])[2]){
        dm[,v] = as.numeric(dm[,v]) - 1
        nlevP[v] = length(table(dm[,v]))
      } 
      
      suffStat = list(dm = dm, nlev  = nlevP, adaptDF = F)
      
      pc.fit <- mod.pc(suffStat = suffStat,
                       indepTest = get(test),
                       alpha = alpha, labels = labels, verbose = TRUE)
      
      # plot(pc.fit, main = " ")
      RANGEG[[i]][[paste("dist",dist)]][[test]]$plot = pc.fit
      
      ev = evaluation(ADJ, pc.fit)
      
      RANGEG[[i]][[paste("dist",dist)]][[test]]$ev = ev
      
      if(is.null((ev))){
        RANGEG[[i]][[paste("dist",dist)]][[test]]$ev = list(NA)
      }
      
      upc.fit <- evaluate::evaluate({
        function(){pc.fit <- mod.pc(suffStat = suffStat,
                                    indepTest = get(test),
                                    alpha = alpha, labels = labels, verbose = TRUE) }
      }, new_device = F)
      
      RANGEG[[i]][[paste("dist",dist)]][[test]]$spc = compute.score(u = upc.fit, alpha = alpha, ns = d)
      
    }
    }
    
  }
  
  for (test in tests){
    par(mar = c(5.1, 4.1, 4.1, 2.1))

    spc_RANGEG = matrix(NA, nrow = N, ncol = length(RANGEG) )  
    truePos_RANGEG = matrix(NA, nrow = N, ncol = length(RANGEG) )  
    trueNeg_RANGEG = matrix(NA, nrow = N, ncol = length(RANGEG) )  
    for (dist in 1:N){
      spc_RANGEG[dist,] <- unlist(lapply(RANGEG, function(x) {x[[paste0("dist ", dist)]][[test]]$spc$ScorePertienceCausale2}))
      truePos_RANGEG[dist,] = unlist(lapply(RANGEG, function(x) {x[[paste0("dist ", dist)]][[test]]$ev$truePosSqlt} ))
      trueNeg_RANGEG[dist,] = unlist(lapply(RANGEG, function(x) {x[[paste0("dist ", dist)]][[test]]$ev$trueNegSqlt} ))
    }
    
    colnames(spc_RANGEG) = names(RANGEG)
    rownames(spc_RANGEG) = paste0("dist ", 1:N)
    colnames(truePos_RANGEG) = names(RANGEG)
    rownames(truePos_RANGEG) = paste0("dist ", 1:N)
    colnames(trueNeg_RANGEG) = names(RANGEG)
    rownames(trueNeg_RANGEG) = paste0("dist ", 1:N)
     
# calcul des quantiles
    
    quantile(spc_RANGEG[,50], probs = c(0.05,0.5,0.95))
    Quantiles_spc_RANGEG = apply(spc_RANGEG,2, quantile, probs = c(0.05,0.5,0.95), na.rm = TRUE)
    Quantiles_truePos_RANGEG = apply(truePos_RANGEG, 2, quantile,probs = c(0.05,0.5,0.95), na.rm = TRUE)
    Quantiles_trueNeg_RANGEG = apply(trueNeg_RANGEG, 2, quantile,probs = c(0.05,0.5,0.95), na.rm = TRUE)
    
    
    str(spc_RANGEG)
    Ex_spc_RANGEG = spc_RANGEG[indexExample,]
    Ex_truePos_RANGEG = truePos_RANGEG[indexExample,]
    Ex_trueNeg_RANGEG = trueNeg_RANGEG[indexExample,]
    
    if (!nosave){
    # postscript(file = paste0("results/degradation/degradationSPC_IPMU",paste0(param, collapse = "_"), "_", test,".eps"), 
    #           pointsize = 18, width = 20, height = 5)
    ##### FLAG -----
    pdf(file =  paste0("results/degradation/degradationSPC_IPMU",paste0(param, collapse = "_"), "_", test,".eps"), pointsize = 18, width = 12, height = 5)
    }
    
    par(mar = c(5.1, 3.1, 1.1, 1.1))
    
    plot(1:length(RANGEG), Quantiles_spc_RANGEG[2,], xaxt = "n", yaxt = "n", type = 'l', ylim = c(0,1), xlab =" ", ylab = " ", lwd = 2) 
   
    maxNotNA = min(which(is.na(Quantiles_spc_RANGEG[3,])))-1
    polygon(c(1:maxNotNA, maxNotNA:1), c(Quantiles_spc_RANGEG[3,1:maxNotNA], rev(Quantiles_spc_RANGEG[1,1:maxNotNA])), 
            col = adjustcolor( 1, alpha.f = 0.1),  border = FALSE)
    
    lines(1:length(RANGEG), Quantiles_truePos_RANGEG[2,], lty = 2, col = 4, lwd = 2)
    polygon(c(1:maxNotNA, maxNotNA:1), 
            c(Quantiles_truePos_RANGEG[3,1:maxNotNA], rev(Quantiles_truePos_RANGEG[1,1:maxNotNA])), 
            col = adjustcolor( 4, alpha.f = 0.2),  border = FALSE)
    
    lines(1:length(RANGEG), Quantiles_trueNeg_RANGEG[2,], lty = 3, col = 2, lwd = 2)
    polygon(c(1:maxNotNA, maxNotNA:1), 
            c(Quantiles_trueNeg_RANGEG[3,1:maxNotNA], rev(Quantiles_trueNeg_RANGEG[1,1:maxNotNA])), 
            col = adjustcolor( 2, alpha.f = 0.2),  border = FALSE)
    
    axis(1, at =1:length(RANGEG),
         labels = names(RANGEG),las = 3,cex.axis=1)
    axis(2, at =c(0,0.25,0.5,0.75,1),cex.axis=1)
    
       
    
    abline(v = nbaddBound + 1, lwd = 2 )

    
    legend('topleft', legend = c("cri", "TP", "TN"), col = c(1,4,2), lty = c(1,2,3), lwd = 2, bg="white", cex=1)
    
    if (!nosave){
    dev.off()
    
   postscript(file = paste0("results/degradation/degradationSPC_IPMUoneExample",paste0(param, collapse = "_"), "_", test,".eps"), 
            pointsize = 18, width = 12, height = 5)
    # pdf(file =  paste0("results/degradation/degradationSPC_IPMUoneExample",paste0(param, collapse = "_"), "_", test,".eps"),
     #    pointsize = 18, width = 12, height = 5)
      
    }
    
    par(mar = c(5.1, 3.1, 1.1, 1.1))
    plot(1:length(RANGEG), Ex_spc_RANGEG, xaxt = "n", yaxt = "n", type = 'o', ylim = c(0,1), xlab =" ", ylab = " ", lwd = 1.5)
    
    lines(1:length(RANGEG), Ex_truePos_RANGEG, lty = 2, col = 4, lwd = 1.5)
    
    lines(1:length(RANGEG), Ex_trueNeg_RANGEG, lty = 3, col = 2, lwd = 1.5)
    
    axis(1, at =1:length(RANGEG),
         labels = names(RANGEG),las = 3, cex.axis=1)
    axis(2, at = c(0,0.25,0.5,0.75,1), cex.axis=1)
    
    
    
    abline(v = nbaddBound + 1 , lwd = 1.5)
    
    
    legend('topleft', legend = c("cri", "TP", "TN"), col = c(1,4,2), lty = c(1,2,3), lwd = 1, bg="white", cex=1)
    
    if (!nosave){
    dev.off()
    }
    
    
    
    
    Graphes_interessants = 1:length(RANGEG)
    #for (g in Graphes_interessants){
    #  plot(RANGEG[[g]][[test]]$plot, main = names(RANGEG)[g]) 
    #}
    
    
    
  }
  
  
  
  
  
  
}

degradationIPMU(SIM, nbaddBound  = nbaddBound, alpha, seed.degr = seed.degr, param, nosave = FALSE)





degradation.etude.edgeIPMU <- function(SIM, indexExample, nbaddBound = nbaddBound, alpha, seed.degr = seed.degr, test = test, param, nosave){
  
  N = length(SIM$dm)
  set.seed(seed.degr)
  d = dim(SIM$ADJ)[1]
  ADJ = SIM$ADJ
  DATAcont = SIM$DATAcont[[indexExample]]
  
  SEP = SIM$SEP
  temp = unlist(SEP)
  nbbound = length(temp[!temp %in% c(0,1)])
  
  
  
  RANGEG = vector(mode = "list", length = nbaddBound + nbbound + 1)
  
  names(RANGEG) = c(paste0("D+", nbaddBound:1,"b"),"D*", paste0("D-", 1:nbbound,"b") )
  
  RANGEG[["D*"]] = list()
  RANGEG[["D*"]]$SEP = SIM$SEP
  RANGEG[["D*"]]$DATAdiscr = SIM$DATAdiscr[[indexExample]]
  
  
  #set.seed(seed.degr)
  #d = dim(EXP$SIM$ADJ)[1]
  #ADJ = EXP$SIM$ADJ
  #DATAcont = EXP$SIM$DATAcont
  
  temp = unlist(SIM$SEP)
  nbbound = length(temp[!temp %in% c(0,1)])
  
  
  
  RANGEG = vector(mode = "list", length = nbaddBound + nbbound + 1)
  
  names(RANGEG) = c(paste0("D+",nbaddBound:1,"b"),"D*", paste0("D-", 1:nbbound,"b") )
  
  RANGEG[["D*"]] = list()
  RANGEG[["D*"]]$SEP = SIM$SEP
  RANGEG[["D*"]]$DATAdiscr = SIM$DATAdiscr[[indexExample]]
  
  SEPminus = SIM$SEP
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
  
  
  SEPplus = SIM$SEP
  
  
  for (i in 1:nbaddBound){
    
    varAl = sample(1:d,1)
    
    newBound = runif(1)
    SEPplus[[varAl]] = c(newBound, SEPplus[[varAl]])[order( c(newBound, SEPplus[[varAl]]))]
    
    RANGEG[[nbaddBound + 1 - i]] = list()
    
    RANGEG[[nbaddBound + 1 - i]]$SEP = SEPplus
    
    RANGEG[[nbaddBound + 1 - i]]$DATAdiscr = discret(DATAcont, SEP = SEPplus)
    
  }
  
  
  for (i in 1:length((RANGEG))){
    # RANGEG[[i]]$DATAdiscr
    
    RANGEG[[i]][[test]] = list()
    
    labels = SIM$NameVars
    nlev = rep(NA, dim(RANGEG[[i]]$DATAdiscr)[2])
    
    dm = RANGEG[[i]]$DATAdiscr
    for (v in 1:dim(RANGEG[[i]]$DATAdiscr)[2]){
      dm[,v] = as.numeric(dm[,v]) - 1
      nlev[v] = length(table(dm[,v]))
    } 
    
    suffStat = list(dm = dm, nlev = nlev, adaptDF = F)
    
    pc.fit <- mod.pc(suffStat = suffStat,
                     indepTest = get(test), 
                     alpha = alpha, labels = labels, verbose = TRUE)
    
    RANGEG[[i]][[test]]$plot = pc.fit
    
    if ( is.null(pc.fit) ){
      RANGEG[[i]][[test]]$plot = "no estimation"
    }
    
    ev = evaluation(ADJ, pc.fit)
    
    RANGEG[[i]][[test]]$ev = ev
    
    
    upc.fit <- evaluate::evaluate({
      function(){pc.fit <- mod.pc(suffStat = suffStat,
                                  indepTest = get(test), 
                                  alpha = alpha, labels = labels, verbose = TRUE) }
    }, new_device = F) 
    
    RANGEG[[i]][[test]]$spc = compute.score(u = upc.fit, alpha = alpha, ns = d)
    
    
    
  }
  
  
  spc_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$ScorePertienceCausale2} ))
  
  plot(1:length(RANGEG), spc_RANGEG, xaxt = "n" , yaxt = "n", type = 'o', ylim = c(0,1), 
       main = " ", xlab =" ", ylab = " ")
  
  lines(1:length(RANGEG), spc_RANGEG, xaxt = "n", type = 'o')
  axis(1, at = 1:length(RANGEG),
       labels = names(RANGEG), las = 3, cex.axis=1.5)
  axis(2, at =c(0,0.25,0.5,0.75,1), cex.axis=1.5)
  abline(v = nbaddBound + 1 )
  # abline(h = alpha)
  # par(xpd = T )
  # text(length(RANGEG) + 5 , alpha, expression(alpha), cex = 1.6)
  
  
  sumAsupp_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$sumAsupp} ))
  sumArest_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$spc$sumArest} ))
  
  truePos_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$ev$truePosSqlt} ))
  trueNeg_RANGEG = unlist(lapply(RANGEG, function(x) {x[[test]]$ev$trueNegSqlt} ))

  
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
    if (!nosave){
      
      #FAG 2 -----
    postscript(file = paste0("results/degradation/degradationEdegeWiseSPCedgeIPMU", edge, paste0(param, collapse = "_"), "_", test,".eps"), 
               pointsize = 18, width = 7, height = 7)
      
      # pdf(file =  paste0("results/degradation/degradationSPC_IPMU",paste0(param, collapse = "_"), "_", test,".eps"), pointsize = 18, width = 12, height = 5)
    }
    
  
  par(mar = c(5.1, 3.1, 1.1, 1.2))
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
    axis(2, at = c(0,0.25,0.5,0.75,1), cex.axis = 1.5)
    abline(v = nbaddBound + 1 )
    abline(h = alpha)
    par(xpd = T )
    text(length(RANGEG) + 5, alpha, expression(alpha), cex = 1.6)
    if (!nosave){
    dev.off()}
    
  }
  
  
  
}


seed.degr = 1 # 6, 15, 36, 44,47,  53, 57, 60, 62, 74, 81
indexExample = 43 # 27 (n0n) # 21#12# 11#7 # 43
degradation.etude.edgeIPMU(SIM, indexExample = indexExample, nbaddBound = nbaddBound, alpha, seed.degr = seed.degr, test = "Mutual.Information.test", param, nosave = FALSE)
  

