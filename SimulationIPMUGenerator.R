

library(pcalg)



mattodm <- function(m){
  dm = m
  d = dim(m)[2]
  nlev = rep(NA,d)
  for (v in 1:d){
    dm[,v] = as.numeric(dm[,v]) - 1
    nlev[v] = length(table(dm[,v]))
  }
  return(list(dm = dm, nlev = nlev))
  
}



simulationIPMU <- function(structure, n, seed, N){
  # ce qui change c'est que j'ai maintenant N disctribtion a generer. et non plus une seule
  # sinon tout parail
  set.seed(seed)
  
  if (structure == "fork"){
    
    V <- c("V1", "V2", "V3")
    edL <- vector("list",length=3)
    names(edL) <- V
    edL[[1]] <- list(edges=c(2,3),weights=c(1,1))
    g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
    
  }
  
  
  if (structure == "vstructure"){
    

    V <- c( "V1", "V2", "V3")
    edL <- vector("list",length=3)
    names(edL) <- V
    edL[[1]] <- list(edges=c(),weights=c())
    edL[[2]] <- list(edges=c(1),weights=c(1))
    edL[[3]] <- list(edges=c(1),weights=c(1))
    g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
    
    
  }
  
  if (structure == "diamond"){
    

    V <- c("V1", "V2", "V3","V4")
    edL <- vector("list",length=4)
    names(edL) <- V
    edL[[1]] <- list(edges=c(2,3),weights=c(1,1))
    edL[[2]] <- list(edges=c(4),weights=c(1))
    edL[[3]] <- list(edges=c(4),weights=c(1))
    g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
    
  }
  
  if (structure == "mediator"){
    V <- c("V1", "V2", "V3")
    edL <- vector("list",length=3)
    names(edL) <- V
    edL[[1]] <- list(edges=c(2,3),weights=c(1,1))
    edL[[2]] <- list(edges=c(3),weights=c(1))
    
    g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
  }
  
  
  
  
  if (structure == "chain"){
    

    V <- c("V1", "V2", "V3")
    edL <- vector("list",length=3)
    names(edL) <- V
    edL[[1]] <- list(edges=c(2),weights=c(1))
    edL[[2]] <- list(edges=c(3),weights=c(1))
    
    g <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")
    
  }
  
  
  myDAG = g
  d = dim(as(myDAG,"matrix"))[1]
  
  
  DATA = vector(mode = 'list', N)
  names(DATA) = paste0("dist", 1:N)
  
  RULES = rep()
  
  SEP = vector(mode = "list", d)
  NameVars = paste0("V", 1:d)
  
  # names(DATA) = NameVars
  
  for (j in 1:N){
    DATA[[j]] = vector(mode = 'list', d)
    names(DATA[[j]]) = paste0("V", 1:d)
  for (i in 1:d){
    DATA[[j]][[i]] = as.list(rep(NA, n))
  }}
  
  names(SEP) = NameVars
  
  
   
  ADJ = as(myDAG,"matrix")
  
  # SUPPORT = ADJ
  SUPPORT = vector(mode = 'list', N)
  names(SUPPORT) = paste0("support dist", 1:N)
  for (j in 1:N){
    SUPPORT[[j]] = ADJ
  }
  
  
  ADJstate = matrix(0, nrow = dim(ADJ)[1], ncol = dim(ADJ)[2])
  
  neT = sum(ADJ)
  
  level1 = which(apply(ADJ, 2, sum) == 0)
  
  for (i in 1:d){
    SEP[[i]] = c( 0 , 1)
  }
  
  
  
  
  
  
  for (l in level1){
    for (j in 1:N){
    for (i in 1:n){
      DATA[[j]][[l]][[i]] =  runif(n = 1, min = 0, max = 1)
    }
    }
  }
  

  for (l in level1){
    for (j in 1:N){
    DATA[[j]][[l]] = DATA[[j]][[l]][sample(1:n)]
    }
  }
  
  
  meanNA <- function(x){
    if (length(x) == 1){
      if (is.na(x)){
        return(NA)
      }
    }
    return(mean(x, na.rm = T))
  }
  
  randomPick <- function(x){
    if (length(x) == 1){
      if (is.na(x)){
        return(NA)
      }
      
      return(x)
    }
    x = x[-which(is.na(x))]
    x = sample(x,size = 1)
    return(x)
  }
  
  presenceNA = T
  
  NodeToTreat = "init" 
  
  
  while (!is.null(NodeToTreat) | presenceNA){
    
    
    nodesC = c()
    for (j in 1:d){
      temp = c()
      for (dist in 1:N){
      temp = c(temp, unlist(lapply(DATA[[dist]][[j]], randomPick)))
      }
      
      if ( sum(!is.finite(temp)) == 0 ){
        nodesC = c(nodesC, j)
      }
    }
    
    
    for (node in nodesC){
      effets = which(ADJ[node,] == 1)
      
      effetNonTraite = effets[!(effets %in% which(ADJstate[node,] == 1))]
      if (length(effetNonTraite) == 0){next}
      for (effet in effetNonTraite){
        
        tempa = runif(2)
        Ia = tempa[order(tempa)]
        SEP[[node]] = c(SEP[[node]], Ia )
        
        
        tempb = runif(2)
        Ib = tempb[order(tempb)]
        Ib
        SEP[[effet]] = c(SEP[[effet]], Ib )
        
        for (dist in 1:N){
        DATAnode = unlist(lapply(DATA[[dist]][[node]],randomPick))
        
        index = which(DATAnode >= Ia[1] & DATAnode < Ia[2])
        for (ind in index){
          DATA[[dist]][[effet]][[ind]] = c(DATA[[dist]][[effet]][[ind]], runif(1, min = Ib[1], max = Ib[2]))
        }
        SUPPORT[[dist]][node, effet] = length(index)
        }
        
        ADJstate[node, effet] = 1
        RULES = c(RULES, paste0("IF X",node, " in [", round(Ia[1],3), " , ", round(Ia[2],3), " ] THEN X",effet, " in [", round(Ib[1],3), " , ", round(Ib[2],3), "] (support = ", length(index), ")" ))
        
        
      }
    }
    
    
     
    NodeToTreat = rep()
    presenceNA = F
    for (node in 1:d){
      if (sum(ADJ[,node] == ADJstate[,node]) == d){ # si tout les effets ont ?t? trait?s ? 
        
        for (dist in 1:N){
        DATA[[dist]][[node]] = lapply(DATA[[dist]][[node]], randomPick)
        
        if (sum(!is.finite(unlist(lapply(DATA[[dist]][[node]], randomPick)))) != 0){ # si il reste des NAs
          if (!(node %in% NodeToTreat)){
          NodeToTreat = c(NodeToTreat, node)}
          presenceNA = T
        }
        }}
    }
    
      
    for (node in NodeToTreat){
      
      for (dist in 1:N){
      DATA[[dist]][[node]] = lapply(DATA[[dist]][[node]], randomPick)
      realTogenerate = which(!is.finite(unlist(DATA[[dist]][[node]])))
      
        
      z = length(realTogenerate)
      
      
      leftones = which(is.na(unlist(DATA[[dist]][[node]])))
      
      
      DATA[[dist]][[node]][leftones] = runif(length(leftones), min = 0, 1)
      }
    }
    
    for( i in 1:d ){
      temp = SEP[[i]]
      SEP[[i]] = temp[order(temp)]
    }
    
  }
  
  
  DATAcont = vector(mode = 'list', N)
  names(DATAcont) = paste0("DATAcontdist", 1:N)
  
  DATAdiscr = vector(mode = 'list', N)
  names(DATAdiscr) = paste0("DATAdiscrdist", 1:N)
  
  mattodm_res =  vector(mode = 'list', N)
  names(mattodm_res) = paste0("mattodm_resdist", 1:N)
  
  dm =  vector(mode = 'list', N)
  names(dm) = paste0("dm_resdist", 1:N)
  
  for (dist in 1:N){
  DATAcont[[dist]] = matrix(nrow = n, ncol = d)
  DATAdiscr[[dist]] = data.frame(DATAcont[[dist]])
  
  for (v in 1:d){
    datacont = unlist(DATA[[dist]][[v]])
    if( length(datacont) != length(DATAcont[[dist]][,v]) ){
      return(list(DATAcont = DATAcont, v= v ,datacont = datacont ))
    }
    
    DATAcont[[dist]][,v] = datacont
    DATAdiscr[[dist]][,v] = cut(datacont, SEP[[v]])
  }

  mattodm_res[[dist]] = mattodm(m = DATAdiscr[[dist]])
  dm[[dist]] = mattodm_res[[dist]]$dm 
  
}
  
  nlev = mattodm_res[[dist]]$nlev # ca c'est la meme info pour toutes les distribution
   
  return(list(DATAcont = DATAcont, DATAdiscr = DATAdiscr,
              SEP = SEP, ADJ = ADJ, myDAG = myDAG, NameVars = NameVars, DATA = DATA, dm = dm, nlev = nlev,
              RULES = RULES, neT = neT, SUPPORT = SUPPORT))
}

