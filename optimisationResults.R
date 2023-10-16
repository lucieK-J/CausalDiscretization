



source(file="tests_for_conditionnalInd.R")
source(file = "SimulationGenerator.R")
source("CausalRelivanceScore.R")
source(file = "usefulFunctions.R")
source(file = "CausalDiscretization.R")


n = 200

m = 10

structures = c( "chain", "fork", "vstructure", "mediator", "diamond" )

alpha = 0.01

nInits = 30

N = 500 # number of simulations

esp = 0.0001

tests = c("gaussCItest", "kci.test", "pSCCI", "Mutual.Information.test")

discretTests = c( "pSCCI", "Mutual.Information.test")




experiences <- function(SIM, alpha, nInits, m, esp, tests, display){
  
  
  invisible(mapply(assign, names(SIM), SIM, MoreArgs = list(envir = globalenv())))
  
  EXP = vector(mode = "list", 5)
  
  
  EXP$SIM = SIM
  
  ADJ = SIM$ADJ
  d = dim(SIM$DATAcont)[2]
  
  dataTypes = c('eqfreq', "cont",  "optimal", 'optimised')
  for (dataType in dataTypes){
    progressInDataTypes = (which(dataType == dataTypes)-1)/length(dataTypes)
    print(paste0("start now ", dataType, " (progress is ", progressInDataTypes, ")"))
    EXP[[dataType]] = list()
    

    
    if (dataType == "cont"){
      data = SIM$DATAcont
      nlev = NULL
      labels = SIM$NameVars
    }
    
    
    if (dataType == "optimised"){
      data = SIM$DATAcont
      nlev = NULL
      labels = SIM$NameVars
      
      data = SIM$DATAcont
      Vtodiscr = NameVars
      # esp = 0.001
      data = data.frame(DATAcont)
      colnames(data) = NameVars
      
      
    }
    
    
    if (dataType == "eqfreq"){
      data = SIM$DATAcont 
      for (var in 1:dim(SIM$DATAcont)[2]){
        data[,var] = as.numeric(arules::discretize(data[,var], method = "frequency", breaks = SIM$nlev[var])) - 1
      }
      nlev = SIM$nlev
      labels = SIM$NameVars
    }
    
    
    
    
    if (dataType == "optimal"){
      data = SIM$dm
      nlev = SIM$nlev
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
      nlev = nlevPS
    }
    
      for (test in tests){
      startingTime = proc.time()
      
      print(paste0("-------- test is ", test ))
      
      if(test == "gaussCItest"){
        if(dataType %in% c("optimised", "eqfreq",  "optimal")){next}
        suffStat = list(C = cor(data), n = n)
      }
      
      if(test == "kci.test"){
        if(dataType %in% c("optimised", "eqfreq",  "optimal")){next}
        suffStat = list(C = cor(data), n = n, dm = data)
      }
      
      if (test %in% c("disCItest", "Jonckheere.Terpstra.test", "Mutual.Information.test","pSCCI")){
        suffStat = list(dm = data, nlev  = nlev, adaptDF = F)
        if (dataType == "cont"){
          next
        }
      }
      
      
      EXP[[dataType]][[test]] = list()
      
      if (dataType == "optimised"){
        
        ResOptimisation = causal.discretization.with.several.inits(nInits, data, m, Vtodiscr, esp, alpha, onlyInit = F, criteria, test, truth, seed = NULL, randomInit = T, display = display, SIM)
        EXP[[dataType]][[test]]$ev = ResOptimisation$evOptimised
        EXP[[dataType]][[test]]$spc =  ResOptimisation$scoreOptimised
        endingTime = proc.time()
        EXP[[dataType]][[test]]$ev$time = endingTime - startingTime
        
      } else {
      
      
      pc.fit <- mod.pc(suffStat = suffStat,
                       indepTest = get(test), 
                       alpha = alpha, labels = labels, verbose = FALSE)
      
      endingTime = proc.time()

      
      EXP[[dataType]][[test]]$plot = pc.fit
      
      ev = evaluation(ADJ, pc.fit)
      EXP[[dataType]][[test]]$ev = ev
      EXP[[dataType]][[test]]$ev$time = endingTime - startingTime
      
      invisible(upc.fit <- evaluate::evaluate({
        function(){pc.fit <- mod.pc(suffStat = suffStat,
                                    indepTest = get(test), 
                                    alpha = alpha, labels = labels, verbose = TRUE) }
      }, new_device = F))
      
      
      EXP[[dataType]][[test]]$spc = compute.score(u = upc.fit, alpha = alpha, ns = d)
      }
    }
    
  }
  
  return(EXP)
  
}



generation.N.sim <-function(N, param, tests, display){
  SIM1..N = vector(mode = "list", N)
  n = as.numeric(param[2])
  structure = param[1]
  alpha = as.numeric(param[3])
  nInits =  as.numeric(param[4])
  m = as.numeric(param[5])
  esp = as.numeric(param[6])
  for (seed in 1:N){
    print(paste0('Iteration ', seed))
    set.seed(seed)
    SIM =  simulation(structure, n, seed)
    
    k=1
    while(sum(SIM$nlev == 1) != 0){ SIM = simulation(structure, n, seed+10000*k) ; k=k+1}
    
    
    EXP = experiences(SIM, alpha, nInits, m , esp, tests, display)

    SIM1..N[[seed]] = EXP
  }
  
  save(SIM1..N, file = paste0("results/comparaison/", N,"sim_", paste(param, collapse = "")))
}






library(xtable)
tableauResum <- function(N, param, tests){
  file = paste0("results/comparaison/", N,"sim_", paste(param, collapse = ""))
  load(file)
  
  LIST_RES = list(param = list(), results = list())
  LIST_RES$param = param
 
  
  LIST_RES$res = list(baselines = list(), 
                      eqFreq = list(), 
                      causal_discretization = list(),
                      optimal = list())
  
  
  
  measuresTodisplay = c("trueNegSqlt", "truePosSqlt", "time")
  measure ="trueNegSqlt"

  for (test in c('gaussCItest', 'kci.test' )){
  for (measure in measuresTodisplay){
    LIST_RES$res$baselines[[test]][[measure]] =  mean(unlist(lapply(SIM1..N[1:N], function(x) {x[['cont']][[test]]$ev[[measure]]})))
  }}

  
  # causal_discretization 
  for (test in discretTests){
  for (measure in measuresTodisplay){
    LIST_RES$res$causal_discretization[[test]][[measure]] =  mean(unlist(lapply(SIM1..N[1:N], function(x) {x[['optimised']][[test]]$ev[[measure]]})))
  }
  }
  SIM1..N[[1]]$optimised$disCItest$ev$truePosSqlt
  SIM1..N[[1]]$optimal$disCItest$ev$truePosSqlt
  
  # optimal 
  for (test in discretTests){
    for (measure in measuresTodisplay){
      LIST_RES$res$optimal[[test]][[measure]] =  mean(unlist(lapply(SIM1..N[1:N], function(x) {x[['optimal']][[test]]$ev[[measure]]})))
    }
  }
  
  # eqFreq
  measure = "truePosSqlt"
  test = "pSCCI"
  for (test in discretTests){
    for (measure in measuresTodisplay){
      LIST_RES$res$eqFreq[[test]][[measure]] =  mean(unlist(lapply(SIM1..N[1:N], function(x) {x$eqfreq[[test]]$ev[[measure]]})))
    }
  }

  
  nblines = 2 + length(discretTests)
  nbcol = 3
  linesNames = rep()
  colNames = c('TP', 'TN', 'time')
  lev1 = names(LIST_RES$res)
  niceResults = c()
  
  for (i in lev1){
    linesNames = c(linesNames, paste0(i, ' ', names(LIST_RES$res[[i]])  )  )
    for (test in names(LIST_RES$res[[i]]) ){
      time = LIST_RES$res[[i]][[test]]$time
      TP = LIST_RES$res[[i]][[test]]$truePosSqlt
      TN = LIST_RES$res[[i]][[test]]$trueNegSqlt
      if (!is.null(c(TP, TN, time))){
        # print(c(TP, TN, time))
        niceResults = rbind(niceResults, c(TP, TN, time))
      } else niceResults = rbind(niceResults, c(NA, NA, NA))
      
    }
  }
  
  colnames(niceResults) = colNames
  rownames(niceResults) = linesNames
  
  niceResults
  
  xtable(niceResults)

  
  return(niceResults)

}

# _________________________________________________ GENERATION START


for (structure in structures){
  param = c(structure,n, alpha, nInits, m, esp) # attention : order is important
  generation.N.sim(N, param, tests, display = F)
  
  niceResults = tableauResum(N, param, tests)
  file 
  write.csv(file = paste0('results/comparaison/Results', N, paste0(param, collapse = "_") ), x = niceResults)
}

# _________________________________________________ GENERATION END



# _________________________________________________ READING RESUTS

niceResultsCombined = matrix()
Times = c()
k = 1

for (structure in structures){
  param = c(structure,n, alpha, nInits, m, esp)  # attention : order is important
  names(param) = c("structure","n", "alpha", "nInits", "m", "esp")
  fileToread = paste0('results/comparaison/Results', N, paste0(param, collapse = "_") )
  niceResults = read.csv(file =  fileToread)
  rowNames = niceResults$X
  niceResults$X  = NULL
  Times =  cbind(Times, niceResults$time)
  niceResults$time  = NULL
  colnames(niceResults) = paste0(structure, ' ', colnames(niceResults))
  niceResultsCombined = cbind(niceResultsCombined, niceResults)
  }

niceResultsCombined[,1] = NULL

niceResultsCombined
Times
niceResultsCombined = cbind(niceResultsCombined, apply(Times, 1, mean))
colnames(niceResultsCombined)[length(colnames(niceResultsCombined))] = "mean execution time"

rownames(niceResultsCombined) = rowNames

library(xtable)

xtable(niceResultsCombined[,c(1,3,5,7,9)])
niceResultsCombined[,c(1,3,5,7,9)]

xtable(niceResultsCombined[,c(2,4,6,8,10)])
niceResultsCombined[,c(2,4,6,10)]

niceResultsCombined[ , c(10,11)]





