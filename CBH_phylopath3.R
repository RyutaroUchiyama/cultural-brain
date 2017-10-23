


install.packages("devtools")
install.packages("ape")
install.packages("gtools")
devtools::install_github("Ax3man/phylopath") #needs to be latest github version
library(phylopath) 
library(ape)
library(gtools)
detach("package:ggm") # Disable ggm package temporarily because it masks 'phylopath::DAG' function. ggm will be reactivated later in the script.

# ================================================================================
# === Basic settings =============================================================
# ================================================================================

# --- read data (need to set path and file name) ---------
primatetree <- read.nexus("consensusTree_10kTrees_Primates_Version3.nex")
primates <- read.csv("cbh_analysis/CBHprimates.csv")    #enter path to data

Nvar <- 5 # number of variables 
dagobj <-  DAG(Inn ~ Inn, Juv~GS, BM~ECV) 
# phylopath::DAG(). Enter the variables that are to be included in model. The relationships between the variables do not matter here, although they need to be entered in x~y format.     
# For an odd-number size matrix, repeat one of the variables (any one) twice, e.g., DAG(V1~V2,V3~V3) for three variables. 
# The order of the variables is inverted as result of DAG(), e.g., DAG(L1~R1,L2~R2) becomes [R2,L2,R1,L1] in resulting matrix -- merely a stylistic point.

# ================================================================================
# ================================================================================

# --- preparing variables -----------
primates <- as.data.frame(primates)   # need to convert from tibble to data frame, in order to use row names as required by phylo_path() 
rownames(primates) <- primates[,1]    # row names need to be species names for phylo_path() function. 
bb <- is.element(rownames(primates), primatetree$tip.label)
primates <- primates[bb,] 

primates$soclearning.corrected = round(primates$sociallearning / log10(primates$articlecount_zr), digits=3)
primates$soclearning.corrected[intersect(which(primates$articlecount_zr < 10), which(primates$sociallearning==0))] <- NA
primates$innov.corrected = round(primates$innovation / log10(primates$articlecount_zr), digits=3)
primates$innov.corrected[intersect(which(primates$articlecount_zr < 10), which(primates$innovation==0))] <- NA

primates$logst.ecv <- as.vector(scale(primates$log_ecv,center=T,scale=T))
primates$logst.bodymass <- as.vector(scale(primates$log_bodymass,center=T,scale=T))
primates$logst.groupsize <- as.vector(scale(primates$log_groupsize,center=T,scale=T))
primates$logst.juvenileper <- as.vector(scale(primates$log_juvenileper,center=T,scale=T))
primates$logst.soclearning.corrected <- as.vector(scale(log10(primates$soclearning.corrected+1),center=T,scale=T))
primates$logst.innov.corrected <- as.vector(scale(log10(primates$innov.corrected+1),center=T,scale=T))
  # as.vector() turns these variables into atomic vectors, necessary for using phylopath::best and phylopath::est_DAG

primates$ECV <- primates$logst.ecv
primates$BM <- primates$logst.bodymass
primates$GS <- primates$logst.groupsize
primates$Juv <- primates$logst.juvenileper
primates$SocL <- primates$logst.soclearning.corrected
primates$Inn <- primates$logst.innov.corrected

dagobj[1:Nvar,1:Nvar] <- matrix(0,Nvar,Nvar) # reset all entries of DAG object to 0
zerosones <- permutations(2,Nvar^2,v=c(0,1),repeats.allowed=TRUE) #gtools::permutations - generate all possible binary strings, to allocate to DAG matrices

# --- Remove models with reflexive edges (because these are not acyclic)
cycles <- matrix(0,dim(zerosones)[1],1)
for (i in 1:Nvar){
  cycles[which(zerosones[,Nvar*(i-1)+i] == 1)] <- 1
}
zerosones_2 <-  zerosones[-which(cycles==1),]
rm(zerosones)

# --- Allocate all possible binary strings to (a list of) DAG objects
models_1 <- rep(list(dagobj),dim(zerosones_2)[1])
for (i in 1:dim(zerosones_2)[1]){
  models_1[[i]][,] <- zerosones_2[i,]
}
rm(zerosones_2)

# --- Further remove models, this time ones that have bidirectional edges
cycles_2 <- matrix(0,length(models_1),1)
for (k in 1:length(models_1)){
  cycles_2[k] <- length(intersect( which(models_1[[k]] == t(models_1[[k]])), which(models_1[[k]]==1))) >0
}
models_2 <- models_1[which(cycles_2==0)]
rm(models_1)

# --- then do an additional check for cyclic models using isAcyclic()
library(ggm)
cycles_3 <- matrix(0,length(models_2),1)
for (i in 1:length(models_2)){
  cycles_3[i] <- isAcyclic(models_2[[i]])==0
}
models_3 <- models_2[which(cycles_3==0)]
rm(models_2)

# --- remove models for which d-separation cannot be evaluated due to full connectedness, using ggm::basiSet()
cycles_basisnull <- matrix(0,length(models_3),1)
for (i in 1:length(models_3)){
  cycles_basisnull[i] <- is.null(basiSet(models_3[[i]])) # ggm::basiSet
}
detach("package:ggm") # Disable ggm package because it masks DAG() function in 'DAG::phylopath'
models <- models_3[which(cycles_basisnull==0)] # 519 with 4 vars, 19 with 3 vars.
rm(models_3)

# --- give names to models
modelnames <- matrix(NA,length(models),1)
for (i in 1:length(models)){
  modelnames[i] <- paste("m",i,sep="")
}
names(models) <- modelnames

# --- results of analysis
### result <- phylo_path( models, rhino, rhino_tree)
result <- phylo_path( models, primates, primatetree)
### resultsummary <- summary(result)[1:100,]
resultsummary <- summary(result)[1:10,]

# --- get names and graphs of N best models
bestmodels_N <- ifelse(dim(resultsummary)[1] > 100, 100, dim(resultsummary)[1]) 
bestmodels_names <- matrix(NA,bestmodels_N,1)
for (i in 1:bestmodels_N){
  bestmodels_names[i] <- summary(result)[i,"model"]
}
bestmodels_causalmodels <-  result$model_set[bestmodels_names]

# --- get path coefficients and SDs for the top N models (phylopath::est_DAG)
bestmodels_coefs <- list()
bestmodels_SEs <- list()
for (i in 1:bestmodels_N){
  bestmodels_coefs[[i]] <- est_DAG( bestmodels_matrices[[i]], result$data, result$tree, model="lambda")$coef
  bestmodels_SEs[[i]] <- est_DAG( bestmodels_matrices[[i]], result$data, result$tree, model="lambda")$se
  names(bestmodels_coefs)[i] <- bestmodels_names[i]
  names(bestmodels_SEs)[i] <- bestmodels_names[i]
}

# --- output
timestamp <- gsub(" ","_", gsub(":", "-", Sys.time()))
write.csv(resultsummary, paste("resultsummary_",timestamp,".csv", sep="") )
write.csv(bestmodels_causalmodels, paste("bestmodels_causalmodels_",timestamp,".csv", sep="") )
write.csv(bestmodels_coefs, paste("bestmodels_coefs_",timestamp,".csv", sep="") )
write.csv(bestmodels_SEs, paste("bestmodels_SEs_",timestamp,".csv", sep="") )




