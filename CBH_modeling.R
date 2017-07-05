
#>> Note:
#>> Very high group size and very low prec score in golden monkey is driving c~e up to the position of 
#>> a dominant model, despite it being an outlier. I have thus deleted precocial score for golden monkey 
#>> to avoid this data point.  


# ===== Set parameters here =====================

input.data = CBHdata
combs.index <- c(0,1,6,16,26,31) # -- indices of rows where the number of predictor variables changes.

# ^ set "these"combs.index" manually according to the number of variables. 
#   Current setting is for combinations of six variables (1 DV, 5 IVs)
#   If more than six, edit code in a few indicated places.
#   Can probably automate this using choose(n,k) but too lazy right now.

# Current variables:
# >>(a) ecv, >>(b) bodymass, >>(c) TT.precscore, >>(d) juv.period, >>(e) pop.groupsize, >>(f) cylinder.score.

# IN THE input.data DATA FRAME, USE LOWER-CASE LETTERS (IN ORDER FROM "a") FOR VARIABLE NAMES TO BE ANALYZED.

# ======= Generate models combinatorially =====================
combs.varindices <- letters[1:length(combs.index)] 
for (v in 1:length(combs.index)){
  list.vessel <- list()
  letternum <- v # for which variable
  for (j in 1:(length(combs.varindices)-2)){
    combs.DVcount <- length(combs.varindices)-j
    combs <- t(combn(combs.varindices[-letternum], combs.DVcount))
    combs.p <- matrix("+",dim(combs)[1],dim(combs)[2]-1)
    combs.x <- matrix(NA,dim(combs)[1],dim(combs)[2]*2-1)
    for (i in 1:(dim(combs)[2]-1)){
      combs.x[,(i-1)*2+1] <- combs[,i]
      combs.x[,(i-1)*2+2] <- combs.p[,i] }
    combs.x[,dim(combs.x)[2]] <- combs[,dim(combs)[2]]
    combs.1 <- matrix(rep(c(letters[letternum]),dim(combs.x)[1]), dim(combs.x)[1],1)
    combs.2 <- matrix(rep(c("~"),dim(combs.x)[1]), dim(combs.x)[1],1)
    combs.range <- (combs.index[j]+1) : combs.index[j+1]
    if (j==1){
      list.vessel[1]<- matrix(apply(cbind(combs.1,combs.2,combs.x), 1, paste, collapse="")) }
    else{ for (k in 1:length(combs.range)){
      row <- (combs.range[1]:combs.range[k])[k]
      list.vessel[row]<- matrix(apply(cbind(combs.1,combs.2,combs.x), 1, paste, collapse=""))[k] 
    }}} 
  combs.singlesinds <-  (combs.index[j+1]+1):combs.index[j+2]
  combs.singles <- combs.varindices[-letternum]
  combs.1 <- matrix(rep(c(letters[letternum]),length(combs.singlesinds)), length(combs.singlesinds),1)
  combs.2 <- matrix(rep(c("~"),length(combs.singlesinds)), length(combs.singlesinds),1)
  list.vessel[combs.singlesinds] <- matrix(apply(cbind(combs.1,combs.2,combs.singles), 1, paste, collapse=""))
  if (v==1){ list.a <- list.vessel }; if (v==2){ list.b <- list.vessel }; if (v==3){ list.c <- list.vessel }
  if (v==4){ list.d <- list.vessel }; if (v==5){ list.e <- list.vessel }; if (v==6){ list.f <- list.vessel }
}   
# ----- ^ edit if more variables.


# ======= Compute model statistics =====================
matrix.aic <- matrix(NA,31,6)
matrix.aic.AIC <- matrix(NA,31,6)
matrix.models <- matrix(NA,31,6)  
matrix.fstatistic <- matrix(NA,31,6)
matrix.rsquared <- matrix(NA,31,6)  
matrix.paramnum <- matrix(NA,31,6) 
model.all =list(list.a,list.b,list.c,list.d,list.e,list.f) # -- edit if more variables.

for (j in 1:length(combs.index)){
list.input <- model.all[[j]]
inputcol <- j
for (i in 1:length(list.input)){
  fit.model <- lm(list.input[[i]], data=input.data)
  lmsummary <- summary(fit.model)
  k <- length(fit.model$coefficients)
  matrix.models[i, inputcol] <-  list.input[[i]]
  matrix.aic[i, inputcol] <- round( 2*k - logLik(fit.model), digits=3) 
  matrix.aic.AIC[i, inputcol] <- round(AIC(fit.model), digits=3) 
  matrix.paramnum[i, inputcol] <-  k
  matrix.fstatistic[i, inputcol] <- round( lmsummary$fstatistic[[1]], digits=3)
  matrix.rsquared[i, inputcol] <- round( lmsummary$r.squared, digits=3) 
}}

matrix.daic <- round(matrix.aic - min(matrix.aic),digits=4)
matrix.weights <- exp(-0.5*matrix.daic)/sum(exp(-0.5*matrix.daic))


# ======= Compute parameter weights =====================
param.weights <- data.frame(matrix(NA,30,2))
param.weights[1,1] <- sum(matrix.weights[grep("a~b",matrix.models)]); param.weights[16,1] <- sum(matrix.weights[grep("b~a",matrix.models)])
param.weights[2,1] <- sum(matrix.weights[grep("a~c",matrix.models)]); param.weights[17,1] <- sum(matrix.weights[grep("c~a",matrix.models)])
param.weights[3,1] <- sum(matrix.weights[grep("a~d",matrix.models)]); param.weights[18,1] <- sum(matrix.weights[grep("d~a",matrix.models)])
param.weights[4,1] <- sum(matrix.weights[grep("a~e",matrix.models)]); param.weights[19,1] <- sum(matrix.weights[grep("e~a",matrix.models)])
param.weights[5,1] <- sum(matrix.weights[grep("a~f",matrix.models)]); param.weights[20,1] <- sum(matrix.weights[grep("f~a",matrix.models)])
param.weights[6,1] <- sum(matrix.weights[grep("b~c",matrix.models)]); param.weights[21,1] <- sum(matrix.weights[grep("c~b",matrix.models)])
param.weights[7,1] <- sum(matrix.weights[grep("b~d",matrix.models)]); param.weights[22,1] <- sum(matrix.weights[grep("d~b",matrix.models)])
param.weights[8,1] <- sum(matrix.weights[grep("b~e",matrix.models)]); param.weights[23,1] <- sum(matrix.weights[grep("e~b",matrix.models)])
param.weights[9,1] <- sum(matrix.weights[grep("b~f",matrix.models)]); param.weights[24,1] <- sum(matrix.weights[grep("f~b",matrix.models)])
param.weights[10,1] <- sum(matrix.weights[grep("c~d",matrix.models)]); param.weights[25,1] <- sum(matrix.weights[grep("d~c",matrix.models)])
param.weights[11,1] <- sum(matrix.weights[grep("c~e",matrix.models)]); param.weights[26,1] <- sum(matrix.weights[grep("e~c",matrix.models)])
param.weights[12,1] <- sum(matrix.weights[grep("c~f",matrix.models)]); param.weights[27,1] <- sum(matrix.weights[grep("f~c",matrix.models)])
param.weights[13,1] <- sum(matrix.weights[grep("d~e",matrix.models)]); param.weights[28,1] <- sum(matrix.weights[grep("e~d",matrix.models)])
param.weights[14,1] <- sum(matrix.weights[grep("d~f",matrix.models)]); param.weights[29,1] <- sum(matrix.weights[grep("f~d",matrix.models)])
param.weights[15,1] <- sum(matrix.weights[grep("e~f",matrix.models)]); param.weights[30,1] <- sum(matrix.weights[grep("f~e",matrix.models)])
param.weights[,1] <- round(param.weights[,1], digits=3)
param.weights[1,2] <- "a~b"; param.weights[2,2] <- "a~c"; param.weights[3,2] <- "a~d"; param.weights[4,2] <- "a~e"
param.weights[5,2] <- "a~f"; param.weights[6,2] <- "b~c"; param.weights[7,2] <- "b~d"; param.weights[8,2] <- "b~e"
param.weights[9,2] <- "b~f"; param.weights[10,2] <- "c~d"; param.weights[11,2] <- "c~e"; param.weights[12,2] <- "c~f"
param.weights[13,2] <- "d~e"; param.weights[14,2] <- "d~f"; param.weights[15,2] <- "e~f";
param.weights[16,2] <- "b~a"; param.weights[17,2] <- "c~a"; param.weights[18,2] <- "d~a"; param.weights[19,2] <- "e~a"
param.weights[20,2] <- "f~a"; param.weights[21,2] <- "c~b"; param.weights[22,2] <- "d~b"; param.weights[23,2] <- "e~b"
param.weights[24,2] <- "f~b"; param.weights[25,2] <- "d~c"; param.weights[26,2] <- "e~c"; param.weights[27,2] <- "f~c"
param.weights[28,2] <- "e~d"; param.weights[29,2] <- "f~d"; param.weights[30,2] <- "f~e"
param.weights.bi <- data.frame("bidirectional sums"=(param.weights[1:15,1] + param.weights[16:30,1]), "param elements"= param.weights[1:15,2])

# ======= Find best-fit models =====================
bestmodels <- data.frame(matrix(NA,18,4))
colnames(bestmodels) <- c("variable","model","row","AIC")
for (i in 1:6){
  mins <- sort(matrix.aic[,i], index.return=T)$ix[1:3]
  bestmodels[((i-1)*3+1):((i-1)*3+3),1] <- rep(letters[i],3)
  bestmodels[((i-1)*3+1):((i-1)*3+3),2] <- matrix.models[mins,i] 
  bestmodels[((i-1)*3+1):((i-1)*3+3),3] <- mins 
  bestmodels[((i-1)*3+1):((i-1)*3+3),4] <- round(matrix.aic[mins,i],digits=2) 
}


write.csv(bestmodels, "best_models.csv")









