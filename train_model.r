library(dplyr)
library(splines)
library(properties)
library(tidyr)
library(FactoMineR)
library(nnet)
library(zoo)
library(lubridate)
library(parallel)
library(MGLM)
source("Models/MGLM_functions.R")
library(PLNmodels)
library(ggplot2)


result_directory = "results"
ifelse(!dir.exists(file.path(result_directory)), dir.create(file.path(result_directory)), FALSE)  


##########################################################################################################
#SCORE FUNCTION FOR COMPARING OBSERVED vs FOUND SEGMENTS
##########################################################################################################

scoreDiff = function(Tis_trouv, Tis_obs){
  groups_found = round(Tis_trouv)*index(tau_obs)
  d_groups_found = as.data.frame(groups_found) %>% 
    mutate(change = ifelse( (V1==0 & lead(V1!=0) & lag(V1==0)) | (V2==0 & lead(V2)!=0 & lag(V2==0)) | (V3==0 & lead(V3)!=0 & lag(V3==0))  , 1, 0)) %>% 
    replace(is.na(.), 0)
  d_groups_obs = as.data.frame(Tis_obs) %>% 
    mutate(change = ifelse( (V1==0 & lead(V1!=0) & lag(V1==0)) | (V2==0 & lead(V2)!=0 & lag(V2==0)) | (V3==0 & lead(V3)!=0 & lag(V3==0))  , 1, 0)) %>% 
    replace(is.na(.), 0)
  
  scoreDiff = sum(d_groups_found$change!=d_groups_obs$change)
  return(scoreDiff)
}



##########################################################################################################
#FUNCTION TO VISUALIZE THE RESULTS
##########################################################################################################

output_segments <- function(Res, try){
  Seg_df = as.data.frame(Res$Tis)
  colnames(Seg_df) = c(paste0('Seg', as.character(seq(ncol(Seg_df)-1))),'j')
  Seg_df = Seg_df %>% pivot_longer(cols = -c(j), names_to = 'Segment') 
  
  
  PLOT_SEGMENTS = ggplot(Seg_df 
                         ,aes(x=j,y=value)) + geom_line() +  facet_grid(Segment~.)  +
    theme( legend.position = "top",legend.direction = "horizontal", strip.text.y = element_text(size=10,face='bold'),
           strip.text.y.right = element_text(angle = 0),
           legend.box = "horizontal", panel.spacing = unit(2, "lines"),
           legend.justification = c(0, 1),
           legend.title=element_text(size=30, face='bold'),
           axis.title=element_text(size=25,face="bold"), axis.text = element_text(size=10),axis.text.x = element_text( vjust = 0.5),
           legend.text = element_text(size = 20,face='bold'),legend.key.width=unit(4, "cm"),legend.key.height=unit(1, "cm"))+
    ylab(expression(tau[js])) +
    scale_y_continuous(breaks = c(0,0.5,1)) +
    xlab('j') + ggtitle(paste0("S = ",as.character(Res$S),", BIC = ",as.character(round(Res$BIC))))
  
  png(paste0("results/Segments_try",as.character(try),".png"))
  print(PLOT_SEGMENTS)
  dev.off()
}



##########################################################################################################
#TRAIN MODEL TO FIND SEGMENTS
##########################################################################################################
myProps <- read.properties("parameters")
Ntries = as.numeric(myProps$`Ntries`)
model_train = eval(myProps$`model_train`)
nb_S = as.numeric(myProps$`nb_S`)
nb_try = as.numeric(myProps$`nb_try`)
max_try = as.numeric(myProps$`max_try`)
suite_jours = eval(parse(text=myProps$`suite_jours`))
prop_jours = eval(parse(text=myProps$`prop_jours`))
diff_dat = eval(parse(text=myProps$`diff_dat`))
model_sim = eval(myProps$`model_sim`)
nombreJ = as.numeric(myProps$`nombreJ`)
nombreT = as.numeric(myProps$`nombreT`)

formule_Seg = ~ bs(j,5) - 1
formule_mod = "~ factor(hr)"

#lists for gathering score metrics
bestS = c()
nbMisclass = c()

#N tries for training model
for(try in seq(1,Ntries)){
  
  #Data generation
  source("preprocess_data.r")
  X = readRDS(file = "data/X.rds")
  Y = readRDS(file = "data/Y.rds")
  tau_obs = readRDS(file = "data/Tis_obs.rds")

  
  #PoiMult
  if(model_train == "PoiMult"){
      source("Models/PoiMult.R")
      list_bic = c()
      list_Res= list()
      for(S in seq(2,nb_S)){
        Res = PoiMult_model(Y = Y, X = X, formule_Seg, formule_mod, formule_mod, S=S, itermax=25,
                      tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH")
        list_Res[[S]] = Res
        if(Res$convergence){
          list_bic = c(list_bic, Res$BIC)
        }
      }
  }
  
  
  #sPoiMult
  if(model_train == "sPoiMult"){
    source("Models/sPoiMult.R")
    list_bic = c()
    list_Res= list()
    for(S in seq(2,nb_S)){
      Res = sPoiMult_model(Y = Y, X = X, formule_Seg, formule_mod, formule_mod, S=S, itermax=25,
                          tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH")
      list_Res[[S]] = Res
      if(Res$convergence){
        list_bic = c(list_bic, Res$BIC)
      }
    }
  }
  
  
  #NegPol
  if(model_train == "NegPol"){
    source("Models/NegPol.R")
    paral = as.logical(myProps$paral)
    nbcores = as.numeric(myProps$nbcores)
    list_bic = c()
    list_Res= list()
    for(S in seq(2,nb_S)){
      Res = NegPol_model(Y = Y, X = X, formule_Seg, formule_mod, formule_mod, S=S, itermax=25,
                           tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH",paral=paral,nbcores=nbcores)
      list_Res[[S]] = Res
      if(Res$convergence){
        list_bic = c(list_bic, Res$BIC)
      }
    }
  }
  
  
  #sNegPol
  if(model_train == "sNegPol"){
    source("Models/sNegPol.R")
    paral = as.logical(myProps$paral)
    nbcores = as.numeric(myProps$nbcores)
    list_bic = c()
    list_Res= list()
    for(S in seq(2,nb_S)){
      Res = sNegPol_model(Y = Y, X = X, formule_Seg, formule_mod, formule_mod, S=S, itermax=25,
                         tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH",paral=paral,nbcores=nbcores)
      list_Res[[S]] = Res
      if(Res$convergence){
        list_bic = c(list_bic, Res$BIC)
      }
    }
  }
  
  
  
  #PLNdiag
  if(model_train == "PLNdiag"){
    source("Models/PLNdiag.R")
    list_bic = c()
    list_Res= list()
    for(S in seq(2,nb_S)){
      Res = PLNdiag_model(Y = Y, X = X, formule_Seg, formule_mod, S=S, itermax=25,
                          tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH")
      list_Res[[S]] = Res
      if(Res$convergence){
        list_bic = c(list_bic, Res$BIC)
      }
    }
  }
  
  
  
  #sPLNdiag
  if(model_train == "sPLNdiag"){
    source("Models/sPLNdiag.R")
    list_bic = c()
    list_Res= list()
    for(S in seq(2,nb_S)){
      Res = sPLNdiag_model(Y = Y, X = X, formule_Seg, formule_mod, S=S, itermax=25,
                          tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH")
      list_Res[[S]] = Res
      if(Res$convergence){
        list_bic = c(list_bic, Res$BIC)
      }
    }
  }
  
  
  
  #PLNfull
  if(model_train == "PLNfull"){
    source("Models/PLNfull.R")
    list_bic = c()
    list_Res= list()
    for(S in seq(2,nb_S)){
      Res = PLNfull_model(Y = Y, X = X, formule_Seg, formule_mod, S=S, itermax=25,
                           tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH")
      list_Res[[S]] = Res
      if(Res$convergence){
        list_bic = c(list_bic, Res$BIC)
      }
    }
  }
  
  
  #sPLNfull
  if(model_train == "sPLNfull"){
    source("Models/sPLNfull.R")
    list_bic = c()
    list_Res= list()
    for(S in seq(2,nb_S)){
      Res = sPLNfull_model(Y = Y, X = X, formule_Seg, formule_mod, S=S, itermax=25,
                          tol=2e-6, nb_try = nb_try, max_try = max_try, min_seg = 0, init = "CAH")
      list_Res[[S]] = Res
      if(Res$convergence){
        list_bic = c(list_bic, Res$BIC)
      }
    }
  }
  
  
  
  #Visualize some results
  if(try==1 | try%%5==0){
    output_segments(list_Res[[which(list_bic==min(unlist(list_bic)))[1] + 1]], try)
  }
  
  #Retrieval of results
  bestS= c(bestS, (which(list_bic==min(unlist(list_bic)))[1])+1) #best number of segments found per try
  
  #number of bad classifications for S=3
  groups_obs = as.data.frame(tau_obs*index(tau_obs))
  nbMisclass = c(nbMisclass, scoreDiff(list_Res[[3]]$Tis[,1:3], groups_obs))
}


#Write results in a text file
sink("results/results_simu.txt")
cat("Results of simulations with :")
cat("\n")
cat(paste0("Model to simulate data : ", eval(myProps$`model_sim`)))
cat(paste0("\nModel to train : ", model_train,"\n"))
cat("\n")
cat(paste("Data simulated : ",paste(paste0(as.character(nombreJ*prop_jours)," days of type ", as.character(suite_jours)  ) , collapse = ', '), "\n")[1] )
cat(paste0("A day is ", as.character(nombreT), " time slots\n"))
cat(paste0("The impacts of segments on each sub-series are : ", paste(as.character(diff_dat), collapse=' ,')))
cat("\n\n")
cat("Results :\n")
cat("Percentage of trials where the simulation found 3 segments : ", as.character(round((length(bestS==3)/length(bestS))*100)), "%\n" )
cat("Mean misclassification rate : ", as.character(round(mean(nbMisclass)/nombreJ)), "%" )
sink()
