library(Matrix)
library(properties)
library(MASS)



##########################################################################################################
#ARGUMENTS
##########################################################################################################

myProps <- read.properties("parameters")
suite_jours = eval(parse(text=myProps$`suite_jours`))
prop_jours = eval(parse(text=myProps$`prop_jours`))
diff_dat = eval(parse(text=myProps$`diff_dat`))
model_sim = eval(myProps$`model_sim`)
nombreJ = as.numeric(myProps$`nombreJ`)
nombreT = as.numeric(myProps$`nombreT`)


##########################################################################################################
#DATA BUILDING
##########################################################################################################

#Covariates data
X_data = data.frame('hr' = rep(seq(1,nombreT), nombreJ),
                      'j' = rep(1:nombreJ, each=nombreT))
  
#Count data, 1 matrix per segment (3 segments)
Y_list = list(matrix(0,nrow=nombreJ*nombreT,ncol=5), 
                    matrix(0,nrow=nombreJ*nombreT,ncol=5),
                    matrix(0,nrow=nombreJ*nombreT,ncol=5))
  

#Segments positions: for days first, for timestamps then
tau_j =  as.matrix(t(sparseMatrix(t(rep(suite_jours,(nombreJ*prop_jours))), seq(1,nombreJ),x = 1)))
tau = as.matrix(t(sparseMatrix(t(rep(suite_jours,nombreT*(nombreJ*prop_jours))), seq(1,nombreJ*nombreT),x = 1)))


#Apply a perturbation of d% to each segment  
pert = matrix(rep(diff_dat, each = (nombreJ*nombreT)), ncol = 5)
pert_seg = list(pert, matrix(0, nrow = (nombreJ*nombreT), ncol=5), -pert)
  
#Elicitation of simulation parameters 

prop_mult <- t(replicate(nombreT, diff(c(0, sort(runif(4)), 1))))
distMoy = sample(300:2000, nombreT, replace = TRUE)
  


#Simulation data from specified model 

#Poisson log normal model
if(model_sim == "pln"){
    log_dist = log(distMoy * prop_mult)
    sigma = as.numeric(myProps$`sigma_pln`)
    
    for(s in seq(1, 3)){
      theta = t(sapply(1:nombreT, function(i)  MASS::mvrnorm(1, mu = log_dist[i,], Sigma = diag(sigma,5)) ))
      Y_list[[s]] =  t(sapply(1:nombreT, function(i)  rpois(5, exp(theta[i,]) ) )) 
      for(j in seq(2,nombreJ)){
        theta = t(sapply(1:nombreT, function(i)  MASS::mvrnorm(1, mu = log_dist[i,], Sigma = diag(sigma,5)) ))
        Y_list[[s]] = rbind(Y_list[[s]],  t(sapply(1:nombreT, function(i)  rpois(5, exp(theta[i,]) ) ))  )
      }
      Y_list[[s]] = Y_list[[s]] + pert_seg[[s]]*Y_list[[s]]/100
    }
  }
  
  
#Poisson sums and Multinomial share model
if(model_sim == "PoiMult"){

    for(s in seq(1, 3)){
      tot_NB =  rpois(nombreT,  lambda = distMoy)
      Y_list[[s]] =  t(sapply(1:nombreT, function(i)  rmultinom(1, tot_NB[i], prop_mult[i,]) )) 
      for(j in seq(2,nombreJ)){
        tot_NB =  rpois(nombreT, lambda = distMoy)
        Y_list[[s]] = rbind(Y_list[[s]],   t(sapply(1:nombreT, function(i)  rmultinom(1, tot_NB[i], prop_mult[i,]) )) )
      }
      Y_list[[s]] = Y_list[[s]] + pert_seg[[s]]*Y_list[[s]]/100
    }
}


#Negative binomial sums and Dirichlet Multinomial share model
if(model_sim == "NegPol"){
  size = as.numeric(myProps$`size_NegPol`)
  for(s in seq(1, 3)){
    tot_NB =  rnbinom(nombreT, size=size, mu = distMoy)
    Y_list[[s]] =  t(sapply(1:nombreT, function(i)  rgamma(5,rmultinom(1, tot_NB[i], prop_mult[i,])) )) 
    for(j in seq(2,nombreJ)){
      tot_NB =  rnbinom(nombreT, size=size, mu = distMoy)
      Y_list[[s]] = rbind(Y_list[[s]],   t(sapply(1:nombreT, function(i)  rgamma(5,rmultinom(1, tot_NB[i], prop_mult[i,])) )) )
    }
    Y_list[[s]] = Y_list[[s]] + pert_seg[[s]]*Y_list[[s]]/100
  }
}
  
  
  
#Associating counts generated from different segments  
Y_data = tau[,1]*Y_list[[1]] + tau[,2]*Y_list[[2]] + 
    tau[,3]*Y_list[[3]] 
  
  
X = X_data
Y = apply(Y_data, 2, as.integer)
Y[Y<0] = 1



##########################################################################################################
#SAVE DATA
##########################################################################################################
data_directory = "data"
ifelse(!dir.exists(file.path(data_directory)), dir.create(file.path(data_directory)), FALSE)  
saveRDS(X, file = "data/X.rds")
saveRDS(Y, file = "data/Y.rds")
saveRDS(tau_j, file = 'data/Tis_obs.rds')  
