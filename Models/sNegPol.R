
sNegPol_model <- function(Y = Y, X = X, formuleSeg = ~bs(j,10) , formulenb =  ~bs(hr,5), formuleMult =   ~bs(hr,5),
                        S=3,itermax=50,tol=10^-3, alphatol = 10^-3, nb_try = 5, max_try = 15, min_seg=3, init = "CAH",paral=TRUE,nbcores=3){
  Y=as.matrix(Y)
  Ytot = rowSums(Y)
  df = as.data.frame(cbind(Ytot, X))
  criterions = list()
  res = list()
  converged = FALSE
  all_tries = FALSE
  i = 0

  
  print(paste0("Model with ",as.character(S), " segments"))
  while((all_tries == FALSE) & i < max_try){
    hotRestart = TRUE
    i = i+1
    print(paste0('Try number : ', as.character(i)))
    
    
    Xmp = model.matrix(as.formula(formulenb),X)
    Xmm = model.matrix(as.formula(formuleMult),X)
    Xms = model.matrix(formuleSeg,as.data.frame(unique(X$j)) %>% `colnames<-`(c("j")))
    
    fin_iter = FALSE
    skip_try = FALSE
    nb_suppr = 0
    nbT = dim(X)[1]
    nbJ = dim(Xms)[1]
    nbH = nbT/nbJ
    
    capts = colnames(Y)
    fo = formula(paste0("cbind(",paste0(capts,collapse = ","),")", as.character(formuleMult)))
    
    ###############################################
    # Initialization segments
    ###############################################
    
    #Initial clustering 
    C_j = matrix(0,nbJ,S)
    
    #Initialization HAC   
    Y_t = as.data.frame(cbind(X$hr,X$j,Y[,sample(seq(1,ncol(Y)), size = .75*ncol(Y))]))
    colnames(Y_t)[1:2] = c('hr','j')
    Y_t = Y_t %>% pivot_longer(-c(hr,j),names_to = "poste", values_to = "nb") %>%
      mutate(poste_hr = paste0(poste,"_",hr)) %>%
      dplyr::select(j,poste_hr,nb) %>%
      spread(key=poste_hr,value=nb)
    
    #  PCA 
    res.pca <- PCA(Y_t, ncp = 3, graph = FALSE)
    # HCPC
    res.hcpc <- HCPC(res.pca,nb.clust=S, graph = FALSE)
    
    
    for(s in seq(1,S)){
      C_j[res.hcpc$data.clust$clust==s,s]=1
      C_j[res.hcpc$data.clust$clust==s,-s]=0/(S-1)
    }
    
    #First initialization of the parameters over a few characteristic days
    grp_jours = list()
    for(s in seq(1,S)){
      grp_jours[[s]] = sample(seq(1,nbJ),3,prob=C_j[,s], replace=TRUE)
    }
    parameters = list(nbinom=list(),dirmultinomial=list())
    
    for(s in seq(1,S)){
      jour_poids = grp_jours[[s]]
      vec_j = rep(0,nbJ)
      vec_j[jour_poids] = 1
      vec_poids = rep(vec_j, each=nbH)
      #Negative binomial
      tryCatch(
        parameters$nbinom[[s]] <- suppressWarnings(glm.nb(as.formula(paste0("Ytot",formulenb)), data = df, weights = vec_poids, control=list(epsilon = 1e-5,maxit=20,trace=F))),
        error = function(e){
          message("one segment has disappeared")
          skip_try <<- TRUE
        })
      if(skip_try) { next }
      
      #Dirichlet Multinomial
      tryCatch(
        parameters$dirmultinomial[[s]] <-  MGLMreg.fit(Y, init=NULL, Xmm, dist = "DM",weight = vec_poids, parallel = paral, cores = nbcores, maxiters = 20),
        error = function(e){
          message(e)
          skip_try <<- TRUE
        })
      if(skip_try) { next }
      
    }
    

    
    #First segments initialization
    Tis_j = matrix(0,nbJ,S)
    lTis_j = matrix(0,nbJ,S)
    
    for(s in seq(1,S)){
      thetapred=exp(model.matrix(formula(formuleMult),X)%*%parameters$dirmultinomial[[s]]@coefficients)
      lTis_j[,s]= as.matrix(as.data.frame(cbind(dnbinom(Ytot,mu=parameters$nbinom[[s]]$fitted.values, size = parameters$nbinom[[s]]$theta ,log = TRUE)+
                                                  ddirmn(Y,alpha=thetapred), X$j))  %>% 
                              `colnames<-`(c("ll","j")) %>% 
                              group_by(j) %>% summarise(ll=sum(ll)) %>% 
                              dplyr::select(ll))
    }
    
    lnorm =  apply(lTis_j,1,max) + log(apply(exp(lTis_j -(apply(lTis_j,1,max)%*%matrix(1,1,S))), 1, sum))
    lTis_j_n = lTis_j - lnorm%*%matrix(1,1,S)
    Tis_j = round(exp(lTis_j_n),5)
    
    
    iter = 0
    fin_iter = FALSE
    criterion = -Inf
    old_criterion = -Inf
    
    
    
    
    #The initialization of the parameters and the EM algo are integrated in a while loop which 
    #allows a hotRestart when a segment has disappeared
    
    nbr_supr = 0
    while(hotRestart != FALSE){
      hotRestart = FALSE
      
      if(iter > 0){
        
        ###############################################
        # Initialization of the parameters associated with the segments
        ###############################################

        #segments
        tryCatch(
          reg_alpha <- multinom(Tis_j ~ Xms, trace=FALSE,maxit=70,MaxNWts = 10000000),
          error = function(e){
            message(e)
            skip_try <<- TRUE
          })
        if(skip_try) { next }
        
        
        
        # Formatting in probas of belonging to each segment by time slot
        Tis_j = as.data.frame(Tis_j)
        Tis_j[,"j"] = unique(X$j)
        Tis_t = X %>%
          dplyr::select(j) %>%
          left_join(Tis_j, by = c('j'='j')) %>%
          dplyr::select(-c(j))
        Tis_t = as.matrix(Tis_t)
        
        
        
        for(s in seq(1,S)){
          Tis_t[,s][Tis_t[,s]==0] = 0
          #Negative binomial
          tryCatch(
            parameters$nbinom[[s]] <- suppressWarnings(glm.nb(as.formula(paste0("Ytot",formulenb)), data = df, weights = as.matrix(Tis_t[,s]), control=list(epsilon = 1e-5,maxit=20,trace=F))),
            error = function(e){
              message("One segment has disappeared")
              skip_try <<- TRUE
            })
          if(skip_try) { next }
          
          #Dirichlet Multinomial
          tryCatch(
            parameters$dirmultinomial[[s]] <- MGLMreg.fit(Y, init=NULL, Xmm, dist = "DM",weight =  Tis_t[,s], parallel = paral, cores = nbcores, maxiters = 20,display=T),
            error = function(e){
              message(e)
              skip_try <<- TRUE
            })
          if(skip_try) { next }
        }
        
      }
      
      
      
      
      
      
      ################################################
      #EM algorithm
      ################################################
      while(fin_iter != TRUE & iter < itermax & skip_try != TRUE){
        iter = iter + 1
        
        ##################
        # E-STEP
        ##################

        if(iter > 1){
          
          pijs = fitted(reg_alpha)
          lpijs = as.data.frame(log(pijs))
          
          for(s in seq(1,S)){
            thetapred=exp(model.matrix(formula(formuleMult),X)%*%parameters$dirmultinomial[[s]]@coefficients)
            lTis_j[,s]= as.matrix(as.data.frame(cbind(dnbinom(Ytot,mu=parameters$nbinom[[s]]$fitted.values, size = parameters$nbinom[[s]]$theta ,log = TRUE)+
                                                        ddirmn(Y,alpha=thetapred), X$j))  %>% 
                                    `colnames<-`(c("ll","j")) %>% 
                                    group_by(j) %>% summarise(ll=sum(ll)) %>% 
                                    dplyr::select(ll)) +
              lpijs[,s]
          }
          
          #Creation of normalized probas for segment membership
          lnorm =  apply(lTis_j,1,max) + log(apply(exp(lTis_j -(apply(lTis_j,1,max)%*%matrix(1,1,S))), 1, sum))
          lTis_j_n = lTis_j - lnorm%*%matrix(1,1,S)
          Tis_j = round(exp(lTis_j_n),5)
          
        }
        
        

        old_criterion=criterion
        criterion=sum(lnorm)
        criterions[[i]] = criterion
        
        #print((criterion-old_criterion)/abs(old_criterion))
        
        if(iter >1  && abs(criterion-old_criterion)/abs(old_criterion)<tol && (criterion-old_criterion) > 0 ){
          fin_iter=TRUE
        }
        
        #Test disappearance of a segment, if a segment has less than min_seg consecutive days, we launch a HotRestart
        max_cons = sapply(as.data.frame(Tis_j[,1:S]), function(x) max(with(rle(x==0), lengths[!values])) ) #largest number of consecutive days per segment
        if((sum(max_cons <= min_seg) > 0)){
          nbr_supr = nbr_supr + 1
          print(nbr_supr)
          li_it = 0
          while((sum(max_cons <= min_seg) > 0)){
            li_it = li_it + 1
            nb_deads = length(which(max_cons <= min_seg))
            list_seg_disp =  which(max_cons <= min_seg)
            mat_worst_probas = Tis_j[,1:S] * lTis_j
            mat_worst_probas[,list_seg_disp] = 0
            mult_elts = 5
            pos_elts_grps = data.frame()
            while(nrow(pos_elts_grps) < (length(list_seg_disp))+1){
              elts <- order(mat_worst_probas)[seq_len(mult_elts*nb_deads)]
              pos <- as.data.frame(arrayInd(elts, dim(mat_worst_probas), useNames = TRUE)) %>%
                arrange(row) %>%
                mutate(group=ifelse(((row == (lead(row)-1))), 0, 1)) %>%
                mutate(group = ifelse(is.na(group), 1, group)) %>%
                mutate(group_c = cumsum(group)) %>%
                group_by(group_c) %>% mutate(nb=n()) %>%
                mutate(nb = ifelse(group == 1, 1, nb)) %>%
                arrange(desc(nb),row,col)
              pos_elts =  pos %>% group_by(group_c,nb) %>% summarise(n=n()) %>% arrange(desc(nb))
              pos_elts_grps = pos_elts %>% filter(nb > (min_seg+1))
              mult_elts = mult_elts + 1
            }
            c = 0
            for(seg_disp in list_seg_disp){
              c = c + 1
              print(paste0("Segment ", as.character(seg_disp), " disappeared, creating another one"))
              
              random_group = sample(seq(1,nrow(pos_elts_grps)),1)
              elts_to_select = pos_elts_grps[random_group,]
              pos_sel = (pos %>% inner_join(elts_to_select, by=c('group_c'='group_c','nb'='nb')))[,c('row','col')]
              pos_elts_grps = pos_elts_grps[-random_group,]
              
              pos_sel = as.matrix(pos_sel)
              Tis_j[pos_sel] = rep(0,nrow(pos_sel))
              Tis_j[pos_sel[,1],seg_disp] = rep(1,nrow(pos_sel))
            }
            max_cons = sapply(as.data.frame(Tis_j[,1:S]), function(x) max(with(rle(x==0), lengths[!values])) ) 
            if(li_it > 10){
              next
            }
          }
          hotRestart = TRUE
          print("Hot Restart")
          if(hotRestart){break} #back to parameters initialization
          if(nbr_supr > 5){
            hotRestart = FALSE
            skip_try = TRUE
            next
          }
        }
        
        
        
        
        #################
        # M-STEP
        #################

        #segments : alpha parameters
        if(iter>1){
          tryCatch(
            #reg_alpha <- multinom(Tis_j ~ ., coefs_alpha,Wts=reg_alpha$wts, maxit=200,reltol=tol,trace=FALSE),
            reg_alpha <- multinom(Tis_j ~ Xms,Wts=reg_alpha$wts,trace=FALSE,maxit=70,MaxNWts = 10000000),
            error = function(e){
              message(e)
              skip_try <<- TRUE
            })
          if(skip_try) { next }
        } else{
          tryCatch(
            reg_alpha <- multinom(Tis_j ~ Xms,trace=FALSE,maxit=70,MaxNWts = 10000000),
            error = function(e){
              message(e)
              skip_try <<- TRUE
            })
          if(skip_try) { next }
          
        }
        
        

        # Formatting in probas of belonging to each segment by time slot
        Tis_j = as.data.frame(Tis_j)
        Tis_j[,"j"] = unique(X$j)
        Tis_t = X %>%
          dplyr::select(j) %>%
          left_join(Tis_j, by = c('j'='j')) %>%
          dplyr::select(-c(j))
        Tis_t = as.matrix(Tis_t)
        
        #Negative binomial and Dirichlet Multinomiale
        for(s in seq(1,S)){
          Tis_t[,s][Tis_t[,s]==0] = 0
          coefs_nbinom = parameters$nbinom[[s]]$coefficients
          coefs_nbinom[is.na(coefs_nbinom)] = 0
          
          #NB
          tryCatch(
            parameters$nbinom[[s]] <- suppressWarnings(glm.nb(as.formula(paste0("Ytot",formulenb)), data = df, start = coefs_nbinom, weights = as.matrix(Tis_t[,s]), control=list(epsilon = 1e-6, maxit=20,trace=F))),
            error = function(e){
              message(e)
              parameters$nbinom[[s]] <- suppressWarnings(glm.nb(as.formula(paste0("Ytot",formulenb)), data = df, weights = as.matrix(Tis_t[,s]), control=list(epsilon = 1e-5, maxit=20,trace=F)))
            })
          if(skip_try) { next }
          
          
          
          #Dirichlet Multinomial
          tryCatch(
            #parameters$dirmultinomial[[s]] <- MGLMtune(fo, data = cbind(Y,X), penalty="sweep", dist = "DM",weight = Tis_t[,s], maxiters=50),
            parameters$dirmultinomial[[s]] <- MGLMreg.fit(Y, init=parameters$dirmultinomial[[s]]@coefficients, Xmm, dist = "DM",weight =  Tis_t[,s], parallel = paral, cores = nbcores, maxiters=50,display = F),
            error = function(e){
              message(e)
              skip_try <<- TRUE
            })
          if(skip_try) { next }
        }
        
 
      }
      
      if(nbr_supr > 5) {
        skip_try=TRUE
        break
      }
    }
    

    
    if(skip_try){
      criterions[[i]] = -Inf
    }
    
    
    all_tries = sum(!is.infinite(unlist(criterions))) == nb_try
    
    parameters$segmentation=reg_alpha
    res[[i]] = list(Tis_j,lTis_j,parameters)
  }
  
  if(all_tries){converged = TRUE}
  
  criterion_max = max(unlist(criterions))
  criterion_max_idx = which(criterions==max(unlist(criterions)))[1]
  parameters_best = res[[criterion_max_idx]][[3]]
  Tis_best = res[[criterion_max_idx]][[1]]
  lTis_best = res[[criterion_max_idx]][[2]]
  
  n = dim(Y)[2]
  nu = (S-1)*dim(Xms)[2] + S * (dim(Xmp)[2]) + S * dim(Xmm)[2] * (n)
  BIC = -2 * criterion_max + nu * log(n*nbT)
 
  
  return(list(convergence = converged, S=S, crit=criterion_max, Tis=Tis_best, lTis=lTis_best, params=parameters_best, BIC = BIC))
  
}

