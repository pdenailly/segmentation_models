



PLNfull_model <- function(Y = Y, X = X, formuleSeg = ~bs(j,10) , formPLN = "~factor(hr)", 
                             S=3,itermax=50,tol=10^-3, alphatol = 10^-3, nb_try = 5, max_try = 15, min_seg=3, init = "CAH"){
  Y=as.matrix(Y)
  rownames(Y)<-NULL
  data_Count = prepare_data(Y,X)
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
    
    Xms = model.matrix(formuleSeg,as.data.frame(unique(X$j)) %>% `colnames<-`(c("j")))
    
    fin_iter = FALSE
    skip_try = FALSE
    nb_suppr = 0
    nbT = dim(X)[1]
    nbJ = dim(Xms)[1]
    nbH = nbT/nbJ
    
    ###############################################
    # Initialisation segments
    ###############################################
    
    #Clustering initial
    C_j = matrix(0,nbJ,S)
    
    #Initialisation CAH   
    Y_t = as.data.frame(cbind(X$hr,X$j,Y[,sample(seq(1,ncol(Y)), size = .75*ncol(Y))]))
    colnames(Y_t)[1:2] = c('hr','j')
    Y_t = Y_t %>% pivot_longer(-c(hr,j),names_to = "poste", values_to = "nb") %>%
      mutate(poste_hr = paste0(poste,"_",hr)) %>%
      dplyr::select(j,poste_hr,nb) %>%
      spread(key=poste_hr,value=nb)
    
    #  ACP 
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
    parameters = list(pln=list())
    
    for(s in seq(1,S)){
      jour_poids = grp_jours[[s]]
      vec_j = rep(0,nbJ)
      vec_j[jour_poids] = 1
      vec_poids = rep(vec_j, each=nbH)
      vec_poids[vec_poids==0] = 1e-5
      tryCatch(
        parameters$pln[[s]] <- suppressWarnings(PLN(as.formula(paste("Abundance",as.character(formPLN))),data_Count, weights  = vec_poids  , control = list(covariance = "full",maxeval=500,trace=0))),
        error = function(e){
          message(e)
          skip_try <<- TRUE
        })
      if(skip_try) { next }
      
    }
    
    #1st initialisation segments
    Tis_j = matrix(0,nbJ,S)
    lTis_j = matrix(0,nbJ,S)
    
    for(s in seq(1,S)){
      lTis_j[,s] =  as.matrix(as.data.frame(cbind(parameters$pln[[s]]$loglik_vec,X$j)) %>%
                                rename(j = names(.)[length(names(.))]) %>%
                                `colnames<-`(c("ll","j")) %>%
                                group_by(j) %>% summarise(ll=sum(ll)) %>% 
                                dplyr::select(ll)
      ) 
    }
    lnorm =  apply(lTis_j,1,max) + log(apply(exp(lTis_j -(apply(lTis_j,1,max)%*%matrix(1,1,S))), 1, sum))
    lTis_j_n = lTis_j - lnorm%*%matrix(1,1,S)
    Tis_j = round(exp(lTis_j_n),5)
    
    
    iter = 0
    fin_iter = FALSE
    criterion = -Inf
    old_criterion = -Inf
    
    
    nbr_supr = 0
    while(hotRestart != FALSE){
      hotRestart = FALSE
      
      if(iter > 0){
        
        ###############################################
        # Initialization of the parameters associated with the segments
        ###############################################

        
        #segments
        pijs = sapply(as.vector(colSums(Tis_j)/nbJ), function(i) rep(i,nbJ))
        
 
        
        # #Formatting in probas of belonging to each segment by slice
        Tis_j = as.data.frame(Tis_j)
        Tis_j[,"j"] = unique(X$j)
        Tis_t = X %>%
          dplyr::select(j) %>%
          left_join(Tis_j, by = c('j'='j')) %>%
          dplyr::select(-c(j))
        Tis_t = as.matrix(Tis_t)
        
        
        
        
        for(s in seq(1,S)){
          Tis_t[,s][Tis_t[,s]==0] = .Machine$double.eps
          tryCatch(
            parameters$pln[[s]] <- suppressWarnings(PLN(as.formula(paste("Abundance",as.character(formPLN))),data_Count, weights =  Tis_t[,s] , control = list(covariance = "full",maxeval=500,trace=0))),
            error = function(e){
              message(e)
              skip_try <<- TRUE
            })
          if(skip_try) { next }
          
        }
        
      }
      
      
      Xmm_PLN = model.matrix(as.formula(formPLN),X)
      
      
      
      ################################################
      #EM algorithm
      ################################################
      while(fin_iter != TRUE & iter < itermax & skip_try != TRUE){
        iter = iter + 1
        
        
        ##################
        # E-STEP
        ##################

        if(iter > 1){
          #Probability of belonging to the segments a posteriori
          lpijs = as.data.frame(log(pijs))
          
          
          for(s in seq(1,S)){
            lTis_j[,s] =  as.matrix(as.data.frame(cbind(parameters$pln[[s]]$loglik_vec,X$j)) %>%
                                      rename(j = names(.)[length(names(.))]) %>%
                                      #pivot_longer(-j) %>%
                                      `colnames<-`(c("ll","j")) %>%
                                      group_by(j) %>% summarise(ll=sum(ll)) %>% 
                                      dplyr::select(ll)
            ) +
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
        
        if(iter >1  && abs(criterion-old_criterion)/abs(old_criterion)<tol){
          fin_iter=TRUE
        }
        
        
        
        #Test disappearance of a segment, if a segment has less than min_seg consecutive days, we launch a HotRestart
        max_cons = sapply(as.data.frame(Tis_j[,1:S]), function(x) max(with(rle(x==0), lengths[!values])) ) 
        if((sum(max_cons <= min_seg) > 0)){
          nbr_supr = nbr_supr + 1
          while((sum(max_cons <= min_seg) > 0)){
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
              print(paste0("Le segment ", as.character(seg_disp), " meurt, formation d'un nouveau segment"))
              
              random_group = sample(seq(1,nrow(pos_elts_grps)),1)
              elts_to_select = pos_elts_grps[random_group,]
              pos_sel = (pos %>% inner_join(elts_to_select, by=c('group_c'='group_c','nb'='nb')))[,c('row','col')]
              pos_elts_grps = pos_elts_grps[-random_group,]
              
              pos_sel = as.matrix(pos_sel)
              Tis_j[pos_sel] = rep(0,nrow(pos_sel))
              Tis_j[pos_sel[,1],seg_disp] = rep(1,nrow(pos_sel))
            }
            max_cons = sapply(as.data.frame(Tis_j[,1:S]), function(x) max(with(rle(x==0), lengths[!values])) ) 
          }
          hotRestart = TRUE
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

        #segments
        pijs = sapply(as.vector(colSums(Tis_j)/nbJ), function(i) rep(i,nbJ))
        
        
        if(S <= 10){
          plot.ts(Tis_j)
        } else{
          plot.ts(Tis_j[,1:10])
        }
        
        
        ## UPDATE THE MIXTURE MODEL VIA OPTIMIZATION OF PLNmixture
        
        # #Formatting in probas of belonging to each segment by slice
        Tis_j = as.data.frame(Tis_j)
        Tis_j[,"j"] = unique(X$j)
        Tis_t = X %>%
          dplyr::select(j) %>%
          left_join(Tis_j, by = c('j'='j')) %>%
          dplyr::select(-c(j))
        Tis_t = as.matrix(Tis_t)
        
        
        
        for(s in seq(1,S)){
          Tis_t[,s][Tis_t[,s]==0] = .Machine$double.eps
          tryCatch(
            parameters$pln[[s]] <- suppressWarnings(PLN(as.formula(paste("Abundance",as.character(formPLN))),data_Count, weights =Tis_t[,s] , 
                                                        control = list(inception = parameters$pln[[s]], covariance = "full",maxeval=500,trace=0))),
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
    
    res[[i]] = list(Tis_j,lTis_j,parameters)
  }
  
  if(all_tries){converged = TRUE}
  
  criterion_max = max(unlist(criterions))
  criterion_max_idx = which(criterions==max(unlist(criterions)))[1]
  parameters_best = res[[criterion_max_idx]][[3]]
  Tis_best = res[[criterion_max_idx]][[1]]
  lTis_best = res[[criterion_max_idx]][[2]]
  
  
  
  n = dim(Y)[2]
  nu =  S * dim(Xmm_PLN)[2] *n  + S * n*(n+1)/2
  BIC = -2 * criterion_max + nu * log(n*nbT)
  
  return(list(convergence = converged, S=S, crit=criterion_max, Tis=Tis_best, lTis=lTis_best, params=parameters_best, BIC = BIC))
  
}
