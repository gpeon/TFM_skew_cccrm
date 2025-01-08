ccc_vc_bayes <- function(dataset,ry, rind, rmet = NULL, rtime = NULL, covar = NULL, cl = 0.95, parallel = FALSE,workers = 1, n_chains = 4, n_burnin = 1000, n_iter = 2000, priors = brms::empty_prior(), distr_family = 'gaussian', correct = FALSE, int = TRUE, log_alpha = FALSE, plot.gel = FALSE){
  
  # remains to implement: parallel = T
  
  # Parallellization
  #-------------------------------------------------------
  ncores <- parallelly::availableCores(omit = 1)
  if(workers >= ncores){
    workers <- ncores
  }
  
  if(parallel){
    oplan <- future::plan("multisession", workers = workers)
    on.exit(future::plan(oplan))
    message("Parallel execution.")
    message(paste("Number of cores available:",ncores,sep=" "))
    message(paste("Number of cores used:",workers,sep=" "))
  }
  
  # Dataset
  #-------------------------------------------------------
  dades<- data.frame(dataset) |> dplyr::select(any_of(c(ry,rind,rtime,rmet,covar)))
  
  dades <- dades |> dplyr::rename(y = all_of(ry),
                                  ind = all_of(rind),
                                  met = all_of(rmet),
                                  time = all_of(rtime))
  
  dades$y<- as.numeric(dades$y)
  dades$ind<- as.factor(dades$ind)
  
  # Formula
  #-------------------------------------------------------
  
  form <- y ~ 1 + (1|ind)
  
  if(!is.null(rmet)){
    dades$met <- as.factor(dades$met)
    if(int){
      form <- update.formula(form,~.+met+(0+met|ind))
    }else{
      form <- update.formula(form,~.+met)
    }
  }
  
  if(!is.null(rtime)){
    dades$time <- as.factor(dades$time)
    if(!is.null(rmet)){
      form <- update.formula(form,~.+time+met*time+(0+time|ind))
    }else{
      form <- update.formula(form,~.+time+(0+time|ind))
    }
  }
  
  if(!is.null(covar)){
    covar_update <- paste(deparse1(form),"+",paste(covar,collapse =" + "))
    
    form <- as.formula(covar_update)
  }
  
  
  dades_stan <- brms::standata(form,data = dades,family = distr_family, prior = priors)
  
  # Necessary information for the CCC (The contrast matrix)
  #-------------------------------------------------------
  
  if(is.null(covar)){
    length_b <- dades_stan$K
  }else{
    length_b <- dades_stan$M_2
    if(!is.null(rtime)){
      length_b <- length_b*dades_stan$M_3
    }
  }
  
  
  nm <- dades_stan$M_2
  nd <- nm*(nm-1)/2
  
  
  if(!is.null(rtime)){
    nt <- dades_stan$M_3
    # All differences design matrix
    C<-array(0,dim=c(length_b,nt*nm))
    k<-0
    for (i in 1:nm){
      for (j in 1:nt){
        k<-k+1
        C[,k]<-c(1,contrasts(dades$met)[i,],contrasts(dades$time)[j,],
                 c(contrasts(dades$met)[i,]%*%t(contrasts(dades$time)[j,])))
      }
    }
    # Difference between methods at each time matrix
    L<-array(0,dim=c(length_b,nt*nd))
    k<-0
    for (i in 1:(nt*(nm-1))){
      for (j in (1:(nm-1))){
        if ((i+nt*j)<=(nt*nm)){
          k<-k+1
          L[,k]=C[,i]-C[,i+nt*j]
        }
      }
    }
  }else{
    aux1 <- diag(nm)
    aux1[1,1:nm] <- 1
    
    L=array(0,dim=c(nm,nd))
    cont=0
    for (i in 1:nm){
      for (j in 1:(nm-1)){
        if (i>j){
          cont=cont+1
          L[,cont]=aux1[,i]-aux1[,j]
        }
      }
    }
  }
  
  
  auxD <- diag(nrow = ncol(L),ncol = ncol(L))
  dades_stan$auxD <- auxD
  dades_stan$L <- L
  
  # Stan code editing
  #-------------------------------------------------------
  
  stanvars_mod <- brms::stanvar(scode = "vector[M_1] sd_1 = rep_vector(sd_1_scalar, M_1);",
                                block = "tparameters", position = "start")
  if((!is.null(rmet)) & int){
    stanvars_mod <- stanvars_mod + brms::stanvar(scode = "vector[M_2] sd_2 = rep_vector(sd_2_scalar, M_2);",
                                                 block = "tparameters", position = "start")
    if(!is.null(rtime)){
      stanvars_mod <- stanvars_mod + brms::stanvar(scode = "vector[M_3] sd_3 = rep_vector(sd_3_scalar, M_3);",
                                                   block = "tparameters", position = "start")
    }
  }
  
  # Add ccc (Longitudinal,Non-longitudinal with subject-method interaction, Non-longitudinal)
  if(!is.null(rtime)){
    
    if(is.null(covar)){
      stanvars_mod <- stanvars_mod + brms::stanvar(scode = "matrix[K, M_3*(M_2-1)] L;",
                                                   block = "data", x = L) + 
        brms::stanvar(scode = "matrix[M_3*(M_2-1), M_3*(M_2-1)] auxD;",
                      block = "data", x = auxD) +
        brms::stanvar(scode = "\n  vector[Kc + 1] fixed_effects;\n  fixed_effects[1] = b_Intercept;\n  for (k in 1:Kc) {\n    fixed_effects[k + 1] = b[k];\n  }", block = "genquant") +
        brms::stanvar(scode = "vector[M_3] difmed = transpose(L) * fixed_effects;", block = "genquant") +
        brms::stanvar(scode = "real aux1 = (transpose(difmed) * auxD * difmed); ", block = "genquant") +
        brms::stanvar(scode = "real var_b = fmax(aux1/(M_2*(M_2-1)*M_3),0);", block = "genquant")
    }else{
      stanvars_mod <- stanvars_mod + brms::stanvar(scode = "matrix[M_3*M_2, M_3*(M_2-1)] L;",
                                                   block = "data", x = L) + 
        brms::stanvar(scode = "matrix[M_3*(M_2-1), M_3*(M_2-1)] auxD;",
                      block = "data", x = auxD) +
        brms::stanvar(scode = "\n  vector[M_3*M_2] fixed_effects;\n  fixed_effects[1] = b_Intercept;\n  for (k in 1:(M_3*M_2-1)) {\n    fixed_effects[k + 1] = b[k];\n  }", block = "genquant") +
        brms::stanvar(scode = "vector[M_3] difmed = transpose(L) * fixed_effects;", block = "genquant") +
        brms::stanvar(scode = "real aux1 = (transpose(difmed) * auxD * difmed); ", block = "genquant") +
        brms::stanvar(scode = "real var_b = fmax(aux1/(M_2*(M_2-1)*M_3),0);", block = "genquant")
    }

    
    
    if(!log_alpha){
      stanvars_mod <- stanvars_mod +
        brms::stanvar(scode = "real ccc_num = sd_1_scalar^2+sd_3_scalar^2;", block = "genquant") + 
        brms::stanvar(scode = "real ccc_denom = sd_1_scalar^2+sd_2_scalar^2+sd_3_scalar^2+var_b+sigma^2;",
                      block = "genquant") + 
        brms::stanvar(scode = "real ccc = ccc_num/ccc_denom;", block = "genquant")
    }else if(distr_family == "skew_normal"){
      stanvars_mod <- stanvars_mod +
        brms::stanvar(scode = "real ccc_num = sd_1_scalar^2+sd_3_scalar^2;", block = "genquant") + 
        brms::stanvar(scode = "real ccc_denom = sd_1_scalar^2+sd_2_scalar^2+sd_3_scalar^2+var_b+sigma^2;",
                      block = "genquant") + 
        brms::stanvar(scode = "real ccc = ccc_num/ccc_denom;", block = "genquant")
    }else{
      stanvars_mod <- stanvars_mod + 
        brms::stanvar(scode = "real term_exp = exp(sd_1_scalar^2);", block = "genquant") +
        brms::stanvar(scode = "real var_a = fmin(exp(2*logn_mu)*term_exp*(term_exp-1),1e10);",
                      block = "genquant") + 
        brms::stanvar(scode = "real ccc_num = var_a+sd_3_scalar^2;", block = "genquant") +  
        brms::stanvar(scode = "real ccc_denom = var_a+sd_2_scalar^2+sd_3_scalar^2+var_b+sigma^2;",
                      block = "genquant") + 
        brms::stanvar(scode = "real ccc = ccc_num/ccc_denom;", block = "genquant")
    }
    
  }else{
    
    if(is.null(covar)){
      stanvars_mod <- stanvars_mod + brms::stanvar(scode = "matrix[K, (M_2-1)] L;",
                                                   block = "data", x = L) +
        brms::stanvar(scode = "matrix[(M_2-1), (M_2-1)] auxD;",
                      block = "data", x = auxD) + 
        brms::stanvar(scode = "\n  vector[Kc + 1] fixed_effects;\n  fixed_effects[1] = b_Intercept;\n  for (k in 1:Kc) {\n    fixed_effects[k + 1] = b[k];\n  }", block = "genquant") + 
        brms::stanvar(scode = "vector[M_2*(M_2-1)/2] difmed = transpose(L) * fixed_effects;", block = "genquant") + 
        brms::stanvar(scode = "real aux1 = (transpose(difmed) * difmed); ", block = "genquant") + 
        brms::stanvar(scode = "real var_b = fmax(aux1/(M_2*(M_2-1)),0);", block = "genquant")
    }else{
      stanvars_mod <- stanvars_mod + brms::stanvar(scode = "matrix[M_2*M_3, (M_2-1)] L;",
                                                   block = "data", x = L) +
        brms::stanvar(scode = "matrix[(M_2-1), (M_2-1)] auxD;",
                      block = "data", x = auxD) + 
        brms::stanvar(scode = "\n  vector[M_2*M_3] fixed_effects;\n  fixed_effects[1] = b_Intercept;\n  for (k in 1:(M_2*M_3)) {\n    fixed_effects[k + 1] = b[k];\n  }", block = "genquant") + 
        brms::stanvar(scode = "vector[M_2*(M_2-1)/2] difmed = transpose(L) * fixed_effects;", block = "genquant") + 
        brms::stanvar(scode = "real aux1 = (transpose(difmed) * difmed); ", block = "genquant") + 
        brms::stanvar(scode = "real var_b = fmax(aux1/(M_2*(M_2-1)),0);", block = "genquant")
    }

    
    if(int){
      if(!log_alpha){
        stanvars_mod <- stanvars_mod + 
          brms::stanvar(scode = "real ccc_num = sd_1_scalar^2;", block = "genquant") +  
          brms::stanvar(scode = "real ccc_denom = sd_1_scalar^2+sd_2_scalar^2+var_b+sigma^2;",
                        block = "genquant") + 
          brms::stanvar(scode = "real ccc = ccc_num/ccc_denom;", block = "genquant")
      }else{
        stanvars_mod <- stanvars_mod + 
          brms::stanvar(scode = "real term_exp = exp(sd_1_scalar^2);", block = "genquant") +
          brms::stanvar(scode = "real var_a = fmin(exp(2*logn_mu)*term_exp*(term_exp-1),1e10);",
                        block = "genquant") + 
          brms::stanvar(scode = "real ccc_num = var_a;", block = "genquant") +  
          brms::stanvar(scode = "real ccc_denom = var_a+sd_2_scalar^2+var_b+sigma^2;",
                        block = "genquant") + 
          brms::stanvar(scode = "real ccc = ccc_num/ccc_denom;", block = "genquant")
      }
    }else{
      if(!log_alpha){
        stanvars_mod <- stanvars_mod + 
          brms::stanvar(scode = "real ccc_num = sd_1_scalar^2;", block = "genquant") +  
          brms::stanvar(scode = "real ccc_denom = sd_1_scalar^2+var_b+sigma^2;",
                        block = "genquant") + 
          brms::stanvar(scode = "real ccc = ccc_num/ccc_denom;", block = "genquant")
      }else{
        stanvars_mod <- stanvars_mod + 
          brms::stanvar(scode = "real term_exp = exp(sd_1_scalar^2);", block = "genquant") +
          brms::stanvar(scode = "real var_a = fmin(exp(2*logn_mu)*term_exp*(term_exp-1),1e10);",
                        block = "genquant") + 
          brms::stanvar(scode = "real ccc_num = var_a;", block = "genquant") +  
          brms::stanvar(scode = "real ccc_denom = var_a+var_b+sigma^2;",
                        block = "genquant") + 
          brms::stanvar(scode = "real ccc = ccc_num/ccc_denom;", block = "genquant")
      }
    }
    
    
  }
  
  # Lognormal case
  
  if(log_alpha){
    # priors = brms::empty_prior()
    # priors <- priors + brms::prior(student_t(3, 0, 1.5),class = "logn", coef = mu)
    priors <- priors + brms::prior(gamma(1.5,2),class = "sd", group = ind,coef = Intercept)
    
    stanvars_mod <- stanvars_mod +
      brms::stanvar(scode = "real logn_mu;",block = "parameters") +
      brms::stanvar(scode = "lprior += student_t_lpdf(logn_mu | 3, 0, 1.5);", block = "tparameters")
  }
  
  # Gen code
  mod_stan <- brms::stancode(form,data = dades,family = distr_family, stanvars = stanvars_mod, prior = priors)
  
  
  # Modify parameters
  mod_stan <- stringr::str_replace(mod_stan,"vector<lower=0>\\[M_1\\] sd_1;","real<lower=0> sd_1_scalar;")
  if((!is.null(rmet)) & int){
    mod_stan <- stringr::str_replace(mod_stan,"vector<lower=0>\\[M_2\\] sd_2;","real<lower=0> sd_2_scalar;")
    if(!is.null(rtime)){
      mod_stan <- stringr::str_replace(mod_stan,"vector<lower=0>\\[M_3\\] sd_3;","real<lower=0> sd_3_scalar;")
    }
  }
  
  if(log_alpha){
    mod_stan <- stringr::str_replace(mod_stan,"r_1_1 = \\(","r_1_1 = exp\\(logn_mu + ")
  }
  
  # Regen code
  #mod_stan <- rstan::stan_model(model_code = mod_stan)
  
  # Model
  #-------------------------------------------------------
  
  model.stan <- rstan::stan(model_code = mod_stan, data = dades_stan,
                            iter = n_iter, chains = n_chains, warmup = n_burnin, cores = workers)
  
  # Estimation
  #-------------------------------------------------------
  
  # Extract the VC
  pars_ccc <- c("var_a","sd_1_scalar","sd_2_scalar","sd_3_scalar","var_b","sigma")
  
  if(!log_alpha){
    pars_ccc <- pars_ccc[pars_ccc != "var_a"]
  }else{
    pars_ccc <- pars_ccc[pars_ccc != "sd_1_scalar"]
  }
  
  if(!int){
    pars_ccc <- pars_ccc[pars_ccc != "sd_2_scalar"]
  }
  
  if(is.null(rtime)){
    pars_ccc <- pars_ccc[pars_ccc != "sd_3_scalar"]
  }
  
  vars_ccc <- rstan::extract(model.stan, pars = pars_ccc)
  
  ccc_bay_vc <- coda:::summary.mcmc.list(rstan::As.mcmc.list(model.stan, 
                                                             pars = pars_ccc))$statistics[,1]
  
  if(plot.gel){
    aux_eff <- coda::effectiveSize(rstan::As.mcmc.list(model.stan, pars = "ccc")) |> round(4)
    coda::gelman.plot(rstan::As.mcmc.list(model.stan, pars = "ccc"))
    legend("topleft", legend = paste0("ESS (n):\n",aux_eff), bty = "n")
  }
  
  if(!correct){
    ccc_bay <- rstan::As.mcmc.list(model.stan, pars = "ccc")
    ccc_bay_ci <-HDInterval::hdi(ccc_bay)
    ccc_bay <- c(bayestestR::map_estimate(ccc_bay)$MAP_Estimate,
                 coda:::summary.mcmc.list(ccc_bay)$statistics[2],
                 coda:::summary.mcmc.list(ccc_bay)$quantiles[3])
    
  }else{
    # Correct the var_b estimate
    Sb <- cov(as.data.frame(rstan::extract(model.stan, pars = c("fixed_effects"))))
    aux_betas <- rstan::summary(model.stan, pars = c("Intercept","b"), probs = c(0.1, 0.9))$summary
    b <- aux_betas[,1]
    difmed <- t(L)%*%b
    A<-L%*%auxD%*%t(L)
    if(!is.null(rtime)){
      trace_varb <- sum(diag((A%*%Sb)))/(nm*(nm-1)*nt)
    }else{
      trace_varb <- sum(diag((A%*%Sb)))/(nm*(nm-1))
    }
    
    vars_ccc$var_b <- sapply(vars_ccc$var_b, function(i) max(i-trace_varb,0))
    
    # Compute corrected ccc
    if(!log_alpha){
      if(!is.null(rtime)){
        aux_vars <- vars_ccc |> as.data.frame() |>
          dplyr::mutate(ccc = (sd_1_scalar^2 + sd_3_scalar^2)/(sd_1_scalar^2 + sd_2_scalar^2 + sd_3_scalar^2 + var_b + sigma^2))
      }else if(int){
        aux_vars <- vars_ccc |> as.data.frame() |> 
          dplyr::mutate(ccc = (sd_1_scalar^2)/(sd_1_scalar^2 + sd_2_scalar^2 + var_b + sigma^2))
      }else{
        aux_vars <- vars_ccc |> as.data.frame() |> 
          dplyr::mutate(ccc = (sd_1_scalar^2)/(sd_1_scalar^2 + sd_2_scalar^2 + var_b + sigma^2))
      }
    }else{
      if(!is.null(rtime)){
        aux_vars <- vars_ccc |> as.data.frame() |>
          dplyr::mutate(ccc = (var_a + sd_3_scalar^2)/(var_a + sd_2_scalar^2 + sd_3_scalar^2 + var_b + sigma^2))
      }else if(int){
        aux_vars <- vars_ccc |> as.data.frame() |> 
          dplyr::mutate(ccc = var_a/(var_a + sd_2_scalar^2 + var_b + sigma^2))
      }else{
        aux_vars <- vars_ccc |> as.data.frame() |> 
          dplyr::mutate(ccc = var_a/(var_a + sd_2_scalar^2 + var_b + sigma^2))
      }
    }
    
    ccc_bay <- c(bayestestR::map_estimate(aux_vars$ccc)$MAP_Estimate,
                 sd(coda::as.mcmc(aux_vars$ccc)),
                 median(coda::as.mcmc(aux_vars$ccc)))
    ccc_bay_ci <- coda::HPDinterval(coda::as.mcmc(aux_vars$ccc))
    
    ccc_bay_vc["var_b"] <- mean(vars_ccc$var_b)
  }
  
  res <- c(ccc_bay[1],ccc_bay_ci,ccc_bay[2],ccc_bay[3])
  conf.lab <- paste(cl*100,"%",sep="")
  names(res) <- c("CCC",paste("LL CI",conf.lab),
                  paste("UL CI",conf.lab),"SE CCC","Median CCC")
  response <- list(ccc=res,vc = ccc_bay_vc, sigma = NULL, model = model.stan,
                   resamples=NULL,transf=NULL,m=NULL,N=NULL)
  
  class(response)<-"ccc"
  return(response)
}