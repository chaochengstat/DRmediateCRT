

#' Estimation of the marginal mediator model
est_m_marginal = function(formula.M,data) {
  m_model = lm(formula.M,data)
  m_coef = m_model$coef
  m_sigma  = sd(m_model$residuals)
  epsilon = qnorm(pnorm(m_model$model[,1],mean=m_model$fitted.values,sd=m_sigma))
  list(m_coef=m_coef,m_sigma=m_sigma,epsilon=epsilon)
}

#' Estimation of the Gaussian copula model with the exchangable correlation matrix
est_copula = function(m_marginal_model,data,id_name="id",N_name="N") {
  data_all = data.frame(data,epsilon=m_marginal_model$epsilon)
  K = unique(data_all[,id_name])
  fn = function (par) {
    loglik = 0;i=1
    while (i<dim(data_all)[1]) {
      j=i-1+data_all[i,N_name]
      N = data_all[i,N_name]
      x = data_all$epsilon[i:j]
      rho = par[1]
      Sigma_inv = 1/(1-rho)*(diag(N)-rho/(1+(N-1)*rho)*(rep(1,N) %*% t(rep(1,N))))
      Sigma_det = ((1-rho)^N)*(1+rho/(1-rho)*N)
      loglik = loglik + c(t(x) %*% Sigma_inv %*% x + log(Sigma_det))
      i=j+1
    }
    loglik
  }
  rho = optim(par=0.5,fn=fn,method="L-BFGS-B",lower=0.005,upper=0.995)
  return(c(rho=rho[[1]]))
}


#' Estimation of the marginal outcome model
est_y_marginal = function(formula.Y,data,A_level=1,A_name="A") {
  data_sub = data[which(data[,A_name]==A_level),]
  y_model = lm(formula.Y,data_sub)
  y_coef = y_model$coef
  list(y_coef=y_coef,se=summary(y_model)$coefficients)
}


#' Estimation of the mediation functional theta(a,a*) based on the semiparametric doubly robust estimator under original EIFs
#'
#' @param data the dataset
#' @param a1 the value *a* in *\theta(a,a^\star)*
#' @param a0 the value *a^\star* in *\theta(a,a^\star)*
#' @param formula.M formula for the mediator model
#' @param formula.Y formula for the outcome model
#' @param m_res regression output from the marginal mediator model
#' @param y1_res regression output from the marginal outcome model under the treated clusters
#' @param y0_res regression output from the marginal outcome model under the control clusters
#' @param copula_res output from the copula model
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns The estimates of theta(a,a*)
MF_theta_par_est1 = function(data,a1,a0,formula.M,formula.Y,m_res,y1_res,y0_res,copula_res,is_stab=1) {
  #a1=0;a0=0;is_stab=1
  data_sum = unique(data[,c(id_name,A_name)])
  a_prob = 0.5
  # number of clusters
  K = length(unique(data[,id_name]))
  # model matrix for mediator
  data0 = data1 = data; data0[,A_name]=0;data1[,A_name]=1;
  data_m0 =model.matrix(formula.M,data0)
  data_m1 =model.matrix(formula.M,data1)
  # model matrix for outcome
  data_y = model.matrix(formula.Y,data)
  k=1
  # calculate the mediator weights (and the stablizing weights)
  Sigma_M.list = mean_m0_marginal.list = mean_m1_marginal.list =
    sd_m0_marginal.list = sd_m1_marginal.list = f_m0c.list = f_m1c.list = as.list(rep(1,K))
  Stab_C = Stab_I = rep(1,dim(data)[1])
  for (k in (1:K)) {
    if (k==1) {ind = 1:data[1,N_name]}
    if (k>1) {ind = (max(ind)+1):(max(ind)+data[max(ind)+1,N_name])}
    N=length(ind)
    A = data[ind[1],A_name]
    # f(M|0,C)
    Sigma_M = (1-copula_res)*diag(N) + copula_res*(rep(1,N) %*% t(rep(1,N)))
    mean_m0_marginal = c(data_m0[ind,,drop=F] %*% m_res$m_coef)
    sd_m0_marginal = m_res$m_sigma
    epsilon_m0 = qnorm(pnorm(data[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal))
    f_m0_marginal = dnorm(data[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal)
    f_epsilon_m0 = dnorm(epsilon_m0)
    copula_m0 = mvtnorm::dmvnorm(epsilon_m0,mean=rep(0,N),sigma=Sigma_M)
    f_m0c = copula_m0*prod(f_m0_marginal/f_epsilon_m0)
    # f(M|1,C)
    mean_m1_marginal = c(data_m1[ind,,drop=F] %*% m_res$m_coef)
    sd_m1_marginal = m_res$m_sigma
    epsilon_m1 = qnorm(pnorm(data[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal))
    f_m1_marginal = dnorm(data[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal)
    f_epsilon_m1 = dnorm(epsilon_m1)
    copula_m1 = mvtnorm::dmvnorm(epsilon_m1,mean=rep(0,N),sigma=Sigma_M)
    f_m1c = copula_m1*prod(f_m1_marginal/f_epsilon_m1)
    # write the lists
    Sigma_M.list[[k]] = Sigma_M
    mean_m0_marginal.list[[k]] = mean_m0_marginal
    mean_m1_marginal.list[[k]] = mean_m1_marginal
    sd_m0_marginal.list[[k]] = sd_m0_marginal
    sd_m1_marginal.list[[k]] = sd_m1_marginal
    f_m0c.list[[k]] = f_m0c
    f_m1c.list[[k]] = f_m1c
    Stab_C[ind] = (1/a_prob)*(f_m0c/f_m1c)*(1/N)
    Stab_I[ind] = (1/a_prob)*(f_m0c/f_m1c)
    #print((1/a_prob)*(f_m0c/f_m1c)*(1/N))
  }
  if (is_stab==0) Stab_C = Stab_I = rep(1,dim(data)[1])
  if (a1==a0) Stab_C = Stab_I = rep(1,dim(data)[1])
  # outcome mean
  y_model_C = y_model_I = list()

  X_y = model.matrix(formula.Y,data[which(data[,A_name]==a1),])
  Y_y = data[which(data[,A_name]==a1),Y_name]
  weights_y = sqrt(Stab_C[which(data[,A_name]==a1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  y_model_C$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)

  X_y = model.matrix(formula.Y,data[which(data[,A_name]==a1),])
  Y_y = data[which(data[,A_name]==a1),Y_name]
  weights_y = sqrt(Stab_I[which(data[,A_name]==a1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  y_model_I$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)

  # calculate the EIF
  k=1
  B_monte_carlo=25
  I_eif = I_eif_ind = 0
  I_g = I_g_ind = 0
  N.list=c()
  for (k in (1:K)) {
    if (k==1) {ind = 1:data[1,N_name]}
    if (k>1) {ind = (max(ind)+1):(max(ind)+data[max(ind)+1,N_name])}
    N=length(ind);N.list=c(N.list,N)
    A = data[ind[1],A_name]
    # ###### 1st row of EIF
    # Y_bar
    Y_bar = mean(data[ind,Y_name])
    # model-based outcome mean E[Y_bar|1,M,C]
    Y_bar_model = mean(c(data_y[ind,,drop=F] %*% y_model_C$y_coef))
    # calculate I1
    if (a1 != a0) {
      I1 = (Y_bar - Y_bar_model)*(f_m0c.list[[k]]/f_m1c.list[[k]])*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
    } else {
      I1 = (Y_bar - Y_bar_model)*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
    }
    ###### 2nd row of EIF
    # eta
    copula_M = pnorm(mvtnorm::rmvnorm(B_monte_carlo,mean=rep(0,N),sigma=Sigma_M.list[[k]]))
    if (a0==0) {
      M_sim  = qnorm(copula_M, mean= rep(1,B_monte_carlo) %*% t(mean_m0_marginal.list[[k]]),sd=sd_m0_marginal.list[[k]])
    } else {
      M_sim  = qnorm(copula_M, mean= rep(1,B_monte_carlo) %*% t(mean_m1_marginal.list[[k]]),sd=sd_m1_marginal.list[[k]])
    }
    eta.list = rep(0,B_monte_carlo)
    data_now = data_y[ind,]
    data_name = colnames(data_now)
    for ( l in (1:B_monte_carlo)) {
      if (M_name %in% data_name ) data_now[,M_name] = M_sim[l,];
      if (M_bar_name %in% data_name ) data_now[,M_bar_name] = (sum(M_sim[l,])-M_sim[l,])/(N-1)
      eta.list[l] = mean(c(data_now %*% y_model_C$y_coef))
    }
    eta = mean(eta.list)
    I2 = (Y_bar_model-eta)*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
    # 3rd row of EIF
    I3 = eta
    # EIF-based estimator
    I_eif = I_eif + (I1+I2+I3)
    I_g = I_g + I3

    ###### 1st row of EIF
    # Y_bar
    Y_bar = mean(data[ind,Y_name])
    # model-based outcome mean E[Y_bar|1,M,C]
    Y_bar_model = mean(c(data_y[ind,,drop=F] %*% y_model_I$y_coef))
    # calculate I1
    if (a1 != a0) {
      I1 = (Y_bar - Y_bar_model)*(f_m0c.list[[k]]/f_m1c.list[[k]])*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
    } else {
      I1 = (Y_bar - Y_bar_model)*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
    }
    ###### 2nd row of EIF
    # eta
    eta.list = rep(0,B_monte_carlo)
    data_now = data_y[ind,]
    for ( l in (1:B_monte_carlo)) {
      if (M_name %in% data_name) data_now[,M_name] = M_sim[l,];
      if (M_bar_name %in% data_name) data_now[,M_bar_name] = (sum(M_sim[l,])-M_sim[l,])/(N-1)
      eta.list[l] = mean(c(data_now %*% y_model_I$y_coef))
    }
    eta = mean(eta.list)
    I2 = (Y_bar_model-eta)*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
    # 3rd row of EIF
    I3 = eta
    # EIF-based estimator
    I_eif_ind = I_eif_ind + N*(I1+I2+I3)
    I_g_ind = I_g_ind + I3*N
  }

  list(Cluster=c(m1=I_g/K,eif=I_eif/K),
       Individual=c(m1=I_g_ind/K,eif=I_eif_ind/K)/mean(N.list))
}


#' Estimation of the mediation functional theta(a,a*) based on the machine learning estimator under original EIFs
#'
#' @param data the dataset
#' @param a1 the value *a* in *\theta(a,a^\star)*
#' @param a0 the value *a^\star* in *\theta(a,a^\star)*
#' @param formula.M formula for the mediator model
#' @param formula.Y formula for the outcome model
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param V number of groups based on the cross-fitting procedure
#'
#' @returns The estimates of theta(a,a*)
MF_theta_nonpar_est1 = function(data,a1=1,a0=0,formula.M,formula.Y,is_stab=1,V=5,is_true_copula=1,myseed=2023) {
  #a1=1;a0=0;V=5;myseed=2023;is_stab=1
  # create a dataset including all covariates (including the interaction terms)
  data1 = data0 = data; data1[,A_name] = 1; data0[,A_name] = 0
  data_m = model.matrix(formula.M,data)[,-1];data_m1 = model.matrix(formula.M,data1)[,-1];data_m0 = model.matrix(formula.M,data0)[,-1]
  additional_columns = setdiff(colnames(data_m),colnames(data))
  if (length(additional_columns)>0) {
    data = cbind(data,data_m[,additional_columns,drop=FALSE])
    data1 = cbind(data1,data_m1[,additional_columns,drop=FALSE])
    data0 = cbind(data0,data_m0[,additional_columns,drop=FALSE])
    colnames(data) = colnames(data1)=colnames(data0)=gsub("[(*) ]","",colnames(data))
  }
  # names of the mediator and outcome models
  M_model_names = gsub("[(*) ]","",colnames(data_m))
  Y_model_names = gsub("[(*) ]","",colnames(model.matrix(formula.Y,data))[-1])
  # treatment probability
  data_sum = aggregate(data, by = list(data[,id_name]), mean)[,-1]
  N.list = data_sum[,N_name]
  a_prob = 0.5
  # number of clusters
  K = dim(data_sum)[1]
  # initialize influence function values
  I_eif_C = I_eif_I =rep(0,K)
  I_g_C = I_g_I =rep(0,K)
  set.seed(myseed)
  folds = createFolds(1:K, k = V)
  v=1
  Sigma_M.list = mean_m0_marginal.list = mean_m1_marginal.list =
    sd_m0_marginal.list = sd_m1_marginal.list = f_m0c.list = f_m1c.list = as.list(rep(1,K))
  for (v in (1:V)) {
    index_main = -which(data[,id_name] %in% folds[[v]])
    index_vali = which(data[,id_name] %in% folds[[v]])
    data_main=data[index_main,]; data_main1=data1[index_main,]; data_main0=data0[index_main,]
    data_vali=data[index_vali,]; data_vali1=data1[index_vali,]; data_vali0=data0[index_vali,]
    data_sum_main = data_sum[-which(data_sum[,id_name] %in% folds[[v]]),]
    data_sum_vali = data_sum[which(data_sum[,id_name] %in% folds[[v]]),]
    # mediator model (coef & sd of error term)
    m_m = SuperLearner(
      Y          = data_main[,M_name],
      X          = data_main[,M_model_names,drop=F],
      family     = gaussian(),
      id = data_main[,id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    #m_m$epsilon = data_main[,M_name] - predict(m_m, onlySL = TRUE)$pred[,1]
    m_m$m_sigma = sqrt(mean((data_main[,M_name] - predict(m_m, onlySL = TRUE)$pred[,1])^2))
    m_m$epsilon = qnorm(pnorm(data_main[,M_name],mean=predict(m_m, onlySL = TRUE)$pred[,1],sd=m_m$m_sigma))
    # copula estimate
    if (is_true_copula==1)  copula_res = est_copula(m_m,data_main,id_name="id",N_name="N")
    if (is_true_copula==0)  copula_res = 0.0001
    #print(copula_res)
    # outcome model
    m_y = SuperLearner(
      Y          = data_main[which(data_main[,A_name]==a1),Y_name],
      X          = data_main[which(data_main[,A_name]==a1),Y_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==a1),id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    # prediction for the validation data
    Stab_C = Stab_I = rep(1,dim(data_vali)[1])
    for (k in folds[[v]]) {
      ind = which(data_vali[,id_name]==k)
      N=length(ind)
      A = data_vali[ind[1],A_name]
      # f(M|0,C)
      Sigma_M = (1-copula_res)*diag(N) + copula_res*(rep(1,N) %*% t(rep(1,N)))
      mean_m0_marginal = predict(m_m, newdata = data_vali0[ind,], onlySL = TRUE)$pred[,1]
      sd_m0_marginal = m_m$m_sigma
      epsilon_m0 = qnorm(pnorm(data_vali[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal))
      f_m0_marginal = dnorm(data_vali[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal)
      f_epsilon_m0 = dnorm(epsilon_m0)
      copula_m0 = mvtnorm::dmvnorm(epsilon_m0,mean=rep(0,N),sigma=Sigma_M)
      f_m0c = copula_m0*prod(f_m0_marginal/f_epsilon_m0)
      # f(M|1,C)
      mean_m1_marginal = predict(m_m, newdata = data_vali1[ind,], onlySL = TRUE)$pred[,1]
      sd_m1_marginal = m_m$m_sigma
      epsilon_m1 = qnorm(pnorm(data_vali[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal))
      f_m1_marginal = dnorm(data_vali[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal)
      f_epsilon_m1 = dnorm(epsilon_m1)
      copula_m1 = mvtnorm::dmvnorm(epsilon_m1,mean=rep(0,N),sigma=Sigma_M)
      f_m1c = copula_m1*prod(f_m1_marginal/f_epsilon_m1)
      # write the lists
      Sigma_M.list[[k]] = Sigma_M
      mean_m0_marginal.list[[k]] = mean_m0_marginal
      mean_m1_marginal.list[[k]] = mean_m1_marginal
      sd_m0_marginal.list[[k]] = sd_m0_marginal
      sd_m1_marginal.list[[k]] = sd_m1_marginal
      f_m0c.list[[k]] = f_m0c
      f_m1c.list[[k]] = f_m1c
      Stab_C[ind] = (1/a_prob)*(f_m0c/f_m1c)*(1/N)
      Stab_I[ind] = (1/a_prob)*(f_m0c/f_m1c)
    }
    # TMLE if (is_stab = 1)
    Y_1step = predict(m_y, newdata = data_vali, onlySL = TRUE)$pred[,1]
    if (is_stab == 1) {
      est_y_tmle_C = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==a1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==a1),Y_name],
                             weights=Stab_C[which(data_vali$A==a1)],
                             offset=Y_1step[which(data_vali$A==a1)],
                             family=gaussian(link = "identity"))$coef
      est_y_tmle_I = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==a1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==a1),Y_name],
                             weights=Stab_I[which(data_vali$A==a1)],
                             offset=Y_1step[which(data_vali$A==a1)],
                             family=gaussian(link = "identity"))$coef
    } else {
      est_y_tmle_C = est_y_tmle_I = 0
    }
    Y_final_C = Y_1step + est_y_tmle_C
    Y_final_I = Y_1step + est_y_tmle_I
    B_monte_carlo=25
    for (k in folds[[v]]) {
      ind = which(data_vali[,id_name]==k)
      N = length(ind)
      A = data_vali[ind[1],A_name]
      # ###### 1st row of EIF
      # Y_bar
      Y_bar = mean(data_vali[ind,Y_name])
      # model-based outcome mean E[Y_bar|1,M,C]
      Y_bar_model_C = mean(Y_final_C[ind])
      Y_bar_model_I = mean(Y_final_I[ind])
      # calculate I1
      if (a1 != a0) {
        I1_C = (Y_bar - Y_bar_model_C)*(f_m0c.list[[k]]/f_m1c.list[[k]])*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
        I1_I = (Y_bar - Y_bar_model_I)*(f_m0c.list[[k]]/f_m1c.list[[k]])*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
      } else {
        I1_C = (Y_bar - Y_bar_model_C)*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
        I1_I = (Y_bar - Y_bar_model_I)*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
      }
      ###### 2nd row of EIF
      # eta
      copula_M = pnorm(mvtnorm::rmvnorm(B_monte_carlo,mean=rep(0,N),sigma=Sigma_M.list[[k]]))
      if (a0==0) {
        M_sim  = qnorm(copula_M, mean= rep(1,B_monte_carlo) %*% t(mean_m0_marginal.list[[k]]),sd=sd_m0_marginal.list[[k]])
      } else {
        M_sim  = qnorm(copula_M, mean= rep(1,B_monte_carlo) %*% t(mean_m1_marginal.list[[k]]),sd=sd_m1_marginal.list[[k]])
      }
      eta.list_C = eta.list_I = rep(0,B_monte_carlo)
      data_now = data_vali[ind,]
      for ( l in (1:B_monte_carlo)) {
        data_now[,M_name] = M_sim[l,]; data_now[,M_bar_name] = (sum(M_sim[l,])-M_sim[l,])/(N-1)
        eta.list_C[l] = mean(est_y_tmle_C+predict(m_y,newdata = data_now, onlySL = TRUE)$pred[,1])
        eta.list_I[l] = mean(est_y_tmle_I+predict(m_y,newdata = data_now, onlySL = TRUE)$pred[,1])
      }
      eta_C = mean(eta.list_C)
      eta_I = mean(eta.list_I)
      I2_C = (Y_bar_model_C-eta_C)*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
      I2_I = (Y_bar_model_I-eta_I)*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
      # 3rd row of EIF
      I3_C = eta_C
      I3_I = eta_I
      # EIF-based estimator
      I_eif_C[k] = (I1_C+I2_C+I3_C)
      I_g_C[k] = I3_C
      I_eif_I[k] = (I1_I+I2_I+I3_I)*N
      I_g_I[k] = I3_I*N
    }
  }
  list(Cluster=c(m1=sum(I_g_C)/K,eif=sum(I_eif_C)/K),
       Individual=c(m1=sum(I_g_I)/K,eif=sum(I_eif_I)/K)/mean(N.list),
       I_eif_C = I_eif_C,
       I_eif_I = I_eif_I)
}



#' Estimation of the mediation functional tau based on the semiparametric doubly robust estimator under original EIFs
#'
#' @param data the dataset
#' @param formula.M formula for the mediator model
#' @param formula.Y formula for the outcome model
#' @param m_res regression output from the marginal mediator model
#' @param y1_res regression output from the marginal outcome model under the treated clusters
#' @param y0_res regression output from the marginal outcome model under the control clusters
#' @param copula_res output from the copula model
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns The estimates of tau
MF_tau_par_est1 = function(data,formula.M,formula.Y,m_res,y1_res,y0_res,copula_res,is_stab=1) {
  # treatment probability
  data_sum = unique(data[,c(id_name,A_name)])
  a_prob = 0.5
  # number of clusters
  K = length(unique(data[,id_name]))
  # model matrix for mediator
  data0 = data1 = data; data0[,A_name]=0;data1[,A_name]=1;
  data_m0 =model.matrix(formula.M,data0)
  data_m1 =model.matrix(formula.M,data1)
  # model matrix for outcome
  data_y = model.matrix(formula.Y,data1)

  k=1
  B_monte_carlo=25
  Sigma_M.list = mean_m0_marginal.list = mean_m1_marginal.list =
    sd_m0_marginal.list = sd_m1_marginal.list = f_m_noj_0c.list = f_m1_marginal.list = f_m1c.list = as.list(rep(1,K))
  Stab_C = Stab_I = rep(1,dim(data)[1])
  for (k in (1:K)) {
    if (k==1) {ind = 1:data[1,N_name]}
    if (k>1) {ind = (max(ind)+1):(max(ind)+data[max(ind)+1,N_name])}
    N=length(ind)
    A = data[ind[1],A_name]

    # f(M|0,C)
    Sigma_M = (1-copula_res)*diag(N) + copula_res*(rep(1,N) %*% t(rep(1,N)))
    mean_m0_marginal = c(data_m0[ind,,drop=F] %*% m_res$m_coef)
    sd_m0_marginal = m_res$m_sigma
    epsilon_m0 = qnorm(pnorm(data[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal))
    f_m0_marginal = dnorm(data[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal)
    f_epsilon_m0 = dnorm(epsilon_m0)
    f_m_noj_0c=c()
    for (j in (1:N)) {
      copula_m_noj_0 = mvtnorm::dmvnorm(epsilon_m0[-j],mean=rep(0,N-1),sigma=Sigma_M[1:(N-1),1:(N-1)])
      f_m_noj_0c[j] = copula_m_noj_0*prod(f_m0_marginal[-j]/f_epsilon_m0[-j])
    }
    # f(M|1,C)
    mean_m1_marginal = c(data_m1[ind,,drop=F] %*% m_res$m_coef)
    sd_m1_marginal = m_res$m_sigma
    epsilon_m1 = qnorm(pnorm(data[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal))
    f_m1_marginal = dnorm(data[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal)
    f_epsilon_m1 = dnorm(epsilon_m1)
    copula_m1 = mvtnorm::dmvnorm(epsilon_m1,mean=rep(0,N),sigma=Sigma_M)
    f_m1c = copula_m1*prod(f_m1_marginal/f_epsilon_m1)

    # write the lists
    Sigma_M.list[[k]] = Sigma_M
    mean_m0_marginal.list[[k]] = mean_m0_marginal
    mean_m1_marginal.list[[k]] = mean_m1_marginal
    sd_m0_marginal.list[[k]] = sd_m0_marginal
    sd_m1_marginal.list[[k]] = sd_m1_marginal
    f_m_noj_0c.list[[k]] = f_m_noj_0c
    f_m1_marginal.list[[k]] = f_m1_marginal
    f_m1c.list[[k]] = f_m1c
    Stab_C[ind] = (1/a_prob)*(f_m1_marginal*f_m_noj_0c/f_m1c)*(1/N)
    Stab_I[ind] = (1/a_prob)*(f_m1_marginal*f_m_noj_0c/f_m1c)
  }
  if (is_stab==0) Stab_C = Stab_I = rep(1,dim(data)[1])

  # outcome mean
  y_model_C = y_model_I = list()
  X_y = model.matrix(formula.Y,data[which(data[,A_name]==1),])
  Y_y = data[which(data[,A_name]==1),Y_name]
  weights_y = sqrt(Stab_C[which(data[,A_name]==1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  y_model_C$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)

  X_y = model.matrix(formula.Y,data[which(data[,A_name]==1),])
  Y_y = data[which(data[,A_name]==1),Y_name]
  weights_y = sqrt(Stab_I[which(data[,A_name]==1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  y_model_I$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)


  k=1
  B_monte_carlo=25
  I_eif_C = I_eif_I = 0
  I_g_C = I_g_I = 0
  N.list=c()
  for (k in (1:K)) {
    if (k==1) {ind = 1:data[1,N_name]}
    if (k>1) {ind = (max(ind)+1):(max(ind)+data[max(ind)+1,N_name])}
    N=length(ind);N.list=c(N.list,N)
    A = data[ind[1],A_name]
    ###### construct the EIF
    # Y_bar
    Y = data[ind,Y_name]
    M = data[ind,M_name]
    # model-based outcome mean E[Y_bar|1,M,C]
    Y_bar_model_C = c(data_y[ind,,drop=F] %*% y_model_C$y_coef)
    Y_bar_model_I = c(data_y[ind,,drop=F] %*% y_model_I$y_coef)
    # calculate the intergals
    copula_M = pnorm(mvtnorm::rmvnorm(B_monte_carlo,mean=rep(0,N),sigma=Sigma_M.list[[k]]))
    M_sim_0  = qnorm(copula_M,mean= rep(1,B_monte_carlo) %*% t(mean_m0_marginal.list[[k]]),sd=sd_m0_marginal.list[[k]])
    copula_M = pnorm(mvtnorm::rmvnorm(B_monte_carlo,mean=rep(0,N),sigma=Sigma_M.list[[k]]))
    M_sim_1  = qnorm(copula_M,mean= rep(1,B_monte_carlo) %*% t(mean_m1_marginal.list[[k]]),sd=sd_m1_marginal.list[[k]])
    eta_C=kappa_C=xi_C = matrix(NA, nrow = N, ncol = B_monte_carlo)
    eta_I=kappa_I=xi_I = matrix(NA, nrow = N, ncol = B_monte_carlo)
    data_now1 = data_now2 = data_now3 = data_y[ind,]
    data_name = colnames(data_y)
    for ( l in (1:B_monte_carlo)) {
      if (M_name %in% data_name ) data_now1[,M_name] = M_sim_1[l,];
      if (M_bar_name %in% data_name ) data_now1[,M_bar_name] = (sum(M_sim_0[l,])-M_sim_0[l,])/(N-1)
      eta_C[,l] = c(data_now1 %*% y_model_C$y_coef)
      eta_I[,l] = c(data_now1 %*% y_model_I$y_coef)
      if (M_bar_name %in% data_name ) data_now2[,M_bar_name] = (sum(M_sim_0[l,])-M_sim_0[l,])/(N-1)
      xi_C[,l] = c(data_now2 %*% y_model_C$y_coef)
      xi_I[,l] = c(data_now2 %*% y_model_I$y_coef)
      if (M_name %in% data_name ) data_now3[,M_name] = M_sim_1[l,]
      kappa_C[,l] = c(data_now3 %*% y_model_C$y_coef)
      kappa_I[,l] = c(data_now3 %*% y_model_I$y_coef)
    }
    eta_C = apply(eta_C,1,mean)
    eta_I = apply(eta_I,1,mean)
    xi_C = apply(xi_C,1,mean)
    xi_I = apply(xi_I,1,mean)
    kappa_C = apply(kappa_C,1,mean)
    kappa_I = apply(kappa_I,1,mean)
    # calculate I1
    I1_C = (Y - Y_bar_model_C)*(f_m_noj_0c.list[[k]]*f_m1_marginal.list[[k]]/f_m1c.list[[k]])*(A/a_prob)
    ###### 2nd row of EIF
    I2_C = (xi_C-eta_C)*(A/a_prob)
    ###### 3rd row of EIF
    I3_C = (kappa_C-eta_C)*((1-A)/(1-a_prob))
    # 4th row of EIF
    I4_C = eta_C
    # EIF-based estimator
    I_eif_C = I_eif_C + mean(I1_C+I2_C+I3_C+I4_C)
    I_g_C = I_g_C + mean(I4_C)

    ##### individual-level estimand
    # calculate I1
    I1_I = (Y - Y_bar_model_I)*(f_m_noj_0c.list[[k]]*f_m1_marginal.list[[k]]/f_m1c.list[[k]])*(A/a_prob)
    ###### 2nd row of EIF
    I2_I = (xi_I-eta_I)*(A/a_prob)
    ###### 3rd row of EIF
    I3_I = (kappa_I-eta_I)*((1-A)/(1-a_prob))
    # 4th row of EIF
    I4_I = eta_I
    # EIF-based estimator
    I_eif_I = I_eif_I + N*mean(I1_I+I2_I+I3_I+I4_I)
    I_g_I = I_g_I + N*mean(I4_I)
  }
  list(Cluster=c(m1=I_g_C/K,eif=I_eif_C/K),
       Individual=c(m1=I_g_I/K,eif=I_eif_I/K)/mean(N.list))
}

#' Estimation of the mediation functional tau based on the machine learning estimator under original EIFs
#'
#' @param data the dataset
#' @param formula.M formula for the mediator model
#' @param formula.Y formula for the outcome model
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param V number of groups based on the cross-fitting procedure
#'
#' @returns The estimates of tau
MF_tau_nonpar_est1 = function(data,formula.M,formula.Y,is_stab=1,V=5,is_true_copula=1,myseed=2023) {
  #is_stab=1;V=5;is_true_copula=1;myseed=2023
  # create a dataset including all covariates (including the interaction terms)
  data1 = data0 = data; data1[,A_name] = 1; data0[,A_name] = 0
  data_m = model.matrix(formula.M,data)[,-1];data_m1 = model.matrix(formula.M,data1)[,-1];data_m0 = model.matrix(formula.M,data0)[,-1]
  additional_columns = setdiff(colnames(data_m),colnames(data))
  if (length(additional_columns)>0) {
    data = cbind(data,data_m[,additional_columns,drop=FALSE])
    data1 = cbind(data1,data_m1[,additional_columns,drop=FALSE])
    data0 = cbind(data0,data_m0[,additional_columns,drop=FALSE])
    colnames(data) = colnames(data1)=colnames(data0)=gsub("[(*) ]","",colnames(data))
  }
  # names of the mediator and outcome models
  M_model_names = gsub("[(*) ]","",colnames(data_m))
  Y_model_names = gsub("[(*) ]","",colnames(model.matrix(formula.Y,data))[-1])
  # treatment probability
  data_sum = unique(data[,c(id_name,N_name,A_name)])
  N.list = data_sum[,N_name]
  a_prob = 0.5
  # number of clusters
  K = length(unique(data[,id_name]))
  I_eif_C = I_eif_I =rep(0,K)
  I_g_C = I_g_I =rep(0,K)
  set.seed(myseed)
  folds = createFolds(1:K, k = V)
  v=1
  B_monte_carlo=25
  Sigma_M.list = mean_m0_marginal.list = mean_m1_marginal.list =
    sd_m0_marginal.list = sd_m1_marginal.list = f_m_noj_0c.list =
    f_m1_marginal.list = f_m1c.list = as.list(rep(1,K))
  for (v in (1:V)) {
    index_main = -which(data[,id_name] %in% folds[[v]])
    index_vali = which(data[,id_name] %in% folds[[v]])
    data_main=data[index_main,]; data_main1=data1[index_main,]; data_main0=data0[index_main,]
    data_vali=data[index_vali,]; data_vali1=data1[index_vali,]; data_vali0=data0[index_vali,]
    data_sum_main = data_sum[-which(data_sum[,id_name] %in% folds[[v]]),]
    data_sum_vali = data_sum[which(data_sum[,id_name] %in% folds[[v]]),]
    # mediator model (coef & sd of error term)
    m_m = SuperLearner(
      Y          = data_main[,M_name],
      X          = data_main[,M_model_names,drop=F],
      family     = gaussian(),
      id = data_main[,id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    m_m$m_sigma = sqrt(mean((data_main[,M_name] - predict(m_m, onlySL = TRUE)$pred[,1])^2))
    m_m$epsilon = qnorm(pnorm(data_main[,M_name],mean=predict(m_m, onlySL = TRUE)$pred[,1],sd=m_m$m_sigma))
    # copula estimate
    if (is_true_copula==1)  copula_res = est_copula(m_m,data_main,id_name="id",N_name="N")
    if (is_true_copula==0)  copula_res = 0.0001
    # outcome model
    m_y = SuperLearner(
      Y          = data_main[which(data_main[,A_name]==1),Y_name],
      X          = data_main[which(data_main[,A_name]==1),Y_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==1),id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    # prediction for the validation data
    Stab_C = Stab_I = rep(1,dim(data_vali)[1])
    for (k in folds[[v]]) {
      ind = which(data_vali[,id_name]==k)
      N=length(ind)
      A = data_vali[ind[1],A_name]
      # f(M|0,C)
      Sigma_M = (1-copula_res)*diag(N) + copula_res*(rep(1,N) %*% t(rep(1,N)))
      mean_m0_marginal = predict(m_m, newdata = data_vali0[ind,], onlySL = TRUE)$pred[,1]
      sd_m0_marginal = m_m$m_sigma
      epsilon_m0 = qnorm(pnorm(data_vali[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal))
      f_m0_marginal = dnorm(data_vali[ind,M_name],mean=mean_m0_marginal,sd = sd_m0_marginal)
      f_epsilon_m0 = dnorm(epsilon_m0)
      f_m_noj_0c=c()
      for (j in (1:N)) {
        copula_m_noj_0 = mvtnorm::dmvnorm(epsilon_m0[-j],mean=rep(0,N-1),sigma=Sigma_M[1:(N-1),1:(N-1)])
        f_m_noj_0c[j] = copula_m_noj_0*prod(f_m0_marginal[-j]/f_epsilon_m0[-j])
      }
      # f(M|1,C)
      mean_m1_marginal = predict(m_m, newdata = data_vali1[ind,], onlySL = TRUE)$pred[,1]
      sd_m1_marginal = m_m$m_sigma
      epsilon_m1 = qnorm(pnorm(data_vali[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal))
      f_m1_marginal = dnorm(data_vali[ind,M_name],mean=mean_m1_marginal,sd = sd_m1_marginal)
      f_epsilon_m1 = dnorm(epsilon_m1)
      copula_m1 = mvtnorm::dmvnorm(epsilon_m1,mean=rep(0,N),sigma=Sigma_M)
      f_m1c = copula_m1*prod(f_m1_marginal/f_epsilon_m1)
      # write the lists
      Sigma_M.list[[k]] = Sigma_M
      mean_m0_marginal.list[[k]] = mean_m0_marginal
      mean_m1_marginal.list[[k]] = mean_m1_marginal
      sd_m0_marginal.list[[k]] = sd_m0_marginal
      sd_m1_marginal.list[[k]] = sd_m1_marginal
      f_m_noj_0c.list[[k]] = f_m_noj_0c
      f_m1_marginal.list[[k]] = f_m1_marginal
      f_m1c.list[[k]] = f_m1c
      Stab_C[ind] = (1/a_prob)*(f_m1_marginal*f_m_noj_0c/f_m1c)*(1/N)
      Stab_I[ind] = (1/a_prob)*(f_m1_marginal*f_m_noj_0c/f_m1c)
    }
    if (is_stab==0) Stab_C = Stab_I = rep(1,dim(data)[1])
    # TMLE if (is_stab = 1)
    Y_1step = predict(m_y, newdata = data_vali, onlySL = TRUE)$pred[,1]
    if (is_stab == 1) {
      est_y_tmle_C = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==1),Y_name],
                             weights=Stab_C[which(data_vali$A==1)],
                             offset=Y_1step[which(data_vali$A==1)],
                             family=gaussian(link = "identity"))$coef
      est_y_tmle_I = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==1),Y_name],
                             weights=Stab_I[which(data_vali$A==1)],
                             offset=Y_1step[which(data_vali$A==1)],
                             family=gaussian(link = "identity"))$coef
    } else {
      est_y_tmle_C = est_y_tmle_I = 0
    }
    Y_final_C = Y_1step + est_y_tmle_C
    Y_final_I = Y_1step + est_y_tmle_I
    B_monte_carlo=25
    for (k in folds[[v]]) {
      ind = which(data_vali[,id_name]==k)
      N=length(ind)
      A = data_vali[ind[1],A_name]
      ###### construct the EIF
      # Y_bar
      Y = data_vali[ind,Y_name]
      M = data_vali[ind,M_name]
      # model-based outcome mean E[Y_bar|1,M,C]
      Y_bar_model_C = Y_final_C[ind]
      Y_bar_model_I = Y_final_I[ind]
      # calculate the intergals
      copula_M = pnorm(mvtnorm::rmvnorm(B_monte_carlo,mean=rep(0,N),sigma=Sigma_M.list[[k]]))
      M_sim_0  = qnorm(copula_M,mean= rep(1,B_monte_carlo) %*% t(mean_m0_marginal.list[[k]]),sd=sd_m0_marginal.list[[k]])
      copula_M = pnorm(mvtnorm::rmvnorm(B_monte_carlo,mean=rep(0,N),sigma=Sigma_M.list[[k]]))
      M_sim_1  = qnorm(copula_M,mean= rep(1,B_monte_carlo) %*% t(mean_m1_marginal.list[[k]]),sd=sd_m1_marginal.list[[k]])
      eta_C=kappa_C=xi_C = matrix(NA, nrow = N, ncol = B_monte_carlo)
      eta_I=kappa_I=xi_I = matrix(NA, nrow = N, ncol = B_monte_carlo)
      data_now1 = data_now2 = data_now3 = data_vali[ind,]
      for ( l in (1:B_monte_carlo)) {
        data_now1[,M_name] = M_sim_1[l,]; data_now1[,M_bar_name] = (sum(M_sim_0[l,])-M_sim_0[l,])/(N-1)
        eta_C[,l] = est_y_tmle_C+predict(m_y,newdata = data_now1, onlySL = TRUE)$pred[,1]
        eta_I[,l] = est_y_tmle_I+predict(m_y,newdata = data_now1, onlySL = TRUE)$pred[,1]
        data_now2[,M_bar_name] = (sum(M_sim_0[l,])-M_sim_0[l,])/(N-1)
        xi_C[,l] = est_y_tmle_C+predict(m_y,newdata = data_now2, onlySL = TRUE)$pred[,1]
        xi_I[,l] = est_y_tmle_I+predict(m_y,newdata = data_now2, onlySL = TRUE)$pred[,1]
        data_now3[,M_name] = M_sim_1[l,]
        kappa_C[,l] = est_y_tmle_C+predict(m_y,newdata = data_now3, onlySL = TRUE)$pred[,1]
        kappa_I[,l] =est_y_tmle_I+predict(m_y,newdata = data_now3, onlySL = TRUE)$pred[,1]
      }
      eta_C = apply(eta_C,1,mean)
      eta_I = apply(eta_I,1,mean)
      xi_C = apply(xi_C,1,mean)
      xi_I = apply(xi_I,1,mean)
      kappa_C = apply(kappa_C,1,mean)
      kappa_I = apply(kappa_I,1,mean)
      # calculate I1
      I1_C = (Y - Y_bar_model_C)*(f_m_noj_0c.list[[k]]*f_m1_marginal.list[[k]]/f_m1c.list[[k]])*(A/a_prob)
      ###### 2nd row of EIF
      I2_C = (xi_C-eta_C)*(A/a_prob)
      ###### 3rd row of EIF
      I3_C = (kappa_C-eta_C)*((1-A)/(1-a_prob))
      # 4th row of EIF
      I4_C = eta_C

      ##### individual-level estimand
      # calculate I1
      I1_I = (Y - Y_bar_model_I)*(f_m_noj_0c.list[[k]]*f_m1_marginal.list[[k]]/f_m1c.list[[k]])*(A/a_prob)
      ###### 2nd row of EIF
      I2_I = (xi_I-eta_I)*(A/a_prob)
      ###### 3rd row of EIF
      I3_I = (kappa_I-eta_I)*((1-A)/(1-a_prob))
      # 4th row of EIF
      I4_I = eta_I

      I_eif_C[k] = mean(I1_C+I2_C+I3_C+I4_C)
      I_g_C[k] = mean(I4_C)
      I_eif_I[k] = mean(I1_I+I2_I+I3_I+I4_I)*N
      I_g_I[k] = mean(I4_I)*N
    }
  }
  list(Cluster=c(m1=sum(I_g_C)/K,eif=sum(I_eif_C)/K),
       Individual=c(m1=sum(I_g_I)/K,eif=sum(I_eif_I)/K)/mean(N.list),
       I_eif_C = I_eif_C,
       I_eif_I = I_eif_I)
}


#' Estimation of the mediation functional theta(a,a*) based on the semiparametric doubly robust estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param a1 the value *a* in *\theta(a,a^\star)*
#' @param a0 the value *a^\star* in *\theta(a,a^\star)*
#' @param formula.A formula for the conditional treatment model s(a,m,c,n)
#' @param formula.Y formula for the outcome model eta_j(a,m,c,n)
#' @param formula.Ya formula for the outcome model eta_j^star(a,a*,c,n)
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns The estimates of theta(a,a*)
MF_theta_par_est2 = function(data,a1=1,a0=0,formula.A,formula.Y,formula.Ya,is_stab=1) {
  #a1=1;a0=0;is_stab=1
  # treatment probability
  data_sum = aggregate(data, by = list(data[,id_name]), mean)[,-1]
  a_prob = 0.5
  # number of clusters
  K = dim(data_sum)[1]
  # model matrix for outcome
  data_y = model.matrix(formula.Y,data)
  # conditional treatment probability (s model)
  m_s = glm(formula.A,data_sum,family=binomial(link="logit"))
  s1_pred    = predict(m_s, newdata = data_sum,type="response")
  s0_pred = 1-s1_pred
  w_sum = ((a0*s1_pred+(1-a0)*s0_pred)/(a1*s1_pred+(1-a1)*s0_pred))*((a_prob*a1+(1-a_prob)*(1-a1))/(a_prob*a0+(1-a_prob)*(1-a0)))
  w_ind = rep(w_sum,times=data_sum[,N_name])
  Stab_C = 1/(a_prob*a1+(1-a_prob)*(1-a1))*w_ind*(1/data[,N_name])
  Stab_I = 1/(a_prob*a1+(1-a_prob)*(1-a1))*w_ind
  if (is_stab==0) Stab_C = Stab_I = rep(1,dim(data)[1])
  if (a1==a0) Stab_C = Stab_I = rep(1,dim(data)[1])
  # outcome model
  y_model_C = y_model_I = list()

  X_y = model.matrix(formula.Y,data[which(data[,A_name]==a1),])
  Y_y = data[which(data[,A_name]==a1),Y_name]
  weights_y = sqrt(Stab_C[which(data[,A_name]==a1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  y_model_C$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)
  eta_a1_C = c(data_y %*% y_model_C$y_coef)

  X_y = model.matrix(formula.Y,data[which(data[,A_name]==a1),])
  Y_y = data[which(data[,A_name]==a1),Y_name]
  weights_y = sqrt(Stab_I[which(data[,A_name]==a1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  y_model_I$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)
  eta_a1_I = c(data_y %*% y_model_I$y_coef)

  # marginal outcome model
  data_y_mar = cbind(data,eta_a1_I,eta_a1_C)
  m_y_mar_C = lm.fit(x=model.matrix(formula.Ya,data=data_y_mar[which(data_y_mar[,A_name]==a0),]),
                     y=data_y_mar[which(data_y_mar[,A_name]==a0),"eta_a1_C"])
  mu_y_mar_C = c(model.matrix(formula.Ya,data=data_y_mar) %*% m_y_mar_C$coefficients)
  m_y_mar_I = lm.fit(x=model.matrix(formula.Ya,data=data_y_mar[which(data_y_mar[,A_name]==a0),]),
                     y=data_y_mar[which(data_y_mar[,A_name]==a0),"eta_a1_I"])
  mu_y_mar_I = c(model.matrix(formula.Ya,data=data_y_mar) %*% m_y_mar_I$coefficients)

  # write EIF
  k=1
  I_eif_C = I_eif_I = rep(0,K)
  I_g_C = I_g_I = rep(1,K)
  N.list=c()
  for (k in (1:K)) {
    if (k==1) {ind = 1:data[1,N_name]}
    if (k>1) {ind = (max(ind)+1):(max(ind)+data[max(ind)+1,N_name])}
    N=length(ind);N.list=c(N.list,N)
    A = data[ind[1],A_name]
    # cluser-average estimand
    Y_bar = data[ind,Y_name]
    Y_bar_model_C = eta_a1_C[ind]
    I1_C = (Y_bar - Y_bar_model_C)*w_ind[ind]*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
    I2_C = (Y_bar_model_C-mu_y_mar_C[ind])*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
    I3_C = mu_y_mar_C[ind]
    I_eif_C[k] = mean(I1_C+I2_C+I3_C)
    I_g_C[k] = mean(I3_C)
    # individual-average estimand
    Y_bar = data[ind,Y_name]
    Y_bar_model_I = eta_a1_I[ind]
    I1_I = (Y_bar - Y_bar_model_I)*w_ind[ind]*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
    I2_I = (Y_bar_model_I-mu_y_mar_I[ind])*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
    I3_I = mu_y_mar_I[ind]
    I_eif_I[k] = mean(I1_I+I2_I+I3_I)*N
    I_g_I[k] = mean(I3_I)*N
  }
  list(Cluster=c(m1=sum(I_g_C)/K,eif=sum(I_eif_C)/K),
       Individual=c(m1=sum(I_g_I)/K,eif=sum(I_eif_I)/K)/mean(N.list))
}

#' Estimation of the mediation functional theta(a,a*) based on the machine learner estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param a1 the value *a* in *\theta(a,a^\star)*
#' @param a0 the value *a^\star* in *\theta(a,a^\star)*
#' @param formula.A formula for the conditional treatment model s(a,m,c,n)
#' @param formula.Y formula for the outcome model eta_j(a,m,c,n)
#' @param formula.Ya formula for the outcome model eta_j^star(a,a*,c,n)
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns The estimates of theta(a,a*)
MF_theta_nonpar_est2 = function(data,a1=1,a0=0,formula.A,formula.Y,formula.Ya,is_stab=1,V=5,myseed=2023) {
  #a1=1;a0=0;is_stab=1;V=5;myseed=2023
  # create a dataset including all covariates (including the interaction terms)
  data1 = data0 = data; data1[,A_name] = 1; data0[,A_name] = 0
  data_sum = aggregate(data, by = list(data[,id_name]), mean)[,-1]
  data_a = model.matrix(formula.A,data_sum)[,-1]
  additional_columns = setdiff(colnames(data_a),colnames(data))
  if (length(additional_columns)>0) {
    data_sum = cbind(data_sum,data_a[,additional_columns,drop=FALSE])
    colnames(data_sum) =gsub("[(*) ]","",colnames(data_sum))
  }
  # names of the mediator and outcome models
  A_model_names = gsub("[(*) ]","",colnames(data_a))
  Y_model_names = colnames(model.matrix(formula.Y,data))[-1]
  Ya_model_names = colnames(model.matrix(formula.Ya,data))[-1]
  # treatment probability
  N.list = data_sum[,N_name]
  a_prob = 0.5
  # number of clusters
  K = dim(data_sum)[1]

  # initialize influence function values
  I_eif_C = I_eif_I =rep(0,K)
  I_g_C = I_g_I =rep(0,K)
  set.seed(myseed)
  folds = createFolds(1:K, k = V)
  v=1
  for (v in (1:V)) {
    index_main = -which(data[,id_name] %in% folds[[v]])
    index_vali = which(data[,id_name] %in% folds[[v]])
    data_main=data[index_main,]; data_main1=data1[index_main,]; data_main0=data0[index_main,]
    data_vali=data[index_vali,]; data_vali1=data1[index_vali,]; data_vali0=data0[index_vali,]
    data_sum_main = data_sum[ -which(data_sum[,id_name] %in% folds[[v]]),]
    data_sum_vali = data_sum[which(data_sum[,id_name] %in% folds[[v]]),]
    # propensity score models (A|M,C,N)
    m_s = SuperLearner(
      Y          = data_sum_main[,A_name],
      X          = data_sum_main[,A_model_names,drop=F],
      family     = binomial(),
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE, trimLogit = 0.01),
      cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
    )
    s1_pred    = predict(m_s, newdata = data_sum_vali[,A_model_names], onlySL = TRUE)$pred[,1]
    s0_pred = 1-s1_pred
    w_sum = ((a0*s1_pred+(1-a0)*s0_pred)/(a1*s1_pred+(1-a1)*s0_pred))*((a_prob*a1+(1-a_prob)*(1-a1))/(a_prob*a0+(1-a_prob)*(1-a0)))
    w_ind = rep(w_sum,times=data_sum_vali$N)
    Stab_C = 1/(a_prob*a1+(1-a_prob)*(1-a1))*w_ind*(1/data_vali[,N_name])
    Stab_I = 1/(a_prob*a1+(1-a_prob)*(1-a1))*w_ind
    if (is_stab==0) Stab_C = Stab_I = rep(1,dim(data_vali)[1])
    if (a1==a0) Stab_C = Stab_I = rep(1,dim(data_vali)[1])

    # outcome model
    m_y = SuperLearner(
      Y          = data_main[which(data_main[,A_name]==a1),Y_name],
      X          = data_main[which(data_main[,A_name]==a1),Y_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==a1),id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    # stablized outcome estimate
    Y_1step_vali = predict(m_y, newdata = data_vali, onlySL = TRUE)$pred[,1]
    if (is_stab == 1) {
      est_y_tmle_C = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==a1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==a1),Y_name],
                             weights=Stab_C[which(data_vali$A==a1)],
                             offset=Y_1step_vali[which(data_vali$A==a1)],
                             family=gaussian(link = "identity"))$coef
      est_y_tmle_I = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==a1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==a1),Y_name],
                             weights=Stab_I[which(data_vali$A==a1)],
                             offset=Y_1step_vali[which(data_vali$A==a1)],
                             family=gaussian(link = "identity"))$coef
    } else {
      est_y_tmle_C = est_y_tmle_I = 0
    }
    Y_final_C_vali = Y_1step_vali + est_y_tmle_C
    Y_final_I_vali = Y_1step_vali + est_y_tmle_I

    Y_1step_main = predict(m_y, newdata = data_main, onlySL = TRUE)$pred[,1]
    Y_final_C_main = Y_1step_main + est_y_tmle_C
    Y_final_I_main = Y_1step_main + est_y_tmle_I

    # marginal outcome model
    m_y_mar_C = SuperLearner(
      Y          = Y_final_C_main[which(data_main[,A_name]==a0)],
      X          = data_main[which(data_main[,A_name]==a0),Ya_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==a0),id_name],
      SL.library = c("SL.glm","SL.ranger"), #"SL.randomForest"
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    m_y_mar_I = SuperLearner(
      Y          = Y_final_I_main[which(data_main[,A_name]==a0)],
      X          = data_main[which(data_main[,A_name]==a0),Ya_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==a0),id_name],
      SL.library = c("SL.glm","SL.ranger"), #"SL.randomForest"
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    mu_C_vali = predict(m_y_mar_C,newdata=data_vali[,Ya_model_names,drop=F],onlySL = TRUE)$pred[,1]
    mu_I_vali = predict(m_y_mar_I,newdata=data_vali[,Ya_model_names,drop=F],onlySL = TRUE)$pred[,1]
    for (k in (folds[[v]])) {
      ind = which(data_vali[,id_name]==k)
      N=length(ind);
      A = data_vali[ind[1],A_name]
      # cluser-average estimand
      Y_bar = data_vali[ind,Y_name]
      Y_bar_model_C = Y_final_C_vali[ind]
      I1_C = (Y_bar - Y_bar_model_C)*w_ind[ind]*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
      I2_C = (Y_bar_model_C-mu_C_vali[ind])*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
      I3_C = mu_C_vali[ind]
      I_eif_C[k] = mean(I1_C+I2_C+I3_C)
      I_g_C[k] = mean(I3_C)
      # individual-average estimand
      Y_bar = data_vali[ind,Y_name]
      Y_bar_model_I = Y_final_I_vali[ind]
      I1_I = (Y_bar - Y_bar_model_I)*w_ind[ind]*((A==a1)/(a_prob*a1 + (1-a_prob)*(1-a1)))
      I2_I = (Y_bar_model_I-mu_I_vali[ind])*((A==a0)/(a_prob*a0 + (1-a_prob)*(1-a0)))
      I3_I = mu_I_vali[ind]
      I_eif_I[k] = mean(I1_I+I2_I+I3_I)*N
      I_g_I[k] = mean(I3_I)*N
    }
  }
  list(Cluster=c(m1=sum(I_g_C)/K,eif=sum(I_eif_C)/K),
       Individual=c(m1=sum(I_g_I)/K,eif=sum(I_eif_I)/K)/mean(N.list),
       I_eif_C = I_eif_C,
       I_eif_I = I_eif_I)
}



#' Estimation of the mediation functional tau based on the semiparametric doubly robust estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param formula.A formula for the treatment model s(a,m,c,n)
#' @param formula.M formula for the mediator model kappa_j(a,mj,c,n)
#' @param formula.Ma formula for the mediator model k_j^star(a,m,c,n)
#' @param formula.Y formula for the outcome model eta_j(a,m,c,n)
#' @param formula.Yb formula for the outcome model eta_j^dagger(a,a*,mj,c,n)
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns The estimates of tau
MF_tau_par_est2 = function(data,formula.A,formula.M,formula.Ma,formula.Y,formula.Yb,is_stab=1) {
  # treatment probability
  data_sum = aggregate(data, by = list(data[,id_name]), mean)[,-1]
  N.list = data_sum[,N_name]
  a_prob = 0.5
  # number of clusters
  K = length(unique(data[,id_name]))
  # model matrix for mediator
  data0 = data1 = data; data0[,A_name]=0;data1[,A_name]=1;
  data_m0 =model.matrix(formula.M,data0)
  data_m1 =model.matrix(formula.M,data1)
  # model matrix for outcome
  data_y = model.matrix(formula.Y,data1)
  # s_model
  m_s = glm(formula.A,data_sum,family=binomial(link="logit"))
  s1_pred    = predict(m_s, newdata = data_sum,type="response")
  s0_pred = 1-s1_pred
  s_sum = c(s0_pred/s1_pred)*c(a_prob/(1-a_prob))
  s_ind = rep(s_sum,times=data_sum$N)
  # marginal mediator model
  m_m_mar = lm(formula.M,data=data)
  M_hat_mar = predict(m_m_mar,newdata=data,type="response")
  m_m_mar_sd = sd(data[,M_name] - M_hat_mar)
  M_hat_mar1 = predict(m_m_mar,newdata=data1,type="response")
  M_hat_mar0 = predict(m_m_mar,newdata=data0,type="response")
  # conditional mediator model M_j|M_{-j},A,C,N
  m_m_con = lm(formula.Ma,data=data)
  M_hat_con = predict(m_m_con,newdata=data,type="response")
  m_m_con_sd = sd(data[,M_name] - M_hat_con)
  w_ind = c(dnorm(data[,M_name],mean=predict(m_m_mar,newdata=data1,type="response"),sd=m_m_mar_sd))/c(
    dnorm(data[,M_name],mean=predict(m_m_con,newdata=data0,type="response"),sd=m_m_con_sd))*(s_ind)
  # stablized weights
  Stab_C = (1/a_prob)*w_ind*(1/data[,N_name])
  Stab_I = (1/a_prob)*w_ind
  if (is_stab==0) {Stab_C=Stab_I = rep(1,dim(data)[1])}
  # outcome model E[Y|X,M,C,N]
  m_y_C = m_y_I = list()
  X_y = model.matrix(formula.Y,data[which(data[,A_name]==1),])
  Y_y = data[which(data[,A_name]==1),Y_name]
  weights_y = sqrt(Stab_C[which(data[,A_name]==1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  m_y_C$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)
  eta1_C = c(model.matrix(formula.Y,data1) %*% m_y_C$y_coef)

  X_y = model.matrix(formula.Y,data[which(data[,A_name]==1),])
  Y_y = data[which(data[,A_name]==1),Y_name]
  weights_y = sqrt(Stab_I[which(data[,A_name]==1)])
  X_weight = (weights_y %*% matrix(1,nrow=1,ncol=dim(X_y)[2])) * X_y
  Y_weight = Y_y * weights_y
  m_y_I$y_coef = c(solve(t(X_weight) %*% X_weight,tol=10^(-27)) %*% t(X_weight) %*% Y_weight)
  eta1_I = c(model.matrix(formula.Y,data1) %*% m_y_I$y_coef)
  # outcome model estimates
  B_monte_carlo=25
  M1_sim = t(rep(1,B_monte_carlo) %*% t(M_hat_mar1))+ matrix(rnorm(B_monte_carlo*dim(data)[1],mean=0,sd=m_m_mar_sd),ncol=B_monte_carlo)
  mu_c_mat_C = mu_c_mat_I = matrix(0,ncol=B_monte_carlo,nrow=dim(data)[1])
  dat_sim = data1
  for (j in (1:B_monte_carlo)) {
    #dat_sim[,M_bar_name] = (dat_sim[,M_bar_name]*dat_sim[,N_name]-dat_sim[,M_name]+M1_sim[,j])/dat_sim[,N_name]
    dat_sim[,M_name] = M1_sim[,j]
    mu_c_mat_C[,j] = c(model.matrix(formula.Y,data=dat_sim) %*% m_y_C$y_coef)
    mu_c_mat_I[,j] = c(model.matrix(formula.Y,data=dat_sim) %*% m_y_I$y_coef)
  }
  mu_c_C = apply(mu_c_mat_C,1,mean)
  mu_c_I = apply(mu_c_mat_I,1,mean)
  # outcome model (mu_b)
  w_ind_b = c(dnorm(data[,M_name],mean=predict(m_m_mar,newdata=data0,type="response"),sd=m_m_mar_sd))/c(
    dnorm(data[,M_name],mean=predict(m_m_con,newdata=data0,type="response"),sd=m_m_con_sd))
  m_y_b_C = lm.fit(y=c(w_ind_b*eta1_C)[which(data[,A_name]==0)],
                   x=model.matrix(formula.Yb,data)[which(data[,A_name]==0),])
  mu_b_C = c(model.matrix(formula.Yb,data0) %*% m_y_b_C$coef)
  m_y_b_I = lm.fit(y=c(w_ind_b*eta1_I)[which(data[,A_name]==0)],
                   x=model.matrix(formula.Yb,data)[which(data[,A_name]==0),])
  mu_b_I = c(model.matrix(formula.Yb,data0) %*% m_y_b_I$coef)

  # calculate mu_d
  mu_d_mat_C = mu_d_mat_I = matrix(0,ncol=B_monte_carlo,nrow=dim(data)[1])
  dat_sim_d = data1
  for (j in (1:B_monte_carlo)) {
    dat_sim_d[,M_name] = M1_sim[,j]
    mu_d_mat_C[,j] = c(model.matrix(formula.Yb,data=dat_sim_d) %*% m_y_b_C$coef)
    mu_d_mat_I[,j] = c(model.matrix(formula.Yb,data=dat_sim_d) %*% m_y_b_I$coef)
  }
  mu_d_C = apply(mu_d_mat_C,1,mean)
  mu_d_I = apply(mu_d_mat_I,1,mean)

  I_g0_C = mu_d_C
  I_g_C = aggregate(x=I_g0_C,by=list(data[,id_name]),mean)[,2]

  I_eif0_C = data[,A_name]/a_prob*w_ind*(data[,Y_name]-eta1_C) +
    data[,A_name]/a_prob*(mu_b_C-mu_d_C) +
    (1-data[,A_name])/(1-a_prob)*(mu_c_C-mu_d_C) +
    mu_d_C
  I_eif_C = aggregate(x=I_eif0_C,by=list(data[,id_name]),mean)[,2]

  I_g0_I = data[,N_name]*mu_d_I
  I_g_I = aggregate(x=I_g0_I,by=list(data[,id_name]),mean)[,2]

  I_eif0_I = data[,N_name]*(data[,A_name]/a_prob*w_ind*(data[,Y_name]-eta1_I) +
                              data[,A_name]/a_prob*(mu_b_I-mu_d_I) +
                              (1-data[,A_name])/(1-a_prob)*(mu_c_I-mu_d_I) +
                              mu_d_I)
  I_eif_I = aggregate(x=I_eif0_I,by=list(data[,id_name]),mean)[,2]

  list(Cluster=c(m1=sum(I_g_C)/K,eif=sum(I_eif_C)/K),
       Individual=c(m1=sum(I_g_I)/K,eif=sum(I_eif_I)/K)/mean(N.list))
}


#' Estimation of the mediation functional tau based on the machine learning estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param formula.A formula for the treatment model s(a,m,c,n)
#' @param formula.M formula for the mediator model kappa_j(a,mj,c,n)
#' @param formula.Ma formula for the mediator model k_j^star(a,m,c,n)
#' @param formula.Y formula for the outcome model eta_j(a,m,c,n)
#' @param formula.Yb formula for the outcome model eta_j^dagger(a,a*,mj,c,n)
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param V number of groups based on the cross-fitting procedure
#'
#' @returns The estimates of tau
MF_tau_nonpar_est2 = function(data,formula.A,formula.M,formula.Ma,formula.Y,formula.Yb,is_stab=1,V=5,myseed=2023) {
  #is_stab=1;V=5;myseed=2023
  B_monte_carlo=25
  data1 = data0 = data; data1[,A_name] = 1; data0[,A_name] = 0
  data_m = model.matrix(formula.M,data)[,-1];data_m1 = model.matrix(formula.M,data1)[,-1];data_m0 = model.matrix(formula.M,data0)[,-1]
  data_ma = model.matrix(formula.Ma,data)[,-1];data_ma1 = model.matrix(formula.Ma,data1)[,-1];data_ma0 = model.matrix(formula.Ma,data0)[,-1]
  additional_columns1 = setdiff(colnames(data_m),colnames(data))
  additional_columns2 = setdiff(colnames(data_ma),colnames(data))
  if (length(additional_columns1)>0) {
    data = cbind(data,data_m[,additional_columns1,drop=FALSE])
    data1 = cbind(data1,data_m1[,additional_columns1,drop=FALSE])
    data0 = cbind(data0,data_m0[,additional_columns1,drop=FALSE])
    colnames(data) = colnames(data1)=colnames(data0)=gsub("[(*) ]","",colnames(data))
  }
  if (length(additional_columns2)>0 & length(additional_columns1)<1 ) {
    data = cbind(data,data_ma[,additional_columns2,drop=FALSE])
    data1 = cbind(data1,data_ma1[,additional_columns2,drop=FALSE])
    data0 = cbind(data0,data_ma0[,additional_columns2,drop=FALSE])
    colnames(data) = colnames(data1)=colnames(data0)=gsub("[(*) ]","",colnames(data))
  }
  data_sum = aggregate(data, by = list(data[,id_name]), mean)[,-1]
  data_a = model.matrix(formula.A,data_sum)[,-1]
  additional_columns = setdiff(colnames(data_a),colnames(data_sum))
  if (length(additional_columns)>0) {
    data_sum = cbind(data_sum,data_a[,additional_columns,drop=FALSE])
    colnames(data_sum) =gsub("[(*) ]","",colnames(data_sum))
  }
  # names of the mediator and outcome models
  A_model_names = gsub("[(*) ]","",colnames(model.matrix(formula.A,data_sum))[-1])
  M_model_names = gsub("[(*) ]","",colnames(model.matrix(formula.M,data))[-1])
  Ma_model_names = gsub("[(*) ]","",colnames(model.matrix(formula.Ma,data))[-1])
  Y_model_names = gsub("[(*) ]","",colnames(model.matrix(formula.Y,data))[-1])
  Yb_model_names = gsub("[(*) ]","",colnames(model.matrix(formula.Yb,data))[-1])
  # treatment probability
  N.list = data_sum[,N_name]
  a_prob = 0.5
  # number of clusters
  K = length(unique(data[,id_name]))

  # the efficient influence function
  I_eif_C = I_eif_I =rep(0,K)
  I_g_C = I_g_I =rep(0,K)

  set.seed(myseed)
  folds = createFolds(1:K, k = V)
  v=1
  for (v in (1:V)) {
    index_main = -which(data[,id_name] %in% folds[[v]])
    index_vali = which(data[,id_name] %in% folds[[v]])
    data_main=data[index_main,]; data_main1=data1[index_main,]; data_main0=data0[index_main,]
    data_vali=data[index_vali,]; data_vali1=data1[index_vali,]; data_vali0=data0[index_vali,]
    data_sum_main = data_sum[-which(data_sum[,id_name] %in% folds[[v]]),]
    data_sum_vali = data_sum[which(data_sum[,id_name] %in% folds[[v]]),]

    # propensity score models (A|M,C,N)
    m_s = SuperLearner(
      Y          = data_sum_main[,A_name],
      X          = data_sum_main[,A_model_names,drop=F],
      family     = binomial(),
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE, trimLogit = 0.01),
      cvControl  = list(V = 5L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL)
    )

    s1_pred    = predict(m_s, newdata = data_sum_vali[,A_model_names,drop=F], onlySL = TRUE)$pred[,1]
    s0_pred = 1-s1_pred
    s_sum = c(s0_pred/s1_pred)*c(a_prob/(1-a_prob))
    s_ind = rep(s_sum,times=data_sum_vali$N)
    # marginal mediator model M_j|A,C,N
    m_m_mar = SuperLearner(
      Y          = data_main[,M_name],
      X          = data_main[,M_model_names,drop=F],
      family     = gaussian(),
      id = data_main[,id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    M_hat_mar_main = predict(m_m_mar,newdata=data_main,onlySL = TRUE)$pred[,1]
    m_m_mar_sd = sd(data_main[,M_name] - M_hat_mar_main)
    M_hat_mar_main1 = predict(m_m_mar,newdata=data_main1,onlySL = TRUE)$pred[,1]
    M_hat_mar_main0 = predict(m_m_mar,newdata=data_main0,onlySL = TRUE)$pred[,1]
    M_hat_mar_vali1 = predict(m_m_mar,newdata=data_vali1,onlySL = TRUE)$pred[,1]
    # conditional mediator model M_j|M_{-j},A,C,N
    m_m_con = SuperLearner(
      Y          = data_main[,M_name],
      X          = data_main[,Ma_model_names,drop=F],
      family     = gaussian(),
      id = data_main[,id_name],
      SL.library = c("SL.glm","SL.ranger"), #,"SL.glmnet","SL.gam","SL.step.interaction"
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )

    M_hat_con = predict(m_m_con,newdata=data_main,onlySL = TRUE)$pred[,1]
    m_m_con_sd = sd(data_main[,M_name] - M_hat_con)
    w_ind = c(dnorm(data_vali[,M_name],mean=predict(m_m_mar,newdata=data_vali1,onlySL = TRUE)$pred[,1],sd=m_m_mar_sd))/c(
      dnorm(data_vali[,M_name],mean=predict(m_m_con,newdata=data_vali0,onlySL = TRUE)$pred[,1],sd=m_m_con_sd))*(s_ind)
    # stablized weights
    Stab_C = (1/a_prob)*w_ind*(1/data_vali[,N_name])
    Stab_I = (1/a_prob)*w_ind
    if (is_stab==0) Stab_C = Stab_I = rep(1,dim(data_vali)[1])

    # outcome model
    m_y = SuperLearner(
      Y          = data_main[which(data_main[,A_name]==1),Y_name],
      X          = data_main[which(data_main[,A_name]==1),Y_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==1),id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    # TMLE if (is_stab = 1)
    Y_1step_vali = predict(m_y, newdata = data_vali, onlySL = TRUE)$pred[,1]
    Y_1step_main = predict(m_y, newdata = data_main, onlySL = TRUE)$pred[,1]
    if (is_stab == 1) {
      est_y_tmle_C = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==1),Y_name],
                             weights=Stab_C[which(data_vali$A==1)],
                             offset=Y_1step_vali[which(data_vali$A==1)],
                             family=gaussian(link = "identity"))$coef
      est_y_tmle_I = glm.fit(x=matrix(rep(1,length(data_vali[which(data_vali$A==1),Y_name])),ncol=1),
                             y=data_vali[which(data_vali$A==1),Y_name],
                             weights=Stab_I[which(data_vali$A==1)],
                             offset=Y_1step_vali[which(data_vali$A==1)],
                             family=gaussian(link = "identity"))$coef
    } else {
      est_y_tmle_C = est_y_tmle_I = 0
    }
    Y_final_C = Y_1step_vali + est_y_tmle_C
    Y_final_I = Y_1step_vali + est_y_tmle_I
    eta1_C = Y_1step_main + est_y_tmle_C
    eta1_I = Y_1step_main + est_y_tmle_I

    M1_main = t(rep(1,B_monte_carlo) %*% t(M_hat_mar_main1))+ matrix(rnorm(B_monte_carlo*dim(data_main)[1],mean=0,sd=m_m_mar_sd),ncol=B_monte_carlo)
    M1_vali = t(rep(1,B_monte_carlo) %*% t(M_hat_mar_vali1))+ matrix(rnorm(B_monte_carlo*dim(data_vali)[1],mean=0,sd=m_m_mar_sd),ncol=B_monte_carlo)
    mu_c_main_mat_C = mu_c_main_mat_I = matrix(0,ncol=B_monte_carlo,nrow=dim(data_main)[1])
    mu_c_vali_mat_C = mu_c_vali_mat_I = matrix(0,ncol=B_monte_carlo,nrow=dim(data_vali)[1])
    dat_sim_main = data_main1
    dat_sim_vali = data_vali1
    for (j in (1:B_monte_carlo)) {
      #dat_sim_main[,M_bar_name] = (dat_sim_main[,M_bar_name]*dat_sim_main[,N_name]-dat_sim_main[,M_name]+M1_main[,j])/dat_sim_main[,N_name]
      dat_sim_main[,M_name] = M1_main[,j]
      y_predict_now = predict(m_y,newdata=dat_sim_main,onlySL = TRUE)$pred[,1]
      mu_c_main_mat_C[,j] = est_y_tmle_C + y_predict_now
      mu_c_main_mat_I[,j] = est_y_tmle_I + y_predict_now

      #dat_sim_vali[,M_bar_name] = (dat_sim_vali[,M_bar_name]*dat_sim_vali[,N_name]-dat_sim_vali[,M_name]+M1_vali[,j])/dat_sim_vali[,N_name]
      dat_sim_vali[,M_name] = M1_vali[,j]
      y_predict_now = predict(m_y,newdata=dat_sim_vali,onlySL = TRUE)$pred[,1]
      mu_c_vali_mat_C[,j] = est_y_tmle_C + y_predict_now
      mu_c_vali_mat_I[,j] = est_y_tmle_I + y_predict_now
    }
    mu_c_main_C = apply(mu_c_main_mat_C,1,mean)
    mu_c_vali_C = apply(mu_c_vali_mat_C,1,mean)
    mu_c_main_I = apply(mu_c_main_mat_I,1,mean)
    mu_c_vali_I = apply(mu_c_vali_mat_I,1,mean)

    # outcome model (mu_b)
    w_ind_b_main = c(dnorm(data_main[,M_name],mean=predict(m_m_mar,newdata=data_main0,onlySL = TRUE)$pred[,1],sd=m_m_mar_sd))/c(
      dnorm(data_main[,M_name],mean=predict(m_m_con,newdata=data_main0,onlySL = TRUE)$pred[,1],sd=m_m_con_sd))
    m_y_b_C = SuperLearner(
      Y          = c(eta1_C*w_ind_b_main)[which(data_main[,A_name]==0)],
      X          = data_main[which(data_main[,A_name]==0),Yb_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==0),id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    m_y_b_I = SuperLearner(
      Y          = c(eta1_I*w_ind_b_main)[which(data_main[,A_name]==0)],
      X          = data_main[which(data_main[,A_name]==0),Yb_model_names,drop=F],
      family     = gaussian(),
      id = data_main[which(data_main[,A_name]==0),id_name],
      SL.library = c("SL.glm","SL.ranger"),
      control    = list(saveFitLibrary = TRUE),
      cvControl  = list(V = 5L, shuffle = TRUE, validRows = NULL)
    )
    mu_b_vali_C = predict(m_y_b_C,newdata=data_vali0,onlySL = TRUE)$pred[,1]
    mu_b_vali_I = predict(m_y_b_I,newdata=data_vali0,onlySL = TRUE)$pred[,1]

    # calculate mu_d
    mu_d_main_mat_C = mu_d_main_mat_I = matrix(0,ncol=B_monte_carlo,nrow=dim(data_main)[1])
    mu_d_vali_mat_C = mu_d_vali_mat_I = matrix(0,ncol=B_monte_carlo,nrow=dim(data_vali)[1])
    dat_sim_d_main = data_main1
    dat_sim_d_vali = data_vali1
    for (j in (1:B_monte_carlo)) {
      dat_sim_d_main[,M_name] = M1_main[,j]
      mu_d_main_mat_C[,j] = predict(m_y_b_C,newdata=dat_sim_d_main,onlySL = TRUE)$pred[,1]
      mu_d_main_mat_I[,j] = predict(m_y_b_I,newdata=dat_sim_d_main,onlySL = TRUE)$pred[,1]

      dat_sim_d_vali[,M_name] = M1_vali[,j]
      mu_d_vali_mat_C[,j] = predict(m_y_b_C,newdata=dat_sim_d_vali,onlySL = TRUE)$pred[,1]
      mu_d_vali_mat_I[,j] = predict(m_y_b_I,newdata=dat_sim_d_vali,onlySL = TRUE)$pred[,1]
    }
    mu_d_main_C = apply(mu_d_main_mat_C,1,mean)
    mu_d_vali_C = apply(mu_d_vali_mat_C,1,mean)
    mu_d_main_I = apply(mu_d_main_mat_I,1,mean)
    mu_d_vali_I = apply(mu_d_vali_mat_I,1,mean)
    for (k in folds[[v]]) {
      ind = which(data_vali[,id_name]==k)
      N=length(ind)
      A = data_vali[ind[1],A_name]
      ###### construct the EIF
      # Y_bar
      Y = data_vali[ind,Y_name]
      M = data_vali[ind,M_name]

      # cluster-average g-formula
      I_g_C[k] = mean(mu_d_vali_C[ind])
      I_eif_C[k] = mean(A/a_prob*w_ind[ind]*(Y-Y_final_C[ind]) +
                          A/a_prob*(mu_b_vali_C[ind]-mu_d_vali_C[ind]) +
                          (1-A)/(1-a_prob)*(mu_c_vali_C[ind]-mu_d_vali_C[ind]) +
                          mu_d_vali_C[ind])

      I_g_I[k] = N*mean(mu_d_vali_I[ind])
      I_eif_I[k] = N*mean(A/a_prob*w_ind[ind]*(Y-Y_final_I[ind]) +
                            A/a_prob*(mu_b_vali_I[ind]-mu_d_vali_I[ind]) +
                            (1-A)/(1-a_prob)*(mu_c_vali_I[ind]-mu_d_vali_I[ind]) +
                            mu_d_vali_I[ind])
    }
  }
  list(Cluster=c(m1=sum(I_g_C)/K,eif=sum(I_eif_C)/K),
       Individual=c(m1=sum(I_g_I)/K,eif=sum(I_eif_I)/K)/mean(N.list),
       I_eif_C = I_eif_C,
       I_eif_I = I_eif_I)
}


#' Point estimation of mediation effects based on the semiparametric doubly robust estimator under original EIFs
#'
#' @param data the dataset
#' @param formula.M formula for the mediator model
#' @param formula.Y formula for the outcome model
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns Point Estimates of NIE, NDE, and SME
mediate_par_est1 = function(data,formula.M,formula.Y,is_stab,is_true_copula=1) {
  m_res = est_m_marginal(formula.M,data)
  if (is_true_copula==1) {
    copula_res = est_copula(m_res,data,id_name="id",N_name="N")
  } else {
    copula_res = 0.0001
  }
  #print(copula_res)
  y1_res = est_y_marginal(formula.Y,data,A_level=1,A_name="A")
  y0_res = est_y_marginal(formula.Y,data,A_level=0,A_name="A")
  # estimate NIE (parametric)
  MF_11_par = MF_theta_par_est1(data,a1=1,a0=1,formula.M,formula.Y,m_res,y1_res,y0_res,copula_res,is_stab)
  MF_10_par = MF_theta_par_est1(data,a1=1,a0=0,formula.M,formula.Y,m_res,y1_res,y0_res,copula_res,is_stab)
  MF_00_par = MF_theta_par_est1(data,a1=0,a0=0,formula.M,formula.Y,m_res,y1_res,y0_res,copula_res,is_stab)
  NIE_C = MF_11_par$Cluster - MF_10_par$Cluster
  NIE_I = MF_11_par$Individual - MF_10_par$Individual
  NDE_C = MF_10_par$Cluster - MF_00_par$Cluster
  NDE_I = MF_10_par$Individual - MF_00_par$Individual
  # estimate SME (parametric)
  MF_110_par = MF_tau_par_est1(data,formula.M,formula.Y,m_res,y1_res,y0_res,copula_res,is_stab)
  SME_C = MF_11_par$Cluster - MF_110_par$Cluster
  SME_I = MF_11_par$Individual - MF_110_par$Individual
  out = c(NIE_C[1],NIE_I[1],
          NDE_C[1],NDE_I[1],
          SME_C[1],SME_I[1],
          NIE_C[2],NIE_I[2],
          NDE_C[2],NDE_I[2],
          SME_C[2],SME_I[2])
  names(out) = c("NIE_C_mf","NIE_I_mf","NDE_C_mf","NDE_I_mf","SME_C_mf","SME_I_mf",
                 "NIE_C_eif1_par","NIE_I_eif1_par","NDE_C_eif1_par","NDE_I_eif1_par","SME_C_eif1_par","SME_I_eif1_par")
  out
}



#' Mediation Analysis based on the semiparametric doubly robust estimator under original EIFs
#'
#' @param data the dataset
#' @param formula.M formula for the mediator model
#' @param formula.Y formula for the outcome model
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param B number of bootstrap samples
#'
#' @returns Point and interval estimates of NIE, NDE, and SME
mediate_par_est1_boot = function(data,formula.M,formula.Y,is_stab,B=100,is_true_copula=1) {
  K = max(data[,id_name])
  for (b in (1:B)) {
    if (b==1) {
      res0 = mediate_par_est1(data,formula.M,formula.Y,is_stab,is_true_copula)
      res_mat = matrix(NA,ncol=length(res0),nrow=B)
      res_mat[b,]=res0
      colnames(res_mat) = names(res0)
    } else {
      k_indices = sample(1:K,size=K,replace=TRUE)
      loc_list = sapply(k_indices, function(x) {which(data[,id_name] %in% x)} )
      i_indices = unlist(loc_list)
      id_indices = rep(1:K,times=sapply(loc_list, function(x) length(x) ))
      data_boot = data[i_indices,]
      data_boot[,id_name] = id_indices
      res_mat[b,] = tryCatch(mediate_par_est1(data_boot,formula.M,formula.Y,is_stab,is_true_copula),error= function(x) return(res0))
    }
  }
  CI_low = apply(res_mat,2,quantile,probs=0.025,na.rm=TRUE)
  CI_up = apply(res_mat,2,quantile,probs=0.975,na.rm=TRUE)
  cbind(res0,CI_low,CI_up)
}


#' Mediation Analysis based on the machine learning estimator under original EIFs
#'
#' @param data the dataset
#' @param formula.M formula for the mediator model
#' @param formula.Y formula for the outcome model
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns Point and interval estimates of NIE, NDE, and SME
mediate_nonpar_est1 = function(data,formula.M,formula.Y,is_stab,is_true_copula=1) {
  data_sum = aggregate(data, by = list(data[,id_name]), mean)[,-1]
  N.list = data_sum[,N_name]
  K=max(data_sum[,id_name])
  # estimate NIE (nonparametric)
  MF_11_nonpar = MF_theta_nonpar_est1(data,a1=1,a0=1,formula.M,formula.Y,is_stab,V=5,is_true_copula,myseed=2023)
  MF_10_nonpar = MF_theta_nonpar_est1(data,a1=1,a0=0,formula.M,formula.Y,is_stab,V=5,is_true_copula,myseed=2023)
  NIE_C = MF_11_nonpar$Cluster[2] - MF_10_nonpar$Cluster[2]
  NIE_I = MF_11_nonpar$Individual[2] - MF_10_nonpar$Individual[2]
  NIE_C_SE = sqrt(mean((MF_11_nonpar$I_eif_C-MF_10_nonpar$I_eif_C-NIE_C)^2)/K)
  NIE_I_SE = sqrt(mean(( (MF_11_nonpar$I_eif_I-MF_10_nonpar$I_eif_I-N.list*NIE_I)/mean(N.list))^2)/K)
  # estimate NDE (nonparametric)
  MF_00_nonpar = MF_theta_nonpar_est1(data,a1=0,a0=0,formula.M,formula.Y,is_stab,V=5,is_true_copula,myseed=2023)
  NDE_C = MF_10_nonpar$Cluster[2] - MF_00_nonpar$Cluster[2]
  NDE_I = MF_10_nonpar$Individual[2] - MF_00_nonpar$Individual[2]
  NDE_C_SE = sqrt(mean((MF_10_nonpar$I_eif_C-MF_00_nonpar$I_eif_C-NDE_C)^2)/K)
  NDE_I_SE = sqrt(mean(( (MF_10_nonpar$I_eif_I-MF_00_nonpar$I_eif_I-N.list*NDE_I)/mean(N.list))^2)/K)
  # estimate SME (nonparametric)
  MF_110_nonpar = MF_tau_nonpar_est1(data,formula.M,formula.Y,is_stab,V=5,is_true_copula,myseed=2023)
  SME_C = MF_11_nonpar$Cluster[2] - MF_110_nonpar$Cluster[2]
  SME_I = MF_11_nonpar$Individual[2] - MF_110_nonpar$Individual[2]
  SME_C_SE = sqrt(mean((MF_11_nonpar$I_eif_C-MF_110_nonpar$I_eif_C-SME_C)^2)/K)
  SME_I_SE = sqrt(mean(( (MF_11_nonpar$I_eif_I-MF_110_nonpar$I_eif_I-N.list*SME_I)/mean(N.list))^2)/K)

  point = c(NIE_C[1],NIE_I[1],NDE_C[1],NDE_I[1],SME_C[1],SME_I[1])
  names(point) = c("NIE_C_eif1_ml","NIE_I_eif1_ml","NDE_C_eif1_ml","NDE_eif1_ml","SME_C_eif1_ml","SME_I_eif1_ml")
  num_cov = length(X_names)+1
  CI_low = point - sqrt(K/(K-num_cov))*qt(0.975,df=K-num_cov)*c(NIE_C_SE,NIE_I_SE,NDE_C_SE,NDE_I_SE,SME_C_SE,SME_I_SE)
  CI_up = point + sqrt(K/(K-num_cov))*qt(0.975,df=K-num_cov)*c(NIE_C_SE,NIE_I_SE,NDE_C_SE,NDE_I_SE,SME_C_SE,SME_I_SE)
  cbind(point,CI_low,CI_up)
}


#' Point estimation of mediation effects based on the semiparametric doubly robust estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param formula.A formula for the treatment model s(a,m,c,n)
#' @param formula.M formula for the mediator model kappa_j(a,mj,c,n)
#' @param formula.Ma formula for the mediator model k_j^star(a,m,c,n)
#' @param formula.Y formula for the outcome model eta_j(a,m,c,n)
#' @param formula.Ya formula for the outcome model eta_j^star(a,a*,c,n)
#' @param formula.Yb formula for the outcome model eta_j^dagger(a,a*,mj,c,n)
#' @param is_stab stabilization or not (1 yes, 0 no)
#'
#' @returns Point Estimates of NIE, NDE, and SME
mediate_par_est2 = function(data,formula.A,formula.M,formula.Ma,
                                formula.Y,formula.Ya,formula.Yb,is_stab) {
  # estimate NIE (parametric)
  MF_11_par = MF_theta_par_est2(data,a1=1,a0=1,formula.A,formula.Y,formula.Ya,is_stab)
  MF_10_par = MF_theta_par_est2(data,a1=1,a0=0,formula.A,formula.Y,formula.Ya,is_stab)
  NIE_C = MF_11_par$Cluster - MF_10_par$Cluster
  NIE_I = MF_11_par$Individual - MF_10_par$Individual
  # estimate NDE (parametric)
  MF_00_par = MF_theta_par_est2(data,a1=0,a0=0,formula.A,formula.Y,formula.Ya,is_stab)
  NDE_C = MF_10_par$Cluster - MF_00_par$Cluster
  NDE_I = MF_10_par$Individual - MF_00_par$Individual
  # estimate SME (parametric)
  MF_110_par = MF_tau_par_est2(data,formula.A,formula.M,formula.Ma,formula.Y,formula.Yb,is_stab)
  SME_C = MF_11_par$Cluster - MF_110_par$Cluster
  SME_I = MF_11_par$Individual - MF_110_par$Individual
  out = c(NIE_C[1],NIE_I[1],NDE_C[1],NDE_I[1],SME_C[1],SME_I[1],
          NIE_C[2],NIE_I[2],NDE_C[2],NDE_I[2],SME_C[2],SME_I[2])
  names(out) = c("NIE_C_g","NIE_I_g","NDE_C_g","NDE_I_g","SME_C_g","SME_I_g",
                 "NIE_C_eif2_par","NIE_I_eif2_par","NDE_C_eif2_par","NDE_I_eif2_par","SME_C_eif2_par","SME_I_eif2_par")
  out
}



#' Mediation analysis based on the semiparametric doubly robust estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param formula.A formula for the treatment model s(a,m,c,n)
#' @param formula.M formula for the mediator model kappa_j(a,mj,c,n)
#' @param formula.Ma formula for the mediator model k_j^star(a,m,c,n)
#' @param formula.Y formula for the outcome model eta_j(a,m,c,n)
#' @param formula.Ya formula for the outcome model eta_j^star(a,a*,c,n)
#' @param formula.Yb formula for the outcome model eta_j^dagger(a,a*,mj,c,n)
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param B number of bootstrap samples
#'
#' @returns Point and interval estimates of NIE, NDE, and SME
mediate_par_est2_boot = function(data,formula.A,formula.M,formula.Ma,
                                     formula.Y,formula.Ya,formula.Yb,is_stab,B=100) {
  K = max(data[,id_name])
  for (b in (1:B)) {
    if (b==1) {
      res0 = mediate_par_est2(data,formula.A,formula.M,formula.Ma,formula.Y,formula.Ya,formula.Yb,is_stab)
      res_mat = matrix(NA,ncol=length(res0),nrow=B)
      res_mat[b,]=res0
      colnames(res_mat) = names(res0)
    } else {
      k_indices = sample(1:K,size=K,replace=TRUE)
      loc_list = sapply(k_indices, function(x) {which(data[,id_name] %in% x)} )
      i_indices = unlist(loc_list)
      id_indices = rep(1:K,times=sapply(loc_list, function(x) length(x) ))
      data_boot = data[i_indices,]
      data_boot[,id_name] = id_indices
      res_mat[b,] = tryCatch(mediate_par_est2(data_boot,formula.A,formula.M,formula.Ma,formula.Y,formula.Ya,formula.Yb,is_stab),error= function(x) return(res0))
    }
  }
  CI_low = apply(res_mat,2,quantile,probs=0.025,na.rm=TRUE)
  CI_up = apply(res_mat,2,quantile,probs=0.975,na.rm=TRUE)
  cbind(res0,CI_low,CI_up)
}




#' Mediation analysis based on the machine learning estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param formula.A formula for the treatment model s(a,m,c,n)
#' @param formula.M formula for the mediator model kappa_j(a,mj,c,n)
#' @param formula.Ma formula for the mediator model k_j^star(a,m,c,n)
#' @param formula.Y formula for the outcome model eta_j(a,m,c,n)
#' @param formula.Ya formula for the outcome model eta_j^star(a,a*,c,n)
#' @param formula.Yb formula for the outcome model eta_j^dagger(a,a*,mj,c,n)
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param B number of bootstrap samples
#'
#' @returns Point and interval estimates of NIE, NDE, and SME
mediate_nonpar_est2 = function(data,formula.A,formula.M,formula.Ma,
                                   formula.Y,formula.Ya,formula.Yb,is_stab) {
  data_sum = aggregate(data, by = list(data[,id_name]), mean)[,-1]
  N.list = data_sum[,N_name]
  K=max(data_sum[,id_name])
  # estimate NIE (nonparametric)
  MF_11_nonpar = MF_theta_nonpar_est2(data,a1=1,a0=1,formula.A,formula.Y,formula.Ya,is_stab,V=5,myseed=2023)
  MF_10_nonpar = MF_theta_nonpar_est2(data,a1=1,a0=0,formula.A,formula.Y,formula.Ya,is_stab,V=5,myseed=2023)
  NIE_C = MF_11_nonpar$Cluster[2] - MF_10_nonpar$Cluster[2]
  NIE_I = MF_11_nonpar$Individual[2] - MF_10_nonpar$Individual[2]
  NIE_C_SE = sqrt(mean((MF_11_nonpar$I_eif_C-MF_10_nonpar$I_eif_C-NIE_C)^2)/K)
  NIE_I_SE = sqrt(mean(( (MF_11_nonpar$I_eif_I-MF_10_nonpar$I_eif_I-N.list*NIE_I)/mean(N.list))^2)/K)
  # estimate NDE (nonparametric)
  MF_00_nonpar = MF_theta_nonpar_est2(data,a1=0,a0=0,formula.A,formula.Y,formula.Ya,is_stab,V=5,myseed=2023)
  NDE_C = MF_10_nonpar$Cluster[2] - MF_00_nonpar$Cluster[2]
  NDE_I = MF_10_nonpar$Individual[2] - MF_00_nonpar$Individual[2]
  NDE_C_SE = sqrt(mean((MF_10_nonpar$I_eif_C-MF_00_nonpar$I_eif_C-NDE_C)^2)/K)
  NDE_I_SE = sqrt(mean(( (MF_10_nonpar$I_eif_I-MF_00_nonpar$I_eif_I-N.list*NDE_I)/mean(N.list))^2)/K)
  # estimate SME (nonparametric)
  MF_110_nonpar = MF_tau_nonpar_est2(data,formula.A,formula.M,formula.Ma,formula.Y,formula.Yb,is_stab,V=5,myseed=2023)
  SME_C = MF_11_nonpar$Cluster[2] - MF_110_nonpar$Cluster[2]
  SME_I = MF_11_nonpar$Individual[2] - MF_110_nonpar$Individual[2]
  SME_C_SE = sqrt(mean((MF_11_nonpar$I_eif_C-MF_110_nonpar$I_eif_C-SME_C)^2)/K)
  SME_I_SE = sqrt(mean(( (MF_11_nonpar$I_eif_I-MF_110_nonpar$I_eif_I-N.list*SME_I)/mean(N.list))^2)/K)

  point = c(NIE_C[1],NIE_I[1],NDE_C[1],NDE_I[1],SME_C[1],SME_I[1])
  names(point) = c("NIE_C_np2","NIE_I_np2","NDE_C_np2","NDE_I_np2","SME_C_np2","SME_I_np2")
  num_cov = length(X_names)+1
  CI_low = point - sqrt(K/(K-num_cov))*qt(0.975,df=K-num_cov)*c(NIE_C_SE,NIE_I_SE,NDE_C_SE,NDE_I_SE,SME_C_SE,SME_I_SE)
  CI_up = point + sqrt(K/(K-num_cov))*qt(0.975,df=K-num_cov)*c(NIE_C_SE,NIE_I_SE,NDE_C_SE,NDE_I_SE,SME_C_SE,SME_I_SE)
  cbind(point,CI_low,CI_up)
}


#' Mediation analysis based on the machine learning estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param id_name cluster id
#' @param N_name cluster size
#' @param A_name treatment
#' @param M_name mediator
#' @param Y_name outcome
#' @param X_names individual-level covariates
#' @param V_names cluster-level covariates (if no cluster covariates, then set V_names=NULL)
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param B number of bootstrap
#'
#' @returns Point and interval estimates of NIE, NDE, and SME
mediateCRT_est1=function(data,id_name,N_name,A_name,M_name,Y_name,X_names,V_names,is_stab,B) {
  # id_name = "id"
  # N_name = "N"
  # A_name = "A"
  # M_name = "M"
  # Y_name = "Y"
  # X_names = c("X","X2")
  # V_names = c("C")
  # is_stab=1
  # B=10

  id_name <<- id_name
  N_name <<- N_name
  A_name <<- A_name
  M_name <<- M_name
  Y_name <<- Y_name
  X_names <<- X_names
  V_names <<- V_names

  # names for M_bar and X_bar
  M_bar_name <<- paste0(M_name,"_bar")
  X_bar_names <<- paste0(X_names,"_bar")

  # create group-average Xs and M
  x_cols <- intersect(c(X_names,M_name), names(data))
  data = data %>% group_by(id) %>% mutate(
    across(all_of(x_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}_allbar")) %>%
    ungroup()
  data = data %>% group_by(id) %>%
    mutate(across(all_of(x_cols),~ (sum(.x) - .x) / (n() - 1),.names = "{.col}_bar")) %>%
    ungroup()
  data = as.data.frame(data)
  # formulas for the regression models
  if (is.null(V_names)) {
    formula.Y = as.formula(paste(Y_name, "~", M_name, "+", paste0(M_name,"_bar"),
                      "+", paste(X_names,collapse = " + "),"+",N_name))
    formula.M = as.formula(paste(M_name,"~",A_name,"+","I(",A_name,"*",N_name,")","+",paste(X_names,collapse = " + "),
          "+",N_name))
  } else {
    formula.Y = as.formula(paste(Y_name, "~", M_name, "+", paste0(M_name,"_bar"),
                                 "+", paste(X_names,collapse = " + "),"+",paste(V_names,collapse = " + "),"+",N_name))
    formula.M = as.formula(paste(M_name,"~",A_name,"+","I(",A_name,"*",N_name,")","+",paste(X_names,collapse = " + "),
                                 "+",paste(V_names,collapse = " + "),"+",N_name))
  }
  res_par = mediate_par_est1_boot(data,formula.M,formula.Y,is_stab,B,is_true_copula=1)
  res_nonpar = mediate_nonpar_est1(data,formula.M,formula.Y,is_stab,is_true_copula=1)

  out = list(`Mediation functional Estimator` = res_par[1:6,],
             `Doubly Robust Estimator` = res_par[7:12,],
             `Machine Learning Estimator` = res_nonpar)
  colnames(out[[1]]) = c("point","CI_lower","CI_upper")
  colnames(out[[2]]) = c("point","CI_lower","CI_upper")
  colnames(out[[3]]) = c("point","CI_lower","CI_upper")
  out
}


#' Mediation analysis based on the machine learning estimator under reparameterized EIFs
#'
#' @param data the dataset
#' @param id_name cluster id
#' @param N_name cluster size
#' @param A_name treatment
#' @param M_name mediator
#' @param Y_name outcome
#' @param X_names individual-level covariates
#' @param V_names cluster-level covariates (if no cluster covariates, then set V_names=NULL)
#' @param is_stab stabilization or not (1 yes, 0 no)
#' @param B number of bootstrap
#'
#' @returns Point and interval estimates of NIE, NDE, and SME
mediateCRT_est2=function(data,id_name,N_name,A_name,M_name,Y_name,X_names,V_names,is_stab,B) {
  # id_name = "id"
  # N_name = "N"
  # A_name = "A"
  # M_name = "M"
  # Y_name = "Y"
  # X_names = c("X","X2")
  # V_names = c("C")
  # is_stab=1
  # B=10
  id_name <<- id_name
  N_name <<- N_name
  A_name <<- A_name
  M_name <<- M_name
  Y_name <<- Y_name
  X_names <<- X_names
  V_names <<- V_names

  # names for M_bar and X_bar
  M_bar_name <<- paste0(M_name,"_bar")
  X_bar_names <<- paste0(X_names,"_bar")

  # create group-average Xs and M
  x_cols <- intersect(c(X_names,M_name), names(data))
  data = data %>% group_by(id) %>% mutate(
    across(all_of(x_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}_allbar")) %>%
    ungroup()
  data = data %>% group_by(id) %>%
    mutate(across(all_of(x_cols),~ (sum(.x) - .x) / (n() - 1),.names = "{.col}_bar")) %>%
    ungroup()
  data = as.data.frame(data)
  # formulas for the regression models
  formula.Y = as.formula(paste(Y_name, "~", M_name, "+", paste0(M_name,"_bar"),
                               "+", paste(X_names,collapse = " + "),"+",paste(V_names,collapse = " + "),"+",N_name))
  formula.M = as.formula(paste(M_name,"~",A_name,"+","I(",A_name,"*",N_name,")","+",paste(X_names,collapse = " + "),
                               "+",paste(V_names,collapse = " + "),"+",N_name))
  formula.A = as.formula(paste(A_name,"~",paste0(M_name,"_allbar"),"+",paste(paste0(X_names,"_allbar"),collapse="+"),
                               "+",paste(V_names,collapse="+"),"+",N_name))
  formula.Ma = as.formula(paste(M_name,"~",A_name,"+","I(",A_name,"*",N_name,")","+",paste0(M_name,"_bar"),"+",paste(X_names,collapse = " + "),
                                "+",paste(V_names,collapse = " + "),"+",N_name))
  formula.Ya = as.formula(paste(Y_name, "~", paste(X_names,collapse = " + "),"+",paste(paste0(X_names,"_bar"),collapse = "+"),"+",
                                paste(V_names,collapse = " + "),"+",N_name))
  formula.Yb = as.formula(paste(Y_name, "~", M_name, "+",paste(X_names,collapse = " + "),"+",
                                paste(paste0(X_names,"_bar"),collapse = "+"),"+",paste(V_names,collapse = " + "),"+",N_name))

  res_par = mediate_par_est2_boot(data,formula.A,formula.M,formula.Ma,
                                      formula.Y,formula.Ya,formula.Yb,is_stab,B)
  res_nonpar = mediate_nonpar_est2(data,formula.A,formula.M,formula.Ma,
                                       formula.Y,formula.Ya,formula.Yb,is_stab)

  out = list(`Doubly Robust Estimator` = res_par[7:12,],
             `Machine Learning Estimator` = res_nonpar)
  colnames(out[[1]]) = c("point","CI_lower","CI_upper")
  colnames(out[[2]]) = c("point","CI_lower","CI_upper")
  out
}







