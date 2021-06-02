#' @export
#'
#'

trs_alasso = function(Z,Y,delta_i,r,lamb_vec,solu = TRUE)
{
  fun_E = function(Y, Z, delta, beta_ini, lamb_ini, r, Q = 60)
  {
    
    gq = gauss.quad(Q, kind="laguerre")
    nodes = gq$nodes
    whts  = gq$weights
    dgamma_Q = dgamma(nodes, shape = 1/r, scale = r)
    
    n = length(Y)  #n个Y值
    Y_mat_nn = matrix(rep(Y,n), nrow = n, byrow = TRUE)
    Lamb_ini = colSums( lamb_ini * (Y <= Y_mat_nn) )
    
    beta_Z_ini = c(Z %*% beta_ini)
    exp_beZ = exp(beta_Z_ini)
    lam_expbeZ_m = matrix( rep(lamb_ini* exp_beZ, Q ), nrow = Q, byrow = TRUE )
    Lam_expbeZ_m = matrix( rep(Lamb_ini* exp_beZ, Q ), nrow = Q, byrow = TRUE )
    E_de_m = whts * t( t(nodes * lam_expbeZ_m)^delta ) *
      exp( -nodes * Lam_expbeZ_m + matrix(rep(log(dgamma_Q)+nodes,n), nrow = Q)   )
    E_nu_m = nodes * E_de_m  #没有whts，因为包含在E_de_m中
    
    E = colSums(E_nu_m) / colSums(E_de_m)
    return(E)
  }
  G = function(x,r)
  {
    if(r>0)
    {return(log(1+r*x)/r)}
    else if(r==0)
    {return(x)}
    else{print("r<0")}
  }
  dG = function(x,r)
  {
    if(r>0)
    {return(1/(1+r*x))}
    else if(r==0)
    {return(1)}
    else{print("r<0")}
  }
  ddG = function(x,r)
  {
    if(r>0)
    {return(-r/(1+r*x)^2)}
    else if(r==0)
    {return(0)}
    else{print("r<0")}
  }
  fun_c = function(Y, Z, delta, beta_ini, lamb_ini, r)
  {
    n = length(Y)
    Y_j.m = matrix(rep(Y,n), nrow = n, byrow = TRUE)
    Y_i.m = matrix(rep(Y,n), nrow = n, byrow = FALSE)
    Y_jleqi = Y_j.m <= Y_i.m
    #Y_jgeqi = Y_j.m >= Y_i.m
    beta_Z_ini = c(Z %*% beta_ini)
    temp1 = Y_jleqi %*% (exp(beta_Z_ini) * lamb_ini)
    til.c = c( (delta == 0) * dG(temp1,r) +
                 (delta == 1) * (-ddG(temp1,r)/dG(temp1,r) + dG(temp1,r)) )
    return(til.c) 
  }
  loglik = function(n,delta,z,beta,Y,E)
  {
    Y_m = matrix(rep(Y,n), nrow = n, byrow = TRUE)
    Ie = (t(Y_m)<=Y_m)*matrix(rep( E*exp(Z%*%beta),n  ),nrow = n, byrow = TRUE)
    return( sum( delta*(Z%*%beta - log(rowSums(Ie))) ) )
  }
  normalize = function(x)
  {
    y = (x-mean(x))/sd(x)
    return(y)
  }
  n = length(Y)
  beta_len = ncol(Z)
  #------定义记录跳过的循环------
  skip_iter = rep(0,beta_len)#每一行是每次迭代的所有lamb值
  skip_para = matrix(0,nrow = 1, ncol = 2)
  colnames(skip_para) = c("lamb","adalasso内循环次数")
  #------初始值-----
  possibleError_ini <- tryCatch(
    {
      res_EM = EM.est(Y, Z, delta_i, alpha = r)
      beta_ini = c(res_EM$beta_new)
      lamb_ini = c(res_EM$lamb_Y)
    },error=function(e)
    {
      print(e)
    })
  if(inherits(possibleError_ini, "error"))
  {
    cat("\n","wrong EM","\n")
    return(list(beta_res = NA, 
                GCV_res = NA,
                lamb_res = NA,
                beta_all = NA,
                CSV_all = NA,
                lamb_all = NA,
                skip_iter = NA,
                skip_para = NA))
  }
  #------定义每个lamb的beta,GCV,loglik------
  beta_all = matrix(0, nrow = length(lamb_vec), ncol = beta_len)
  GCV = c()
  loglik_vec = c()
  #------开始循环lamb------
  for (lamb_ind in 1:length(lamb_vec))
  {
    lamb = lamb_vec[lamb_ind]
    beta_latmp = beta_ini
    skip_loop = FALSE #每个lamb开始默认为FALSE
    #------开始adalasso
    for (ada_rep in 1:200)
    {
      possibleError <- tryCatch(
        {
          til_c = fun_E(Y = Y, Z = Z, beta_ini = beta_latmp, delta = delta_i, lamb_ini = lamb_ini ,r = r )
          #til_c = fun_c(Y = Y, Z = Z, beta_ini = beta_latmp, delta = delta_i, lamb_ini = lamb_ini ,r = r )
          Exp_betaZ_tmp = c(exp(Z%*%beta_latmp))*til_c
          mat_Y = matrix(rep(Y,n), nrow = n, byrow =FALSE)#每一行相同
          Del_m = (Exp_betaZ_tmp*( mat_Y >= t(mat_Y)))
          Del = (t(Z) - ((t(Z) %*% (Del_m))/matrix(rep(colSums(Del_m),beta_len), nrow = beta_len, byrow = TRUE) )) %*% (-delta_i)
          Del2 = ddloglik_transmdl(n,delta_i[order(Y)],Z[order(Y),],beta_latmp,til_c[order(Y)])
          X = chol(Del2)
          W = solve(t(X))%*%(Del2%*%beta_latmp-Del)
          ada_new = shoot(X = X, y = W, beta_ini_s = beta_ini, lamb_s = lamb)
        },error=function(e)
        {
          print(e)
        })
      if(inherits(possibleError, "error"))
      {
        cat("\n","lamb:",lamb,"\n")
        skip_loop = TRUE
        print("skip赋值TRUE")
        print(skip_loop)
        #记录循环终止时的参数
        skip_para = rbind(skip_para,c(lamb,ada_rep))
        skip_iter[lamb_ind] = 1
        #终止循环，给ada_new和beta_latmp赋值成相等的大向量
        ada_new=beta_latmp=rep(10^(5),beta_len)
      }
      
      if(max(abs(ada_new-beta_latmp)) < 10^(-5)){break}
      beta_latmp = ada_new
      #print(beta_latmp)
    }
    #------GCV------
    #如果跳过循环那么GCV被赋值成100000，确保不会被选中
    if(skip_loop)
    {
      beta_all[lamb_ind,] = beta_latmp
      GCV[lamb_ind] = 10^(5)
    }else
    {
      beta_all[lamb_ind,] = beta_latmp
      A = diag(1/abs(beta_latmp*beta_ini))
      diag(A)[which(beta_latmp == 0)] = 0
      diag(A)[which(beta_ini == 0)] = 0
      p =sum(diag(  solve(Del2 + lamb*A) %*% Del2  ) )
      GCV[lamb_ind] = (-loglik(n,delta_i,Z,beta_latmp,Y,til_c)) / (n*(1-p/n)^2)
      cat("\n","lamb_ind:",lamb_ind)
      #print(beta_all[lamb_ind,])
      #loglik_vec[lamb_ind] = -loglik(n,delta_i,Z,beta_latmp,Y,til_c)
    }
    
  }
  
  
  #------solution path------
  #save.image("trans gcv n=100 r=0.5 rep=200.RData")
  if(solu)
  {
    solutionpath = t(beta_all)
    col_c = c(2:beta_len)
    plot(x = lamb_vec, y = solutionpath[1,], type = "l", col = col_c[1],
         xlim = c(0,max(lamb_vec)+20), ylim = c(min(solutionpath), max(solutionpath)+0.2))
    for (i in 2:beta_len) 
    {
      lines(x = lamb_vec, y = solutionpath[i,], col = col_c[i])
    }
    
    legend("right",                                    #图例位置为右上角
           legend=paste("beta",c(1:beta_len), sep = ""),        #图例内容
           col=col_c,                 #图例颜色
           lty=1,lwd=2) 
    
    abline(v=lamb_vec[which(GCV == min(GCV))],lwd=1,col="black")
  }
  return(list(beta_res = beta_all[which(GCV == min(GCV))[1],], 
              GCV_res = GCV[which(GCV == min(GCV))[1]],
              lamb_res = lamb_vec[which(GCV == min(GCV))[1]],
              beta_all = beta_all,
              GCV_all = GCV,
              skip_para = skip_para))
}
