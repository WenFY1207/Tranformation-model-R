#' @export


Simu = function(M, n, alpha.true, rho, cu.res = TRUE, cu.rep = TRUE, plot.Lamb = TRUE,
               alpha.know = TRUE,  alpha.scale = seq(0.5,1.0,by=0.1)){

  library(survival)
  library(survminer)
  library(statmod)

  num_of_repeat = M

  #true value: beta1, beta2, Lambda1, Lambda2
  para_true = c(-0.5, 1, rho*log(1+1), rho*log(1+2))

  result = matrix(nrow = 1, ncol = 5) #??¼beta1.2,lambda1.2,censor
  result_SE = matrix(nrow = 1, ncol = 4)
  result_CP = matrix(nrow = 1, ncol = 4)

  result_X1 = matrix(nrow = 1, ncol = n)
  result_X2 = matrix(nrow = 1, ncol = n)
  result_Y = matrix(nrow = 1, ncol = n)
  result_Ydot = matrix(nrow = 1, ncol = n)
  result_m = c()
  result_delta = matrix(nrow = 1, ncol = n)
  #ǰm??Ϊ?��?ֵ????n-m???? -1 ????
  result_lambda_Ydot = matrix(nrow = 1, ncol = n)

  result_till_now = matrix(nrow = 4, ncol = 4)#ʵʱ????
  colnames(result_till_now) = c("Bias","SD","SE","CP")
  rownames(result_till_now) = c("beta1","beta2","Lambda1","Lambda2")
  bias_till_now = rep(0,4)
  bias2_till_now = rep(0,4)
  SD_till_now = rep(0,4)
  SE_till_now = rep(0,4)
  CP_till_now = rep(0,4)

  if(alpha.know)
  {
    alpha.c = alpha.true
  }else
  {
    alpha.c = alpha.scale
  }



  length_alpha = length(alpha.c)
  result_alpha = c()
  result_lik = matrix(0,nrow = 1, ncol = length_alpha)

  H.fun = function(x,alpha)
  {
    if(alpha == 0)
    {
      x
    }else if(alpha > 0)
    {
      log(1+alpha * x) / alpha
    }else
    {
      print("H:wrong alpha < 0")
    }
  }

  Hdot.fun = function(x,alpha)
  {
    if(alpha == 0)
    {
      1
    }else if(alpha > 0)
    {
      1 / (1+alpha*x)
    }else
    {
      print("Hdot:wrong alpha < 0")
    }
  }

  Q = 60
  gq = gauss.quad(Q, kind="laguerre")
  nodes = gq$nodes
  whts  = gq$weights
  #??E????ESS??EB?õ?

  plot_x = seq(0.1,2.5,by=0.1)
  n_x = length(plot_x)
  plot_Lam_Ydot = matrix(0, nrow = 1, ncol= n_x)
  plot_upper = matrix(0, nrow = 1, ncol= n_x)
  plot_lower = matrix(0, nrow = 1, ncol= n_x)
  u_975 = qnorm(0.975)

  for(now_repeat in 1:num_of_repeat)
  {
    gen_data = generate_data(n, alpha.true, rho, c(-0.5,1), now_repeat)

    X = gen_data$X
    X1 = X[,1]
    X2 = X[,2]
    delta_i = gen_data$delta
    beta_X = gen_data$beta_X
    Y = gen_data$Y

    #??¼X1,X2,Y,delta_i,m??ֵ???ڵ?һ??
    result_X1 = rbind(X1, result_X1)
    result_X2 = rbind(X2, result_X2)
    result_Y = rbind(Y, result_Y)
    result_delta = rbind(delta_i, result_delta)

    ##?Ҿ???1,2?????ı?1С??Y_i?ı???i##
    ##????1????I_1??2??I_2
    I_1 = which(Y ==  1- min(1-Y[Y<=1]))
    I_2 = which(Y ==  2- min(2-Y[Y<=2]))



    ##############??ʼ????####################
    #-------????Y?Ķ???--------
    Y_hat = sort(unique(Y[delta_i==1]))
    l = length(Y_hat)

    Y_dot = sort(unique(Y))
    n = length(Y)  #n??Yֵ
    m = length(Y_dot)   # m?????ظ???Yֵ??Y'
    result_m = c(m,result_m)
    result_Ydot = rbind(c(Y_dot, rep(0,(n-m))), result_Ydot)

    #Y1 Y2 ... Yn
    #Y1 Y2 ... Yn  l*n
    Y_mat = matrix(rep(Y,l), nrow = l, byrow = TRUE)

    #Y1 Y2 ... Yn
    #Y1 Y2 ... Yn  m*n
    Y_mat_mn = matrix(rep(Y,m), nrow = m, byrow = TRUE)

    #Yhat 1 ... Yhat 1
    #Yhat 2 ... Yhat 2
    #...
    #Yhat l ... Yhat l   l*n
    Y_hat_mat = matrix(rep(Y_hat,n), nrow = l, byrow = FALSE)

    #Ydot 1 ... Ydot 1
    #Ydot 2 ... Ydot 2
    #...
    #Ydot m ... Ydot m   m*n
    Y_dot_mat = matrix(rep(Y_dot,n), nrow = m, byrow = FALSE)

    #SE???õ?
    Y_eq_Yhat = (Y_mat == Y_hat_mat)
    Y_geq_Yhat = (Y_mat >= Y_hat_mat)

    #?? beta^t+1 ?õ?
    m_Y_k = matrix(Y, ncol = n, nrow = n, byrow = TRUE)
    m_Y_i = matrix(Y, ncol = n, nrow = n, byrow = FALSE)
    m_Y_kgeqi = (m_Y_k>=m_Y_i)
    m_X1 = matrix(X1, nrow = n, ncol = n, byrow = TRUE)
    m_X2 = matrix(X2, nrow = n, ncol = n, byrow = TRUE)

    if(alpha.true > 0)
    {
      lik.fun = c()#??¼ÿ??alpha????Ȼ????ֵ
      beta_best = c()
      LambdaY_best = c()
      lambdaY_best = c()
      lambdaYd_best = c()
      for(al_i in 1:length(alpha.c))
      {
        alpha.i = alpha.c[al_i]
        # for(al_i in 1)
        # {
        #   alpha.i = 0.06
        EM.res = EM.est(Y = Y, alpha = alpha.i, Q = Q, X, delta = delta_i )
        beta_new = EM.res$beta_new
        Lambda_Y_new = EM.res$Lamb_Y
        lambda_Y_new = EM.res$lamb_Y
        lambda_Ydot_new = EM.res$lamb_Ydot

        beX = c(X %*% c(beta_new))
        Lam.ExpbeX = Lambda_Y_new * exp(c(X %*% c(beta_new)))
        lam.ExpbeX = lambda_Y_new * exp(c(X %*% c(beta_new)))
        #log(prod((lam.ExpbeX*(Hdot.fun(Lam.ExpbeX,alpha.i)))^delta_i * exp(-(H.fun(Lam.ExpbeX,alpha.i))))  )
        aa = delta_i*(log(lambda_Y_new) + beX + log(Hdot.fun(Lam.ExpbeX,alpha.i)))
        aa[is.na(aa)] = 0
        lik.new =  sum(aa - H.fun(Lam.ExpbeX,alpha.i))
        lik.fun = c( lik.fun, lik.new)


        if(al_i == 1)
        {
          beta_best = beta_new
          LambdaY_best = Lambda_Y_new
          lambdaY_best = lambda_Y_new
          lambdaYd_best = lambda_Ydot_new
          lik_best = lik.new
          alpha_best = alpha.i
          #print(beta_best)
        }else if((al_i>1)&(lik.fun[al_i-1]<lik.new))
        {
          beta_best = beta_new
          LambdaY_best = Lambda_Y_new
          lambdaY_best = lambda_Y_new
          lambdaYd_best = lambda_Ydot_new
          lik_best = lik.new
          alpha_best = alpha.i
          #print(beta_best)
          #print(alpha.i)
        }

      }
      result_lik = rbind(lik.fun, result_lik)
      result_alpha = c(alpha_best,result_alpha)
      result_lambda_Ydot = rbind(c(lambdaYd_best, rep(-1,(n-m))), result_lambda_Ydot)
      result = rbind(c(beta_best,LambdaY_best[I_1],LambdaY_best[I_2],1-sum(delta_i)/n),result)
      #-------------------------SE-------------------------
      #EM.res??EM.est()????һ?ν???,Yһ??ʱY_eq_Yhat??Y_geq_Yhat??ͬ
      Y_eq_Yhat = EM.res$Y_eq_Yhat
      Y_geq_Yhat =  EM.res$Y_geq_Yhat

      #??��Ĭ?϶???Ϊ????��
      #??E(S)E(B)E(SS^T)
      dgamma_Q = dgamma(nodes, shape = 1/alpha_best, scale = alpha_best)

      beta_X = beta_best[1]*X1 + beta_best[2]*X2
      #-------------------------E(S)----------------------
      ES_nu = matrix(0, nrow = l+2, ncol = n)
      #ES_nuÿһ????n??Ԫ??,??ʾi=1...n,????��????ΪES?ķ???
      #ES_nu
      #L1/beta1         L2/beta1     ...   Ln/beta1
      #L1/beta2         L2/beta2     ...   Ln/beta2
      #L1/Lambda{Yhat 1} L2/Lambda{Yhat 1} ... Ln/Lambda{Yhat 1}
      #...
      #L1/Lambda{Yhat l} L2/Lambda{Yhat l} ... Ln/Lambda{Yhat l}
      #
      #??1,2?зֱ????㣬??3:(m+2)??ֱ??????????

      Y_mat = matrix(rep(Y,l), nrow = l, byrow = TRUE)
      #Y1 Y2 ... Yn
      #Y1 Y2 ... Yn  l*n

      Y_hat_mat = matrix(rep(Y_hat,n), nrow = l, byrow = FALSE)
      #Yhat 1 ... Yhat 1
      #Yhat 2 ... Yhat 2
      #...
      #Yhat l ... Yhat l   l*n
      ES = c()
      ES_g = matrix(0, nrow = l+2, ncol = n)
      ES_de_c = rep(0,n)
      #ES????ÿһ?з?ĸ??ͬ????��?ۼ?ES_de_each??nά
      Lambda_b_nozero = lambdaY_best
      Lambda_b_nozero[which(Lambda_b_nozero==0)]=10^(-10) #??Lambda{Y}??0??Ϊ1,????Lambda{Y'}
      for(g in 1:Q) #??????:??��1-n
      {

        xi_i = nodes[g]
        ES_de_each = whts[g] * (xi_i* lambdaY_best* exp(beta_X))^(delta_i) *
          exp(-xi_i* LambdaY_best* exp(beta_X)) *
          dgamma_Q[g] *
          exp(xi_i)

        ES_nu[1:2,] =  ES_nu[1:2,] + t((X * (delta_i-xi_i*LambdaY_best*exp(beta_X)))*
                                         ES_de_each)
        ES_nu[3:(l+2),] = ES_nu[3:(l+2),] +
          t( ES_de_each *
               t(
                 matrix(rep(delta_i,l), nrow = l, byrow = TRUE)*
                   (Y_eq_Yhat)/
                   matrix(rep(Lambda_b_nozero,l), nrow = l, byrow = TRUE)-
                   matrix(rep(xi_i*exp(beta_X),l), nrow = l, byrow = TRUE)*
                   (Y_geq_Yhat)

               )
          )
        ES_de_c = ES_de_c + ES_de_each
      }

      ES_g = ES_nu / matrix(rep(ES_de_c,(l+2)), nrow = (l+2), byrow = TRUE)
      ESES = ES_g%*%t(ES_g)

      #---------------------E(SST),E(B)---------------

      list_ESSEB = ESSEB(l, n, Q, nodes, whts, lambdaY_best, LambdaY_best,
                         beta_X, delta_i, dgamma_Q, Lambda_b_nozero,
                         Y_eq_Yhat, Y_geq_Yhat, X)

      I = (list_ESSEB$EB-list_ESSEB$ESS+ESES)
      SE_square = diag(solve(I))
      SE = sqrt(SE_square)

      e1 = c(0,0,(Y_hat <= 1))
      e2 = c(0,0,(Y_hat <= 2))

      #??ά
      SE_Lambda = sqrt( diag( rbind(e1,e2) %*% solve(I) %*% cbind(e1,e2) ) )

      plot_x_m = matrix(rep(plot_x,m), nrow = n_x ,byrow = FALSE)
      plot_I = (plot_x_m)>=( matrix(rep(Y_dot,n_x), nrow = n_x, byrow = TRUE) )
      plot_Lam_Ydot_c = rowSums(plot_I * matrix(rep(lambdaYd_best,n_x), nrow = n_x, byrow = TRUE))
      plot_Lam_Ydot = rbind( plot_Lam_Ydot_c ,plot_Lam_Ydot)

      e_nx = cbind(matrix(0, nrow = n_x, ncol = 2),
                   (matrix(rep(Y_hat,n_x),nrow = n_x, byrow = TRUE )<= matrix(rep(plot_x, l), nrow = n_x, byrow = FALSE) ) )
      plot_SE = sqrt(diag(e_nx%*% solve(I) %*%t(e_nx)))




    }else if(alpha.true==0)#alpha=0
    {
      lamI =  rowSums(matrix(rep(c(1:2), m), nrow = 2, byrow = FALSE) >
                        matrix(rep(Y_dot,2), nrow = 2, byrow = TRUE) )

      dataframe = data.frame(Y = Y, delta = delta_i, X1 = X1, X2 = X2)

      res.cox = coxph(Surv(time = Y ,event = delta ) ~ X1 + X2, data =  dataframe)
      su.res.cox = summary(res.cox)

      sur = survfit(res.cox, newdata = data.frame(X1=0,X2=0))

      result_lambda_Ydot = rbind(c(sur$cumhaz, rep(-1,(n-m))), result_lambda_Ydot)
      result = rbind(c(su.res.cox$coefficients[,1],sur$cumhaz[lamI] ,1-sum(delta_i)/n),result)
      SE_Lambda = sur$std.err[lamI]
      SE = su.res.cox$coefficients[,3]

      plot_x_m = matrix(rep(plot_x,m), nrow = n_x ,byrow = FALSE)
      plot_I = (plot_x_m)>=( matrix(rep(Y_dot,n_x), nrow = n_x, byrow = TRUE) )
      plot_Lam_Ydot_c = rowSums(plot_I * matrix(rep(lambda_Ydot_new,n_x), nrow = n_x, byrow = TRUE))
      plot_Lam_Ydot = rbind( plot_Lam_Ydot_c ,plot_Lam_Ydot)

      plot_SE = sur$std.err[rowSums(plot_I)]

    }else
    {
      print("wrong alpha<0")
    }

    para_now = result[1,1:4]
    SE_now = c(SE[1:2], SE_Lambda)

    u_975 = qnorm(0.975)
    CP_now = (para_true >= (para_now-u_975*SE_now))&
      (para_true <= (para_now+u_975*SE_now))

    result_SE = rbind(SE_now, result_SE)
    result_CP = rbind(CP_now, result_CP)


    plot_upper = rbind(plot_Lam_Ydot_c + u_975*plot_SE, plot_upper)
    plot_lower = rbind(plot_Lam_Ydot_c - u_975*plot_SE, plot_lower)


    if(cu.rep){print(now_repeat)}
    if(cu.res)
    {
      result_till_now[,1] = bias_till_now = ((now_repeat-1)*bias_till_now + para_now - para_true )/now_repeat
      bias2_till_now = ((now_repeat-1)*bias2_till_now + (para_now - para_true)^2 )/now_repeat
      result_till_now[,2] = SD_till_now = sqrt( (now_repeat/(now_repeat-1)) * (bias2_till_now - bias_till_now^2) )
      # print(SD_till_now[1])
      # print(sd(result[1:now_repeat,1]))
      result_till_now[,3] = SE_till_now = ((now_repeat-1)*SE_till_now + SE_now )/now_repeat
      result_till_now[,4] = CP_till_now = ((now_repeat-1)*CP_till_now + CP_now )/now_repeat
      print(result_till_now)
    }



  }
  colnames(result) = c("beta1","beta2","lambda1","lambda2","censor")

  #-------------------------??????????Ҫ????---------------------------
  bias_beta1 = mean(result[1:num_of_repeat,1])+0.5
  bias_beta2 = mean(result[1:num_of_repeat,2])-1
  bias_Lambda1 = mean(result[1:num_of_repeat,3])-(rho*(log(1+1)))
  bias_Lambda2 = mean(result[1:num_of_repeat,4])-(rho*(log(1+2)))
  sd_beta1 = sd(result[1:num_of_repeat,1])
  sd_beta2 = sd(result[1:num_of_repeat,2])
  sd_Lambda1 = sd(result[1:num_of_repeat,3])
  sd_Lambda2 = sd(result[1:num_of_repeat,4])
  se_beta1 = mean(result_SE[1:num_of_repeat,1])
  se_beta2 = mean(result_SE[1:num_of_repeat,2])
  se_Lambda1 = mean(result_SE[1:num_of_repeat,3])
  se_Lambda2 = mean(result_SE[1:num_of_repeat,4])
  cp_beta1 = sum(result_CP[1:num_of_repeat,1])/num_of_repeat
  cp_beta2 = sum(result_CP[1:num_of_repeat,2])/num_of_repeat
  cp_Lambda1 = sum(result_CP[1:num_of_repeat,3])/num_of_repeat
  cp_Lambda2 = sum(result_CP[1:num_of_repeat,4])/num_of_repeat
  result_table = matrix(c(bias_beta1  , sd_beta1 , se_beta1 , cp_beta1 ,
                          bias_beta2  , sd_beta2 , se_beta2 , cp_beta2 ,
                          bias_Lambda1,sd_Lambda1,se_Lambda1,cp_Lambda1,
                          bias_Lambda2,sd_Lambda2,se_Lambda2,cp_Lambda2),nrow = 4,byrow = TRUE)
  colnames(result_table) = c("Bias","SD","SE","CP")
  rownames(result_table) = c("beta1","beta2","Lambda1","Lambda2")
  #????????
  show_result = round(result_table,3)
  show_result[,4] = show_result[,4]*100
  show_result
  paste(c("alpha=","rho=","n="),c(alpha.true,rho,n))

  #latex
  # for(p in 1:4)
  # {print(paste0("$",show_result[p,1],"$&$",show_result[p,2],"$&$",show_result[p,3],"$&$",show_result[p,4],"$"))}

  if(plot.Lamb)
  {
    #?Ƚ?Lambda?��?ֵ????ʵֵ
    plot(x = plot_x, y = colMeans(plot_Lam_Ydot[1:num_of_repeat,]), type = "l",lty=2, col = "blue", ylim = c(0,0.7),
         ylab = "Lambda(Y)", xlab = "Y", main = paste0("n=",n,", rho=",rho,", alpha=",alpha.true,", ",num_of_repeat," replicates"))
    lines(x = plot_x, y = rho*log(1+plot_x), col="black")
    lines(x = plot_x, y = colMeans(plot_upper[1:num_of_repeat,]), type = "l",lty=2, col = "red")
    lines(x = plot_x, y = colMeans(plot_lower[1:num_of_repeat,]), type = "l",lty=2, col = "red")
    legend("topleft",legend = c("true","estimate","confidence interval"),lty = c(2,1,2), col = c("blue", "black","red"))
  }

  # length(lik.fun)
  # plot(alpha.c, lik.fun)
  # lines(alpha.c, result_lik[1,])
  # length(alpha.c)
  # alpha.c[which(lik.fun == max(lik.fun))]

  return(list( main.est = result[1:num_of_repeat,], SE = result_SE[1:num_of_repeat,], CP = result_CP[1:num_of_repeat,],
               X1 = result_X1[1:num_of_repeat,], X2 = result_X2[1:num_of_repeat,],
               Y = result_Y[1:num_of_repeat,],  Ydot = result_Ydot[1:num_of_repeat,],
               m = result_m,  delta = result_delta[1:num_of_repeat,],
               lambda_Ydot = result_lambda_Ydot[1:num_of_repeat,], alpha = result_alpha))

}
