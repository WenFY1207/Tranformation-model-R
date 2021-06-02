#' @export
#'
generate_data = function(n, alpha, rho, beta_true, now_repeat=1)
{

  set.seed(7+700*(now_repeat+1))
  X1 = rbinom(n,1,0.5)
  epsil = rnorm(n)
  X2 = X1 + epsil * (abs(epsil) < 1 ) + (abs(epsil) >= 1 )
  X = matrix(c(X1,X2) , ncol = 2, byrow = FALSE)
  ##生成T,  代码t##
  U = runif(n)
  beta_X = c(X %*% beta_true) #向量beta*向量X ，所有对象的值组成向量

  if(alpha == 0)
  {
    t = (exp(-log(U)/(rho * exp(beta_X))) - 1)
  }else
  {
    t=(exp((U^(-alpha) - 1)/(alpha * rho * exp(beta_X))) - 1)
  }

  ##生成截断时间C， 观测数据Y ##
  C = runif(n,2,6)
  Y = pmin(C,t)
  delta_i = ifelse( C >= t, 1, 0)
  return(list(X = X, delta = delta_i, beta_X = beta_X, Y = Y))
}
