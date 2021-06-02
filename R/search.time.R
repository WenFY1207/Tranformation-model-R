search.Lamb= function(res.tran, pred)
{
  len.t = length(pred)
  len.Ydot = length(res.tran$Ydot)
  index.t = rowSums( matrix(rep(res.tran$Ydot,len.t),nrow = len.t ,byrow = TRUE)
                     <=matrix(rep(pred, len.Ydot), nrow=len.t,  byrow = FALSE) )
  Lambda = res.tran$Lamb.est[index.t]
  pvalue = res.tran$p.Lambda[index.t]
  upper = res.tran$Lamb.upper[index.t]
  lower = res.tran$Lamb.lower[index.t]
  return(list(Lambda = Lambda, pvalue = pvalue,
              upper = upper, lower = lower))
}
