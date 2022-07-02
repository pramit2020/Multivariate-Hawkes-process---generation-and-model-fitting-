#gradients 

#dl_dmu
dl_dmu <- function(arrivals, mu_0,alpha, beta){
  
  grad = numeric(N)
  for(i in 1:N){
    grad[i] = -T+sum(1/logarithmic_part(i,arrivals,mu_0,alpha,beta))
  }
  
  return(grad)
}


#dl_d(alpha)

dl_dalpha<-function(arrivals, mu_0,alpha, beta){
  
  grad_alpha = matrix(0,N,N)
  for(i in 1:N){
    n_arr = as.numeric(arrivals[[i]])
    for(j in 1:N){
      a = term(i,j,arrivals,mu_0,alpha, beta)/alpha[i,j]
      b_vec = sapply(c(1:length(n_arr)), FUN = function(x){R(i,j,x,arrivals,mu_0,alpha,beta)})
      grad_alpha[i,j] <- -a+sum(b_vec/logarithmic_part(i,arrivals,mu_0,alpha,beta))
    }
  }
  
  return(grad_alpha)
}

#auxiliary function for dl_dbeta
dR_dbeta <- function(m,n,k, arrivals, mu_0, alpha, beta){
  
  if(k==0) return(0)
  else{
    
    t2=arrivals[[m]][[k]] #node m arrivals
    if(k==1) t1=0
    else t1=arrivals[[m]][[k-1]]
    node_n = as.numeric(arrivals[[n]])
    node_n_arrivals = node_n[node_n > t1 & node_n <t2]
    #non null
    if(length(node_n_arrivals)>0){
    r1 = sum(censored_exp(-beta[m,n]*(t2-node_n_arrivals))*(node_n_arrivals-t2))
    r2=censored_exp(-beta[i,j]*(t2-t1))*(t1-t2)*R(m,n,k,arrivals,mu_0,alpha,beta)
    r3= censored_exp(-beta[i,j]*(t2-t1))*dR_dbeta(m,n,k-1,arrivals, mu_0, alpha, beta)
    return(r1+r2+r3)
    }
    
    else{
      r2=censored_exp(-beta[i,j]*(t2-t1))*(t1-t2)*R(m,n,k,arrivals,mu_0,alpha,beta)
      r3= censored_exp(-beta[i,j]*(t2-t1))*dR_dbeta(m,n,k-1,arrivals, mu_0, alpha, beta)
      return(r2+r3)
    }
    
    
  }
  
}


dl_dbeta<-function(arrivals, mu_0, alpha, beta){
  
  grad_beta = matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:N){
    n_arr = as.numeric(arrivals[[i]])
    a = term(i,j,arrivals,mu_0,alpha, beta)/alpha[i,j]
    b_vec = sapply(c(1:length(n_arr)), FUN = function(x){dR_dbeta(i,j,x,arrivals,mu_0,alpha,beta)})
    grad_beta[i,j] <- a/beta[i,j]+sum(b_vec/logarithmic_part(i,arrivals,mu_0,alpha,beta))*alpha[i,j]
    
    ## mimicking the term thingy
    v = T-as.numeric(arrivals[[j]])
    t = censored_exp(-beta[i,j]*v) #vector
    grad_beta[i,j]= grad_beta[i,j] + alpha[i,j]*sum(t*v)/beta[i,j]
     
    }
  }
  
  return(grad_beta)
  
}
