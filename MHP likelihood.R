#MHP likelihood
#R function, as defined in @Yuenda Chen et al ( uploaded in Dropbox) 
#R works with tensors, checked

R <- function(m,n,k, arrivals, mu_0, alpha, beta){
  
  if(k==0) return(0)
  else{
    
    t2=arrivals[[m]][[k]] #node m arrivals
   # print(c(m,n,k))
    if(k==1) t1=0
    else t1=arrivals[[m]][[k-1]]
    node_n = as.numeric(arrivals[[n]])
   # print(node_n)
    node_n_arrivals = node_n[node_n > t1 & node_n <t2]
    r = exp(-beta[m,n]*(t2-node_n_arrivals))
    r_dash <- sum(r)
   # print(paste0("r_dash is ",r_dash))
    return(r_dash+exp(-beta[m,n]*(t2-t1))*R(m,n,k-1,arrivals, mu_0, alpha, beta))
  }
  
  
}

# terms of type alpha_{ij}*sum(1-exp(-beta_mn*()))
term <- function(m,n,arrivals, mu_0, alpha, beta){
  v = T-as.numeric(arrivals[[n]])
  t = 1-exp(-beta[m,n]*v)
  return(alpha[m,n]*sum(t)/beta[m,n])
}

#for node m
#returns a vector, of the length same as number of arrivals in node m
logarithmic_part <- function(m,arrivals, mu_0, alpha, beta){
    #log term
    a = as.numeric(arrivals[[m]])
    array=torch_rand(length(a))
    for(k in 1:length(a)){
      
    arg = 0
    for(j in 1:N){
      
      arg = arg+ torch_mul(alpha[m,j],R(m,j,k,arrivals, mu_0, alpha, beta))
        
      }
   # r_vec = sapply(1:N, FUN = function(x){return(R(m,x,k, arrivals, mu_0, alpha, beta))}, simplify = ) #gives a list
  #  arg = mu_0[m]+torch_sum(torch_mul(alpha[m,],r_vec))
    array[k]= arg+mu_0[m]
    }
    return(array)
}


node_loglik<-function(m,arrivals, mu_0, alpha, beta){
#exp terms
  p=0
  for(j in 1:N){
    p=p+term(m,j,arrivals, mu_0, alpha, beta)
  }
  
#add log terms and -mu*T
  return(-mu_0[m]*T-p+sum(log(logarithmic_part(m,arrivals, mu_0, alpha, beta))))
  
}


#alpha, beta, mu 
# we know n = number of nodes and arrival data 
#wanna ket l(alpha, beta, mu)

loglikelihood_full <-function(arrivals,mu_0,alpha, beta){
  l = torch_rand(N)
  for(i in 1:N){
    l[i] = node_loglik(i,arrivals, mu_0, alpha, beta)
    cat("node ", i, "done")
  }
  return(sum(l))
}

# 
# t1<-Sys.time()
# # lieklihood
# L = loglikelihood_full(arrivals,mu_0,alpha, beta)
# L$backward()
# L$grad_fn
# 
# alpha$grad
# beta$grad
# mu_0$grad
# 
# t2<-Sys.time()
# beepr::beep()
# t2-t1
