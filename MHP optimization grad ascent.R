#MHP gradient ascent 

#initialize
init<-list() #mu, alpha, beta
init[[1]] <-c(0.2,0.6)
init[[2]]<-matrix(c(0.05,0.4,0.65,0.3),2,2, byrow = T)
init[[3]]<-matrix(c(1,1,1,0.7),2,2, byrow = T)

#projection

proj_matrix<-function(x){
  
  y=x #dummy
  p=dim(x)[1]
  q=dim(x)[2]
  for(i in 1:p){
    for(j in 1:q){
      y[i,j] = max(x[i,j],0)
    }
  }
  return(y)
}


proj_vector <-function(x){
  
  y=x
  for(i in 1:length(x)){
    y[i] = max(x[i],0)
  }
  
  return(y)
}


#projected grad ascent
t1= Sys.time()
proj_grad_asc<-function(arrivals,init, max_iter=100, step_size=10^{-3}, tolerance = 10^{-4}){
  vals = numeric()
  c=0
  curr = init
  del_norm = sqrt(sum(dl_dmu(arrivals, curr[[1]], curr[[2]], curr[[3]])^2))+norm(dl_dalpha(arrivals, curr[[1]], curr[[2]], curr[[3]]))+dl_dbeta(arrivals, curr[[1]], curr[[2]], curr[[3]])
  
  
  while(del_norm>tolerance && c < max_iter){
  
    update = list()
    update[[1]] = proj_vector(curr[[1]]+step_size*dl_dmu(arrivals, curr[[1]],curr[[2]],curr[[3]]))
    update[[2]] = proj_matrix(curr[[2]]+step_size*dl_dalpha(arrivals, curr[[1]],curr[[2]],curr[[3]]))
    update[[3]] = proj_matrix(curr[[3]]+step_size*dl_dbeta(arrivals, curr[[1]],curr[[2]],curr[[3]])) 
    curr = update
    c = c+1
   
    del_norm = sqrt(sum(dl_dmu(arrivals, curr[[1]], curr[[2]], curr[[3]])^2))+norm(dl_dalpha(arrivals, curr[[1]], curr[[2]], curr[[3]]))+dl_dbeta(arrivals, curr[[1]], curr[[2]], curr[[3]])
    vals <- append(vals,loglikelihood_full(arrivals, curr[[1]],curr[[2]],curr[[3]]))
    print(c)
    print(paste0("Time elapsed is ", Sys.time()-t1))
  }
  
  return(list(curr, vals))
}


##running
optimisation<- proj_grad_asc(arrivals, init, step_size = 10^{-4})
plot(optimisation[[2]], xlab = "iteration", ylab = "log-likleihood", ty="l", lwd=2, col="blue")
optimisation[[1]] #parameters
#
beep()
t2= Sys.time()
t2-t1
  