
par(mfrow=c(1,2))

# likelihood
NUM_tweets = length(arrivals[[1]])
print(NUM_tweets)
T

mu_start<-torch_tensor(c(0.4),requires_grad = TRUE)
alpha_start<- torch_tensor(matrix(c(4.2),1,1), requires_grad = TRUE)
beta_start<-torch_tensor(matrix(c(4.7),1,1), requires_grad = TRUE)

t1=Sys.time()

#list of learning rates
num_iterations <- 100
#lr_list <- c(10^{-3}, 10^{-2}, 10^{-2}*5,10^{-1})
lr_list <- c( 1)


for(lr in lr_list){

#lr <- 10^{-3}
tolerance = 10^{-3}


optimizer <- optim_adadelta(list(mu_start,alpha_start, beta_start), lr)
my_plot = numeric()
parameter_error = as.numeric(torch_norm(mu_0$sub(mu_start))+torch_norm(alpha$sub(alpha_start))+torch_norm(beta$sub(beta_start)))

for (i in 1:num_iterations) {
  
  optimizer$zero_grad() #zero out any existing gradients 
  
  value <- -loglikelihood_full(arrivals,torch_clamp(mu_start, min=0),torch_clamp(alpha_start, min=0),torch_clamp(beta_start,min=0)) # we need to minimize neg log likleihood
 # print(paste0("Loglikelihood is: ", as.numeric(-value)))
  my_plot <- append(my_plot,  as.numeric(-value))
  
  
  if (i %% 10 == 0)  cat("Value is: ", as.numeric(value), "\n")
 
  
  value$backward()
  optimizer$step()
  

  # if(1){
  #   cat("Gradient for mu is", as.numeric(mu_start$grad),"\n\n")
  #   cat("Gradient for alpha is: ", as.matrix(alpha_start$grad), "\n\n")
  #   cat("Gradient for beta is: ", as.matrix(beta_start$grad), "\n\n")
  # }

  if (i %% 10 == 0) {
  cat("Iteration: ", i, "\n")
  print(mu_start)
  print(alpha_start)
  print(beta_start)
  }
 
  print(paste0("Time elapsed is ", Sys.time()-t1, " after iteration : ", i ))
  
  error = (as.numeric(torch_norm(mu_0$sub(mu_start))+torch_norm(alpha$sub(alpha_start))+torch_norm(beta$sub(beta_start))))/3
  error = sqrt(error) #RMSised 
  parameter_error <- append(parameter_error, error)
  
  grad_norm = torch_norm(alpha_start)+torch_norm(beta_start)+torch_norm(mu_start)
  if(as.numeric(grad_norm) < tolerance) print("gradient very small - tolerance reached")
    
  if(as.numeric(grad_norm) < tolerance)  break # early exit
    
    
  }
  
plot(1:length(my_plot), my_plot, lwd=2, ty="l", col="blue", xlab ="Iteration", ylab = "Loglikelihood", main = paste0("lr = ", lr, " and (mu. alpha. beta) = (0.5, 4, 5)"))
plot(1:length(parameter_error), parameter_error, lwd=2, ty="l", col="red", xlab ="Iteration", ylab = "Error", main = paste0("lr = ", lr,  " and (mu. alpha. beta) = (0.5, 4, 5)"))
beep()
}






