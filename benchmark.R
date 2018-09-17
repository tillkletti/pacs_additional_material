setwd('/mnt/Partage/Politecnico/PACS/project/')

rm(list=ls())
graphics.off()

library(fdaPDE)

sigmanoise <- 0.1 # deviazione standard dell'errore gaussiano

# sequenza di rho per GCV
rhoseq <- 10^seq(-7,-0.01,by=0.1)

# sequenza di rho per ottimizzazione iterativa
rhovec <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
rhovec_l <- length(rhovec)

# dominio
lato <- 5
x_L <- y_L <- -lato
x_U <- y_U <- lato
limite <- 1
len_bordo <- length(seq(from=x_L,to=x_U,by=limite))
boundary <- cbind(c(seq(from=x_L,to=x_U,by=limite)[-1],rep(x_U,len_bordo)[-1], seq(from=x_U,to=x_L,by=-limite)[-1],rep(x_L,len_bordo)[-1]), 
                  c(rep(y_L,len_bordo)[-1],seq(from=y_L,to=y_U,by=limite)[-1],rep(y_U,len_bordo)[-1],seq(from=y_U,to=y_L,by=-limite)[-1]))

n_iterazioni <- 4
# array dove salver? i risultati

result_optim_LBFGSB <- array(NA,dim=c(rhovec_l,2,n_iterazioni)) # angolo e intensit?

GCV_optim <- array(NA,dim=c(rhovec_l,n_iterazioni))
GCV_smooth <- array(NA,dim=c(rhovec_l,n_iterazioni))

rho_smooth <- array(NA,dim=c(rhovec_l,n_iterazioni))

index <- rep(NA,n_iterazioni)
angle_final <- rep(NA,n_iterazioni)
intensity_final <- rep(NA,n_iterazioni)
rho_smooth_final <- rep(NA,n_iterazioni)

# times
times <- matrix(NA_real_, n_iterazioni, 2, dimnames = list(c(), c("R", "CPP")))

for(iterazione in 1:n_iterazioni){
  #####
  data_train <- read.table(paste(getwd(),"/random_fields/random_field_train_",iterazione,".txt",sep=''), header=FALSE)
  data_test  <- read.table(paste(getwd(),"/random_fields/random_field_test_",iterazione,".txt",sep=''),  header=FALSE)
  
  set.seed(iterazione)
  data <- as.numeric(data_train[,3]) + rnorm(length(data_train[,3]), 0, sigmanoise*diff(range((as.numeric(data_train[,3]))))) # aggiungo rumore gauusiano ai dati
  
  data_locations <- matrix(NA,nrow=dim(data_train)[1],ncol=2)
  data_locations[,1] <- as.numeric(data_train[,1])
  data_locations[,2] <- as.numeric(data_train[,2])
  
  p <- matrix(data=NA,nrow=length(data)+dim(boundary)[1],ncol=2)
  p[,1] <- c(data_locations[,1],boundary[,1])
  p[,2] <- c(data_locations[,2],boundary[,2])
  
  isboundary <- matrix(data=NA,nrow=dim(boundary)[1],ncol=2)
  isboundary[,1] <- (length(data)+1):(length(data)+dim(boundary)[1])
  isboundary[,2] <- c((length(data)+2):(length(data)+dim(boundary)[1]),(length(data)+1))
  
  mesh_1 <- create.MESH.2D(p, order = 1, segments = isboundary)
  mesh <- refine.MESH.2D(mesh_1, maximum_area=0.2, delaunay=TRUE)
  
  basisobj <- create.FEM.basis(mesh)
  
  n <- length(data) # numero di dati
  omega <- (lato*2)^2 # area dominio
  
  ##### 
  ##### CPP anisotropic smoothing
  kappa <- matrix(0, 2, 2)
  PDE_parameters <- list(K = kappa, b = c(0, 0), c = 0)
  times[iterazione,2] <- system.time(smoothing_aniso_cpp <- aniso.smooth.FEM.PDE.basis(observations=c(data,rep(NA,dim(mesh$nodes)[1]-length(data))), FEMbasis=basisobj, lambda=rhovec/(1-rhovec)*n/omega, GCVmethod = 1, PDE_parameters = PDE_parameters))[["elapsed"]]

  #####
  ##### R anisotropic smoothing
  times[iterazione, 1] <- system.time({
    function_to_optimize = function(x){

      angle <- x[1]
      intensity <- x[2]

      Q <- cbind(c(cos(angle),sin(angle)),c(-sin(angle),cos(angle)))
      sigma <- matrix(c(1,0,0,intensity)/sqrt(intensity),nrow=2,ncol=2,byrow=TRUE)
      kappa <- Q%*%sigma%*%solve(Q)

      PDE_parameters <- list(K = kappa, b = c(0, 0), c = 0)

      smoothing <- smooth.FEM.PDE.basis(observations=c(data,rep(NA,dim(mesh$nodes)[1]-length(data))), FEMbasis=basisobj, lambda= rho/(1-rho)*n/omega, PDE_parameters = PDE_parameters)

      H <- sum(( eval.FEM(smoothing$fit.FEM, data_locations) - data )^2) # MSE sul training set

      return(H)
    }

    for(i in 1:rhovec_l){

      rho <- rhovec[i]

      if(i==1) start <- c(pi/2,5)
      if(i!=1) start <- result_optim_LBFGSB[i-1,,iterazione]

      # optime <- system.time({
        result_optim_LBFGSB[i,,iterazione] <- optim(par = start, fn = function_to_optimize, method = "L-BFGS-B", lower = c(0,1), upper = c(pi,1000))$par

        # tengo conto della periodicit? dell'angolo:
        if(result_optim_LBFGSB[i,1,iterazione]==pi) result_optim_LBFGSB[i,,iterazione] <- optim(par = c(0 ,result_optim_LBFGSB[i,2,iterazione]), fn = function_to_optimize, method = "L-BFGS-B", lower = c(0,1), upper = c(pi,1000))$par
        if(result_optim_LBFGSB[i,1,iterazione]==0)  result_optim_LBFGSB[i,,iterazione] <- optim(par = c(pi,result_optim_LBFGSB[i,2,iterazione]), fn = function_to_optimize, method = "L-BFGS-B", lower = c(0,1), upper = c(pi,1000))$par
      # })["elapsed"]
      # print(paste0("(R) Time to optimize [", i, "]: ", optime))
      ## GCV
      angle <- result_optim_LBFGSB[i,1,iterazione]
      intensity <- result_optim_LBFGSB[i,2,iterazione]

      # print(paste0('(angle[', i, '], intensity[', i, '] = (', angle, ', ', intensity, ')'))

      Q <- cbind(c(cos(angle),sin(angle)),c(-sin(angle),cos(angle)))
      sigma <- matrix(c(1,0,0,intensity)/sqrt(intensity),nrow=2,ncol=2,byrow=TRUE)
      kappa <- Q%*%sigma%*%solve(Q)

      PDE_parameters <- list(K = kappa, b = c(0, 0), c = 0)
      #smoothing

      # time_to_GCV <- system.time({
        smoothing_aniso <- smooth.FEM.PDE.basis(observations=c(data,rep(NA,dim(mesh$nodes)[1]-length(data))), FEMbasis=basisobj, lambda= rhoseq/(1-rhoseq)*n/omega, PDE_parameters = PDE_parameters, GCV=TRUE, GCVmethod = 1)
      GCVseq <- smoothing_aniso$GCV
      minGCVid <- which.min(GCVseq)
      # cat(GCVseq, sep = ", ")
      # })["elapsed"]
      # print(paste0("(R) Time to compute GCV [", i, "]: ", time_to_GCV))

      rho_smooth[i,iterazione] <- rhoseq[minGCVid]
      GCV_smooth[i,iterazione] <- GCVseq[minGCVid]
    }
    # time_to_final <- system.time({
    index[iterazione] <- which.min(GCV_smooth[,iterazione])

    angle_final[iterazione] <- result_optim_LBFGSB[index[iterazione],1,iterazione]
    intensity_final[iterazione] <- result_optim_LBFGSB[index[iterazione],2,iterazione]

    Q <- cbind(c(cos(angle_final[iterazione]),sin(angle_final[iterazione])),c(-sin(angle_final[iterazione]),cos(angle_final[iterazione])))
    sigma <- matrix(c(1,0,0,intensity_final[iterazione])/sqrt(intensity_final[iterazione]),nrow=2,ncol=2,byrow=TRUE)
    kappa <- Q%*%sigma%*%solve(Q)

    # print(paste0('Final Kappa = (', angle_final[iterazione], ', ', intensity_final[iterazione]))
    # print(paste0('Final Lambda = ', rho_smooth[index[iterazione],iterazione]/(1-rho_smooth[index[iterazione],iterazione])*n/omega))

    PDE_parameters <- list(K = kappa, b = c(0, 0), c = 0)

    smoothing_final_aniso <- smooth.FEM.PDE.basis(observations=c(data,rep(NA,dim(mesh$nodes)[1]-length(data))), FEMbasis=basisobj, lambda= rho_smooth[index[iterazione],iterazione]/(1-rho_smooth[index[iterazione],iterazione])*n/omega, PDE_parameters = PDE_parameters)
    # })["elapsed"]
    # print(paste0("(R) Time to compute final regression: ", time_to_final))
  })["elapsed"]
  
  save(smoothing_aniso_cpp, file = paste("smoothing_aniso_cpp",iterazione,".RData",sep=''))
  save(data, smoothing_final_aniso, file = paste("smoothing_iso_aniso",iterazione,".RData",sep=''))
}
save(times, file = "time_benchmark")