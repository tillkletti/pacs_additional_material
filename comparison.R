setwd('/mnt/Partage/Politecnico/PACS/project/')
library(parallel)
library(fdaPDE)
graphics.off()

n_iterazioni <- 30

res <- mclapply(1:n_iterazioni, function (iterazione, RMSE_r, RMSE_cpp, RMSE_star_r, RMSE_star_cpp) {
# res <- lapply(1:n_iterazioni, function (iterazione, RMSE_r, RMSE_cpp, RMSE_star_r, RMSE_star_cpp) {
  load(paste0('smoothing_aniso_cpp', iterazione, '.RData'))
  load(paste0('smoothing_iso_aniso', iterazione, '.RData'))
  
  data_train <- read.table(paste(getwd(),"/random_fields/random_field_train_",iterazione,".txt",sep=''), header=FALSE)
  data_locations <- matrix(NA,nrow=dim(data_train)[1],ncol=2)
  data_locations[,1] <- as.numeric(data_train[,1])
  data_locations[,2] <- as.numeric(data_train[,2])
  
  data_test  <- read.table(paste(getwd(),"/random_fields/random_field_test_",iterazione,".txt",sep=''),  header=FALSE)
  test_data <- as.numeric(data_test[,3])
  test_locations <- matrix(NA,nrow=dim(data_test)[1],ncol=2)
  test_locations[,1] <- as.numeric(data_test[,1])
  test_locations[,2] <- as.numeric(data_test[,2])
  
  fitFEM_r <- smoothing_final_aniso$fit.FEM
  fitFEM_cpp <- smoothing_aniso_cpp$fit.FEM
  
  RMSE_r <- sqrt(sum(( eval.FEM(fitFEM_r, data_locations) - data )^2)/length(data))
  RMSE_cpp <- sqrt(sum(( eval.FEM(fitFEM_cpp, data_locations) - data )^2)/length(data))
  
  RMSE_star_r <- sqrt(sum(( eval.FEM(fitFEM_r, test_locations) - test_data )^2)/length(test_data))
  RMSE_star_cpp <- sqrt(sum(( eval.FEM(fitFEM_cpp, test_locations) - test_data )^2)/length(test_data))
  
  return(c(RMSE_r, RMSE_cpp, RMSE_star_r, RMSE_star_cpp))
}, mc.cores = getOption("mc.cores", 4L))
# })

RMSE <- as.data.frame(t(simplify2array(res, higher = FALSE)))
colnames(RMSE) <- paste('RMSE', c('(R script)', '(C++ code)', 'star (R script)', 'star (C++ code)'))
rm(res)

X11()
boxplot(RMSE, names = colnames(RMSE), ylim = c(0, max(RMSE)))
summary(RMSE)
