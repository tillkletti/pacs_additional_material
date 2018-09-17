# Analysis of the Meuse dataset
# This is a demo, you can call it via the command demo(Meuse)

rm(list = ls())
graphics.off()

library(fdaPDE)
library(car)
data("MeuseData")
data("MeuseBorder")

head(MeuseData)
head(MeuseBorder)

# Choose which variable we wish to predict
zinc <- MeuseData$zinc
# Its plot is
# plot3d(MeuseData$x,MeuseData$y,zinc, type = "s")

# Create a mesh from the data
mesh <- create.MESH.2D(nodes = MeuseData[,c(2,3)], segments = MeuseBorder, order = 1)
plot(mesh)

# Refine the mesh specifying the maximum area
mesh <- refine.MESH.2D(mesh, maximum_area=6000, delaunay=TRUE)
plot(mesh)

# Create the FEM basis object
FEM_basis <- create.FEM.basis(mesh)

# Choose a set of regularization coefficients
lambda = 10^seq(1,9,3)

# Do an isotropic smoothing
# smoothing_iso <- smooth.FEM.basis(observations=c(zinc,rep(NA,dim(mesh$nodes)[1]-length(zinc))), FEMbasis=FEM_basis, lambda = lambda, GCV=TRUE)
# 
# plot(smoothing_iso$fit.FEM)
# plot(smoothing_iso$PDEmisfit.FEM)
# Which lambda is the good one ?

# Do an anisotropic smoothing using the CPP code after defining the regularization PDE parameters
kappa <- matrix(0, 2, 2) # Initial value of K
PDE_parameters <- list(K = kappa, b = c(0, 0), c = 0)

smoothing_aniso <- aniso.smooth.FEM.PDE.basis(observations = c(zinc,rep(NA,dim(mesh$nodes)[1]-length(zinc))),FEMbasis = FEM_basis, lambda = lambda, PDE_parameters = PDE_parameters, covariates = c(MeuseData$copper,rep(NA,dim(mesh$nodes)[1]-length(MeuseData$copper))))

# Plot the result
plot(smoothing_aniso$fit.FEM)
# plot(smoothing_aniso$PDEmisfit.FEM)

# display the anisotropy numerically and graphically
anisotropyMatrix <- function(angle, intensity){
  Q <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)), byrow = TRUE, ncol = 2, nrow = 2)
  Sigma <- matrix(c(1,0,0,intensity), byrow = TRUE, nrow = 2, ncol = 2)/sqrt(intensity)
  return(Q%*%Sigma%*%solve(Q))
}

smoothing_aniso$kappa
anisotropyMatrix(smoothing_aniso$kappa[1],smoothing_aniso$kappa[2])
plot(mesh)
# The ellipse is a graphical representation of the estimated anisotropy
ellipse(c(mean(MeuseData$x),mean(MeuseData$y)), shape = anisotropyMatrix(smoothing_aniso$kappa[1],smoothing_aniso$kappa[2]), radius=200, center.pch = 20, col = "blue", center.cex = 1, lwd = 2, lty = 1)
# This result is very satisfying from the intuitive point of view. Indeed it seems obvious that the concentrations of some metal in a river should be very similar in two points that are connected by the current.

