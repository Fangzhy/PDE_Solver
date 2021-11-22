
#This script is used to solve the chemical diffusion problem in skin layers involving different diffusivities in different direction.
#Fang, 12/31/2012
#(1) dc/dt = D_l*(d^2C/dr^2) + D_t*(d^2c/dz^2)
#IC: c(r,z,0) = c0;
#BC: 1) c(r,-L/2,t)=0, 2) c(r,L/2,t)=0, 3) D_l*(dc/dr|r=Rin)=-kC(Rin,z,t), 4)D_l*(dc/dr|r=Rout)=0  

library("ReacTran") #https://cran.r-project.org/web/packages/ReacTran/vignettes/PDE.pdf
# #### Method 1, assume zeros boundary conditions don't use boundary conditions of flux, 
# rm(list=ls())#clear the memory
# #treat x as radius of sweat duct r and y as skin depth Z
# pde2D <- function (t, C, parms) {
#   CONC <- matrix(nr = Nx, nc = Ny, C)
#   Tran <- tran.2D(CONC, D.x = Dx, D.y = Dy, C.x.up = 0,
#                   C.y.up = 0, C.y.down = 0, dx = Gridx, dy = Gridy)
#   dCONC <- Tran$dC - r * CONC
#   #dCONC[ii]<- dCONC[ii] + p
#   return(list(as.vector(dCONC)))
# }
# #
# Nx <- 100;
# Ny <- 100;
# dy <- dx <- 1
# Gridx <- setup.grid.1D(x.up = 10, x.down = 100, N = Nx)
# Gridy <- setup.grid.1D(x.up = -20, x.down = 20, N = Ny)
# Dy <- 0.5
# Dx <- 2.5
# r <- 0.1
# 
# #
# Cini <- rep(100, Nx*Ny)
# times <- 0:15
# print(system.time(
#   out <- ode.2D (y = Cini, times = times,
#                  parms = NULL, func = pde2D,
#                  dimens = c(Nx, Ny), method = "ode45")
# ))
# 
# #
# par(oma = c(0,0,1,0))
# image(out, xlab = "Sweat Radius", ylab = "Skin Depth Z",zlim=range(out),add.contour = TRUE,
#       mfrow = c(4, 4), ask = FALSE,
#       main = paste("t = ", times),
#       grid = list(x = Gridx$x.mid, y = Gridy$x.mid))
# mtext(side = 3, outer = TRUE, cex = 1.25, line = -1,
#       "PDE for chemical transport in skin layers and sweat duct")
# 
# 
# image(out, main = "PDE for chemical transport in skin layers and sweat duct",grid = list(x = Gridx$x.mid, y = Gridy$x.mid),add.contour = TRUE,
#       xlab = "Sweat Radius", ylab = "Skin Depth Z",zlim=range(out))

#### Method 2; use boundary conditions of flux on radius direction (i.e., x) 
rm(list=ls())#clear the memory
#treat x as radius of sweat duct r and y as skin depth Z
pde2D <- function (t, C, parms) {
  CONC <- matrix(nr = Nx, nc = Ny, C)
  Tran <- tran.2D(CONC, D.x = Dx, D.y = Dy, C.y.up = 0, C.y.down = 0, 
                  flux.x.down=0, flux.x.up=-k*CONC[1,], dx = Gridx, dy = Gridy)
  dCONC <- Tran$dC
  return(list(as.vector(dCONC)))
}
#
Nx <- 100;
Ny <- 50;
#dy <- dx <- 1
Gridx <- setup.grid.1D(x.up = 10, x.down = 100, N = Nx)
Gridy <- setup.grid.1D(x.up = -20, x.down = 20, N = Ny)
Dy <- 2
Dx <- 6
k <- 4

#
Cini <- rep(100, Nx*Ny)
times <- 0:31
print(system.time(
  out <- ode.2D (y = Cini, times = times,
                 parms = NULL, func = pde2D,
                 dimens = c(Nx, Ny), method = "ode45")
))

#
par(oma = c(0,0,1,0))
image(out, xlab = "Sweat Radius", ylab = "Skin Depth Z",zlim=range(out),add.contour = TRUE,
      mfrow = c(4, 4), ask = FALSE,
      main = paste("t = ", times),
      grid = list(x = Gridx$x.mid, y = Gridy$x.mid))
mtext(side = 3, outer = TRUE, cex = 1.25, line = -1,
      "PDE for chemical transport in skin layers and sweat duct")
#
image(out, main = "PDE for chemical transport in skin layers and sweat duct",grid = list(x = Gridx$x.mid, y = Gridy$x.mid),add.contour = TRUE,
      xlab = "Sweat Radius", ylab = "Skin Depth Z",zlim=range(out))
#Ratio of Wt to W0
dydx = as.matrix(Gridy$dx)%*%t(Gridx$dx)
dydx=as.vector(dydx);
W0=sum(Cini*dydx)
Wt=rep(0,1,length(times))
for (i in 1:length(times)){
  Wt[i]=sum(out[i,2:ncol(out)]*dydx)
}
plot(Wt/W0,xlab='Time', type= 'b', ylab='Wt / W0', main = 'Amount of chemicals remained in skin at each time step')
plot(1-Wt/W0,xlab='Time', type= 'b', ylab='1 - Wt / W0', main = 'Amount of chemicals desorped into water at each time step')
