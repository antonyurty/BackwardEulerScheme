################ Simulating mBm #################
D = function(x,y){
  num = sqrt(gamma(2*x+1)*gamma(2*y+1)*sin(pi*x)*sin(pi*y))
  denom = 2*gamma(x+y+1)*sin(0.5*pi*(x+y))
  
  num/denom
}

C = function(t,s, H){
  D(H(t), H(s))*( t^(H(t)+H(s)) + s^(H(t)+H(s)) - abs(t-s)^(H(t) + H(s)) )
}

plot( (0:1000)/1000, 0.2*sin(2*pi*((0:1000)/1000))+0.5, type = "l", ylab = "H(t)", xlab = "t")

H = function(t){
  0.2*sin(2*pi*t) +0.5
}


######## Computing covariance

Time = (0:10000)/10000
Covariance = matrix(nrow = length(Time)-1, ncol = length(Time)-1)

for(i in 1:(length(Time)-1)){
  for(j in 1:(length(Time)-1)){
    Covariance[i,j] = C(Time[i+1], Time[j+1], H)
  }
}

L = chol(Covariance)

save(L, file="Lmatrix.Rdata")

gen.BH = function(L, Time){
  xi = rnorm(length(Time)-1)
  c(0, t(L) %*% xi)
}

L = load("~/Lmatrix.Rdata")

L[1,1]

############### Simulating sandwich ###############

require(nleqslv)

BH1 = gen.BH(L, Time)
BH2 = gen.BH(L, Time)
BH3 = gen.BH(L, Time)
BH4 = gen.BH(L, Time)
BH5 = gen.BH(L, Time)
BH6 = gen.BH(L, Time)
BH7 = gen.BH(L, Time)
BH8 = gen.BH(L, Time)
BH9 = gen.BH(L, Time)
BH10 = gen.BH(L, Time)

plot(Time, BH3, type = "l")


b = function(x,t,z, dt){
  x - dt*(x-sin(10*t))^(-4) + dt*(sin(10*t)+2 - x)^(-4) - z
}

Euler = function(Time, BH, b, Y0){
  N = length(Time)-1
  dt = Time[2] - Time[1]
  Y = numeric(N+1)
  Y[1] = Y0
  
  dBH = BH[-1] - BH[-length(BH)]
  
  for(i in 1:N){
    Y[i+1] = nleqslv(Y[i], b, t = Time[i+1], z = Y[i] + dBH[i], dt = dt)$x
  }
  cbind(Time, Y)
}

start_time <- Sys.time()
Y1 = Euler(Time, BH1, b, 1)
Y2 = Euler(Time, BH2, b, 1)
Y3 = Euler(Time, BH3, b, 1)
Y4 = Euler(Time, BH4, b, 1)
Y5 = Euler(Time, BH5, b, 1)
Y6 = Euler(Time, BH6, b, 1)
Y7 = Euler(Time, BH7, b, 1)
Y8 = Euler(Time, BH8, b, 1)
Y9 = Euler(Time, BH9, b, 1)
Y10 = Euler(Time, BH10, b, 1)
end_time <- Sys.time()
(end_time - start_time)/10

plot(Y1, type = "l", ylim = c(-1,3))
lines(Y2, type = "l", col = "red")
lines(Y3, type = "l", col = "blue")
lines(Y4, type = "l", col = "purple")
lines(Y5, type = "l", col = "green")
lines(Y6, type = "l", col = "brown")
lines(Y7, type = "l", col = "darkslategray")
lines(Y8, type = "l", col = "indianred4")
lines(Y9, type = "l", col = "lightsteelblue4")
lines(Y10, type = "l", col = "peru")


lines(Time, sin(10*Time), col = "gray", lty = 2)
lines(Time, sin(10*Time)+2, col = "gray", lty = 2)



