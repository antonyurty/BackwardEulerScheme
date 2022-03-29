require(somebm)

Time = (0:10000)/10000
N = length(Time)-1
H = 0.7
k1 = 1
k2 = 1

Euler = function(Time, k1, k2, BH, Y0){
  N = length(Time)-1
  dt = Time[2] - Time[1]
  Y = numeric(N+1)
  Y[1] = Y0
  dBH = BH[-1] - BH[-length(BH)]
  
  for(i in 1:N){
    B0 = Y[i] + dBH[i] + dt*(k1 - k2)
    B1 = -1 - dt*( k1 + k2 )
    B2 = - Y[i] - dBH[i]
    
    p = B1 - B2^2/3
    q = 2*B2^3/27 - B2*B1/3 + B0
    
    Q = (p/3)^3 + (q/2)^2
    
    alpha = ( -q/2 + sqrt(as.complex(Q)) )^(1/3)
    beta =  ( -q/2 - sqrt(as.complex(Q)) )^(1/3)
    x3 = -(alpha+beta)/2 - rt3*(alpha - beta)/2 - B2/3
    
    Y[i+1] = as.numeric(x3)
  }
  cbind(Time, Y)
}

BH1 = c(0, fbm(H, N-1))
BH2 = c(0, fbm(H, N-1))
BH3 = c(0, fbm(H, N-1))
BH4 = c(0, fbm(H, N-1))
BH5 = c(0, fbm(H, N-1))
BH6 = c(0, fbm(H, N-1))
BH7 = c(0, fbm(H, N-1))
BH8 = c(0, fbm(H, N-1))
BH9 = c(0, fbm(H, N-1))
BH10 = c(0, fbm(H, N-1))

start_time <- Sys.time()
Y1 = Euler(Time, 0.5, 0.5, BH1, 0)
Y2 = Euler(Time, 0.5, 0.5, BH2, 0)
Y3 = Euler(Time, 0.5, 0.5, BH3, 0)
Y4 = Euler(Time, 0.5, 0.5, BH4, 0)
Y5 = Euler(Time, 0.5, 0.5, BH5, 0)
Y6 = Euler(Time, 0.5, 0.5, BH6, 0)
Y7 = Euler(Time, 0.5, 0.5, BH7, 0)
Y8 = Euler(Time, 0.5, 0.5, BH8, 0)
Y9 = Euler(Time, 0.5, 0.5, BH9, 0)
Y10 = Euler(Time, 0.5, 0.5, BH10, 0)
end_time <- Sys.time()
(end_time - start_time)/10

plot(Y1, type = "l", col = "black", ylim = c(-1,1))
lines(Y2, type = "l", col = "red")
lines(Y3, type = "l", col = "blue")
lines(Y4, type = "l", col = "purple")
lines(Y5, type = "l", col = "green")
lines(Y6, type = "l", col = "brown")
lines(Y7, type = "l", col = "darkslategray")
lines(Y8, type = "l", col = "indianred4")
lines(Y9, type = "l", col = "lightsteelblue4")
lines(Y10, type = "l", col = "peru")

abline(-1,0, col = "gray", lty = 2)
abline(1,0, col = "gray", lty = 2)
