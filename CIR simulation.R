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
  denom = 2*(1+k2*dt)
  for(i in 1:N){
    num1 = Y[i] + dBH[i]
    num2 = num1^2 + 4*k1*dt*(1+k2*dt)
    num = num1 + sqrt( num2  )
    
    Y[i+1] = num/denom
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
Y1 = Euler(Time, 1, 1, BH1, 1)
Y2 = Euler(Time, 1, 1, BH2, 1)
Y3 = Euler(Time, 1, 1, BH3, 1)
Y4 = Euler(Time, 1, 1, BH4, 1)
Y5 = Euler(Time, 1, 1, BH5, 1)
Y6 = Euler(Time, 1, 1, BH6, 1)
Y7 = Euler(Time, 1, 1, BH7, 1)
Y8 = Euler(Time, 1, 1, BH8, 1)
Y9 = Euler(Time, 1, 1, BH9, 1)
Y10 = Euler(Time, 1, 1, BH10, 1)
end_time <- Sys.time()
(end_time - start_time)/10


plot(Y1, type = "l", col = "black", ylim = c(0,2.5))
lines(Y2, type = "l", col = "red")
lines(Y3, type = "l", col = "blue")
lines(Y4, type = "l", col = "purple")
lines(Y5, type = "l", col = "green")
lines(Y6, type = "l", col = "brown")
lines(Y7, type = "l", col = "darkslategray")
lines(Y8, type = "l", col = "indianred4")
lines(Y9, type = "l", col = "lightsteelblue4")
lines(Y10, type = "l", col = "peru")

abline(0,0, col = "gray", lty = 2)
