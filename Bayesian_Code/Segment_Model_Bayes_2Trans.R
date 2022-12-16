###############################
###   Bayesian Model    ####
###        For NPQ Transition model          ####
###############################
Seg_model_2Trans <-
  "model {
### Data Loop
for (t in 1:N){
NPQ[t] ~ dnorm( mu.NPQ[t] , tau )
}
 for(i in 1:N)
  {
mu.NPQ[i] <- Intercept[K[i]] + (slope[J[i]+K[i]] * (Time[i]- Trans[K[i]]))
J[i] <- step(Time[i] - Trans[1])
K[i] <- 1 + step(Time[i] - Trans[2])
}
Intercept[1] ~ dnorm(0.0,10)T(0,)#~ dnorm(0.5,4)#T(0,)
Intercept[2] ~ dnorm(0.0,10)T(0,)#~ dnorm(0.08,4)#T(0,)

slope[1] ~ dnorm(0.0,1e2)T(,0)#dnorm(-0.01,1000)#T(,0)
slope[2]  <- (Intercept[2] - Intercept[1])/(Trans[2]-Trans[1])
slope[3] ~  dnorm(0.0,1e8)T(,0)#dnorm(-4.33333E-05,7e9)#T(,0)


### Trans Est full Pulse sequence 
TransT[1] ~ dnorm(1000,1.6e-05)
TransT[2] ~ dnorm(2000,1.6e-05)
#TransT[1] ~ dunif(minX,maxX-1)
#TransT[2] ~ dunif(minX+1,maxX)

Trans[1:2] <- sort(TransT)

### calcs
yi = (slope[2]*Trans[1]) + Intercept[2]
Intercepta = yi - (slope[1]*Trans[1])

QE = (Intercepta - Intercept[1]) / Intercepta
QT = Intercept[1] / Intercepta


tau ~ dgamma(1 ,1)
}
" 
