BlackScholesFormula  <- function (
    spot,
    timetomat,
    strike,
    r, 
    q=0, 
    sigma, 
    opttype=1, 
    greektype=1)
{ 
  
  d1<-(log(spot/strike)+ ((r-q)+0.5*sigma^2)*timetomat)/(sigma*sqrt(timetomat))
  d2<-d1-sigma*sqrt(timetomat)
  
  if (opttype==1 && greektype==1) result<-spot*exp(-q*timetomat)*pnorm(d1)-strike*exp(-r*timetomat)*pnorm(d2)
  
  if (opttype==2 && greektype==1) result<-spot*exp(-q*timetomat)*pnorm(d1)-strike*exp(-r*timetomat)*pnorm(d2)-spot*exp(-q*timetomat)+strike*exp(-r*timetomat)
  
  if (opttype==1 && greektype==2) result<-exp(-q*timetomat)*pnorm(d1)
  
  if (opttype==2 && greektype==2) result<-exp(-q*timetomat)*(pnorm(d1)-1)
  
  if (greektype==3) result<-exp(-q*timetomat)*dnorm(d1)/(spot*sigma*sqrt(timetomat))
  
  if (greektype==4) result<-exp(-q*timetomat)*spot*dnorm(d1)*sqrt(timetomat)
  
  BlackScholesFormula<-result
  
}

# INITIALIZE
# ==========

# JP
S0<-30000
muf<-0.05
rf<-0.05
sigmaf<-0.2

# Exchange rate 
X0 <- 1/100
sigmaX<-0.2

# US 
rd <- 0.05
q <- rd - rf + sigmaX * sigmaf

sigma_hedge <- 0.2

T<-2
K<-S0
capT <- T
Nhedge<-252; Nrep<-1000

opttype <- 2

# HEDGE
# =====

St<-rep(S0, length=Nrep)
Xt<-rep(X0, length=Nrep)
dt<-T/Nhedge

# V(0) 
initialoutlay <- X0 * BlackScholesFormula(S0,T,K,rd,q,sigmaf,opttype,1)
initialoutlay / X0
V <- rep(initialoutlay,length=Nrep)

# h_Bf(0)
h_Bf <- X0 * BlackScholesFormula(St,T,K,rd,q,sigma_hedge,opttype,2) / Xt * St

# h_Bd(0)
h_Bd <- V - h_Bf * Xt

for(i in 2:Nhedge){
  
  # Update the processes
  BM = rnorm(Nrep)
  St <- St * exp( (muf-0.5 * sigmaf^2) * dt + sigmaf * sqrt(dt) * BM)	
  Xt <- Xt * exp( (rd-rf-0.5 * sigmaX^2) * dt + sigmaX * sqrt(dt) * BM)
  
  # V(t)
  V <- h_Bf * exp(rf * dt) * Xt + h_Bd * exp(rd * dt)

  # h_Bf(t)
  h_Bf <- X0 * BlackScholesFormula(St,T,K,rd,q,sigma_hedge,opttype,2) / Xt * St

  # h_Bd(t)
  h_Bd <- V - h_Bf * Xt
}

BM <- rnorm(Nrep)
ST <- St * exp((muf-0.5*sigmaf^2) * dt + sigmaf * sqrt(dt) * BM)	
XT <- Xt * exp((rd-rf-0.5*sigmaX^2) * dt + sigmaX * sqrt(dt) * BM)

V <- h_Bf * exp(rf * dt) * XT + h_Bd * exp(rd * dt)

hedgeerror <- V-X0*pmax(K-ST,0)
optionpayoff <- X0*pmax(K-ST,0)

# SUMMARY STATS & GRAPHS
# ======================

print(paste("Initial investment =",round(initialoutlay,4)))
#print(paste("Average discounted option payoff =",round(exp(-r*T)*mean(optionpayoff),4)))
#print(paste("Average discounted portfolio value =",round(exp(-r*T)*mean(Vpf),4)))

x11()
plot(ST,V,col="blue",xlab="S(T)",ylab="Value of hedge portfolio(T)",main="Discrete hedging of a Quanto-put option")
points(sort(ST),X0*pmax(K-sort(ST),0),type='l',lwd=3,col='black')

