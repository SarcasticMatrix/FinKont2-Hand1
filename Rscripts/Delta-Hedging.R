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

BlackScholesImpVol  <- function (obsprice,spot,timetomat,strike,r, q=0, opttype=1)
{ difference<- function(sigBS, obsprice,spot,timetomat,strike,r,q,opttype)
{BlackScholesFormula (spot,timetomat,strike,r,q,sigBS, opttype,1)-obsprice
}
uniroot(difference, c(10^-3,0.5),obsprice=obsprice,spot=spot,timetomat=timetomat,strike=strike,r=r,q=q,opttype=opttype,tol = .Machine$double.eps)$root
}

# INITIALIZE
# ==========

S0<-100
mu<-0.05; r<-0.05; sigma<-0.2; 

sigma_hedge<-0.2

T<-1; K<-102
capT <- T
Nhedge<-252; Nrep<-1000

opttype <- 2

# HEDGE
# =====

St<-rep(S0, length=Nrep)
dt<-T/Nhedge
#initialoutlay<-BlackScholesFormula(S0,capT,strike, r,0,sigma_hedge,1,1)
initialoutlay<-BlackScholesFormula(S0,T,K, r,0,sigma,opttype,1)

Vpf<-rep(initialoutlay,length=Nrep)

a<-BlackScholesFormula(St,T,K,r,0,sigma_hedge,opttype,2)
b<-Vpf-a*St

for(i in 2:Nhedge){
  St<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))	
  Vpf<-a*St+b*exp(dt*r)    
  a<- BlackScholesFormula(St,(capT-(i-1)*dt),K, r,0,sigma_hedge,opttype,2)
  b<-(Vpf-a*St)
}

ST<-St*exp((mu-0.5*sigma^2)*dt +sigma*sqrt(dt)*rnorm(Nrep))
Vpf<-a*ST+b*exp(dt*r) 
hedgeerror<-(Vpf-pmax(ST-K,0))
optionpayoff<-pmax(ST-K,0)

# SUMMARY STATS & GRAPHS
# ======================

print(paste("Initial investment =",round(initialoutlay,4)))
print(paste("Average discounted option payoff =",round(exp(-r*T)*mean(optionpayoff),4)))
print(paste("Average discounted portfolio value =",round(exp(-r*T)*mean(Vpf),4)))

x11()
plot(ST,Vpf,col="blue",xlab="S(T)",ylab="Value of hedge portfolio(T)",main="Discrete hedging of a call option",xlim=range(ST),ylim=range(optionpayoff))
text(min(ST),max(optionpayoff),paste("#hegde points =",Nhedge),adj=0)
points(sort(ST),pmax(K-sort(ST),0),type='l',lwd=3)