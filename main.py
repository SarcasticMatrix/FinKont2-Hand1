from scripts.BlackScholes import BlackScholes

import numpy as np

########################################
############## Parameters ##############
########################################
# US Stock
rUS = 0.05

# FX
X0 = 1/100
muX = 0.05
rX = 0.05
sigmaX = 0.2

# Japanese Stock 
S0 = 30000
mu = 0.05
rJ = 0.05
sigma = 0.2

# Hhedging 
sigma_hedge = 0.2
Nhedge = 24

# Option
T = 1
K = 102

# Simulation
Nrep = 1000

Xt = np.array([X0]*Nrep)
St = np.array([S0]*Nrep)
dt = T/Nhedge

initialoutlay = BlackScholes(T=T,strike=K,r=rJ,q=0,sigma=sigma)
initialoutlay = X0 * initialoutlay.valuation(spot=S0)
Vpf = np.array([initialoutlay]* Nrep)

model = BlackScholes(T=T,strike=K,r=rJ,q=0,sigma=sigma_hedge,option_type='call',greek_type='delta')
a = np.array([model.valuation(s) for s in St.tolist()])
b = Vpf - a*St

from tqdm import tqdm
for i in tqdm(range(1,Nhedge)):
    # Begining of time t 
    Xt = Xt * np.exp((muX-0.5*sigmaX**2)*dt +sigmaX*np.sqrt(dt)*np.random.normal(size=Nrep))	
    St = St * np.exp((mu-0.5*sigma**2)*dt +sigma*np.sqrt(dt)*np.random.normal(size=Nrep))
    Vpf = a * St + b * np.exp(dt*rJ)  

    # End of time t
    model = BlackScholes(
        T=T-(i-1)*dt, 
        strike=K, 
        r=rJ-sigmaX*sigma, 
        q=0, 
        sigma=sigma,
        greek_type='delta')
    a = np.array([model.valuation(s) for s in St.tolist()])
    alpha = X0 * np.exp((rJ - sigmaX * sigma - rUS)*(T-(i-1)*dt)) / Xt
    b = Vpf - a * alpha * St

ST = St*np.exp((mu-0.5*sigma**2)*dt + sigma * np.sqrt(dt)*np.random.normal(size=Nrep))
XT = Xt * np.exp((muX-0.5*sigmaX**2)*dt +sigmaX*np.sqrt(dt)*np.random.normal(size=Nrep))
alpha = X0 * np.exp((rJ - sigmaX * sigma - rUS)*(T-(i-1)*dt)) / Xt
b = Vpf - a * alpha * St
Vpf = a * ST + b * np.exp(dt*rJ) 
hedgeerror = (Vpf-alpha*np.maximum(ST-K,0))
optionpayoff = alpha*np.maximum(ST-K,0)
  
print("Initial investment =",round(initialoutlay,4))
print("Average discounted option payoff =",round(np.exp(-rJ*T)*np.mean(optionpayoff),4))
print("Average discounted portfolio value =",round(np.exp(-rJ*T)*np.mean(Vpf),4))


import matplotlib.pyplot as plt 
plt.figure(figsize=(20,10))
plt.scatter(ST,Vpf,label="Discrete hedging of a call option",color='black')
plt.ylabel("Value of hedge portfolio(T)")
plt.xlabel("S(T)")
x = np.sort(ST)
plt.plot(x,alpha*np.maximum(x-K,0),label='Payoff',color='red')
plt.title("Replication of vanila option")
plt.show()