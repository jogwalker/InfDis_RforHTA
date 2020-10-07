#Load libraries
library(ggplot2)
library(deSolve)
library(tidyverse)

# Initial conditions

initial_values=c(SP=999999,IP=1,TIP=0,CP=0,TCP=0,RP=0,RCP=0,SX=0,IX=0,TIX=0,CX=0,TCX=0,RX=0,RCX=0,D=0) # assume one PWID is infected

# Time points

time=seq(from=1,to=100,by=1)

# Baseline parameters
parameters=c(
  gamma=0, # baseline rate of treatment 
  beta=0.5, # transmission rate per contact per year 
  tl=1/0.25, # rate of leaving treatment 
  svr=0.95, # cure rate
  muC=1/25, # death rate with cirrhosis
  muX=1/70, # death rate for exPWID
  muP=1/40, # death rate for PWID
  zeta= 1/20, # progression to cirrhosis
  tau=1/100 # average injecting duration of 10 years
)

# SIR model function 

sir_model <- function(time,state,parameters){
  with(as.list(c(state,parameters)),{
    
    # force of infection assuming only transmission in PWID
    NP=SP+IP+TIP+CP+TCP+RP+RCP
    lambda=beta*((IP+TIP+CP+TCP)/NP) 
    # B = muP*NP + muX*(SX+IX+CX+TIX+TCX+RX+RCX)
    dSP=  -lambda*SP - tau*SP - muP*SP # birth (non HCV deaths), - infection, cessation, death
    dIP= lambda*SP + (1-svr)*tl*TIP -gamma*IP  - zeta*IP - muP*IP - tau*IP# +infection, failed treatment, - treatment, -progression, -death, -cessation
    dCP= -gamma*CP -(muC+muP)*CP + zeta*IP + (1-svr)*tl*TCP -tau*CP# +progression, failed treatment, -death, - treatment, - cessation
    dTIP= gamma*IP - tl*TIP - muP*TIP - tau*TIP # treatment, - leaving treatment, - death, cessation
    dTCP= gamma*CP - tl*TCP - muC*TCP - muP*TCP -tau*TCP
    dRP= svr*tl*TIP - muP*RP - tau*RP# recovered, death, cessation
    dRCP = svr*tl*TCP - muP*RCP -muC*RCP - tau*RCP
      
    dSX= tau*SP - muX*SX
    dIX= tau*IP - muX*IX - gamma*IX + (1-svr)*tl*TIX
    dCX= tau*CP - muX*CX - muC*CX - gamma*CX + (1-svr)*tl*TCX
    dTIX= tau*TIP - muX*TIX + gamma*IX - tl*TIX
    dTCX= tau*TCP - muX*TCX - muC*TCX + gamma*CX - tl*TCX
    dRX= tau*RP - muX*RX + svr*tl*TIX
    dRCX = tau*RCP - muX*RCX - muC*RCX + svr*tl*TCX
      
    dD= muP*NP + muX*(SX+IX+CX+TIX+TCX+RX+RCX) +muC*(CP+RCP+TCP+CX+TCX+RCX)
      
    return(list(c(dSP,dIP,dCP,dTIP,dTCP,dRP,dRCP,dSX,dIX,dCX,dTIX,dTCX,dRX,dRCX,dD)))
  }
  )
}





#Solving the differential equations
output<-as.data.frame(ode(y=initial_values,func = sir_model,parms=parameters,times = time))


out_long <- output %>% pivot_longer(2:6)

# To plot number in each compartment over time
ggplot(data = out_long,          
       aes(x = time, y = value/1000000, colour = name, group = name)) +  
  geom_line() +xlab("Time (years)")+ylab("Proportion of the population")+scale_color_discrete(name="State")



