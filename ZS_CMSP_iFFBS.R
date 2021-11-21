#ZSCMSP model from Table 1 with iFFBS sampler for the hidden states
rm(list=ls())
library(nimble)

memory.limit(size=1E10)
setwd("~/GitHub/ZS_CMSP_code")

#bring in cleaned data
load("final_data_cleaned.Rdata")

#cases_it=y_it
cases <-dnew2[,10:93]
cases <- as.matrix(cases)
plot(cases[1,])

#Rain_t 
Rain <- d$TotPrec[1:84]

#Temp_t
Temp <- d$MeanMaxTemp[1:84]

#HDI_i 
HDI <- dnew2$HDI

#N_i 
N <- dnew2$Npop

#population in thousands 
pop <- N/1000

#construct neighborhood matrix 
N_matrix <- is_connect_nei
colnames(N_matrix) <- c("CodDistrict.y",paste0("V",N_matrix$CodDistrict))
N_matrix <- N_matrix[order(N_matrix$CodDistrict.y),]
N_matrix <- N_matrix[,c("CodDistrict.y",paste0("V",1:12),paste0("V",14:161))]


#check symmetric
sum(N_matrix[1:160,2:161] != t(N_matrix[1:160,2:161]))
N_matrix <- N_matrix[1:160,2:161]
#its good 
sum(N_matrix[1:160,1:160] != t(N_matrix[1:160,1:160]))
min(rowSums(N_matrix))
#symmetric and every location has at least one neighbor

NM <- as.matrix(N_matrix)

#make a count vector 
count <- rep(NA,160)
num <- as.numeric(as.carAdjacency(N_matrix)$num)
count[1] <- 1 
for(i in 1:159){
  count[i+1]  <- count[i]+num[i]
}

#also need S_it=NA if y_it=0 and S_it=1 if y_it>0
S <- matrix(nrow=160,ncol=84)
for(i in 1:160){
  for(t in 1:84){
    if(cases[i,t]>0){
      S[i,t] <- 1
    }
  }
}

dengeeConsts <- list(N=160,T=84,adj=as.numeric(as.carAdjacency(N_matrix)$adj),
                     num=as.numeric(as.carAdjacency(N_matrix)$num),
                     pop=pop,Temp=Temp,Rain=Rain,HDI=HDI,
                     count=count,psi =cases,mpsi=mean(cases),lpsi=log(cases+1),
                     mlpsi=mean(log(cases+1)))

dengeeData <- list(y=cases,S=S) 

dengeeCode <- nimbleCode({
  #priors
  beta0 ~ dnorm(0,sd=100)
  beta1 ~dnorm(0,sd=100)
  beta2 ~ dnorm(0, sd=100)
  rho~dnorm(0,sd=100)
  for(j in 1:14){
    alpha[j]~dnorm(0,sd=100)
  }
  for(i in 1:N){
    b0[i]~dnorm(beta0,sd=sigma_b0)
    b[i]~dnorm(rho,sd=sigma_b)
  }
  sigma_b ~ dgamma(.001,.001)
  sigma_b0 ~ dgamma(.001,.001)
  #likelihood
  for(i in 1:N) {
    for(t in 2:T){
      mup[i,t] <- exp(b0[i]+beta1*(Rain[t-1]-mean(Rain[]))+
                        beta2*(Temp[t-1]-mean(Temp[])))*psi[i,t-1]+exp(b[i])
      mu[i,t] <- mup[i,t]*(S[i,t])+.00001
      y[i,t] ~ dpois(mu[i,t])
      like[i,t] <- dpois(x=y[i,t],lambda=mu[i,t])
    }
  }
  #markov chain
  for(i in 1:N) {
    S[i,1] ~ dbern(.5)
    for(t in 2:T){
      
      #calc Neighbor Sum for t=t-1 
      Snei[num[i],i,(t-1)] <- S[adj[count[i]],(t-1)]
      for(j in 2:num[i]){
        Snei[num[i]-j+1,i,(t-1)] <- Snei[num[i]-j+2,i,(t-1)] + S[adj[count[i]+j-1],(t-1)]
      }
      
      lp01[i,t] <- alpha[1]+alpha[2]*(pop[i]-mean(pop[]))+alpha[3]*(HDI[i]-mean(HDI[]))+
        alpha[4]*(Temp[t-1]-mean(Temp[]))+alpha[5]*(Rain[t-1]-mean(Rain[]))+
        alpha[6]*Snei[1,i,(t-1)]
      lp11[i,t] <- alpha[7]+alpha[8]*(pop[i]-mean(pop[]))+
        alpha[9]*(HDI[i]-mean(HDI[]))+alpha[10]*(Temp[t-1]-mean(Temp[]))+
        alpha[11]*(Rain[t-1]-mean(Rain[]))+alpha[12]*Snei[1,i,(t-1)]+
        alpha[13]*(lpsi[i,t-1]-mlpsi)
      logit(p1[i,t]) <- lp01[i,t]*(1-S[i,t-1])+lp11[i,t]*S[i,t-1]
      S[i,t] ~ dbern(p1[i,t])
    }
  }
})

inits <- list(beta0 = rnorm(n=1,mean=0,sd=1),beta1 = rnorm(n=1,mean=0,sd=.01),
              beta2 = rnorm(n=1,mean=0,sd=.01) ,
              sigma_b = .01,sigma_b0 = .01,rho=rnorm(n=1,mean=0,sd=1) ,
              S=matrix(rep(0,160*84),nrow = 160,ncol=84),
              alpha=rnorm(n=14,mean=0,sd=c(1,.01,1,.1,.01,.1,.1,.01,1,.1,.01,.1,.01,.01)),
              b=rnorm(n=160,mean=0,sd=.01),b0=rnorm(n=160,mean=0,sd=.01))

dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)

Cdengee <- compileNimble(dengeemodel)

#bring in the custom samplers
source("samplers.R")

dengeeConf <- configureMCMC(dengeemodel, print = TRUE)

#add iFFBS sampler, this might not work in new versions of nimble
for(loc in 1:160){
  
  
  dq <- !dengeemodel$isData(paste0("S[",loc,", ]"))
  if(sum(dq)>0){
    dengeeConf$removeSampler(paste0("S[",loc,", dq]"))
    dengeeConf$addSampler(target = dengeemodel$expandNodeNames(paste0("S[",loc,", dq]")),
                          type = "iZIPFFBS4")
  }
  
}

dengeeConf

#try a RW block sampler
dengeeConf$removeSampler(c("alpha[1]","alpha[6]","alpha[7]","alpha[12]"))
dengeeConf$addSampler(target = c("alpha[1]","alpha[6]"),type = "RW_block")
dengeeConf$addSampler(target = c("alpha[7]","alpha[12]"),type = "RW_block")

#sample sd on log scale
dengeeConf$removeSampler(c("sigma_b","sigma_b0"))
dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))
dengeeConf$addSampler(target=c("sigma_b0"),type="RW",control=list(log=TRUE))

print(dengeeConf)

dengeeConf$addMonitors(c("beta0","rho" ,"b","sigma_b","alpha","like","beta1","beta2","b0"))

dengeeMCMC <- buildMCMC(dengeeConf)

CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)

start_time <- Sys.time()
niter <- 40000
#when not running in parrelel, random inits must be specified here as well
samples <- runMCMC(CdengeeMCMC, niter = niter,nchains = 3,nburnin=10000
                   ,samplesAsCodaMCMC = TRUE,thin=3,
                   inits=list(beta0 = rnorm(n=1,mean=0,sd=1),beta1 = rnorm(n=1,mean=0,sd=.01),
                              beta2 = rnorm(n=1,mean=0,sd=.01) ,
                              sigma_b = .01,sigma_b0 = .01,rho=rnorm(n=1,mean=0,sd=1) ,
                              S=matrix(rep(0,160*84),nrow = 160,ncol=84),
                              alpha=rnorm(n=14,mean=0,sd=c(1,.01,1,.1,.01,.1,.1,.01,1,.1,.01,.1,.01,.01)),
                              b=rnorm(n=160,mean=0,sd=.01),b0=rnorm(n=160,mean=0,sd=.01)))
end_time <- Sys.time()

end_time - start_time

run_time <- end_time-start_time

#check convergence
library(coda)
effectiveSize(samples[,c("beta0","beta1","beta2","sigma_b","rho","alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]",
                         "alpha[7]","alpha[8]","alpha[9]","alpha[10]","alpha[11]","alpha[12]","alpha[13]","alpha[14]")])
min(effectiveSize(samples[,c("beta0","beta1","beta2","sigma_b","rho","alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]",
                             "alpha[7]","alpha[8]","alpha[9]","alpha[10]","alpha[11]","alpha[12]","alpha[13]","alpha[14]")])/as.numeric(run_time))
gelman.diag(samples[,c("beta0","beta1","beta2","sigma_b","rho","alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]",
                       "alpha[7]","alpha[8]","alpha[9]","alpha[10]","alpha[11]","alpha[12]","alpha[13]","alpha[14]")])


#can also examine trace plots
plot(samples[,c("beta0","rho","b[1]","sigma_b")])
plot(samples[,c("beta0","rho","b0[1]","sigma_b0")])
plot(samples[,c("beta1","beta2")])
plot(samples[,c("b[4]","b[143]","b[50]","b[100]")])
plot(samples[,c("alpha[1]","alpha[2]","alpha[3]","alpha[4]")])
plot(samples[,c("alpha[1]","alpha[2]")])
plot(samples[,c("alpha[5]","alpha[6]","alpha[7]","alpha[8]")])
plot(samples[,c("alpha[2]","alpha[8]")])
plot(samples[,c("alpha[9]","alpha[10]","alpha[11]","alpha[12]")])
plot(samples[,c("alpha[13]","alpha[12]")])
plot(samples[,c("alpha[6]","alpha[12]")])
plot(samples[,c("alpha[4]","alpha[10]")])
plot(samples[,c("alpha[5]","alpha[11]")])
plot(samples[,c("alpha[13]","alpha[7]","alpha[1]")])


#WAIC
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
like   <- samps[,startsWith(names(samps),"like")] 
#0s and NAs from t=1 which is saved this way for some reason (we start likelihood at t=2)
#it is simply removed below and has no impact
fbar   <- colMeans(like+.000001)
Pw     <- sum(apply(log(like+.000001),2,var),na.rm = TRUE)
WAIC   <- -2*sum(log(fbar),na.rm = TRUE)+2*Pw
remove(like)
remove(fbar)





