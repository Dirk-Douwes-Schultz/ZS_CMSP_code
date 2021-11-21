## best fitting (bolded in Table 1) ZS-CMSNB model from Table 1 of the main text with 
#binary samplers for the hidden states
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

#now calc prior mean/sd of log(r_i), see notes (from Bauer 2018 paper)
theta1 <- .5
theta2 <- 4

lru <- rep(NA,160)
lrsd <- rep(NA,160)
for(i in 1:160){
  lru[i] <- log(mean(cases[i,]))-log(theta1)
  lrsd[i] <- log(theta1/theta2)/(-1.64)
}

lru[lru==-Inf] <- min(lru[lru!=-Inf])

#sum cases in neighboring areas
log_sum_nei_cases <- matrix(nrow=160,ncol=84)
for(i in 1:160){
  for(t in 2:84){
    nei <- which(N_matrix[i,]==1)
    log_sum_nei_cases[i,t] <- 0
    for(j in nei){
      log_sum_nei_cases[i,t] <- log_sum_nei_cases[i,t]+log(cases[j,t-1]+1)
    }
  }
}
hist(log_sum_nei_cases)
mean_log_sum_nei_cases <- mean(log_sum_nei_cases,na.rm = TRUE)

#sum of prevalence in neighboring areas
log_sum_nei_prev <- matrix(nrow=160,ncol=84)
for(i in 1:160){
  for(t in 2:84){
    nei <- which(N_matrix[i,]==1)
    log_sum_nei_prev[i,t] <- 0
    for(j in nei){
      log_sum_nei_prev[i,t] <- log_sum_nei_prev[i,t]+log((cases[j,t-1]/pop[j])+1)
    }
  }
}
hist(log_sum_nei_prev)
mean_log_sum_nei_prev <- mean(log_sum_nei_prev,na.rm = TRUE)

#need to calculate log of prevalence
prev <- matrix(nrow = 160,ncol=84)
for(i in 1:160){
  for(t in 1:84){
    prev[i,t] <- cases[i,t]/pop[i]
  }
}

hist(log(prev+1))


dengeeConsts <- list(N=160,T=84,adj=as.numeric(as.carAdjacency(N_matrix)$adj),
                     num=as.numeric(as.carAdjacency(N_matrix)$num),
                     pop=pop,Temp=Temp,Rain=Rain,HDI=HDI,
                     count=count,psi =cases,mpsi=mean(cases),lpsi=log(cases+1),
                     mlpsi=mean(log(cases+1)),lprev = log(prev+1),mlprev=mean(log(prev+1)),
                     lru=lru,lrsd=lrsd,
                     log_sum_nei_cases=log_sum_nei_cases,
                     mean_log_sum_nei_cases=mean_log_sum_nei_cases,
                     log_sum_nei_prev=log_sum_nei_prev,
                     mean_log_sum_nei_prev=mean_log_sum_nei_prev)

dengeeData <- list(y=cases,S=S) 

dengeeCode <- nimbleCode({
  #priors
  beta0 ~ dnorm(0,sd=100)
  beta1 ~dnorm(0,sd=100)
  beta2 ~ dnorm(0, sd=100)
  rho~dnorm(0,sd=100)
  for(j in 1:15){
    alpha[j]~dnorm(0,sd=100)
  }
  for(i in 1:N){
    b0[i]~dnorm(beta0,prec_b0)
    b[i]~dnorm(rho,sd=sigma_b)
  }
  sigma_b ~ dunif(0,10)
  prec_b0 ~ dgamma(.1,.1)
  sigma_b0 <- 1/sqrt(prec_b0)
  for(i in 1:160){
    r[i] <- exp(lr[i])
    lr[i] ~ dnorm(lru[i],sd=lrsd[i])
  }
  #likelihood
  for(i in 1:N) {
    for(t in 2:T){
      mup[i,t] <- exp(b0[i]+beta1*(Rain[t-1]-mean(Rain[]))+
                        beta2*(Temp[t-1]-mean(Temp[])))*psi[i,t-1]+exp(b[i])
      
      p[i,t] <- r[i]/(r[i]+(S[i,t])*mup[i,t]) - 1e-10*(1-S[i,t])
      y[i,t] ~ dnegbin(p[i,t],r[i])
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
      
      #calc Neighbor Sum pop for t=t-1 
      Snei_pop[num[i],i,(t-1)] <- S[adj[count[i]],(t-1)]*log(pop[adj[count[i]]]*pop[i])
      for(j in 2:num[i]){
        Snei_pop[num[i]-j+1,i,(t-1)] <- Snei_pop[num[i]-j+2,i,(t-1)] + S[adj[count[i]+j-1],(t-1)]*log(pop[adj[count[i]+j-1]]*pop[i])
      }
      
      lp01[i,t] <- alpha[1]+alpha[2]*(pop[i]-mean(pop[]))+alpha[3]*(HDI[i]-mean(HDI[]))+
        alpha[4]*(Temp[t-1]-mean(Temp[]))+alpha[5]*(Rain[t-1]-mean(Rain[]))+
        alpha[6]*Snei[1,i,(t-1)]+alpha[14]*(log_sum_nei_prev[i,t]-mean_log_sum_nei_prev)+
        alpha[15]*Snei_pop[1,i,(t-1)]
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
              sigma_b = .01,prec_b0 = 1,rho=rnorm(n=1,mean=0,sd=1) ,
              S=matrix(rep(0,160*84),nrow = 160,ncol=84),
              alpha=rnorm(n=15,mean=0,sd=c(1,.01,1,.1,.01,.1,.1,.01,1,.1,.01,.1,.01,.01,.01)),
              b=rnorm(n=160,mean=0,sd=.01),b0=rnorm(n=160,mean=0,sd=.01),
              lr=runif(n=160,min=-10,max=10))

dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)

Cdengee <- compileNimble(dengeemodel)

dengeeConf <- configureMCMC(dengeemodel, print = TRUE)

#try a RW block sampler
dengeeConf$removeSampler(c("alpha[1]","alpha[6]","alpha[7]","alpha[12]","alpha[2]","alpha[14]","alpha[15]"))
dengeeConf$addSampler(target = c("alpha[1]","alpha[6]","alpha[2]","alpha[14]","alpha[15]"),type = "AF_slice")
dengeeConf$addSampler(target = c("alpha[7]","alpha[12]"),type = "AF_slice")

#sample sd on log scale
dengeeConf$removeSampler(c("sigma_b"))
dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))

print(dengeeConf)

dengeeConf$addMonitors(c("beta0","rho" ,"b","sigma_b","alpha","beta1","beta2","r","b0","sigma_b0"))

dengeeMCMC <- buildMCMC(dengeeConf)

CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)

initsFunction <- function() list(beta0 = rnorm(n=1,mean=0,sd=1),beta1 = rnorm(n=1,mean=0,sd=.01),
                                 beta2 = rnorm(n=1,mean=0,sd=.01) ,
                                 sigma_b = .01,prec_b0 = 1,rho=rnorm(n=1,mean=0,sd=1) ,
                                 S=matrix(rep(0,160*84),nrow = 160,ncol=84),
                                 alpha=rnorm(n=15,mean=0,sd=c(1,.01,1,.1,.01,.1,.1,.01,1,.1,.01,.1,.01,.01,.01)),
                                 b=rnorm(n=160,mean=0,sd=.01),b0=rnorm(n=160,mean=0,sd=.01),
                                 lr=runif(n=160,min=-10,max=10))




start_time <- Sys.time()
samples <- runMCMC(CdengeeMCMC, niter = 80000,nchains = 3,nburnin=30000
                   ,samplesAsCodaMCMC = TRUE,thin=5,inits = initsFunction)




end_time <- Sys.time()
end_time - start_time

run_time <- end_time-start_time


#check convergence
library(coda)
effectiveSize(samples[,c("r[1]","beta0","beta1","beta2","sigma_b","sigma_b0","rho","alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]",
                         "alpha[7]","alpha[8]","alpha[9]","alpha[10]","alpha[11]","alpha[12]","alpha[13]","alpha[14]","alpha[15]")])
min(effectiveSize(samples[,c("r[1]","beta0","beta1","beta2","sigma_b","rho","alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]",
                             "alpha[7]","alpha[8]","alpha[9]","alpha[10]","alpha[11]","alpha[12]","alpha[13]","alpha[14]","alpha[15]")]))
gelman.diag(samples[,c("r[1]","beta0","beta1","beta2","sigma_b","sigma_b0","rho","alpha[1]","alpha[2]","alpha[3]","alpha[4]","alpha[5]","alpha[6]",
                       "alpha[7]","alpha[8]","alpha[9]","alpha[10]","alpha[11]","alpha[12]","alpha[13]","alpha[14]","alpha[15]")])


#test a few location specific effects
gelman.diag(samples[,c("b[1]","b[2]","b[3]","b[4]",
                       "b[21]","b[39]","b[112]","b[55]",
                       "b[23]","b[29]","b[15]","b[74]",
                       "b[41]","b[123]","b[16]","b[59]")])

gelman.diag(samples[,c("b0[1]","b0[2]","b0[3]","b0[4]",
                       "b0[21]","b0[39]","b0[112]","b0[55]",
                       "b0[23]","b0[29]","b0[15]","b0[74]",
                       "b0[41]","b0[123]","b0[16]","b0[59]")])

gelman.diag(samples[,c("r[1]","r[2]","r[3]","r[4]",
                       "r[21]","r[39]","r[112]","r[55]",
                       "r[23]","r[29]","r[15]","r[74]",
                       "r[41]","r[123]","r[16]","r[59]",
                       "r[35]","r[120]","r[127]","r[135]")])

#can also examine trace plots
plot(samples[,c("beta0","rho")])
plot(samples[,c("sigma_b0","sigma_b")])
plot(samples[,c("beta1","beta2")])
plot(samples[,c("b[4]","b[143]","b[50]","b[100]")])
plot(samples[,c("b0[4]","b0[143]","b0[50]","b0[100]")])
plot(samples[,c("r[35]","r[120]","r[127]","r[135]")])
plot(samples[,c("alpha[1]","alpha[2]","alpha[3]","alpha[4]")])
plot(samples[,c("alpha[5]","alpha[6]","alpha[7]","alpha[8]")])
plot(samples[,c("alpha[9]","alpha[10]","alpha[11]","alpha[12]")])
plot(samples[,c("alpha[13]","alpha[14]","alpha[15]")])

#WAIC
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
lppd <- 0
pwaic <- 0
for(i in 1:160){
  print(i)
  for(t in 2:84){
    
    mup <-  as.numeric(unlist((exp(samps[paste0("b0.",i,".")]+samps[paste0("beta1")]*(Rain[t-1]-mean(Rain[]))+
                                     samps[paste0("beta2")]*(Temp[t-1]-mean(Temp[])))*cases[i,t-1]+exp(samps[paste0("b.",i,".")]))))
    
    pnegbin<- as.numeric(unlist(samps[paste0(paste0("r.",i,"."))]/(samps[paste0(paste0("r.",i,"."))]+samps[paste0("S.",i,"..",t,".")]*mup)-.000001*(1-samps[paste0("S.",i,"..",t,".")])))
    
    lppd <- lppd + log(mean(dnbinom(x = cases[i,t], size=as.numeric(unlist(samps[paste0(paste0("r.",i,"."))])), prob=pnegbin)))
    pwaic <- pwaic + var(log(dnbinom(x = cases[i,t], size=as.numeric(unlist(samps[paste0(paste0("r.",i,"."))])), prob=pnegbin)))
  }
  
}

waic <- -2*(lppd-pwaic)

#WAIC = 67,632

