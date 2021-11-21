#endemic/epidemic model from Table 1 of the main text
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

#now calc prior mean/sd of log(r_i) (from Bauer 2018 paper cited in main text)
theta1 <- .5
theta2 <- 4

lru <- rep(NA,160)
lrsd <- rep(NA,160)
for(i in 1:160){
  lru[i] <- log(mean(cases[i,]))-log(theta1)
  lrsd[i] <- log(theta1/theta2)/(-1.64)
}

lru[lru==-Inf] <- min(lru[lru!=-Inf])



dengeeConsts <- list(N=160,T=84,
                     pop=pop,Temp=Temp,Rain=Rain,HDI=HDI,
                     count=count,psi =cases,mpsi=mean(cases),lpsi=log(cases+1),
                     mlpsi=mean(log(cases+1)),lru=lru,lrsd=lrsd)

dengeeData <- list(y=cases) 


dengeeCode <- nimbleCode({
  #priors
  beta0 ~ dnorm(0,sd=100)
  beta1 ~dnorm(0,sd=100)
  beta2 ~ dnorm(0, sd=100)
  rho~dnorm(0,sd=100)
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
      
      p[i,t] <- r[i]/(r[i]+mup[i,t]) 
      y[i,t] ~ dnegbin(p[i,t],r[i])
    }
  }
})

inits <- list(beta0 = rnorm(n=1,mean=0,sd=1),beta1 = rnorm(n=1,mean=0,sd=.01),
              beta2 = rnorm(n=1,mean=0,sd=.01) ,
              sigma_b = .01,prec_b0 = 1,rho=rnorm(n=1,mean=0,sd=1),
              b=rnorm(n=160,mean=0,sd=.01),b0=rnorm(n=160,mean=0,sd=.01),
              lr=runif(n=160,min=-10,max=10))

dengeemodel <- nimbleModel(dengeeCode, dengeeConsts, dengeeData, inits)

Cdengee <- compileNimble(dengeemodel)

dengeeConf <- configureMCMC(dengeemodel, print = TRUE)

#sample sd on log scale
dengeeConf$removeSampler(c("sigma_b"))
dengeeConf$addSampler(target=c("sigma_b"),type="RW",control=list(log=TRUE))

print(dengeeConf)

dengeeConf$addMonitors(c("beta0","rho" ,"b","sigma_b","beta1","beta2","r","b0","sigma_b0"))

dengeeMCMC <- buildMCMC(dengeeConf)

CdengeeMCMC <- compileNimble(dengeeMCMC, project = dengeemodel,resetFunctions = TRUE)

initsFunction <- function() list(beta0 = rnorm(n=1,mean=0,sd=1),beta1 = rnorm(n=1,mean=0,sd=.01),
                                 beta2 = rnorm(n=1,mean=0,sd=.01) ,
                                 sigma_b = .01,prec_b0 = 1,rho=rnorm(n=1,mean=0,sd=1),
                                 b=rnorm(n=160,mean=0,sd=.01),b0=rnorm(n=160,mean=0,sd=.01),
                                 lr=runif(n=160,min=-10,max=10))

start_time <- Sys.time()
samples <- runMCMC(CdengeeMCMC, niter = 80000,nchains = 3,nburnin=30000
                   ,samplesAsCodaMCMC = TRUE,thin=5,inits = initsFunction)

end_time <- Sys.time()

end_time - start_time

run_time <- end_time-start_time

##check convergence
library(coda)
effectiveSize(samples[,c("r[1]","beta0","beta1","beta2","sigma_b","sigma_b0","rho")])
min(effectiveSize(samples[,c("r[1]","beta0","beta1","beta2","sigma_b","sigma_b0","rho")]))
gelman.diag(samples[,c("r[1]","beta0","beta1","beta2","sigma_b","sigma_b0","rho")])

#randomly check location specific parameters
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

#look at trace plots
plot(samples[,c("beta0","rho","b[1]")])
plot(samples[,c("sigma_b0","sigma_b")])
plot(samples[,c("beta1","beta2","r[1]")])
plot(samples[,c("b[4]","b[143]","b[50]","b[100]")])
plot(samples[,c("b0[4]","b0[143]","b0[50]","b0[100]")])
plot(samples[,c("r[35]","r[120]","r[127]","r[135]")])

##WAIC
samps <- data.frame(rbind(samples[[1]],samples[[2]],samples[[3]]))
lppd <- 0
pwaic <- 0
for(i in 1:160){
  print(i)
  for(t in 2:84){
    
    mup <-  as.numeric(unlist(exp(samps[paste0("b0.",i,".")]+samps[paste0("beta1")]*(Rain[t-1]-mean(Rain[]))+
                                                                      samps[paste0("beta2")]*(Temp[t-1]-mean(Temp[])))*cases[i,t-1]+exp(samps[paste0("b.",i,".")])))
    
    pnegbin<- as.numeric(unlist(samps[paste0(paste0("r.",i,"."))]/(samps[paste0(paste0("r.",i,"."))]+mup)))
    
    lppd <- lppd + log(mean(dnbinom(x = cases[i,t], size=as.numeric(unlist(samps[paste0(paste0("r.",i,"."))])), prob=pnegbin)))
    pwaic <- pwaic + var(log(dnbinom(x = cases[i,t], size=as.numeric(unlist(samps[paste0(paste0("r.",i,"."))])), prob=pnegbin)))
  }
  
}

waic <- -2*(lppd-pwaic)

#WAIC = 69,939