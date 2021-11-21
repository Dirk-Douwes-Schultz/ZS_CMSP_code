#samplers are here

#functions for debugging
write_output <- function(x){
  
  ttt <- as.numeric(Sys.time())
  writeLines(paste0(x), paste0("~/GitHub/ZS_CMSP_code/error_logs/",ttt,"_output.txt"))
}

#this will break sampler which is needed
write_output_broken <- function(x){
  
  ttt <- Sys.time()
  writeLines(paste0(x), paste0("~/GitHub/ZS_CMSP_code/error_logs/",ttt,"_output.txt"))
}



R_write_output <- nimbleRcall(function(x = double(1)){}, Rfun = 'write_output',
                              returnType = void())

R_write_output_broken <- nimbleRcall(function(x = double(1)){}, Rfun = 'write_output_broken',
                                     returnType = void())

#iFFBS sampler
iZIPFFBS4 <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    #setup
    nnames <- model$expandNodeNames(target)
    times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames))
    loc <- as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames[[1]]))
    numnodes <- length(nnames)
    startbs <- max(times)
    startfs <- min(times)
    Mt <- length(model$y[loc, ])
    startmidfzt <- ifelse(startfs==1,2,1)
    endmidfzt <- ifelse(startbs==Mt,numnodes-1,numnodes)
    startizt <- ifelse(startbs==Mt,2,1)
    #all dependencies together, not same as "dependencies" because that deletes duplicates
    calcNodes <- model$getDependencies(target)
    #grab the dependencies of each target node
    dependencies <- NULL
    start_depend <- rep(NA,numnodes)
    start_depend[1] <- 1
    end_depend <- rep(NA,numnodes)
    index <- 1
    for (n in nnames){
      d <- model$getDependencies(n)
      dependencies <- c(dependencies,d)
      end_depend[index] <- length(d)+start_depend[index]-1
      start_depend[index+1] <- end_depend[index]+1
      index <- index+1
    }
    #need forward dependencies, this may have to be modified if you change model specification
    f_dependencies <- NULL
    f_start_depend <- rep(NA,numnodes)
    f_start_depend[1] <- 1
    f_end_depend <- rep(NA,numnodes)
    index <- 1
    for(n in nnames){
      d <- model$getDependencies(n)
      d <- d[grepl(paste0(".*(S)\\[(?!",loc,").*,.*\\].*"),d,perl=TRUE)]
      f_dependencies <- c(f_dependencies,d)
      f_end_depend[index] <- length(d)+f_start_depend[index]-1
      f_start_depend[index+1] <- f_end_depend[index]+1
      index <- index+1
    }
    #be careful because S[loc, Mt] has no forward dependencies
    #Now grab the neighbors
    neighbors <- NULL 
    start_nei <- rep(NA,numnodes)
    start_nei[1] <- 1
    end_nei <- rep(NA,numnodes)
    index <- 1
    for(zt in 1:numnodes){
      d <- unique(as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", 
                                  f_dependencies[f_start_depend[zt]:f_end_depend[zt]])))
      neighbors <- c(neighbors,d)
      end_nei[index] <- length(d)+start_nei[index]-1
      start_nei[index+1] <- end_nei[index]+1
      index <- index+1 
    }
    #again watch out for Mt dont use this for final timepoint as no forward neighbors
    # all matrrixs have collumn order; 0 state, 1 state
    #adjusted filtered probs 
    q <- matrix(nrow=numnodes,ncol=2)
    q[1,1:2] <- c(-.99,-.99)
    #log forward adjustment
    lf <- matrix(nrow=numnodes,ncol=2)
    lf[1,1:2] <- c(-.99,-.99)
    #transition matrix
    tm <- matrix(nrow=2,ncol=2)
    tm[1,1:2] <- c(-.99,-.99)
    
  },
  
  run = function() {
    
    #now start filter
    #if startfs is 1 then have to run inital update
    if(startfs==1){
      lf[1,model$S[loc,1]+1] <<- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      oS <- 1-model$S[loc,1]
      #if oS is 1 then going from 0->1 so need to add spread parameter
      if(oS==1){
        ff <- 1
      }else{
        ff <- -1
      }
      loS <- 0
      for(j in start_nei[1]:end_nei[1]){
        c_nei <- neighbors[j]
        #if previous neighbor value is 1 then need to consider 0 state spread parameter
        if(model$S[c_nei,1]==1){
          add <- ff*model$alpha[12]
        }else{
          add <- ff*(model$alpha[6]+model$alpha[15]*log(model$pop[loc]*model$pop[c_nei]))
        }
        newp <- expit(logit(model$p1[c_nei,2])+add)
        if(model$S[c_nei,2]==1){
          loS <- loS +log(newp)
        }else{
          loS <-loS +log(1-newp)
        }
        #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
        if(is.na(loS)|is.nan(loS)){
          R_write_output(x=c(loS,loc,lf[1,2],oS,c_nei,add,model$p1[c_nei,2],newp,
                             model$S[c_nei,2],(1-model$S[c_nei,2])*log(1-newp),
                             model$S[c_nei,2]*log(newp)))
          R_write_output_broken(x=c(loS))
        }
      }
      
      lf[1,oS+1] <<- loS
      
      q2 <- 1/(1+exp(lf[1,1]-lf[1,2]))
      q[1,1:2] <<- c(1-q2,q2)
      
    }
    
    #now normal filter
    if(startmidfzt<=endmidfzt){
      for(zt in startmidfzt:endmidfzt){
        #zt <- 2
        ct <- times[zt]
        #check if previous time is data or not
        if((zt-1) != 0){
          if(times[zt-1]==ct-1){
            q_tm1 <- q[zt-1,1:2]
          }else{
            q_tm1 <- c(0,1)
          }
        }else{
          q_tm1 <- c(0,1)
        }
        
        tm[1,1:2] <<- c(1-expit(model$lp01[loc,ct]),expit(model$lp01[loc,ct]))
        tm[2,1:2] <<- c(1-expit(model$lp11[loc,ct]),expit(model$lp11[loc,ct]))
        
        p <- t(tm) %*% asCol(q_tm1[1:2])
        
        lf[zt,model$S[loc,ct]+1] <<- model$getLogProb(f_dependencies[f_start_depend[zt]:f_end_depend[zt]])
        
        oS <- 1-model$S[loc,ct]
        #if oS is 1 then going from 0->1 so need to add spread parameter
        if(oS==1){
          ff <- 1
        }else{
          ff <- -1
        }
        loS <- 0
        
        #calculate forward cross product of dependent chains
        for(j in start_nei[zt]:end_nei[zt]){
          c_nei <- neighbors[j]
          #if previous neighbor value is 1 then need to consider persistence state spread parameter
          if(model$S[c_nei,ct]==1){
            add <- ff*model$alpha[12]
          }else{
            add <- ff*(model$alpha[6]+model$alpha[15]*log(model$pop[loc]*model$pop[c_nei]))
            #if(ct==78 & loc==3){
              #R_write_output(x=c(c_nei,model$pop[loc],model$pop[c_nei],add,model$alpha[6],model$alpha[15]))
              #R_write_output_broken(x=c(loS))
            #}
          }
          newp <- expit(logit(model$p1[c_nei,ct+1])+add)
          if(model$S[c_nei,ct+1]==1){
            loS <-loS+log(newp)
          }else{
            loS <- loS+log(1-newp)
          }
          # loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
        }
        
        lf[zt,oS+1] <<- loS
        
        #calc filtered prob
        lp0 <- model$r[loc]*(log(model$r[loc])-log(model$r[loc]+model$mup[loc,ct]))
        #lp0 <- -model$mup[loc,ct]
        q2 <- 1/(1+exp(lf[zt,1]+log(p[1,1])-lf[zt,2]-log(p[2,1])-lp0))
        q[zt,1:2] <<- c(1-q2,q2)
        
      }
    }
    
    
    #now final filter 
    if(startbs==Mt){
      if(times[numnodes-1]==Mt-1){
        q_tm1 <- q[numnodes-1,1:2]
      }else{
        q_tm1 <- c(0,1)
      }
      tm[1,1:2] <<- c(1-expit(model$lp01[loc,Mt]),expit(model$lp01[loc,Mt]))
      tm[2,1:2] <<- c(1-expit(model$lp11[loc,Mt]),expit(model$lp11[loc,Mt]))
      p <- t(tm) %*% asCol(q_tm1[1:2])
      lp0 <- model$r[loc]*(log(model$r[loc])-log(model$r[loc]+model$mup[loc,Mt]))
      #lp0 <- -model$mup[loc,Mt]
      q2 <- 1/(1+exp(log(p[1,1])-log(p[2,1])-lp0))
      q[numnodes,1:2] <<- c(1-q2,q2)
    }
    
    #backwards sampling time.
    
    #start with Mt 
    if(startbs==Mt){
      prev <- model$S[loc,Mt]
      model$S[loc,Mt] <<- rbinom(n=1, size=1, prob=q[numnodes,2])
      if(model$S[loc,Mt]!=prev){
        model$calculate(nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]])
        #copy(from = model, to = mvSaved, row = 1, 
        #nodes = dependencies[start_depend[numnodes]:end_depend[numnodes]], logProb = TRUE)
      }
    }
    
    for(izt in startizt:(numnodes)){
      zt <- numnodes-izt+1
      ct <- times[zt]
      prev <- model$S[loc,ct]
      tm[1,1:2] <<- c(1-expit(model$lp01[loc,ct+1]),expit(model$lp01[loc,ct+1]))
      tm[2,1:2] <<- c(1-expit(model$lp11[loc,ct+1]),expit(model$lp11[loc,ct+1]))
      trans <- tm[1:2,model$S[loc,ct+1]+1]
      #    b2 <- trans[2]*q[zt,2]/(trans[2]*q[zt,2]+trans[1]*q[zt,1])
      #    b2 <- 1/(1+exp(log(trans[2])+log(q[zt,2])-log(trans[1])-log(q[zt,1])))
      b2 <- 1/(1+exp(log(trans[1])+log(q[zt,1])-log(trans[2])-log(q[zt,2])))
      if(is.na(b2)|is.nan(b2)){
        R_write_output(x=c(b2,trans[1],q[zt,1],trans[2],q[zt,2],zt,loc,startfs,lf[1,1],lf[1,2],model$S[loc,1]))
        R_write_output_broken(x=c(b2))
        # b2 <- .5
      }
      model$S[loc,ct] <<- rbinom(n=1, size=1,prob=b2)
      if( model$S[loc,ct]!=prev){
        model$calculate(nodes = dependencies[start_depend[zt]:end_depend[zt]])
        #copy(from = model, to = mvSaved, row = 1, 
        #nodes = dependencies[start_depend[zt]:end_depend[zt]], logProb = TRUE)
      }
    }
    
    copy(from = model, to = mvSaved, row = 1, 
         nodes = calcNodes, logProb = TRUE)
    
  },
  
  methods = list(   reset = function () {}   )
  
)


#bFFBS2 sampler
bZIPFFBS2 <- nimbleFunction(
  
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    
    #setup
    nnames <- model$expandNodeNames(target)
    times <- as.numeric(gsub(".*\\[.*,(.*)\\].*", "\\1", nnames))
    times <- unique(sort(times))
    locs <- unique(as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", nnames)))
    loc1 <- locs[1]
    loc2 <- locs[2]
    numnodes <- length(nnames)
    numtimes <- length(times)
    startbs <- max(times)
    startfs <- min(times)
    Mt <- length(model$y[loc1, ])
    startmidfzt <- ifelse(startfs==1,2,1)
    endmidfzt <- ifelse(startbs==Mt,numtimes-1,numtimes)
    startizt <- ifelse(startbs==Mt,2,1)
    calcNodes <- model$getDependencies(target)
    #dependencie matrix 
    dependencies1 <- NULL
    dependencies2 <- NULL
    start_depend1 <- rep(NA,numtimes)
    start_depend1[1] <- 1
    end_depend1 <- rep(NA,numtimes)
    start_depend2 <- rep(NA,numtimes)
    start_depend2[1] <- 1
    end_depend2 <- rep(NA,numtimes)
    index <- 1
    for (ct in times){
      d1 <- model$getDependencies(paste0("S[",loc1," ,",ct,"]"))
      d2 <- model$getDependencies(paste0("S[",loc2," ,",ct,"]"))
      dependencies1 <- c(dependencies1,d1)
      dependencies2 <- c(dependencies2,d2)
      end_depend1[index] <- length(d1)+start_depend1[index]-1
      start_depend1[index+1] <- end_depend1[index]+1
      end_depend2[index] <- length(d2)+start_depend2[index]-1
      start_depend2[index+1] <- end_depend2[index]+1
      index <- index+1
    }
    # it has one extra length dont worry
    
    #now need to calculate forward nodes
    #first calculate the neighbors
    d1 <- model$getDependencies(paste0("S[",loc1," ,",times[1],"]"))
    nei1 <- d1[grepl(paste0(".*(S)\\[(?!",loc1,",).*\\].*"),d1,perl=TRUE)]
    nei1 <- unique(as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", 
                                   nei1)))
    
    d2 <- model$getDependencies(paste0("S[",loc2," ,",times[1],"]"))
    nei2 <- d2[grepl(paste0(".*(S)\\[(?!",loc2,",).*\\].*"),d2,perl=TRUE)]
    nei2 <- unique(as.numeric(gsub(".*\\[(.*),.*\\].*", "\\1", 
                                   nei2)))
    
    #are they neighbors? 
    neilocs <- ifelse(loc2 %in% nei1,1,0)
    
    nei_combine <- unique(c(nei1,nei2))
    nei_combine <- sort(nei_combine[!nei_combine %in% c(loc1,loc2)])
    
    nei1_only <- nei1[(!nei1 %in% nei2) &(!nei1 %in% c(loc1,loc2))]
    nei2_only <- nei2[(!nei2 %in% nei1) &(!nei2 %in% c(loc1,loc2))]
    nei_both <- nei_combine[(!nei_combine %in% c(nei1_only,nei2_only))&(!nei_combine %in% c(loc1,loc2))]
    
    lnei1_only <- length(nei1_only)
    lnei2_only <- length(nei2_only)
    lnei_both <- length(nei_both)
    
    #add 1 on to it
    if(lnei1_only==1){
      nei1_only <- c(nei1_only,1)
    }
    if(lnei2_only==1){
      nei2_only <- c(nei2_only,1)
    }
    if(lnei_both==1){
      nei_both <- c(nei_both,1)
    }
    
    #need f_depencies 
    f_dependencies <- NULL
    f_start_depend <- rep(NA,numtimes)
    f_start_depend[1] <- 1
    f_end_depend <- rep(NA,numtimes)
    index <- 1
    for(ct in times){
      d <- model$getDependencies(c(paste0("S[",loc1,", ",ct,"]"),paste0("S[",loc2,", ",ct,"]")))
      d <- d[grepl(paste0(".*(S)\\[(?!(",loc1,",|",loc2,",)).*\\].*"),d,perl=TRUE)]
      f_dependencies <- c(f_dependencies,d)
      f_end_depend[index] <- length(d)+f_start_depend[index]-1
      f_start_depend[index+1] <- f_end_depend[index]+1
      index <- index+1
    }
    
    #filter probs 
    q <- matrix(nrow=numtimes,ncol = 4)
    q[1,1:4] <- c(-.99,-.99,-.99,-.99)
    #blocked transition matrix
    tmB <- matrix(nrow=4,ncol=4)
    tmB[1,1:4] <- c(-.99,-.99,-.99,-.99)
    #log forward adjustment
    lf <- matrix(nrow=numtimes,ncol=4)
    lf[1,1:4] <- c(-.99,-.99,-.99,-.99)
  },
  
  run = function() {
    
    #now start filter 
    if(startfs==1){
      s1 <- model$S[loc1, 1]
      s2 <- model$S[loc2, 1]
      
      #calc (0,0) 
      if(s1==0 & s2==0){
        lf[1,1] <<- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      }else{
        if(model$y[loc1, 1]==0 & model$y[loc2, 1]==0){
          ff1 <- 0-s1
          ff2 <- 0-s2
          loS <- 0
          if(lnei1_only>0){
            for(j in 1:lnei1_only){
              c_nei <- nei1_only[j]
              if(model$S[c_nei,1]==1){
                add <- ff1*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          if(lnei2_only>0){
            for(j in 1:lnei2_only){
              c_nei <- nei2_only[j]
              if(model$S[c_nei,1]==1){
                add <- ff2*model$alpha[12]
              }else{
                add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          if(lnei_both>0){
            for(j in 1:lnei_both){
              c_nei <- nei_both[j]
              if(model$S[c_nei,1]==1){
                add <- (ff1+ff2)*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          lf[1,1] <<- loS
        }else{
          lf[1,1] <<- -Inf
        }
      }
      
      #calc (1,0) 
      if(s1==1 & s2==0){
        lf[1,2] <<- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      }else{
        if(model$y[loc2, 1]==0){
          ff1 <- 1-s1
          ff2 <- 0-s2
          loS <- 0
          if(lnei1_only>0){
            for(j in 1:lnei1_only){
              c_nei <- nei1_only[j]
              if(model$S[c_nei,1]==1){
                add <- ff1*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          if(lnei2_only>0){
            for(j in 1:lnei2_only){
              c_nei <- nei2_only[j]
              if(model$S[c_nei,1]==1){
                add <- ff2*model$alpha[12]
              }else{
                add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          if(lnei_both>0){
            for(j in 1:lnei_both){
              c_nei <- nei_both[j]
              if(model$S[c_nei,1]==1){
                add <- (ff1+ff2)*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          lf[1,2] <<- loS
        }else{
          lf[1,2] <<- -Inf
        }
      }
      
      #calc (0,1) 
      if(s1==0 & s2==1){
        lf[1,3] <<- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      }else{
        if(model$y[loc1, 1]==0){
          ff1 <- 0-s1
          ff2 <- 1-s2
          loS <- 0
          if(lnei1_only>0){
            for(j in 1:lnei1_only){
              c_nei <- nei1_only[j]
              if(model$S[c_nei,1]==1){
                add <- ff1*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          if(lnei2_only>0){
            for(j in 1:lnei2_only){
              c_nei <- nei2_only[j]
              if(model$S[c_nei,1]==1){
                add <- ff2*model$alpha[12]
              }else{
                add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          if(lnei_both>0){
            for(j in 1:lnei_both){
              c_nei <- nei_both[j]
              if(model$S[c_nei,1]==1){
                add <- (ff1+ff2)*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,2])+add)
              if(model$S[c_nei,2]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
            }
          }
          lf[1,3] <<- loS
        }else{
          lf[1,3] <<- -Inf
        }
      }
      
      #calc (1,1) 
      if(s1==1 & s2==1){
        lf[1,4] <<- model$getLogProb(f_dependencies[f_start_depend[1]:f_end_depend[1]])
      }else{
        
        ff1 <- 1-s1
        ff2 <- 1-s2
        loS <- 0
        if(lnei1_only>0){
          for(j in 1:lnei1_only){
            c_nei <- nei1_only[j]
            if(model$S[c_nei,1]==1){
              add <- ff1*model$alpha[12]
            }else{
              add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
            }
            newp <- expit(logit(model$p1[c_nei,2])+add)
            if(model$S[c_nei,2]==1){
              loS <- loS+log(newp)
            }else{
              loS <- loS+log(1-newp)
            }
            #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
          }
        }
        if(lnei2_only>0){
          for(j in 1:lnei2_only){
            c_nei <- nei2_only[j]
            if(model$S[c_nei,1]==1){
              add <- ff2*model$alpha[12]
            }else{
              add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
            }
            newp <- expit(logit(model$p1[c_nei,2])+add)
            if(model$S[c_nei,2]==1){
              loS <- loS+log(newp)
            }else{
              loS <- loS+log(1-newp)
            }
            #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
          }
        }
        if(lnei_both>0){
          for(j in 1:lnei_both){
            c_nei <- nei_both[j]
            if(model$S[c_nei,1]==1){
              add <- (ff1+ff2)*model$alpha[12]
            }else{
              add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
            }
            newp <- expit(logit(model$p1[c_nei,2])+add)
            if(model$S[c_nei,2]==1){
              loS <- loS+log(newp)
            }else{
              loS <- loS+log(1-newp)
            }
            #loS <- loS+model$S[c_nei,2]*log(newp)+(1-model$S[c_nei,2])*log(1-newp)
          }
        }
        lf[1,4] <<- loS
      }
      
      nl <- lf[1,1:4]-max(lf[1,1:4])
      q[1,1:4] <<- exp(nl)/sum(exp(nl))
      
    }
    
    ##now between filter
    if(startmidfzt<=endmidfzt){
      for(zt in startmidfzt:endmidfzt){
        ct <- times[zt]
        if((zt-1) != 0){
          if(times[zt-1]==ct-1){
            q_tm1 <- q[zt-1,1:4]
          }else{
            q_tm1 <- c(0,0,0,1)
          }
        }else{
          q_tm1 <- c(0,0,0,1)
        }
        
        #construct transition matrix
        s1 <- model$S[loc1, ct]
        s2 <- model$S[loc2, ct]
        s1m <- model$S[loc1, ct-1]
        s2m <- model$S[loc2, ct-1]
        alpha6_new <- model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[loc2])
        #p(loc)(value)(condloc1loc2)
        if(neilocs==1){
          #(0,0)
          p1100 <- expit(model$lp01[loc1,ct]+(0-s2m)*alpha6_new)
          p2100 <- expit(model$lp01[loc2,ct]+(0-s1m)*alpha6_new)
          #(1,0)
          p1110 <- expit(model$lp11[loc1,ct]+(0-s2m)*model$alpha[12])
          p2110 <- expit(model$lp01[loc2,ct]+(1-s1m)*alpha6_new)
          #(0,1)
          p1101 <- expit(model$lp01[loc1,ct]+(1-s2m)*alpha6_new)
          p2101 <- expit(model$lp11[loc2,ct]+(0-s1m)*model$alpha[12])
          #(1,1)
          p1111 <- expit(model$lp11[loc1,ct]+(1-s2m)*model$alpha[12])
          p2111 <- expit(model$lp11[loc2,ct]+(1-s1m)*model$alpha[12])
          
        }else{
          p1100 <- expit(model$lp01[loc1,ct])
          p2100 <- expit(model$lp01[loc2,ct])
          p1110 <- expit(model$lp11[loc1,ct])
          p2110 <- expit(model$lp01[loc2,ct])
          p1101 <- expit(model$lp01[loc1,ct])
          p2101 <- expit(model$lp11[loc2,ct])
          p1111 <- expit(model$lp11[loc1,ct])
          p2111 <- expit(model$lp11[loc2,ct])
        }
        tmB[1,1:4] <<- c((1-p1100)*(1-p2100),
                         p1100*(1-p2100),
                         (1-p1100)*p2100,
                         p1100*p2100)
        tmB[2,1:4] <<- c((1-p1110)*(1-p2110),
                         p1110*(1-p2110),
                         (1-p1110)*p2110,
                         p1110*p2110)
        tmB[3,1:4] <<- c((1-p1101)*(1-p2101),
                         p1101*(1-p2101),
                         (1-p1101)*p2101,
                         p1101*p2101)
        tmB[4,1:4] <<- c((1-p1111)*(1-p2111),
                         p1111*(1-p2111),
                         (1-p1111)*p2111,
                         p1111*p2111)
        
        p <- t(tmB[1:4,1:4]) %*% asCol(q_tm1[1:4])
        
        # #testing
        # model$S[loc1, ct-1] <- 1
        # model$S[loc2, ct-1] <- 0
        # model$calculate(model$getDependencies(c("S[13, 7]","S[18, 7]")))
        
        #now calculate forward probs+liklihood (note p(y|s) cancels if y>0 so only add in if y=0)
        s1 <- model$S[loc1, ct]
        s2 <- model$S[loc2, ct]
        
        #calc (0,0) 
        if(s1==0 & s2==0){
          #log likelikihood is 0 since p(y=0|s=0)=1 -> ll=log(1*1)
          lf[zt,1] <<- model$getLogProb(f_dependencies[f_start_depend[zt]:f_end_depend[zt]])+0
        }else{
          if(model$y[loc1, ct]==0 & model$y[loc2, ct]==0){
            ff1 <- 0-s1
            ff2 <- 0-s2
            loS <- 0
            if(lnei1_only>0){
              for(j in 1:lnei1_only){
                c_nei <- nei1_only[j]
                if(model$S[c_nei,ct]==1){
                  add <- ff1*model$alpha[12]
                }else{
                  add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            if(lnei2_only>0){
              for(j in 1:lnei2_only){
                c_nei <- nei2_only[j]
                if(model$S[c_nei,ct]==1){
                  add <- ff2*model$alpha[12]
                }else{
                  add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            if(lnei_both>0){
              for(j in 1:lnei_both){
                c_nei <- nei_both[j]
                if(model$S[c_nei,ct]==1){
                  add <- (ff1+ff2)*model$alpha[12]
                }else{
                  add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            #log likelikihood is 0 since p(y=0|s=0)=1 -> ll=log(1*1)
            lf[zt,1] <<- loS+0
          }else{
            lf[zt,1] <<- -Inf
          }
        }
        
        #calc (1,0) 
        #log likelihood depends on y1
        if(model$y[loc1, ct]==0){
          #ll <- -model$mup[loc1,ct]
          ll <- model$r[loc1]*(log(model$r[loc1])-log(model$r[loc1]+model$mup[loc1,ct]))
        }else{
          ll <- 0
        }
        if(s1==1 & s2==0){
          lf[zt,2] <<- model$getLogProb(f_dependencies[f_start_depend[zt]:f_end_depend[zt]])+ll
        }else{
          if(model$y[loc2, ct]==0){
            ff1 <- 1-s1
            ff2 <- 0-s2
            loS <- 0
            if(lnei1_only>0){
              for(j in 1:lnei1_only){
                c_nei <- nei1_only[j]
                if(model$S[c_nei,ct]==1){
                  add <- ff1*model$alpha[12]
                }else{
                  add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            if(lnei2_only>0){
              for(j in 1:lnei2_only){
                c_nei <- nei2_only[j]
                if(model$S[c_nei,ct]==1){
                  add <- ff2*model$alpha[12]
                }else{
                  add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            if(lnei_both>0){
              for(j in 1:lnei_both){
                c_nei <- nei_both[j]
                if(model$S[c_nei,ct]==1){
                  add <- (ff1+ff2)*model$alpha[12]
                }else{
                  add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            lf[zt,2] <<- loS+ll
          }else{
            lf[zt,2] <<- -Inf
          }
        }
        
        #calc (0,1) 
        #log likelihood depends on y1
        if(model$y[loc2, ct]==0){
          #ll <- -model$mup[loc2,ct]
          ll <- model$r[loc2]*(log(model$r[loc2])-log(model$r[loc2]+model$mup[loc2,ct]))
        }else{
          ll <- 0
        }
        if(s1==0 & s2==1){
          lf[zt,3] <<- model$getLogProb(f_dependencies[f_start_depend[zt]:f_end_depend[zt]])+ll
        }else{
          if(model$y[loc1, ct]==0){
            ff1 <- 0-s1
            ff2 <- 1-s2
            loS <- 0
            if(lnei1_only>0){
              for(j in 1:lnei1_only){
                c_nei <- nei1_only[j]
                if(model$S[c_nei,ct]==1){
                  add <- ff1*model$alpha[12]
                }else{
                  add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            if(lnei2_only>0){
              for(j in 1:lnei2_only){
                c_nei <- nei2_only[j]
                if(model$S[c_nei,ct]==1){
                  add <- ff2*model$alpha[12]
                }else{
                  add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            if(lnei_both>0){
              for(j in 1:lnei_both){
                c_nei <- nei_both[j]
                if(model$S[c_nei,ct]==1){
                  add <- (ff1+ff2)*model$alpha[12]
                }else{
                  add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
                }
                newp <- expit(logit(model$p1[c_nei,ct+1])+add)
                if(model$S[c_nei,ct+1]==1){
                  loS <- loS+log(newp)
                }else{
                  loS <- loS+log(1-newp)
                }
                #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
              }
            }
            lf[zt,3] <<- loS+ll
          }else{
            lf[zt,3] <<- -Inf
          }
        }
        
        #calc (1,1) 
        #log likelihood depends on y1 and y2
        if(model$y[loc1, ct]==0 & model$y[loc2, ct]==0){
          #ll <- -model$mup[loc1,ct]-model$mup[loc2,ct]
          ll <- model$r[loc1]*(log(model$r[loc1])-log(model$r[loc1]+model$mup[loc1,ct]))+model$r[loc2]*(log(model$r[loc2])-log(model$r[loc2]+model$mup[loc2,ct]))
        }else if(model$y[loc1, ct]==0){
          #ll <- -model$mup[loc1,ct]
          ll <- model$r[loc1]*(log(model$r[loc1])-log(model$r[loc1]+model$mup[loc1,ct]))
        }else if(model$y[loc2, ct]==0){
          #ll <- -model$mup[loc2,ct]
          ll <-  model$r[loc2]*(log(model$r[loc2])-log(model$r[loc2]+model$mup[loc2,ct]))
        }
        if(s1==1 & s2==1){
          lf[zt,4] <<- model$getLogProb(f_dependencies[f_start_depend[zt]:f_end_depend[zt]])+ll
        }else{
          
          ff1 <- 1-s1
          ff2 <- 1-s2
          loS <- 0
          if(lnei1_only>0){
            for(j in 1:lnei1_only){
              c_nei <- nei1_only[j]
              if(model$S[c_nei,ct]==1){
                add <- ff1*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,ct+1])+add)
              if(model$S[c_nei,ct+1]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
            }
          }
          if(lnei2_only>0){
            for(j in 1:lnei2_only){
              c_nei <- nei2_only[j]
              if(model$S[c_nei,ct]==1){
                add <- ff2*model$alpha[12]
              }else{
                add <- ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,ct+1])+add)
              if(model$S[c_nei,ct+1]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
            }
          }
          if(lnei_both>0){
            for(j in 1:lnei_both){
              c_nei <- nei_both[j]
              if(model$S[c_nei,ct]==1){
                add <- (ff1+ff2)*model$alpha[12]
              }else{
                add <- ff1*(model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[c_nei]))+ff2*(model$alpha[6]+model$alpha[15]*log(model$pop[loc2]*model$pop[c_nei]))
              }
              newp <- expit(logit(model$p1[c_nei,ct+1])+add)
              if(model$S[c_nei,ct+1]==1){
                loS <- loS+log(newp)
              }else{
                loS <- loS+log(1-newp)
              }
              #loS <- loS+model$S[c_nei,ct+1]*log(newp)+(1-model$S[c_nei,ct+1])*log(1-newp)
            }
          }
          lf[zt,4] <<- loS+ll
        }
        
        nl <- lf[zt,1:4]+log(p[1:4,1])
        nls <- nl-max(nl)
        q[zt,1:4] <<- exp(nls)/sum(exp(nls))
      }
    }
    
    
    if(startbs==Mt){
      #now final filter
      
      if(times[numtimes-1]==Mt-1){
        q_tm1 <- q[numtimes-1,1:4]
      }else{
        q_tm1 <- c(0,0,0,1)
      }
      
      #calc final transition matrix
      s1 <- model$S[loc1, Mt]
      s2 <- model$S[loc2, Mt]
      s1m <- model$S[loc1, Mt-1]
      s2m <- model$S[loc2, Mt-1]
      alpha6_new <- model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[loc2])
      #p(loc)(value)(condloc1loc2)
      if(neilocs==1){
        #(0,0)
        p1100 <- expit(model$lp01[loc1,Mt]+(0-s2m)*alpha6_new)
        p2100 <- expit(model$lp01[loc2,Mt]+(0-s1m)*alpha6_new)
        #(1,0)
        p1110 <- expit(model$lp11[loc1,Mt]+(0-s2m)*model$alpha[12])
        p2110 <- expit(model$lp01[loc2,Mt]+(1-s1m)*alpha6_new)
        #(0,1)
        p1101 <- expit(model$lp01[loc1,Mt]+(1-s2m)*alpha6_new)
        p2101 <- expit(model$lp11[loc2,Mt]+(0-s1m)*model$alpha[12])
        #(1,1)
        p1111 <- expit(model$lp11[loc1,Mt]+(1-s2m)*model$alpha[12])
        p2111 <- expit(model$lp11[loc2,Mt]+(1-s1m)*model$alpha[12])
        
      }else{
        p1100 <- expit(model$lp01[loc1,Mt])
        p2100 <- expit(model$lp01[loc2,Mt])
        p1110 <- expit(model$lp11[loc1,Mt])
        p2110 <- expit(model$lp01[loc2,Mt])
        p1101 <- expit(model$lp01[loc1,Mt])
        p2101 <- expit(model$lp11[loc2,Mt])
        p1111 <- expit(model$lp11[loc1,Mt])
        p2111 <- expit(model$lp11[loc2,Mt])
      }
      tmB[1,1:4] <<- c((1-p1100)*(1-p2100),
                       p1100*(1-p2100),
                       (1-p1100)*p2100,
                       p1100*p2100)
      tmB[2,1:4] <<- c((1-p1110)*(1-p2110),
                       p1110*(1-p2110),
                       (1-p1110)*p2110,
                       p1110*p2110)
      tmB[3,1:4] <<- c((1-p1101)*(1-p2101),
                       p1101*(1-p2101),
                       (1-p1101)*p2101,
                       p1101*p2101)
      tmB[4,1:4] <<- c((1-p1111)*(1-p2111),
                       p1111*(1-p2111),
                       (1-p1111)*p2111,
                       p1111*p2111)
      
      p <- t(tmB[1:4,1:4]) %*% asCol(q_tm1[1:4])
      
      ##calc likelihood
      # (0,0)
      if(model$y[loc1, Mt]==0 & model$y[loc2, Mt]==0){
        lf[numtimes,1] <<- 0
      }else{
        lf[numtimes,1] <<- -Inf
      }
      # (1,0)
      if(model$y[loc2, Mt]==0 ){
        if(model$y[loc1, Mt]==0){
          #lf[numtimes,2] <<- -model$mup[loc1,Mt]
          lf[numtimes,2] <<- model$r[loc1]*(log(model$r[loc1])-log(model$r[loc1]+model$mup[loc1,Mt]))
        }else{
          lf[numtimes,2] <<- 0
        }
      }else{
        lf[numtimes,2] <<- -Inf
      }
      # (0,1)
      if(model$y[loc1, Mt]==0 ){
        if(model$y[loc2, Mt]==0){
          #lf[numtimes,3] <<- -model$mup[loc2,Mt]
          lf[numtimes,3] <<- model$r[loc2]*(log(model$r[loc2])-log(model$r[loc2]+model$mup[loc2,Mt]))
        }else{
          lf[numtimes,3] <<- 0
        }
      }else{
        lf[numtimes,3] <<- -Inf
      }
      # (1,1)
      if(model$y[loc1, Mt]==0 & model$y[loc2, Mt]==0){
        #lf[numtimes,4] <<- -model$mup[loc1,Mt]-model$mup[loc2,Mt]
        lf[numtimes,4] <<- model$r[loc1]*(log(model$r[loc1])-log(model$r[loc1]+model$mup[loc1,Mt]))+model$r[loc2]*(log(model$r[loc2])-log(model$r[loc2]+model$mup[loc2,Mt]))
      }else if(model$y[loc1, Mt]==0){
        #lf[numtimes,4] <<- -model$mup[loc1,Mt]
        lf[numtimes,4] <<- model$r[loc1]*(log(model$r[loc1])-log(model$r[loc1]+model$mup[loc1,Mt]))
      }else if(model$y[loc2, Mt]==0){
        #lf[numtimes,4] <<- -model$mup[loc2,Mt]
        lf[numtimes,4] <<- model$r[loc2]*(log(model$r[loc2])-log(model$r[loc2]+model$mup[loc2,Mt]))
      }
      
      nl <- lf[numtimes,1:4]+log(p[1:4,1])
      nls <- nl-max(nl)
      q[numtimes,1:4] <<- exp(nls)/sum(exp(nls))
    }
    
    #first backward sample with Mt
    if(startbs==Mt){
      prev <- c(model$S[loc1,Mt],model$S[loc2,Mt])
      news <- rcat(1,q[numtimes,1:4])
      
      if(news==1){
        model$S[loc1,Mt] <<- 0
        model$S[loc2,Mt] <<- 0
      }else if(news==2){
        model$S[loc1,Mt] <<- 1
        model$S[loc2,Mt] <<- 0
      }else if(news==3){
        model$S[loc1,Mt] <<- 0
        model$S[loc2,Mt] <<- 1
      }else{
        model$S[loc1,Mt] <<- 1
        model$S[loc2,Mt] <<- 1
      }
      
      if(model$S[loc1,Mt]!= prev[1]){
        model$calculate(nodes = dependencies1[start_depend1[numtimes]:end_depend1[numtimes]])
      }
      if(model$S[loc2,Mt]!= prev[2]){
        model$calculate(nodes = dependencies2[start_depend2[numtimes]:end_depend2[numtimes]])
      }
    }
    
    #now rest of backward sample]
    for(izt in startizt:(numtimes)){
      zt <- numtimes-izt+1
      ct <- times[zt]
      
      prev <- c(model$S[loc1,ct],model$S[loc2,ct])
      
      #calculate the forward state
      if(model$S[loc1,ct+1]==1 & model$S[loc2,ct+1]==1){
        fs <- 4
      }else if(model$S[loc1,ct+1]==1 & model$S[loc2,ct+1]==0){
        fs <- 2
      }else if(model$S[loc1,ct+1]==0 & model$S[loc2,ct+1]==1){
        fs <- 3
      }else{
        fs <- 1
      }
      
      #calculate transition matrix for ct+1
      s1 <- model$S[loc1, ct+1]
      s2 <- model$S[loc2, ct+1]
      s1m <- model$S[loc1, ct]
      s2m <- model$S[loc2, ct]
      alpha6_new <- model$alpha[6]+model$alpha[15]*log(model$pop[loc1]*model$pop[loc2])
      #p(loc)(value)(condloc1loc2)
      if(neilocs==1){
        #(0,0)
        p1100 <- expit(model$lp01[loc1,ct+1]+(0-s2m)*alpha6_new)
        p2100 <- expit(model$lp01[loc2,ct+1]+(0-s1m)*alpha6_new)
        #(1,0)
        p1110 <- expit(model$lp11[loc1,ct+1]+(0-s2m)*model$alpha[12])
        p2110 <- expit(model$lp01[loc2,ct+1]+(1-s1m)*alpha6_new)
        #(0,1)
        p1101 <- expit(model$lp01[loc1,ct+1]+(1-s2m)*alpha6_new)
        p2101 <- expit(model$lp11[loc2,ct+1]+(0-s1m)*model$alpha[12])
        #(1,1)
        p1111 <- expit(model$lp11[loc1,ct+1]+(1-s2m)*model$alpha[12])
        p2111 <- expit(model$lp11[loc2,ct+1]+(1-s1m)*model$alpha[12])
        
      }else{
        p1100 <- expit(model$lp01[loc1,ct+1])
        p2100 <- expit(model$lp01[loc2,ct+1])
        p1110 <- expit(model$lp11[loc1,ct+1])
        p2110 <- expit(model$lp01[loc2,ct+1])
        p1101 <- expit(model$lp01[loc1,ct+1])
        p2101 <- expit(model$lp11[loc2,ct+1])
        p1111 <- expit(model$lp11[loc1,ct+1])
        p2111 <- expit(model$lp11[loc2,ct+1])
      }
      tmB[1,1:4] <<- c((1-p1100)*(1-p2100),
                       p1100*(1-p2100),
                       (1-p1100)*p2100,
                       p1100*p2100)
      tmB[2,1:4] <<- c((1-p1110)*(1-p2110),
                       p1110*(1-p2110),
                       (1-p1110)*p2110,
                       p1110*p2110)
      tmB[3,1:4] <<- c((1-p1101)*(1-p2101),
                       p1101*(1-p2101),
                       (1-p1101)*p2101,
                       p1101*p2101)
      tmB[4,1:4] <<- c((1-p1111)*(1-p2111),
                       p1111*(1-p2111),
                       (1-p1111)*p2111,
                       p1111*p2111)
      
      trans <- tmB[1:4,fs]
      lp <- log(trans)+log(q[zt,1:4])
      lp <- lp-max(lp)
      b <- exp(lp)/sum(exp(lp))
      news <- rcat(1,b)
      
      if(news==1){
        model$S[loc1,ct] <<- 0
        model$S[loc2,ct] <<- 0
      }else if(news==2){
        model$S[loc1,ct] <<- 1
        model$S[loc2,ct] <<- 0
      }else if(news==3){
        model$S[loc1,ct] <<- 0
        model$S[loc2,ct] <<- 1
      }else{
        model$S[loc1,ct] <<- 1
        model$S[loc2,ct] <<- 1
      }
      
      if(model$S[loc1,ct]!= prev[1]){
        model$calculate(nodes = dependencies1[start_depend1[zt]:end_depend1[zt]])
      }
      if(model$S[loc2,ct]!= prev[2]){
        model$calculate(nodes = dependencies2[start_depend2[zt]:end_depend2[zt]])
      }
      
    }
    
    copy(from = model, to = mvSaved, row = 1, 
         nodes = calcNodes, logProb = TRUE)
    
  },
  
  methods = list(   reset = function () {}   )
  
)
