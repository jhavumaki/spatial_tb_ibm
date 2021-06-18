rm(list=ls())
rm(.Random.seed, envir=globalenv())
library(plyr)
library(dplyr)
library(profvis)
library(fastmatch)
library(purrr)
dir<-"" #point to where the household setup files are
source("state_transitionsv3.R") #state transitions
#load files
paramLHS<-read.csv('lhs_10K_210419.csv') #read in sampled parameters
paramLHS<-paramLHS[,-1]
steady_state=11000#6750
numSamples =40000000#30000000
timestep=1
tspan = 12000#8250#2500# 6800# 6840 
totalrandseed=1

#determine range of parameter sets that will be used, can be informed by array_id from a slurm script, designed for running one at a time
minLHS=1#1*(ARRAYID-1)+1
maxLHS=1#ARRAYID*1
numLHS=maxLHS-minLHS+1

incidence=c();
ARTI=c();
AhhRTI=c();
AnnualCInf=c();
AnnualHHInf=c();
newTx=matrix(NA, ncol= numLHS*totalrandseed, nrow=tspan)

#annual checks
incidence=matrix(NA, ncol= numLHS*totalrandseed, nrow=tspan/12)
ARTI=matrix(NA, ncol= numLHS*totalrandseed, nrow=tspan/12)
AhhRTI=matrix(NA, ncol= numLHS*totalrandseed, nrow=tspan/12)
AnnualCInf=matrix(NA, ncol= numLHS*totalrandseed, nrow=tspan/12)
AnnualHHInf=matrix(NA, ncol= numLHS*totalrandseed, nrow=tspan/12)

HH_ATB=c();
PropHH=matrix(NA, nrow=numLHS,ncol=totalrandseed)

#set persistent variables in transmission function
kk<<-0
gg<<-0
cc<<-0
Sys.time()

#load spatial area
hID<-sample(1:10, 1)
hh_data<-read.csv(paste(dir,"hh_df", hID, ".csv", sep="")) 
hh_data<-hh_data[,-1]
dist1<-read.csv(paste(dir,"distance_table", hID, ".csv", sep="")) 
print(dist1[1:10,])
dist1<-dist1[,-1]

  for (l in minLHS:maxLHS){
    
    
    print(l)
    rm(.Random.seed)
    rseed=round(runif(1)*10^9)
    set.seed(rseed)
    print(rseed)
    
    #set params from LHS
    theta = paramLHS[1,l]
    epses = paramLHS[2:6,l]
    tau = paramLHS[7,l]
    gamma =paramLHS[8,l]
    kappa =paramLHS[9,l]
    p_ipt =paramLHS[10,l]
    CDR=paramLHS[11,l] 
    ATB_treat = (CDR*kappa+CDR*theta+CDR*gamma)/(1-CDR) 
    treatment_duration = 12/paramLHS[12,l] # convert to per year
    ipt_duration = 12/paramLHS[13,l] # convert to per year 
    spatialRange<-paramLHS[24,l]#in KM
    dist2<-dist1[dist1$d<=spatialRange,] #remove all entries > spatial range
    centroid_totals<-aggregate(hh_data$size, by=list(Category=hh_data$centroid), FUN=sum)
    betaHH = paramLHS[14,l]/(mean(hh_data$size)-1)
    dist2[dist2$d==0,"d"]<-0.05
    meanG_dist<-mean(aggregate(dist2$d, by=list(Category=dist2$a), FUN= function(x) sum(x^(-paramLHS[28,l])))[,2]) 
    betaComm= paramLHS[26,l]*(paramLHS[14,l]/(meanG_dist*mean(centroid_totals$x)))
    inf_multiplier = paramLHS[15:23,l]
    init_N= sum(hh_data[,3])
    maxID = init_N
    IDs = 1:maxID 
    hh_indiv_start=cumsum(hh_data$size)

    #sample times to next state
      monthly_rates=rep(epses/(12/timestep), each=(12/timestep))
      fast_progression_prob = 1-sum(monthly_rates) #risk of progression to LL
      fltimes<-sample(1:((12/timestep)*5), numSamples, prob=c(monthly_rates), replace=TRUE)
      lltimes = rgeom(numSamples,prob=1.0-exp(-tau/(12/timestep)))+1
      recoverytimes =rgeom(numSamples,prob=1.0-exp(-gamma/(12/timestep)))+1
      ATBdeathtimes=rgeom(numSamples,prob=1.0-exp(-kappa/(12/timestep)))+ 1
      txtimes =rgeom(numSamples,prob=1.0-exp(-ATB_treat/(12/timestep))) + 1
      deathtimes = NA+matrix(0,numSamples,1)
      
      hh_id<-rep(1:nrow(hh_data), hh_data$size)
      centroid_ID<-rep(hh_data$centroid, hh_data$size)
      
   
      
      #Initial States
      numC=round(paramLHS[25,l]);#367
      curr_state<-matrix(5, nrow=maxID, ncol=1)
      eligibleCs<-dist1[dist1$a==paramLHS[29,l] & dist1$d<paramLHS[27,l], "b"] #i.e., these are within seed range from randomly selected centroid
       numCs<-sample(which(centroid_ID %in% eligibleCs), numC)
      curr_state[numCs]<-2
      
      
      
      S_count=sum(curr_state==5);
      E_count=sum(curr_state==1);
      L_count=sum(curr_state==0);
      C_count=sum(curr_state==2);
      T_count=sum(curr_state==3);
      D_count=sum(curr_state==4);
      R_count=sum(curr_state==6);
      PTs_count=sum(curr_state==7);
      PTl_count=sum(curr_state==8);
      
      state_ts=data.frame(L_count, E_count, C_count, T_count,  D_count, S_count, R_count,  PTs_count,  PTl_count)
      colnames(state_ts)<-c('L', 'E','C','T','D','S','R','PTs','PTl')
      
      state_ts1<-state_ts
      
      #draw next state
      draw.next.state <- function(curr_state,curr_t, timestep,fltimes, recoverytimes, ATBdeathtimes, txtimes, lltimes,
                                  deathtimes, fast_progression_prob, treatment_duration, ipt_duration) {
        if(curr_state==5){
          z<-matrix(c(5,Inf), ncol=2)#s_trans(1)
        } else if (curr_state==1){
          z<-fl_trans(1, timestep, fltimes, fast_progression_prob)
        } else if (curr_state==0){
          z<-ll_trans(lltimes, 1, timestep)
        } else if (curr_state==2){
          z<-c_trans(recoverytimes, deathtimes, txtimes, 1, timestep)
        } else if (curr_state==3){
          z<-tx_trans(treatment_duration,1, timestep)
        }else if (curr_state==4){
          z<-matrix(c(4,Inf), ncol=2)#d_trans(1)
        } else if (curr_state==6){
          z<-matrix(c(6,Inf), ncol=2)#r_trans(1)
        }
        return(z)
      }
      
      #because curr_state is only 5 and 2:
      init_states<-matrix(c(5, Inf), ncol=2, nrow=init_N, byrow=T)
      
      
      update_states<-c_trans(recoverytimes, deathtimes, txtimes, sum(curr_state==2), timestep) 
      init_states[which(curr_state==2),1]<- update_states[,1]#update_states$states
      init_states[which(curr_state==2),2]<- update_states[,2]#update_states$delays
      
      next_state = init_states[,1]
      next_event_t =init_states[,2]
      
      #calculate death time
      death_time= d_transT(theta,sum(state_ts), timestep)
      
      state_df = data.frame(IDs, curr_state, next_state, next_event_t, death_time, hh_id, centroid_ID)
      
      curr_state =c()
      init_states =c()
      next_state =c()
      next_event_t=c()
      hh_id = c()
      incATB=c()
      atb_pos=c()
      atb_all<-c()
      deaths<-c()
      hhInf<-c()
      commInf<-c()
      cumI=0
      curr_t=c()
      y=1
      output<-matrix(0, nrow=nrow(centroid_totals), ncol=tspan-steady_state )
      outputInc<-matrix(0, nrow=nrow(centroid_totals), ncol=tspan-steady_state )
      detected<-matrix(0,nrow=tspan-steady_state, ncol=3  )
      HH_ATB<-matrix(0, nrow=tspan-steady_state, ncol=2 )
      ind_foi1=matrix(0, nrow=maxID,ncol=1)
      foi_hh1=matrix(0,init_N,1)
      foi_comm1=matrix(0,init_N,1)
      
      #start simulation
      for (curr_t in 1:tspan){
       
        a<-state_df[,2] #time trim
        infectedsC=which(a==2)#time trim
        ind_foi=ind_foi1
        
        if (length(infectedsC)>0){
          #efficient Household FOI  
          foi_hh=foi_hh1
          numC<-data.frame(table(state_df[infectedsC,"hh_id"]))
          numC$Var1<-as.numeric(as.character(numC$Var1))
          numC$Freq<-as.numeric(as.character(numC$Freq))
          hh<-matrix(F,nrow=max(state_df$hh_id), ncol=2)
          hh[numC$Var1, 1]<-T
          hh[numC$Var1, 2]<-numC$Freq*betaHH
          
          foi_hh<-hh[state_df$hh_id,2]
         
          #Community FOI
          foi_comm=foi_comm1
          numC<-data.frame(table(state_df[infectedsC, "centroid_ID"]))
          numC$Var1<-as.numeric(as.character(numC$Var1))
          numC$Freq<-as.numeric(as.character(numC$Freq))
          
          #force of infection is from column a to b
          distC<-matrix(F,nrow=max(state_df$centroid_ID), ncol=1)
          distC[numC$Var1]<-T
          distC<-dist2[distC[dist2$a],]
          distC$numCa<-numC[fmatch(distC$a,numC$Var1),"Freq"] #count number of C per centroid in a, 
          distC$foi<-(distC$d^(-paramLHS[28,l]))*betaComm*distC$numCa # calculate Guassian kernel min>0 since already imposed spatial range
          centroidFoi<-data.frame(distC %>% group_by(b) %>% summarise(foi = sum(foi))) #sum FOI imposed upon each centroid in column B, checked
           
          
          
          Cent<-matrix(F,nrow=max(state_df$centroid_ID), ncol=2)
          Cent[centroidFoi$b,1]<-T
          Cent[centroidFoi$b,2]<-centroidFoi$foi
          foi_comm<-Cent[state_df$centroid_ID,2]

          #add FOIs
          foi= foi_hh+foi_comm 
          
          #incorporte immunity
          statearray=state_df[,2]
          
          
          #L, E, C, T, D, S, R, PTs, Ptl
          ind_foi<-(1.0-exp(-1*foi*inf_multiplier[statearray+1])) 
          
        }
        
       
        newLatent =which(runif(init_N)<ind_foi) #note bool
        

        
        state_df[newLatent,c(3, 4)]= data.frame(1, curr_t)
        
        # Get the group of individuals who have a transition on this time step
        transition_ids = which(state_df[,4] == curr_t); #new latent, but also anyone transitioning at this time step
        num_transitions = length(transition_ids);
        
        incATB = sum(state_df[transition_ids,3] ==2) #new ATB

        
        #keep track of how many incident cases  (state 2) and how many case notifications (state 3)
        if (curr_t>steady_state){  
      
          curr_cases_centroid<-table(state_df[transition_ids,"centroid_ID"], state_df[transition_ids,"next_state"])
          if (sum(colnames(curr_cases_centroid)=="3")==1){
            curr_cases_centroid3<- data.frame(rownames(curr_cases_centroid), curr_cases_centroid[,which(colnames(curr_cases_centroid)=="3")])
            colnames(curr_cases_centroid3)<-c("centroid_ID", "n")
            curr_cases_centroid3$centroid_ID<-as.numeric(as.character(curr_cases_centroid3$centroid_ID))
            curr_cases_centroid3$n<-as.numeric(as.character(curr_cases_centroid3$n))
            output[curr_cases_centroid3$centroid_ID,curr_t-steady_state]<-curr_cases_centroid3$n
          }
          
          
          if (sum(colnames(curr_cases_centroid)=="2")==1){
            curr_cases_centroid2<- data.frame(rownames(curr_cases_centroid), curr_cases_centroid[,which(colnames(curr_cases_centroid)=="2")])
            colnames(curr_cases_centroid2)<-c("centroid_ID", "n")
            curr_cases_centroid2$centroid_ID<-as.numeric(as.character(curr_cases_centroid2$centroid_ID))
            curr_cases_centroid2$n<-as.numeric(as.character(curr_cases_centroid2$n))
            outputInc[curr_cases_centroid2$centroid_ID,curr_t-steady_state]<-curr_cases_centroid2$n
          }
          
          
          
        }
        
        
        
        if (num_transitions > 0) {
          #Now draw a set of next states for these individuals
          state_df[transition_ids,2]=state_df[transition_ids,3]
          
          states1=state_df[transition_ids,2]
          curr_time=rep(curr_t, length(states1))
          
          next_states<-t(data.frame(sapply(states1, function(x) draw.next.state(x,curr_t, timestep,fltimes, recoverytimes, ATBdeathtimes, txtimes, lltimes, deathtimes, fast_progression_prob, treatment_duration, ipt_duration))))
          state_df[transition_ids,3]<-unlist(next_states[,1])
          state_df[transition_ids,4]<-curr_t+unlist(next_states[,2])
        }
        
        #calculate how many cases are detected (vs. not detected)
        if (curr_t>steady_state){
          detected[curr_t-steady_state,1]<-curr_t
          detected[curr_t-steady_state,2]<-sum(state_df[transition_ids,"curr_state"]==2 & state_df[transition_ids,"next_event_t"]<=tspan)
          detected[curr_t-steady_state,3]<-sum(state_df[transition_ids,"curr_state"]==3)
        }
        
        #track co-prevalence household case
        if (curr_t>steady_state & sum(state_df[transition_ids,2]==3)>0){
          
          Newly_detected<-state_df[transition_ids[state_df[transition_ids,2]==3],]
          Newly_detected_hh<-state_df[which(state_df$hh_id %in% Newly_detected$hh_id),]
          
          if (sum(Newly_detected_hh[,2]==2)>0 | sum(Newly_detected_hh[,2]==3)>0){
            a<-table(Newly_detected_hh[,"hh_id"], Newly_detected_hh[,"curr_state"])
            aATB<-a[,colnames(a)=="2" |  colnames(a)=="3"]
            
            
            if (is.null(dim(aATB))){
              HH_ATB[curr_t-steady_state,1]<-sum(aATB-1)
              
            } else  if (!is.null(dim(aATB))){
            HH_ATB[curr_t-steady_state,1]<-sum(rowSums(aATB)-1)
            }
            HH_ATB[curr_t-steady_state,2]<-sum(rowSums(a)-1) #total number of cases in household with case already
              
             
          }
        }
        
        
        
        #vital dynamics
        disease_death_ids = state_df[,2] == 4;
        natural_death_ids = state_df[,5] == curr_t;
        
        if (sum(natural_death_ids) > 0) {
         
          
          aa<-curr_t+d_transT(theta,sum(natural_death_ids), timestep)
          ab<-which(natural_death_ids)
          state_df[ab,2:5] <-data.frame(5,5,Inf, aa)
        }
        
        #births
        if (sum(disease_death_ids)>0){
          
          state_df[which(disease_death_ids),2]<-sample(c(5,1,0,6), sum(disease_death_ids),prob=c(state_ts[curr_t,"S"],state_ts[curr_t,"E"], state_ts[curr_t,"L"], state_ts[curr_t,"R"])/init_N, replace=T)
          
          next_states<-t(data.frame(sapply(state_df[which(disease_death_ids),2], function(x) draw.next.state(x,curr_t, timestep,fltimes, recoverytimes, ATBdeathtimes, txtimes, lltimes, deathtimes, fast_progression_prob, treatment_duration, ipt_duration))))
          state_df[disease_death_ids,3]<-unlist(next_states[,1])
          state_df[disease_death_ids,4]<-curr_t+unlist(next_states[,2])
          state_df[disease_death_ids,5]<- curr_t+d_transT(theta,sum(disease_death_ids), timestep)
        }
        

      
        state_ts[curr_t+1,]<-hist(state_df[,2], breaks=seq(-1,8,1), plot=F)$counts
        
        
        atb_all[curr_t]<-incATB
        
      
          
          #reset monthly counts
          atb_pos<-c()
          
          commInf<-c()
          hhInf<-c()
          deaths<-c()
          incATB<-c()
        
        uid<-c() 
        a<-c()
        infTemp<-c()
        infP<-c()
      }
      kk<<-0
      gg<<-0
      cc<<-0
      
      #save files
      filename_ts= paste('state_ts_state_none','lhs', l, 'rng' ,rseed,'seed','HH', hID,'_none','.csv', sep="")
      write.csv(state_ts,filename_ts)
      
      filename_ts= paste('hh_atb_none','lhs', l, 'rng' ,rseed,'seed','HH', hID,'_none','.csv', sep="")
      write.csv(HH_ATB,filename_ts)
      
      
        filename_ts= paste('output_state_none', 'lhs', l, 'rng',rseed,'seed','HH', hID,'_none','.csv', sep="")
        output<-data.frame(output)
        colnames(output)<- (steady_state+1):tspan
        write.csv(output,filename_ts)
        
        filename_ts= paste('out_Inc_put_state_none', 'lhs', l, 'rng',rseed,'seed','HH', hID,'_none','.csv', sep="")
        outputInc<-data.frame(outputInc)
        colnames(outputInc)<- (steady_state+1):tspan
        write.csv(outputInc,filename_ts)
        
        
        filename_ts= paste('detected', 'lhs', l, 'rng',rseed,'seed','HH', hID,'_none','.csv', sep="")
        detected<-data.frame(detected)
        colnames(detected)<- c("curr_t", "atb", "tx")
        write.csv(detected,filename_ts)

    }
Sys.time()


