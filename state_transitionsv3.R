
# Susceptible -----------------------------------------------------------------------
s_trans <- function(n) {

states = 5#(matrix(n,ncol=1)* 5);
delays = Inf#Inf*matrix(n,ncol=1)
#s_out=data.frame(states, delays)
s_out=matrix(c(states, delays), ncol=2)
return(s_out)
}


# Fast Latent ----------------------------------------------------------------------
fl_trans <- function(n, timestep, fltimes, fast_progression_prob){
states=c()
times=c()
fast_prog = as.numeric(runif(n) > fast_progression_prob) # for each n see if random number drawn from uniform is less than fast prog pro 
states = 0#matrix(0,nrow=n,ncol=1)
delays = 60#5*(12/timestep)*matrix(n,ncol=1)
if (fast_prog==1){
  states=2
  delays=fltimes[kk+1]
  kk<<-kk+1
} 

#fl_out=data.frame(states, delays)
fl_out=matrix(c(states, delays), ncol=2)
return(fl_out)
}

# LL ----------------------------------------------------------------------
ll_trans<-function(lltimes, n, timestep){
states = 2#(matrix(n,ncol=1) * 2)
delays = lltimes[(gg+1):(gg+n)]
gg<<-gg+n
#ll_out = data.frame(states, delay)
ll_out = matrix(c(states, delays), ncol=2)

return(ll_out)
}

# cpos --------------------------------------------------------------------
c_trans<- function(recoverytimes, deathtimes, txtimes, n, timestep){
states=c()
delays=c()
state_labels = c(6, 4, 3)
recovery_time = recoverytimes[(cc+1):(cc+n)]
death_time= deathtimes[(cc+1):(cc+n)]
treatment_time = txtimes[(cc+1):(cc+n)]
cc<<-cc+n;

if (sum(is.na(death_time))>0){
  state_labels = c(6, 3)
  #df=data.frame(recovery_time,treatment_time)
  df =matrix(c(recovery_time, treatment_time), ncol=2)
  for (i in 1:n){
    state = which(df[i,] == min(df[i,])) #which indices have min delay time per row
    if (length(state) > 1) {# if there are multiple mins, take sample of only 1
      state = sample(state,1);
    } 
    time = df[i,state] #take actual time
    states[i] = state_labels[state]
    delays[i] = time
  }  
} else if (sum(!is.na(death_time))>0){
df=matrix(c(recovery_time, death_time, treatment_time), ncol=3) #data.frame(recovery_time,death_time,treatment_time)
for (i in 1:n){
stateTemp = which(df[i,] == min(df[i,])) #which indices have min delay time per row
if (length(stateTemp) == 1) {
states[i] =stateTemp
} else if (length(stateTemp) > 1) {# if there are multiple mins, take sample of only 1
states[i] = sample(states[i],1);
} 
}
a<-data.frame(1:n, states)
delays = df[cbind(a[,1],a[,2])] #take actual time
states = state_labels[states]
delays[i] = time
}  
#c_out =data.frame(states, delays)
c_out =matrix(c(states, delays), ncol=2)
return(c_out)
}

# tx ----------------------------------------------------------------------
tx_trans<-function(treatment_duration,n, timestep){
states = (matrix(n,ncol=1) * 6)
delays = 6#(treatment_duration*12/timestep)*ones(n,1);
#tx_out =data.frame(states, delays)
tx_out =matrix(c(states, delays), ncol=2)
return(tx_out)
}

# D -----------------------------------------------------------------------
d_trans<-function(n){
  states = (matrix(n,ncol=1) * 4)
  delays = Inf*matrix(n,ncol=1)#(treatment_duration*12/timestep)*ones(n,1);
  #d_out =data.frame(states, delays)
  d_out =matrix(c(states, delays), ncol=2)
  return(d_out)
}


# D times -----------------------------------------------------------------

d_transT<-function(life_expectancy,n, timestep){
death_p = 1.0-exp(-(life_expectancy)/(12/timestep));
delays = (rgeom(n,death_p) + 1);
d_out = delays
return(d_out)
}

# R -----------------------------------------------------------------------
r_trans<-function(n){
  states = (matrix(n,ncol=1) * 6)
  delays = Inf*matrix(n,ncol=1)#(treatment_duration*12/timestep)*ones(n,1);
  #r_out =data.frame(states, delays)
  r_out =matrix(c(states, delays), ncol=2)
  return(r_out)
}


