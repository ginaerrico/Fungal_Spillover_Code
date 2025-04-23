library(deSolve)
library(ReacTran)

library(abc)
library(abcrf)
library(BiocManager)
library(ggplot2)
library(dplyr)
library(lime)
library(keras)
library(tensorflow)
library(EBImage)
library(tfdatasets)
library(RColorBrewer)
library(lattice)
library(vegan)

##To create the turbulence vector
#turb<-runif(n=51, min=0.1, max=2)
#turb<-as.data.frame(turb)
#write.csv(turb,"turb.csv")
#setwd("/Users/bebache/Desktop/Spillover/Revision")
setwd("C:/Users/Zach/Downloads/GinasRStuff")
turb<-read.csv("turb.csv")
turb<-turb[,2]

############################################
# FOUR MODELS
#No spillover ------------------------------------------------------------

model0<- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx1 = x1*(alpha1 - a11*x1 - a12*x2 - (beta11*y1)/(1+e11*x1))-((beta12*y2)/(1+e12*x1))
    dx2 = x2*(alpha2 - a22*x2 - a21*x1 - (beta21*y1)/(1+e21*x2))-((beta22*y2)/(1+e22*x2))
    dy1 = y1*((delta1*beta11*x1)/(1+e11*x1)+(delta2*beta21*x2)/(1+e21*x2)-m)
    return(list(c(dx1,dx2, dy1)))
  })
}

# model1: No velocity
model1 <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    x1 <- State[1:N]
    x2 <- State[(N+1):(2*N)]
    y1<-State[((2*N)+1):(3*N)]
    y2<-State[((3*N)+1):(4*N)]
    
    dx1 = x1*(alpha1 - a11*x1 - a12*x2 - ((beta11*y1)/(1+e11*x1))-((beta12*y2)/(1+e12*x1)))
    dx2 = x2*(alpha2 - a22*x2 - a21*x1 - ((beta21*y1)/(1+e21*x2))-((beta22*y2)/(1+e22*x2)))
    dy1 = y1*((delta1*beta11*x1)/(1+e11*x1)+(delta2*beta21*x2)/(1+e21*x2)-m)
    dy2 = tran.1D (C = y2, C.up = 5, C.down = 0,D = 10, dx = Grid, v=0)$dC 
    
    return(list(c(dx1,dx2, dy1,dy2)))
  })
}

# model2: 0.25 velocity
model2 <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    x1 <- State[1:N]
    x2 <- State[(N+1):(2*N)]
    y1<-State[((2*N)+1):(3*N)]
    y2<-State[((3*N)+1):(4*N)]
    
    dx1 = x1*(alpha1 - a11*x1 - a12*x2 - ((beta11*y1)/(1+e11*x1))-((beta12*y2)/(1+e12*x1)))
    dx2 = x2*(alpha2 - a22*x2 - a21*x1 - ((beta21*y1)/(1+e21*x2))-((beta22*y2)/(1+e22*x2)))
    dy1 = y1*((delta1*beta11*x1)/(1+e11*x1)+(delta2*beta21*x2)/(1+e21*x2)-m)
    dy2 = tran.1D (C = y2, C.up = 5, C.down = 0,D = 10, dx = Grid, v=0.25)$dC 
    
    return(list(c(dx1,dx2, dy1,dy2)))
  })
}

# model3: turbulence
model3 <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    x1 <- State[1:N]
    x2 <- State[(N+1):(2*N)]
    y1<-State[((2*N)+1):(3*N)]
    y2<-State[((3*N)+1):(4*N)]
    
    dx1 = x1*(alpha1 - a11*x1 - a12*x2 - ((beta11*y1)/(1+e11*x1))-((beta12*y2)/(1+e12*x1)))
    dx2 = x2*(alpha2 - a22*x2 - a21*x1 - ((beta21*y1)/(1+e21*x2))-((beta22*y2)/(1+e22*x2)))
    dy1 = y1*((delta1*beta11*x1)/(1+e11*x1)+(delta2*beta21*x2)/(1+e21*x2)-m)
    dy2 = tran.1D (C = y2, C.up = 5, C.down = 0,D = 10, dx = Grid, v=turb)$dC 
    
    return(list(c(dx1,dx2, dy1,dy2)))
  })
}

### 0] no diffusion
#testing outcomes at different intersp. comp. levels - No diff
A21=seq(0.01,1, by=0.01) #intersp. competition 
A12=seq(0.01,1, by=0.01) 
outcome0<- c() #empty outcome vector
a_21<-c()
a_12<-c()

for(i in 1:length(A21))
{
  for(j in 1:length(A12))
  {
    Time <- seq(0, 100, by = 1)
    State <- c(x1=5, x2=5, y1=5) 
    Pars <- c(
      alpha1= 3,
      alpha2= 3,
      a11= 0.5,
      a22= 0.5,
      a21= A21[i],
      a12= A12[j],
      beta11= 0.5, #or 0.3  for specialist resident
      beta12= 0.3,
      beta21= 0.5, # or 0.7 for specialist resident
      beta22=0,
      e11= 0.1,
      e21= 0.1,
      e12= 0.1,
      e22=0.1,
      delta1= 0.1,
      delta2= 0.1,
      m=0.2,
      y2=0
    )
    out0 <- ode(y = State, func = model0,
                times = Time, parms = Pars)
    outcome0 <- rbind(outcome0, out0[nrow(out0),-1])
    a_12<-c(a_12,A12[j])
    a_21<-c(a_21,A21[i])
    
  }
}

outcome_P1<-(outcome0[,1])
outcome_P2<-(outcome0[,2])
outcome_Y1<-(outcome0[,3])


outcome_P1[is.na(outcome_P1)==TRUE]<-0
outcome_P2[is.na(outcome_P2)==TRUE]<-0
outcome_Y1[is.na(outcome_Y1)==TRUE]<-0

#printing the outcome with new cut off at 0.00001
cc<-c() #vector of color
for(i in 1:length(outcome_P2)) {
  
  if((outcome_P2[i]>0.00001 & outcome_P1[i]<=0.00001) & outcome_Y1[i]>0)
  {
    cc<-c(cc,"palegreen2") # x2 is winning
  }
  if((outcome_P2[i]<=0.00001 & outcome_P1[i]>0.00001) & outcome_Y1[i]>0)
  {
    cc<-c(cc,"darkgreen") # x1 winning
  }
  if((outcome_P2[i]>0.00001 & outcome_P1[i]>0.00001) & outcome_Y1[i]>0)
  {
    cc<-c(cc,"orange1") # coexistence
  }
  if((outcome_P2[i]<=0.00001 & outcome_P1[i]<=0.00001)| (outcome_Y1[i]<0))
  {cc<-c(cc,"lightgrey")} #crash
}

plot(a_12,a_21,col=cc,pch=20)

    Pars <- c(
      alpha1= 3,
      alpha2= 3,
      a11= 0.5,
      a22= 0.5,
      a21= 0.41,
      a12= 0.39,
      beta11= 0.5,
      beta12= 0.3,
      beta21= 0.5,
      beta22=0,
      e11= 0.1,
      e21= 0.1,
      e12= 0.1,
      e22=0.1,
      delta1= 0.1,
      delta2= 0.1,
      m=0.2, # too high. Drop to 0.2
      y2=0
    )
    out0 <- ode(y = State, func = model0,
                times = Time, parms = Pars)

### 1] Difussion without wind


N<-50
Grid <- setup.grid.1D(x.up = 0, x.down = 100, N = N)

x1ini <- rep(x = 5, times = N)
x2ini <- rep(x = 5, times = N)
y1ini <- rep(x = 5, times = N)
y2ini<-c(2,rep(x = 0, times = (N-1)))

State <- c(x1ini, x2ini,y1ini,y2ini)
Time <- seq(0, 100, by = 1)

A21=seq(0.01,1, by=0.01) #varying comp of x1 on x2
A12=seq(0.01,1, by=0.01) #
outcome01<- c() #empty outcome vector
a_21<-c()
a_12<-c()

for(i in 1:length(A21))
{
  for(j in 1:length(A12))
  {
    Pars <- c(
      alpha1= 3,
      alpha2= 3,
      a11= 0.5,
      a22= 0.5,
      a21= A21[i],
      a12= A12[j],
      beta11= 0.5,
      beta12= 0.3, #effect of spillover on P1
      beta21= 0.5,
      e11= 0.1,
      e21= 0.1,
      e12= 0.3, 
      delta1= 0.1,
      delta2= 0.1,
      m=0.2, #mortality too high. Drop it to 0.2..
      e22=0.3,
      beta22=0.3 #effect of spillover on P2 (specialist is 0)
    )
    out01 <- ode.1D(y = State, func = model1,
                    times = Time, parms = Pars, nspec = 4,
                    names = c("x1","x2","y1","y2"), dimens = N)
    outcome01 <- rbind(outcome01, out01[nrow(out01),-1])
    a_12<-c(a_12,A12[j])
    a_21<-c(a_21,A21[i])
    
  }
}

outcome_P1<-(outcome01[,1:50])
outcome_P2<-(outcome01[,51:100])
outcome_Y1<-(outcome01[,101:150])

outcome_P1[is.na(outcome_P1)==TRUE]<-0
outcome_P2[is.na(outcome_P2)==TRUE]<-0
outcome_Y1[is.na(outcome_Y1)==TRUE]<-0


for (k in 1:50) { #looping through the space
  cc<-c() #vector of color
  for(i in 1:nrow(outcome_P2))
    
  {
    if(outcome_P2[i,k]<0.00001 & outcome_P1[i,k]>0.00001 & outcome_Y1[i,k]>0)
      
    {
      cc<-c(cc,"darkgreen") # x1 is winning
    }
    if(outcome_P2[i,k]>0.00001 & outcome_P1[i,k]<0.00001 & outcome_Y1[i,k]>0)
    {
      cc<-c(cc,"palegreen2") # x2 winning
    }
    if(outcome_P2[i,k]>0.00001 & outcome_P1[i,k]>0.00001 & outcome_Y1[i,k]>0)
    {
      cc<-c(cc,"orange1") # coexistence
    }
    if((outcome_P2[i,k]<0.00001 & outcome_P1[i,k]<0.00001) | (outcome_Y1[i,k]<0))
    {cc<-c(cc,"lightgrey")} #crash
  }
  plot(a_12,a_21,pch=20,col=cc,xlab="Plant 2",ylab="Plant 1") # competition
  dev.copy(png,paste("A",k,".png",sep="")) #saving on computer
  dev.off()
}

### 2] Repeat with wind


A21=seq(0.01,1, by=0.01) #varying comp of x1 on x2
A12=seq(0.01,1, by=0.01) #
outcome01<- c() #empty outcome vector
a_21<-c()
a_12<-c()

for(i in 1:length(A21))
{
  for(j in 1:length(A12))
  {
    Pars <- c(
      alpha1= 3,
      alpha2= 3,
      a11= 0.5,
      a22= 0.5,
      a21= A21[i],
      a12= A12[j],
      beta11= 0.5,
      beta12= 0.3, #effect of spillover on P1
      beta21= 0.5,
      e11= 0.1,
      e21= 0.1,
      e12= 0.3, 
      delta1= 0.1,
      delta2= 0.1,
      m=0.2,
      e22=0.3,
      beta22=0.3 #effect of spillover on P2 (specialist is 0)
    )
    out01 <- ode.1D(y = State, func = model2,
                    times = Time, parms = Pars, nspec = 4,
                    names = c("x1","x2","y1","y2"), dimens = N)
    outcome01 <- rbind(outcome01, out01[nrow(out01),-1])
    a_12<-c(a_12,A12[j])
    a_21<-c(a_21,A21[i])
    
  }
}

outcome_P1<-(outcome01[,1:50])
outcome_P2<-(outcome01[,51:100])
outcome_Y1<-(outcome01[,101:150])

outcome_P1[is.na(outcome_P1)==TRUE]<-0
outcome_P2[is.na(outcome_P2)==TRUE]<-0
outcome_Y1[is.na(outcome_Y1)==TRUE]<-0


for (k in 1:50) { #looping through the space
  cc<-c() #vector of color
  for(i in 1:nrow(outcome_P2))
    
  {
    if(outcome_P2[i,k]<0.00001 & outcome_P1[i,k]>0.00001 & outcome_Y1[i,k]>0)
      
    {
      cc<-c(cc,"darkgreen") # x1 is winning
    }
    if(outcome_P2[i,k]>0.00001 & outcome_P1[i,k]<0.00001 & outcome_Y1[i,k]>0)
    {
      cc<-c(cc,"palegreen2") # x2 winning
    }
    if(outcome_P2[i,k]>0.00001 & outcome_P1[i,k]>0.00001 & outcome_Y1[i,k]>0)
    {
      cc<-c(cc,"orange1") # coexistence
    }
    if((outcome_P2[i,k]<0.00001 & outcome_P1[i,k]<0.00001) | (outcome_Y1[i,k]<0))
    {cc<-c(cc,"lightgrey")} #crash
  }
  plot(a_12,a_21,pch=20,col=cc,xlab="Plant 2",ylab="Plant 1") # competition
  dev.copy(png,paste("A",k,".png",sep="")) #saving on computer
  dev.off()
}
 
###  3] Repeat with turbulence


A21=seq(0.01,1, by=0.01) #varying comp of x1 on x2
A12=seq(0.01,1, by=0.01) #
outcome01<- c() #empty outcome vector
a_21<-c()
a_12<-c()

for(i in 1:length(A21))
{
  for(j in 1:length(A12))
  {
    Pars <- c(
      alpha1= 3,
      alpha2= 3,
      a11= 0.5,
      a22= 0.5,
      a21= A21[i],
      a12= A12[j],
      beta11= 0.5,
      beta12= 0.3, #effect of spillover on P1
      beta21= 0.5,
      e11= 0.1,
      e21= 0.1,
      e12= 0.3, 
      delta1= 0.1,
      delta2= 0.1,
      m=0.2,
      e22=0.3,
      beta22=0.3 #effect of spillover on P2
    )
    out01 <- ode.1D(y = State, func = model3,
                    times = Time, parms = Pars, nspec = 4,
                    names = c("x1","x2","y1","y2"), dimens = N)
    outcome01 <- rbind(outcome01, out01[nrow(out01),-1])
    a_12<-c(a_12,A12[j])
    a_21<-c(a_21,A21[i])
    
  }
}

outcome_P1<-(outcome01[,1:50])
outcome_P2<-(outcome01[,51:100])
outcome_Y1<-(outcome01[,101:150])

outcome_P1[is.na(outcome_P1)==TRUE]<-0
outcome_P2[is.na(outcome_P2)==TRUE]<-0
outcome_Y1[is.na(outcome_Y1)==TRUE]<-0


for (k in 1:50) { #looping through the space
  cc<-c() #vector of color
  for(i in 1:nrow(outcome_P2))
    
  {
    if(outcome_P2[i,k]<0.00001 & outcome_P1[i,k]>0.00001 & outcome_Y1[i,k]>0)
      
    {
      cc<-c(cc,"darkgreen") # x1 is winning
    }
    if(outcome_P2[i,k]>0.00001 & outcome_P1[i,k]<0.00001 & outcome_Y1[i,k]>0)
    {
      cc<-c(cc,"palegreen2") # x2 winning
    }
    if(outcome_P2[i,k]>0.00001 & outcome_P1[i,k]>0.00001 & outcome_Y1[i,k]>0)
    {
      cc<-c(cc,"orange1") # coexistence
    }
    if((outcome_P2[i,k]<0.00001 & outcome_P1[i,k]<0.00001) | (outcome_Y1[i,k]<0))
    {cc<-c(cc,"lightgrey")} #crash
  }
  plot(a_12,a_21,pch=20,col=cc,xlab="Plant 2",ylab="Plant 1") # competition
  dev.copy(png,paste("B",k,".png",sep="")) #saving on computer
  dev.off()
}

## repeat 1-3 with specialist spillover with beta22=0 and then 0-3 with specialist resident and both generalist or specialist spillover


### Simulating results for a large number of parameters



outcome0<- c() #empty outcome vector
Pars_save<-c() #empty parameter vector

for(i in 1:10000)
{
  
    Time <- seq(0, 100, by = 1)
    State <- c(x1=5, x2=5, y1=5) 
    Pars <- c(
      alpha1= runif(1,1,10),
      alpha2= runif(1,1,10),
      a11= runif(1,0.01,1),
      a22= runif(1,0.01,1),
      a21= runif(1,0.01,1),
      a12= runif(1,0.01,1),
      beta11= runif(1,0.01,1), # Resident pathogen on P1
      beta12= runif(1,0.01,1),
      beta21= runif(1,0.01,1), # Resident pathogen on P2
      beta22=runif(1,0.01,1),
      e11= runif(1,0.01,1),
      e21= runif(1,0.01,1),
      e12= runif(1,0.01,1),
      e22=runif(1,0.01,1),
      delta1= runif(1,0.01,1),
      delta2= runif(1,0.01,1),
      m=runif(1,0.01,1),
      y2=0
    )
    out0 <- ode(y = State, func = model0,
                times = Time, parms = Pars)
    outcome0 <- rbind(outcome0, out0[nrow(out0),-1])
    Pars_save<-rbind(Pars_save,Pars)
}

outcome_P1<-(outcome0[,1])
outcome_P2<-(outcome0[,2])
outcome_Y1<-(outcome0[,3])

outcome_P1[is.na(outcome_P1)==TRUE]<-0
outcome_P2[is.na(outcome_P2)==TRUE]<-0
outcome_Y1[is.na(outcome_Y1)==TRUE]<-0

#printing the outcome
cc<-c() #vector of color
for(i in 1:length(outcome_P2)) {
  
  if((outcome_P2[i]>0.00001 & outcome_P1[i]<=0.00001) & outcome_Y1[i]>0)
  {
    cc<-c(cc,"palegreen2") # x2 is winning
  }
  if((outcome_P2[i]<=0.00001 & outcome_P1[i]>0.00001) & outcome_Y1[i]>0)
  {
    cc<-c(cc,"darkgreen") # x1 winning
  }
  if((outcome_P2[i]>0.00001 & outcome_P1[i]>0.00001) & outcome_Y1[i]>0)
  {
    cc<-c(cc,"orange1") # coexistence
  }
  if((outcome_P2[i]<=0.00001 & outcome_P1[i]<=0.00001)| (outcome_Y1[i]<0))
  {cc<-c(cc,"lightgrey")} #crash
}

cc_0<-cc


Pars_rda<-rda(Pars_save[cc!="lightgrey",1:6]~cc[cc!="lightgrey"])
anova(Pars_rda,by="terms")

Pars_ord<-metaMDS(Pars_save[cc!="lightgrey",1:6])
plot(Pars_ord,type="n")
points(Pars_ord,col=cc[cc!="lightgrey"],pch=20)
ordiellipse(Pars_ord,cc[cc!="lightgrey"],col=cc[cc!="lightgrey"],draw="polygon")
text(Pars_ord, display = c("species"), labels=dimnames(Pars_save)[[2]][1:6])

fitness_ratio<-Pars_save[,1]/Pars_save[,2]
competition_ratio<-(Pars_save[,5]*Pars_save[,4])/(Pars_save[,6]*Pars_save[,3])
resident_ratio<-(Pars_save[,7]*Pars_save[,12])/(Pars_save[,11]*Pars_save[,9])
resident_growth_ratio<-Pars_save[,15]/Pars_save[,16]
resident_death<-Pars_save[,17]
spillover_ratio<-(Pars_save[,8]*Pars_save[,14])/(Pars_save[,13]*Pars_save[,10])
spillover_resident_ratio<-spillover_ratio/resident_ratio
parameters<-cbind(fitness_ratio,competition_ratio,resident_ratio,resident_growth_ratio,resident_death,spillover_ratio,spillover_resident_ratio)

parameters_ord<-metaMDS(parameters[cc!="lightgrey",])
plot(parameters_ord,type="n")
points(parameters_ord,col=cc[cc!="lightgrey"],pch=20)
ordiellipse(parameters_ord,cc[cc!="lightgrey"],col=cc[cc!="lightgrey"],draw="polygon")
text(parameters_ord, display = c("species"), labels=dimnames(parameters)[[2]])

# set up spatial grid and initial conditions for diffusion
N<-50
Grid <- setup.grid.1D(x.up = 0, x.down = 100, N = N)

x1ini <- rep(x = 5, times = N)
x2ini <- rep(x = 5, times = N)
y1ini <- rep(x = 5, times = N)
y2ini<-c(5,rep(x = 0, times = (N-1)))

State <- c(x1ini, x2ini,y1ini,y2ini)
Time <- seq(0, 100, by = 1)

outcome01<- c() #empty outcome vector model1
outcome02<-c() #empty outcome vector model2
outcome03<-c() #empty outcome vector model3

for(i in 1:nrow(Pars_save))
{
  
    Pars <- Pars_save[i,-18]
    out01 <- ode.1D(y = State, func = model1,
                    times = Time, parms = Pars, nspec = 4,
                    names = c("x1","x2","y1","y2"), dimens = N)
    outcome01 <- rbind(outcome01, out01[nrow(out01),-1])
    out02 <- ode.1D(y = State, func = model2,
                    times = Time, parms = Pars, nspec = 4,
                    names = c("x1","x2","y1","y2"), dimens = N)
    outcome02 <- rbind(outcome02, out02[nrow(out02),-1])
    out03 <- ode.1D(y = State, func = model3,
                    times = Time, parms = Pars, nspec = 4,
                    names = c("x1","x2","y1","y2"), dimens = N)
    outcome03 <- rbind(outcome03, out03[nrow(out03),-1])
}

image(out01,mfrow=NULL,grid=Grid$x.mid,main="y2 - Diffusion",which="y2")
image(out02,mfrow=NULL,grid=Grid$x.mid,main="y2 - Diffusion & Wind",which="y2")
image(out03,mfrow=NULL,grid=Grid$x.mid,main="y2 - Diffusion & Turbulences",which="y2")

outcome1_P1<-(outcome01[,1:50])
outcome1_P2<-(outcome01[,51:100])
outcome1_Y1<-(outcome01[,101:150])
outcome2_P1<-(outcome02[,1:50])
outcome2_P2<-(outcome02[,51:100])
outcome2_Y1<-(outcome02[,101:150])
outcome3_P1<-(outcome03[,1:50])
outcome3_P2<-(outcome03[,51:100])
outcome3_Y1<-(outcome03[,101:150])

outcome1_P1[is.na(outcome1_P1)==TRUE]<-0
outcome1_P2[is.na(outcome1_P2)==TRUE]<-0
outcome1_Y1[is.na(outcome1_Y1)==TRUE]<-0
outcome2_P1[is.na(outcome2_P1)==TRUE]<-0
outcome2_P2[is.na(outcome2_P2)==TRUE]<-0
outcome2_Y1[is.na(outcome2_Y1)==TRUE]<-0
outcome3_P1[is.na(outcome3_P1)==TRUE]<-0
outcome3_P2[is.na(outcome3_P2)==TRUE]<-0
outcome3_Y1[is.na(outcome3_Y1)==TRUE]<-0

cc_1<-c()
cc_2<-c()
cc_3<-c()

for (k in 1:50) { #looping through the space
  cc1<-c() #vector of color
  cc2<-c() #vector of color
  cc3<-c() #vector of color
  for(i in 1:nrow(outcome1_P2))
    
  {
    if(outcome1_P2[i,k]<=0.00001 & outcome1_P1[i,k]>0.00001 & outcome1_Y1[i,k]>0)
      
    {
      cc1<-c(cc1,"darkgreen") # x1 is winning
    }
    if(outcome1_P2[i,k]>0.00001 & outcome1_P1[i,k]<=0.00001 & outcome1_Y1[i,k]>0)
    {
      cc1<-c(cc1,"palegreen2") # x2 winning
    }
    if(outcome1_P2[i,k]>0.00001 & outcome1_P1[i,k]>0.00001 & outcome1_Y1[i,k]>0)
    {
      cc1<-c(cc1,"orange1") # coexistence
    }
    if((outcome1_P2[i,k]<=0.00001 & outcome1_P1[i,k]<=0.00001) | (outcome1_Y1[i,k]<0))
    {cc1<-c(cc1,"lightgrey")} #crash
 
    if(outcome2_P2[i,k]<=0.00001 & outcome2_P1[i,k]>0.00001 & outcome2_Y1[i,k]>0)
    {
      cc2<-c(cc2,"darkgreen") # x1 is winning
    }
    if(outcome2_P2[i,k]>0.00001 & outcome2_P1[i,k]<=0.00001 & outcome2_Y1[i,k]>0)
    {
      cc2<-c(cc2,"palegreen2") # x2 winning
    }
    if(outcome2_P2[i,k]>0.00001 & outcome2_P1[i,k]>0.00001 & outcome2_Y1[i,k]>0)
    {
      cc2<-c(cc2,"orange1") # coexistence
    }
    if((outcome2_P2[i,k]<=0.00001 & outcome2_P1[i,k]<=0.00001) | (outcome2_Y1[i,k]<0))
    {cc2<-c(cc2,"lightgrey")} #crash
    
    if(outcome3_P2[i,k]<=0.00001 & outcome3_P1[i,k]>0.00001 & outcome3_Y1[i,k]>0)
    {
      cc3<-c(cc3,"darkgreen") # x1 is winning
    }
    if(outcome3_P2[i,k]>0.00001 & outcome3_P1[i,k]<=0.00001 & outcome3_Y1[i,k]>0)
    {
      cc3<-c(cc3,"palegreen2") # x2 winning
    }
    if(outcome3_P2[i,k]>0.00001 & outcome3_P1[i,k]>0.00001 & outcome3_Y1[i,k]>0)
    {
      cc3<-c(cc3,"orange1") # coexistence
    }
    if((outcome3_P2[i,k]<=0.00001 & outcome3_P1[i,k]<=0.00001) | (outcome3_Y1[i,k]<0))
    {cc3<-c(cc3,"lightgrey")} #crash
    
     }
  cc_1<-rbind(cc_1,cc1)
  cc_3<-rbind(cc_3,cc3)
  cc_2<-rbind(cc_2,cc2)
}

#plot result at the border: l=1 you can change this location from 1 to 50.
l<-1
parameters_ord<-metaMDS(parameters[cc_1[l,]!="lightgrey",])
plot(parameters_ord,type="n")
points(parameters_ord,col=cc_1[l,][cc_1[l,]!="lightgrey"],pch=20)
ordiellipse(parameters_ord,cc_1[l,][cc_1[l,]!="lightgrey"],col=cc_1[l,][cc_1[l,]!="lightgrey"],draw="polygon")
text(parameters_ord, display = c("species"), labels=dimnames(parameters)[[2]])

parameters_ord<-metaMDS(parameters[cc_2[l,]!="lightgrey",])
plot(parameters_ord,type="n")
points(parameters_ord,col=cc_2[l,][cc_2[l,]!="lightgrey"],pch=20)
ordiellipse(parameters_ord,cc_2[l,][cc_2[l,]!="lightgrey"],col=cc_2[l,][cc_2[l,]!="lightgrey"],draw="polygon")
text(parameters_ord, display = c("species"), labels=dimnames(parameters)[[2]])

parameters_ord<-metaMDS(parameters[cc_3[l,]!="lightgrey",])
plot(parameters_ord,type="n")
points(parameters_ord,col=cc_3[l,][cc_3[l,]!="lightgrey"],pch=20)
ordiellipse(parameters_ord,cc_3[l,][cc_3[l,]!="lightgrey"],col=cc_3[l,][cc_3[l,]!="lightgrey"],draw="polygon")
text(parameters_ord, display = c("species"), labels=dimnames(parameters)[[2]])


cc<-cc_0

statistics_coexistence_nodiffusion<-cbind(parameters[cc=="orange1",],Pars_save[cc=="orange1",-18])
statistics_P1_nodiffusion<-cbind(parameters[cc=="darkgreen",],Pars_save[cc=="darkgreen",-18])
statistics_P2_nodiffusion<-cbind(parameters[cc=="palegreen2",],Pars_save[cc=="palegreen2",-18])

modindex<-c(cc[cc=="orange1"],cc[cc=="darkgreen"],cc[cc=="palegreen2"])
x_train<-rbind(statistics_coexistence_nodiffusion,statistics_P1_nodiffusion,statistics_P2_nodiffusion)

shuf<-sample(c(1:nrow(x_train)),1,replace=FALSE)
data.statistics<-t(x_train[shuf,])
y_test<-modindex[shuf]

simulations<-cbind(modindex[-shuf],x_train[-shuf,])
simulations<-data.frame(simulations)
simulations[,1]<-as.factor(simulations[,1])
colnames(simulations)<-"modindex"
colnames(simulations)[2:ncol(simulations)]<-c(dimnames(parameters)[[2]],dimnames(Pars_save)[[2]][-18])
for(i in 2:ncol(simulations))
{
  simulations[,i]<-as.numeric(simulations[,i])
}

model.rf1 <- abcrf(modindex~., data = simulations, ntree=250)
model.rf1
plot(model.rf1,simulations,obs=data.frame(data.statistics),col=simulations[,1])
predict(model.rf1, data.frame(data.statistics), simulations, ntree=250)

ppx<-as.matrix(simulations[,-1])%*%as.numeric(model.rf1$model.lda$scaling[,1])
ppy<-as.matrix(simulations[,-1])%*%as.numeric(model.rf1$model.lda$scaling[,2])
plot(ppx,ppy,pch=20,col=as.character(simulations[,1]),ylab="LD2",xlab="LD1")

disagree<-cc_0
disagree[disagree=="orange1"]<-2
disagree[disagree=="darkgreen"]<-1
disagree[disagree=="palegreen2"]<--1
disagree[disagree=="lightgrey"]<-NA
disagree<-as.numeric(disagree)

disagree1<-cc_1[1,]
disagree1[disagree1=="orange1"]<-2
disagree1[disagree1=="darkgreen"]<-1
disagree1[disagree1=="palegreen2"]<--1
disagree1[disagree1=="lightgrey"]<-NA
disagree1<-as.numeric(disagree1)

total<-disagree1-disagree # -3 = P1 win instead of coexistence, -2 = P1 wins instead of P2 , -1 = P2 wins instead of coexistence, 1 = Coexistence instead of P2, 2 = P2 wins instead of P1, 3 = coexistence instead of P1

simulations_disagree<-cbind(rep(1,1000),parameters,Pars_save[,-18])
simulations_disagree<-data.frame(simulations_disagree)
colnames(simulations_disagree)[1]<-"modindex"
simulations_disagree<-simulations_disagree[total!=0 & is.na(total)==FALSE,]
simulations_disagree[,1]<-as.factor(total[total!=0 & is.na(total)==FALSE])
model.rfdisagree <- abcrf(modindex~., data = simulations_disagree, ntree=250)
model.rfdisagree
plot(model.rfdisagree,simulations_disagree)

new<-as.numeric(as.character(simulations_disagree[,1]))
new[new%in%c(-3,-1)]<-"Disrupt"
new[new%in%c(1,3)]<-"Promote"
new[new==-2]<-"P1"
new[new==2]<-"P2"

simulations_disagree2<-simulations_disagree
simulations_disagree2[,1]<-as.factor(as.character(new))
model.rfdisagree2 <- abcrf(modindex~., data = simulations_disagree2, ntree=250)
model.rfdisagree2
plot(model.rfdisagree2,simulations_disagree2)


####
# P1
test<-rbind(Pars_save[cc_0=="darkgreen",-18],Pars_save[cc_1[1,]=="darkgreen",-18],Pars_save[cc_2[1,]=="darkgreen",-18],Pars_save[cc_3[1,]=="darkgreen",-18])
group<-c(rep("cc_0",length(cc_0[cc_0=="darkgreen"])),rep("cc_1",length(cc_1[1,][cc_1[1,]=="darkgreen"])),rep("cc_2",length(cc_2[1,][cc_2[1,]=="darkgreen"])),rep("cc_3",length(cc_3[1,][cc_3[1,]=="darkgreen"])))

test_rda<-rda(test~group)
anova(test_rda) #significant
#Model: rda(formula = test ~ group)
#Df Variance      F Pr(>F)   
#Model        3   0.0076 2.8555  0.005 **
#Residual 13242  11.7719  

test_ord<-metaMDS(test)
plot(test_ord,type="n",main="P1")
points(test_ord,pch=20,col=brewer.pal(11,"BrBG")[as.numeric(as.factor(group))])
ordiellipse(test_ord,group,col=brewer.pal(11,"BrBG")[1:4],draw="polygon",border=FALSE)
#text(test_ord, display = c("species"), labels=dimnames(test)[[2]])

test<-data.frame(test)
env<- cbind( (test$a11*test$a12)/(test$a21*test$a22),(test$beta11*test$beta12)/(test$beta21*test$beta22),test$alpha1/test$alpha2)
tt<-envfit(test_ord, env)
tt
plot(tt,labels=c("Competition","Predation","Fitness"),col="Black",cex=2)

legend("topright",col=brewer.pal(11,"BrBG")[1:4],legend=c("No spillover","Slow spillover","Fast spillover","Turbulences"),pch=20)


ww<-table(cut(test$a21[group=="cc_0"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)),cut(test$a11[group=="cc_0"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)))
dimnames(ww)[[1]]<-seq(0.025,0.975,length.out=20)
dimnames(ww)[[2]]<-seq(0.025,0.975,length.out=20)
levelplot(ww,xlab="P1 Interspecific Competition",ylab="P1 Intraspecific Competition",main="No diffusion")
l<-levelplot(ww,xlab="P1 Interspecific Competition",ylab="P1 Intraspecific Competition",main="No diffusion")


ww<-table(cut(test$a21[group=="cc_2"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)),cut(test$a11[group=="cc_2"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)))
dimnames(ww)[[1]]<-seq(0.025,0.975,length.out=20)
dimnames(ww)[[2]]<-seq(0.025,0.975,length.out=20)
levelplot(ww,xlab="P1 Interspecific Competition",ylab="P1 Intraspecific Competition",main="Diffusion",at=l$panel.args.common$at)


ww<-table(cut(test$a21[group=="cc_3"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)),cut(test$a11[group=="cc_3"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)))
dimnames(ww)[[1]]<-seq(0.025,0.975,length.out=20)
dimnames(ww)[[2]]<-seq(0.025,0.975,length.out=20)
levelplot(ww,xlab="P1 Interspecific Competition",ylab="P1 Intraspecific Competition",main="Turbulence",at=l$panel.args.common$at)

# P2

test2<-rbind(Pars_save[cc_0=="palegreen2",-18],Pars_save[cc_1[1,]=="palegreen2",-18],Pars_save[cc_2[1,]=="palegreen2",-18],Pars_save[cc_3[1,]=="palegreen2",-18])
group2<-c(rep("cc_0",length(cc_0[cc_0=="palegreen2"])),rep("cc_1",length(cc_1[1,][cc_1[1,]=="palegreen2"])),rep("cc_2",length(cc_2[1,][cc_2[1,]=="palegreen2"])),rep("cc_3",length(cc_3[1,][cc_3[1,]=="palegreen2"])))


test_rda2<-rda(test2~group2)
anova(test_rda2) #significant
#Model: rda(formula = test2 ~ group2)
#Df Variance      F Pr(>F)   
#Model        3   0.0077 2.8314  0.008 **
#Residual 13043  11.7958 

test_ord2<-metaMDS(test2)
plot(test_ord2,type="n",main="P2")
points(test_ord2,pch=20,col=brewer.pal(11,"BrBG")[as.numeric(as.factor(group2))])
ordiellipse(test_ord2,group2,col=brewer.pal(11,"BrBG")[1:4],draw="polygon",border=FALSE)
#text(test_ord2, display = c("species"), labels=dimnames(test2)[[2]])

test2<-data.frame(test2)
env2<- cbind( (test2$a11*test2$a12)/(test2$a21*test2$a22),(test2$beta11*test2$beta12)/(test2$beta21*test2$beta22),test2$alpha1/test2$alpha2)
tt2<-envfit(test_ord2, env2)
tt2
plot(tt2,labels=c("Competition","Predation","Fitness"),col="Black",cex=2)

legend("topright",col=brewer.pal(11,"BrBG")[1:4],legend=c("No spillover","Slow spillover","Fast spillover","Turbulences"),pch=20)


ww<-table(cut(test2$a12[group2=="cc_0"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)),cut(test2$a22[group2=="cc_0"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)))
dimnames(ww)[[1]]<-seq(0.025,0.975,length.out=20)
dimnames(ww)[[2]]<-seq(0.025,0.975,length.out=20)
levelplot(ww,xlab="P2 Interspecific Competition",ylab="P2 Intraspecific Competition",main="No diffusion")
l<-levelplot(ww,xlab="P2 Interspecific Competition",ylab="P2 Intraspecific Competition",main="No diffusion")
  
ww<-table(cut(test2$a12[group2=="cc_2"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)),cut(test2$a22[group2=="cc_2"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)))
dimnames(ww)[[1]]<-seq(0.025,0.975,length.out=20)
dimnames(ww)[[2]]<-seq(0.025,0.975,length.out=20)
levelplot(ww,xlab="P2 Interspecific Competition",ylab="P2 Intraspecific Competition",main="Diffusion",at=l$panel.args.common$at)


ww<-table(cut(test2$a12[group2=="cc_3"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)),cut(test2$a22[group2=="cc_3"],breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1)))
dimnames(ww)[[1]]<-seq(0.025,0.975,length.out=20)
dimnames(ww)[[2]]<-seq(0.025,0.975,length.out=20)
levelplot(ww,xlab="P2 Interspecific Competition",ylab="P2 Intraspecific Competition",main="Turbulence",at=l$panel.args.common$at)



#Coexistence

test3<-rbind(Pars_save[cc_0=="orange1",-18],Pars_save[cc_1[1,]=="orange1",-18],Pars_save[cc_2[1,]=="orange1",-18],Pars_save[cc_3[1,]=="orange1",-18])
group3<-c(rep("cc_0",length(cc_0[cc_0=="orange1"])),rep("cc_1",length(cc_1[1,][cc_1[1,]=="orange1"])),rep("cc_2",length(cc_2[1,][cc_2[1,]=="orange1"])),rep("cc_3",length(cc_3[1,][cc_3[1,]=="orange1"])))

test_rda3<-rda(test3~group3)
anova(test_rda3) #significant
#Model: rda(formula = test3 ~ group3)
#Df Variance      F Pr(>F)    
#Model       3   0.4176 32.465  0.001 ***
#Residual 2787  11.9499     

test_ord3<-metaMDS(test3)
plot(test_ord3,type="n",main="Coexistence")
points(test_ord3,pch=20,col=brewer.pal(11,"BrBG")[c(1,4,7,10)][as.numeric(as.factor(group3))])
ordiellipse(test_ord3,group3,col=brewer.pal(11,"BrBG")[c(1,4,7,10)],draw="polygon",border=FALSE)
#text(test_ord3, display = c("species"), labels=dimnames(test3)[[2]])


test3<-data.frame(test3)
env3<- cbind( (test3$a11*test3$a12)/(test3$a21*test3$a22),(test3$beta11*test3$beta12)/(test3$beta21*test3$beta22),test3$alpha1/test3$alpha2)
tt3<-envfit(test_ord3, env3)
tt3
plot(tt3,labels=c("Competition","Predation","Fitness"),col="Black",cex=2)

legend("topright",col=brewer.pal(11,"BrBG")[c(1,4,7,10)],legend=c("No spillover","Spillover","Spillover with wind","Turbulences"),pch=20)

step_mean<-function(x){
  m<-c()
  for(i in 2:length(x))
  {m<-c(m,(x[i-1]+(x[i]-x[i-1])/2))}
    return(as.numeric(m))
  }

ww<-table(cut(test3$a21[group3=="cc_0"]/test3$a12[group3=="cc_0"],breaks= seq(0,2,length.out=21)),cut(test3$a11[group3=="cc_0"]/test3$a22[group3=="cc_0"],breaks=seq(0,2,length.out=21)))
dimnames(ww)[[1]]<-step_mean(seq(0,2,length.out=21))
dimnames(ww)[[2]]<-step_mean(seq(0,2,length.out=21))
levelplot(ww,xlab="P1/P2 Interspecific Competition",ylab="P1/P2 Intraspecific Competition",main="No diffusion")
l<-levelplot(ww,xlab="P1/P2 Interspecific Competition",ylab="P1/P2 Intraspecific Competition",main="No diffusion")

ww<-table(cut(test3$a21[group3=="cc_2"]/test3$a12[group3=="cc_2"],breaks= seq(0,2,length.out=21)),cut(test3$a11[group3=="cc_2"]/test3$a22[group3=="cc_2"],breaks=seq(0,2,length.out=21)))
dimnames(ww)[[1]]<-step_mean(seq(0,2,length.out=21))
dimnames(ww)[[2]]<-step_mean(seq(0,2,length.out=21))
levelplot(ww,xlab="P1/P2 Interspecific Competition",ylab="P1/P2 Intraspecific Competition",main="Diffusion",at=l$panel.args.common$at)

ww<-table(cut(test3$a21[group3=="cc_3"]/test3$a12[group3=="cc_3"],breaks= seq(0,2,length.out=21)),cut(test3$a11[group3=="cc_3"]/test3$a22[group3=="cc_3"],breaks=seq(0,2,length.out=21)))
dimnames(ww)[[1]]<-step_mean(seq(0,2,length.out=21))
dimnames(ww)[[2]]<-step_mean(seq(0,2,length.out=21))
levelplot(ww,xlab="P1/P2 Interspecific Competition",ylab="P1/P2 Intraspecific Competition",main="Turbulence",at=l$panel.args.common$at)




