#Mathias 07-Jan-14
#Estimating seed dispersal kernels for the megafauna
#=================================================

#=================
#1. seed retention 
#================

#Gamma distribution for seed retention time - Guttal et al. 2011
mean=50 #Chosen according to seed retention time of large mammals
var=400 #Chosen to reproduce a similar range
k=(mean^2)/var
theta=mean/var #Wrong in Guttal et al. 2011

a=rgamma(10000,shape=k,rate=theta)


#======================
#2. Simulating movement
#======================

require(adehabitat)

R=100 #Replicates
steps=250 #number of movement bouts
#I'm assuming we register positions each hour (250 hours)


dist=rep(NA,R)
dist.final=rep(NA,R)
cumdist=matrix(NA,steps,R)
st.dist=matrix(NA,steps,R)

for (i in 1:R){ #Simulating R levy walks
	
	w <- simm.levy(1:steps, mu = 3,l0=50, burst = "mu = 3") #Simulating levy walk

	mov=w[[1]]

	dist[i]=sum(mov$dist,na.rm=TRUE) #total travel distance
	cumdist[,i]=cumsum(mov$dist) #cumulative travel distance per timestep
	dist.final[i]=sqrt(mov[steps,1]^2+mov[steps,2]^2) # farthest point reached after all timesteps
	st.dist[,i]=sqrt(mov[,1]^2+mov[,2]^2) # distance (straight line) of each point relative to starting point

}

#======================
#3. seed dispersal Kernel
#======================
source('~/Documents/Colaboracoes/M.Galetti - Oxford/Data/code/mykernel.R', chdir = TRUE)

#3.1. Substituing the time axis
time=sort(round(a)) #discrete hours
dist.time=st.dist[time,] #distance each seed was dispersed
mykernel(dist.time, bw=40, h= 5)


f=splinefun(density(dist.time)) #to compute the probability density associated with any given value

sum(f(0:5000)) #Computes the probability for a given interval (integrating) - from 0 to maximum should yield ~1

#======================
#4. Total dispersal kernel
#======================

#Gamma distribution for seed retention time - Guttal et al. 2011
mean=50 #Chosen according to seed retention time of large mammals
var=400 #Chosen to reproduce a similar range
k=(mean^2)/var
theta=mean/var #Wrong in Guttal et al. 2011

#-----------------------------------
#4.1.seed retention for three species
k1=k
k2=k*1.2
k3=k*0.8

d1=rgamma(1000,shape=k1,rate=theta)
d2=rgamma(1000,shape=k2,rate=theta)
d3=rgamma(1000,shape=k3,rate=theta)


#-----------------------------------
#4.2.simulating levy walks

R=100 #Replicates
steps=250 #number of movement bouts
#I'm assuming we register positions each hour (250 hours)


st.dist1=matrix(NA,steps,R)
st.dist2=st.dist1
st.dist3=st.dist1

require(adehabitat)

for (i in 1:R){ #Simulating R levy walks
	w1 <- simm.levy(1:steps, mu = 3,l0=50, burst = "mu = 3") #Simulating levy walk
	w2 <- simm.levy(1:steps, mu = 3,l0=50, burst = "mu = 3") #Simulating levy walk
	w3 <- simm.levy(1:steps, mu = 3,l0=20, burst = "mu = 3") #Simulating levy walk
	
	mov1=w1[[1]]
	mov2=w2[[1]]
	mov3=w3[[1]]
	
	st.dist1[,i]=sqrt(mov1[,1]^2+mov3[,2]^2) # distances at each timestep 
	st.dist2[,i]=sqrt(mov2[,1]^2+mov2[,2]^2)
	st.dist3[,i]=sqrt(mov3[,1]^2+mov1[,2]^2)
}

#-----------------------------------
#4.3.Combining retention time and movement

#discrete hours
time1=sort(round(d1))
time2=sort(round(d2)) 
time3=sort(round(d3))  

#distance each seed was dispersed
dist.time1=st.dist1[time1,]
dist.time2=st.dist2[time2,]
dist.time3=st.dist3[time3,]

totd=c(dist.time1,dist.time2,dist.time3)

#=============================
#5. Plotting
#============================
#setwd("~/Documents/Colaboracoes/M.Galetti - Oxford")
#pdf("DKernel.pdf")

#-----------------------------------
# 5.1 Total kernel (histogram and rug)
require(MASS)
truehist(totd, xlim= c(0,6000),
         # ylim=c(0,0.012),
         prob= T, h= 5, xlab= "Distance (m)",
         ylab= "Probability", col= rgb(0.5, 0.5, 0.5, 0.3),
         lty= 0)
rug(totd, side= 1,col= "tomato")

#------------------------------
# 5.2 Species-specific kernels
dens1<-density(dist.time1, bw= 50,from= 0,to= max(totd)) # add density estimate
lines(dens1,xlim= c(0, 10000), col=rgb(0,104,139,max=255), lwd=2)

dens2<-density(dist.time2, bw= 50, from= 0, to= max(totd)) # add density estimate
lines(dens2, xlim= c(0, 10000), col=rgb(205,102,29,200,max=255),lwd=2)

dens3<-density(dist.time3, bw= 50, from= 0, to= max(totd)) # add density estimate

lines(dens3, xlim=c(0, 10000),col=rgb(107,142,35,max=255),lwd=2)

dens.tot<-density(totd, bw= 50,from= 0,to= max(totd))
lines(dens.tot,xlim= c(0, 10000), col= "black", lwd=3)

#------------------
#5.3 plotting histograms
par(new=TRUE)
par(oma=c(6,4,2,2))
#par(mar=c(1,2,1,2))
par(mfcol=c(4,2), mfg=c(1,2))

#histograms for each species
#par(mfrow=c(3,1), lwd=1.3, cex.axis=1.2, cex.lab=1.3)

par(mar=c(1,2,4,5))
hist(dist.time1, xlim=c(0,6000),col=rgb(0,104,139,200,max=255),breaks=20,xlab="",main="",ylab="", freq=TRUE)

par(mar=c(2,2,2,5))
hist(dist.time2, xlim=c(0,6000), col=rgb(205,102,29,200,max=255),breaks=20, xlab="", main="",freq=TRUE)
mtext("frequency", 2,cex=0.8,padj=-3)

par(mar=c(4,2,1,5))
hist(dist.time3, xlim=c(0,6000), col=rgb(107,142,35,200,max=255), xlab="Distance(m)",main="",ylab="",freq=TRUE,breaks=12)
#dev.off()

#computing probabilities
# f=splinefun(dens.tot) #to compute the probability density associated with any given value
# sum(f(1000:3000))



#===================================
#5.4 comparing Kernels
#setwd("~/Documents/Colaboracoes/M.Galetti - Oxford")
#pdf("DKernel2.pdf")

#----------------------
#5.4.1 with Pleistocene megafauna
# Total kernel (histogram and rug)

par(mfrow=c(2,1))
require(MASS)
truehist(totd, xlim= c(0,6000),
          ylim=c(0,0.0020),
         prob= T, h= 5, xlab= "Distance (m)",
         ylab= "Probability", col= rgb(0.5, 0.5, 0.5, 0.3),
         lty= 0)
rug(totd, side= 1,col= "tomato")

# Species-specific kernels
dens1<-density(dist.time1, bw= 50,from= 0,to= max(totd)) # add density estimate
lines(dens1,xlim= c(0, 10000), col=rgb(0,104,139,max=255), lwd=2)

dens2<-density(dist.time2, bw= 50, from= 0, to= max(totd)) # add density estimate
lines(dens2, xlim= c(0, 10000), col=rgb(205,102,29,200,max=255),lwd=2)

dens3<-density(dist.time3, bw= 50, from= 0, to= max(totd)) # add density estimate
lines(dens3, xlim=c(0, 10000),col=rgb(107,142,35,max=255),lwd=2)

dens.tot<-density(totd, bw= 50,from= 0,to= max(totd))
lines(dens.tot,xlim= c(0, 10000), col= "black", lwd=3)

#------------------------
#computing probabilities
f=splinefun(dens.tot) #to compute the probability density associated with any given value
P.tot=round(sum(f(1000:6000)),2)
text(x=5000,y=.001,substitute(paste(italic(P)[1000], " = ",a),list(a=P.tot)))


#----------------------
#5.4.2 after extinctions
truehist(dist.time3, xlim= c(0,6000),
          ylim=c(0,0.0020),
         prob= T, h= 5, xlab= "Distance (m)",
         ylab= "Probability", col= rgb(0.5, 0.5, 0.5, 0.3),
         lty= 0)
rug(dist.time3, side= 1,col= "tomato")



dens3<-density(dist.time3, bw= 50, from= 0, to= max(totd)) # add density estimate
lines(dens3, xlim=c(0, 10000),col=rgb(107,142,35,max=255),lwd=2)

#------------------------
#computing probabilities
f=splinefun(dens3) #to compute the probability density associated with any given value
P.def=round(sum(f(1000:6000)),2)

text(x=5000,y=.001,substitute(paste(italic(P)[1000], " = ",a),list(a=P.def)))

#dev.off()

