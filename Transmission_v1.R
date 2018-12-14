###
### Alison C. Ketz 12/10/2018
###
### Transmission following Heisey et al (2010) rcode
### for fitting FOI model given the cross-sectional harvest data
###

###
### Preliminaries
###

rm(list=ls())

# setwd("C:/Users/aketz/Documents/Transmission/Transmission_v1")
setwd("~/Documents/Transmission/Transmission_v1")

library(coda)
# library(Hmisc)
library(lubridate)
library(xlsx)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(xtable)
library(nimble)
library(dplyr)
library(tidyr)
library(lattice)
library(beepr)

###
### Load previous runs
###

#load()

###
### Load data and clean
###
#
# #Load CWD aging data
# age.cwd.df = read.xlsx("C:/Users/aketz/Documents/Data/Harvest/CWDagingIowaDaneGrant_01-24-2018.xlsx",sheetName="CWDagingIowaDaneGrant_01-24-201")
# names(age.cwd.df)
#
# table(age.cwd.df$age)
# cwd.df = read.csv("C:/Users/aketz/Documents/Data/WDNR_surveillance/Surveillance_data.csv")
cwd.df = read.csv("~/Documents/Data/WDNR_surveillance/Surveillance_data.csv")
levels(cwd.df$age) = c("","1","9+","9+","2","3","4-5","6-8","ADULT","0")
cwd.df$kill_date=as.Date(cwd.df$kill_date)
cwd.df=cwd.df[cwd.df$age!="ADULT",]
cwd.df$age=as.factor(as.character(cwd.df$age))

#proportion of tested deer without results
# table(cwd.df$has_results)[1]/sum(table(cwd.df$has_results)) #.0061

#removing cases where there were no cwd results
cwd.df=cwd.df[cwd.df$has_results=="Y",]

#proportion of cases without sex
# table(cwd.df$sex)[1]/dim(cwd.df)[1]#=.004

#removing all cases w/o sex
cwd.df=cwd.df[cwd.df$sex!="",]

#checking the stats of individuals without ages assigned
# plot(cwd.df$kill_date[cwd.df$age==""])

###assessing the missing age data from the surveillance dataset
#missing age data from surveillance data by year
# hist(year(cwd.df$kill_date[cwd.df$age==""]))
# table(cwd.df$age)[1]/dim(cwd.df)[1]#.011  = 1 percent of data are not aged
# table(cwd.df$positive[cwd.df$age==""])[2]/length(cwd.df$positive[cwd.df$age==""])#.0176 raw positive w/o age classification
# table(cwd.df$positive)[2]/length(cwd.df$positive)#.031 overall raw positive

cwd.df=cwd.df[cwd.df$age!="",]
cwd.df$age=as.factor(cwd.df$age)
cwd.df$positive = as.numeric(cwd.df$positive)-1

#removing the individuals with kill date that are NA
cwd.df = cwd.df[!is.na(cwd.df$kill_date),]

### Removing data where no section nor quarter section was recorded
# head(cwd.df[is.na(cwd.df$sect),],20)
# plot(cwd.df$kill_date[is.na(cwd.df$sect)]) #more than half of data missing space was in 2002, early....
# hist(cwd.df$sect)
# sum(is.na(cwd.df$sect))/dim(cwd.df)[1] #.0187

#removing NAs of section level
cwd.df = cwd.df[!is.na(cwd.df$sect),]
# sum(is.na(cwd.df$quar_sect))/dim(cwd.df)[1] #10% of quarter sections not recorded (although section level available)

###################################################################################################################################
#redoing Table 1 from Dennis's Ecological Monographs (2010) paper
###################################################################################################################################

apparent.prev = cwd.df %>% group_by(.dots=c("sex","age")) %>% summarise(n=n(),positive = sum(positive), prevalence = sum(positive)/n())
# apparent.prev = apparent.prev[c(7,1,3,4,5,6,2,14,8,10,11,12,13,9),]
apparent.prev
###################################################################################################################################

#convert sex to numeric integer values, Females = 1, Males = 0
levels(cwd.df$sex)=c(NA,1,0)


### convert age class to numeric
cwd.df$age=as.factor(as.character(cwd.df$age))
levels(cwd.df$age)=c(0,1,2,3,4,7,10)
cwd.df$age.num = as.numeric(as.character(cwd.df$age))

#calculate approximate birthdate based on age at kill_date
ageclass=as.numeric(levels(as.factor(cwd.df$age.num)))

#add birth_date to dataframe and calculate values for it
# assuming birth date is May 15, depending on kill_date and ageclass at kill_date
cwd.df$birth_date = NA
for(j in 1:dim(cwd.df)[1]){
  if(format(cwd.df$kill_date[j], format="%m-%d")<="05-15"){
    cwd.df$birth_date[j] = paste((year(cwd.df$kill_date[j])-cwd.df$age.num[j]-1),"-05-15",sep="")
  }else{
    cwd.df$birth_date[j] = paste(year(cwd.df$kill_date[j])-cwd.df$age.num[j],"-05-15",sep="")
  }
}
cwd.df$birth_date=as.Date(cwd.df$birth_date)

#calculating approximate age in days and months
cwd.df$agedays = cwd.df$kill_date - cwd.df$birth_date
cwd.df$agemonths = as.integer(ceiling(cwd.df$agedays/(365.25/12)))

# cwd.df[which(cwd.df$agemonths==min(cwd.df$agemonths)),]



#eda of ages in months
cwd.df$age.num[cwd.df$age.num==0]=0.5

hist(cwd.df$age.num)
boxplot(cwd.df$age.num)
plot(cwd.df$age.num)
summary(cwd.df$age.num)

###
### Format data for running in the model
###

cwd.df = cwd.df %>% arrange(cwd.df$sect,cwd.df$age,cwd.df$kill_date)

#constants for model
N.sect = max(cwd.df$sect) #number of spatial units = 36 sections
runsum.dat = c(0,cumsum(as.numeric(table(cwd.df$sect))))
agem.dat = cwd.df$agemonths

class(cwd.df$agemonths)



max(runsum.dat)
agem.dat[dim(cwd.df)[1]]
agem.dat[runsum.dat[N.sect+1]]
max(agem.dat)



###
### simulation for indexing matrices gamma and dayhaz
###

beta0 = -6
gam =  matrix(NA,nr=dim(cwd.df)[1],nc=max(agem.dat))
dayhaz =  matrix(NA,nr=dim(cwd.df)[1],nc=max(agem.dat))

for (i in 1:N.sect) {
  for (j in (runsum.dat[i]+1):runsum.dat[i+1]) {
    for (k in 1:agem.dat[j]) {
      gam[j,k]<-(beta0)
      dayhaz[j,k] <- exp(gam[j,k])
    } #k
  } #j
} #i

head(dayhaz)

dayhaz.in=dayhaz
gam.in = gam
gam.in[!is.na(gam)]=NA
gam.in[is.na(gam)]=0
dayhaz.in[!is.na(dayhaz)]=NA
dayhaz.in[is.na(dayhaz)]=0




##################################################################################################################################
###
### Model code
###
##################################################################################################################################
#Data variable definitions:
#N = number of cells.
#runsum[i] = the count of all subjects in cells with indices less than i.
#runsum[i]+1 to runsum[i+1] are all subjects in cell i.
#agem[j] = the age of subject j in months.
#female[j] = indicates whether subject j is a female (1) or male (0).
#Tage = oldest age bin.
#j is rows of individs
#k is age column
##################################################################################################################################

#Table 2 - overall intercept, nonspatial
modelcode = nimbleCode({
    beta0 ~ dflat()
    for (i in 1:N) {
        for (j in (runsum[i]+1):runsum[i+1]) {
            for (k in 1:agem[j]) {
                gamma[j,k]<-(beta0)
                dayhaz[j,k] <- exp(gamma[j,k])
            } #k
            icumhaz[j]<-sum(dayhaz[j,1:agem[j]])
            p[j] <-1-exp(-icumhaz[j])
            pos[j] ~ dbern(p[j])
        } #j
    } #i
})#end model


nimData = list(pos=cwd.df$positive,
               dayhaz = dayhaz.in,
               gamma = gam.in,
               icumhaz=rep(NA,dim(cwd.df)[1]),
               p = rep(NA,dim(cwd.df)[1])
               )

nimConsts = list(runsum=runsum.dat,
                 N = N.sect,
                 agem=cwd.df$agemonths)
#female
names(cwd.df)

nimInits = list(beta0 = rnorm(1,-6,1))


modelout<-nimbleModel(code= modelcode,
                      name="Transmission_Int_nospace",
                      constants = nimConsts,
                      data = nimData,
                      inits = nimInits,
                      check=FALSE,calculate=FALSE)



#number of MCMCr iterations, Chains, and Burn-in
reps<-10000
bin<-0.5*reps
n.chains<-1

##One-line MCMC call
starttime<-Sys.time()
mcmcout<-nimbleMCMC(model=modelout,
                    nchains=n.chains,
                    nburnin = bin,
                    niter=reps,
                    monitors=c("beta0"))
Sys.time()-starttime
beep(sound=5)

# out<-mcmc.list(lapply(mcmcout$samples, mcmc))
out = mcmc(mcmcout)
fit.sum=summary(mcmcout)

pdf("figures/param_plot_v6.pdf")
densityplot(out[,"tau"],ylab="tau")
traceplot(out[,"tau"],ylab="tau")
densityplot(out[,"beta0"],ylab="beta0")
traceplot(out[,"beta0"],ylab="beta0")
dev.off()

##################################################################################################################################
###
### Spatial Models
###
##################################################################################################################################

#Data variable definitions:
#N = number of cells.
#runsum[i] = the count of all subjects in cells with indices less than i.
#runsum[i]+1 to runsum[i+1] are all subjects in cell i.
#agem[j] = the age of subject j in months.
#female[j] = indicates whether subject j is a female (1) or male (0).
#Tage = oldest age bin.
#j is rows of individs
#k is age column


nT = max(cwd.df$agemonths)

#create num vector
num.yr1 = c(1,rep(2,nT-2),1)

#create adjacency vector yr1
temp<-as.matrix(bandSparse(n=nT,k=c(1),symmetric = T))
temp2<-matrix(0,nT,nT)
for(i in 1:nrow(temp2)){
  temp2[i,] = temp[i,]*c(1:nT)
}
adj = t(temp2)[which(t(temp2)!=0)]






#Table 2 - overall intercept, spatial
modelcode = nimbleCode({

    hprec ~ dgamma(0.5,0.0005)
    sprec ~ dgamma(0.5,0.0005)

    hprec1<-0.2*hprec

    for(j in 1:N){
        hetero[j]~dnorm(0,hprec1)
    }

    sdspace<-sd(space[1:N])
    sdhetero<-sd(hetero[1:N])
    psi<-sdspace/(sdhetero+sdspace)

### weighting for the car normal spatial effect - will have to figure this one out.
for(k in 1:sumNumNeigh){
    sweights[k] <- 1
}
    space[1:N] ~ dcar_normal(sadj[1:N], sweights[1:N], snum[1:N], sprec)
    beta0 ~ dflat()

    for (i in 1:N) {
        for (j in (runsum[i]+1):runsum[i+1]) {
            for (k in 1:agem[j]) {
                gamma[j,k]<-(space[i]+hetero[i]+beta0)
                dayhaz[j,k] <- exp(gamma[j,k])
            } #k
            icumhaz[j]<-sum(dayhaz[j,1:agem[j]])
            p[j] <-1-exp(-icumhaz[j])
            pos[j] ~ dbern(p[j])
        } #j #individual
    } #i #spatial location i
})#end model



nimData = list(pos=cwd.df$positive)

nimConsts = list(runsum=runsum.dat,
                 N = N.sect,
                 agem=cwd.df$agemonths)
#female


nimInits = list(beta0 = runif(1,0,1))





modelout<-nimbleModel(code= modelcode,
                      name="Transmission_Int_Space",
                      constants = nimConsts,
                      data = nimData,
                      inits = nimInits)



#number of MCMCr iterations, Chains, and Burn-in
reps<-50000
bin<-0.5*reps
n.chains<-1
