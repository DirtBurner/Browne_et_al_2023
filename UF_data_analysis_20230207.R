## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----load_packages---------------------------------------------------------------------------------------------------------------------------------
library(zoo)
library(MASS)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(nlstools)
library(ezknitr)



## ----user defined values---------------------------------------------------------------------------------------------------------------------------

core_name <- "LMG1311_JKC1"
collectyear=2013 #year in CE that core was collected
countyear=2016 #year in CE that core was counted for radioisotopes
Depth_M=c(0.4,1.2,2,2.8,3.6,5.2,6.8,8.5,10.5,13.5,18.5,23.5,33.5)
Pb210T=c(66.15303314,62.89575353,60.16993251,56.1229667,55.26586492,47.39009405,31.40510612,25.86412078,17.92583009,8.888739124,7.011996368,4.371368536,2.879387116)
Pb210u=c(0.923975511,1.043404646,1.043855765,1.077380321,0.871133708,0.899180449,0.629698332,0.731284593,0.554923413,0.379968221,0.413155963,0.313522917,0.248186141)
Ra226=c(6.821091269,5.016832664,5.003393066,5.380786404,4.610007074,5.353152765,5.295144045,7.390493346,5.153360507,5.420436271,4.578162811,4.411351507,4.888541848)
Ra226U=c(3.326330511,2.170383149,1.85628713,2.386004097,1.694741448,2.654557505,2.557715672,2.11977753,2.053059948,2.398919635,1.810716536,1.487010557,1.910846372)
DBD=c(0.254175,0.256045,0.26,0.274595,0.36452,0.36572,0.312005,0.28579,0.349225,0.372895,0.20087,0.1829,0.343255)


## ----Code initialization---------------------------------------------------------------------------------------------------------------------------
a=length(Pb210T)
nmeasures=a[1]
ageresults=matrix(0,nmeasures,4)
ageresults[,1]=Depth_M


## ----Establish radioisotopic values, include=FALSE-------------------------------------------------------------------------------------------------
Pb <- Pb210T #total Pb-210 activity (dpm/g)
Pbunc <- Pb210u #total Pb-210 activity uncertainty (dpm/g)
Ra <- Ra226 #total Ra-226 activity (dpm/g)
Raunc <- Ra226U #total Ra-226 activity uncertainty (dpm/g)
xsPb <- Pb-Ra
lambdaPb=log(2)/22.3
replace(xsPb, xsPb<0, 0.01) #replaces potential negative xsPb values with very small values
xsPbcor0=xsPb/(exp(-0.000085332*(countyear-collectyear))) #decay corrects activity to time of core collection
xsPbcor <- replace(xsPbcor0, xsPbcor0<0, 0.01) 


## ----regression depth range------------------------------------------------------------------------------------------------------------------------
# Establish depth range for regression
ztop=2 # number of first sample at base of surface mixed layer; MAKE SURE THIS NUMBER IS ACCURATE FOR EACH PROFILE

A0=xsPbcor[ztop]#first sample at base of surface mixed layer
AZ=xsPbcor/A0#normalized activity
zbot <- nmeasures



## ----Calculate Cumulative Mass---------------------------------------------------------------------------------------------------------------------
# Calculate Cumulative Mass
cmval <- matrix(0,nmeasures)
DBDval <- DBD
depth <- Depth_M
dbdmean <- rollmean(DBDval,2)
for (v in 2:nmeasures ) {
  cmval[v]  <-  ((depth[v]-depth[v-1])*dbdmean[v-1])+cmval[v-1]
}
DBDave <- mean(DBDval[3:nmeasures])# average of dry bulk density values excluding top two samples; adjust as necessary to average only the vertically stable part of the profile



## ----Calculate nls---------------------------------------------------------------------------------------------------------------------------------

activitydata <- data.frame(AZ[ztop:zbot], cmval[ztop:zbot])
colnames(activitydata) <- c("AZ", "cmval")
modeltop <-ztop# number of first sample at base of surface mixed layer; MAKE SURE THIS NUMBER IS ACCURATE FOR EACH PROFILE
modelbot <- 11 #sample number for base of excess activity; avoids using any data points that are at supported levels of Pb-210 activity
modeldata <- activitydata[modeltop:modelbot,]
expo.der <- deriv3(~exp(-LS * cmval),
                   c("LS"),
                   function(cmval, LS) NULL)
str(expo.der)

LS<- lambdaPb/0.1
start <- list(LS=LS)

expmod <- nls(AZ~expo.der(cmval, LS), data=modeldata, start=start)
summary(expmod)
confint(expmod)

pred <- data.frame(cmval=seq(cmval[modeltop],cmval[modelbot],l=25))
der <- do.call(expo.der, args=c(list(cmval=pred$cmval), as.list(coef(expmod))))

F <- attr(der, "gradient") 
U <- chol(vcov(expmod))
se <- sqrt(apply(F%*%t(U), 1, function(x) sum(x^2))) #standard error on regression


plot(AZ~cmval, data=modeldata,xlab="Cumulative mass depth in core (g/cm2)",
     ylab="Normalized Excess Pb-210 Activity (dpm/g)",
     #xlim=c(0,cmval[nmeasures]), ylim=c(0,1))
    xlim=c(0,10), ylim=c(0,1))
matlines(pred$cmval, c(der)+
           outer(se, qt(c(.5, .025,.975), df=df.residual(expmod))),
         type="l", col=c(1,2,2), lty=c(1,2,2))
legend("topright",
       legend=c("observed values", "predicted values",
                "confidence interval (95%)"),
       lty=c(NA,1,2), col=c(1,1,2), pch=c(1,NA,NA), bty="n")




## ----calculate ages, include=FALSE-----------------------------------------------------------------------------------------------------------------
exp_s <- summary(expmod)
exp_c <- confint(expmod)#95% confidence intervals on regression
massratemid <- lambdaPb/exp_s$coefficients[1]# midpoint mass acc. rate, g/cm2/yr
expmmodfitci <- confint(expmod, level=0.95)
massratelow <- lambdaPb/expmmodfitci[1]# lower mass acc. rate, g/cm2/yr
massratehigh <- lambdaPb/expmmodfitci[2]# upper mass acc. rate, g/cm2/yr
fitpval <- exp_s$coefficients[4]
fitpvalr <- round(fitpval, digits = 5)


## ----establish range of sed rates------------------------------------------------------------------------------------------------------------------
#establish range of sed rates
sedratemid <- massratemid/DBDave #95% CI linear sedimentation rate, data-point estimate (cm/y)
sedratelow <- massratelow/DBDave #95% CI linear sedimentation rate, low estimate (cm/y)
sedratehigh <- massratehigh/DBDave #95% CI linear sedimentation rate, high estimate (cm/y)
ageresults[,2] <- round(collectyear-ageresults[,1]/sedratemid)
ageresults[,3] <- round(collectyear-ageresults[,1]/sedratelow)
ageresults[,4] <- round(collectyear-ageresults[,1]/sedratehigh)
agemodel_out <- as.data.frame(ageresults)
colnames(agemodel_out) <- c("DIC", "agemid","agehigh","agelow")


## ----plot age-depth data---------------------------------------------------------------------------------------------------------------------------
# plot age-depth data
agedepth <- ageresults[,1]
agemid <- ageresults[,2]
agelow <- ageresults[,4]
agehigh <- ageresults[,3]
ageplot <- data.frame(agedepth, agemid, agelow, agehigh)
ageplotmlt <- melt(ageplot,id.vars="agedepth",measure.vars = c("agemid","agehigh", "agelow"))

#depthmax <-ageplot[nmeasures,1]
depthmax <-20
agemax <- 1900
p1 <- ggplot(ageplotmlt, aes(x = agedepth, y = value, group = variable)) +
  geom_line(aes(colour=variable, group=variable), size=1) +
  scale_color_manual(values=c("#000000","#CC6666","#0072B2"))+
  scale_x_continuous(breaks = seq(0, depthmax, 5), limits = c(0, depthmax)) +
  scale_y_continuous(breaks = seq(agemax, 2020, 20), limits = c(agemax, 2020)) +
  labs(x = "Depth in Core (cm)", y = "Age (Common Era)") + 
  ggtitle(core_name) +
  theme(legend.position = c(0.08, 0.35))
  theme_bw() +
  theme(panel.background = element_rect(fill = 'grey95'))+
  theme(panel.grid.major = element_line(colour = "black", size=0.5),
        panel.grid.minor = element_line(colour = "grey75"),
        plot.title = element_text(size = rel(1.1), face = "bold"))

p1



agemodel_out


## ----save plot-------------------------------------------------------------------------------------------------------------------------------------
ggsave(paste0(core_name,".pdf"), 
       p1, width = 6, height = 5)

write.csv(agemodel_out, file=(paste0(core_name,"modelresults.csv")))

