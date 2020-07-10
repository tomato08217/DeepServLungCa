library(Hmisc)
library(survival)
setwd("D:/杨斌/977/0609")
train<-read.csv("train_clinical.csv")
#read the data
d<-train

#we decided we are interested in calibration at 36 mo and 60 mo,
#so censor after 36 and 60
#already done in excel


#calculate predicted probability at pre-specified time (adm.cens)
survcox_d<-coxph(data=train, Surv(GND36_time,GND36_status)~risk_score+Age+TTF_1+Stage+Ki_67)
summary(survcox_d)
survfit_d=survfit(survcox_d, newdata=d, se.fit=FALSE)
survcox_e<-coxph(data=train, Surv(GND60_time,GND60_status)~risk_score+Age+TTF_1+Stage+Ki_67)
summary(survcox_e)
survfit_e=survfit(survcox_e, newdata=d, se.fit=FALSE)

survpr36=survfit_d$surv[215,]
estsurv36=survpr36
estinc36=1-survpr36

survpr60=survfit_e$surv[312,]
estsurv60=survpr60
estinc60=1-survpr60

#split into deciles
d$dec36=as.numeric(cut2(estinc36, g=10))
d$dec60=as.numeric(cut2(estinc60, g=10))

#check that there are 5 or more events in each group
#if not then collapse groups
table(d$dec36, d$GND36_status)
table(d$dec60, d$GND60_status)

source("GND_test.v2.r")
#calculate the GND test
GND.result1=GND.calib(pred=estinc36, tvar=d$GND36_time, out=d$GND36_status, 
                     cens.t=36, groups=d$dec36, adm.cens=36)
GND.result2=GND.calib(pred=estinc60, tvar=d$GND60_time, out=d$GND60_status, 
                     cens.t=60, groups=d$dec60, adm.cens=60)

GND.result1
GND.result2
