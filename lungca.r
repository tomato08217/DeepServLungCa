setwd("D:/杨斌/977/0609")
train<-read.csv("train_clinical.csv")
test<-read.csv("test_clinical.csv")
clinical<-read.csv("976_clinical.csv")
fvars<-c("Sex","Smoking.status","Family.History","TTF_1","Ki_67","Histologic.subtype","Stage","T","N","Chemotherapy","ICI","EGFR.gene.mutation","Lymphatic.metastasis")
clinical[fvars] <- lapply(clinical[fvars],factor)

library(tableone)
vars<-c("Sex","Age","Smoking.status","Family.History","TTF_1","Ki_67","Histologic.subtype","Stage","T","N","Chemotherapy","ICI","EGFR.gene.mutation","Lymphatic.metastasis","risk_score","PLT","M.1","N.1","L","CPR","LMR","NLR")
nonnormal<-c("risk_score","Age","PLT","M.1","N.1","L","CPR","LMR","NLR")
tabletwo<-CreateTableOne(vars=vars,strata=c("Cohort"),data=clinical)
table2<-print(tabletwo,nonnormal=nonnormal,cramVars=fVars,showAllLevels=T)
write.csv(table2,file="table1.csv")

clinical$Sex<-as.factor(clinical$Sex)
clinical$Age<-as.numeric(as.character(clinical$Age))
clinical$Smoking.status<-as.factor(clinical$Smoking.status)
clinical$Grade.of.Differentiation<-as.factor(clinical$Grade.of.Differentiation)
clinical$FamilyHistory<-as.factor(clinical$FamilyHistory)
clinical$Histologic.subtype<-as.factor(clinical$Histologic.subtype)
clinical$T<-as.factor(clinical$T)
clinical$N<-as.factor(clinical$N)
clinical$M<-as.factor(clinical$M)
clinical$Stage<-as.factor(clinical$Stage)
clinical$Chemotherapy<-as.factor(clinical$Chemotherapy)
clinical$EGFR.gene.mutation<-as.factor(clinical$EGFR.gene.mutation)
clinical$TTF_1<-as.factor(clinical$TTF_1)
clinical$Ki_67<-as.factor(clinical$Ki_67)
y_clinical<- Surv(clinical$time,clinical$Event) #bundle y
library(survival)
a<-coxph(y_clinical~Age,data=clinical)
summary(a)
b<-coxph(y_clinical~Sex,data=clinical)
summary(b)
c<-coxph(y_clinical~Smoking.status,data=clinical)
summary(c)
d<-coxph(y_clinical~Family.History,data=clinical)
summary(d)
e<-coxph(y_clinical~TTF_1,data=clinical)
summary(e)
f<-coxph(y_clinical~Ki_67,data=clinical)
summary(f)
g<-coxph(y_clinical~Histologic.subtype,data=clinical)
summary(g)
h<-coxph(y_clinical~Grade.of.Differentiation,data=clinical)
summary(h)
i<-coxph(y_clinical~Stage,data=clinical)
summary(i)
j<-coxph(y_clinical~T,data=clinical)
summary(j)
k<-coxph(y_clinical~N,data=clinical)
summary(k)
l<-coxph(y_clinical~M,data=clinical)
summary(l)
m<-coxph(y_clinical~ICI,data=clinical)
summary(m)
n<-coxph(y_clinical~PLT,data=clinical)
summary(n)
o<-coxph(y_clinical~N.1,data=clinical)
summary(o)
p<-coxph(y_clinical~M.1,data=clinical)
summary(p)
q<-coxph(y_clinical~L,data=clinical)
summary(q)
r<-coxph(y_clinical~CPR,data=clinical)
summary(r)
s<-coxph(y_clinical~NLR,data=clinical)
summary(s)
t<-coxph(y_clinical~LMR,data=clinical)
summary(t)
u<-coxph(y_clinical~EGFR.gene.mutation,data=clinical)
summary(u)

train$Sex<-as.factor(train$Sex)
train$FamilyHistory<-as.factor(train$FamilyHistory)
train$Histologic.subtype<-as.factor(train$Histologic.subtype)
train$Stage<-as.factor(train$Stage)
train$Chemotherapy<-as.factor(train$Chemotherapy)
train$Target<-as.factor(train$Target)
train$EGFR.gene.mutation<-as.factor(train$EGFR.gene.mutation)
train$Lymphatic.metastasis<-as.factor(train$Lymphatic.metastasis)
train$TTF_1<-as.factor(train$TTF_1)
train$Ki_67<-as.factor(train$Ki_67)
train$T<-as.factor(train$T)
train$N<-as.factor(train$N)
train$Smoking.status<-as.factor(train$Smoking.status)

test$Sex<-as.factor(test$Sex)
test$FamilyHistory<-as.factor(test$FamilyHistory)
test$Histologic.subtype<-as.factor(test$Histologic.subtype)
test$Stage<-as.factor(test$Stage)
test$Chemotherapy<-as.factor(test$Chemotherapy)
test$Target<-as.factor(test$Target)
test$EGFR.gene.mutation<-as.factor(test$EGFR.gene.mutation)
test$Lymphatic.metastasis<-as.factor(test$Lymphatic.metastasis)
test$TTF_1<-as.factor(test$TTF_1)
test$Ki_67<-as.factor(test$Ki_67)
test$T<-as.factor(test$T)
test$N<-as.factor(test$N)


library(survival)

y_train<- Surv(train$time,train$Event) #bundle y
y_test<-Surv(test$time,test$Event)

#cox 
Rmodel <- coxph(y_train~risk_score, data = train)
Cmodel1<-coxph(y_train ~ Age + Sex + Smoking.status +  Histologic.subtype + Stage  + TTF_1 + Ki_67+T+N+N.1+M.1+L+LMR,  data = train)
TNMmodel <- coxph(y_train~T+N, data = train)
summary(Cmodel1)
#AIC backward
Cmodel<- step(Cmodel1, direction = "backward")
#or
Cmodel<-coxph(y_train~Age+TTF_1+Stage+Ki_67+Histologic.subtype+M.1,train)
CRmodel<-coxph(y_train~risk_score+Age+TTF_1+Stage+Ki_67,train)


#test with 95%CI

method<- survConcordance(Surv(test$time,test$Event) ~ predict(Rmodel, test))
cat(method$concordance, "95%CI:", method$concordance-method$std.err*1.96, method$concordance+method$std.err*1.96)

method<- survConcordance(Surv(test$time,test$Event) ~ predict(Cmodel, test))
cat(method$concordance, "95%CI:", method$concordance-method$std.err*1.96, method$concordance+method$std.err*1.96)

method<- survConcordance(Surv(test$time,test$Event) ~ predict(CRmodel, test))
cat(method$concordance, "95%CI:", method$concordance-method$std.err*1.96, method$concordance+method$std.err*1.96)

#forest map
library(survminer)
ggforest(CRmodel,train)

#nomogram
library(rms)
library(openxlsx)
train_model<-read.xlsx("model.xlsx",sheet=1)
test_model<-read.xlsx("model.xlsx",sheet=2)
dd<-datadist(train_model)
options(datadist="dd")
#build with cph
CRmodel2<-cph(y_train~risk_score+Age+TTF_1+Stage+Ki_67,x=T,y=T,surv=T,train_model,time.inc=100)#define predicion time point
surv<- Survival(CRmodel2)
surv1<-function(x)surv(36,lp=x)
surv2<-function(x)surv(60,lp=x)
nom<-nomogram(CRmodel2, fun=list(surv1,surv2),funlabel=c("3 year surivival","5 year survival"),lp=F,fun.at=c('0.90','0.70','0.5','0.3','0.1'),maxscale=10)
plot(nom,xfrac = 0.6)
#calculate total points
library(nomogramFormula)
results <- formula_rd( nom )
points_train<- points_cal(formula = results$formula,rd=train_model )
points_test<- points_cal(formula = results$formula,rd=test_model )
#calculate probablities
points2<-prob_cal(reg = CRmodel2,times = c(36,60))

#calibration curve,build with cph,x=T,y=T
library(rms)
CRmodel36<-cph(y_train~risk_score+Age+TTF_1+Stage+Ki_67,train,x=T,y=T,surv=TRUE,time.inc=36)
cal36<- calibrate(CRmodel36,cmethod='KM', method='boot', u=36, m=100,B=50)
CRmodel60<-cph(y_train~risk_score+Age+TTF_1+Stage+Ki_67,train,x=T,y=T,surv=TRUE,time.inc=60)
cal60<- calibrate(CRmodel60,cmethod='KM', method='boot', u=60, m=100,B=50)#u=time.inc,m=how many groups
plot(cal36,lwd = 2,lty = 1,errbar.col = c("blue"),xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-predicted probability (%)",ylab = "Observed probability (%)",col = c("blue"),cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal36[,c('mean.predicted',"KM")],type= 'b',lwd = 2,col = c("blue"),pch = 16)
plot(cal60,lwd = 2,lty = 1,errbar.col = c("orange"),xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-predicted probability (%)",ylab = "Observed probability (%)",col = c("blue"),cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal60[,c('mean.predicted',"KM")],type= 'b',lwd = 2,col = c("orange"),pch = 16)

#decision curve
library(devtools) 
library(rmda)
crmodelforDCA<-decision_curve(Event~risk_score+Age+TTF_1+Stage+Ki_67,data =train,thresholds = seq(0,1, by = .01),bootstraps = 10)
rmodelforDCA<-decision_curve(Event~risk_score,data=train,thresholds = seq(0,1, by = .01),bootstraps = 10)
cmodelforDCA<-decision_curve(Event~Age+TTF_1+Stage+Ki_67,data=train,thresholds = seq(0,1, by = .01),bootstraps = 10)
tnmmodelforDCA<-decision_curve(Event~T+N,data=train,thresholds = seq(0,1, by = .01),bootstraps = 10)
plot_decision_curve( list(cmodelforDCA,rmodelforDCA,crmodelforDCA,tnmmodelforDCA), curve.names = c("cmodel","rmodel", "crmodel","tnmmodel"),col = c("blue","black","red","orange"), cost.benefit.axis = FALSE,confidence.intervals =FALSE,legend.position = "none")
legend("bottomright", inset=.05, , c("clinical model","deep learning model","clinical plus deep learning model","TNM model"),lty=c(1),cex=0.3, col=c("blue","black", "red","orange"))

#NRI IRI
library(survIDINRI)
train_model<-read.xlsx("model.xlsx",sheet=1)
test_model<-read.xlsx("model.xlsx",sheet=2)
t0=36
t1=60
outcome=train_model[,c(1,2)]
covs1<-as.matrix(train_model[,c(-1:-2)])
covs0<-as.matrix(train_model[,c(-1:-3)])
x1<-IDI.INF(outcome, covs0, covs1, t0, npert=100)
x2<-IDI.INF(outcome, covs0, covs1, t1, npert=100)
IDI.INF.OUT(x1) 
IDI.INF.OUT(x2) 
#m1 Result of IDI. Point and corresponding (1-alpha/2) confidence interval are given
#m2 Result of continuous-NRI. Point and corresponding (1-alpha/2) confidence interval are given.
#m3 Result of median improvement in risk score

#compare CINDEX
library(survcomp)
cindex_cr <- concordance.index(predict(CRmodel),surv.time = train$time,  surv.event = train$Event,method = "noether") 
cindex_c <- concordance.index(predict(Cmodel),surv.time = train$time,  surv.event = train$Event,method = "noether") 
cindex_r<-concordance.index(predict(Rmodel),surv.time = train$time,  surv.event = train$Event,method = "noether") 
cindex_tnm<-concordance.index(predict(TNMmodel),surv.time = train$time,  surv.event = train$Event,method = "noether") 

cindex.comp(cindex_cr, cindex_tnm)

#ph test
cox.zph(CRmodel)
library(survminer)
p<-cox.zph(CRmodel)
ggcoxzph(p)


# Integrated Brier score and prediction error curve
library(pec)
library(riskRegression)
library(rms)
library(Hmisc)
library(survival)
library(survminer)
library(survcomp)
dd=datadist(train)
options(datadist="dd")
Models <- list("CRmodel"= coxph(Surv(time,Event)~risk_score+Age+TTF_1+Stage+Ki_67+Histologic.subtype+M, data=train,x=TRUE,y=TRUE),
               "Cmodel" = coxph(Surv(time,Event)~ Age+TTF_1+Stage+Ki_67+Histologic.subtype+M, data=train,x=TRUE,y=TRUE),
               "Rmodel" = coxph(Surv(time,Event)~ risk_score, data=train,x=TRUE,y=TRUE),
               "TNMmodel"=coxph(Surv(time,Event)~T+N,data=train,x=T,y=T))
p <- pec(object = Models,cens.model = "cox", data=train, splitMethod="Boot632plus", B=100,reference = FALSE)
print(p)
par(mai=c(1,1,1,1))
plot(p,type="l",smooth=TRUE,legend = FALSE,xlim=c(0,24),axis1.at=seq(0,24,2), xlab="Survival weeks", ylab="Prediction error",col = c("red", "blue","black"),lwd = c(3,3,3),lty = c(1,1,1))



#KM
library(survival)
library(survminer)
res.cut <- surv_cutpoint(train, time = "time", event = "Event",variables = c("risk_score", "Age"))
res.cat <- surv_categorize(res.cut)
risk_score_fitforKM <- survfit(Surv(time,Event) ~risk_score, data = res.cat)
age_fitforKM <- survfit(Surv(time,Event) ~Age, data = res.cat)
ggsurvplot(risk_score_fitforKM , data = res.cat, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom",legend.title = "risk",legend.labs = c("Deep Learning high risk","Deep Learning low risk"))
ggsurvplot(age_fitforKM , data = res.cat, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom",legend.title = "risk",legend.labs = c("older than 67", "younger than 67"))

Stage_fitforKM <- survfit(Surv(time,Event)~ Stage, data = train)
ggsurvplot(Stage_fitforKM , data = train, risk.table = TRUE,pval=T,conf.int = FALSE,xlab='Time in Months',legend = "top",legend.title = "risk")
Ki67_fitforKM <- survfit(Surv(time,Event)~ Ki_67, data = train)
ggsurvplot(Ki67_fitforKM , data = train, risk.table = TRUE,pval=T,conf.int = FALSE,xlab='Time in Months',legend = "top",legend.title = "risk",legend.labs = c("Ki67-", "Ki67+"))
TTF_fitforKM <- survfit(Surv(time,Event)~TTF_1, data = train)
ggsurvplot(TTF_fitforKM , data = train, risk.table = TRUE,pval=T,conf.int = FALSE,xlab='Time in Months',legend = "top",legend.title = "risk",legend.labs = c("TTF-1-", "TTF-1+"))

#median

library("sqldf")
risk_score_highrisk <- subset(clinical,risk_score>median(clinical$risk_score))
risk_score_lowrisk <- subset(clinical,risk_score<median(clinical$risk_score))
riskscorehigh_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = risk_score_highrisk)
riskscorelow_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = risk_score_lowrisk)
ggsurvplot(riskscorehigh_chemo_fitforKM, data = risk_score_highrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")
ggsurvplot(riskscorelow_chemo_fitforKM, data = risk_score_lowrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")

age_highrisk <- subset(clinical,Age>median(clinical$Age))
age_lowrisk <- subset(clinical,Age<median(clinical$Age))
agehigh_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = age_highrisk)
agelow_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data =age_lowrisk)
ggsurvplot(agehigh_chemo_fitforKM, data = age_highrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")
ggsurvplot(agelow_chemo_fitforKM, data = age_lowrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")

TTF_highrisk <- subset(clinical,TTF_1==0)
TTF_lowrisk <- subset(clinical,TTF_1==1)
TTFhigh_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = TTF_highrisk)
TTFlow_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data =TTF_lowrisk)
ggsurvplot(TTFhigh_chemo_fitforKM, data = TTF_highrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")
ggsurvplot(TTFlow_chemo_fitforKM, data = TTF_lowrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")

Ki67_highrisk <- subset(clinical,Ki_67==1)
Ki67_lowrisk <- subset(clinical,Ki_67==0)
Ki67high_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = Ki67_highrisk)
Ki67low_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data =Ki67_lowrisk)
ggsurvplot(Ki67high_chemo_fitforKM, data = Ki67_highrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")
ggsurvplot(Ki67low_chemo_fitforKM, data = Ki67_lowrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")

stage_highrisk <- subset(clinical,Stage>2)
stage_lowrisk <- subset(clinical,Stage<3)
stage_high_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = stage_highrisk)
stage_low_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data =stage_lowrisk)
ggsurvplot(stage_high_chemo_fitforKM, data = stage_highrisk , risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")
ggsurvplot(stage_low_chemo_fitforKM, data = stage_lowrisk, risk.table = TRUE, pval=T,conf.int = F,xlab='Time in Months',legend = "bottom")

#risk score/new score risk stratification
risk.group<-ifelse(clinical$risk_score>median(clinical$risk_score),'high','low')
forKM_riskscore_rate<-survfit(Surv(time,Event)~ risk.group, data =clinical)
ggsurvplot(forKM_riskscore_rate, data = clinical, risk.table = TRUE,pval=T,pval.method=TRUE,conf.int = FALSE,xlab='Time in Months',legend = "bottom",legend.title = "risk",legend.labs = c("Deep learning high risk","Deep learning low risk")) 
newscore.group<-ifelse(clinical$new_score>median(clinical$new_score),'high','low')
forKM_newscore_rate<-survfit(Surv(time,Event)~ newscore.group, data =clinical)
ggsurvplot(forKM_newscore_rate, data = clinical, risk.table = TRUE,pval=T,pval.method=TRUE,conf.int = FALSE,xlab='Time in Months',legend = "bottom",legend.title = "risk",legend.labs = c("New score high risk","New score low risk")) 

#newscore
clinical<-read.csv("976_clinical.csv")
library("sqldf")
low1 <- subset(clinical, new_score_cutoff3==0)
middle1 <- subset(clinical, new_score_cutoff3==1)
high1 <- subset(clinical, new_score_cutoff3==2)

low2<-subset(clinical, new_score_median_7.2==0)
high2<-subset(clinical, new_score_median_7.2==1)

low3<-subset(clinical, new_score_cutoff==0)
high3<-subset(clinical, new_score_cutoff==1)

low4 <- subset(clinical, new_score_Q3==0)
middle4 <- subset(clinical, new_score_Q3==1)
high4 <- subset(clinical, new_score_Q3==2)

high1_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = high1)
ggsurvplot(high1_chemo_fitforKM, data = high1, risk.table = TRUE, pval=T,conf.int = F)
middle1_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = middle1)
ggsurvplot(middle1_chemo_fitforKM, data = middle1, risk.table = TRUE, pval=T,conf.int = F) 
low1_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = low1)
ggsurvplot(low1_chemo_fitforKM, data = low1, risk.table = TRUE, pval=T,conf.int = F)  
high2_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = high2)
ggsurvplot(high2_chemo_fitforKM, data = high2, risk.table = TRUE, pval=T,conf.int = F)
low2_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = low2)
ggsurvplot(low2_chemo_fitforKM, data = low2, risk.table = TRUE, pval=T,conf.int = F) 
high3_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = high3)
ggsurvplot(high3_chemo_fitforKM, data = high3, risk.table = TRUE, pval=T,conf.int = F)
low3_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = low3)
ggsurvplot(low3_chemo_fitforKM, data = low3, risk.table = TRUE, pval=T,conf.int = F)  
high4_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = high4)
ggsurvplot(high4_chemo_fitforKM, data = high4, risk.table = TRUE, pval=T,conf.int = F)
middle4_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = middle4)
ggsurvplot(middle4_chemo_fitforKM, data = middle4, risk.table = TRUE, pval=T,conf.int = F) 
low4_chemo_fitforKM <- survfit(Surv(time,Event) ~Chemotherapy, data = low4)
ggsurvplot(low4_chemo_fitforKM, data = low4, risk.table = TRUE, pval=T,conf.int = F)  


#pearson
radiomics<-read.csv("976_radiomics_1198.csv")
y<-radiomics[,5:1204]
c<-0
d<-0
for (i in 1:1200) c[i]<-cor.test(y[,i],y$NLR)$estimate
for (i in 1:1200) d[i]<-cor.test(y[,i],y$LMR)$estimate
pearson <- read.csv("pearson.csv",sep=",",header=T)
x<-pearson[,2:3]
rownames(x)<-pearson$X
x = t(x)
pheatmap(x,cluster_rows = FALSE, cluster_cols = FALSE,main = "Features correlation",show_colnames = F,show_rownames = T,annotation_colors = ann_colors,border_color = "NA")
