library(chron) 
library(lubridate)
library(scales)
library(dplyr) 
library(tidyr)
library(rlist)
library(stringr)
library(tidyverse)
library(anytime)
library(reshape2)
library(zoo)
library(magrittr)
library(reshape)
library(arsenal)
library("dplyr", character.only = TRUE)
library(data.table)
library(imputeTS)
library(ggplot2)
library(sfsmisc)
library(reshape)
library(e1071)
library(splines)
options(scipen=999)
library(nlme)
library(segmented)
library(tictoc)
library(lmtest)
library(investr)

All_Data_BUS<-read.csv(file="C:\\Users\\gaslani\\OneDrive - Massey University\\FirstProject\\NewTasks\\LocationTwo\\BT_BUS.csv")
All_Data_BUS<-All_Data_BUS %>% select(-c(X,Bus.prop))
colnames(All_Data_BUS)<-c("BT", "Bus.count","ATC","day","time")

# Adding another variable as time of day and Type of day ( weekdays vs weekend)
f =as.POSIXct(All_Data_BUS$time,tz="Europe/London", origin="1970-01-01")
Hour<-hour(f)
mon<-month(f)

All_Data_BUS=All_Data_BUS  %>% mutate(Hour = Hour)%>% mutate(Month= mon) %>% mutate(DayWeek = weekdays(as.Date(day)))

test<-All_Data_BUS %>%
  mutate(TypeOfDay = case_when(DayWeek=="Monday" |DayWeek=="Tuesday"| DayWeek=="Wednesday" | DayWeek=="Thursday"| DayWeek=="Friday" ~ "Weekday", 
                               DayWeek=="Saturday" |DayWeek=="Sunday"  ~ "Weekend"))%>%
  mutate(TypeOfHour = case_when(7<=Hour & Hour<10  |15<=Hour & Hour<19 ~ "Busy-time", 
                                TRUE  ~ "Non-Busy-time")) %>%
  mutate(season = case_when(Month==3 |Month==4| Month==5  ~ "Spring", 
                            Month==6 |Month==7| Month==8  ~ "Summer", 
                            Month==9 |Month==10| Month==11  ~ "Autumn",
                            Month==12 |Month==1| Month==2 ~ "Winter")) 


test$TypeOfDay<-as.factor(test$TypeOfDay)
test$Hour<-as.factor(test$Hour)
test$TypeOfHour<-as.factor(test$TypeOfHour)
test$Month<-as.factor(test$Month)
test$season<-as.factor(test$season)
#test$season<-	factor(c("Spring", "Summer","Autumn","Winter"), levels=c("Spring", "Summer","Autumn","Winter"))
test$DayWeek<-as.factor(test$DayWeek)
#test$DayWeek<-factor(c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"), levels=c("Monday", "Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday" ))

#*****************************************************
#%%%%%%%%%%% Type of Hour %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#=================================================================================
#*# @@@@@@@@@@@@@@@@@@@@@@@@@@@ ALL Days seperately @@@@@@@@@@@@@@@@@@@@@@@@
dd<-rep(1:7, by=1) # different days
dh<-rep(1:24, by=1) # different hours
DAY<-c("Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday")




#*****************************************************
#*****************************************************
#*#*****************************************************
#*****************************************************v

# Linear Mixed Effect Models Part
lme.test<-test %>% select(c(BT,ATC,Hour,DayWeek))
lme.test$DayWeek<-ordered(lme.test$DayWeek, levels = c("Monday", "Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday" ))

# Range for x0=[l,u]
l<-20  
u<-120
yzero=30
#*******************************************************************
#*#***************************************************************

# ****************** Profile log-likelihood via optim() **********************
#***************************************************************
#*#***************************************************************
Profile.LOG=function(par,data){
  t<-as.integer(c(30,par))
  t1<-c(dh[10],DAY[1])
  extra.row.lme<-rbind(data,c(t,t1))
  extra.row.lme$BT<-as.integer(extra.row.lme$BT)
  extra.row.lme$ATC<-as.integer(extra.row.lme$ATC)
  regoutput.lme<-lme(BT~ATC, data =extra.row.lme, random = ~ ATC | DayWeek/Hour,control = lmeControl(opt = "optim") )
  return(logLik(regoutput.lme))
}


optim.res<-optim(par =c(70), fn =  Profile.LOG, data = lme.test, method =c("Brent") ,
                 lower = l, upper = u,control=list(fnscale=-1))


f <- function(x,maxloglik){
  t<-as.integer(c(30,x))
  t1<-c(dh[10],DAY[1])
  extra.row.lme<-rbind(lme.test,c(t,t1))
  extra.row.lme$BT<-as.integer(extra.row.lme$BT)
  extra.row.lme$ATC<-as.integer(extra.row.lme$ATC)
  mm=lme(BT~ATC, data =extra.row.lme,
         random = ~ ATC | DayWeek/Hour,control = lmeControl(opt = "optim"))
  logLik(mm) - maxloglik + qchisq(.80,1)/2 
}

# Computing CI 
r1.opt<-uniroot(f,c(20,optim.res$par),maxloglik=optim.res$value)
r2.opt<-uniroot(f,c(optim.res$par,120),maxloglik=optim.res$value)


#*******************************************************************
#*#***************************************************************

# ****************** Profile log-likelihood via loop **********************
#***************************************************************
#*#***************************************************************

vec.lme=seq.int(l,u, 1)
res.lme=rep(0,length(vec.lme))
tic("model fitting")
for(i in 1:length(vec.lme)){
  t<-as.integer(c(yzero,vec.lme[i]))
  t1<-c(dh[10],DAY[1])
  extra.row.lme<-rbind(lme.test,c(t,t1))
  extra.row.lme$BT<-as.integer(extra.row.lme$BT)
  extra.row.lme$ATC<-as.integer(extra.row.lme$ATC)
  regoutput.lme<-lme(BT~ATC, data =extra.row.lme,
                     random = ~ ATC | DayWeek/Hour,
                     control = lmeControl(opt = "optim"))
 
  res.lme[i]<-(logLik(regoutput.lme))
}
xxx<-toc()

result<-data.frame(Estimate.ATC=numeric(),LogLike=numeric())
result<-as.data.frame(cbind(vec.lme,res.lme))

max_point<-result[which(result$res.lme==max(result$res.lme)),]

# Computing CI 
r1<-uniroot(f,c(l,max_point[1,1]),maxloglik=max_point[1,2])
r2<-uniroot(f,c(max_point[1,1],200),maxloglik=max_point[1,2])



ggplot(data=result, aes(x=vec.lme,y=res.lme)) +
  geom_point()+
  labs(x="X",y="logLikelihood")+
  geom_vline(aes(xintercept=round(max_point[1,1],1)), size=1.5, color="red")+
  geom_hline(yintercept=(max_point[1,2]-0.8211872 ))+
  geom_vline(aes(xintercept=r1$root), size=1,linetype='dotted', color="blue")+
  geom_vline(aes(xintercept=r2$root), size=1,linetype='dotted', color="blue")+
  annotate(x=max_point[1,1],y=+Inf,label=paste0("ATC=",round(max_point[1,1],1)),vjust=2,geom="label")+
  annotate(x=round(r1$root,1),y=+Inf,label=paste0("low=",round(r1$root,1)),vjust=2,geom="label")+
  annotate(x=round(r2$root,1),y=+Inf,label=paste0("up=",round(r2$root,1)),vjust=2,geom="label")+
  ggtitle(paste("Profile LogLikelihood \n Assuming Hour & Day-",dh[10],DAY[1]))#+










#*#****************************************






