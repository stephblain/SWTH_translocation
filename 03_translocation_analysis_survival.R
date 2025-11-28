#######################################################
# Set up packages and data ############################
#######################################################

library(tidyverse)
library(Matrix)
library(marked)
library(NatParksPalettes)
library(ggpubr)

setwd("C:/Users/Steph/GitHub_data/translocation/")

thrush1<-read.csv("survival.translocatedBirds.ch_10days.20250914.csv")
pemberton1<-read.csv("survival.PembertonBirds.2020_2023.ch_10days.20250914.csv")

today<-gsub("-","",Sys.Date())

thrush<-thrush1%>%
  mutate(ch=gsub("days_","",ch))%>%
  #mutate(sex_f=as.factor(paste("f",sex_binary,sep="")))%>%
  mutate(release_site=as.factor(release_site))%>%
  mutate(ancestry=hiest_ancestry*(-1)+1)%>%
  select(ch,ancestry,release_site)
table(thrush$release_site)

#######################################################
# Run CJS analysis ####################################
#######################################################


model1.parameters=list(Phi=list(formula=~release_site*ancestry),
                      p=list(formula=~time+release_site+ancestry))

model1="CJS"

cjsfit=crm(thrush,model=model1,hessian=T,
            model.parameters=model1.parameters)

sink(paste(model1,"beta.translocation",today,"txt",sep="."))
cjsfit$results
sink(file = NULL)

#Only significant effect is the elevated detection probability for Squamish vs Lilloet


fit1<-predict(cjsfit)
phi1<-fit1$Phi
phi1$est.adj<-phi1$estimate^29
p1<-fit1$p

write.csv(phi1,paste(model1,"fit.phi.translocation",today,"csv",sep="."),row.names=F)
write.csv(p1,paste(model1,"fit.p.translocation",today,"csv",sep="."),row.names=F)

theme_set(theme_classic())

ggplot(phi1,aes(x=est.adj,fill=release_site))+
  geom_histogram(colour="white")+
  scale_fill_manual(values=natparks.pals("Banff")[c(4,1)])

ggplot(phi1,aes(x=release_site,y=ancestry,colour=est.adj))+
  geom_jitter(size=4,shape=19,alpha=0.8,width=0.1)+
  viridis::scale_colour_viridis(option="plasma",end=0.9)
  

#Fit a CJS model with detections binned into 10 day chunks to create 30 time points
#No evidence for the predicted interaction between release site and ancestry
#None of release site, ancestry, or the interaction between release site and 
#ancestry had a significant effect on survival


#######################################################
# Add pemberton into the analysis #####################
#######################################################

all(colnames(thrush1)==colnames(pemberton1))

thrush<-rbind(thrush1,pemberton1)%>%
  mutate(ch=gsub("days_","",ch))%>%
  mutate(release_site=factor(release_site,levels=c("Pemberton","Lilloet","Squamish")))%>%
  mutate(release_year=as.factor(release_year))%>%
  mutate(ancestry=hiest_ancestry*(-1)+1)%>%
  select(ch,ancestry,release_site,release_year)

model1.parameters=list(Phi=list(formula=~release_site*ancestry),
                       p=list(formula=~time+release_site+release_year+ancestry))


model1="CJS"

cjsfit=crm(thrush,model=model1,hessian=T,
           model.parameters=model1.parameters)

sink(paste(model1,"beta.translocation_Pemberton",today,"txt",sep="."))
cjsfit$results
sink(file = NULL)

#Only significant effect is the elevated detection probability for Squamish vs Lilloet

fit1<-predict(cjsfit)
phi1<-fit1$Phi
phi1$est.adj<-phi1$estimate^29
p1<-fit1$p

write.csv(phi1,paste(model1,"fit.phi.translocation_Pemberton",today,"csv",sep="."),row.names=F)
write.csv(p1,paste(model1,"fit.p.translocation_Pemberton",today,"csv",sep="."),row.names=F)


#######################################################
# Add pemberton into the analysis #####################
#######################################################

pemberton<-thrush1%>%filter(release_site=="Pemberton")

model1.parameters=list(Phi=list(formula=~release_year*ancestry),
                       p=list(formula=~time+release_year+ancestry))


model1="CJS"

cjsfit=crm(thrush,model=model1,hessian=T,
           model.parameters=model1.parameters)

#sink(paste(model1,"beta.translocation_PembertonOnly",today,"txt",sep="."))
cjsfit$results
#sink(file = NULL)