library(xlsx)
library(readxl)
library(gemtc)
library(dplyr)
library(igraph)


setwd("~/OneDrive/NMA")


antibiotics <- data.frame(read_xlsx("Data_example.xlsx",
                                    na=c("","NA","999","."),
                                    skip=0),
                          stringsAsFactors = F)
colnames(antibiotics)[colnames(antibiotics)=="X__1"] <- "trt"

d1<-unique(antibiotics$IntPrim)
temp1<-data.frame(matrix(,nrow = 27,ncol = 2))
temp1$X1<-d1
for (i in 1:27){
  temp1[i,2]<-LETTERS[i]
}

temp1[27,2]<-"AA"


antibiotics <- merge(antibiotics, temp1, by.x=c("IntPrim"), by.y=c("X1"),all = T)
colnames(antibiotics)[colnames(antibiotics)=="X2"] <- "Class"


bacteremia<-antibiotics[,c("ID","BacPtNum","BacPtDen","Class","trt","Author","IntPrim")] #quitar intprim
bacteremia<-na.omit(bacteremia)
length(unique(bacteremia$ID))
bacteremia<-bacteremia %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
bacteremia<-bacteremia[!bacteremia$ID %in% names(which(table(bacteremia$ID)<2)),]
colnames(bacteremia)[colnames(bacteremia)=="ID"] <- "study"
colnames(bacteremia)[colnames(bacteremia)=="BacPtNum"] <- "responders"
colnames(bacteremia)[colnames(bacteremia)=="BacPtDen"] <- "sampleSize"
d<-unique(bacteremia$trt)
temp_bac<-data.frame(matrix(,nrow = 21,ncol = 2))
temp_bac$X1<-d
for (i in 1:21){
  temp_bac[i,2]<-letters[i]
}
bacteremia <- merge(bacteremia, temp_bac, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(bacteremia)[colnames(bacteremia)=="X2"] <- "treatment"


bacteremia_class<-bacteremia %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
bacteremia_class<-na.omit(bacteremia_class)

bacteremia<-bacteremia %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
bacteremia<-na.omit(bacteremia)
#####

bacteremia_class %>% 
  group_by(study) %>% 
  summarise(n = n()) %>% filter(n==3 | n==4)

bacteremia_class %>% 
  group_by(study) %>% filter(study==18 | study==23 | study==84 | study==27202)

bacteremia_class[21:23,]
bacteremia_class[22,3]<-bacteremia_class[22,3]+bacteremia_class[21,3]
bacteremia_class[22,2]<-bacteremia_class[22,2]+bacteremia_class[21,2]
bacteremia_class<-bacteremia_class[-21,]
bacteremia_class[27:29,]
bacteremia_class[28,3]<-bacteremia_class[28,3]+bacteremia_class[27,3]
bacteremia_class[28,2]<-bacteremia_class[28,2]+bacteremia_class[27,2]
bacteremia_class<-bacteremia_class[-27,]
bacteremia_class[85:87,]
bacteremia_class[86,3]<-bacteremia_class[86,3]+bacteremia_class[87,3]
bacteremia_class[86,2]<-bacteremia_class[86,2]+bacteremia_class[87,2]
bacteremia_class<-bacteremia_class[-87,]
bacteremia_class[111:113,]
bacteremia_class[111,3]<-bacteremia_class[111,3]+bacteremia_class[112,3]
bacteremia_class[111,2]<-bacteremia_class[111,2]+bacteremia_class[112,2]
bacteremia_class<-bacteremia_class[-112,]

######


bacteremia_class[11:13,]
bacteremia_class[12,3]<-bacteremia_class[12,3]+bacteremia_class[11,3]
bacteremia_class[12,4]<-bacteremia_class[12,4]+bacteremia_class[11,4]
bacteremia_class<-bacteremia_class[-11,]
bacteremia_class[17:19,]
bacteremia_class[18,3]<-bacteremia_class[18,3]+bacteremia_class[17,3]
bacteremia_class[18,4]<-bacteremia_class[18,4]+bacteremia_class[17,4]
bacteremia_class<-bacteremia_class[-17,]
bacteremia_class[57:59,]
bacteremia_class[58,3]<-bacteremia_class[58,3]+bacteremia_class[59,3]
bacteremia_class[58,4]<-bacteremia_class[58,4]+bacteremia_class[59,4]
bacteremia_class<-bacteremia_class[-59,]
bacteremia_class<-bacteremia_class[,1:5]
length(unique(bacteremia_class$study))

bacteremia<-bacteremia[,c(1:4,7)]

colnames(bacteremia_class)[colnames(bacteremia_class)=="Class"] <- "treatment"
network_bacteremia_class<-mtc.network(as.data.frame(bacteremia_class))
network_bacteremia<-mtc.network(as.data.frame(bacteremia))
plot(network_bacteremia_class,vertex.color="coral")
plot(network_bacteremia,vertex.color="coral")
legend("topleft", legend=unique(paste(bacteremia$treatment,bacteremia$trt,sep = " - ")),
       cex=.8,bty = "n",text.font = 2)


legend("topleft", legend=c("B - No antibiotic", "C - Fluoroquinolone and Non-Anti-Pseudomonal Beta Lactam","E - Fluoroquinolone", "F - Placebo","H - Cephalosporin",
                           "J - Fluoroquinolone and Macrolide","K - Fluoroquinolone and Rifamycin","L - TMP-SMX","N - Quinolone_not_fluoro"),
       cex=.8,bty = "n",text.font = 2)

modelfixed_bacteremia_cl<-mtc.model(network_bacteremia_class,linearModel="fixed",likelihood = "binom",link = "logit",n.chain = 3)
modelrandom_bacteremia_cl<-mtc.model(network_bacteremia_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

modelrandom_bacteremia_trt<-mtc.model(network_bacteremia,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

resultsfixed_bacteremia_cl<-mtc.run(modelfixed_bacteremia_cl,n.adapt = 10000,n.iter = 100000)
Fluofixed<-relative.effect(resultsfixed_bacteremia_cl,t1="E")

forest(Fluofixed,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
legend("topleft", legend=c("B - No antibiotic", "D - Nonabsorbable","E - Fluoroquinolone", "F- Placebo","H - Cephalosporin",
                           "L - TMP-SMX","V - TMP-SMX_Not_Daily"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

resultsfixed_bacteremia_cl$deviance$DIC

resultsrandom_bacteremia_cl<-mtc.run(modelrandom_bacteremia_cl,n.adapt = 10000,n.iter = 100000)
resultsrandom_bacteremia_trt<-mtc.run(modelrandom_bacteremia_trt,n.adapt = 10000,n.iter = 100000)

resultsrandom_bacteremia_cl$deviance$DIC

#random
Fluorandom<-relative.effect(resultsrandom_bacteremia_cl,"E")
Ciprorandom<-relative.effect(resultsrandom_bacteremia_trt,"a")

forest(Fluorandom,digits=2,left.label="Favor classes",right.label="Favor fluoroquinolones")
legend("topleft", legend=c("B - No antibiotic", "C - Fluoroquinolone and Non-Anti-Pseudomonal Beta Lactam","E - Fluoroquinolone", "F - Placebo","H - Cephalosporin",
                           "J - Fluoroquinolone and Macrolide","K - Fluoroquinolone and Rifamycin","L - TMP-SMX","N - Quinolone_not_fluoro"),
       cex=.8,bty = "n",text.font = 2)

forest(Ciprorandom,digits=2,left.label="Favor others",right.label="Favors ciprofloxacin")
legend("topleft", legend=unique(paste(bacteremia$treatment,bacteremia$trt,sep = " - ")),
       cex=.7,bty = "n",text.font = 2)
gelman.diag(resultsrandom_bacteremia_cl)
mtc.nodesplit.comparisons(network_bacteremia_class)
mtc.nodesplit.comparisons(network_bacteremia)


noderandom_bacteremia_cl<-mtc.nodesplit(network_bacteremia_class,linearModel="random",likelihood =
                                          "binom",link = "log",n.chain = 3)

noderandom_bacteremia_trt<-mtc.nodesplit(network_bacteremia,linearModel="random",likelihood =
                                          "binom",link = "log",n.chain = 3)


plot(summary(noderandom_bacteremia_cl))
plot(summary(noderandom_bacteremia_trt))
legend("topleft", legend=c("B - No antibiotic","E - Fluoroquinolone", "F- Placebo","H - Cephalosporin"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

ranks_bacteremia<-rank.probability(resultsrandom_bacteremia_cl, preferredDirection=-1)
plot(ranks_bacteremia,beside=T,col=1:9,main="Probability of being the best treatment",
     names.arg=c("No antibiotic", "Fluoroquinolone and NAPB","Fluoroquinolone", "Placebo","Cephalosporin",
                 "Fluoroquinolone and Macrolide","Fluoroquinolone and Rifamycin","TMP-SMX","Quinolone_not_fluoro"),cex.names=.8)
legend("topleft", c("1","2","3","4","5","6","7","8","9"), cex=.8, bty="n",fill=1:9)

summary(network_bacteremia_class)

heterorandom_bacteremia<-mtc.anohe(network_bacteremia_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_bacteremia))
########################### fever  

fever<-antibiotics[,c("ID","FevPtNum" ,"FevPtDen","Class","trt","Author")]
fever<-na.omit(fever)
length(unique(fever$ID))
fever<-fever %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
fever<-fever[!fever$ID %in% names(which(table(fever$ID)<2)),]
colnames(fever)[colnames(fever)=="ID"] <- "study"
colnames(fever)[colnames(fever)=="FevPtNum"] <- "responders"
colnames(fever)[colnames(fever)=="FevPtDen"] <- "sampleSize"
d<-unique(fever$trt)
temp_fev<-data.frame(matrix(,nrow = 15,ncol = 2))
temp_fev$X1<-d
for (i in 1:15){
  temp_fev[i,2]<-letters[i]
}
fever <- merge(fever, temp_fev, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(fever)[colnames(fever)=="X2"] <- "treatment"


fever_class<-fever %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
fever_class<-na.omit(fever_class)

fever<-fever %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
fever<-na.omit(fever)

fever_class[25:27,]
fever_class[25,3]<-fever_class[25,3]+fever_class[26,3]
fever_class[25,4]<-fever_class[25,4]+fever_class[26,4]
fever_class<-fever_class[-26,]

fever_class<-fever_class[,1:5]

fever<-fever[,c(1:4,7)]

colnames(fever_class)[colnames(fever_class)=="Class"] <- "treatment"
network_fever_class1<-mtc.network(as.data.frame(fever_class))
network_fever<-mtc.network(as.data.frame(fever))


summary(network_fever_class1)
plot(network_fever_class1,vertex.color="coral")
plot(network_fever,vertex.color="coral")

legend("topleft", legend=c("B - No antibiotic", "C - Fluoroquinolone and Non-Anti-Pseudomonal Beta Lactam",
                           "D - Nonabsorbable","E - Fluoroquinolone", "F - Placebo","H - Cephalosporin",
                           "K - Fluoroquinolone and Rifamycin","L - TMP-SMX","N - Quinolone_not_fluoro","V - Tetracycline"),
       cex=.8,bty = "n",text.font = 2)

modelfixed_fever_cl<-mtc.model(network_fever_class,linearModel="fixed",likelihood = "binom",link = "logit",n.chain = 3)
modelrandom_fever_cl1<-mtc.model(network_fever_class1,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

modelrandom_fever_trt<-mtc.model(network_fever,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

resultsfixed_fever_cl<-mtc.run(modelfixed_fever_cl,n.adapt = 10000,n.iter = 100000)
Fluofixed<-relative.effect(resultsfixed_fever_cl,t1="E")


forest(Fluofixed,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
legend("topleft", legend=c("B - No antibiotic", "D - Nonabsorbable","E - Fluoroquinolone", "F- Placebo","H - Cephalosporin",
                           "L - TMP-SMX","V - TMP-SMX_Not_Daily"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

resultsfixed_fever_cl$deviance$DIC

##### random

resultsrandom_fever_cl1<-mtc.run(modelrandom_fever_cl1,n.adapt = 10000,n.iter = 100000)
resultsrandom_fever_trt<-mtc.run(modelrandom_fever_trt,n.adapt = 10000,n.iter = 100000)
resultsrandom_fever_cl$deviance$DIC

Fluorandom_fev1<-relative.effect(resultsrandom_fever_cl1,"E")
Ciprorandom_fev<-relative.effect(resultsrandom_fever_trt,"a")
forest(Fluorandom_fev1,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
forest(Ciprorandom_fev,digits=2,left.label="Favor others",right.label="Favor ciprofloxacin")

legend("topleft", legend=c("B - No antibiotic", "C - Fluoroquinolone and Non-Anti-Pseudomonal Beta Lactam",
                           "D - Nonabsorbable","E - Fluoroquinolone", "F - Placebo","H - Cephalosporin",
                           "K - Fluoroquinolone and Rifamycin","L - TMP-SMX","N - Quinolone_not_fluoro","V - Tetracycline"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

gelman.diag(resultsrandom_fever_cl)
mtc.nodesplit.comparisons(network_fever_class)
noderandom_fever_cl1<-mtc.nodesplit(network_fever_class1,linearModel="random",likelihood =
                                     "binom",link = "log",n.chain = 3)

noderandom_fever_trt<-mtc.nodesplit(network_fever,linearModel="random",likelihood =
                                     "binom",link = "log",n.chain = 3)

plot(summary(noderandom_fever_cl))
plot(summary(noderandom_fever_trt))

legend("topleft", legend=c("B - No antibiotic","E - Fluoroquinolone", "F- Placebo","H - Cephalosporin"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

ranks_fever<-rank.probability(resultsrandom_fever_cl, preferredDirection=-1)
plot(ranks_fever,beside=T,col=1:10,main="Probability of being the best treatment",
     names.arg=c("No antibiotic", "Fluoroquinolone and NAPB","Nonabsorbable","Fluoroquinolone", "Placebo","Cephalosporin",
                 "Fluoroquinolone and Rifamycin","TMP-SMX","Quinolone_not_fluoro","Tetracycline"),cex.names=.8)
legend("topleft", c("1","2","3","4","5","6","7","8","9","10"), cex=.7, bty="n",fill=1:10)


heterorandom_fever_trt<-mtc.anohe(network_fever,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_fever_trt))

################### mortality

mort<-antibiotics[,c("ID","AllMortNum" ,"AllMortDen","Class","trt","Author")]
mort<-na.omit(mort)
length(unique(mort$ID))
mort<-mort %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
mort<-mort[!mort$ID %in% names(which(table(mort$ID)<2)),]
colnames(mort)[colnames(mort)=="ID"] <- "study"
colnames(mort)[colnames(mort)=="AllMortNum"] <- "responders"
colnames(mort)[colnames(mort)=="AllMortDen"] <- "sampleSize"
d<-unique(mort$trt)
temp_mort<-data.frame(matrix(,nrow = 26,ncol = 2))
temp_mort$X1<-d
for (i in 1:26){
  temp_mort[i,2]<-letters[i]
}

mort <- merge(mort, temp_mort, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(mort)[colnames(mort)=="X2"] <- "treatment"


mort_class<-mort %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
mort_class<-na.omit(mort_class)

mort<-mort %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
mort<-na.omit(mort)

mort_class[21:23,]
mort_class[22,3]<-mort_class[22,3]+mort_class[21,3]
mort_class[22,4]<-mort_class[22,4]+mort_class[21,4]
mort_class<-mort_class[-21,]

mort_class<-mort_class[,1:5]

mort<-mort[,c(1:4,7)]

colnames(mort_class)[colnames(mort_class)=="Class"] <- "treatment"
network_mort_class<-mtc.network(as.data.frame(mort_class))

network_mort<-mtc.network(as.data.frame(mort))

summary(network_mort_class)
plot(network_mort_class,vertex.color="coral")
plot(network_mort,vertex.color="coral")

legend("topleft", legend=c("B - No antibiotic", "C - Fluoroquinolone and Non-Anti-Pseudomonal Beta Lactam",
                           "D - Nonabsorbable","E - Fluoroquinolone", "F - Placebo","H - Cephalosporin",
                           "J - Fluoroquinolone and Macrolide","K - Fluoroquinolone and Rifamycin","L - TMP-SMX","N - Quinolone_not_fluoro","AA - Fluoroquinolone and Mitroimidazole"),
       cex=.8,bty = "n",text.font = 2)

modelfixed_mort_cl<-mtc.model(network_mort_class,linearModel="fixed",likelihood = "binom",link = "logit",n.chain = 3)
modelrandom_mort_cl<-mtc.model(network_mort_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

modelrandom_mort_trt<-mtc.model(network_mort,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

resultsfixed_mort_cl<-mtc.run(modelfixed_mort_cl,n.adapt = 10000,n.iter = 100000)

summary(network_mort_class)
Fluofixed_mort<-relative.effect(resultsfixed_mort_cl,t1="E")


forest(Fluofixed_mort,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
legend("topleft", legend=c("B - No antibiotic", "D - Nonabsorbable","E - Fluoroquinolone", "F- Placebo","H - Cephalosporin",
                           "L - TMP-SMX"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

resultsfixed_mort_cl$deviance$DIC

####random

resultsrandom_mort_cl<-mtc.run(modelrandom_mort_cl,n.adapt = 10000,n.iter = 100000)
resultsrandom_mort_trt<-mtc.run(modelrandom_mort_trt,n.adapt = 10000,n.iter = 100000)
resultsrandom_mort_cl$deviance$DIC

Fluorandom_mort<-relative.effect(resultsrandom_mort_cl,"E")
Ciprorandom_mort<-relative.effect(resultsrandom_mort_trt,"h")
forest(Fluorandom_mort,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
forest(Ciprorandom_mort,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")

legend("topleft", legend=c("B - No antibiotic", "C - Fluoroquinolone and Non-Anti-Pseudomonal Beta Lactam",
                           "D - Nonabsorbable","E - Fluoroquinolone", "F - Placebo","H - Cephalosporin",
                           "J - Fluoroquinolone and Macrolide","K - Fluoroquinolone and Rifamycin","L - TMP-SMX","N - Quinolone_not_fluoro","AA - Fluoroquinolone and Mitroimidazole"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

gelman.diag(resultsrandom_mort_cl)
mtc.nodesplit.comparisons(network_mort_class)
noderandom_mort_cl<-mtc.nodesplit(network_mort_class,linearModel="random",likelihood =
                                    "binom",link = "log",n.chain = 3)

noderandom_mort_trt<-mtc.nodesplit(network_mort,linearModel="random",likelihood =
                                    "binom",link = "log",n.chain = 3)
plot(summary(noderandom_mort_cl))
legend("topleft", legend=c("B - No antibiotic","L - TMP-SMX"),
       cex=.8,bty = "n",ncol=1,text.font = 2)

ranks_mort<-rank.probability(resultsrandom_mort_cl, preferredDirection=-1)
plot(ranks_mort,beside=T,col=1:11,main="Probability of being the best treatment",
     names.arg=c("Fluoroquinolone and Mitroimidazole","No antibiotic", "Fluoroquinolone and NAPB","Nonabsorbable","Fluoroquinolone", "Placebo","Cephalosporin",
                 "Fluoroquinolone and Macrolide","Fluoroquinolone and Rifamycin","TMP-SMX","Quinolone_not_fluoro"),cex.names=.7)
legend("topleft", c("1","2","3","4","5","6","7","8","9","10","11"), cex=.7, bty="n",fill=1:11)


heterorandom_mort<-mtc.anohe(network_mort_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_mort))

save(network_bacteremia_class,Fluorandom,ranks_bacteremia,noderandom_bacteremia_cl,network_fever_class,
     Fluorandom_fev,ranks_fever,noderandom_fever_cl,network_mort_class,
     Fluorandom_mort,ranks_mort,noderandom_mort_cl,file="antib_lump.RData")

save(network_bacteremia_class,Fluorandom,ranks_bacteremia,noderandom_bacteremia_cl,network_fever_class,
     Fluorandom_fev,ranks_fever,noderandom_fever_cl,network_mort_class,
     Fluorandom_mort,ranks_mort,noderandom_mort_cl,file="antib3.RData")

#####################
network_bacteremia<-mtc.network(bacteremia)

bacteremia<-bacteremia %>% arrange(treatment)
colours_network_bac<-c("coral","darkseagreen","coral","coral","goldenrod","chartreuse4","coral","coral","coral","coral","coral","hotpink",
                       "coral","coral","hotpink","goldenrod","cornflowerblue","mediumslateblue","coral","coral")
plot(network_bacteremia,vertex.color=colours_network_bac)
bacteremia<-bacteremia %>% arrange(treatment)
legend("topleft", legend=unique(paste(bacteremia$treatment,bacteremia$trt,sep = " - ")),
       cex=.8,bty = "n")

colours_bac<-c("coral", "darkseagreen", "goldenrod", "chartreuse4", "hotpink",
               "cornflowerblue", "mediumslateblue")


legend("topright", legend=c("Fluoroquinolone", "Placebo","Cephalosporin","No antibiotic",
                            "TMP-SMX","Nonabsorbable","TMP-SMX_Not_Daily"),
       col="#777777",pch=21,
       pt.bg=colours_bac, pt.cex=2, cex=.8, bty="n", ncol=1)
noderandom_bacteremia_trt<-mtc.nodesplit(network_bacteremia,linearModel="random",likelihood =
                                           "binom",link = "logit",n.chain = 3)
plot(summary(noderandom_bacteremia_trt)) 

save(network_bacteremia,bacteremia,noderandom_bacteremia_trt,file="antib3.RData")

#### extra

heterorandom_bacteremia<-mtc.anohe(network_bacteremia_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
heterorandom_fever<-mtc.anohe(network_fever_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
heterorandom_mort<-mtc.anohe(network_mort_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)

save(heterorandom_bacteremia,heterorandom_fever,heterorandom_mort,file="extra.RData")

save(bacteremia,network_bacteremia,Ciprorandom,noderandom_bacteremia_trt,network_fever_class1,Fluorandom_fev1,
     network_bacteremia_class,Fluorandom,noderandom_bacteremia_cl,network_mort_class,Fluorandom_mort,noderandom_mort_cl,
     Ciprorandom_fev,network_fever,noderandom_fever_cl1,noderandom_fever_trt,fever,
     mort,network_mort,Ciprorandom_mort,noderandom_mort_trt,file="antib4.RData")

save(network_bacteremia_class,noderandom_bacteremia_cl,bacteremia_class,file="antib5.RData") 

#####################

inf<-antibiotics[,c("ID","InfMortNum","InfMortDen","Class","trt","Author")] #quitar intprim
inf<-na.omit(inf)
length(unique(inf$ID))
inf<-inf %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
inf<-inf[!inf$ID %in% names(which(table(inf$ID)<2)),]
colnames(inf)[colnames(inf)=="ID"] <- "study"
colnames(inf)[colnames(inf)=="InfMortNum"] <- "responders"
colnames(inf)[colnames(inf)=="InfMortDen"] <- "sampleSize"
d<-unique(inf$trt)
temp_inf<-data.frame(matrix(,nrow = 26,ncol = 2))
temp_inf$X1<-d
for (i in 1:26){
  temp_inf[i,2]<-letters[i]
}
inf <- merge(inf, temp_inf, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(inf)[colnames(inf)=="X2"] <- "treatment"


inf_class<-inf %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
inf_class<-na.omit(inf_class)

inf<-inf %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
inf<-na.omit(inf)
#####

inf_class %>% 
  group_by(study) %>% 
  summarise(n = n()) %>% filter(n==3 | n==4)

inf_class %>% 
  group_by(study) %>% filter(study==18 | study==23 | study==1295)



######


inf_class[15:17,]
inf_class[15,3]<-inf_class[15,3]+inf_class[16,3]
inf_class[15,4]<-inf_class[15,4]+inf_class[16,4]
inf_class<-inf_class[-16,]
inf_class[23:25,]
inf_class[23,3]<-inf_class[23,3]+inf_class[24,3]
inf_class[23,4]<-inf_class[23,4]+inf_class[24,4]
inf_class<-inf_class[-24,]

inf_class<-inf_class[,1:5]
length(unique(inf_class$study))

inf<-inf[,c(1:4,7)]

colnames(inf_class)[colnames(inf_class)=="Class"] <- "treatment"
network_inf_class<-mtc.network(as.data.frame(inf_class))
network_inf<-mtc.network(as.data.frame(inf))
plot(network_inf_class,vertex.color="coral")
plot(network_inf,vertex.color="coral")

modelrandom_inf_cl<-mtc.model(network_inf_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

modelrandom_inf_trt<-mtc.model(network_inf,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)


resultsrandom_inf_cl<-mtc.run(modelrandom_inf_cl,n.adapt = 10000,n.iter = 100000)
resultsrandom_inf_trt<-mtc.run(modelrandom_inf_trt,n.adapt = 10000,n.iter = 100000)

resultsrandom_inf_cl$deviance$DIC

Fluorandom_inf<-relative.effect(resultsrandom_inf_cl,"E")
Ciprorandom_inf<-relative.effect(resultsrandom_inf_trt,"h")
forest(Fluorandom_inf,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
forest(Ciprorandom_inf,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")

mtc.nodesplit.comparisons(network_inf_class)
noderandom_inf_cl<-mtc.nodesplit(network_inf_class,linearModel="random",likelihood =
                                   "binom",link = "log",n.chain = 3)

noderandom_inf_trt<-mtc.nodesplit(network_inf,linearModel="random",likelihood =
                                    "binom",link = "log",n.chain = 3)
plot(summary(noderandom_inf_cl))


heterorandom_inf<-mtc.anohe(network_inf_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_inf))

####################

BSI<-antibiotics[,c("ID","BSIResNum","BSIResDen","Class","trt","Author")] #quitar intprim
BSI<-na.omit(BSI)
length(unique(BSI$ID))
BSI<-BSI %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
BSI<-BSI[!BSI$ID %in% names(which(table(BSI$ID)<2)),]
colnames(BSI)[colnames(BSI)=="ID"] <- "study"
colnames(BSI)[colnames(BSI)=="BSIResNum"] <- "responders"
colnames(BSI)[colnames(BSI)=="BSIResDen"] <- "sampleSize"
d<-unique(BSI$trt)
temp_BSI<-data.frame(matrix(,nrow = 10,ncol = 2))
temp_BSI$X1<-d
for (i in 1:10){
  temp_BSI[i,2]<-letters[i]
}
BSI <- merge(BSI, temp_BSI, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(BSI)[colnames(BSI)=="X2"] <- "treatment"


BSI_class<-BSI %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
BSI_class<-na.omit(BSI_class)

BSI<-BSI %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
BSI<-na.omit(BSI)
#####

BSI_class %>% 
  group_by(study) %>% 
  summarise(n = n()) %>% filter(n==3 | n==4)




######


BSI_class<-BSI_class[,1:5]
length(unique(BSI_class$study))

BSI<-BSI[,c(1:4,7)]

colnames(BSI_class)[colnames(BSI_class)=="Class"] <- "treatment"
network_BSI_class<-mtc.network(as.data.frame(BSI_class))
network_BSI<-mtc.network(as.data.frame(BSI))
plot(network_BSI_class,vertex.color="coral")
plot(network_BSI,vertex.color="coral")

modelrandom_BSI_cl<-mtc.model(network_BSI_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

modelrandom_BSI_trt<-mtc.model(network_BSI,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)


resultsrandom_BSI_cl<-mtc.run(modelrandom_BSI_cl,n.adapt = 10000,n.iter = 100000)
resultsrandom_BSI_trt<-mtc.run(modelrandom_BSI_trt,n.adapt = 10000,n.iter = 100000)

resultsrandom_BSI_cl$deviance$DIC

Fluorandom_BSI<-relative.effect(resultsrandom_BSI_cl,"E")
Ciprorandom_BSI<-relative.effect(resultsrandom_BSI_trt,"h")
forest(Fluorandom_BSI,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
forest(Ciprorandom_BSI,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")

mtc.nodesplit.comparisons(network_BSI_class)
noderandom_BSI_cl<-mtc.nodesplit(network_BSI_class,linearModel="random",likelihood =
                                   "binom",link = "log",n.chain = 3)

noderandom_BSI_trt<-mtc.nodesplit(network_BSI,linearModel="random",likelihood =
                                    "binom",link = "log",n.chain = 3)
plot(summary(noderandom_BSI_trt))


heterorandom_BSI<-mtc.anohe(network_BSI_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_BSI))

#################
Cdiff<-antibiotics[,c("ID","CdiffNum","CdiffDen","Class","trt","Author")] #quitar intprim
Cdiff<-na.omit(Cdiff)
length(unique(Cdiff$ID))
Cdiff<-Cdiff %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
Cdiff<-Cdiff[!Cdiff$ID %in% names(which(table(Cdiff$ID)<2)),]
colnames(Cdiff)[colnames(Cdiff)=="ID"] <- "study"
colnames(Cdiff)[colnames(Cdiff)=="CdiffNum"] <- "responders"
colnames(Cdiff)[colnames(Cdiff)=="CdiffDen"] <- "sampleSize"
d<-unique(Cdiff$trt)
temp_Cdiff<-data.frame(matrix(,nrow =5,ncol = 2))
temp_Cdiff$X1<-d
for (i in 1:5){
  temp_Cdiff[i,2]<-letters[i]
}
Cdiff <- merge(Cdiff, temp_Cdiff, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(Cdiff)[colnames(Cdiff)=="X2"] <- "treatment"


Cdiff_class<-Cdiff %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
Cdiff_class<-na.omit(Cdiff_class)

Cdiff<-Cdiff %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
Cdiff<-na.omit(Cdiff)
#####

Cdiff_class %>% 
  group_by(study) %>% 
  summarise(n = n()) %>% filter(n==3 | n==4)




######


Cdiff_class<-Cdiff_class[,1:5]
length(unique(Cdiff_class$study))

Cdiff<-Cdiff[,c(1:4,7)]

colnames(Cdiff_class)[colnames(Cdiff_class)=="Class"] <- "treatment"
network_Cdiff_class<-mtc.network(as.data.frame(Cdiff_class))
network_Cdiff<-mtc.network(as.data.frame(Cdiff))
plot(network_Cdiff_class,vertex.color="coral")
plot(network_Cdiff,vertex.color="coral")

modelrandom_Cdiff_cl<-mtc.model(network_Cdiff_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)



resultsrandom_Cdiff_cl<-mtc.run(modelrandom_Cdiff_cl,n.adapt = 10000,n.iter = 100000)


resultsrandom_Cdiff_cl$deviance$DIC

Fluorandom_Cdiff<-relative.effect(resultsrandom_Cdiff_cl,"E")
forest(Fluorandom_Cdiff,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")


mtc.nodesplit.comparisons(network_Cdiff_class)


heterorandom_Cdiff<-mtc.anohe(network_Cdiff_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_Cdiff))

#############

IFD<-antibiotics[,c("ID","IFDNum","IFDDem","Class","trt","Author")] #quitar intprim
IFD<-na.omit(IFD)
length(unique(IFD$ID))
IFD<-IFD %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
IFD<-IFD[!IFD$ID %in% names(which(table(IFD$ID)<2)),]
colnames(IFD)[colnames(IFD)=="ID"] <- "study"
colnames(IFD)[colnames(IFD)=="IFDNum"] <- "responders"
colnames(IFD)[colnames(IFD)=="IFDDem"] <- "sampleSize"
d<-unique(IFD$trt)
temp_IFD<-data.frame(matrix(,nrow =20,ncol = 2))
temp_IFD$X1<-d
for (i in 1:20){
  temp_IFD[i,2]<-letters[i]
}
IFD <- merge(IFD, temp_IFD, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(IFD)[colnames(IFD)=="X2"] <- "treatment"


IFD_class<-IFD %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
IFD_class<-na.omit(IFD_class)

IFD<-IFD %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
IFD<-na.omit(IFD)
#####

IFD_class %>% 
  group_by(study) %>% 
  summarise(n = n()) %>% filter(n==3 | n==4)

IFD_class[9:11,]
IFD_class[9,3]<-IFD_class[9,3]+IFD_class[10,3]
IFD_class[9,4]<-IFD_class[9,4]+IFD_class[10,4]
IFD_class<-IFD_class[-10,]
IFD_class[37:39,]
IFD_class[38,3]<-IFD_class[38,3]+IFD_class[39,3]
IFD_class[38,4]<-IFD_class[38,4]+IFD_class[39,4]
IFD_class<-IFD_class[-39,]


######


IFD_class<-IFD_class[,1:5]
length(unique(IFD_class$study))

IFD<-IFD[,c(1:4,7)]

colnames(IFD_class)[colnames(IFD_class)=="Class"] <- "treatment"
network_IFD_class<-mtc.network(as.data.frame(IFD_class))
network_IFD<-mtc.network(as.data.frame(IFD))
plot(network_IFD_class,vertex.color="coral")
plot(network_IFD,vertex.color="coral")

modelrandom_IFD_cl<-mtc.model(network_IFD_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

modelrandom_IFD_trt<-mtc.model(network_IFD,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)


resultsrandom_IFD_cl<-mtc.run(modelrandom_IFD_cl,n.adapt = 10000,n.iter = 100000)
resultsrandom_IFD_trt<-mtc.run(modelrandom_IFD_trt,n.adapt = 10000,n.iter = 100000)

resultsrandom_IFD_cl$deviance$DIC

Fluorandom_IFD<-relative.effect(resultsrandom_IFD_cl,"E")
Ciprorandom_IFD<-relative.effect(resultsrandom_IFD_trt,"h")
forest(Fluorandom_IFD,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
forest(Ciprorandom_IFD,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")

mtc.nodesplit.comparisons(network_IFD_class)
noderandom_IFD_cl<-mtc.nodesplit(network_IFD_class,linearModel="random",likelihood =
                                   "binom",link = "log",n.chain = 3)

noderandom_IFD_trt<-mtc.nodesplit(network_IFD,linearModel="random",likelihood =
                                    "binom",link = "log",n.chain = 3)
plot(summary(noderandom_IFD_trt))


heterorandom_IFD<-mtc.anohe(network_IFD_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_IFD))

###############

MSK<-antibiotics[,c("ID","MSKNum","MSKDem","Class","trt","Author")] #quitar intprim
MSK<-na.omit(MSK)
length(unique(MSK$ID))
MSK<-MSK %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
MSK<-MSK[!MSK$ID %in% names(which(table(MSK$ID)<2)),]
colnames(MSK)[colnames(MSK)=="ID"] <- "study"
colnames(MSK)[colnames(MSK)=="MSKNum"] <- "responders"
colnames(MSK)[colnames(MSK)=="MSKDem"] <- "sampleSize"
d<-unique(MSK$trt)
temp_MSK<-data.frame(matrix(,nrow =6,ncol = 2))
temp_MSK$X1<-d
for (i in 1:6){
  temp_MSK[i,2]<-letters[i]
}
MSK <- merge(MSK, temp_MSK, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(MSK)[colnames(MSK)=="X2"] <- "treatment"


MSK_class<-MSK %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
MSK_class<-na.omit(MSK_class)

MSK<-MSK %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
MSK<-na.omit(MSK)
#####

MSK_class %>% 
  group_by(study) %>% 
  summarise(n = n()) %>% filter(n==3 | n==4)




######


MSK_class<-MSK_class[,1:5]
length(unique(MSK_class$study))

MSK<-MSK[,c(1:4,7)]

colnames(MSK_class)[colnames(MSK_class)=="Class"] <- "treatment"
network_MSK_class<-mtc.network(as.data.frame(MSK_class))
network_MSK<-mtc.network(as.data.frame(MSK))
plot(network_MSK_class,vertex.color="coral")
plot(network_MSK,vertex.color="coral")

modelrandom_MSK_cl<-mtc.model(network_MSK_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)




resultsrandom_MSK_cl<-mtc.run(modelrandom_MSK_cl,n.adapt = 10000,n.iter = 100000)


resultsrandom_MSK_cl$deviance$DIC

Fluorandom_MSK<-relative.effect(resultsrandom_MSK_cl,"E")

forest(Fluorandom_MSK,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")

heterorandom_MSK<-mtc.anohe(network_MSK_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_MSK))


##############

FNP<-antibiotics[,c("ID","FNPtNum","FNPtDen","Class","trt","Author")] #quitar intprim
FNP<-na.omit(FNP)
length(unique(FNP$ID))
FNP<-FNP %>% arrange(ID)  %>%
  group_by(ID) %>%
  filter(any(Class=="E",Class=="H"))
FNP<-FNP[!FNP$ID %in% names(which(table(FNP$ID)<2)),]
colnames(FNP)[colnames(FNP)=="ID"] <- "study"
colnames(FNP)[colnames(FNP)=="FNPtNum"] <- "responders"
colnames(FNP)[colnames(FNP)=="FNPtDen"] <- "sampleSize"
d<-unique(FNP$trt)
temp_FNP<-data.frame(matrix(,nrow =10,ncol = 2))
temp_FNP$X1<-d
for (i in 1:10){
  temp_FNP[i,2]<-letters[i]
}
FNP <- merge(FNP, temp_FNP, by.x=c("trt"), by.y=c("X1"),all = T)
colnames(FNP)[colnames(FNP)=="X2"] <- "treatment"


FNP_class<-FNP %>% arrange(study) %>%
  group_by(study) %>%
  mutate(Class = replace(Class, n_distinct(Class)==1, NA) )
FNP_class<-na.omit(FNP_class)

FNP<-FNP %>% arrange(study) %>%
  group_by(study) %>%
  mutate(treatment = replace(treatment, n_distinct(treatment)==1, NA) )
FNP<-na.omit(FNP)
#####

FNP_class %>% 
  group_by(study) %>% 
  summarise(n = n()) %>% filter(n==3 | n==4)




######


FNP_class<-FNP_class[,1:5]
length(unique(FNP_class$study))

FNP<-FNP[,c(1:4,7)]

colnames(FNP_class)[colnames(FNP_class)=="Class"] <- "treatment"
network_FNP_class<-mtc.network(as.data.frame(FNP_class))
network_FNP<-mtc.network(as.data.frame(FNP))
plot(network_FNP_class,vertex.color="coral")
plot(network_FNP,vertex.color="coral")

modelrandom_FNP_cl<-mtc.model(network_FNP_class,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)

modelrandom_FNP_trt<-mtc.model(network_FNP,linearModel="random",likelihood = "binom",link = "log",n.chain = 3)


resultsrandom_FNP_cl<-mtc.run(modelrandom_FNP_cl,n.adapt = 10000,n.iter = 100000)
resultsrandom_FNP_trt<-mtc.run(modelrandom_FNP_trt,n.adapt = 10000,n.iter = 100000)

resultsrandom_FNP_cl$deviance$DIC

Fluorandom_FNP<-relative.effect(resultsrandom_FNP_cl,"E")
Ciprorandom_FNP<-relative.effect(resultsrandom_FNP_trt,"h")
forest(Fluorandom_FNP,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")
forest(Ciprorandom_FNP,digits=2,left.label="Favor classes",right.label="Favor Fluoroquinolones")

mtc.nodesplit.comparisons(network_FNP)

noderandom_FNP_trt<-mtc.nodesplit(network_FNP,linearModel="random",likelihood =
                                    "binom",link = "log",n.chain = 3)
plot(summary(noderandom_FNP_trt))


heterorandom_FNP<-mtc.anohe(network_FNP_class,linearModel="random",likelihood = "binom",link = "logit",n.chain = 3)
plot(summary(heterorandom_FNP))


#####

save(inf,network_inf,Ciprorandom_inf,noderandom_inf_trt,network_BSI_class,Fluorandom_BSI,
     network_inf_class,Fluorandom_inf,noderandom_inf_cl,network_IFD_class,Fluorandom_IFD,noderandom_IFD_cl,
     Ciprorandom_BSI,network_BSI,noderandom_BSI_cl,noderandom_BSI_trt,BSI,
     IFD,network_IFD,Ciprorandom_IFD,noderandom_IFD_trt,Cdiff,network_Cdiff_class,Fluorandom_Cdiff,
     MSK,network_MSK_class,Fluorandom_MSK,FNP,network_FNP,Ciprorandom_FNP,noderandom_FNP_trt,network_FNP_class,Fluorandom_FNP,
     file="antib6.RData")

