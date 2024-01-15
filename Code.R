library("nlme")
library(foreign)
library(pROC)
library(ggplot2)
library(ggrepel)
library("survival")
library("survminer")


#### Linear modeling Roc curves

combdata <- read.csv("C:/Users/wppwjw/OneDrive - Cardiff University/Work/Wioleta/Jan 24 rebut/Data_zscored.csv")
combdata <- as.data.frame(combdata)

mydata = combdata

control = subset(combdata, MS == 0)
case  = subset(combdata, MS == 1)

#list of all new genes to add
vnames_c <- names(mydata[,2:21])  #just CSF
len_c = length(vnames_c)
vnames_s <- names(mydata[,22:45])   #just Serum
len_s = length(vnames_s)
vnames_cs <- names(mydata[,2:45])    #Both
len_cs = length(vnames_cs)

new_dat = cbind(mydata[,2:45],MS = mydata$MS,Sex = mydata$Sex,age = mydata$age)
new_dat_backup = new_dat


# Create random integer from 1 to 4 to allow Split Data into Training and Testing
set.seed(888) #1234 or 888 or whatever
ran = floor(runif(nrow(new_dat),1,4.999999))
new_dat = cbind(new_dat,ran)

#first - one biomarker
#create dataframe
out = NULL
out <- data.frame(Integers=integer(), Characters=character(), Doubles=double(), Integers=integer(), 
                  Integers=integer(), Integers=integer(), Doubles=double(), Doubles=double(), Doubles=double(), 
                  Doubles=double(), Doubles=double(), Doubles=double(), Doubles=double(), Doubles=double(), 
                  Doubles=double(), Doubles=double())
colnames(out) <- c('Num', 'Biomarkers', 'AUC', 'End','Patients','Variables','Train_AUC1','Test_AUC1',
                   'Train_AUC2','Test_AUC2','Train_AUC3','Test_AUC3','Train_AUC4','Test_AUC4',
                   'Train_Mean','Test_Mean')
out$Biomarkers = as.character(out$Biomarkers)
out_start = out
out_tot = out

#set the below depending on  CSF or Serum or both

#vnames = vnames_s
#ng = len_s
#vnames = vnames_c
#ng = len_c
vnames = vnames_cs
ng = len_cs

num_var = 1
for(i in 1:ng) #ng = number of biomarkers
{ 
  #construct regression equation initially using all data for train and test
  stv = paste(vnames[i])
  V = as.formula(paste("MS ~ ",stv))
#  V = as.formula(paste("MS ~ Sex + age + ",stv))
  fit = glm(V, family = binomial, data = new_dat) 
  p = predict(fit, newdata = new_dat, type = "response")
  roc.info = roc(new_dat$MS, p, quiet = TRUE)
  out[i,1] = i
  out[i,2] = stv
  out[i,3] = roc.info$auc
  out[i,4] = i
  out[i,5] = fit$df.null + num_var #how many samples used in the regression
  out[i,6] = num_var
  
  #Loop through the 4 combinations of 75% Train and 25% Test
  k = 7
  for(j in 1:4){
    Train = subset(new_dat, new_dat$ran != j)
    Test = subset(new_dat, new_dat$ran == j)
    fit = glm(V, family = binomial, data = Train)
    p = predict(fit, newdata = Train, type = "response")
    roc.info = roc(Train$MS, p, quiet = TRUE)
    out[i,k] = roc.info$auc
    p = predict(fit, newdata = Test, type = "response")
    roc.info = roc(Test$MS, p, quiet = TRUE)
    out[i,k+1] = roc.info$auc
    k = k + 2
  }
  
  #Finf the mean AUC's for Train and Test
  out[i,15] = (out[i,7] + out[i,9] + out[i,11] + out[i,13])/4
  out[i,16] = (out[i,8] + out[i,10] + out[i,12] + out[i,14])/4
  print(i)
} 

out_order <- out[order(out$Test_Mean, decreasing = TRUE),]
out_order



#out_order1clog = out_order  #c complete csf, s complete serum
#out_order1slog = out_order  #c complete csf, s complete serum
out_order1cslog = out  #c complete csf, s complete serum


#multiple biomarkers

out_start = out_order1cslog
out_tot = out_order1cslog
out = out_order1cslog

#vnames = vnames_s
#ng = len_s
#vnames = vnames_c
#ng = len_c
vnames = vnames_cs
ng = len_cs
 

b_max = 3    #adjustable max number of biomarkers
for(k in 2:b_max)   #
{
  print(paste("iteration with  ",k," biomarkers"))
  nrow = nrow(out)
  new_row = out[,2]
  roc = out[,3]
  done = out[,4]   #this ensures previous biomarkers are not re-added in this iteration 
  num_patient = out[,5]
  
  icount = 0
  out = out_start
  for(j in 1:nrow)
  {
      stvg = new_row[j]
      i = done[j]+1
      while(i <= ng)
      {
        stv = paste(stvg," + ",vnames[i])
 #       V = as.formula(paste("MS ~ ",stv))    #with or without sex and age
        V = as.formula(paste("MS ~ Sex + age + ",stv))
        
        #first use all data for both Train and Test
        fit = glm(V, family = binomial, data = new_dat) 
        p = predict(fit, newdata = new_dat, type = "response")
        roc.info = NULL
        roc.info = roc(new_dat$MS, p, quiet = TRUE)
        icount = icount + 1
        out[icount,1] = icount
        out[icount,2] = stv
        out[icount,3] = roc.info$auc
        out[icount,4] = i
        out[icount,5] = fit$df.null + 1 #how many samples used in the regression
        out[icount,6] = k
          
        l = 7
        for(j in 1:4){
          Train = subset(new_dat, new_dat$ran != j)
          Test = subset(new_dat, new_dat$ran == j)
          fit = glm(V, family = binomial, data = Train) 
          p = predict(fit, newdata = Train, type = "response")
          roc.info = roc(Train$MS, p, quiet = TRUE)
          out[icount,l] = roc.info$auc
          p = predict(fit, newdata = Test, type = "response")
          roc.info = roc(Test$MS, p, quiet = TRUE)
          out[icount,l+1] = roc.info$auc
          l = l + 2
        }
        out[icount,15] = (out[icount,7] + out[icount,9] + out[icount,11] + out[icount,13])/4
        out[icount,16] = (out[icount,8] + out[icount,10] + out[icount,12] + out[icount,14])/4

        i = i + 1
      }
  }
  out_order <- out[order(out$Test_Mean, decreasing = TRUE),]
  out_tot = rbind(out_tot,out_order)
####  
  nout_order = min(nrow(out_order),200)  #these 2 lines allow choice of how many ordered biomarkers combinations from
  out = out_order[1:nout_order,]         #this iteration are used in the next iteration with 1 more biomarker added
####
}
out_order <- out_tot[order(out_tot$Test_Mean, decreasing = TRUE),]
head(out_order,30)

#write out_order to file

out_str = paste0("C:/Users/wppwjw/OneDrive - Cardiff University/Work/Wioleta/Jan 24 rebut/new test.csv") #
write.csv(out_order, out_str,col.names = TRUE)




### Survival analysis

# read in whole datafile

combdata <- read.csv("C:/Users/wppwjw/OneDrive - Cardiff University/Work/Wioleta/Jan 24 rebut/Data_zscored.csv")
combdata <- as.data.frame(combdata)

mydata = combdata




#list of all new genes to add
vnames_c <- names(mydata[,2:21])  #just CSF
len_c = length(vnames_c)
vnames_s <- names(mydata[,22:45])   #just Serum
len_s = length(vnames_s)
vnames_cs <- names(mydata[,2:45])    #Both
len_cs = length(vnames_cs)

nc = ncol(mydata)


bio_dat = mydata 
bio_dat = subset(bio_dat, TimetoEDSS6 >= 0)
data = bio_dat
vnames = vnames_cs

out = NULL
out <- data.frame(Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),
                  Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),
                  Integers=integer())
num_bio = 1
len = length(vnames)
ng = len
maxc = 0

for(i in 1:len)  
{
  stv = paste(vnames[i])
  
#  EDSS6
  V = as.formula(paste("Surv(EDSS6, ReachedEDSS6) ~ ",stv, "+ Sex + DMT + age" )) 
  
#  Time to Relapse
#   V = as.formula(paste("Surv(TimeToNewRelapseCensor, NewRelapseCensored) ~ ",stv, "+ Sex + age" ))
  
  
  res.cox <- coxph(V, data = bio_dat)
  x = summary(res.cox)  
  out_tmp = c(num_bio, x$coefficients[1,c(1,3,5)],x$conf.int[1,c(1,3,4)],x$n,x$nevent,x$concordance[1],x$rsq[1],extractAIC(res.cox)[2],i)
  out = rbind(out,out_tmp)
}

out = cbind(vnames[], out)

colnames(out) = c("Biomarkers","Num_Biomarkers","Beta","SE","p","HR","HR_CI_Low","HR_CI_High","N","Nevents","Concordance","Rsq","AIC","Pos")
out = data.frame(out)
out <- out[order(out$Concordance, decreasing = TRUE),] 

out = round_df(out, 6, rf = "round")

out

##### Add additional biomarkers
out_tmp = out
out_hld = out


for(k in 2:3)   #2:6
{
  num_bio = k
  nrow = nrow(out_tmp)
  stvg = out_tmp[,1]
  po = out_tmp[,14]
  out_tmp <- data.frame(Character = character(), Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),
                        Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),Doubles=double(),
                        Doubles=double(),Doubles=double(),Doubles=double(),Integers=integer(),
                        stringsAsFactors = FALSE)
  
  for(j in 1:nrow)
  {
    print(paste(k,j,nrow))
    
    i = as.numeric(as.character(po[j])) + 1  #ensures biomarkers aren't duplicated in the same iteration

    while(i <= ng)
    {
      stv = paste(stvg[j]," + ",vnames[i])
   
#  EDSS6
      V = as.formula(paste("Surv(EDSS6, ReachedEDSS6) ~ ",stv, "+ Sex + DMT + age" )) 
      
#  Time to Relapse
      #   V = as.formula(paste("Surv(TimeToNewRelapseCensor, NewRelapseCensored) ~ ",stv, "+ Sex + age" ))
      
      res.cox <- coxph(V, data = bio_dat)
      x = summary(res.cox)

      out_tmp1 = c(stv, num_bio, x$coefficients[1,c(1,3,5)],x$conf.int[1,c(1,3,4)],x$n,x$nevent,x$concordance[1],x$rsq[1],extractAIC(res.cox)[2],i)
      out_tmp = rbind(out_tmp,out_tmp1)
      
      i = i + 1
      
    }
  }
  colnames(out_tmp) = c("Biomarkers","Num_Biomarkers","Beta","SE","p","HR","HR_CI_Low","HR_CI_High","N","Nevents","Concordance","Rsq","AIC","Pos")
  
  #  out_tmp <- out_tmp[order(out_tmp$p, decreasing = TRUE),]
  out_tmp <- out_tmp[order(out_tmp$Concordance, decreasing = TRUE),]  
  out = rbind(out,out_tmp)
  
###  the 4 lines below allow a cut down (the best performing) number of biomarkers passed onto the next iteration
  cut_down = 1500
  ll = nrow(out_tmp)
  ll = min(ll,cut_down)
  out_tmp = out_tmp[1:ll,]
###
}

head(out,20)

#out_str = paste0("    path   ") #
#write.csv(out, out_str,col.names = TRUE)