######
###### Code for model building using age and sex adjusted protein level data
#-----------------------------------------------------------------------------------------------------------------------


##### Adjust the candidate protein levels by the effects of age and sex
#---------------------------------------------------------------------------------------------
library(ggplot2)
library(plotly)

setwd('/PATH')

Merge_data = read.csv('./Merge_raw_protein_levels.csv', header = T)
Merge_Pheno_data=Merge_data[c(1:4)]
Merge_data_control<-subset(Merge_data,Merge_data$Diagnosis=="pTau-CN")

### Examine the effects of age and gender on each protein level and generate the adjusted protein levels
library(robustbase)

Adjusted_Merge_data=data.frame(matrix(ncol = 18,nrow = length(Merge_data[,1])))
colnames(Adjusted_Merge_data)=names(Merge_data)[5:22]
Adjusted_Merge_data=cbind(Merge_Pheno_data,Adjusted_Merge_data)

Age_Sex_Effects_on_proteins=data.frame(matrix(ncol = 8,nrow = 18))
colnames(Age_Sex_Effects_on_proteins)=c("Age_Estimate","Age_SE","Age_T_value","Age_P_value","Gender_Estimate","Gender_SE","Gender_T_value","Gender_P_value")
rownames(Age_Sex_Effects_on_proteins)=names(Merge_data)[5:22]

for (k in 1:18){
  Protein_age_sex_test<-lmrob(Merge_data_control[,k+4]~Age+Gender,data=Merge_data_CN,k.max=900000)
  
  Age_Sex_Effects_on_proteins[k,1:4]=as.vector(summary(Protein_age_sex_test)$coefficients[2,])
  Age_Sex_Effects_on_proteins[k,5:8]=as.vector(summary(Protein_age_sex_test)$coefficients[3,])
  
  Adjusted_Merge_data[,k+4]=Merge_data[,k+4]-Merge_data[,2]*as.numeric(Age_Sex_Effects_on_proteins[k,1])-Merge_data[,3]*as.numeric(Age_Sex_Effects_on_proteins[k,5])
}

write.csv(Age_Sex_Effects_on_proteins,'./Summary_AgeSex_Effects.csv',row.names = T)

rownames(Adjusted_Merge_data)=Merge_data$ID
write.csv(Adjusted_Merge_data,'./Merge_AgeSex_adjusted_protein_levels.csv',row.names = T)


#-----------------------------------------------------------------------------------------------------------------------
## Model building based on the adjusted protein levels
##### Input overall file
setwd('/PATH')

Merge_data = read.csv('./Merge_AgeSex_adjusted_protein_levels.csv', header = T)
Merge_data_control<-subset(Merge_data,Merge_data$Diagnosis=="pTau-CN")
Merge_data_AD<-subset(Merge_data,Merge_data$Diagnosis=="pTau+AD")
Merge_data_for_model=rbind(Merge_data_control,Merge_data_AD)
Merge_data_for_model$Diagnosis<-relevel(factor(Merge_data_for_model$Diagnosis),ref = "pTau-CN")

#----------------------------------------------------------------------
##### Model building
#----------------------------------------------------------------------
library(pROC)

Data_for_model=Merge_data_for_model[c(4,5:22)]
Data_for_model=na.omit(Data_for_model)

linmod<-glm(Diagnosis~.,data= Data_for_model,family = binomial)
pred=predict(linmod, type=c("response"))
Data_for_model$AD_probability_score=100*pred
graph<-plot.roc(Data_for_model$Diagnosis,Data_for_model$AD_probability_score,percent=TRUE,col="#008600")
plot(graph)
auc(graph)

Model_estimate=data.frame(matrix(ncol = 6,nrow = 19))
colnames(Model_estimate)=c("Protein_ID","Estimate","SE","Z_value","P_value")
Model_estimate[1,1]="Intercept"
Model_estimate[1,2:5]=as.vector(summary(linmod)$coefficients[1,1:4])
Model_estimate[2:19,1]=names(Data_for_model)[2:19]

for (q in (2:19)){
  Model_estimate[q,2:5]=as.vector(summary(linmod)$coefficients[q,1:4])
}

write.csv(Model_estimate,"./Summary_model_estimate.csv",row.names = F)
