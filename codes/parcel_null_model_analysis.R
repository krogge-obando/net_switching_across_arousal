#This code will compute the null model results for the parcelation FIND atlas

#Load packaages
#install.packages("rlist")
library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(lsr)
library(reshape2)
library(rlist)


#writefunction

FIND_compute_null_pvalue<-function(df_null) {
  
  #df_null<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_parcels_100_null_switching_rate_results.csv",header=0)
  
  
  iterations<-rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
                    24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
                    44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
                    64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
                    84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100),times=90)
  
  iterations<-sort(iterations)

  parcels<-c("d_dmn1","d_dmn2","d_dmn3","d_dmn4","d_dmn5","d_dmn6","d_dmn7","d_dmn8","d_dmn9",
             "v_dmn1","v_dmn2","v_dmn3","v_dmn4","v_dmn5","v_dmn6","v_dmn7","v_dmn8","v_dmn9","v_dmn10",
             "a_sal1","a_sal2","a_sal3","a_sal4","a_sal5","a_sal6","a_sal7",
             "p_sal1","p_sal2","p_sal3","p_sal4","p_sal5","p_sal6","p_sal7","p_sal8","p_sal9","p_sal10","p_sal11","p_sal12",
             "lcen_1","lcen_2","lcen_3","lcen_4","lcen_5","lcen_6",
             "rcen_1","rcen_2","rcen_3","rcen_4","rcen_5","rcen_6",
             "aud_1","aud_2","aud_3",
             "basg_1","basg_2","basg_3","basg_4","basg_5",
             "pVis_1","pVis_2",
             "hVis_1","hVis_2",
             "lang_1","lang_2","lang_3","lang_4","lang_5","lang_6","lang_7",
             "prec_1","prec_2","prec_3","prec_4",
             "sensmotor_1","sensmotor_2","sensmotor_3","sensmotor_4","sensmotor_5","sensmotor_6",
             "vis_spat1","vis_spat2","vis_spat3","vis_spat4","vis_spat5","vis_spat6","vis_spat7","vis_spat8","vis_spat9","vis_spat10", "vis_spat11")
 
 parcels<-as.data.frame(parcels)
  
 parcels<-parcels[order(factor(parcels$parcels,levels=c(c("d_dmn","v_dmn","a_sal","p_sal","l_cen","r_cen","aud",
                                                          "bas_g","pVis","hVis","lang","prec","sensmotor","vis_spat")))),]
 df_null_labled<-rbind(iterations,df_null)
 df_null_labled<-rbind(parcels,df_null_labled)
 
 rownames(df_null_labled)<-c('parcels','iterations','vcon02_scan01','vcon02_scan02','vcon03_scan01','vcon04_scan02','vcon05_scan01','vcon07_scan01','vcon08_scan01',
                             'vcon08_scan02', 'vcon09_scan01','vcon09_scan02','vcon10_scan01','vcon10_scan02','vcon11_scan02','vcon12_scan01',
                             'vcon13_scan01','vcon14_scan01','vcon14_scan02','vcon15_scan01','vcon15_scan02','vcon17_scan01','vcon17_scan02',
                             'vcon18_scan01','vcon18_scan02','vcon19_scan01','vcon19_scan02','vcon20_scan01','vcon20_scan02','vcon22_scan01',
                             'vcon22_scan02','vcon23_scan01', 'vcon24_scan01','vcon25_scan01','vcon25_scan02', 'vcon26_scan01','vcon27_scan01',
                             'vcon27_scan02','vcon29_scan01','vcon29_scan02','vcon30_scan01','vcon32_scan01','vcon33_scan01','vcon33_scan02',
                             'vcon35_scan01','vcon36_scan01','vcon37_scan01')
 
 alert_scans<-df_null_labled[c("parcels","iterations","vcon02_scan01","vcon07_scan01","vcon15_scan01","vcon15_scan02","vcon17_scan02","vcon19_scan02","vcon20_scan01", "vcon20_scan02","vcon29_scan01","vcon29_scan02"),]
 drowsy_scans<-df_null_labled[c("parcels","iterations", "vcon05_scan01", "vcon08_scan01", "vcon08_scan02", "vcon09_scan01", "vcon09_scan02", "vcon12_scan01", "vcon14_scan01",
                                "vcon22_scan01", "vcon23_scan01", "vcon24_scan01", "vcon25_scan02", "vcon32_scan01", "vcon36_scan01", "vcon37_scan01"),]
 
 #generating df_of each network for alert and drowsy
 
 alert_scans<-as.data.frame(t(alert_scans))
 drowsy_scans<-as.data.frame(t(drowsy_scans))
 
 alert_scans<-melt(alert_scans,id=c("parcels","iterations"))
 drowsy_scans<-melt(drowsy_scans,id=c("parcels","iterations"))
 


#forloop to generate 100 cohen D values of psal
df_psal_cohenD_list <- data.frame(matrix(ncol=100, nrow=12))

for(j in 1:12) {
  output_name_a<-unlist(paste0('df_','psal',j,'_alert','_inter'),)
  output_name_d<-unlist(paste0('df_','psal',j,'_drowsy','_inter'),)
  input<-unlist(paste0('p_sal',j),)
  new_alert_df<-assign(output_name_a,subset(alert_scans,alert_scans$parcels==input))
  new_drowsy_df<-assign(output_name_d,subset(drowsy_scans,drowsy_scans$parcels==input))
  
for (i in 1:100) {
  output_a<-subset(new_alert_df,new_alert_df$iterations==i)
  output_d<-subset(new_drowsy_df,new_drowsy_df$iterations==i)
  
  df_psal_cohenD<-cohensD(as.numeric(output_a$value),as.numeric(output_d$value))
  df_psal_cohenD_list[j,i]<-df_psal_cohenD}
  df_psal_cohenD_list[is.na(df_psal_cohenD_list)] <- 0}

df_psal_cohenD_list<-t(df_psal_cohenD_list)
colnames(df_psal_cohenD_list)<-c("p_sal1","p_sal2","p_sal3","p_sal4","p_sal5","p_sal6","p_sal7","p_sal8","p_sal9","p_sal10","p_sal11","p_sal12")


df_asal_cohenD_list <- data.frame(matrix(ncol=100, nrow=7))

for(j in 1:7) {
  output_name_a<-unlist(paste0('df_','asal',j,'_alert','_inter'),)
  output_name_d<-unlist(paste0('df_','asal',j,'_drowsy','_inter'),)
  input<-unlist(paste0('a_sal',j),)
  new_alert_df<-assign(output_name_a,subset(alert_scans,alert_scans$parcels==input))
  new_drowsy_df<-assign(output_name_d,subset(drowsy_scans,drowsy_scans$parcels==input))
  
  for (i in 1:100) {
    output_a<-subset(new_alert_df,new_alert_df$iterations==i)
    output_d<-subset(new_drowsy_df,new_drowsy_df$iterations==i)
    output_a[is.na(output_a)] <- 0
    output_d[is.na(output_d)] <- 0
    df_asal_cohenD<-cohensD(as.numeric(output_a$value),as.numeric(output_d$value))
    df_asal_cohenD_list[j,i]<-df_asal_cohenD
    df_asal_cohenD_list[is.na(df_asal_cohenD_list)] <- 0}}

df_asal_cohenD_list<-t(df_asal_cohenD_list)
df_asal1_drowsy_inter
colnames(df_asal_cohenD_list)<-c("a_sal1","a_sal2","a_sal3","a_sal4","a_sal5","a_sal6","a_sal7")

df_lcen_cohenD_list <- data.frame(matrix(ncol=100, nrow=6))

for(j in 1:6) {
  output_name_a<-unlist(paste0('df_','lcen',j,'_alert','_inter'),)
  output_name_d<-unlist(paste0('df_','lcen',j,'_drowsy','_inter'),)
  input<-unlist(paste0('lcen_',j),)
  new_alert_df<-assign(output_name_a,subset(alert_scans,alert_scans$parcels==input))
  new_drowsy_df<-assign(output_name_d,subset(drowsy_scans,drowsy_scans$parcels==input))
  
  for (i in 1:100) {
    output_a<-subset(new_alert_df,new_alert_df$iterations==i)
    output_d<-subset(new_drowsy_df,new_drowsy_df$iterations==i)
    df_lcen_cohenD<-cohensD(as.numeric(output_a$value),as.numeric(output_d$value))
    df_lcen_cohenD_list[j,i]<-df_lcen_cohenD
    df_lcen_cohenD_list[is.na(df_lcen_cohenD_list)] <- 0}}

df_lcen_cohenD_list<-t(df_lcen_cohenD_list)
colnames(df_lcen_cohenD_list)<-c("lcen_1","lcen_2","lcen_3","lcen_4","lcen_5","lcen_6")

df_rcen_cohenD_list <- data.frame(matrix(ncol=100, nrow=6))

for(j in 1:6) {
  output_name_a<-unlist(paste0('df_','rcen',j,'_alert','_inter'),)
  output_name_d<-unlist(paste0('df_','rcen',j,'_drowsy','_inter'),)
  input<-unlist(paste0('rcen_',j),)
  new_alert_df<-assign(output_name_a,subset(alert_scans,alert_scans$parcels==input))
  new_drowsy_df<-assign(output_name_d,subset(drowsy_scans,drowsy_scans$parcels==input))
  
  for (i in 1:100) {
    output_a<-subset(new_alert_df,new_alert_df$iterations==i)
    output_d<-subset(new_drowsy_df,new_drowsy_df$iterations==i)
    output_a[is.na(output_a)] <- 0
    output_d[is.na(output_d)] <- 0
    df_rcen_cohenD<-cohensD(as.numeric(output_a$value),as.numeric(output_d$value))
    df_rcen_cohenD_list[j,i]<-df_rcen_cohenD
    df_rcen_cohenD_list[is.na(df_rcen_cohenD_list)] <- 0}}

df_rcen_cohenD_list<-t(df_rcen_cohenD_list)
colnames(df_rcen_cohenD_list)<-c("rcen_1","rcen_2","rcen_3","rcen_4","rcen_5","rcen_6")

df_ddmn_cohenD_list <- data.frame(matrix(ncol=100, nrow=9))

for(j in 1:9) {
  output_name_a<-unlist(paste0('df_','ddmn',j,'_alert','_inter'),)
  output_name_d<-unlist(paste0('df_','ddmn',j,'_drowsy','_inter'),)
  input<-unlist(paste0('d_dmn',j),)
  new_alert_df<-assign(output_name_a,subset(alert_scans,alert_scans$parcels==input))
  new_drowsy_df<-assign(output_name_d,subset(drowsy_scans,drowsy_scans$parcels==input))
  
  for (i in 1:100) {
    output_a<-subset(new_alert_df,new_alert_df$iterations==i)
    output_d<-subset(new_drowsy_df,new_drowsy_df$iterations==i)
    
    df_ddmn_cohenD<-cohensD(as.numeric(output_a$value),as.numeric(output_d$value))
    df_ddmn_cohenD_list[j,i]<-df_ddmn_cohenD
    df_ddmn_cohenD_list[is.na(df_ddmn_cohenD_list)] <- 0}}

df_ddmn_cohenD_list<-t(df_ddmn_cohenD_list)
colnames(df_ddmn_cohenD_list)<-c("d_dmn1","d_dmn2","d_dmn3","d_dmn4","d_dmn5","d_dmn6","d_dmn7","d_dmn8","d_dmn9")

df_vdmn_cohenD_list <- data.frame(matrix(ncol=100, nrow=10))

for(j in 1:10) {
  output_name_a<-unlist(paste0('df_','vdmn',j,'_alert','_inter'),)
  output_name_d<-unlist(paste0('df_','vdmn',j,'_drowsy','_inter'),)
  input<-unlist(paste0('v_dmn',j),)
  new_alert_df<-assign(output_name_a,subset(alert_scans,alert_scans$parcels==input))
  new_drowsy_df<-assign(output_name_d,subset(drowsy_scans,drowsy_scans$parcels==input))
  
  for (i in 1:100) {
    output_a<-subset(new_alert_df,new_alert_df$iterations==i)
    output_d<-subset(new_drowsy_df,new_drowsy_df$iterations==i)
    df_vdmn_cohenD<-cohensD(as.numeric(output_a$value),as.numeric(output_d$value))
    df_vdmn_cohenD_list[j,i]<-df_vdmn_cohenD
    df_vdmn_cohenD_list[is.na(df_vdmn_cohenD_list)] <- 0}}

df_vdmn_cohenD_list<-t(df_vdmn_cohenD_list)
colnames(df_vdmn_cohenD_list)<-c("v_dmn1","v_dmn2","v_dmn3","v_dmn4","v_dmn5","v_dmn6","v_dmn7","v_dmn8","v_dmn9","v_dmn10")


#Compute null pvalue from computing number of experimental cohenD values with null cohenD values

#df_raw<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_switching_rate_labeled.csv")
df_raw<-read.csv('/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_parcels_switching_rate_labeled.csv')

df_alert_raw<-subset(df_raw,df_raw$Arousal_State=="alert")
df_drowsy_raw<-subset(df_raw,df_raw$Arousal_State=="drowsy")

psal_pvalue_list <- data.frame(matrix(ncol=12, nrow=1))
for (i in 1:12) {
  input<-paste0('p_sal',i)
  alert<-select(df_alert_raw,input)
  drowsy<-select(df_drowsy_raw,input)
  output_name<-paste0('experimental_df_psal',i,'_cohenD')
  experiment_cohenD<-assign(output_name,cohensD(alert[,1],drowsy[,1]))
  psal_pvalue_list[i]<-sum(as.integer(as.logical(df_psal_cohenD_list[,i]>=experiment_cohenD)))/100
}

colnames(psal_pvalue_list)<-c("p_sal1","p_sal2","p_sal3","p_sal4","p_sal5","p_sal6","p_sal7","p_sal8","p_sal9","p_sal10","p_sal11","p_sal12")

asal_pvalue_list<-data.frame(matrix(ncol=7,nrow=1))
for (i in 1:7) {
  input<-paste0('a_sal',i)
  alert<-select(df_alert_raw,input)
  drowsy<-select(df_drowsy_raw,input)
  output_name<-paste0('experimental_df_asal',i,'_cohenD')
  assign(output_name,cohensD(alert[,1],drowsy[,1]))
  experiment_cohenD<-assign(output_name,cohensD(alert[,1],drowsy[,1]))
  asal_pvalue_list[i]<-sum(as.integer(as.logical(df_asal_cohenD_list[,i]>=experiment_cohenD)))/100
}

colnames(asal_pvalue_list)<-c("a_sal1","a_sal2","a_sal3","a_sal4","a_sal5","a_sal6","a_sal7")
rownames(asal_pvalue_list)<-c("pvalue")

lcen_pvalue_list<-data.frame(matrix(ncol=6,nrow=1))
for (i in 1:6) {
  input<-paste0('lcen_',i)
  alert<-select(df_alert_raw,input)
  drowsy<-select(df_drowsy_raw,input)
  output_name<-paste0('experimental_df_lcen',i,'_cohenD')
  assign(output_name,cohensD(alert[,1],drowsy[,1]))
  experiment_cohenD<-assign(output_name,cohensD(alert[,1],drowsy[,1]))
  lcen_pvalue_list[i]<-sum(as.integer(as.logical(df_lcen_cohenD_list[,i]>=experiment_cohenD)))/100
}
colnames(lcen_pvalue_list)<-c("lcen_1","lcen_2","lcen_3","lcen_4","lcen_5","lcen_6")
rownames(lcen_pvalue_list)<-c("pvalue")

rcen_pvalue_list<-data.frame(matrix(ncol=6,nrow=1))
for (i in 1:6) {
  input<-paste0('rcen_',i)
  alert<-select(df_alert_raw,input)
  drowsy<-select(df_drowsy_raw,input)
  output_name<-paste0('experimental_df_rcen',i,'_cohenD')
  assign(output_name,cohensD(alert[,1],drowsy[,1]))
  experiment_cohenD<-assign(output_name,cohensD(alert[,1],drowsy[,1]))
  rcen_pvalue_list[i]<-sum(as.integer(as.logical(df_rcen_cohenD_list[,i]>=experiment_cohenD)))/100
}
colnames(rcen_pvalue_list)<-c("rcen_1","rcen_2","rcen_3","rcen_4","rcen_5","rcen_6")
rownames(rcen_pvalue_list)<-c("pvalue")

ddmn_pvalue_list<-data.frame(matrix(ncol=9,nrow=1))
for (i in 1:9) {
  input<-paste0('d_dmn',i)
  alert<-select(df_alert_raw,input)
  drowsy<-select(df_drowsy_raw,input)
  output_name<-paste0('experimental_df_ddmn',i,'_cohenD')
  assign(output_name,cohensD(alert[,1],drowsy[,1]))
  experiment_cohenD<-assign(output_name,cohensD(alert[,1],drowsy[,1]))
  ddmn_pvalue_list[i]<-sum(as.integer(as.logical(df_ddmn_cohenD_list[,i]>=experiment_cohenD)))/100
}
colnames(ddmn_pvalue_list)<-c("d_dmn1","d_dmn2","d_dmn3","d_dmn4","d_dmn5","d_dmn6","d_dmn7","d_dmn8","d_dmn9")
rownames(ddmn_pvalue_list)<-c("pvalue")

vdmn_pvalue_list<-data.frame(matrix(ncol=10,nrow=1))
for (i in 1:10) {
  input<-paste0('v_dmn',i)
  alert<-select(df_alert_raw,input)
  drowsy<-select(df_drowsy_raw,input)
  output_name<-paste0('experimental_df_vdmn',i,'_cohenD')
  assign(output_name,cohensD(alert[,1],drowsy[,1]))
  experiment_cohenD<-assign(output_name,cohensD(alert[,1],drowsy[,1]))
  vdmn_pvalue_list[i]<-sum(as.integer(as.logical(df_vdmn_cohenD_list[,i]>=experiment_cohenD)))/100}
colnames(vdmn_pvalue_list)<-c("v_dmn1","v_dmn2","v_dmn3","v_dmn4","v_dmn5","v_dmn6","v_dmn7","v_dmn8","v_dmn9","v_dmn10")
rownames(vdmn_pvalue_list)<-c("pvalue")

return(list(asal_pvalue_list=asal_pvalue_list,psal_pvalue_list=psal_pvalue_list,lcen_pvalue_list=lcen_pvalue_list,
            rcen_pvalue_list=rcen_pvalue_list,ddmn_pvalue_list=ddmn_pvalue_list,vdmn_pvalue_list=vdmn_pvalue_list))
}

df_null_1<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_parcels_100_null_switching_rate_results.csv",0)
FIND_null_1_pval<-FIND_compute_null_pvalue(df_null_1)
FIND_null_1_pval<-as.data.frame((FIND_null_1_pval))
FIND_null_1_pval<-as.data.frame(t(FIND_null_1_pval))

#list.save(FIND_null_1_pval,'/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_parcels_null.rdata')

write.csv(FIND_null_1_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_model_1_results.csv")

#df_null_2<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_parcels_100_null_gs_switching_rate_results.csv",0)

df_null_2<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_parcels_gs_corrected_switching_rate_null.csv",0)
FIND_null_2_pval<-FIND_compute_null_pvalue(df_null_2)
FIND_null_2_pval<-as.data.frame((FIND_null_2_pval))
FIND_null_2_pval<-as.data.frame(t(FIND_null_2_pval))

list.save(FIND_null_2_pval,'/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_parcels_null_gs.rdata')

write.csv(FIND_null_2_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_gs_model_results.csv")

df_null_3<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_parcels_100_null_var_switching_rate_results.csv",0)
FIND_null_3_pval<-FIND_compute_null_pvalue(df_null_3)
FIND_null_3_pval<-as.data.frame((FIND_null_3_pval))
FIND_null_3_pval<-as.data.frame(t(FIND_null_3_pval))
#list.save(FIND_null_3_pval,'/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_parcels_null_var.rdata')

write.csv(FIND_null_3_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_var_model_results.csv")

df_null_4<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_parcels_stat_corr_corrected_switching_rate_null.csv",0)

FIND_null_4_pval<-FIND_compute_null_pvalue(df_null_4)

FIND_null_4_pval<-as.data.frame((FIND_null_4_pval))
FIND_null_4_pval<-as.data.frame(t(FIND_null_4_pval))

write.csv(FIND_null_4_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_stat_cor_model_results.csv")

  


