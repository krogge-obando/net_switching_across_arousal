#This code will recompute the FIND networks switching rate null model p-value

#Author Kim Rogge-Obando
# Updated: 04/16/2024

install.packages("lsr")
library(lsr)
library(reshape2)



#writefunction

FIND_compute_null_pvalue<-function(df_null){
  
  rownames(df_null)<-c('vcon02_scan01','vcon02_scan02','vcon03_scan01','vcon04_scan02','vcon05_scan01','vcon07_scan01','vcon08_scan01',
                            'vcon08_scan02', 'vcon09_scan01','vcon09_scan02','vcon10_scan01','vcon10_scan02','vcon11_scan02','vcon12_scan01',
                            'vcon13_scan01','vcon14_scan01','vcon14_scan02','vcon15_scan01','vcon15_scan02','vcon17_scan01','vcon17_scan02',
                            'vcon18_scan01','vcon18_scan02','vcon19_scan01','vcon19_scan02','vcon20_scan01','vcon20_scan02','vcon22_scan01',
                            'vcon22_scan02','vcon23_scan01', 'vcon24_scan01','vcon25_scan01','vcon25_scan02', 'vcon26_scan01','vcon27_scan01',
                            'vcon27_scan02','vcon29_scan01','vcon29_scan02','vcon30_scan01','vcon32_scan01','vcon33_scan01','vcon33_scan02',
                            'vcon35_scan01','vcon36_scan01','vcon37_scan01')

#write.csv(df_null,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_null_data_raw.csv")

#df_null <- read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_null_data_raw.csv")

iterations<-rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
                  24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
                  44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,
                  64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
                  84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100),times=14)


networks<-rep(c("a_sal","p_sal","d_dmn","v_dmn","l_cen","r_cen","aud","bas_g",
                "h_vis","lang","prec","p_vis","sens","vis_s"),times=100);

networks<-as.data.frame(networks)

networks<-networks[order(factor(networks$networks,levels=c(c("a_sal","p_sal","d_dmn","v_dmn","l_cen","r_cen","aud","bas_g",
                                                             "h_vis","lang","prec","p_vis","sens","vis_s")))),]
df_null_labled<-rbind(iterations,df_null)
df_null_labled<-rbind(networks,df_null_labled)


rownames(df_null_labled)<-c('networks','iterations','vcon02_scan01','vcon02_scan02','vcon03_scan01','vcon04_scan02','vcon05_scan01','vcon07_scan01','vcon08_scan01',
                     'vcon08_scan02', 'vcon09_scan01','vcon09_scan02','vcon10_scan01','vcon10_scan02','vcon11_scan02','vcon12_scan01',
                     'vcon13_scan01','vcon14_scan01','vcon14_scan02','vcon15_scan01','vcon15_scan02','vcon17_scan01','vcon17_scan02',
                     'vcon18_scan01','vcon18_scan02','vcon19_scan01','vcon19_scan02','vcon20_scan01','vcon20_scan02','vcon22_scan01',
                     'vcon22_scan02','vcon23_scan01', 'vcon24_scan01','vcon25_scan01','vcon25_scan02', 'vcon26_scan01','vcon27_scan01',
                     'vcon27_scan02','vcon29_scan01','vcon29_scan02','vcon30_scan01','vcon32_scan01','vcon33_scan01','vcon33_scan02',
                     'vcon35_scan01','vcon36_scan01','vcon37_scan01')

#colnames(df_null_labled)<-networks
alert_scans<-df_null_labled[c("networks","iterations","vcon02_scan01","vcon07_scan01","vcon15_scan01","vcon15_scan02","vcon17_scan02","vcon19_scan02","vcon20_scan01", "vcon20_scan02","vcon29_scan01","vcon29_scan02"),]
drowsy_scans<-df_null_labled[c("networks","iterations", "vcon05_scan01", "vcon08_scan01", "vcon08_scan02", "vcon09_scan01", "vcon09_scan02", "vcon12_scan01", "vcon14_scan01",
                              "vcon22_scan01", "vcon23_scan01", "vcon24_scan01", "vcon25_scan02", "vcon32_scan01", "vcon36_scan01", "vcon37_scan01"),]

#generating df_of each network for alert and drowsy

alert_scans<-as.data.frame(t(alert_scans))
drowsy_scans<-as.data.frame(t(drowsy_scans))

#

alert_scans<-melt(alert_scans,id=c("networks","iterations"))
drowsy_scans<-melt(drowsy_scans,id=c("networks","iterations"))

#generating df_of each network for alert and drowsy
df_psal_alert_iter<-subset(alert_scans,alert_scans$networks=="p_sal")
df_asal_alert_iter<-subset(alert_scans,alert_scans$networks=="a_sal")
df_lcen_alert_iter<-subset(alert_scans,alert_scans$networks=="l_cen")
df_rcen_alert_iter<-subset(alert_scans,alert_scans$networks=="r_cen")
df_ddmn_alert_iter<-subset(alert_scans,alert_scans$networks=="d_dmn")
df_vdmn_alert_iter<-subset(alert_scans,alert_scans$networks=="v_dmn")

df_psal_drowsy_iter<-subset(drowsy_scans,drowsy_scans$networks=="p_sal")
df_asal_drowsy_iter<-subset(drowsy_scans,drowsy_scans$networks=="a_sal")
df_lcen_drowsy_iter<-subset(drowsy_scans,drowsy_scans$networks=="l_cen")
df_rcen_drowsy_iter<-subset(drowsy_scans,drowsy_scans$networks=="r_cen")
df_ddmn_drowsy_iter<-subset(drowsy_scans,drowsy_scans$networks=="d_dmn")
df_vdmn_drowsy_iter<-subset(drowsy_scans,drowsy_scans$networks=="v_dmn")

#forloop to generate 100 cohen D values of psal
df_psal_cohenD_list <- data.frame(matrix(ncol=100, nrow=1))

for (i in 1:100) {
  df_psal_alert_each_iter<-subset(df_psal_alert_iter,df_psal_alert_iter$iteration==i)
  df_psal_drowsy_each_iter<-subset(df_psal_drowsy_iter,df_psal_drowsy_iter$iteration==i)
  df_psal_cohenD<-cohensD(as.numeric(df_psal_alert_each_iter$value),as.numeric(df_psal_drowsy_each_iter$value))
  df_psal_cohenD_list[i]<-df_psal_cohenD}
#forloop to generate 100 cohen D values of asal
df_asal_cohenD_list <- data.frame(matrix(ncol=100, nrow=1))
for (i in 1:100) {
  df_asal_alert_each_iter<-subset(df_asal_alert_iter,df_asal_alert_iter$iteration==i)
  df_asal_drowsy_each_iter<-subset(df_asal_drowsy_iter,df_asal_drowsy_iter$iteration==i)
  df_asal_cohenD<-cohensD(as.numeric(df_asal_alert_each_iter$value),as.numeric(df_asal_drowsy_each_iter$value))
  df_asal_cohenD_list[i]<-df_asal_cohenD}

#forloop to generate 100 cohen D values of lcen
df_lcen_cohenD_list <- data.frame(matrix(ncol=100, nrow=1))
for (i in 1:100) {
  df_lcen_alert_each_iter<-subset(df_lcen_alert_iter,df_lcen_alert_iter$iteration==i)
  df_lcen_drowsy_each_iter<-subset(df_lcen_drowsy_iter,df_lcen_drowsy_iter$iteration==i)
  df_lcen_cohenD<-cohensD(as.numeric(df_lcen_alert_each_iter$value),as.numeric(df_lcen_drowsy_each_iter$value))
  df_lcen_cohenD_list[i]<-df_lcen_cohenD}

#forloop to generate 100 cohen D values of rcen
df_rcen_cohenD_list <- data.frame(matrix(ncol=100, nrow=1))
for (i in 1:100) {
  df_rcen_alert_each_iter<-subset(df_rcen_alert_iter,df_rcen_alert_iter$iteration==i)
  df_rcen_drowsy_each_iter<-subset(df_rcen_drowsy_iter,df_rcen_drowsy_iter$iteration==i)
  df_rcen_cohenD<-cohensD(as.numeric(df_rcen_alert_each_iter$value),as.numeric(df_rcen_drowsy_each_iter$value))
  df_rcen_cohenD_list[i]<-df_rcen_cohenD}

#forloop to generate 100 cohen D values of ddmn
df_ddmn_cohenD_list <- data.frame(matrix(ncol=100, nrow=1))
for (i in 1:100) {
  df_ddmn_alert_each_iter<-subset(df_ddmn_alert_iter,df_ddmn_alert_iter$iteration==i)
  df_ddmn_drowsy_each_iter<-subset(df_ddmn_drowsy_iter,df_ddmn_drowsy_iter$iteration==i)
  df_ddmn_cohenD<-cohensD(as.numeric(df_ddmn_alert_each_iter$value),as.numeric(df_ddmn_drowsy_each_iter$value))
  df_ddmn_cohenD_list[i]<-df_ddmn_cohenD}

#forloop to generate 100 cohen D values of vdmn
df_vdmn_cohenD_list <- data.frame(matrix(ncol=100, nrow=1))
for (i in 1:100) {
  df_vdmn_alert_each_iter<-subset(df_vdmn_alert_iter,df_vdmn_alert_iter$iteration==i)
  df_vdmn_drowsy_each_iter<-subset(df_vdmn_drowsy_iter,df_vdmn_drowsy_iter$iteration==i)
  df_vdmn_cohenD<-cohensD(as.numeric(df_vdmn_alert_each_iter$value),as.numeric(df_vdmn_drowsy_each_iter$value))
  df_vdmn_cohenD_list[i]<-df_vdmn_cohenD}

#Compute experimental cohenD values

df_raw<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_switching_rate_labeled.csv")

df_alert_raw<-subset(df_raw,df_raw$Vsleep_group2=="alert")
df_drowsy_raw<-subset(df_raw,df_raw$Vsleep_group2=="drowsy")

experimental_df_psal_cohenD <-cohensD(df_alert_raw$p_sal,df_drowsy_raw$p_sal)
experimental_df_asal_cohenD <-cohensD(df_alert_raw$a_sal,df_drowsy_raw$a_sal)
experimental_df_lcen_cohenD <-cohensD(df_alert_raw$l_cen,df_drowsy_raw$l_cen)
experimental_df_rcen_cohenD <-cohensD(df_alert_raw$p_sal,df_drowsy_raw$r_cen)
experimental_df_ddmn_cohenD <-cohensD(df_alert_raw$d_dmn,df_drowsy_raw$d_dmn)
experimental_df_vdmn_cohenD <-cohensD(df_alert_raw$v_dmn,df_drowsy_raw$v_dmn)

psal_p_value<-sum(as.integer(as.logical(df_psal_cohenD_list>=experimental_df_psal_cohenD)))/100
asal_p_value<-sum(as.integer(as.logical(df_asal_cohenD_list>=experimental_df_asal_cohenD)))/100
lcen_p_value<-sum(as.integer(as.logical(df_lcen_cohenD_list>=experimental_df_lcen_cohenD)))/100
rcen_p_value<-sum(as.integer(as.logical(df_rcen_cohenD_list>=experimental_df_rcen_cohenD)))/100
ddmn_p_value<-sum(as.integer(as.logical(df_ddmn_cohenD_list>=experimental_df_ddmn_cohenD)))/100
vdmn_p_value<-sum(as.integer(as.logical(df_vdmn_cohenD_list>=experimental_df_vdmn_cohenD)))/100

null_model_p_value_results <- data.frame(network=c("psal_p_value","asal_p_value","lcen_p_value","rcen_p_value","ddmn_p_value","vdmn_p_value"),
                                         value=c(psal_p_value,
                                                 asal_p_value,
                                                 lcen_p_value,
                                                 rcen_p_value,
                                                 ddmn_p_value,
                                                 vdmn_p_value)) }



df_null_1<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_100_null_switching_rate_results.csv",0)
FIND_null_1_pval<-FIND_compute_null_pvalue(df_null_1)
write.csv(FIND_null_1_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_model.csv")

#df_null_2<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_100_null_gs_switching_rate_results.csv",0)

#This was done with the corrected values
df_null_2<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_NET_gs_corrected_switching_rate_null.csv",0)

FIND_null_2_pval<-FIND_compute_null_pvalue(df_null_2)
write.csv(FIND_null_2_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_gs_model.csv")

df_null_3<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_100_null_var_switching_rate_results.csv",0)
FIND_null_3_pval<-FIND_compute_null_pvalue(df_null_3)
write.csv(FIND_null_3_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_var_model.csv")

df_null_4<-read.csv("/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_networks_stat_corr_switching_rate_null.csv",0)
FIND_null_4_pval<-FIND_compute_null_pvalue(df_null_4)
write.csv(FIND_null_4_pval,"/Users/roggeokk/Desktop/Projects/Vigilance_SwitchingRate_2/data/FIND_atlas_networks_null_stat_corr_model.csv")

