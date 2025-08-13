%This code will generate the null models for the VU-EEG fMRI data

addpath '/brain_benchmark_toolbox-master'

%write to clear the code 
clear; clc; close all; 

%paths

%this is the text file that will load in all the subjects name
scan= readtable('scan_list.txt', 'readvariablenames', 0);

allsub_corr_list=[];
%%

%% Generating null model 1 maintaining temporal and nodal mean

for scan_id = 1:height(scan)

    subject_name = scan{scan_id,1}{1};
    scan_number  = scan{scan_id,2}{1};
    this_subject =[subject_name,'_',scan_number]
    scan_ts_path=append(main_path, this_subject,'_FIND_parcel_ts.mat')
    load(scan_ts_path)


for i = 1:100
[ts_shuffled,conv_hist,Args,Cons]=bbt(networks_ts);
it_ts_shuffle(i,:,:)=[ts_shuffled];
it_conv_hist{i}=[conv_hist];
it_Args{i}=[Args];
it_Cons{i}=[Cons];
end
subj_100it_shuffle_list(scan_id,:,:,:)=[it_ts_shuffle];
subj_100it_conv_hist{scan_id}=[it_conv_hist];
subj_100it_Args{scan_id}=[it_Args];
subj_100it_Cons{scan_id}=[it_Cons];
end

save("FIND_parcels_sub_null_data_no_hold","subj_100it_shuffle_list")
save("FIND_parcels_sub_null_data_no_hold_parameters","subj_100it_shuffle_list","subj_100it_conv_hist","subj_100it_Args","subj_100it_Cons")
%%

%% Generate null model 2  while maintaining  nodal variation or time variation with "varnode" and "vartime"

for scan_id = 1:height(scan)
    subject_name = scan{scan_id,1}{1};
    scan_number  = scan{scan_id,2}{1};
    this_subject =[subject_name,'_',scan_number]
    scan_ts_path=append(main_path, this_subject,'_FIND_parcel_ts.mat')
    load(scan_ts_path)
for i = 1:100

[ts_shuffled, conv_hist, Args, Cons]=bbt(networks_ts,'varnode', true, 'vartime', true);
it_ts_shuffle(i,:,:)=[ts_shuffled];
it_conv_hist{i}=[conv_hist];
it_Args{i}=[Args];
it_Cons{i}=[Cons];
end
subj_100it_shuffle_list(scan_id,:,:,:)=[it_ts_shuffle];
subj_100it_conv_hist{scan_id}=[it_conv_hist];
subj_100it_Args{scan_id}=[it_Args];
subj_100it_Cons{scan_id}=[it_Cons];
end

save("FIND_parcels_sub_null_data_vartime_varnode_hold","subj_100it_shuffle_list")
save("FIND_parcels_sub_null_data_vartime_varnode_hold_parameters","subj_100it_shuffle_list","subj_100it_conv_hist","subj_100it_Args","subj_100it_Cons")

%% Null model 3 maining mean, variation and correlation with global signal 

for scan_id = 1:height(scan)
    subject_name = scan{scan_id,1}{1};
    scan_number  = scan{scan_id,2}{1};
    this_subject =[subject_name,'_',scan_number]
    scan_ts_path=append(main_path, this_subject,'_FIND_parcel_ts.mat')
    load(scan_ts_path)
for i = 1:100
    i

[ts_shuffled conv_hist,Args, Cons]=bbt(networks_ts,'covnodemode', true,'varnode', true, 'vartime', true);
it_ts_shuffle(i,:,:)=[ts_shuffled];
it_conv_hist{i}=[conv_hist];
it_Args{i}=[Args];
it_Cons{i}=[Cons];
end
subj_100it_shuffle_list(scan_id,:,:,:)=[it_ts_shuffle];
subj_100it_conv_hist{scan_id}=[it_conv_hist];
subj_100it_Args{scan_id}=[it_Args];
subj_100it_Cons{scan_id}=[it_Cons];
end

save("FIND_parcels_sub_null_data_global_signal_hold","subj_100it_shuffle_list")
save("FIND_parcels_sub_null_data_global_signal_hold_parameters.mat","subj_100it_conv_hist","subj_100it_Args","subj_100it_Cons")

%% Here I will generate null model 4 which preserves static correlations.

% FIND_NET_ORDER = [1 1 2 2 3 3 4 5 6 7 8 6 9 6]; %Note this was made with
% accordance to the FIND network time series matrix
% 
%     k = max(FIND_NET_ORDER);
%     networks_of_interest = [1 2 3];
%     network_constraint = false(k);
%     network_constraint(networks_of_interest, networks_of_interest) = 1;
% 

FIND_PARCEL_ORDER = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 5 5 5 5 5 6 6 6 6 7 7 7 7 7 7 7 8 8 8 8 9 9 9 9 9 9 6 6 6 6 6 6 6 6 6 6 6];

    k = max(FIND_PARCEL_ORDER);
    networks_of_interest = [1 2 3];
    network_constraint = false(k);
    network_constraint(networks_of_interest, networks_of_interest) = 1;

 

    for scan_id = 1:height(scan)
     
 subject_name = scan{scan_id,1}{1};
    scan_number  = scan{scan_id,2}{1};
    this_subject =[subject_name,'_',scan_number]
    scan_ts_path=append(main_path, this_subject,'_FIND_parcel_ts.mat')
    load(scan_ts_path)

        for i = 1:100           
            
            [ts_shuffled, conv_hist,Args, Cons]= bbt(networks_ts, 'covnodemode', true, 'varnode', ...
                true, 'vartime', true, 'partition', FIND_PARCEL_ORDER, 'covsystem', network_constraint);
            it_ts_shuffle(i,:,:)=[ts_shuffled];
            it_conv_hist{i}=[conv_hist];
            it_Args{i}=[Args];
            it_Cons{i}=[Cons];
        end
        subj_100it_shuffle_list(scan_id,:,:,:)=[it_ts_shuffle];
        subj_100it_conv_hist{scan_id}=[it_conv_hist];
        subj_100it_Args{scan_id}=[it_Args];
        subj_100it_Cons{scan_id}=[it_Cons];
    end

save("FIND_parcels_sub_null_data_stat_corr_hold","subj_100it_shuffle_list")
save("FIND_parcels_sub_null_data_stat_corr_hold_parameters.mat","subj_100it_conv_hist","subj_100it_Args","subj_100it_Cons")

