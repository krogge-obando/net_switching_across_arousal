%This code is an optimized version to compute FIND parcels null model

%Author:Kim Rogge-Obando k.rogge.obando@gmail.com
%add paths

addpath GenLouvain-master/HelperFunctions/
addpath GenLouvain-master/

save_path=["/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/net_flex/"];

scan_list= readtable('/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/scan_list.txt', 'readvariablenames', 0);

order_of_scan_list=[];


%%
for scan = 1:height(scan_list)

    subject_name = scan_list{scan,1}{1};
    scan_number  = scan_list{scan,2}{1};
    this_subject =[subject_name,'_',scan_number]

    scan_ts_file=['/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/FIND_atlas_ts/',this_subject,'_FIND_net_ts.mat'];
 
    load(scan_ts_file)
   
    [node_assignment_by_window_sub,number_of_communities_new...
       ,number_of_communities_sub,node_switching_sub] = compute_switching_rate_5_23_25(scan_ts,72,1);
    
    order_of_scan_list(scan,:)=[this_subject];

    full_file_path=strcat(save_path,this_subject,"FIND_net_flex.mat");

    save(full_file_path,"node_switching_sub","number_of_communities_new","node_assignment_by_window_sub","number_of_communities_sub")

    node_switching_list(scan,:)=[node_switching_sub];

end

eeg_fmri_vu_flexibility_list=strcat(save_path,"eeg_fmri_vu_HC_net_flex_list.mat");

save(eeg_fmri_vu_flexibility_list,"node_switching_list","order_of_scan_list")

%Complete for parcels
%% To run this section clear all and run first section
for scan = 1:height(scan_list)

    subject_name = scan_list{scan,1}{1};
    scan_number  = scan_list{scan,2}{1};
    this_subject =[subject_name,'_',scan_number]

    scan_ts_file=['/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/FIND_atlas_ts/',this_subject,'_FIND_parcel_ts.mat'];
 
    load(scan_ts_file)

    networks_ts=networks_ts';
   
    [node_assignment_by_window_sub,number_of_communities_new...
       ,number_of_communities_sub,node_switching_sub] = compute_switching_rate_5_23_25(networks_ts,72,1);
    
    order_of_scan_list(scan,:)=[this_subject];

    full_file_path=strcat(save_path,this_subject,"FIND_parcel_flex.mat");

    save(full_file_path,"node_switching_sub","number_of_communities_new","node_assignment_by_window_sub","number_of_communities_sub")

    node_switching_list(scan,:)=[node_switching_sub];

end

eeg_fmri_vu_flexibility_list=strcat(save_path,"eeg_fmri_vu_HC_parcels_flex_list.mat");

save(eeg_fmri_vu_flexibility_list,"node_switching_list","order_of_scan_list")
