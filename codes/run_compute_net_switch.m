%This code is made to run a parloop to compute the null models 

%Author:Kim Rogge-Obando k.rogge.obando@gmail.com
%add paths
%%
addpath GenLouvain-master/HelperFunctions/
addpath GenLouvain-master/

%load data and prepare data to be something we can extract

load("FIND_parcels_sub_null_data_global_signal_hold.mat")


parfor scan=1:dims(1)
    
   scan
   scan_ts= scan_null_ts(scan,:,:,:);
   scan_ts=squeeze(scan_ts);
   
   [node_assignment_by_window_sub(scan,:,:),number_of_communities_new(scan,:,:)...
       ,number_of_communities_sub(scan,:,:),node_switching_sub(scan,:,:)] = compute_net_flex(scan_ts,72,100);

end


save_path="/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/null_net_flex/"
file_name="FIND_Parcels_net_flex_null_model_3_gs.mat"
file_path=append(save_path,file_name)

save(file_path,"node_switching_sub","number_of_communities_sub","number_of_communities_new","node_assignment_by_window_sub",'-v7.3')
