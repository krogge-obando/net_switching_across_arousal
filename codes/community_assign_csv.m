%this code will derive a csv file for all the community assignments for the
%data

%load alert scan names
%%
scan_list= readtable('/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/scan_list.txt', 'readvariablenames', 0);


data_path=['/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/net_flex/']
save_path=['/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/community_assignment/']
%%
for scan=1:height(scan_list)

    subject_name = scan_list{scan,1}{1};
    scan_number  = scan_list{scan,2}{1};
    this_subject =[subject_name,'_',scan_number]

    file_path=append(data_path,this_subject,"FIND_net_flex.mat")
    
    load(file_path)
    
    community_assignment=node_assignment_by_window_sub{1,1};

    
    save_file_name=append(save_path,this_subject,'_community.csv');

%set names of csv with 

this_subject_list(scan,:)=[this_subject];

writematrix(community_assignment,save_file_name)

end
%%
%save alert and drowsy name
save_list=append(save_path,'scan_list.csv')
writematrix(this_subject_list,save_list)


%save csv alert or drowsy

