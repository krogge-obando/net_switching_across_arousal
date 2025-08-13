%This code will derive vu-eeg fMRI HC subject specific ts and map of FIND
%atlas parcels

%Author Kim Kundert-Obando

%write to clear the code 
clear; clc; close all; 

%write out all paths
main_path = '/data1/neurdylab/datasets/eegfmri_vu/PROC/';
save_dir='/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/FIND_atlas_ts';

%this is the text file that will load in all the subjects name
scan= readtable('/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/scan_list.txt', 'readvariablenames', 0);

allsub_corr_list=[];

for subject = 1:height(scan)

    subject_name = scan{subject,1}{1};
    scan_number  = scan{subject,2}{1};
    this_subject =[subject_name,'_',scan_number]
  
    % load fMRI file names
    nii_name = [subject_name,'-',scan_number,'_EPI2MNI_sm_nr.nii.gz'];
    subj_dir=[main_path,subject_name,'/meica_proc_',scan_number,'/meica_out/',nii_name];
    subj_fmri = niftiread(subj_dir); 

    %reformat the fMRI file
    dims = size(subj_fmri);
    subj_V2 = reshape(subj_fmri,prod(dims(1:3)),size(subj_fmri,4));
    Y = subj_V2';

    mean_4D=mean(Y(:));

    Y_scaled=1000*Y/mean_4D;
   % V_mask_norm = zscore(subj_V_mask')';

    %Load motion parameters
    par_name = [subject_name,'-',scan_number,'_ecr_e2.volreg_par'];
    par_dir=[main_path,subject_name,'/meica_proc_',scan_number,'/echo2/',par_name];

    par_file  =load(par_dir);
    dMP=diff(par_file);
    dMP = [dMP(1,:); dMP];

    %remove motion parameters %hm is head motion
    X_hm = [ones(size(par_file,1),1), zscore(par_file),zscore(dMP)];
    betas_hm = pinv(X_hm) * Y_scaled;
    hm_confound = X_hm(:,2:end) * betas_hm(2:end,:);
    Y_hm_removed = Y_scaled - hm_confound; % this will be a new fmri with head motion removed

    Y_hm_removed = zscore(Y_hm_removed);
    Y_hm_removed = Y_hm_removed';
    
     for i=1:9
        d_dmn_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/dorsal_DMN/','0',num2str(i),'/',num2str(i),'.nii'];
        d_dmn_mask=niftiread(d_dmn_mask_path);
        d_dmn_reshape = reshape(d_dmn_mask,[],1);
        d_dmn_parcel=find(d_dmn_reshape>0);
        d_dmn_ts = mean(Y_hm_removed(d_dmn_parcel,:),1);
        d_dmn_ts=zscore(d_dmn_ts);
        d_dmn_ts_parcel(i,:)=[d_dmn_ts];
    end

    for i=1:10
        if i<10
            v_dmn_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/ventral_DMN/','0',num2str(i),'/',num2str(i),'.nii'];

        else
            v_dmn_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/ventral_DMN/',num2str(i),'/',num2str(i),'.nii'];
        end
        v_dmn_mask=niftiread(v_dmn_mask_path);
        v_dmn_reshape = reshape(v_dmn_mask,[],1);
        v_dmn_parcel=find(v_dmn_reshape>0);
        v_dmn_ts = mean(Y_hm_removed(v_dmn_parcel,:),1);
        v_dmn_ts=zscore(v_dmn_ts);
        v_dmn_ts_parcel(i,:)=[v_dmn_ts];
    end

    for i=1:7
        a_sal_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/anterior_Salience/','0',num2str(i),'/',num2str(i),'.nii'];
        a_sal_mask=niftiread(a_sal_mask_path);
        a_sal_reshape = reshape(a_sal_mask,[],1);
        a_sal_parcel=find(a_sal_reshape>0);
        a_sal_ts = mean(Y_hm_removed(a_sal_parcel,:),1);
        a_sal_ts = zscore(a_sal_ts);
        a_sal_ts_parcel(i,:)=[a_sal_ts];
    end

    for i=1:12
        if i<10
            p_sal_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/post_Salience/','0',num2str(i),'/',num2str(i),'.nii'];

        else
            p_sal_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/post_Salience/',num2str(i),'/',num2str(i),'.nii'];
        end

        p_sal_mask=niftiread(p_sal_mask_path);
        p_sal_reshape = reshape(p_sal_mask,[],1);
        p_sal_parcel=find(p_sal_reshape>0);
        p_sal_ts = mean(Y_hm_removed(p_sal_parcel,:),1);
        p_sal_ts = zscore(p_sal_ts);
        p_sal_ts_parcel(i,:)=[p_sal_ts];
    end

    for i=1:6
        lcen_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/LECN/','0',num2str(i),'/',num2str(i),'.nii'];
        lcen_mask=niftiread(lcen_mask_path);
        lcen_reshape = reshape(lcen_mask,[],1);
        lcen_parcel=find(lcen_reshape>0);
        lcen_ts = mean(Y_hm_removed(lcen_parcel,:),1);
        lcen_ts = zscore(lcen_ts);
        lcen_ts_parcel(i,:)=[lcen_ts];

    end

    for i=1:6
        rcen_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/RECN/','0',num2str(i),'/',num2str(i),'.nii'];
        rcen_mask=niftiread(rcen_mask_path);
        rcen_reshape = reshape(rcen_mask,[],1);
        rcen_parcel=find(rcen_reshape>0);
        rcen_ts = mean(Y_hm_removed(rcen_parcel,:),1);
        rcen_ts = zscore(rcen_ts);
        rcen_ts_parcel(i,:)=[rcen_ts];
    end

    for i=1:3
        aud_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Auditory/','0',num2str(i),'/',num2str(i),'.nii'];
        aud_mask=niftiread(aud_mask_path);
        aud_reshape = reshape(aud_mask,[],1);
        aud_parcel=find(aud_reshape>0);
        aud_ts = mean(Y_hm_removed(aud_parcel,:),1);
        aud_ts = zscore(aud_ts);
        aud_ts_parcel(i,:)=[aud_ts];
    end

    for i=1:5
        basg_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Basal_Ganglia/','0',num2str(i),'/',num2str(i),'.nii'];
        basg_mask=niftiread(basg_mask_path);
        basg_reshape = reshape(basg_mask,[],1);
        basg_parcel=find(basg_reshape>0);
        basg_ts = mean(Y_hm_removed(basg_parcel,:),1);
        basg_ts = zscore(basg_ts);
        basg_ts_parcel(i,:)=[basg_ts];

    end

    for i=1:2
        pVis_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/prim_Visual/','0',num2str(i),'/',num2str(i),'.nii'];
        pVis_mask=niftiread(pVis_mask_path);
        pVis_reshape = reshape(pVis_mask,[],1);
        pVis_parcel=find(pVis_reshape>0);
        pVis_ts = mean(Y_hm_removed(pVis_parcel,:),1);
        pVis_ts = zscore(pVis_ts);
        pVis_ts_parcel(i,:)=[pVis_ts];
    end

    for i=1:2
        hVis_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/high_Visual/','0',num2str(i),'/',num2str(i),'.nii'];
        hVis_mask=niftiread(hVis_mask_path);
        hVis_reshape = reshape(hVis_mask,[],1);
        hVis_parcel=find(hVis_reshape>0);
        hVis_ts = mean(Y_hm_removed(hVis_parcel,:),1);
        hVis_ts = zscore(hVis_ts);
        hVis_ts_parcel(i,:)=[hVis_ts];
    end

    for i=1:7
        lang_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Language/','0',num2str(i),'/',num2str(i),'.nii'];
        lang_mask=niftiread(lang_mask_path);
        lang_reshape = reshape(lang_mask,[],1);
        lang_parcel=find(lang_reshape>0);
        lang_ts = mean(Y_hm_removed(lang_parcel,:),1);
        lang_ts = zscore(lang_ts);
        lang_ts_parcel(i,:)=[lang_ts];
    end

    for i=1:4
        prec_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Precuneus/','0',num2str(i),'/',num2str(i),'.nii'];
        prec_mask=niftiread(prec_mask_path);
        prec_reshape = reshape(prec_mask,[],1);
        prec_parcel=find(prec_reshape>0);
        prec_ts = mean(Y_hm_removed(prec_parcel,:),1);
        prec_ts = zscore(prec_ts);
        prec_ts_parcel(i,:)=[prec_ts];
    end

    for i=1:6
        sensmotor_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Sensorimotor/','0',num2str(i),'/',num2str(i),'.nii'];
        sensmotor_mask=niftiread(sensmotor_mask_path);
        sensmotor_reshape = reshape(sensmotor_mask,[],1);
        sensmotor_parcel=find(sensmotor_reshape>0);
        sensmotor_ts = mean(Y_hm_removed(sensmotor_parcel,:),1);
        sensmotor_ts = zscore(sensmotor_ts);
        sensmotor_ts_parcel(i,:)=[sensmotor_ts];
    end

    for i=1:11
        if i<10
            vis_spat_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Visuospatial/','0',num2str(i),'/',num2str(i),'.nii'];

        else
            vis_spat_mask_path =['/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Visuospatial/',num2str(i),'/',num2str(i),'.nii'];
        end
        vis_spat_mask=niftiread(vis_spat_mask_path);
        vis_spat_reshape = reshape(vis_spat_mask,[],1);
        vis_spat_parcel=find(vis_spat_reshape>0);
        vis_spat_ts = mean(Y_hm_removed(vis_spat_parcel,:),1);
        vis_spat_ts = zscore(vis_spat_ts);
        vis_spat_ts_parcel(i,:)=[vis_spat_ts];
    end

    networks_ts=[d_dmn_ts_parcel;v_dmn_ts_parcel;a_sal_ts_parcel;...
        p_sal_ts_parcel;lcen_ts_parcel;rcen_ts_parcel; aud_ts_parcel;...
        basg_ts_parcel; pVis_ts_parcel;hVis_ts_parcel; lang_ts_parcel;...
        prec_ts_parcel; sensmotor_ts_parcel; vis_spat_ts_parcel];
    
    save_path=[save_dir,'/',this_subject,'_FIND_parcel_ts.mat'];

    save(save_path,"networks_ts")
    
end
