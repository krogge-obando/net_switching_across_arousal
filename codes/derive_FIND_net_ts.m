%This code will derive vu-eeg fMRI HC subject specific ts and map of FIND atlas

%Author Kim Kundert-Obando

%write to clear the code 
clear; clc; close all; 

%paths
main_path = '/data1/neurdylab/datasets/eegfmri_vu/PROC/';
save_dir='/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/FIND_atlas_ts'

% mask out the skull (nonzero voxels)
mni_mask = niftiread('/data1/neurdylab/MNI152_T1_2mm_brain_mask_filled.nii.gz');
brainVox = find(mni_mask>0);

%Calling in the FIND-LAB network maps
a_sal_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/anterior_Salience/anterior_Salience.nii.gz');
p_sal_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/post_Salience/post_Salience.nii.gz');
d_dmn_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/dorsal_DMN/dDMN.nii.gz');
v_dmn_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/ventral_DMN/vDMN.nii.gz');
l_cen_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/LECN/LECN.nii.gz');
r_cen_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/RECN/RECN.nii.gz');
aud_map   = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Auditory/Auditory.nii.gz');
bas_g_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Basal_Ganglia/Basal_Ganglia.nii.gz');
h_vis_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/high_Visual/high_Visual.nii.gz');
lang_map  = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Language/Language.nii.gz');
prec_map  = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Precuneus/Precuneus.nii.gz');
p_vis_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/prim_Visual/prim_Visual.nii.gz');
sens_map  = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Sensorimotor/Sensorimotor.nii.gz');
vis_s_map = niftiread('/data1/neurdylab/kim_ro/tri_net_proj/Functional_ROIs/Visuospatial/Visuospatial.nii.gz');

%this is the text file that will load in all the subjects name
scan= readtable('/data1/neurdylab/kim_ro/tri_net_proj/redo_proj/scan_list.txt', 'readvariablenames', 0);

allsub_corr_list=[];
%%
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
    subj_V_mask = subj_V2(brainVox,:);
    Y = subj_V_mask';

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

  % line up networks in a matrix:
    Xnets_0 = [a_sal_map(brainVox), p_sal_map(brainVox), d_dmn_map(brainVox),...
    v_dmn_map(brainVox),l_cen_map(brainVox),r_cen_map(brainVox),aud_map(brainVox),...
    bas_g_map(brainVox),h_vis_map(brainVox),lang_map(brainVox),prec_map(brainVox),...
    p_vis_map(brainVox),sens_map(brainVox),vis_s_map(brainVox)];
    
    Xnets = [ones(size(Xnets_0,1),1),  zscore(single(Xnets_0))];
    % spatial regression against this subject's fMRI data
    %  -> results in network time courses
    
    betas_r1 = pinv(Xnets)*Y_hm_removed';
    % the rows of betas_r1 correspond to the network time courses
    % the first row should be ignored (mean offset)
    betas_r1=zscore(betas_r1');

    scan_ts=betas_r1(:,2:15);

    save_path=[save_dir,'/',this_subject,'_FIND_net_ts.mat']

    save(save_path,"scan_ts")
   
end
