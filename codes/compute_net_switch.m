%This code can accomidate running the switching rate for a null ts or experiemental ts
    %for running the experimental ts you need to input iteration as 1.

    %EXTERNAL FUNCTIONS NEEDED TO RUN THIS CODE 

   % 'multilord', 'postprocess_ordinal_multilayer and 'iterated_genlouvain'
   % at https://github.com/GenLourvain/HenLouvain

   %'flexiblity' at http://commdetect.weebly.com/

    %INPUT INFORMATION
%scan ts: is the matrix of time series for null models this is a 3D matrix
%(null by node by TR) for experimental dat this is (node x TR)

%window_l : is the window length for the dynamic connectivity matrix. It is
%suggested to have a window length greater than one minute in length.

%iteration :is the number of nulls in your 3D matrix, put 1 if you are
%running an experimental ts.

%

function [node_assignment_by_window_sub,number_of_communities_new,number_of_communities_sub,...
    node_switching_sub] = compute_switching_rate(scan_ts,window_l,iteration)

for ts=1:iteration
    
    %Setting up the ts in its correct orientation
    null_network_ts=scan_ts(ts,:,:); 
    null_network_ts=squeeze(null_network_ts);
    node_ts=null_network_ts; %want node x TR
    
    %Setting up the information to conduct a sliding window analysis
    window_l= window_l;
    dims=size(nodes_ts);
    num_window=(dims(2)-window_l+1);
    window_l_t=window_l;
 
%Sliding window analysis
    for j = 1:num_window
        window_ts=nodes_ts(:,j:j+window_l_t-1);
        C0=corr(window_ts');
        C0(C0<0) = 0;
        C0(isnan(C0))=0;
        C0_w(:,:,j)=[C0];
    end % end window loop

    %Iteration for consensus matrix, we will conduct multilayer modularity
    %10 times and store the values of community assignment

    for s_i = 1:10

        for w=1:num_window
            A_sub =C0_w;
            A_n=A_sub(:,:,w);
            A_sub_n{w}=A_n;
        end
        
        %parameters for multilayer modularity typically left at 1 for both
        
        gamma=1;
        omega=1;


        [B1,mm] = multiord(A_sub_n',gamma,omega);
        % multilayer modularity set-up
        PP = @(S) postprocess_ordinal_multilayer(S,num_window);
        % S = modularity assignments; Q1 = modularity value; mod_iter_tmp = iterations
        [S,Q1,mod_iter_tmp] = iterated_genlouvain(B1,10000,0,1,'moverandw',[],PP);
        % 2D representation of modularity assignments
        S_new = reshape(S,dims(1),num_window);
        % saving the node assignment by window
        node_assignment_by_window{s_i} =S_new;
    end %end of consensus matrix iteration s_i

    %Computing consensus matrix
    s=10;

    for w_ss=1:num_window
        for i_ss=1:s
            K1_test=node_assignment_by_window{i_ss};
            K1_test2=K1_test(:,w_ss);
            K1(i_ss,:) =[K1_test2];
        end % end of i_ss creating K1 matrix

        W1 = zeros(dims(1),dims(1));
        for j_s = 1:s
            % extract partition similarity matrix
            K1_new=K1(j_s,:);
            KK1 = (K1_new' == K1_new);
            W1 = W1 + KK1;
        end% end of j_s creating W1 similarity matrix
        W1 = W1 / s;
        W1_w{w_ss}=W1;
    end % end of window loop

    %Final input for multilayer modularity after conducting multilayer
    %consensus

    A= W1_w';
  

    gamma =1;%
    omega =1;%
    [B1_new,mm_new] = multiord(A,gamma,omega);
    % multilayer modularity set-up
    PP_new = @(S) postprocess_ordinal_multilayer(S,num_window);
    % S = modularity assignments; Q1 = modularity value; mod_iter_tmp = iterations
    [S_new2,Q1_new,mod_iter_tmp_new] = iterated_genlouvain(B1_new,10000,0,1,'moverandw',[],PP_new);
    % 2D representation of modularity assignments
    S_new3 = reshape(S_new2,dims(1),num_window);
    % saving the node assignment by window
    node_assignment_by_window_sub{ts} =S_new3;
    % total number of community in network
    number_of_communities_new = sum(diff(sort(S_new3)) ~= 0) +1;
    % 14 node switching rate
    node_switching_fin_sub = flexibility(S_new3');
    % total number of community in network
    number_of_communities_sub{ts,:}= [number_of_communities_new];
    % 14 node switching rate
    node_switching_sub(ts,:) = node_switching_fin_sub;
   
end %end loop for ts
end
