%% demo
%   Implementation of paper: Multidimensional Fractional Programming for Normalized Cuts
%   SPDX-FileCopyrightText: 2024 Beichen Huang <polarishuang0.0@gmail.com>
%   SPDX-License-Identifier: Apache-2.0
clear 
clc

dataset_name = {'Breast','Thyroid','Office+Caltech10','Splice','Rice','Landsat','USPS','Epileptic'};
data_num = size(dataset_name,2);
Max_round = 10;
for data_index = 1:data_num
    load(char(dataset_name(data_index)));
    K = length(unique(label)); % number of clusters
    N = size(data,1);          % number of instances  
    fprintf("\n##Dataset: %s##\n",char(dataset_name(data_index)));
    [W,D] = gen_W(data);       %generate similarity matrix W
    
    % mineig = min(eig(W));
    % if mineig<0
    %     lambda = abs(mineig)/min(diag(D));
    % else
    %     lambda = 1e-8;
    % end

    lambda = 1e-8;
    F = W+lambda*D;

    all_obj = zeros(1,Max_round);
    all_clu = zeros(7,Max_round);
    all_iter = zeros(1,Max_round);
    all_time = zeros(1,Max_round);
    for round = 1:Max_round
        tic;

        % generate random starting point
        X0 = zeros(N,K);
        for i = 1:N
            X0(i,randi(K)) = 1;
        end
        
        % do FPC
        fprintf("=== FPC algo Round %d begin ===\n",round);
        [X_out] = FPC_algo(D,W,K,X0,N,F);

        pre =  zeros(N, 1);
        for i = 1:N
            [~, idx] = max(X_out(i, :));
            pre(i) = idx;
        end
        all_obj(1,round) = NCut_obj_orig(X_out,D,W); %calculate the orignal NCut obj
        all_clu(:,round) = ClusteringMeasure_All(label, pre);
        all_time(1,round) = toc;
    end
    fprintf("All rounds done\n")
    avgT = mean(all_time);
    [~,smallest_obj_index] = min(all_obj);
    best_obj = all_obj(smallest_obj_index);
    best_clu = all_clu(:,smallest_obj_index);
    fprintf("obj: %.6f\nACC: %.4f\nNMI: %.4f\nARI: %.4f\naverage T: %.4fs\n",best_obj,best_clu(1),best_clu(2),best_clu(7),avgT);
end








