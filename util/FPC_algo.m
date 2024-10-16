%% FPC algo
%   Implementation of paper: Multidimensional Fractional Programming for Normalized Cuts
%   SPDX-FileCopyrightText: 2024 Beichen Huang <polarishuang0.0@gmail.com>
%   SPDX-License-Identifier: Apache-2.0
function [X] = FPC_algo(D,W,K,X0,N,F)
d = diag(D);
iter_max = 1e4;
iter_obj = zeros(iter_max+1,1);
X = X0;
valid_flag = 1;
for iter = 1:iter_max
    iter_obj(iter) = NCut_obj(X,D,W);
    
    if iter>2 && (iter_obj(iter)-iter_obj(iter-1))/iter_obj(iter-1)<1e-6
        fprintf("=> Converge at iter %d\n",iter);  
        break
    end
    if(isnan(iter_obj(iter)))
        valid_flag = 0;
        fprintf("!!! NaN at iter %d !!!\n",iter); 
        break;
    end
    % update H
    dx = repmat(d'*X,N,1);
    H = (2*F*X)./dx - repmat(diag(X'*F*X)',N,1) .* repmat(d,1,K)./ (dx.^2);
    % update X for the next iter
    [~,X_index] = max(H,[],2);    
    X = zeros(N,K);
    for i = 1:N
        X(i,X_index(i)) = 1;
    end
end
if iter == iter_max
    fprintf("reach the max iteration");
end
  
end



