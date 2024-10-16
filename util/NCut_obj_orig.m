%% Original NCut objective function
function [obj] = NCut_obj_orig(X,D,W)
% X: N * K
% D: N * N
% W: N * N

K = size(X,2);
results = zeros(1,K);
for i = 1:K
    results(i) = (X(:,i)'*(D-W)*X(:,i))./(X(:,i)'*D*X(:,i));
end
obj = sum(results)/2;

end