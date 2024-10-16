%% Maximize version of NCut objective function
function [obj] = NCut_obj(X,D,W)
% X: N * K
% D: N * N
% W: N * N

K = size(X,2);
results = zeros(1,K);
for i = 1:K
    results(i) = (X(:,i)'*(W)*X(:,i))./(diag(D)'*X(:,i));
end
obj = sum(results);

end