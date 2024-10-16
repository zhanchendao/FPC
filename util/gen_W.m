%% generate similarity matrix W
function [W,D] = gen_W(data)
    %data normalization
    for i = 1:size(data,2)
            if(norm(data(:,i)) == 0)
                continue;
            end
            data(:,i) = data(:,i)/norm(data(:,i));
    end
    %Gaussian Kernel
    d = squareform(pdist(data));
    W = exp(-d);
    W = (W+W')/2;
    D = diag(sum(W));
end

