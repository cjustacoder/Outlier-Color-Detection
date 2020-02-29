function [label, repre]= myKmean(X, itertime)
% can only do 2 means
[m,~] = size(X);
dist=zeros(m,m);
for i = 1:m
    dist(:,i) = fitinKLDiv(X(i,:),X);
end
repre = randi(m,1,2);
label = zeros(m,1);
for i = 1:itertime
    % fix mean, find label
    for j = 1:m
        if dist(repre(1),j) < dist(repre(2),j)
            label(j) = 0;
        else
            label(j) = 1;
        end
    end
    % fix label, find mean
    
end

end