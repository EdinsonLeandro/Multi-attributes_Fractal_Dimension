function X=normalization(V,W)
% X = normalization (V, W) performs data normalization, according to
% function proposed by Rafael Banchs in 2001.
%
% Input parameters:
% V: Array with the data to be normalized.
% W: Array with the data in the area of ??interest.
%
% Output Parameters:
% X: Array with each normalized data

% dev: Data Deviation. Difference with its average in the area of ??interest.
desv=zeros(size(V));
for i=1:size(V,2)
    desv(:,i)=V(:,i)-mean(W(:,i));
end

maxdesv=max(abs(desv));

% If in some data the samples are the same, its deviation would be zero.
% This generates error in normalization function. Therefore, it is assumed
% to be constant and equal to 1.

maxdesv(maxdesv == 0) = 1;

% X: Each column has data normalized.

X=zeros(size(V));
for i=1:size(V,2)
    X(:,i)=((V(:,i)-mean(W(:,i)))/maxdesv(i))+1;
end
