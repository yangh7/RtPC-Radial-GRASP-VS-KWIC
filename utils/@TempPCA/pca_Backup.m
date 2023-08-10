function [result,PC,mean_sub] = pca(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2D images
% % % % [M,N] = size(data);
% % % % % subtract off the mean for each dimension
% % % % mn = mean(data,2);
% % % % mean_sub=repmat(mn,1,N);
% % % % data = data - repmat(mn,1,N);
% % % % % calculate the covariance matrix
% % % % covariance=cov(data');
% % % % % find the eigenvectors and eigenvalues
% % % % [PC, V] = eig(covariance);
% % % % % extract diagonal of matrix as vector
% % % % V = diag(V);
% % % % % sort the variances in decreasing order
% % % % [junk, rindices] = sort(-1*V);
% % % % V = V(rindices);
% % % % PC = PC(:,rindices);
% % % % % project the original data set
% % % % signals = PC' * data;

[M,N,S] = size(data);
for i=1:M
data1=squeeze(data(i,:,:));
data1=data1';
mn = mean(data1,2);
mean_sub=repmat(mn,1,N);
data1 = data1 - mean_sub;
covariance=cov(data1');
[PC, V] = eig(covariance);
V = diag(V);
[junk, rindices] = sort(-1*V);
V = V(rindices);
PC = PC(:,rindices);
signal = PC' * data1;
datanew(i,:,:)=signal'/sqrt(size(data,3));
clear signal data1 covariance PC
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3D images
[M,N,S] = size(data);
temp=zeros(N,M*S);
for i=1:M
    temp(:,(i-1)*S+1:i*S)=squeeze(data(i,:,:));
end
data=temp';clear temp
mn = mean(data,2);
mean_sub=repmat(mn,1,M);
data = data - repmat(mn,1,M);
covariance=cov(data');
[PC, V] = eig(covariance);
V = diag(V);
[junk, rindices] = sort(-1*V);
V = V(rindices);
PC = PC(:,rindices);
signals = PC' * data;
signals=signals';
result=zeros(M,N,S);
for i=1:M
    result(i,:,:)=signals(:,(i-1)*S+1:i*S);
end