% function [datanew,PC_Transform,Mean_Transform] = pca(data)

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
PC_Transform(:,:,i)=PC;
Mean_Transform(:,i)=mean_sub;
signal = PC' * data1;
datanew(i,:,:)=signal'/sqrt(size(data,3));
clear data1 covariance PC mean_sub signal
end