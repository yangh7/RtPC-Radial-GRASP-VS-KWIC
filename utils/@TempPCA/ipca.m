function [datanew] = ipca(data,PC_Transform,Mean_Transform)

data=(PC*signals)+mean_sub;

[M,N,S] = size(data);
for i=1:M
data1=squeeze(data(i,:,:))*sqrt(size(data,3));
data1=data1';
PC=PC_Transform(:,:,i);
mean_sub=Mean_Transform(:,i);
signal=(PC*data1)+mean_sub;
datanew(i,:,:)=signal';
clear signal data1
end