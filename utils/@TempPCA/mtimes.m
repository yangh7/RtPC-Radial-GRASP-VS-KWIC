function res = mtimes(a,b)
% (c) Ricardo Otazo 2008
if isa(a,'TempPCA') == 0
    error('In  A*B only A can be TempPCA operator');
end
if a.adjoint
      [nx,ny,nt]=size(b.data);  
%       b.data=b.data*sqrt(size(b.data,3));
      for ii=1:nt
          temp=b.data(:,:,ii);
          Data_Seq(:,ii)=reshape(temp,nx*ny,1);
      end
      Data_Seq=Data_Seq';
      Data_PCA=(b.PC*Data_Seq)';
      for ii=1:nt
          temp=Data_PCA(:,ii);
          Data_PCA_Back(:,:,ii)=reshape(temp,nx,ny);
      end
      res=Data_PCA_Back;
 else  
    [nx,ny,nt]=size(b);  
    for ii=1:nt
        temp=b(:,:,ii);
        Data_Seq(:,ii)=reshape(temp,nx*ny,1);
    end
    covariance=cov(Data_Seq);
    [PC, V] = eig(covariance);
    V = diag(V);
    [junk, rindices] = sort(-1*V);
    V = V(rindices);
    PC = PC(:,rindices);
    Data_PCA = (PC' * Data_Seq')';
    for ii=1:nt
        temp=Data_PCA(:,ii);
        Data_PCA_New(:,:,ii)=reshape(temp,nx,ny);
    end
    res.data=Data_PCA_New;
    res.PC=PC;
end