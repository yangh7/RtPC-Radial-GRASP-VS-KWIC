function  res = MCNUFFT_2D_GPU_kxyt_single(k,w,b1)
% addpath C:\Users\2014GPU\Documents\MATLAB\gpuNUFFT-2.0.8\gpuNUFFT
over_sample = 2;
Grid_kern = 5;
Region = 8;

% Nd = [size(b1,1),size(b1,1),size(zz,2)];
Nd = size(b1(:,:,1));
Jd = [6,6];
% Jd = [6,6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;

ww = w.^2;
res.st = gpuNUFFT([real(col(k)), imag(col(k))]',(col(ww)),over_sample,Grid_kern,Region,[Nd],[]);

size(res.st)
res.adjoint = 0;
res.imSize = size(b1(:,:,1));
res.dataSize = size(k);
res.b1 = b1;

res = class(res,'MCNUFFT_2D_GPU_kxyt_single');

