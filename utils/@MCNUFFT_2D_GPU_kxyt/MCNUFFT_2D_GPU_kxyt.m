function  res = MCNUFFT_2D_GPU_kxyt(kk,w,b1)
% % addpath C:\Users\2014GPU\Documents\MATLAB\gpuNUFFT-2.0.8\gpuNUFFT
% % addpath C:\Users\exx\Desktop\gpuNUFFT-master\gpuNUFFT
over_sample = 2;
Grid_kern = 3;
Region = 8;

Nd = size(b1(:,:,1));
Jd = [6,6];
% Jd = [6,6,6];
Kd = floor([Nd*1.5]);
n_shift = Nd/2;

for iii = 1:size(kk,3)
% %     fprintf( 'Time dimension: %d\n', iii );
    kkk = kk(:,:,iii);
    ww = w(:,:,iii).^2;
    res.st(:,:,iii) = gpuNUFFT([real(col(kkk)), imag(col(kkk))]',(col(ww)),over_sample,Grid_kern,Region,[Nd],[b1]);
end

res.adjoint = 0;
res.dataSize = size(kk);
res.imSize = size(b1);

res = class(res,'MCNUFFT_2D_GPU_kxyt');

