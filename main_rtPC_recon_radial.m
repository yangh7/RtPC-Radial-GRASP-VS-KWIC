addpath(genpath('utils'));
vs_kwic = 3; %1-cs alone; 2-cs+vs; 3-cs+vs+kwic
nrays = 3; %native spokes
vs_layers = ones(5,1); %VS spokes 
kwic_wid = 3:15:64; %KWIC shape
center_ray = 31; %skip first 30 spokes (not necessary)

filename = '.dat'; %your filename

[recon,zf,NUFFT_CoilSens,b1sens,dummy] = RTPC_Recon_CS_VSKWIC(filename,nrays,vs_kwic,vs_layers,kwic_wid,center_ray);
