function [recon, zf, NUFFT_CoilSens, b1sens, dummy] = RTPC_Recon_CS_VSKWIC(filename, nrays, vs_kwic, vs_layers, kwic_wid, center_ray)
% Using Siemens vb13 program to extract kpace data as an m file
% RING trajectory correction
% ============Inputs============
% filename
% nrays: native spokes
% vs_kwic: 1 = no view sharing, no kwic
%          2 = view sharing, no kwic
%          3 = view sharing, with kwic
% vs_layers: shared spokes (symmetrical)
% kwic_wid: KWIC shape
% center_ray: middle spoke of the 1st frame
% ============Outputs============
% recon: reconstructed images
% zf: zero-filled images
% NUFFT_CoilSens: single-coil time-average images
% b1sens: coil sensitivity maps
% dummy: coil-combined time-average images

N_vs = sum(vs_layers);
skip_ray = center_ray-ceil(nrays/2); %to match middle spoke of all methods

% select lambda for iterative CS recon
all_lambda = [0.0100 0.0075 0.0050 0.0025 0.0010 0.0005];
sel_lambda = 5; %4 for CS, 5 for KWIC

% # of rays used for reconstruction per frame
sel_NR = nrays
if vs_kwic>1
    sel_NR = nrays + 2*N_vs + (mod(nrays,2)-1)
end

% set the length of phase correction data
% fixed number
sel_PC_size = 20

%% Read data
dataraid = ReadSiemensMeasVD13_idea( filename ); %[path_data filename] );
tmp_ind = size(dataraid,2);
rawdata2 = dataraid(tmp_ind).rawdata;
mdh = dataraid(tmp_ind).loopcounters;
nMeasurements = length( dataraid(tmp_ind).loopcounters );
clear dataraid
clear lin cha slc phs set angle_mdh EchoSpacing  info3 info4
for iM=1:nMeasurements
    cha(iM)        = mdh(iM,15);    % coil
    lin(iM)        = mdh(iM,1);     % rays in each coil
    % %         rep(iM)        = mdh(iM,7);     % # of images
    slc(iM)        = mdh(iM,3);     % # of images
    phs(iM)        = mdh(iM,6);     % # of contrasts
    set(iM)        = mdh(iM,8);
    angle_mdh(iM)  = mdh(iM,43); %39);
    EchoSpacing(iM) = mdh(iM,40);   % TE
    info3(iM) = mdh(iM,41);   % TE
    info4(iM) = mdh(iM,42);   % TE
end
angle_mdh(1) = angle_mdh(2);
angle_mdh = angle_mdh./10000;
EchoSpacing = EchoSpacing(2)/1000;

clear kSpace angle_info
for iM=1:nMeasurements
    kSpace(:,lin(iM),phs(iM),cha(iM),set(iM)) = (rawdata2(iM,:));
    angle_info( lin(iM),phs(iM),set(iM) ) = angle_mdh(iM);
end
clear rawdata2  cha lin slc  rep  angle_mdh phs set
n1 = min(find(kSpace(1,:,1,1,1)));
n2 = max(find(kSpace(1,:,1,1,1)));
kSpace = kSpace(:,n1:n2,:,:,:);

size(kSpace)
clear tmp_kSpace  Angles
kSpace_PC = kSpace(:,:,1:sel_PC_size,:,:);
tmp_kSpace = kSpace(:,:,sel_PC_size+1:end,:,:);
PC_angle_inc = pi/(sel_PC_size/2*size(kSpace,2));
tmp_PC_angle = (0:1:size(kSpace,2)*sel_PC_size/2-1)*PC_angle_inc;
PC_angles = zeros( 1, sel_PC_size*size(kSpace,2) );
PC_angles(1:2:end) = tmp_PC_angle;
PC_angles(2:2:end) = tmp_PC_angle+pi;
clear kSpace
kSpace = tmp_kSpace;
clear tmp_kSpace

[NM,NR,NP,NC,NS] = size( kSpace )
% Seg=360 in total
tmp_kSpace = reshape( kSpace, [NM, NR*NP, NC, NS] );
clear kSpace

%% View sharing
if vs_kwic==1
    tmp_kSpace = tmp_kSpace(:,1+skip_ray:floor((NR*NP-N_vs-skip_ray)/sel_NR)*sel_NR+skip_ray,:,:);
    kSpace = reshape( tmp_kSpace, [NM, sel_NR, floor((NR*NP-N_vs-skip_ray)/sel_NR), NC, NS] );
else %VS
    kSpace = zeros(NM,sel_NR,floor((NR*NP-N_vs-skip_ray)/nrays),NC,NS);
    for tt=1:size(kSpace,3)
        kSpace(:,:,tt,:,:) = tmp_kSpace(:,(tt-1)*nrays+1+skip_ray-N_vs:tt*nrays+N_vs+skip_ray +(mod(nrays,2)-1),:,:);
    end
end
clear tmp_kSpace

N = 1;
tau = (1.0+sqrt(5.0))/2.0;
golden_angle = pi / (tau+N-1);
if vs_kwic==1
    Angles = mod( golden_angle*(skip_ray:floor((NR*NP-N_vs-skip_ray)/sel_NR)*sel_NR+skip_ray-1), 2*pi );
else %VS
    Angles = zeros(1,sel_NR*(floor((NR*NP-N_vs-skip_ray)/nrays)));
    for tt=1:size(kSpace,3)
        Angles((tt-1)*sel_NR+1:tt*sel_NR) = mod(golden_angle*((tt-1)*nrays+skip_ray-N_vs:tt*nrays+N_vs-1+skip_ray +(mod(nrays,2)-1)),2*pi);
    end
end
Option.Angles = Angles;
[NM,NR,NP,NC,NS] = size( kSpace )

%% Phase correction (RING)
Option.Indexes = 1:NR*NP;
Option.N = N;
Option.AngleRange = '-G'
disp( 'Trajectory correction (RING)...' )
Option.PC_coeff = TrajectoryCorrection_RING( reshape(kSpace(:,:,:,:,1),[NM,NR*NP,NC]), Option )
Option.PhaseCorrection = 'UseTC';
for ReconType = 1:2
    PCcorrected = RadialPhaseCorrection( kSpace(:,:,:,:,ReconType), kSpace(:,:,:,:,ReconType), Option );
    kSpace(:,:,:,:,ReconType) = PCcorrected.kSpace;
end
NumOfComp = NC;

%% DCF + NUFFT (time-average) + Coil sensitivity
clear tmp_kSpace  tmp_angleinfo  ray_index  recon_nufft  ku wu
tmp_angleinfo = Option.Angles(:);

fprintf( 'Calculate coil sensitivity....\n' );
ray_info = -0.5:1/NM:0.5-1/NM;
ray_grid = [];
for aaa = 1:NP*NR
    ray_grid(:,aaa) = ray_info.*exp(1i*tmp_angleinfo(aaa));
end
ray_grid = -imag(ray_grid) - 1i*real(ray_grid);

Option.Nsamples = NM;
Option.AnglePrecision = 0;
Option.Display = 0;
Option.WeightedContrast = 0;
DCF_b1 = AdvancedDCF_2DRadial( tmp_angleinfo, Option );

clear b1sens
NUFFT_b1 = MCNUFFT_2D_GPU_kxyt_single(ray_grid,DCF_b1,ones(NM,NM,NumOfComp) );
for recon_type = 1:NS
    kspace_b1 = reshape( kSpace(:,:,:,:,recon_type), [NM, NR*NP, NumOfComp] );
    NUFFT_CoilSens = NUFFT_b1'*kspace_b1;
    [dummy,tmp_b1sens] = adapt_array_2d( squeeze(NUFFT_CoilSens) );
    b1sens(:,:,:,recon_type) = tmp_b1sens./max( abs(tmp_b1sens(:)) );
    MaxIntRef_inB1sens = max(abs(dummy(:)));
    clear tmp_b1sens %dummy  NUFFT_CoilSens
end
clear  param  NUFFT_CoilSens_Type1
b1sens(:,:,:,2) = abs(b1sens(:,:,:,2)).*exp( 1i.*angle(b1sens(:,:,:,1)) ); %synthetic coil sens
tmp_combined = [ b1sens(:,:,:,1)   b1sens(:,:,:,2) ];

clear  ray_info  ray_grid
ray_info = -0.5:1/NM:0.5-1/NM;
NR2 = sel_NR
NP2 = size(kSpace,3)
clear kSpace_reorgan   tmp_data
for recon_type = 1:NS
    tmp_data = reshape( kSpace(:,:,:,:,recon_type), [NM  NR2*NP2  NumOfComp] );
    kSpace_reorgan(:,:,:,:,recon_type) = reshape( tmp_data(:,1:NR2*NP2,:), [NM  NR2  NP2  NumOfComp] );
end
clear tmp_data
tmp_ang = Option.Angles(1:NR2*NP2);
Option.Angles_reogran = reshape( tmp_ang, [NR2,  NP2] );
for aaa = 1:NP2
    for bbb = 1:NR2
        ray_grid(:,bbb,aaa) = ray_info.*exp(1i*Option.Angles_reogran(bbb,aaa));
    end
end
ray_grid = -imag(ray_grid) - 1i*real(ray_grid);

%% KWIC filter
DCF_CS = ones(NM,NR2,NP2);
if vs_kwic == 3 %with kwic
    pos = 0;
    for ll=length(vs_layers):-1:1
        tmp = floor(kwic_wid(ll)/2);
        DCF_CS(NM/2-tmp:NM/2+tmp, [pos+1:pos+vs_layers(ll),NR2-pos-vs_layers(ll)+1-(mod(nrays,2)-1):min(NR2,NR2-pos-(mod(nrays,2)-1))],:)=0;
        pos=pos+vs_layers(ll);
    end
end

%% NUFFT + CS recon (TTV)
param.TPCA = TempPCA();
param.TTV = TV_Temp();
param.STV = TV_2DPERF();
param.Wavelet = Wavelet();
param.nite = 8;
param.display = 1;
param.VW = 'on';

lambda = all_lambda(sel_lambda);

clear  Recon_CS_TTV_TPCA   Recon_RadRTPC_grid  Recon_iterCS_Times Recon_RadRTPC_grid_all
for ReconType = 1:NS
    param.E = MCNUFFT_2D_GPU_kxyt( ray_grid, DCF_CS, b1sens(:,:,:,ReconType) );
    param.y = permute( kSpace_reorgan(:,:,:,:,ReconType), [1,2,4,3] );
    Recon_RadRTPC_grid = param.E'*param.y;
    Recon_RadRTPC_grid_all(:,:,:,ReconType) = Recon_RadRTPC_grid;
    
    if strcmp( param.VW, 'on' )
        param.STVw = 0;
        param.TTVw = lambda;
        param.Waveletw = 0;
        param.TPCAw = 0;
    else
        param.STVw = 0;
        param.TTVw = MaxIntRef_inB1sens*lambda;
        param.Waveletw = 0;
        param.TPCAw = 0;
    end
    param
    
    tmp_Recon_CS_TTV_TPCA = Recon_RadRTPC_grid;
    ITER = 3
    for ite = 1:ITER,
        disp( [ sprintf( 'iteration = %2d', ite )] )
        tmp_Recon_CS_TTV_TPCA = CSL1NlCg_TTV_TPCA_STV_VW( tmp_Recon_CS_TTV_TPCA, param );
    end
    Recon_CS_TTV_TPCA(:,:,:,ReconType) = tmp_Recon_CS_TTV_TPCA;
end

zf = Recon_RadRTPC_grid_all;
recon = Recon_CS_TTV_TPCA;

end