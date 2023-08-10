
close all; clear all; clc

NM = 192
Nray = 10
Nrepeat = 3
sel_WC = 0
sel_angle = 1

switch sel_angle
    case 1   % golden angle
        GA = pi/((sqrt(5)+1)/2);
        AngleInfo = mod( (0:Nray*Nrepeat-1)*GA, 2*pi )';
    case 2   % tiny golden angle
        Ntiny = 5;
        TGA = pi/((sqrt(5)+1)/2 + (Ntiny-1) )
        AngleInfo = mod( (0:Nray*Nrepeat-1)*TGA, 2*pi )';
    case 3   % regular angle
        RA = pi/(Nray*Nrepeat)
        AngleInfo = ((0:Nray*Nrepeat-1)*RA)';
end
AngleInfo = reshape( AngleInfo, [Nray Nrepeat] );

Option.Nsamples = NM;
Option.AnglePrecision = 0;   % full precision
Option.Display = 1;
if sel_WC == 0
    Option.WeightedContrast = 0
else
    LocOfMainContrastRay = [3 5 7];   % range from 3 to Nray-2 in this DEMO setting
    Option.WeightedContrast = ones( NM, Nray, Nrepeat );
    Option.WeightedContrast( NM/2+1-2:NM/2+1+2, :, : ) = 0;
    for ii = 1:Nrepeat
        Option.WeightedContrast(  NM/2+1,              LocOfMainContrastRay(ii),                              ii ) = 1;
        Option.WeightedContrast( [NM/2+1-1  NM/2+1+1], LocOfMainContrastRay(ii)-1:LocOfMainContrastRay(ii)+1, ii ) = 1;
        Option.WeightedContrast( [NM/2+1-2  NM/2+1+2], LocOfMainContrastRay(ii)-2:LocOfMainContrastRay(ii)+2, ii ) = 1;
    end
end

DCF = AdvancedDCF( AngleInfo, Option );
figure; imagesc( DCF(:,:,1), [0 1] ); colormap jet
