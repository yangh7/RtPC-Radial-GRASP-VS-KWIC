function traj_delay_Coefficients = TrajectoryCorrection_RING( kSpace, param )
% RING method works in Linux system
%
% kSpace: input kSpace data with (nx,nray,ncoil), 
%         nray should be time-combined
% param.Indexes: kSpace ray number that to use for gradient delay estimation,
%                better pick rays that are at steady state: e.g., 1:Nray 
% param.N : tiny golden angle N (i.e.5).
% param.AngleRange :  1) '-H' for [0,180)
%                     2) '-G' for [0,360)
addpath(strcat(getenv('TOOLBOX_PATH'), '/matlab')); %for BART
kSpace = squeeze(kSpace);
[nx,np,nc] = size(kSpace);

indexes     = param.Indexes;
N           = param.N;
AngleRange  = param.AngleRange;

% convert to ring traj
% % % traj = bart(strcat('traj -c -r -H -x',num2str(nx),' -y',num2str(np),' -s',num2str(N)));% put -H if [0-180), -G [0-360)
Instruction = sprintf( 'traj -c -r %s -x%d  -y%d  -s%d', AngleRange, nx, np, N ); 
traj = bart( Instruction );% put -H if [0-180), -G [0-360)
k1 = reshape(kSpace,[1 nx np nc]);
tR = traj(:,:,indexes);
k1R = k1(:,:,indexes,:);


% estimate gradient delay
dfile = 'ring_estdelay.txt';
if exist(dfile, 'file'); delete(dfile);end
diary(dfile)
bart('estdelay -R ',tR,k1R);
diary off
fid = fopen(dfile);
traj_delay = fgetl(fid);
tmp = strfind( traj_delay, ':' );
if length(tmp)~=2 % Warning exists, read the next line
    traj_delay = fgetl(fid);
    tmp = strfind( traj_delay, ':' );
    traj_delay_Coefficients = [str2num(traj_delay(5:tmp(1)-1))  str2num(traj_delay(tmp(1)+1:tmp(2)-1))   str2num(traj_delay(tmp(2)+1:end))]*(-1);
else
    traj_delay_Coefficients = [str2num(traj_delay(1:tmp(1)-1))  str2num(traj_delay(tmp(1)+1:tmp(2)-1))   str2num(traj_delay(tmp(2)+1:end))]*(-1);
end
fclose(fid);

% tmp = strfind( traj_delay, ':' );
% traj_delay_Coefficients = [str2num(traj_delay(1:tmp(1)-1))  str2num(traj_delay(tmp(1)+1:tmp(2)-1))   str2num(traj_delay(tmp(2)+1:end))]*(-1);

end
