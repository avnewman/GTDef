function [ edgrn,edgrnfcts ] = GTdef_read_edgrn_output(edgrn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GTdef_read_edgrn_output.m					%
%											%
% Read green function library output by Fortran code EDGRN in folder 'edgrnfcts' locally%
%											%
% INPUT											%
% Point source library should cover the whole model region 				%
% 3 files 'izmhs.ss'  'izmhs.ds'  'izmhs.cl' in folder './edgrnfcts/'  			%
% ss - strike-slip; ds - dip-slip; cl - compensated linear vector dipole (CLVD)		%
%---------------------------------------------------------------------------------------%
%       GREEN'S FUNNCTIONN PARAMETERS							%
%       =============================                                                   %
%       Green's function source types:                                                  %
%         1 = strike-slip (m12=m21=1)                                                   %
%         2 = dip-slip (m13=m31=1)                                                      %
%         3 = compensated linear vector dipole (CLVD)                                   %
%             (m11=m22=-1/2, m33=1) (no tangential component)                           %
%       Green's function coordinate system:						%
%         (z,r,t) = cylindrical with z being downward(!)				%
%---------------------------------------------------------------------------------------%
%											%
% OUTPUT										%
% (1) edgrn structure                                                                   %
%											%
% (2) edgrnfcts structure								%
% 											%
% REFERENCE  										%
% Wang, R., Martin, F. L., & Roth, F. (2003)						%
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		%
%		                                                                	%
% first created by lfeng Mon Feb 27 01:57:34 SGT 2012					%
% last modified by lfeng Wed Feb 29 02:03:22 SGT 2012					%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderName = './edgrnfcts/';
fssName = [ folderName 'izmhs.ss' ];
fdsName = [ folderName 'izmhs.ds' ];
fclName = [ folderName 'izmhs.cl' ];

if ~exist(folderName,'dir') || ~exist(fssName,'file') || ~exist(fdsName,'file') || ~exist(fclName,'file')
    error('GTdef_read_edgrn_output ERROR: need to calculate discretized point sources using EDGRN first!');
end

% read in edgrn structure locally
obsz = edgrn.obsz;
nr   = edgrn.nr;
minr = edgrn.minr;
maxr = edgrn.maxr;
nz   = edgrn.nz;
minz = edgrn.minz;
maxz = edgrn.maxz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nr->nr+2 nz->nz+2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extend another two records
ssdisp0  = zeros(3,nr+2,nz+2);	% to avoid conflict naming with Matlab
ssstrn0  = zeros(6,nr+2,nz+2);
ssuzr0   = zeros(nr+2,nz+2);
dsdisp0  = zeros(3,nr+2,nz+2);
dsstrn0  = zeros(6,nr+2,nz+2);
dsuzr0   = zeros(nr+2,nz+2);
cldisp0  = zeros(2,nr+2,nz+2);
clstrn0  = zeros(4,nr+2,nz+2);
cluzr0   = zeros(nr+2,nz+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in strike-slip file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fss = fopen(fssName,'r');
param_cell = textscan(fss,'%f',9,'CommentStyle','#');
param = cell2mat(param_cell);
% check if ss file is consistent with edgrn structure
if param(1)~=nr || param(2)~=minr || param(3)~=maxr || ...
   param(4)~=nz || param(5)~=minz || param(6)~=maxz || ...
   param(7)~=obsz
      error('GTdef_read_edgrn_output ERROR: %s is not consistent with GTdef input file!',fssName);
end
% 10 data values in each line
data_cell  = textscan(fss,'%f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
data_trans = cell2mat(data_cell)';
ssdisp0(:,1:end-2,1:end-2) = reshape(data_trans(1:3,:),3,nr,nz);
ssstrn0(:,1:end-2,1:end-2) = reshape(data_trans(4:9,:),6,nr,nz);
ssuzr0(1:end-2,1:end-2)    = reshape(data_trans(10,:),nr,nz);
fclose(fss);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in dip-slip file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fds = fopen(fdsName,'r');
param_cell = textscan(fds,'%f',9,'CommentStyle','#');
param = cell2mat(param_cell);
% check if ss file is consistent with edgrn structure
if param(1)~=nr || param(2)~=minr || param(3)~=maxr || ...
   param(4)~=nz || param(5)~=minz || param(6)~=maxz || ...
   param(7)~=obsz
      error('GTdef_read_edgrn_output ERROR: %s is not consistent with GTdef input file!',fdsName);
end
% 10 data values in each line
data_cell  = textscan(fds,'%f %f %f %f %f %f %f %f %f %f','CommentStyle','#');
data_trans = cell2mat(data_cell)';
dsdisp0(:,1:end-2,1:end-2) = reshape(data_trans(1:3,:),3,nr,nz);
dsstrn0(:,1:end-2,1:end-2) = reshape(data_trans(4:9,:),6,nr,nz);
dsuzr0(1:end-2,1:end-2)    = reshape(data_trans(10,:),nr,nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in clvd file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fcl = fopen(fclName,'r');
param_cell = textscan(fcl,'%f',9,'CommentStyle','#');
param = cell2mat(param_cell);
% check if ss file is consistent with edgrn structure
if param(1)~=nr || param(2)~=minr || param(3)~=maxr || ...
   param(4)~=nz || param(5)~=minz || param(6)~=maxz || ...
   param(7)~=obsz
      error('GTdef_read_edgrn_output ERROR: %s is not consistent with GTdef input file!',fclName);
end
% 7 data values in each line
data_cell  = textscan(fcl,'%f %f %f %f %f %f %f','CommentStyle','#');
data_trans = cell2mat(data_cell)';
cldisp0(:,1:end-2,1:end-2) = reshape(data_trans(1:2,:),2,nr,nz);
clstrn0(:,1:end-2,1:end-2) = reshape(data_trans(3:6,:),4,nr,nz);
cluzr0(1:end-2,1:end-2)    = reshape(data_trans(7,:),nr,nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if observation point coincides with source point, output is zero
% so use the nearest to replace the zero output
dr=(maxr-minr)/(nr-1);
dz=(maxz-minz)/(nz-1);
if  obsz>=minz && obsz<=maxz && minr<1.0e-6
   iz  = round((obsz-minz)/dz);
   dzs = abs(obsz-(minz+dz*iz))/dz;
   if dzs<0.02
      iz = iz+1;
      ssdisp0(:,1,iz) = ssdisp0(:,2,iz);
      dsdisp0(:,1,iz) = dsdisp0(:,2,iz);
      cldisp0(:,1,iz) = cldisp0(:,2,iz);
      ssstrn0(:,1,iz) = ssstrn0(:,2,iz);
      dsstrn0(:,1,iz) = dsstrn0(:,2,iz);
      clstrn0(:,1,iz) = clstrn0(:,2,iz);
      ssuzr0(1,iz)    = ssuzr0(2,iz);
      dsuzr0(1,iz)    = dsuzr0(2,iz);
      cluzr0(1,iz)    = cluzr0(2,iz);
   end
end

% assign to structure
edgrnfcts.ssdisp = ssdisp0;
edgrnfcts.dsdisp = dsdisp0;
edgrnfcts.cldisp = cldisp0;
edgrnfcts.ssstrn = ssstrn0;
edgrnfcts.dsstrn = dsstrn0;
edgrnfcts.clstrn = clstrn0;
edgrnfcts.ssuzr  = ssuzr0;
edgrnfcts.dsuzr  = dsuzr0;
edgrnfcts.cluzr  = cluzr0;
edgrn.dr = dr;
edgrn.dz = dz;
