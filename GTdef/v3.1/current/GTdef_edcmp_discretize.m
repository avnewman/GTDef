function [ pntsrc ] = GTdef_edcmp_discretize(M,edgrn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GTdef_edcmp_discretize.m					%
%											%
% Distretize one rectangular plane interface into point sources				%
%											%
% INPUT											%
% (1) M = [slip;north;east;depth;length;width;str;dip;rake] {row vectors}		%
%     in EDGRN: slip - dislocation; north - xs; east - ys; depth - zs			%
%     M represents only one plane source						%
%											%
% (2) edgrn structure									%
%---------------------------------------------------------------------------------------%
%											%
% OUTPUT										%
% pntsrc = [ pxs pys pzs pmoment ]							%
% pxs pys pzs - coordinates of point sources						%
% pmoment {5*pntsrc.num}								%
%	1 = weight for strike-slip: m12=m21=1;						%
%           poloidal*sin(2 * theta), toroidal*cos(2 * theta)				%
%	2 = weight for dip-slip: m13=m31=1                                              %
%           poloidal * cos(theta), toroidal * sin(theta)				%
%	3 = weight for clvd: m33=-m11=-m22=1                                            %
%           axisymmetric								%
%	4 = weight for 45 deg strike-slip: m11=-m22=1                                   %
%           greenfct4(theta) = green1(theta + 45 deg)					%
%	5 = weight for 45 deg dip-slip: m23=m32=1                                       %
%           greenfct5(theta) = green2(theta - 90 deg)					%
%											%
% REFERENCE  										%
% Wang, R., Martin, F. L., & Roth, F. (2003)						%
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		%
% based on edcmpf77_2.0/edcdisc.F							%
%		                                                                	%
% first created by lfeng Mon Feb 27 11:20:42 SGT 2012					%
% last modified by lfeng Tue Feb 28 20:02:16 SGT 2012					%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% discretise one rectangular plane sources %%%%%%%%%%%%%%%%%%%%%%%
% M = [slip;north;east;depth;length;width;str;dip;rake]
slip   = M(1,1);	% usually unit
xs     = M(2,1);	% north	
ys     = M(3,1);	% east
zs     = M(4,1);	% depth
length = M(5,1);
width  = M(6,1); 
str    = M(7,1);	% st
dip    = M(8,1);	% di
rake   = M(9,1);	% ra

dlength = edgrn.dr;
if dip>0
   dwidth = min([ edgrn.dr edgrn.dz/sind(dip) ]);
else
   dwidth=edgrn.dr;		% dip==0
end

% moment tensor
sm = zeros(3,3);
sm(1,1) = -sind(dip)*cosd(rake)*sind(2.0*str) - sind(2.0*dip)*sind(rake)*(sind(str))^2;
sm(2,2) =  sind(dip)*cosd(rake)*sind(2.0*str) - sind(2.0*dip)*sind(rake)*(cosd(str))^2;
sm(3,3) = -(sm(1,1)+sm(2,2));
sm(1,2) =  sind(dip)*cosd(rake)*cosd(2.0*str) + 0.5*sind(2.0*dip)*sind(rake)*sind(2.0*str);
sm(2,1) = sm(1,2);
sm(1,3) = -cosd(dip)*cosd(rake)*cosd(str) - cosd(2.0*dip)*sind(rake)*sind(str);
sm(3,1) = sm(1,3);
sm(2,3) = -cosd(dip)*cosd(rake)*sind(str) + cosd(2.0*dip)*sind(rake)*cosd(str);
sm(3,2) = sm(2,3);

nx = max([ 1 round(length/dlength) ]);
ny = max([ 1 round(width/dwidth)   ]);
dx = length*1.0/nx;
dy = width*1.0/ny;
pnum = nx*ny;

% IMPORTANT: slip needs to be multiply by area!
if dx>0.0
   slip = slip*dx;
end
if dy>0.0
   slip = slip*dy;
end

% create discretized point sources
xlin = dx*([1:1:nx]-0.5);	% columnwise
ylin = dy*([1:1:ny]-0.5);	% rowwise
[ xxmat,yymat ] = meshgrid(xlin,ylin);
xx  = reshape(xxmat,[],1); yy = reshape(yymat,[],1);
pxs = xs + xx.*cosd(str) - yy.*cosd(dip)*sind(str);
pys = ys + xx.*sind(str) + yy.*cosd(dip)*cosd(str);
pzs = zs + yy.*sind(dip);

ind_z = pzs<edgrn.minz-edgrn.dz;
if any(ind_z)
   error('GTdef_edcmp_discretize ERROR: parts of rectangular source are shallower than green function grids!');
end
ind_z = pzs>edgrn.maxz+edgrn.dz;
if any(ind_z)
   error('GTdef_edcmp_discretize ERROR: parts of rectangular source are deeper than green function grids!');
end

% assign the same pmoment weights to all point sources
pmoment = zeros(pnum,5);
pmoment(:,1) = sm(1,2)*slip;
pmoment(:,2) = sm(1,3)*slip;
pmoment(:,3) = sm(3,3)*slip;
pmoment(:,4) = 0.5*(sm(1,1)-sm(2,2))*slip;
pmoment(:,5) = sm(2,3)*slip;

pntsrc = [ pxs pys pzs pmoment ];
