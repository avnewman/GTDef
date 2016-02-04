function [ disp,strain,stress,tilt ] = GTdef_edcmp(edgrn,layer,edgrnfcts,pntsrc,Xin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           	     	GTdef_edcmp.m					%
%											%
% Add up contributions from all point sources						%
%											%
% INPUT											%
% (1) edgrn structure									%
%		edgrn.nl        	(scalar)                                	%
% 		edgrn.obsz     		(scalar)                                	%
% 		edgrn.nr	        (scalar)                                	%
% 		edgrn.minr,edgrn.maxr   (scalar)                                	%
% 		edgrn.nz                (scalar)                                	%
% 		edgrn.minz,edgrn.maxz   (scalar)                                	%
%		edgrn.srate		(scalar)					%	
%   		edgrn.dr		(scalar)					%
%   		edgrn.dz		(scalar)       					%
%---------------------------------------------------------------------------------------%
% (2) layer = [ id depth vp vs ro ]	(nn*5)						%
%---------------------------------------------------------------------------------------%
% (3) edgrnfcts structure								%
% point strike-slip source								%
%       ssdisp0(1-3): Uz, Ur, Ut                                                        %
%       ssstrn(1-6): Ezz,Err,Ett,Ezr=Erz,Ert=Etr,Etz=Ezt                                %
%	ssuzr(1)									%
% point dip-slip source                                                                 %
%       dsdisp(1-3): Uz, Ur, Ut                                                         %
%       dsstrn(1-6): Ezz,Err,Ett,Ezr=Erz,Ert=Etr,Etz=Ezt                                %
%	dsuzr(1)									%
% point clvd source                                                                     %
%       cldisp(1-2): Uz, Ur (Ut=0)                                                      %
%       clstrn(1-4): Ezz,Err,Ett,Ezr=Erz (Ert=Etr=Etz=Ezt=0)                            %
%	cluzr(1)									%
% Note ssdisp is a Matlab function, so ssdisp0 is used instead				%
%---------------------------------------------------------------------------------------%
% (4) discretized point source matrix							%
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
% Note weights here are azimuthal factors, because Green's functions are stored without %
% considering azimuth angles                                                            %       
%---------------------------------------------------------------------------------------%
% (5) observation locations								%
%  Xin - point site locations in the local cartesian system 	  	  		%
%        [3*n] [ xx;yy;zz ]						  		%
% Note: need to convert to EDCMP coordinate						%
% GTdef: X = east,  Y = north, Z = downward						%
% EDCMP: X = north, Y = east,  Z = downward	                                        %
%											%
% OUTPUT (rowwise to be consistent with Okada)                                          %
% disp   = [ Ux;Uy;Uz ] (3 row vectors)							%
% strain = [exx;eyy;ezz;eyz;exz;exy] (6 row vectors)                                    %
% stress = [sxx;syy;szz;syz;sxz;sxy] (6 row vectors)                                    %
% tilt (2 row vectors) 									%
%											%
% REFERENCE  										%
% Wang, R., Martin, F. L., & Roth, F. (2003)						%
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		%
% based on edcmpf77_2.0/edcgrn.F							%
%		                                                                	%
% first created by lfeng Mon Feb 27 16:37:00 SGT 2012					%
% added stress calculation lfeng Fri May 18 17:29:13 SGT 2012                           %
% last modified by lfeng Sat May 19 16:10:02 SGT 2012                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% edgrn structure
obsz = edgrn.obsz;
minr = edgrn.minr;
maxr = edgrn.maxr;
minz = edgrn.minz;
maxz = edgrn.maxz;
dr   = edgrn.dr;
dz   = edgrn.dz;
% edgrnfcts strucutre
ssdisp0 = edgrnfcts.ssdisp;	% ssdisp0 is an old matlab function
ssstrn = edgrnfcts.ssstrn;
ssuzr  = edgrnfcts.ssuzr;
dsdisp = edgrnfcts.dsdisp;
dsstrn = edgrnfcts.dsstrn;
dsuzr  = edgrnfcts.dsuzr;
cldisp = edgrnfcts.cldisp;
clstrn = edgrnfcts.clstrn;
cluzr  = edgrnfcts.cluzr;
% convert Xin to EDCMP coordinate
% uniform depth for all observation points has been used in EDGRN to calculate green function library
% thus depth is ignored in EDCMP 
xrec = Xin(2,:);	% row vectors
yrec = Xin(1,:);	% row vectors
nrec = size(Xin,2);	

% discretized point sources
pxs_all     = pntsrc(:,1);
pys_all     = pntsrc(:,2);
pzs_all     = pntsrc(:,3);
pmoment_all = pntsrc(:,4:end);
nps = size(pntsrc,1);

% layer = [ id depth vp vs ro ]	(nn*5)
if ~isempty(layer)
    layer_num = size(layer,1);
    % one layer or obsz is within the bottom layer
    if layer_num==1 || obsz>layer(end,2)
        obs_vp = layer(end,3);
        obs_vs = layer(end,4);
        obs_ro = layer(end,5);
    else
        for ii=1:2:layer_num-1
            dp1 = layer(ii,2); dp2 = layer(ii+1,2);
            if obsz>=dp1 && obsz<=dp2 && dp1~=dp2
                vp1 = layer(ii,3); vp2 = layer(ii+1,3);
                vs1 = layer(ii,4); vs2 = layer(ii+1,4);
                ro1 = layer(ii,5); ro2 = layer(ii+1,5);
                obs_vp = vp1+(vp2-vp1)*(obsz-dp1)/(dp2-dp1);
        	obs_vs = vs1+(vs2-vs1)*(obsz-dp1)/(dp2-dp1); 
        	obs_ro = ro1+(ro2-ro1)*(obsz-dp1)/(dp2-dp1);
        	break
            end
        end
    end
    [mu,lambda,nu] = GTdef_vel2mod(obs_vp,obs_vs,obs_ro);
end

% initialize
disp   = zeros(nrec,3);
strain = zeros(nrec,6);
stress = zeros(nrec,6);
tilt   = zeros(nrec,2);

%%%%%%%%%%%%%%%%%%%%% SUPERPOSITION OF ALL DISCRETE POINT SOURCES %%%%%%%%%%%%%%%%%%%%%
for ii = 1:nrec
    % check distance 
    obsx = xrec(ii); obsy = yrec(ii);
    dis_all = sqrt((obsx-pxs_all).^2 +(obsy-pys_all).^2);
    ind_dis = dis_all>maxr;
    if any(ind_dis)
        error('GTdef_edcmp ERROR: green function grids do not cover all source-observation pairs!');
    end
    ind_dis = dis_all<minr;
    if any(ind_dis)
        error('GTdef_edcmp ERROR: green function grids do not cover all source-observation pairs!');
    end
    % azimuth measured from NORTH
    ind_dis = dis_all>0.0;
    azi_all(ind_dis) = atan2(obsy-pys_all(ind_dis),obsx-pxs_all(ind_dis)).*180/pi;
    ind_dis = dis_all==0.0;
    azi_all(ind_dis) = 0.0;

    % loop through all discretized point sources
    uzx=0.0; uzy=0.0;
    for jj = 1:nps
        pzs = pzs_all(jj);
        pmoment = pmoment_all(jj,:)';
        dis = dis_all(jj);	
        azi = azi_all(jj);
        % depth
        iz  = floor((pzs-minz)/dz);
        dzs = (pzs-(minz+dz*iz))/dz;
        iz  = iz+1;
        % distance
        idis = floor((dis-minr)/dr);
        ddis = (dis-(minr+dr*idis))/dr;
        idis = idis+1;
        % weighting factors for the interpolation
        w00 = (1.0-ddis)*(1.0-dzs);
        w10 = ddis*(1.0-dzs);
        w01 = (1.0-ddis)*dzs;
        w11 = ddis*dzs;
        co  = cosd(azi);
        si  = sind(azi);
        co2 = cosd(2.0*azi);
        si2 = sind(2.0*azi);
        % contributions from the strike-slip components
        ps  = pmoment(1)*si2+pmoment(4)*co2;
        sh  = pmoment(1)*co2-pmoment(4)*si2;
        uz  = ps*(w00*ssdisp0(1,idis,iz)+w10*ssdisp0(1,idis+1,iz)+w01*ssdisp0(1,idis,iz+1)+w11*ssdisp0(1,idis+1,iz+1));
        ur  = ps*(w00*ssdisp0(2,idis,iz)+w10*ssdisp0(2,idis+1,iz)+w01*ssdisp0(2,idis,iz+1)+w11*ssdisp0(2,idis+1,iz+1));
        ut  = sh*(w00*ssdisp0(3,idis,iz)+w10*ssdisp0(3,idis+1,iz)+w01*ssdisp0(3,idis,iz+1)+w11*ssdisp0(3,idis+1,iz+1));
        ezz = ps*(w00*ssstrn(1,idis,iz) +w10*ssstrn(1,idis+1,iz)+w01*ssstrn(1,idis,iz+1)+w11*ssstrn(1,idis+1,iz+1));
        err = ps*(w00*ssstrn(2,idis,iz) +w10*ssstrn(2,idis+1,iz)+w01*ssstrn(2,idis,iz+1)+w11*ssstrn(2,idis+1,iz+1));
        ett = ps*(w00*ssstrn(3,idis,iz) +w10*ssstrn(3,idis+1,iz)+w01*ssstrn(3,idis,iz+1)+w11*ssstrn(3,idis+1,iz+1));
        ezr = ps*(w00*ssstrn(4,idis,iz) +w10*ssstrn(4,idis+1,iz)+w01*ssstrn(4,idis,iz+1)+w11*ssstrn(4,idis+1,iz+1));
        ert = sh*(w00*ssstrn(5,idis,iz) +w10*ssstrn(5,idis+1,iz)+w01*ssstrn(5,idis,iz+1)+w11*ssstrn(5,idis+1,iz+1));
        etz = sh*(w00*ssstrn(6,idis,iz) +w10*ssstrn(6,idis+1,iz)+w01*ssstrn(6,idis,iz+1)+w11*ssstrn(6,idis+1,iz+1));
        uzr = ps*(w00*ssuzr(idis,iz)+w10*ssuzr(idis+1,iz)+w01*ssuzr(idis,iz+1)+w11*ssuzr(idis+1,iz+1));
        uzt = sh*(w00*ssdisp0(1,idis,iz)+w10*ssdisp0(1,idis+1,iz) ...
            + w01*ssdisp0(1,idis,iz+1)+w11*ssdisp0(1,idis+1,iz+1))*2.0/dis;
        % contributions from the dip-slip components
        ps  = pmoment(2)*co+pmoment(5)*si;
        sh  = pmoment(2)*si-pmoment(5)*co;
        uz  = uz +ps*(w00*dsdisp(1,idis,iz)+w10*dsdisp(1,idis+1,iz)+w01*dsdisp(1,idis,iz+1)+w11*dsdisp(1,idis+1,iz+1));
        ur  = ur +ps*(w00*dsdisp(2,idis,iz)+w10*dsdisp(2,idis+1,iz)+w01*dsdisp(2,idis,iz+1)+w11*dsdisp(2,idis+1,iz+1));
        ut  = ut +sh*(w00*dsdisp(3,idis,iz)+w10*dsdisp(3,idis+1,iz)+w01*dsdisp(3,idis,iz+1)+w11*dsdisp(3,idis+1,iz+1));
        ezz = ezz+ps*(w00*dsstrn(1,idis,iz)+w10*dsstrn(1,idis+1,iz)+w01*dsstrn(1,idis,iz+1)+w11*dsstrn(1,idis+1,iz+1));
        err = err+ps*(w00*dsstrn(2,idis,iz)+w10*dsstrn(2,idis+1,iz)+w01*dsstrn(2,idis,iz+1)+w11*dsstrn(2,idis+1,iz+1));
        ett = ett+ps*(w00*dsstrn(3,idis,iz)+w10*dsstrn(3,idis+1,iz)+w01*dsstrn(3,idis,iz+1)+w11*dsstrn(3,idis+1,iz+1));
        ezr = ezr+ps*(w00*dsstrn(4,idis,iz)+w10*dsstrn(4,idis+1,iz)+w01*dsstrn(4,idis,iz+1)+w11*dsstrn(4,idis+1,iz+1));
        ert = ert+sh*(w00*dsstrn(5,idis,iz)+w10*dsstrn(5,idis+1,iz)+w01*dsstrn(5,idis,iz+1)+w11*dsstrn(5,idis+1,iz+1));
        etz = etz+sh*(w00*dsstrn(6,idis,iz)+w10*dsstrn(6,idis+1,iz)+w01*dsstrn(6,idis,iz+1)+w11*dsstrn(6,idis+1,iz+1));
        uzr = uzr+ps*(w00*dsuzr(idis,iz)+w10*dsuzr(idis+1,iz)+w01*dsuzr(idis,iz+1)+w11*dsuzr(idis+1,iz+1));
        uzt = uzt+sh*(w00*dsdisp(1,idis,iz)+w10*dsdisp(1,idis+1,iz) ...
            + w01*dsdisp(1,idis,iz+1)+w11*dsdisp(1,idis+1,iz+1))*(-1.0/dis);
        % contributions from the clvd components
        ps  = pmoment(3);
        uz  = uz +ps*(w00*cldisp(1,idis,iz)+w10*cldisp(1,idis+1,iz)+w01*cldisp(1,idis,iz+1)+w11*cldisp(1,idis+1,iz+1));
        ur  = ur +ps*(w00*cldisp(2,idis,iz)+w10*cldisp(2,idis+1,iz)+w01*cldisp(2,idis,iz+1)+w11*cldisp(2,idis+1,iz+1));
        ezz = ezz+ps*(w00*clstrn(1,idis,iz)+w10*clstrn(1,idis+1,iz)+w01*clstrn(1,idis,iz+1)+w11*clstrn(1,idis+1,iz+1));
        err = err+ps*(w00*clstrn(2,idis,iz)+w10*clstrn(2,idis+1,iz)+w01*clstrn(2,idis,iz+1)+w11*clstrn(2,idis+1,iz+1));
        ett = ett+ps*(w00*clstrn(3,idis,iz)+w10*clstrn(3,idis+1,iz)+w01*clstrn(3,idis,iz+1)+w11*clstrn(3,idis+1,iz+1));
        ezr = ezr+ps*(w00*clstrn(4,idis,iz)+w10*clstrn(4,idis+1,iz)+w01*clstrn(4,idis,iz+1)+w11*clstrn(4,idis+1,iz+1));
        uzr = uzr+ps*(w00*cluzr(idis,iz)+w10*cluzr(idis+1,iz)+w01*cluzr(idis,iz+1)+w11*cluzr(idis+1,iz+1));

        % transform to cartesian coordinates
        disp(ii,1) = disp(ii,1)+ur*co-ut*si;
        disp(ii,2) = disp(ii,2)+ur*si+ut*co;
        disp(ii,3) = disp(ii,3)+uz;
 
        % strain order in edcgrn.F [exx eyy ezz exy eyz ezx]
        % adjust strain order      [exx eyy ezz eyz exz exy]
        strain(ii,1) = strain(ii,1)+err*co*co+ett*si*si-ert*si2;
        strain(ii,2) = strain(ii,2)+err*si*si+ett*co*co+ert*si2;
        strain(ii,3) = strain(ii,3)+ezz;	
        strain(ii,4) = strain(ii,4)+ezr*si+etz*co;
        strain(ii,5) = strain(ii,5)+ezr*co-etz*si;
        strain(ii,6) = strain(ii,6)+0.5*(err-ett)*si2+ert*co2;

        uzx = uzx+uzr*co-uzt*si;
        uzy = uzy+uzr*si+uzt*co;
	% transform of hrizontal tilts to vertical tilts
        tilt(ii,1) = 2.0*strain(ii,5)-uzx;
        tilt(ii,2) = 2.0*strain(ii,4)-uzy;
    end
end

% switch east & north
disp = disp(:,[2 1 3]);
tilt = tilt(:,[2 1]);
strain = strain(:,[2 1 3 5 4 6]);
stress = stress(:,[2 1 3 5 4 6]);
% reverse downside (no need for strain or stress)
disp(:,3) = -disp(:,3);

% from strain to stress
if ~isempty(layer)
    theta = strain(:,1)+strain(:,2)+strain(:,3);
    stress(:,1) = lambda.*theta+2*mu.*strain(:,1);
    stress(:,2) = lambda.*theta+2*mu.*strain(:,2);
    stress(:,3) = lambda.*theta+2*mu.*strain(:,3);
    stress(:,4) = 2*mu.*strain(:,4);
    stress(:,5) = 2*mu.*strain(:,5);
    stress(:,6) = 2*mu.*strain(:,6);
end
% from columnwise to rowwise, to be consistent with Okada disloc3d_mod2.m
disp = disp'; strain = strain'; stress = stress'; tilt = tilt';
