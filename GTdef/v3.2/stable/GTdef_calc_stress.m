function [ sspnt,ssflt1,ssflt2 ] = GTdef_calc_stress(sspnt,ssflt1,ssflt2,Min,earth) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_cal_stress                               %
% calculate stresses for points and fault planes                          %
%                                                                         %
% INPUT:                                                                  %
% structures: sspnt, ssflt1, ssflt2                                       %
% ----------------------------------------------------------------------- %
%  Min - subflt info in cartesian coordinate                              %
%    subflt order: along dip first, then along strike                     %
%    for Okada                                                            %
%       = [len width depth dip str east north ss ds ts]     [flt_num*10]  %
%    for layered model                                                    %
%       = [slip north east depth length width str dip rake] [flt_num*9]   %
% ----------------------------------------------------------------------- %
% Earth Structure:							  %
% Either of the two types of earth structure can be used 		  %
%  earth.type = 'homogeneous'                                             %
%      earth.rigidity		        (scalar)	{30e9 Pa}         %
%      earth.poisson		        (scalar)	{0.25}		  %
%  earth.type = 'layered'	        	        		  %
%      earth.edgrn.nl        	        (scalar)                          %
%      earth.edgrn.obsz     	        (scalar)                          %
%      earth.edgrn.nr	                (scalar)                          %
%      earth.edgrn.minr,edgrn.maxr      (scalar)                          %
%      earth.edgrn.nz                   (scalar)                          %
%      earth.edgrn.minz,edgrn.maxz      (scalar)                          %
%      earth.edgrn.srate                (scalar)			  %
%      earth.layer - [ id depth vp vs ro ]	(nn*5)			  %
%      earth.edgrnfcts - green function library                           %
%                                                                         %
% INTERMEDIATE: (rowwise to be consistent with Okada)                     %
% disp   = [ Ux;Uy;Uz ] (3 row vectors)					  %
% strain = [exx;eyy;ezz;eyz;exz;exy] (6 row vectors)                      %
% stress = [sxx;syy;szz;syz;sxz;sxy] (6 row vectors)                      %
% tilt (2 row vectors) 							  %
%                                                                         % 
% OUTPUT:                                                                 %
% structures: sspnt, ssflt1, ssflt2                                       %
%                                                                         %
% first created by Lujia Feng Mon Jun 11 15:09:22 SGT 2012                %
% added check outputs lfeng Wed Jun 13 18:54:07 SGT 2012                  %
% added Min0 check lfeng Wed Jun 13 19:18:34 SGT 2012                     %
% added earth structure lfeng Fri Mar 20 11:39:54 SGT 2015                %
% last modified by Lujia Feng Fri Mar 20 12:27:59 SGT 2015                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% redcue the size of Min0
%                   1   2     3     4   5   6    7     8  9  10
% for Okada Min0 = [len width depth dip str east north ss ds ts]
slip0 = 1e-5;      % threshold for slip calculation [m]
if strcmpi(earth.type,'homogeneous')
    keep_ind = Min0(:,8)>slip0 | Min0(:,9)>slip0 | Min0(:,10)>slip0;
    Min0 = Min0(keep_ind,:);
else
%                   1    2     3    4     5      6     7   8   9 
% for edcmp Min0 = [slip north east depth length width str dip rake]
    keep_ind = Min0(:,1)>slip0;
    Min0 = Min0(keep_ind,:);
end

% stress point
if sspnt.num~=0
    [ disp,strain,stress,tilt ] = ...
      GTdef_calc(Min0',sspnt.crt,earth);
    [ sspnt.shear,sspnt.normal,sspnt.coulomb ] = ...
      GTdef_calc_coulomb(sspnt.str,sspnt.dip,sspnt.rake,sspnt.fric,stress);
    %----------------------------- Debug -----------------------------
    % output strain to check
    fout = fopen('debug_strain','w');
    fprintf(fout,'#stress point name lon lat depth exx eyy ezz eyz exz exy\n');
    for ii=1:sspnt.num
        fprintf(fout,'stress point %s    %14.8f  %-12.8f %-6.4e  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',...
                sspnt.name{ii},sspnt.loc(ii,:),strain(:,ii));
    end
    fclose(fout);
    % output stress to check
    fout = fopen('debug_stress','w');
    fprintf(fout,'#stress point name lon lat depth sxx syy szz syz sxz sxy\n');
    for ii=1:sspnt.num
        fprintf(fout,'stress point %s    %14.8f  %-12.8f %-6.4e  %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',...
                sspnt.name{ii},sspnt.loc(ii,:),stress(:,ii));
    end
    fclose(fout);
    %-----------------------------------------------------------------
end

% stress fault1
if ssflt1.fltnum~=0
    [ disp,strain,stress,tilt ] = ...
      GTdef_calc(Min0',ssflt1.crt,earth);
    [ ssflt1.shear,ssflt1.normal,ssflt1.coulomb ] = ...
      GTdef_calc_coulomb(ssflt1.str,ssflt1.dip,ssflt1.rake,ssflt1.fric,stress);
end

% stress fault2
if ssflt2.fltnum~=0
    [ disp,strain,stress,tilt ] = ...
      GTdef_calc(Min0',ssflt2.crt,earth);
    [ ssflt2.shear,ssflt2.normal,ssflt2.coulomb ] = ...
      GTdef_calc_coulomb(ssflt2.str,ssflt2.dip,ssflt2.rake,ssflt2.fric,stress);
end
