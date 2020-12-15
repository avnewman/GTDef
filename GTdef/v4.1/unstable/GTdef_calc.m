function [ disp,strain,stress,tilt ] = GTdef_calc(earth,Min,Xin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                GTdef_calc.m				        %
%									  	%
% INPUT:                                                                        %
% earth structure                                                               %
% ----------------------------------------------------------------------------- %
%  Min - subflt info in cartesian coordinate                                    %
%  earth.type = homogeneous or layered                                          %
%    a. homogeneous                                                             %
%      Min = [len width depth dip str east north ss ds ts]     [10*flt_num]     %
%    b. layered                                                                 %
%      Min = [slip north east depth length width str dip rake] [9*flt_num]      %
% ----------------------------------------------------------------------------- %
%  Xin - point site locations in the local cartesian system 	  	        %
%        [3*n] [ xx yy zz ]						        %
%                                                                               %
% OUTPUT:                                                                       %
% disp   = [ Ux Uy Uz ]                   (3*nn)                                %
% strain = [ EXX EYY EZZ EYZ EXZ EXY ]    (6*nn)                                %
% stress = [ SXX SYY SZZ SYZ SXZ SXY ]    (6*nn)                                %
% tilt (2 row vectors) 							        %
%                                                                               %
% first created by Lujia Feng Fri May 18 16:23:56 SGT 2012                      %
% added earth structure lfeng Fri Mar 20 11:39:54 SGT 2015                      %
% changed input Xin lfeng Tue Mar 24 13:50:26 SGT 2015                          %
% last modified by Lujia Feng Wed Mar 25 00:55:21 SGT 2015                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

etype     = earth.type;
mu        = earth.rigidity;
nu        = earth.poisson;
edgrn     = earth.edgrn;
layer     = earth.layer;
edgrnfcts = earth.edgrnfcts;

% Xin needs to be [xx;yy;zz] for disloc3d_mod2.m & GTdef_edcmp.m
% Min needs to be [len;width;depth;dip;str;east;north;ss;ds;ts] for Okada or
% [slip;north;east;depth;length;width;str;dip;rake] for layered
Min = Min';
Xin = Xin';

if strcmpi(etype,'homogeneous')
    tilt = [];
    [ disp,~,strain,stress,~,~ ] = disloc3d_mod2(Min,Xin,mu,nu);
else
    pntsrc = [];
    Mnum = size(Min,2);
    for ii=1:Mnum
        [ pntsrc0 ] = GTdef_edcmp_discretize(Min(:,ii),edgrn);
        pntsrc = [ pntsrc; pntsrc0 ];
    end
    [ disp,strain,stress,tilt ] = GTdef_edcmp(edgrn,layer,edgrnfcts,pntsrc,Xin);
end

% convert disp   output from [ Ux;Uy;Uz ] (3*nn) to [ Ux Uy Uz ] (nn*3)
% convert strain output from [ EXX;EYY;EZZ;EYZ;EXZ;EXY ] (6*nn) to [ EXX EYY EZZ EYZ EXZ EXY ] (nn*6)
% convert stress output from [ SXX;SYY;SZZ;SYZ;SXZ;SXY ] (6*nn) to [ SXX SYY SZZ SYZ SXZ SXY ] (nn*6)
disp   = disp';
strain = strain';
stress = stress';
