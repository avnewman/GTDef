function [ disp,strain,stress,tilt ] = GTdef_calc(M,Xin,earth,mu,nu,edgrn,layer,edgrnfcts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_calc.m				        %
%									  	%
% INPUT:                                                                        %
%  M - model parameters                                                         %
%  earth = homogeneous or layered                                               %
%    a. homogeneous                                                             %
%      M = [len;width;depth;dip;str;east;north;ss;ds;ts]                        %
%    b. layered                                                                 %
%      M = [slip;north;east;depth;length;width;str;dip;rake]                    %
%  Xin - point site locations in the local cartesian system 	  	        %
%        [n*3] [ xx; yy; zz ]						        %
%  (earth = homogeneous)                                                        %
%  mu  - shear modulus/rigidity							%
%  nu  - poisson's ratio                                                        %
%  (earth = layered)		        					%
%  edgrn, layer, and edgrnfcts                                                  %
%                                                                               %
% OUTPUT: disp, strain, stress, and tilt                                        %
% stress = [ SXX;SYY;SZZ;SYZ;SXZ;SXY ]    (6*n matrix)                          %
%                                                                               %
% first created by lfeng Fri May 18 16:23:56 SGT 2012                           %
% last modified by lfeng Mon Jun 11 15:53:20 SGT 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(earth,'homogeneous')
    tilt = [];
    [ disp,~,strain,stress,~,~ ] = disloc3d_mod2(M,Xin,mu,nu);
else
    pntsrc = [];
    Mnum = size(M,2);
    for ii=1:Mnum
        [ pntsrc0 ] = GTdef_edcmp_discretize(M(:,ii),edgrn);
        pntsrc = [ pntsrc; pntsrc0 ];
    end
    [ disp,strain,stress,tilt ] = GTdef_edcmp(edgrn,layer,edgrnfcts,pntsrc,Xin);
end
