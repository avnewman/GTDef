function [ source ] = GTdef_edgmoment(srctype,vp,vs,ro)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 GTdef_edgmoment.m                                     %
%                                                                                       %
% Calculate source functions based on fortran code EDGRN edgmoment.F                    %
% change exponential ** in Fortran to ^ in MATLAB                                       %
%                                                                                       %
% INPUT										        %
% srctype - type of source                                                              %
%    'es' - explosion source                                                            %
%    'ss' - strike-slip                                                                 %
%    'ds' - dip-slip                                                                    %
%    'cl' - compensated linear vector dipole (CLVD)                                     %
%    'vs' - vertical-single-force (fz=F0)                                               %
%    'hs' - horizontal-single-force (fx=F0)                                             %
%    'gn' - airgun in small pool (m11=m22=M0, fz=M0/r_source)                           %
% vp      - P-wave velocity of the source                                               %
% vs      - S-wave velocity of the source                                               %
% ro      - density of the source                                                       %
%										        %
% OUTPUT                                                                                %
% source structure                                                                      %
% r0        - point source scale [m] (scalar)                                           %
% ms        - (scalar)                                                                  %
% cs        - (scalar) is ics (integer) in edgmoment.F                                  %
%	    = 1  the azmuth-factor is cos(ms*theta) for poloidal mode (P-SV waves)      %
%                and sin(ms*theta) for toroidal mode (SH wave)                          %
%	    = -1 otherwise                                                              %
% sfct      - six depth-dependent coefficients                                          %
%             the first four (Um,Em,Vm,Fm) are in the poloidal mode                     %
%             the last two (Wm,Gm) are in the toroidal mode                             %
% kpower    -                     (6 integers)                                          %
%										        %
% REFERENCE  									        %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5	        %
%		                                                                        %
% first created by Lujia Feng Wed Dec 10 15:44:04 SGT 2014                              %
% last modified by Lujia Feng Tue Sep  1 13:18:59 SGT 2015                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strength = ro*vs*vs;  % moment of the source

sfct   = zeros(1,6);
kpower = zeros(1,6);

switch srctype
   case 'es'
      % explosion source (m11=m22=m33=M0)
      ms        = 0;
      cs        = 1.0;
      sfct(1)   = -strength/(2.0*pi*ro*vp*vp);
      sfct(4)   = -strength*(vs/vp)^2/pi;
      kpower(4) = 1;
   case 'ss'
      % strike-slip (m12=m21=M0)
      ms        = 2;
      cs        = -1.0;
      sfct(4)   = strength/(2.0*pi);
      sfct(6)   = -sfct(4);
      kpower(4) = 1;
      kpower(6) = 1;
   case 'ds'
      % dip-slip (m13=m31=M0)
      ms        = 1;
      cs        = 1.0;
      sfct(3)   = -strength/(2.0*pi*ro*vs*vs);
      sfct(5)   = sfct(3);
   case 'cl'
      % compensated linear vector dipole (CLVD) (m11=m22=-M0/2, M33=M0)
      ms        = 0;
      cs        = 1.0;
      sfct(1)   = -strength/(2.0*pi*ro*vp*vp);
      sfct(4)   = strength*(3.0-4.0*(vs/vp)^2)/(4.0*pi);
      kpower(4) = 1;
   case 'vs'
      % vertical-single-force (fz=F0)
      ms        = 0;
      cs        = 1.0;
      sfct(2)   = strength/(2.0*pi);
   case 'hs'
      % horizontal-single-force (fx=F0)
      ms        = 1;
      cs        = 1.0;
      sfct(4)   = strength/(2.0*pi);
      sfct(6)   = sfct(4);
%  case 'gn'
%     % airgun in small pool (m11=m22=M0, fz=M0/r_source)
%     ms=0
%     cs=1
%     sfct(2)=-strength/(2.d0*pi*r0)  r0 not provided
%     sfct(4)=-strength/(2.d0*pi)
%     kpower(4)=1
   otherwise
      error('GTdef_edgmoment ERROR: source type must be es, ss, ds, cl, vs, or hs!');
end

source.ms     = ms;
source.cs     = cs;
source.sfct   = sfct;
source.kpower = kpower;
