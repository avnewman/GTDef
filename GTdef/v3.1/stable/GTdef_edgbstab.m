function [ bess ] = GTdef_edgbstab(srctype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_edgbstab.m                                     %
%                                                                                       %
% Calculate Bessel Function table based on fortran code EDGRN edgbstab.F                %
%                                                                                       %
% INPUT										        %
% srctype  - type of source                                                             %
%   'es' 0 = explosion point source (m11=m22=m33=M0)                                    %
%   'ss' 1 = strike-slip (m12=m21=M0)                                                   %
%   'ds' 2 = dip-slip (m13=m31=M0)                                                      %
%   'cl' 3 = compensated linear vector dipole (CLVD) (m11=m22=-M0/2, M33=M0)            %
%   'vs' 4 = vertical-single-force (fz=F0)                                              %
%   'hs' 5 = horizontal-single-force (fx=F0)                                            %
%   'gn' 6 = airgun in small pool (m11=m22=M0, fz=M0/r_source)                          %
%                                                                                       %
% in edgrn code, input parameter is nn = 3 - istype                                     %
%										        %
% OUTPUT                                                                                %
% bess structure                                                                        %
% nbess   = 2048                                                                        %
% ndbess  = 128                                                                         %
% nnbess  = nbess*ndbess                                                                %
% nnbess1 = nnbess+ndbess                                                               %
% bsdx    - step along x in besselj functions                                           %
% bsfct   - table of J_n(x), dJ_n(x)/dx and n*J_n(x)/x                                  %
%             all multiplied by sqrt(x)                                                 %
% Note: bsfct index starts from zero and length is nnbess1+1 because bsfct(0:nnbess1,3) %
%                                                                                       %
% REFERENCE  									        %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5	        %
%		                                                                        %
% first created by Lujia Feng Mon Aug 31 19:31:14 SGT 2015                              %
% used consistent srctype as GTdef_edgmoment.m lfeng Tue Sep  1 10:12:59 SGT 2015       %
% added bess structure lfeng Fri Jan 15 18:29:09 SGT 2016                               %
% changed bsfct because its index starts from zero in Fortran lfeng Fri Jan 22 SGT 2016 %
% last modified by Lujia Feng Fri Jan 22 11:52:13 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% double precision bsdx,bsfct(0:nnbess1,3)
% common /bessels/ bsdx,bsfct

nbess   = 2048;
ndbess  = 128;
nnbess  = nbess*ndbess;
nnbess1 = nnbess+ndbess;

bsfct   = zeros(nnbess1+1,3);

% Nyquist wavenumber
bsdx = pi*2.0/ndbess;

switch srctype
   case 'ss'
      % strike-slip (m12=m21=M0) istype = 1
      for ii=1:nnbess1
         xx            = bsdx*ii;
         xsqrt         = sqrt(xx);
         bsfct(ii+1,1) = xsqrt*besselj(2,xx);
         aa            = xsqrt*besselj(1,xx);
         bb            = xsqrt*besselj(3,xx);
         bsfct(ii+1,2) = 0.5*(aa-bb);
         bsfct(ii+1,3) = 0.5*(aa+bb);
      end

   case 'ds'
      % dip-slip (m13=m31=M0) istype = 2
      for ii=1:nnbess1
         xx            = bsdx*ii;
         xsqrt         = sqrt(xx);
         bsfct(ii+1,1) = xsqrt*besselj(1,xx);
         aa            = xsqrt*besselj(0,xx);
         bb            = xsqrt*besselj(2,xx);
         bsfct(ii+1,2) = 0.5*(aa-bb);
         bsfct(ii+1,3) = 0.5*(aa+bb);
      end

   case 'cl'
      % compensated linear vector dipole (CLVD) (m11=m22=-M0/2, M33=M0) istype = 3
      for ii=1:nnbess1
         xx            = bsdx*ii;
         xsqrt         = sqrt(xx);
         bsfct(ii+1,1) =  xsqrt*besselj(0,xx);
         bsfct(ii+1,2) = -xsqrt*besselj(1,xx);
         bsfct(ii+1,3) = 0.0;
      end

   otherwise
      error('GTdef_edgbstab ERROR: source type must be ss, ds or cl!');
end

bess.nbess   = nbess;
bess.ndbess  = ndbess;
bess.nnbess  = nnbess;
bess.nnbess1 = nnbess1;
bess.bsdx    = bsdx;
bess.bsfct   = bsfct;
