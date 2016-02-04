function [ yy ] = GTdef_edgsh(yy,kk,eps,sublayer,depths,source) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    GTdef_edgsh.m                                      %
%                                                                                       %
% Toroidal mode Haskell propagator that calculates responses to SH sources              %
%                                                                                       %
% INPUT										        %
% (1) kk  - wavenumber                                                                  %
% (2) eps - relative accuracy                                                           %
% (3) earth.sublayer structure (model sublayers)                                        %
%     corresponds to /model/ h,ro,vp,vs,n0 [h->hh,n0->nl]                               %
% sublayer.nl    - number of sublayers                                       (scalar)   %
% sublayer.topz  - top depth of each sublayer [m]                (nl*1 column vector)   %
% sublayer.botz  - bottom depth of each sublayer [m]             (nl*1 column vector)   %
% sublayer.hh    - thickness of each sublayer [m]                (nl*1 column vector)   %
% sublayer.vp    - P-wave velocity for each sublayer [m/s]       (nl*1 column vector)   %
% sublayer.vs    - S-wave velocity for each sublayer [m/s]       (nl*1 column vector)   %
% sublayer.ro    - density for each sublayer [kg/m^3]            (nl*1 column vector)   %
%											%
% (4) depths structure (changing layers including source and receiver)                  %
%     corresponds to /sublayer/ hp,lp,nno    [hp->hh,lp->nl]                            %
%                    /receiver/ zrec,lzrec [zrec->recz,lzrec->lrec]                     %
%                     /source/  zs,ls      [zs->srcz,ls->lsrc]                          %
% depths.nl        - the number of depths including source depth and receiver depth     %
% depths.hh        - thickness of each depth                                            %
% depths.nno       - layer number for each depth                                        %
% depths.srcz      - depth of source                                                    %
% depths.recz      - depth of receiver                                                  %
% depths.lsrc      - depth number for the source                                        %   
% depths.lrec      - depth number for the receiver                                      %   
%											%
% (5) source structure (changing source)                                                %
%     corresponds to /source/ r0,ms,ics,sfct,kpower [ics->cs]                           %
% source.r0      - point source scale [m] (scalar)                                      %
% source.ms      - (scalar)                                                             %
% source.cs      - (scalar) is ics (integer) in edgmoment.F                             %
%                = 1  the azmuth-factor is cos(ms*theta) for poloidal mode (P-SV waves) %
%                     and sin(ms*theta) for toroidal mode (SH wave)                     %
%                = -1 otherwise                                                         %
% source.sfct    - six depth-dependent coefficients                                     %
%                  the first four (Um,Em,Vm,Fm) are in the poloidal mode                %
%                  the last two (Wm,Gm) are in the toroidal mode                        %
% source.kpower  -                     (6 integers)                                     %
%---------------------------------------------------------------------------------------%
% OUTPUT                                                                                %
% yy      - solution vector for P-SV waves                                              %
%         = (U,E,V,F,W,G)^T                                                             %
%   U,V,W - two displacement components in the frequency-wavenumber domain              %
%   P,S,G - two components of elastic surface force on a horizontal plane               %
% see Equation 15 & 20 of Wang et al. 2003                                              %
%										        %
% This function only calculates W & G                                                   %
%										        %
% REFERENCE  									        %
% Wang, R. (1999), A simple orthonormalization method for stable and efficient          %
% computation of Green’s functions, Bull. Seismol. Soc. Am., 89(3), 733–741.            %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5	        %
%		                                                                        %
% first created by Lujia Feng Tue Sep  1 18:21:45 SGT 2015                              %
% found typos associated with ylw lfeng Tue Jan 26 11:21:35 SGT 2016                    %
% last modified by Lujia Feng Tue Jan 26 18:25:25 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps0 = 1e-3;
%g0   = 0.0; % earth gravity in cm/s^2, no gravity is considered, otherwise g0 = 9.82e2 cm/s^2
y0   = zeros(2,1); % solution vectors at receiver depth
%yup  = zeros(2,1); % vector bases from surface to source
ylw  = zeros(2,1); % vector bases from halfspace to source

% layer parameters
%nl   = sublayer.nl;
%topz = sublayer.topz;
%botz = sublayer.botz;
%hh   = sublayer.hh;
vp   = sublayer.vp;
vs   = sublayer.vs;
ro   = sublayer.ro;

% depth parameters
nl   = depths.nl;
hh   = depths.hh;
nno  = depths.nno;
%srcz = depths.srcz; % ~ srcz
%recz = depths.recz; % ~ recz
lsrc = depths.lsrc;
lrec = depths.lrec;

% source parameters
%r0     = source.r0;
%ms     = source.ms;   
%cs     = source.cs;   
sfct   = source.sfct;  
kpower = source.kpower;

%===============================================================================
%       matrix propagation preparation
%===============================================================================
% to avoid numerical instability when the product between the wavenumber and layer thickness is large
% problem less serious for SH waves than P-SV waves

% determine the shallowest layer close to surface
lup = 1; % assume waves can reach the surface
exponent = 0.0;
for ll = lsrc-1:-1:1
   % decreasing waves propagating as e^(-kz) 
   % see Section 2.4 of Wang et al. 2003
   % if kz is too large (==-kz too small), then the evanescent reflection cannot be observed
   exponent = exponent - kk*hh(ll);
   if (exponent<=log(eps))
      lup = ll+1;
      if (lup>lrec), return; end % negligible effect at receiver
   end
end

% determine the deepest layer in the half-space
llw = nl; % assume waves can reach the deepest layer
exponent = 0.0;
for ll = lsrc:nl-1
   exponent = exponent - kk*hh(ll);
   if (exponent<=log(eps))
     llw = ll;
     if (llw<lrec), return; end % negligible effect at receiver
   end
end

%===============================================================================
%       matrix propagation from surface to source
%===============================================================================
yup = [ 1.0; 0.0 ]; % initial vector base at surface
if lup>1
   nn     = nno(lup-1);
   yup(2) = ro(nn)*vs(nn)^2*kk;
end

if (lup==lrec)  
   y0 = yup; 
end

for ll=lup+1:lsrc
   h0  = hh(ll-1);
   nn  = nno(ll-1);
   vp0 = vp(nn);
   vs0 = vs(nn);
   ro0 = ro(nn);
   % determine propagation matrix
   % call edghask(hk,2,k,hh(ll-1),n)
   [ hk ] = GTdef_edghask(2,kk,h0,vp0,vs0,ro0);
   y1     = hk*yup;
   wave   = exp(-kk*h0);
   yup    = y1*wave;
   if (ll>lrec)
      y0  = y0*wave;
   elseif (ll==lrec)
      y0  = yup;
   end
end

%===============================================================================
%       matrix propagation from half-space to source
%===============================================================================

nn = nno(llw);
ylw(1) = 1.0;
% the lowest layer is fluid
if (vs(nn)<eps0*vp(nn))
   ylw(2) = 0.0;
% the lowest layer is solid
else
   ylw(2) = -ro(nn)*vs(nn)^2*kk;
end

if (llw>lsrc && llw==lrec) 
   y0 = ylw;
end

for ll=llw-1:-1:lsrc
   h0  = hh(ll);
   nn  = nno(ll);
   vp0 = vp(nn);
   vs0 = vs(nn);
   ro0 = ro(nn);
   % determine propagation matrix
   % call edghask(hk,2,k,-hh(l),n)
   [ hk ] = GTdef_edghask(2,kk,-h0,vp0,vs0,ro0);
   y1     = hk*ylw;
   wave   = exp(-kk*h0);
   ylw    = y1*wave;
   if (ll<lrec)
      y0  = y0*wave;
   elseif (ll>lsrc && ll==lrec)
      y0  = ylw;
   end
end

%===============================================================================
%       conditions on the source surface
%===============================================================================
% source function
bb   = zeros(2,1);
coef = zeros(2,2);
for ii=1:2
   if (kpower(ii+4)==1)
      fac = kk;
   else
      fac = 1.0;
   end
   bb(ii) = sfct(ii+4)*fac; % 2-element vector
   coef(ii,1) =  yup(ii);
   coef(ii,2) = -ylw(ii);
end

% solve the two-point boundary problem
% linear equation system
% coef: coefficient matrix [nxn]
% bb:   right-hand matrix  [nx1]
% coef * xx = bb
xx = coef\bb;

% receiver is above source
if lrec<=lsrc
  for ii=1:2
    yy(ii+4) = xx(1)*y0(ii);
  end
% receiver is below source
else
  for ii=1:2
    yy(ii+4) = xx(2)*y0(ii);
  end
end
