function [ yy ] = GTdef_edgpsv(yy,kk,eps,sublayer,depths,source) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   GTdef_edgpsv.m                                      %
%                                                                                       %
% Spheroidal mode orthonormal propagator that calculates responses to P-SV sources      %
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
% This function only calculates U,V,E & F                                               %
%										        %
% REFERENCE  									        %
% Wang, R. (1999), A simple orthonormalization method for stable and efficient          %
% computation of Green’s functions, Bull. Seismol. Soc. Am., 89(3), 733–741.            %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5.           %
%		                                                                        %
% first created by Lujia Feng Tue Sep  1 15:24:03 SGT 2015                              %
% last modified by Lujia Feng Wed Jan 27 03:25:21 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps0 = 1e-3;
g0   = 0.0; % earth gravity in cm/s^2, no gravity is considered, otherwise g0 = 9.82e2 cm/s^2
y0   = zeros(4,2); % solution vectors at receiver depth
c0   = zeros(4,2); % 2 four-dimensional constant vector bases
yup  = zeros(4,2); % vector bases from surface to source
ylw  = zeros(4,2); % vector bases from halfspace to source

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
% problem more serious for P-SV waves than SH waves

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
      if (lup>lrec)  
         return; % negligible effect at receiver
      else
         break;
      end
   end
end

% determine the deepest layer in the half-space
llw = nl; % assume waves can reach the deepest layer
exponent = 0.0;
for ll = lsrc:nl-1
   exponent = exponent - kk*hh(ll);
   if(exponent<=log(eps))
     llw = ll;
     if(llw<lrec)  
        return; % negligible effect at receiver
     else
        break;
     end
   end
end

%===============================================================================
%       matrix propagation from surface to source
%===============================================================================
% define two orthonormal vector bases at surface
%      | U |   | 1 |    | 0 |
%  y - | E | = | 0 | or | 0 | 
%      | V |   | 0 |    | 1 |
%      | F |   | 0 |    | 0 |
if lup==1
   yup(1,1) = 1.0;
   yup(2,1) = ro(1)*g0*yup(1,1);
   yup(3,2) = 1.0;
else
   c0(1,1) = 1.0;
   c0(3,2) = 1.0;
   nn = nno(lup-1);
   % call edgmatrix(ma,4,k,0.0,n)
   [ ma ] = GTdef_edgmatrix(kk,0.0,vp(nn),vs(nn),ro(nn)); % L(0)
   yup = ma*c0; % equation 21 in Wang 1999
   % consider gravity
   dro = ro(nno(lup))-ro(nno(lup-1));
   yup(2,:) = yup(2,:)+dro*g0*yup(1,:);
end

if (lup==lrec) 
  y0 = yup;
end

for ll = lup+1:lsrc
   h0  = hh(ll-1);
   nn  = nno(ll-1);
   vp0 = vp(nn);
   vs0 = vs(nn);
   ro0 = ro(nn);
   % thin layer
   if (kk*h0<=eps0)
      % call edghask(hk,4,k,h0,n)
      [ hk ] = GTdef_edghask(4,kk,h0,vp0,vs0,ro0);
      yup    = hk*yup;
   % thick layer
   else
      % determine propagation matrix
      % call edgmatinv(mai,4,k,0.0,n)
      [ mai ] = GTdef_edgmatinv(kk,0.0,vp0,vs0,ro0);
      % call edgmatrix(ma,4,k,h0,n)
      [ ma  ] = GTdef_edgmatrix(kk,h0,vp0,vs0,ro0);
      c0      = mai*yup; % equation 21 in Wang 1999
      wave    = exp(-kk*h0);
      % normalize all modes
      mag1 = norm(c0(:,1));
      mag2 = norm(c0(:,2));
      fac  = 1.0/(mag1*mag2);
      % orthogonalize the p-sv modes
      orth(1,1) =  c0(3,2)*fac;
      orth(1,2) = -c0(1,2)*fac;
      orth(2,1) = -c0(3,1)*fac;
      orth(2,2) =  c0(1,1)*fac;
      c1        = c0*orth;
      if (ll>lrec)
         orth = orth*wave;
         y0   = y0*orth;
      end
%     c1(1,1) = c1(1,1)
      c1(2,1) = c1(2,1)*wave*wave;
      c1(3,1) = 0.0;
      c1(4,1) = c1(4,1)*wave*wave;
%
      c1(1,2) = 0.0;
      c1(2,2) = c1(2,2)*wave*wave;
%     c1(3,2)=c1(3,2)
      c1(4,2) = c1(4,2)*wave*wave;
      yup     = ma*c1;
   end
   dro = ro(nno(ll))-ro(nno(ll-1));
   yup(2,:)=yup(2,:)+dro*g0*yup(1,:);
   if (ll==lrec), y0 = yup; end
end

%===============================================================================
%       matrix propagation from half-space to source
%===============================================================================
% define two orthonormal vector bases at surface
%      | U |   | 0 |    | 0 |
%  y - | E | = | 1 | or | 0 | 
%      | V |   | 0 |    | 0 |
%      | F |   | 0 |    | 1 |
c0 = zeros(4,2); % 2 four-dimensional constant vector bases
nn = nno(llw);
% the lowest layer is fluid
if (vs(nn)<eps0*vp(nn))
   ylw(1,1) = 1.0;
   ylw(3,2) = 1.0;
% the lowest layer is solid
else
   c0(2,1) = 1.0;
   c0(4,2) = 1.0;
   % call edgmatrix(ma,4,k,0.d0,n)
   [ ma  ] = GTdef_edgmatrix(kk,0.0,vp(nn),vs(nn),ro(nn));
   ylw     = ma*c0;
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
   dro = ro(nno(ll+1))-ro(nno(ll));
   ylw(2,:) = ylw(2,:)-dro*g0*ylw(1,:);
   % thin layer
   if (kk*h0<=eps0)
      % call edghask(hk,4,k,-h0,n)
      [ hk ] = GTdef_edghask(4,kk,-h0,vp0,vs0,ro0);
      ylw = hk*ylw;
   % thick layer
   else
      % determine propagation matrix
      % call edgmatinv(mai,4,k,0.d0,n)
      [ mai ] = GTdef_edgmatinv(kk,0.0,vp0,vs0,ro0);
      % call edgmatrix(ma,4,k,-h0,n)
      [ ma  ] = GTdef_edgmatrix(kk,-h0,vp0,vs0,ro0);
      c0      = mai*ylw;
      wave    = exp(-kk*h0);
      % normalize all modes
      mag1 = norm(c0(:,1));
      mag2 = norm(c0(:,2));
      fac  = 1.0/(mag1*mag2);
      % orthogonalize the p-sv modes
      orth(1,1) =  c0(4,2)*fac;
      orth(1,2) = -c0(2,2)*fac;
      orth(2,1) = -c0(4,1)*fac;
      orth(2,2) =  c0(2,1)*fac;
      c1        = c0*orth;
      if (ll<lrec)
        orth = orth*wave;
	y0   = y0*orth;
      end
      c1(1,1) = c1(1,1)*wave*wave;
%     c1(2,1) = c1(2,1)
      c1(3,1) = c1(3,1)*wave*wave;
      c1(4,1) = 0.0;
      c1(1,2) = c1(1,2)*wave*wave;
      c1(2,2) = 0.0;
      c1(3,2) = c1(3,2)*wave*wave;
%     c1(4,2) = c1(4,2)
      ylw     = ma*c1;
   end
   if (ll>lsrc && ll==lrec), y0 = ylw; end
end

%===============================================================================
%       conditions on the source surface
%===============================================================================
%  source function
bb   = zeros(4,1);
coef = zeros(4,4);
for ii=1:4
   if (kpower(ii)==1)
      fac = kk;
   else
      fac = 1.0;
   end
   bb(ii) = sfct(ii)*fac; % 4-element vector
   for jj=1:2
      coef(ii,jj)   =  yup(ii,jj); % 4x4 matrix
      coef(ii,jj+2) = -ylw(ii,jj);
   end
end

% solve the two-point boundary problem
% linear equation system
% coef: coefficient matrix [nxn]
% bb:   right-hand matrix  [nx1]
% coef * xx = bb
xx = coef\bb;

% receiver is above source
if lrec<=lsrc
   for ii=1:4
      yy(ii) = 0.0;
      for jj=1:2
        yy(ii) = yy(ii) + xx(jj)*y0(ii,jj);
      end
   end
% receiver is below source
else
   for ii=1:4
      yy(ii) = 0.0;
      for jj=1:2
         yy(ii) = yy(ii) + xx(jj+2)*y0(ii,jj);
      end
   end
end
