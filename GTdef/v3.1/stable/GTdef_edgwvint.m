function [ uu ] = GTdef_edgwvint(rr,nr,srate,lambda,mu,sublayer,depths,source,bess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_edgwvint.m                                     %
%                                                                                       %
% Calculate point source Greens functions                                               %
% This code is converted from fortran code EDGRN edgwvint.F                             %
%                                                                                       %
% INPUT										        %
% (1) non-structure parameters                                                          %
% rr     - radial distances [m]                                      (1*nr row vector)  %  
% nr     - number of equidistant radial distances [m]	                      (scalar)  %
% srate  - sampling rate for wavenumber integration [10-128]                            %
%     Note: the larger the value is, the more accurate the results are                  %
% lambda - Lamé's first parameter                                                       %
% mu     - shear modulus/rigidity/Lamé's second parameter                               %
%                                                                                       %
% (2) earth.sublayer structure (model sublayers)                                        %
%     corresponds to /model/ h,ro,vp,vs,n0 [h->hh,n0->nl]                               %
% sublayer.nl    - number of sublayers                                       (scalar)   %
% sublayer.topz  - top depth of each sublayer [m]                (nl*1 column vector)   %
% sublayer.botz  - bottom depth of each sublayer [m]             (nl*1 column vector)   %
% sublayer.hh    - thickness of each sublayer [m]                (nl*1 column vector)   %
% sublayer.vp    - P-wave velocity for each sublayer [m/s]       (nl*1 column vector)   %
% sublayer.vs    - S-wave velocity for each sublayer [m/s]       (nl*1 column vector)   %
% sublayer.ro    - density for each sublayer [kg/m^3]            (nl*1 column vector)   %
%                                                                                       %
% (3) depths structure (changing layers including source and receiver)                  %
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
% (4) source structure (changing source)                                                %
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
%											%
% (5) bess structure                                                                    %
%     corresponds to /bessels/ bsdx,bsfct                                               %
% bess.nbess     = 2048                                                                 %
% bess.ndbess    = 128                                                                  %
% bess.nnbess    = nbess*ndbess                                                         %
% bess.nnbess1   = nnbess+ndbess                                                        %
% bess.bsdx      - step along x in besselj functions                                    %
% bess.bsfct     - table of J_n(x), dJ_n(x)/dx and n*J_n(x)/x                           %
%                  all multiplied by sqrt(x)                                            %
% Note: bsfct index starts from zero and length is nnbess1+1 because bsfct(0:nnbess1,3) %
%---------------------------------------------------------------------------------------%
% OUTPUT                                                                                %
% uu: 1=uz, 2=ur, 3=ut, 4=ezz, 5=err, 6=ett, 7=ezr, 8=ert, 9=etz, 10=duz/dr             %
% note:                                                                                 %
% uz, ur, ezz, err, ett, ezr, and duz/dr have the same azimuth-factor                   %
% as the poloidal mode (p-sv)                                                           %
% ut, ert, and etz have the same azimuth-factor as the toroidal mode (sh)               %
%										        %
% REFERENCE                                                                             %
% Wang, R., Martin, F. L., & Roth, F. (2003)                                            %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5            %
%		                                                                        %
% first created by Lujia Feng Tue Sep  1 12:55:48 SGT 2015                              %
% corrected an error that misplaced .^2 for yabs lfeng Tue Jan 26 18:20:12 SGT 2016     %
% corrected idint to floor instead of round lfeng Wed Jan 27 12:37:20 SGT 2016          %
% last modified by Lujia Feng Wed Jan 27 15:07:04 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps  = 1e-8;
eps0 = 1e-3;

% depth parameters
%nl   = depths.nl;
%hh   = depths.hh;
%nno  = depths.nno;
srcz = depths.srcz; % ~ srcz
recz = depths.recz; % ~ recz
%lsrc = depths.lsrc;
%lrec = depths.lrec;

% source parameters
r0     = source.r0;
ms     = source.ms;   
cs     = source.cs;
sfct   = source.sfct;  
%kpower = source.kpower;

% bessel function parameters
%nbess   = bess.nbess;
ndbess  = bess.ndbess;
nnbess  = bess.nnbess;
%nnbess1 = bess.nnbess1;
bsdx    = bess.bsdx;
bsfct   = bess.bsfct;

% check the importance of ps (poloidal) mode
fps = sum(abs(sfct(1:4)));

if (fps>0.0)
   psflag = true;
else
   psflag = false;
end

% check the importance of sh (toroidal) mode
fsh = sum(abs(sfct(5:6)));
if (fsh>0.0)
   shflag = true;
else
   shflag = false;
end

% initialization
r00 = zeros(nr,1);
for ii=1:nr
   r00(ii) = max([r0 1.0e-02*abs(srcz-recz) 1.0e-02*rr(ii)]);
end
uu  = zeros(10,nr);
u0  = zeros(6,1);
uk0 = zeros(6,1);
y0  = zeros(6,1);

% determine wavenumber limit
% determine limits of y(i), i=1,...,6
ncall = 0;
if (srcz==recz) % source is at the same depth as the receiver
   k0   = eps0*2*pi/(r0+abs(srcz-recz)+rr(nr));
   ymax = 0.0;
   yabs = 1.0;
   % calculate k0
   while (yabs>eps0*ymax)
      % call edgkern(y,k0,ps,sh,eps)
      [ yy ] = GTdef_edgkern(psflag,shflag,k0,eps,sublayer,depths,source);
      ncall  = ncall + 1;
      yabs   = sum(yy([1 3 5]).^2);
      yabs   = k0*sqrt(k0*yabs)*exp(-(k0*r0)^2);
      ymax   = max(ymax,yabs);
      k0     = 1.25*k0;
   end
   analflag  = true;
   kk        = eps0*2*pi/(r0+abs(srcz-recz)+rr(nr));
   yabs      = 0.0;
   dyabs     = 1.0;
   while (dyabs>eps*yabs)
      % call edgkern(y,kk,ps,sh,eps)
      [ yy ] = GTdef_edgkern(psflag,shflag,kk,eps,sublayer,depths,source);
      ncall  = ncall + 1;
      yy([2 4 6]) = yy([2 4 6])./(kk*mu);
      yabs   = sum(yy.^2);
      dyabs  = sum((yy-y0).^2);
      y0     = yy;
      if (kk>=k0)
         analflag = false;
         y0  = zeros(6,1);
         break;
      end
      kk = 1.25*kk;
   end
   y0([2 4 6]) = y0([2 4 6])*mu;
else
   analflag = false;
end

% find out klimit
klimit = eps*2*pi/(r0+abs(srcz-recz)+rr(nr));
ymax = 0.0;
while 1
   % call edgkern(y,klimit,ps,sh,eps)
   [ yy ] = GTdef_edgkern(psflag,shflag,klimit,eps,sublayer,depths,source);
   ncall  = ncall + 1;
   yabs   = sum((yy([1 3 5])-y0([1 3 5])).^2); % corrected .^2 error
   yabs   = klimit*sqrt(klimit*yabs)*exp(-(klimit*r00(1))^2);
   ymax   = max(ymax,yabs);
   if (yabs>eps0*ymax)
      klimit = 1.2*klimit;
   else
      break;
   end
end

% determine wavenumber sampling rate
dk = 2*pi/(srate*(r00(1)+rr(nr))+abs(srcz-recz));
nk = 500 + round(klimit/dk); % idnint in fortran = round in matlab
dk = klimit/nk;

% too small distances will be treated as r = 0!
if (rr(1)>0.0)
  ir1 = 1;
else
  ir1 = 2;
end

% Hanker Transformation integrate through wavenumber
for ik=1:nk
   kk = ik*dk;
   % call edgkern(y,k,ps,sh,eps)
   [ yy ] = GTdef_edgkern(psflag,shflag,kk,eps,sublayer,depths,source);
   if(analflag)
      yy([1 3 5]) = yy([1 3 5]) - y0([1 3 5]);
      yy([2 4 6]) = yy([2 4 6]) - y0([2 4 6])*kk;
   % for r=0
   elseif (ir1==2)
      fac = kk*exp(-(kk*r00(1))^2)*dk;
      u0  = u0  + yy*fac;
      uk0 = uk0 + yy*kk*fac;
   end
   for ir=ir1:nr
      fac = sqrt(kk)*dk*exp(-(kk*r00(ir))^2)/sqrt(rr(ir));
      xx  = kk*rr(ir);
      % bessels functions from pre-calculated tables
      nx  = floor(xx/bsdx); % idint in fortran = floor in matlab
      wr  = xx/bsdx-nx;
      wl  = 1.0-wr;
      if (nx>nnbess)
         nx    = nnbess + mod(nx-nnbess,ndbess);
         bs    = fac*(wl*bsfct(nx+1,:)+wr*bsfct(nx+2,:)); % bsfct index starts from 0 for Fortran bsfct(0:nnbess1,3) but 1 for Matlab 
         bs(3) = bs(3)*(nx+wr)*bsdx/xx;
      else
         bs    = fac*(wl*bsfct(nx+1,:)+wr*bsfct(nx+2,:));
      end
      % u1  = uz, displacement component
      % u2  = ur, displacement component
      % u3  = ut, displacement component
      % u4  = normal stress: szz
      % u5  = surface strain: err+ett
      % u6    will be derived later
      % u7  = shear stress: szr
      % u8  = strain component: dut/dr - (dur/dt)/r + ut/r
      % u9  = shear stress: szt
      % u10 = tilt: duz/dr

      uu(1,ir)  = uu(1,ir)  + yy(1)*bs(1);
      uu(2,ir)  = uu(2,ir)  + yy(3)*bs(2)    + cs*yy(5)*bs(3);
      uu(3,ir)  = uu(3,ir)  - cs*yy(3)*bs(3) - yy(5)*bs(2);
      uu(4,ir)  = uu(4,ir)  + yy(2)*bs(1);
      uu(5,ir)  = uu(5,ir)  - yy(3)*kk*bs(1);
      uu(7,ir)  = uu(7,ir)  + yy(4)*bs(2)    + cs*yy(6)*bs(3);
      uu(8,ir)  = uu(8,ir)  + yy(5)*kk*bs(1);
      uu(9,ir)  = uu(9,ir)  - cs*yy(4)*bs(3) - yy(6)*bs(2);
      uu(10,ir) = uu(10,ir) + yy(1)*kk*bs(2);
   end
end
% end of total integral
fprintf('wavenumber samples: %d, really used %d\n',ncall+nk,nk);

% for very small rr including rr=0
if (ir1==2 && ~analflag)
   if (ms==0)
      uu(1,1) = u0(1);
      uu(4,1) = u0(2);
      uu(5,1) = -0.50*uk0(3);
      uu(6,1) = uu(5,1);
   elseif (ms==1)
      uu(2,1)  =  0.5*(u0(3)+cs*u0(5));
      uu(3,1)  = -0.5*(cs*u0(3)+u0(5));
      uu(7,1)  =  0.5*(u0(4)+cs*u0(6));
      uu(9,1)  = -0.5*(cs*u0(4)+u0(6));
      uu(10,1) =  0.5*uk0(1);
   elseif (ms==2)
      uu(5,1)  =  0.25*(uk0(3)+cs*uk0(5));
      uu(6,1)  = -uu(5,1);
      uu(8,1)  = -0.25*(cs*uk0(3)+uk0(5));
   end
end

% ir1 might be 1 or 2
for ir=ir1:nr
   if(analflag)
      if (ms==0)
         bs(1)  =  0.0;
         bs(2)  = -1.0/rr(ir)^2;
         bs(3)  =  0.0;
         bs1(1) = -1.0/rr(ir)^3;
         bs1(2) =  0.0;
         bs1(3) =  0.0;
      elseif (ms==1)
         bs(1)  =  1.0/rr(ir)^2;
         bs(2)  = -1.0/rr(ir)^2;
         bs(3)  =  1.0/rr(ir)^2;
         bs1(1) =  0.0;
         bs1(2) = -2.0/rr(ir)^3;
         bs1(3) =  1.0/rr(ir)^3;
      elseif (ms==2)
         bs(1)  =  2.0/rr(ir)^2;
         bs(2)  = -1.0/rr(ir)^2;
         bs(3)  =  2.0/rr(ir)^2;
         bs1(1) =  3.0/rr(ir)^3;
         bs1(2) = -4.0/rr(ir)^3;
         bs1(3) =  4.0/rr(ir)^3;
      end
      uu(1,ir)  = uu(1,ir)  + y0(1)*bs(1);
      uu(2,ir)  = uu(2,ir)  + y0(3)*bs(2)  + cs*y0(5)*bs(3);
      uu(3,ir)  = uu(3,ir)  - y0(5)*bs(2)  - cs*y0(3)*bs(3);
      uu(4,ir)  = uu(4,ir)  + y0(2)*bs1(1);
      uu(5,ir)  = uu(5,ir)  - y0(3)*bs1(1);
      uu(7,ir)  = uu(7,ir)  + y0(4)*bs1(2) + cs*y0(6)*bs1(3);
      uu(8,ir)  = uu(8,ir)  + y0(5)*bs1(1);
      uu(9,ir)  = uu(9,ir)  - y0(6)*bs1(2) - cs*y0(4)*bs1(3);
      uu(10,ir) = uu(10,ir) + y0(1)*bs1(2);
   end
   % u6 is ett = ur/r + (dut/dt)/r
   uu(6,ir) = (uu(2,ir)+cs*ms*uu(3,ir))/rr(ir);
   
   % u5 now is err = u5(before) - ett
   uu(5,ir) = uu(5,ir) - uu(6,ir);
   
   % u8 now is ert = 0.5 * u8(before) + (dur/dt)/r - ut/r
   %               = 0.5 * (dut/dr + (dur/dt)/r - ut/r)
   uu(8,ir) = 0.5*uu(8,ir) - (cs*ms*uu(2,ir) + uu(3,ir))/rr(ir);
end
% ir1 might be 1 or 2
% u6 is ett = ur/r + (dut/dt)/r
%hehe uu(6,[ir1 end]) = (uu(2,[ir1 end]) + cs*ms*uu(3,[ir1 end]))./rr([ir1 end])';

% u5 now is err = u5(before) - ett
%hehe uu(5,[ir1 end]) = uu(5,[ir1 end]) - uu(6,[ir1 end]);

% u8 now is ert = 0.5 * u8(before) + (dur/dt)/r - ut/r
%               = 0.5 * (dut/dr + (dur/dt)/r - ut/r)
%hehe uu(8,[ir1 end]) = 0.5*uu(8,[ir1 end])-(cs*ms*uu(2,[ir1 end])+uu(3,[ir1 end]))./rr([ir1 end])';

% from stress to strain
uu(4,:) = (uu(4,:)-lambda.*(uu(5,:)+uu(6,:)))./(lambda+2.0*mu);
uu(7,:) = uu(7,:)/(2.0*mu);
uu(9,:) = uu(9,:)/(2.0*mu);
