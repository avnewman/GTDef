function [ edgrnfcts ] = GTdef_edgrn(edgrn,layer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           	     	GTdef_edgrn.m					%
%											%
% Calculate point sources using EDGRN						        %
%											%
% INPUT											%
% (1) edgrn structure									%
% edgrn.obsz  - uniform observation depth                                    (scalar)   %
% edgrn.nr    - number of equidistant radial distances [m]	             (scalar)   %
% edgrn.minr  - minimum radial distance [m]                                  (scalar)   %
% edgrn.maxr  - maximum radial distance [m]                                  (scalar)   %
% edgrn.nz    - number of equidistant source depths [m]                      (scalar)   %
% edgrn.minz  - minimum source depths [m]                                    (scalar)   %
% edgrn.maxz  - maximum source depths [m]                                    (scalar)   %
% edgrn.srate - sampling rate for wavenumber integration [10-128]                       %
%     Note: the larger the value is, the more accurate the results are                  %
% edgrn.dr    - radial distance step [m]                                     (scalar)   %
% edgrn.dz    - source depth step [m]                                        (scalar)   %
% edgrn.rr    - radial distances [m]                                (1*nr row vector)   %  
% edgrn.zz    - source depths [m]                                (nz*1 column vector)   %  
%											%
% edgrn.nl    - number of layers                                             (scalar)   %
% edgrn.z1    - top depth of each layer [m]                      (nl*1 column vector)   %
% edgrn.z2    - bottom depth of each layer [m]                   (nl*1 column vector)   %
% edgrn.ro    - density for each source depth [kg/m^3]           (nl*1 column vector)   %
% edgrn.vp    - P-wave velocity for each source depth [m/s]      (nl*1 column vector)   %
% edgrn.vs    - S-wave velocity for each source depth [m/s]      (nl*1 column vector)   %
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
%		                                                                	%
% first created by Lujia Feng Wed Dec 10 15:31:46 SGT 2014                              %
% last modified by Lujia Feng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine edgwvint(u,r,nr,srate,lambda,mu,itty)
        integer nr,itty
        double precision srate,lambda,mu
        double precision r(nrmax)
        double precision u(10,nrmax)
%
% 	table of J_n(x), dJ_n(x)/dx and n*J_n(x)/x
%	all multiplied by sqrt(x)
%
	double precision bsdx,bsfct(0:nnbess1,3)
	common /bessels/ bsdx,bsfct
%
        double precision eps,eps0
        data eps,eps0/1.0d-08,1.0d-03/
%
        integer lp,nno(nzmax)
        double precision hp(nzmax)
        common /sublayer/ hp,lp,nno

        integer lzrec
        double precision zrec
        common /receiver/ zrec,lzrec
%       n0: number of model layers
        integer n0
        double precision h(lmax),ro(lmax),vp(lmax),vs(lmax)
        common /model/ h,ro,vp,vs,n0
%       source parameters
        integer ls,ms,ics
        integer kpower(6)
        double precision zs,r0
        double precision sfct(6)
        common /source/ zs,r0,sfct,ls,ms,ics,kpower
c
%       parameters for hankel integrations
c
        integer i,ir,ir1,ncall,ik,nk,nx
        double precision k,k0,klimit,dk,x,wl,wr
        double precision cs,fps,fsh,fac,yabs,dyabs,ymax
        double precision y(6),y0(6),u0(6),uk0(6),bs(3),bs1(3)
	double precision r00(nrmax)
        logical ps,sh,analytic
% u: 1=uz, 2=ur, 3=ut, 4=ezz, 5=err, 6=ett, 7=ezr, 8=ert, 9=etz, 10=duz/dr
% NOTE: uz, ur, ezz, err, ett, ezr duz/dr have the same azimuth-factor as the poloidal mode (p-sv)
%	ut, ert and etz have the same azimuth-factor as the toroidal mode (sh)

% edgrn structure
nl    = edgrn.nl;
obsz  = edgrn.obsz; % zrec in EDGRN
nr    = edgrn.nr;
minr  = edgrn.minr; 
maxr  = edgrn.maxr;
nz    = edgrn.nz;  
minz  = edgrn.minz; 
maxz  = edgrn.maxz;
srate = edgrn.srate;
dr    = edgrn.dr;
dz    = edgrn.dz;

% test if poloidal and toroidal modes are significant
if sum(abs(sfct(1:4)))>0.0
  ps = true; 
else
  ps = false;
end
if sum(abs(sfct(5:6)))>0.0
  sh = true;
else
  sh = false;
end

% loop through depth
for ii=1:nz
% initialization
do ir=1,nr
  r00(ir)=dmax1(r0,0.01*abs(zs-zrec),0.01*r(ir))
end
uu  = zeros(10,nr); % 10 observables to output
u0  = zeros(1,6);
uk0 = zeros(1,6);
y0  = zeros(1,6);

% determine wavenumber limit
% determine limits of y(i), i=1,...,6
r0   = 0.5*max(dz,dr);
eps0 = 0.003;
ncall=0
if zs==zrec
   k0   = 2*eps0*pi/(r0+abs(zs-zrec)+r(nr));
   ymax = 0.0;
   yabs = 0.0;
   while(yabs>eps0*ymax)
      call edgkern(y,k0,ps,sh,eps)
      ncall = ncall+1
      do i=1,5,2
        yabs=yabs+y(i)*y(i)
      enddo
      yabs=k0*sqrt(k0*yabs)*exp(-(k0*r0)**2);
      ymax=max(ymax,yabs);
      k0=1.25d0*k0;
   end

   analytic = true;

   k=eps0*pi2/(r0+dabs(zs-zrec)+r(nr))
   call edgkern(y,k,ps,sh,eps)
   ncall=ncall+1
   do i=2,6,2
     y(i)=y(i)/(k*mu)
   enddo
   yabs=0.d0
   dyabs=0.d0
   do i=1,6
     yabs=yabs+y(i)**2
     dyabs=dyabs+(y(i)-y0(i))**2
     y0(i)=y(i)
   enddo
   if(dyabs.gt.eps*yabs)then
     if(k.ge.k0)then
       analytic=.false.
       do i=1,6
         y0(i)=0.d0
       enddo
     else
       k=1.25d0*k
       goto 20
     endif
   endif
   do i=2,6,2
     y0(i)=y0(i)*mu
   enddo
else
   analytic = false;
end
%
	klimit=eps*pi2/(r0+dabs(zs-zrec)+r(nr))
        ymax=0.d0
30      yabs=0.d0
        call edgkern(y,klimit,ps,sh,eps)
        ncall=ncall+1
        do i=1,5,2
          yabs=yabs+(y(i)-y0(i))**2
        enddo
	yabs=klimit*dsqrt(klimit*yabs)*dexp(-(klimit*r00(1))**2)
	ymax=dmax1(ymax,yabs)
        if(yabs.gt.eps0*ymax)then
          klimit=1.2d0*klimit
          goto 30
        endif
%
%	determine wavenumber sampling rate
%
	dk=pi2/(srate*(r00(1)+r(nr))+dabs(zs-zrec))
	nk=500+idnint(klimit/dk)
	dk=klimit/dble(nk)
%
%	too small distances will be treated as r = 0!
%
	if(r(1).gt.0.d0)then
	  ir1=1
	else
	  ir1=2
	endif
%
        do ik=1,nk
          k=dble(ik)*dk
          call edgkern(y,k,ps,sh,eps)
	  if(analytic)then
	    do i=1,5,2
	      y(i)=y(i)-y0(i)
	    enddo
	    do i=2,6,2
	      y(i)=y(i)-y0(i)*k
	    enddo
	  else if(ir1.eq.2)then
c
c	  for r=0
c
	    fac=k*dexp(-(k*r00(1))**2)*dk
	    do i=1,6
	      u0(i)=u0(i)+y(i)*fac
	      uk0(i)=uk0(i)+y(i)*k*fac
	    enddo
	  endif
          do ir=ir1,nr
	    fac=dsqrt(k)*dk*dexp(-(k*r00(ir))**2)/dsqrt(r(ir))
	    x=k*r(ir)
c
c	    bessels functions from pre-calculated tables
c
	    nx=idint(x/bsdx)
	    wr=x/bsdx-dble(nx)
	    wl=1.d0-wr
	    if(nx.gt.nnbess)then
	      nx=nnbess+mod(nx-nnbess,ndbess)
	      do i=1,3
	        bs(i)=fac*(wl*bsfct(nx,i)+wr*bsfct(nx+1,i))
	      enddo
	      bs(3)=bs(3)*(dble(nx)+wr)*bsdx/x
	    else
	      do i=1,3
	        bs(i)=fac*(wl*bsfct(nx,i)+wr*bsfct(nx+1,i))
	      enddo
	    endif
%
%	    u1-3 are displacement components:
%	    u4 = normal stress: szz
%	    u5 = surface strain: err+ett
%	    u6 will be derived later
%	    u7 = shear stress: szr
%	    u8 = strain component: dut/dr - (dur/dt)/r + ut/r
%	    u9 = shear stress: szt
%	    u10 = tilt: duz/dr
%
	    u(1,ir)=u(1,ir)+y(1)*bs(1)
	    u(2,ir)=u(2,ir)+y(3)*bs(2)+cs*y(5)*bs(3)
	    u(3,ir)=u(3,ir)-cs*y(3)*bs(3)-y(5)*bs(2)
	    u(4,ir)=u(4,ir)+y(2)*bs(1)
	    u(5,ir)=u(5,ir)-y(3)*k*bs(1)
	    u(7,ir)=u(7,ir)+y(4)*bs(2)+cs*y(6)*bs(3)
	    u(8,ir)=u(8,ir)+y(5)*k*bs(1)
	    u(9,ir)=u(9,ir)-cs*y(4)*bs(3)-y(6)*bs(2)
	    u(10,ir)=u(10,ir)+y(1)*k*bs(2)
	  enddo
        enddo
c
c       end of total integral
c
        if(itty.eq.1)then
          write(*,'(a,i7,a,i7)')'   wavenumber samples: ',ncall+nk,
     &                          ', really used: ',nk
        endif
c
	if(ir1.eq.2.and..not.analytic)then
c
c	  for very small r including r=0
c
	  if(ms.eq.0)then
	    u(1,1)=u0(1)
	    u(4,1)=u0(2)
	    u(5,1)=-0.5d0*uk0(3)
	    u(6,1)=u(5,1)
	  else if(ms.eq.1)then
	    u(2,1)=0.5d0*(u0(3)+cs*u0(5))
	    u(3,1)=-0.5d0*(cs*u0(3)+u0(5))
	    u(7,1)=0.5d0*(u0(4)+cs*u0(6))
	    u(9,1)=-0.5d0*(cs*u0(4)+u0(6))
	    u(10,1)=0.5d0*uk0(1)
	  else if(ms.eq.2)then
	    u(5,1)=0.25d0*(uk0(3)+cs*uk0(5))
	    u(6,1)=-u(5,1)
	    u(8,1)=-0.25d0*(cs*uk0(3)+uk0(5))
	  endif
	endif
	do ir=ir1,nr
	  if(analytic)then
	    if(ms.eq.0)then
	      bs(1)=0.d0
	      bs(2)=-1.d0/r(ir)**2
	      bs(3)=0.d0
	      bs1(1)=-1.d0/r(ir)**3
	      bs1(2)=0.d0
	      bs1(3)=0.d0
	    else if(ms.eq.1)then
	      bs(1)=1.d0/r(ir)**2
	      bs(2)=-1.d0/r(ir)**2
	      bs(3)=1.d0/r(ir)**2
	      bs1(1)=0.d0
	      bs1(2)=-2.d0/r(ir)**3
	      bs1(3)=1.d0/r(ir)**3
	    else if(ms.eq.2)then
	      bs(1)=2.d0/r(ir)**2
	      bs(2)=-1.d0/r(ir)**2
	      bs(3)=2.d0/r(ir)**2
	      bs1(1)=3.d0/r(ir)**3
	      bs1(2)=-4.d0/r(ir)**3
	      bs1(3)=4.d0/r(ir)**3
	    endif
	    u(1,ir)=u(1,ir)+y0(1)*bs(1)
	    u(2,ir)=u(2,ir)+y0(3)*bs(2)+cs*y0(5)*bs(3)
	    u(3,ir)=u(3,ir)-y0(5)*bs(2)-cs*y0(3)*bs(3)
	    u(4,ir)=u(4,ir)+y0(2)*bs1(1)
	    u(5,ir)=u(5,ir)-y0(3)*bs1(1)
	    u(7,ir)=u(7,ir)+y0(4)*bs1(2)+cs*y0(6)*bs1(3)
	    u(8,ir)=u(8,ir)+y0(5)*bs1(1)
	    u(9,ir)=u(9,ir)-y0(6)*bs1(2)-cs*y0(4)*bs1(3)
	    u(10,ir)=u(10,ir)+y0(1)*bs1(2)
	  endif
%
%	  u6 is ett = ur/r + (dut/dt)/r
%
	  u(6,ir)=(u(2,ir)+cs*dble(ms)*u(3,ir))/r(ir)
%
%	  u5 now is err = u5(before) - ett
%
	  u(5,ir)=u(5,ir)-u(6,ir)
%
%	  u8 now is ert = 0.5 * u8(before) + (dur/dt)/r - ut/r
%	                = 0.5 * (dut/dr + (dur/dt)/r - ut/r)
%
	  u(8,ir)=0.5d0*u(8,ir)-(cs*dble(ms)*u(2,ir)+u(3,ir))/r(ir)
	enddo
do ir=1,nr
	  u(4,ir)=(u(4,ir)-lambda*(u(5,ir)+u(6,ir)))/(lambda+2.d0*mu)
	  u(7,ir)=u(7,ir)/(2.d0*mu)
	  u(9,ir)=u(9,ir)/(2.d0*mu)
end
