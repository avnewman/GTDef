function [ yy ] = GTdef_edgpsv(edgrn,srcz,recz,kk,eps,eps0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    GTdef_edgpsv.m                                     %
%                                                                                       %
% Calculate responses to P-SV source based on EDGRN edgpsv.F                            %
% EDGRN uses cylindrical coordinate system, depth positive down                         %
%                                                                                       %
% INPUT										        %
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
% (2) srcz    - source depth [m]                                                        %
% (3) recz    - receiver depth [m]                                                      %
% (2) kk      - wave number                                                             %
% (3) eps,eps0- relative accuracy  eps = 1.0e-8; eps0 = 1.0e-3                          %
%										        %
%										        %
% OUTPUT                                                                                %
% yy          - solution vector for the poloidal mode [1x4]                             %
%                                                                                       %
% REFERENCE  									        %
% Wang, R., Martin, F. L., & Roth, F. (2003)					        %
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5	        %
%		                                                                        %
% first created by Lujia Feng Wed Dec 10 17:57:53 SGT 2014                              %
% last modified by Lujia Feng 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g0 = 0.0; % or 0.982d+03

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
rr    = edgrn.rr;

zz    = edgrn.zz; 
z1    = edgrn.z1; 
z2    = edgrn.z2; 
ro    = edgrn.ro;
vp    = edgrn.vp;
vs    = edgrn.vs;

zp    = ; % top depth of propagation layers
hp ; % thickness of propagation layers
lp    =
nno   =
srcInd = % ~ ls
recInd = % ~ lzrec

yy  = zeros(1,4);
c0  = zeros(4,2);
yup = zeros(4,2);
y0  = zeros(4,2);

%-------------------- matrix propagation from surface to source --------------------
% y0(4,2,1): 2 starting solution vectors

% find the uppermost layer (upInd) which has significant deformation
upInd    = 1; % start from the surface layer
exponent = 0.0;
for ii = srcInd-1:-1:1
   exponent = exponent - kk*hp(ii);
   if (exponent<=log(eps))
      upInd = ii+1;
      if(upInd>recInd) return, end; % wave vanishes before reaching the receiver -> yy = [0 0 0 0]
   end
end

% find the lowermost layer (lwInd) which has significant deformation
lwInd    = lp; % start from the halfspace layer
exponent = 0.0;
for ii = srcInd:lp-1
   exponent = exponent - kk*hp(ii);
   if (exponent<log(eps))
      lwInd = ii;
      if(lwInd<recInd) return, end; % wave vanishes before reaching the receiver -> yy = [0 0 0 0]
   end
end

if upInd==1 % if the uppermost layer is free surface
   yup(1,1) = 1.0;
   dro      = ro(1);
   yup(2,1) = dro*g0*yup(1,1);
   yup(3,2) = 1.0;
else
   c0(1,1)  = 1.0;
   c0(3,2)  = 1.0;
   nn       = nno(upInd-1);
   AA       = GTdef_edgmatrix(kk,0.0,vp(nn),vs(nn),ro(nn));
   yup      = AA*c0;
   nn1      = nno(upInd);
   dro      = ro(nn1)-ro(nn);
   yup(2,:) = yup(2,:) + dro*g0*yup(1,:);
end

if(upInd==recInd) y0 = yup, end; % if the uppermost layer is receiver

for ii = upInd+1:srcInd
   h0 = hp(ii-1);
   nn = nno(ii-1);
   if (kk*h0<=eps0)
      [ hk ] = GTdef_edghask(kk,h0,vp(nn),vs(nn),ro(nn),eps0,4);
      y1  = hk*yup;
      yup = y1;
   else
      % determine propagation matrix
      AAi   = GTdef_edgmatinv(kk,0.0,vp(nn),vs(nn),ro(nn));
      AA    = GTdef_edgmatrix(kk,h0,vp(nn),vs(nn),ro(nn));
      c0    = AAi*yup;
      wave  = exp(-kk*h0);
      % normalization of all modes
      norm1 = sum(c0(:,1).*c0(:,1));
      norm2 = sum(c0(:,2).*c0(:,2)); 
      fac   = 1.0/sqrt(norm1*norm2);
      % orthogonalization of the p-sv modes
      orth(1,1) =  c0(3,2)*fac;
      orth(1,2) = -c0(1,2)*fac;
      orth(2,1) = -c0(3,1)*fac;
      orth(2,2) =  c0(1,1)*fac;
      c1        = c0*orth;
      if(ii>recInd)
         orth = orth*wave;
	 y1   = y0*orth;
	 y0   = y1;
      end
%     c1(1,1)=c1(1,1)
      c1(2,1) = c1(2,1)*wave*wave;
      c1(3,1) = (0.0,0.0);
      c1(4,1) = c1(4,1)*wave*wave;

      c1(1,2) = (0.0,0.0);
      c1(2,2) = c1(2,2)*wave*wave;
%     c1(3,2) = c1(3,2)
      c1(4,2) = c1(4,2)*wave*wave;
      yup     = AA*c1;
   end
   dro = ro(nno(ii))-ro(nno(ii-1));
   yup(2,:) = yup(2,:) + dro*g0*yup(1,:);
   if(ii==recInd) y0 = yup, end;
end


%-------------------- matrix propagation from half-space to source --------------------
        do i=1,4
          do j=1,2
            c0(i,j)=0.d0
            ylw(i,j)=0.d0
          enddo
        enddo
c
c
c       c0(4,2): 2 coefficient vectors in the half-space
c
	n=nno(lwInd)
        if(vs(n).lt.eps0*vp(n))then
c
c         the lowest layer is fluid
c
	  ylw(1,1)=1.d0
          ylw(3,2)=1.d0
        else
c
c         the lowest layer is solid
c
          c0(2,1)=1.d0
          c0(4,2)=1.d0
          call edgmatrix(ma,4,k,0.d0,n)
          call axb(ma,c0,4,4,2,ylw)
        endif
        if(lwInd.gt.ls.and.lwInd.eq.lzrec)call memcpy(ylw,y0,8)
c
        do l=lwInd-1,ls,-1
          h0=hp(l)
          n=nno(l)
	  dro=ro(nno(l+1))-ro(nno(l))
	  do j=1,2
	    ylw(2,j)=ylw(2,j)-dro*g0*ylw(1,j)
	  enddo
	  if(k*h0.le.eps0)then
	    call edghask(hk,4,k,-h0,n)
	    call axb(hk,ylw,4,4,2,y1)
	    call memcpy(y1,ylw,8)
	  else
c
c
c           determination of propagation matrix
c
            call edgmatinv(mai,4,k,0.d0,n)
            call edgmatrix(ma,4,k,-h0,n)
            call axb(mai,ylw,4,4,2,c0)
	    wave=dexp(-k*h0)
c
c           normalization of all modes
c
	    do j=1,2
              norm(j)=0.d0
              do i=1,4
                norm(j)=norm(j)+c0(i,j)*c0(i,j)
              enddo
	    enddo
	    fac=1.d0/dsqrt(norm(1)*norm(2))
c
c           orthogonalization of the p-sv modes
c
            orth(1,1)=c0(4,2)*fac
            orth(1,2)=-c0(2,2)*fac
            orth(2,1)=-c0(4,1)*fac
            orth(2,2)=c0(2,1)*fac
            call axb(c0,orth,4,2,2,c1)
            if(l.lt.lzrec)then
	      do j=1,2
	        do i=1,2
	          orth(i,j)=orth(i,j)*wave
	        enddo
	      enddo
              call axb(y0,orth,4,2,2,y1)
              call memcpy(y1,y0,8)
            endif
c
	    c1(1,1)=c1(1,1)*wave*wave
c	    c1(2,1)=c1(2,1)
	    c1(3,1)=c1(3,1)*wave*wave
            c1(4,1)=(0.d0,0.d0)
c
            c1(1,2)=c1(1,2)*wave*wave
            c1(2,2)=(0.d0,0.d0)
            c1(3,2)=c1(3,2)*wave*wave
c           c1(4,2)=c1(4,2)
c
            call axb(ma,c1,4,4,2,ylw)
	  endif
          if(l.gt.ls.and.l.eq.lzrec)call memcpy(ylw,y0,8)
        enddo
c
%-------------------- conditions on the source surface --------------------
        do i=1,4
          if(kpower(i).eq.1)then
            fac=k
          else
            fac=1.d0
          endif
          b(i)=sfct(i)*fac
          do j=1,2
            coef(i,j)=yup(i,j)
            coef(i,j+2)=-ylw(i,j)
          enddo
        enddo
        key=0
	call gemp(coef,b,4,1,1.d-99,key)
        if(key.eq.0)then
          print *,'warning in edgpsv: anormal exit from cgemp!'
          return
        endif
        if(lzrec.le.ls)then
          do i=1,4
            y(i)=0.d0
            do j=1,2
              y(i)=y(i)+b(j)*y0(i,j)
            enddo
          enddo
        else
          do i=1,4
            y(i)=0.d0
            do j=1,2
              y(i)=y(i)+b(j+2)*y0(i,j)
            enddo
          enddo
        endif
c
        return
        end
