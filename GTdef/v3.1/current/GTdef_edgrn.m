function [ earth ] = GTdef_edgrn(earth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   GTdef_edgrn.m					%
%											%
% Create point source greens functions                                                  %
% This code is converted from fortran EDGRN code                                        %
%											%
% INPUT											%
% (1) earth.edgrn structure                                                             %
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
% edgrn.r0    - point source scale [m]                                        (scale)   %  
%											%
% (2) earth.layer = [ id depth vp vs ro ]	(nn*5)                                  %
%											%
%---------------------------------------------------------------------------------------%
% INTERMEDIATE                                                                          %
% (1) earth.sublayer structure (model sublayers)                                        %
%     corresponds to /model/ h,ro,vp,vs,n0 [h->hh,n0->nl]                               %
% sublayer.nl    - number of sublayers                                       (scalar)   %
% sublayer.topz  - top depth of each sublayer [m]                (nl*1 column vector)   %
% sublayer.botz  - bottom depth of each sublayer [m]             (nl*1 column vector)   %
% sublayer.hh    - thickness of each sublayer [m]                (nl*1 column vector)   %
% sublayer.vp    - P-wave velocity for each sublayer [m/s]       (nl*1 column vector)   %
% sublayer.vs    - S-wave velocity for each sublayer [m/s]       (nl*1 column vector)   %
% sublayer.ro    - density for each sublayer [kg/m^3]            (nl*1 column vector)   %
%											%
% (2) depths structure (changing layers including source and receiver)                  %
%     corresponds to /sublayer/ hp,lp,nno  [hp->hh,lp->nl]                              %
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
% (3) source structure (changing source)                                                %
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
% (4) bess structure                                                                    %
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
% OUTPUT (rowwise to be consistent with Okada)                                          %
% (1) earth.edgrnfcts structure								%
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
% last modified by Lujia Feng Thu Feb  4 15:19:44 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folderName = 'edgrnfcts';
if ~exist(folderName,'dir')
   mkdir(folderName);
end
fssName = [ folderName '/izmhs.ss' ];
fdsName = [ folderName '/izmhs.ds' ];
fclName = [ folderName '/izmhs.cl' ];

% read in edgrn structure locally
obsz  = earth.edgrn.obsz;
nr    = earth.edgrn.nr;
minr  = earth.edgrn.minr;
maxr  = earth.edgrn.maxr;
nz    = earth.edgrn.nz;
minz  = earth.edgrn.minz;
maxz  = earth.edgrn.maxz;
srate = earth.edgrn.srate;

if nr<=1 || nz<=1
   error('GTdef_edgrn ERROR: nz & nr must be > 1!');
end
if minr >= maxr
   error('GTdef_edgrn ERROR: minr must be < maxr!');
end
if minz >= maxz
   error('GTdef_edgrn ERROR: minz must be < maxz!');
end
if srate<=10
   srate = 10;
end

% calculate other edgrn parameters
dr   = (maxr-minr)/(nr-1);
dz   = (maxz-minz)/(nz-1);
rr   = (minr:dr:maxr)';
zz   = (minz:dz:maxz)';
r0   = 0.5*max(dz,dr);

% check insignificant depth difference between source and observation depth
zdiff = abs(zz-obsz);
ind   = find(zdiff<1e-2*dz);
if isscalar(ind)
   obsz = zz(ind);
elseif length(ind)>=2
   error('GTdef_edgrn ERROR: source depth wrong!');
end

earth.edgrn.obsz = obsz;
earth.edgrn.dr   = dr;
earth.edgrn.dz   = dz;
earth.edgrn.rr   = rr;
earth.edgrn.zz   = zz;
earth.edgrn.r0   = r0;

% divide earth.layer into earth.sublayers
[ sublayer ] = GTdef_edgsublay(earth.layer);

nl   = sublayer.nl;
topz = sublayer.topz;
botz = sublayer.botz;
hh   = sublayer.hh;
vp   = sublayer.vp;
vs   = sublayer.vs;
ro   = sublayer.ro;

% calcuate properties at receiver depth
% call edglayer(zs1,ls,zrec,lzrec,h,n0)
if obsz == 0
   ind = 1;
else
   ind = find(obsz>topz & obsz<=botz);
end
mu     = ro(ind)*vs(ind)*vs(ind);        % shear modulus/rigidity/Lamé's second parameter
lambda = ro(ind)*vp(ind)*vp(ind)-2.0*mu; % Lamé's first parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nr->nr+2 nz->nz+2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extend another two records
ssdisp0  = zeros(3,nr+2,nz+2);	% to avoid conflict naming with Matlab
ssstrn0  = zeros(6,nr+2,nz+2);
ssuzr0   = zeros(nr+2,nz+2);
dsdisp0  = zeros(3,nr+2,nz+2);
dsstrn0  = zeros(6,nr+2,nz+2);
dsuzr0   = zeros(nr+2,nz+2);
cldisp0  = zeros(2,nr+2,nz+2);
clstrn0  = zeros(4,nr+2,nz+2);
cluzr0   = zeros(nr+2,nz+2);

% 1='Source type: strike-slip'
% 2='Source type: dip-slip'
% 3='Source type: compensated linear vector dipole'
for istype=1:3
   switch istype
      case 1
         fprintf(1,'Source type: strike-slip\n');
         fout = fopen(fssName,'w');
         fprintf(fout,'# Green functions calculated with GTdef_edgrn.m\n');
         fprintf(fout,'# Source type: strike-slip;  Dislocation: 1 [m]\n');
         fprintf(fout,'# The first line: nr r1[m] r2[m] nzs zs1[m] zs2[m] obs.depth[m] obs.lambda[Pa] obs.mu[Pa]\n');
         fprintf(fout,'# The following lines are\n');
         fprintf(fout,'#   uz[m]         ur[m]         ut[m]         ezz           err           ett           ezr           ert           etz           duz/dr\n');
         srctype = 'ss';
      case 2
         fprintf(1,'Source type: dip-slip\n');
         fout = fopen(fdsName,'w');
         fprintf(fout,'# Green functions calculated with GTdef_edgrn.m\n');
         fprintf(fout,'# Source type: dip-slip;  Dislocation: 1 [m]\n');
         fprintf(fout,'# The first line: nr r1[m] r2[m] nzs zs1[m] zs2[m] obs.depth[m] obs.lambda[Pa] obs.mu[Pa]\n');
         fprintf(fout,'# The following lines are\n');
         fprintf(fout,'#   uz[m]         ur[m]         ut[m]         ezz           err           ett           ezr           ert           etz           duz/dr\n');
         srctype = 'ds';
      case 3
         fprintf(1,'Source type: compensated linear vector dipole\n');
         fout = fopen(fclName,'w');
         fprintf(fout,'# Green functions calculated with GTdef_edgrn.m\n');
         fprintf(fout,'# Source type: compensated linear vector dipole;  Dislocation: 1 [m]\n');
         fprintf(fout,'# The first line: nr r1[m] r2[m] nzs zs1[m] zs2[m] obs.depth[m] obs.lambda[Pa] obs.mu[Pa]\n');
         fprintf(fout,'# The following lines are\n');
         fprintf(fout,'#   uz[cm]        ur[cm]        ezz           err           ett           ezr           duz/dr\n');
         srctype = 'cl';
   end
   fprintf(fout,'%.0f %e %e %.0f %e %e %e %e %e\n',nr,minr,maxr,nz,minz,maxz,obsz,lambda,mu);

   % calculate Bessel function tables
   [ bess ] = GTdef_edgbstab(srctype);

   % loop through source at different depths
   for ii=1:nz
      srcz = zz(ii); % source depth
      [ depths ] = GTdef_edgdepth(topz,botz,srcz,obsz);
      ind    = depths.nno(depths.lsrc);
      srcvp  = vp(ind);
      srcvs  = vs(ind); 
      srcro  = ro(ind);
      [ source ] = GTdef_edgmoment(srctype,srcvp,srcvs,srcro);
      source.r0  = r0;
      [ uu ] = GTdef_edgwvint(rr,nr,srate,lambda,mu,sublayer,depths,source,bess);

      % write out
     if istype<=2
        for jj=1:nr
                       %     1      2      3      4      5      6      7      8      9     10
           fprintf(fout,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n',uu(:,jj));
        end
     else
        for jj=1:nr
                       %     1      2      3      4      5      6      7 
           fprintf(fout,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n',uu([1 2 4 5 6 7 10],jj));
        end
     end

     % add to greens function
     if istype==1
        ssdisp0(:,1:nr,ii) = uu(1:3,:);
        ssstrn0(:,1:nr,ii) = uu(4:9,:);
        ssuzr0(1:nr,ii)    = uu(10,:);
     end
     if istype==2
        dsdisp0(:,1:nr,ii) = uu(1:3,:);
        dsstrn0(:,1:nr,ii) = uu(4:9,:);
        dsuzr0(1:nr,ii)    = uu(10,:);
     end
     if istype==3
        cldisp0(:,1:nr,ii) = uu(1:2,:);
        clstrn0(:,1:nr,ii) = uu(3:6,:);
        cluzr0(1:nr,ii)    = uu(7,:);
     end
   end

   fclose(fout);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if observation point coincides with source point, output is zero
% so use the nearest to replace the zero output
if  obsz>=minz && obsz<=maxz && minr<1.0e-6
   iz  = round((obsz-minz)/dz);
   dzs = abs(obsz-(minz+dz*iz))/dz;
   if dzs<0.02
      iz = iz+1;
      ssdisp0(:,1,iz) = ssdisp0(:,2,iz);
      dsdisp0(:,1,iz) = dsdisp0(:,2,iz);
      cldisp0(:,1,iz) = cldisp0(:,2,iz);
      ssstrn0(:,1,iz) = ssstrn0(:,2,iz);
      dsstrn0(:,1,iz) = dsstrn0(:,2,iz);
      clstrn0(:,1,iz) = clstrn0(:,2,iz);
      ssuzr0(1,iz)    = ssuzr0(2,iz);
      dsuzr0(1,iz)    = dsuzr0(2,iz);
      cluzr0(1,iz)    = cluzr0(2,iz);
   end
end

% assign to structure
earth.sublayer         = sublayer;
earth.edgrnfcts.ssdisp = ssdisp0;
earth.edgrnfcts.dsdisp = dsdisp0;
earth.edgrnfcts.cldisp = cldisp0;
earth.edgrnfcts.ssstrn = ssstrn0;
earth.edgrnfcts.dsstrn = dsstrn0;
earth.edgrnfcts.clstrn = clstrn0;
earth.edgrnfcts.ssuzr  = ssuzr0;
earth.edgrnfcts.dsuzr  = dsuzr0;
earth.edgrnfcts.cluzr  = cluzr0;
