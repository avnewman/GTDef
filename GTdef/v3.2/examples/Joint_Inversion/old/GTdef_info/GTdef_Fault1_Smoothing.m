
function []=GTdef_Fault1_Smoothing(matfile_in,kappa)
%code runs an inversion through GTdef with a trick for non-regularly
%gridded data.  While this code is broad enough to run in various
%scenarios, it was designed with the 2015 Illapel earthquake in mind
%(published as Willamson et al 2017, JGR). In that study, I used 25km x 25
%km subfaults that varied with strike and dip along the study region. The
%variable strike posed a major problem with GTdef's smoothing algorithm.

% steps before this, run GTdef for a kappa value of zero with the matfile
% flag on.  This saves the matfile with all relevant data. This original
% input treats each subfault as its own master fault.


load(matfile_in)

% what is in matfile_in? Well it is a kappa=0 inversion using all the
% geodetic data that one may want to use in the inversion. The reasoning
% for this parameter is that GTdef will not smooth non-regularly gridded
% subfaults. In the case for the Illapel earthquake subfaults, there was a
% variable strike (although for most cases this was quite small). So I need
% to trick GTdef by just inputting a kappa=0 output and then smoothing
% that, with the code thinking it is a regular grid.

%kappa-  I need a vector of kappa values to explore. ex: kappa=linspace(0,60000,25); % set up kappa range

modspace.modinfo =[];    % I will modify some of the model info so I need to open the modspace file
modspace.kappa=kappa;    % which is part of the matfile_in
modspace.beta=kappa.^2;

% variables to declare:
dd=25000; ds=25000; Nd=8; Ns=23;     % values needed to smooth main faults, lenth down dip, length along strike, number down dip, number along strike
figure
for zz=1:length(kappa)
    % grab all data from .mat file produced by GTdef
beta = kappa(zz)*kappa(zz);

Lgrn  = modspace.Lgrn;   Xgrn  = modspace.Xgrn;   % LOS and GPS data if used
lb    = modspace.lb;     ub    = modspace.ub;     % lower and upper bounds of slip
x0    = modspace.x0;     Aeq   = modspace.Aeq;    beq   = modspace.beq;
Aineq = modspace.Aineq;  bineq = modspace.bineq;  % various inverse parameters
C=modspace.C;  % Green's function matrix
d=modspace.d;  % data vector

%% code is modified from GTdef_fault1dif.m
% this will generate my two-dimensional smoothing matrix (sm)
[ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(dd,ds,Nd,Ns);   % currently only care about sm_2d
sm=sm_2d;
ind_fixed = find(lb==-Inf);     % index for fixed slips
sm(ind_fixed) = 0;              % don't do smoothing for them
modspace.sm=sm;
sm = GTdef_condense(modspace.sm);

%% code is modified from GTdef_invert.m

sm = sm.*beta;                 % smoothing matrix is modified by beta which is a function of kappa
d_sm = zeros(size(sm,1),1);    % append zeros to data vector
C = [ C;sm ]; d = [ d;d_sm ];  % and add smothing matrix to GF matrix and zeros to the data matrix

options = optimset('MaxIter',2000,'TolFun',1e-30);   
[xx,resnorm] = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0,options);   % same inversion scheme as GTdef
modspace.xx=xx;

%T=xx(184:(184*2)-1);
T=zeros(Nd*Ns,1);    % I am only interested in the thrust component, rescaling output 
for n = 1:(Nd*Ns)
    a(n) = 2 + 3*(n - 1) ;
    T(n)=xx(a(n));
end

% quick plot of different kappas
TT=reshape(T,[Nd,Ns]);   % plotting each Kappa for a quick visual check
subplot(6,6,zz)
image(TT,'CDataMapping','scaled'); title( ['Kp=' num2str(kappa(zz))]); caxis([0,15]); colorbar;


%% quantify residuals, roughness, etc

foutName = strcat(basename,'_kp',num2str(kappa(zz),'%-.0f'),'.out');
[ modspace,pnt,los,bsl,nod ] = GTdef_forward(modspace,pnt,los,bsl,nod);

[ flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt ] = GTdef_update_slips(earth,modspace,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt);

GTdef_output(foutName,earth,modspace,bt,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt,addon,pnt,los,bsl,prf,grd,nod);

%% build files needed for GMT.  currently the multisegment fault slip files
Parse_GTdef_out(foutName)
end
[ ~,basename,~ ] = fileparts(finName);
fsumName = [ basename '_man_inv.out' ];
GTdef_summary(fsumName,modspace);

%% quick look at tradeoff curve
 figure
plot(modspace.modinfo(:,11),modspace.modinfo(:,5),'o-'); ylabel('RMS [m]'); xlabel('roughness');
end
