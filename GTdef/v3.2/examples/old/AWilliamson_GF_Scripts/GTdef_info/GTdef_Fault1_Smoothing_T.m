
function []=GTdef_Fault1_Smoothing_T(matfile_in,kappa,run_type)



%code runs an inversion through GTdef with a trick for non-regularly
%gridded data.  While this code is broad enough to run in various
%scenarios, it was designed with the 2015 Illapel earthquake in mind
%(published as Willamson et al 2017, JGR). In that study, I used 25km x 25
%km subfaults that varied with strike and dip along the study region. The
%variable strike posed a major problem with GTdef's smoothing algorithm.

% steps before this, run GTdef for a kappa value of zero with the matfile
% flag on.  This saves the matfile with all relevant data. This original
% input treats each subfault as its own master fault.

% RUN_TYPE FLAG:
% if run_type== 1, then only tsunami data will be used, still
%                  need a geodetic file to get the fault geometry, but GF 
%                  and data will be replaced.
% if run_type==2,  then only geodetic data will be used.
% if run_type==3,  then a joint inversion 
%
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
modspace.Tgrn=[]; Tgrn=modspace.Tgrn;
%% build Tgrn dataset
if run_type ~= 2
    
if exist('TsuGF.mat', 'file')
load TsuGF.mat % tsunami waveform Green's function matrix. derived from python
Tgrn=zeros(size(GF,1),size(GF,2)*3);   % set up the empty Tgrn matrix
                                      % has variables GF and data
tsu.obs=data';    
tsu.wgt=ones(length(data),1)*.75;
tsu.obs_err=ones(length(data),1)*(0.001603);  % adding error to the tsunami waveform
                                              % based on range in
                                              % measurements a few hours
                                              % before the earthquake
tsu.coef=sqrt(tsu.wgt)./tsu.obs_err;
    for n = 1:size(GF,2)  % number of subfaults
         a(n) = 2 + 3*(n - 1) ;
         Tgrn(:,a(n))=GF(:,n);
    end
 % formatting tsunami GF data
 CC = []; dd = []; mod_Tgrn=Tgrn; mod_tsu.obs=tsu.obs;
if ~isempty(Tgrn)
    eqNum = size(Tgrn,1); % equation number
    for ii = 1:eqNum
        mod_Tgrn(ii,:)  = Tgrn(ii,:).*tsu.coef(ii);
        mod_tsu.obs(ii) = tsu.obs(ii).*tsu.coef(ii);
    end
    ind = find(~isnan(tsu.obs));                % exclude nan values
    CC = [ CC;mod_Tgrn(ind,:) ]; dd = [ dd;mod_tsu.obs(ind) ];
end
C_tsu=CC; d_tsu=dd;
modspace.Tgrn=Tgrn;

C=(modspace.C);
d=(modspace.d);

d=[d;  d_tsu]; modspace.d=d;
C=[C;  C_tsu];  modspace.C=C;
end
end
%%
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

if run_type == 1
C=C_tsu;
d=d_tsu;
Lgrn=[];
Xgrn=[];
end
    
%% code is modified from GTdef_fault1dif.m
% this will generate my two-dimensional smoothing matrix (sm)
[ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(dd,ds,Nd,Ns);   % currently only care about sm_2d
sm=sm_2d;
ind_fixed = find(lb==-Inf);     % index for fixed slips
sm(ind_fixed) = 0;              % don't do smoothing for them
modspace.sm=sm;
sm = GTdef_condense(modspace.sm);
sm2=sm;  %%%%%%
modspace.sm   = sm;
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

%%    Quantify the residuals

data_num = 0;
rss = 0; wrss = 0; chi2 = 0;
rms = 0; wrms = 0; rchi2 = 0;
%%%%%

sm=sm2; %%%%%
% find slips that are included in the model
total = length(xx);
num = 0;
for ii = 1:total
    if xx(ii)==0&&isinf(lb(ii))&&isinf(ub(ii))                  % number of slips that are not included
         num = num+1;
    end
end
slip_num = total-num;

%%%%%
if ~isempty(Xgrn)
    Xmod = Xgrn*xx;                                             % model prediction for point data
    pnt_mod = reshape(Xmod,[],3);
    pnt_mod_err = nan(size(pnt_mod));
    pnt.out = [ pnt.loc pnt_mod pnt_mod_err pnt.wgt ];          % reuse pnt.loc and pnt.wgt

    data_ind = find(~isnan(pnt.obs));                           % exclude nan data
    data_num = data_num+length(data_ind);
    pnt_dif2 = (Xmod(data_ind)-pnt.obs(data_ind)).^2;           % squared residuals/differences
    pnt_err2 = pnt.obs_err(data_ind).^2;                        % squared errors
    rss = rss+sum(pnt_dif2);                                    % residual sum of squares
    wrss = wrss+sum(pnt.obs_wgt(data_ind).*pnt_dif2./pnt_err2); % weighted rss
    chi2 = chi2+sum(pnt_dif2./pnt_err2);                        % chi-square
end

if ~isempty(Lgrn)  % if LOS GF matrix is populated
    los_mod = Lgrn*xx;                                          % model prediction for los
    los_mod_err = nan(size(los_mod));
    los.out = [ los.loc los_mod los_mod_err los.dir los.wgt ];  % reuse los.loc, los.dir and los.wgt

    data_ind = find(~isnan(los.obs));                           % exclude nan data
    data_num = data_num+length(data_ind);
    los_dif2 = (los_mod(data_ind)-los.obs(data_ind)).^2;        % squared residuals/differences
    los_err2 = los.obs_err(data_ind).^2;                        % squared errors
    rss  = rss+sum(los_dif2);                                   % residual sum of squares
    wrss = wrss+sum(los.obs_wgt(data_ind).*los_dif2./los_err2); % weighted rss
    chi2 = chi2+sum(los_dif2./los_err2);                        % chi-square
end
if ~isempty(Tgrn)

    tsu_mod=Tgrn*xx;

    tsu_mod_err=nan(size(tsu_mod));
    tsu.out=[tsu_mod tsu_mod_err tsu.obs tsu.wgt];
    data_ind = find(~isnan(tsu.obs));
    data_num = data_num+length(data_ind);
    tsu_dif2 = (tsu_mod(data_ind)-tsu.obs(data_ind)).^2;
    tsu_err2 = tsu.obs_err(data_ind).^2;
    rss  = rss+sum(tsu_dif2);
    wrss = wrss+sum(tsu.wgt(data_ind).*tsu_dif2./tsu_err2); % weighted rss
    chi2 = chi2+sum(tsu_dif2./tsu_err2);
end

if ~isempty(sm) && size(sm,2) == size(xx,1)
    sm_all = sum(abs(sm*xx));

    sm_num = size(sm,1);

        r_1d = nan;
        r_2d = 1e8*0.5*sm_all/sm_num;                           % NOTE: roughness [cm/km^2]; 0.5 accounts for both ds and ss

else %
    r_1d = nan; r_2d = nan;
end


ndf = data_num-slip_num;                                        % nominal number of degrees of freedom
rms = sqrt(rss/data_num);                                       % root mean square
wrms = sqrt(wrss/data_num);
rchi2 = chi2/ndf;                                               % reduced chi-square
strain = NaN;
modinfo          = [ data_num slip_num ndf rss rms wrss wrms chi2 rchi2 r_1d r_2d strain];
modspace.modinfo = [ modspace.modinfo; modinfo ];

modspace.kappa= [ modspace.kappa, kappa(zz)];
modspace.beta=[modspace.beta, beta];



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
