clc; clear all; close all;

% steps before this, run GTdef for a kappa value of zero with the matfile
% flag on.  This saves the matfile with all relevant data. 
load TsuGF.mat    % tsunami waveform Green's function matrix. derived from python
Tgrn=zeros(size(GF,1),size(GF,2)*3);   % set up the empty Tgrn matrix
                                       % has variables GF and data
tsu.obs=data';    
tsu.wgt=ones(length(data),1);
tsu.obs_err=ones(length(data),1)*(0.001603);  % adding error to the tsunami waveform
                                              % based on range in
                                              % measurements a few hours
                                              % before the earthquake
tsu.coef=sqrt(tsu.wgt)./tsu.obs_err;
for n = 1:size(GF,2)  % number of subfaults
    a(n) = 2 + 3*(n - 1) ;
    Tgrn(:,a(n))=GF(:,n);
end

C = []; d = []; mod_Tgrn=Tgrn; mod_tsu.obs=tsu.obs;
if ~isempty(Tgrn)
    eqNum = size(Tgrn,1); % equation number
    for ii = 1:eqNum
        mod_Tgrn(ii,:)  = Tgrn(ii,:).*tsu.coef(ii);
        mod_tsu.obs(ii) = tsu.obs(ii).*tsu.coef(ii);
    end
    ind = find(~isnan(tsu.obs));                % exclude nan values
    C = [ C;mod_Tgrn(ind,:) ]; d = [ d;mod_tsu.obs(ind) ];
end
C_tsu=C; d_tsu=d;
%%%% NOTE THAT TSUNAMI DOES NOT HAVE WEIGHT OR ERROR-IT IS ALL 1 for WT

load Illapel_los_gps.mat
kappa=linspace(10000,200000,10); % set up kappa range
%kappa=100;  % or just rough version

C_land=modspace.C;
d_land=modspace.d;

modspace.modinfo =[];
modspace.beta=[];
modspace.kappa=[];
modspace.Tgrn=Tgrn;

figure
for zz=1:length(kappa)
beta = kappa(zz)*kappa(zz);

Lgrn  = modspace.Lgrn;  % grab all data from .mat file produced by GTdef
Xgrn  = modspace.Xgrn;
lb    = modspace.lb;
ub    = modspace.ub;
x0    = modspace.x0;
Aeq   = modspace.Aeq;
beq   = modspace.beq;
Aineq = modspace.Aineq;
bineq = modspace.bineq;
C=(modspace.C);
d=(modspace.d);
C = []; d = [];
Lgrn=[];
Xgrn=[];
%d=[  d_tsu];
%d=[ d_land];
d=[  d_tsu;d_land];
C=[  C_tsu;C_land];
%C=[C_land];
%C=[  C_tsu];


%% GTdef_fault1dif.m
dd=25000; ds=25000; Nd=8; Ns=23;     % values needed to smooth main faults
[ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(dd,ds,Nd,Ns);   % currently only care about sm_2d

sm=sm_2d;
ind_fixed = find(lb==-Inf);     % index for fixed slips
sm(ind_fixed) = 0;              % don't do smoothing for them

modspace.sm=sm;

sm     = GTdef_condense(modspace.sm);
sm2=sm;  %%%%%%
sm = sm.*beta;
modspace.sm   = sm;

%% GTdef_invert.m

d_sm = zeros(size(sm,1),1);   % append zeros to data vector
C = [ C;sm ]; d = [ d;d_sm ];  % and add smothing matrix to GF matrix

options = optimset('MaxIter',2000,'TolFun',1e-30);
[xx,resnorm] = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0,options);   % same inversion scheme as GTdef
modspace.xx=xx;
T=zeros(184,1);    % I am only interested in the thrust component, rescaling output 
for n = 1:184
    a(n) = 2 + 3*(n - 1) ;
    T(n)=xx(a(n));
end


%% quick plot of different kappas

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


foutName = strcat(basename,'_kp',num2str(kappa(zz),'%-.0f'),'.out');
%[ modspace,pnt,los,bsl,nod ] = GTdef_forward(modspace,pnt,los,bsl,nod);

[ flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt ] = GTdef_update_slips(earth,modspace,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt);

GTdef_output(foutName,earth,modspace,bt,flt1,flt2,flt3,flt4,flt5,flt6,flt7,subflt,addon,pnt,los,bsl,prf,grd,nod);



Parse_GTdef_out(foutName)


%basename='manualGT';
fsumName = [ basename '_inv.out' ];
GTdef_summary(fsumName,modspace);
    


TT=reshape(T,[Nd,Ns]);   % plotting each Kappa for a quick visual check

subplot(5,5,zz)
image(TT,'CDataMapping','scaled'); title( ['Kp=' num2str(kappa(zz))]); caxis([0,8]); colorbar;



%subplot(5,5,zz)
%plot(tsu.obs,'k'); hold on; plot(tsu_mod,'r'); title( ['Kp=' num2str(kappa(zz))]);

end
figure
 plot(modspace.modinfo(:,11),modspace.modinfo(:,5),'o-'); ylabel('RMS [m]'); xlabel('roughness');
