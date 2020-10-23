
clc; clear all; close all;
CROPTO=144;  % 228
%% Collect all GF files to make Observation 
M=dir('mat_file/SS_*.mat');
M=M(1:CROPTO);
cd mat_file
% % % for i = 1:length(M)   % ALL SLIP IS ONE
% % %   fname=M(i).name
% % %   load(fname);
% % %    if i==1
% % %       Observed= DART_wave;
% % %    else
% % %        Observed=[Observed + DART_wave];
% % %    end
% % % end

%               Checkerboard function syntax                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I = checkerboard(n,p,q)                                                 %
% n ? Side length in pixels of each square in the checkerboard pattern    %
% p ? Number of rows of tiles in the checkerboard pattern                 %
% q ? Number of columns of tiles in the checkerboard pattern              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 2; p = 6; q = CROPTO/p;   % changed from 38 to 15
 I = checkerboard(n,p,q);
K=I(1:p,1:q);
%imagesc(K)
ID=reshape(K,[1,p*q]); 
Wavespace=[];
 for i = 1:size(M,1)
      if size(M,1) ~= size(ID,2)
          %printf('error, checkerboard and GF mismatch');
      elseif i==1
          fname=M(i).name;
          load(fname);
          Wavespace = DART_wave*0;
         Wavespace = [Wavespace + DART_wave*ID(i)];  
      else
       fname=M(i).name;
       load(fname);
       Wavespace = [Wavespace + DART_wave*ID(i)];
       
      end
 end
Observed=Wavespace;
cd ..
clear n p q Wavespace fname
% Observed_noisy = awgn(Observed,30,'measured');
% plot(t,Observed(:,1),'k'); hold on; plot(t,Observed_noisy(:,1),'r'); ...
%     plot(t,Observed(:,2)+10,'k'); plot(t,Observed_noisy(:,2)+10,'r'); ...
%     plot(t,Observed(:,3)+20,'k'); plot(t,Observed_noisy(:,3)+20,'r');
%     

save('Observed_Tsunami.mat','Observed','t');  % SAVE

%% WINDOWING OBSERVED DATA
WINDOW=[1500, 2190; 765, 1680; 150, 600];
FNAME=['Dart32401.txt';'Dart32402.txt';'Dart32403.txt'];
for i = 1: length(WINDOW)
Window_Observed_Tsunami(WINDOW(i,1),WINDOW(i,2), FNAME(i,:))
end

%% WINDOWING EACH GF
cd GFs
k=1;
while 1

M=dir(sprintf('%s%03.f%s','SS_',k,'_*.txt'));
    if ~isempty(M)
        for i = 1:length(M)
            M(i).name
            Window_Observed_Tsunami(WINDOW(i,1),WINDOW(i,2), M(i).name); 
        end
    k=k+1;
    else
        fprintf('%s\n','end of windowing');
        break
    end
end
 cd ..
clear k i M
%% BEGIN BUILDING MATRICES FOR INVERSION

%% IMPORT PRJFLT and XYZFLT FROM GTDEF RUN
load prjflt.mat
prjflt=prjflt(1:CROPTO,:);
%%
kappa=[1000:1000:15000]; 

modinfo=[];
for z=1:length(kappa)
%% DATA VECTOR (DART)
% lines up the detieded, windowed tsunami data
load('WIN_Observed_Tsunami.mat')
Wavespace=[WIN_Observed_32401; WIN_Observed_32402; WIN_Observed_32403];%WIN_Observed_32402];     %%%%%%%%%%%%%% OBSERVED
    
    tsu.WinDart32401=[WIN_Observed_32401 , WIN_t_32401];
    tsu.WinDart32402=[WIN_Observed_32402 , WIN_t_32402];
    tsu.WinDart32403=[WIN_Observed_32403 , WIN_t_32403];
    
    load Observed_Tsunami.mat
    tsu.Dart32401=[Observed(:,1),t];
    tsu.Dart32402=[Observed(:,2),t];
    tsu.Dart32403=[Observed(:,3),t];

tsu.obs=Wavespace; 
tsu.wgt=ones(length(tsu.obs),1);                                   % add weight
tsu.obs_err=ones(length(tsu.obs),1)*.1;                            % add error
tsu.coef=sqrt(tsu.wgt)./tsu.obs_err;

clear  t Wavespace WIN_Observed_32401 WIN_Observed_32402 WIN_Observed_32403 WIN_t_32401 WIN_t_32402 WIN_t_32403

%% Green's Functions Matrix
% set up Green's Function matrix. It contains data for all 50 patches 

DIR=dir('mat_file/WIN_SS*.mat'); 
DIR=DIR(1:CROPTO);
slip.coeff_loc=zeros(length(DIR),1);
num_patch=size(DIR,1);
Tgrn=zeros(length(tsu.obs),num_patch); 
    cd mat_file/
        for i=1:length(DIR)
            load(DIR(i).name);
            slip.coeff_loc(i)=str2num((DIR(i).name(8:10))); 
            win_timeseries=[WIN_DART_32401; WIN_DART_32402; WIN_DART_32403];% WIN_DART_32402];  %%%%% GF
            Tgrn(:,slip.coeff_loc(i))= win_timeseries;
        end
    cd ..
clear DIR  win_timeseries i num_gauges num_points WIN_DART_32401 WIN_DART_32402 WIN_DART_32403 WIN_t_32401 WIN_t_32402 WIN_t_32403

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
clear eqNum ii ind mod_Tgrn mod_tsu

lb=zeros(1,num_patch);
ub=ones(1,num_patch);  
x0=1*ones(1,num_patch);

slip_num = num_patch; % number of patches	
xunit = x0; xunit(xunit~=0) = 1;		% set unit slips for nonzero initial slips
Aineq = zeros(slip_num*2,slip_num); bineq = zeros(slip_num*2,1);   % useful only for fault type 3 & 4
Aeq   = zeros(slip_num,slip_num);  beq   = zeros(slip_num,1);

%%  SMOOTHING
beta = kappa(z)*kappa(z);
dd=20000; ds=20000; Nd=6; Ns=CROPTO/Nd;     % values needed to smooth main faults
[ sm ] = GTdef_sm2d_free(dd,ds,Nd,Ns);

sm  = GTdef_condense(sm);
sm = sm.*beta;
d_sm = zeros(size(sm,1),1);   % append zeros to data vector

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them

nrealdata=size(C,1);
C = [ C;sm ]; d = [ d;d_sm ];  % and add smothing matrix to GF matrix
clear dd ds  i ind_fixed 

%% WAVEFORM INVERSION
fprintf(1,'\n............. doing inversion .............\t');
fprintf(1,'\n............. ............... .............\t');
fprintf(1,'\n............. kappa = %10d .............\t',kappa(z));
options = optimset('MaxIter',2000,'TolFun',1e-30);
d=double(d);
C=full(C);
%[xx,resnorm] = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
[xx,resnorm] = lsqlin(C,d,[],[],[],[],lb,ub,x0,options);

%%  Output waveforms

wave_out=C*xx;
subplot(3,1,1)
 plot(tsu.Dart32401(:,2),tsu.Dart32401(:,1),'k'); hold on; plot(tsu.WinDart32401(:,2),tsu.WinDart32401(:,1),'r'); title('DART 32401');
subplot(3,1,2)
 plot(tsu.Dart32402(:,2),tsu.Dart32402(:,1),'k'); hold on; plot(tsu.WinDart32402(:,2),tsu.WinDart32402(:,1),'r'); title('DART 32402');
 subplot(3,1,3)
 plot(tsu.Dart32403(:,2),tsu.Dart32403(:,1),'k'); hold on; plot(tsu.WinDart32403(:,2),tsu.WinDart32403(:,1),'r');title('DART 32403');


%%
% in this case we only have ds
ss=zeros(length(xx),1);  ds=xx; ts=zeros(length(xx),1); 
subfault_slips=[ss, ds, ts];

%% Modinfo
 %find slips that are included in the model
 total = length(xx);
 num = 0;
 for ii = 1:total
      if xx(ii)==0&&isinf(lb(ii))&&isinf(ub(ii))                  % number of slips that are not included
           num = num+1;
      end
 end
 slip_num = total-num;
tsu.out=[]; data_num=0;
rss = 0; wrss = 0; chi2 = 0;
rms = 0; wrms = 0; rchi2 = 0;

 if ~isempty(Tgrn)
     tsu_mod = Tgrn*xx;                                          % model prediction for los
      tsu_mod_err = nan(size(tsu_mod));
      tsu.out = [tsu_mod tsu_mod_err  tsu.wgt ];  % reuse los.loc, los.dir and los.wgt
  
      data_ind = find(~isnan(tsu.obs));                           % exclude nan data
      data_num = data_num+length(data_ind);
      tsu_dif2 = (tsu_mod(data_ind)-tsu.obs(data_ind)).^2;        % squared residuals/differences
      tsu_err2 = tsu.obs_err(data_ind).^2;                        % squared errors
      rss  = rss+sum(tsu_dif2);                                   % residual sum of squares
      wrss = wrss+sum(tsu.wgt(data_ind).*tsu_dif2./tsu_err2); % weighted rss 
      chi2 = chi2+sum(tsu_dif2./tsu_err2);                        % chi-square
  end


if ~isempty(sm) && size(sm,2) == size(xx,1)
     sm_all = sum(abs(sm*xx)); 
     sm_num = size(sm,1);
         r_1d = nan;
         r_2d = 1e8*0.5*sm_all/sm_num;                           % NOTE: roughness [cm/km^2]; 0.5 accounts for both ds and ss
else % lfeng Sat Jun 19 04:55:22 EDT 2010
     r_1d = nan; r_2d = nan;
end

ndf = data_num-slip_num;                                        % nominal number of degrees of freedom
rms = sqrt(rss/data_num);                                       % root mean square
wrms = sqrt(wrss/data_num);
rchi2 = chi2/ndf;                                               % reduced chi-square
modinfo = [ modinfo ;beta kappa(z) data_num slip_num ndf rss rms wrss wrms chi2 rchi2 r_1d r_2d  ];
clear beta data_num slip_num ndf rss rms wrss wrms chi2 rchi2 r_1d r_2d ...
    sm_all sm_num tsu_dif2 tsu_err2 data_ind tsu_mod tsu_mod_err num ii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GTdef Resolution

ix=find(isfinite(lb'));
id=[[1:nrealdata]';ix+nrealdata];
 C0=C(id,:);                        % operator matrix (includes smoothing)
 d=d(1:nrealdata);
 x0=x0(ix);
 
 % Develop weighted and damped over-determined generalized inverse matrix, Gg
 % for the model resolution matrix  (See Menke, 2012 -- Chapter 4)
 G=C0(1:nrealdata,:);                 % strip off the smoothing part of the matrix
     
 e2I=sm.^2;                          % weighted smoothing matrix
 Gg=inv(G'*G+e2I)*G';                % overdetermined.

 R=Gg*G ;                            % Model Resolution Matrix 
 N=G*Gg ;                            % Data Resolution Matrix
 
 R_ds=diag(R);
 R_ss=zeros(length(R_ds),1);   % TEMPORARY FIX WHEN SS AND TS ARE ADDED TO INVERSION
 R_ts=zeros(length(R_ds),1);
 R_out=[R_ss,R_ds,R_ts];
 
 % resolution spread parameter
 Li=30;  % km PATCH SIZE for a 30x30 subfault size
 r_i=Li./sqrt(R_out(:,2)); % resolution spread parameter (in km)
%%
 basename='INV_tsun'; 
foutName = strcat(basename,'_kp',num2str(kappa(z)),'_patches_R.out');
fout     = fopen(foutName,'w');
fprintf(fout,'#(1)dnum   (2)snum     (3)xtop1       (4)ytop1      (5)ztop1       (6)xbot1       (7)ybot1      (8)zbot1      (9)xbot2 (10)ybot2 (11)zbot2 (12)xtop2 (13)ytop2 (14)ztop2 (15)xctr (16)yctr (17)zcrt (18)ss[m] (19)ds[m] (20)ts[m]  (21)Rss (22)Rds  (23)Rts (24)R_i\n');
for ii = 1:length(prjflt)
fprintf(fout,'    %4d    %4d     %12.5f  %11.5f   %12.3e    %12.5f   %11.5f  %12.3e    %12.5f   %11.5f    %12.3e     %12.5f    %11.5f    %12.3e   %12.5f    %11.5f  %12.3e   %10.5f   %10.5f     %10.5f  %6.4f %6.4f  %6.4f %6.4f\n',prjflt(ii,1:17),subfault_slips(ii,1:3), R_out(ii,1:3), r_i(ii)); 
end
fclose(fout);

end
%%
basename='INV_tsun'; 
foutNameINV = strcat(basename,'_INV.out');
fout1 = fopen(foutNameINV,'w');
fprintf(fout1,'#(1)beta (2)kappa (3)data_num (4)slip_num (5)ndf (6)rss [m^2] (7)rms [m] (8)wrrs [m^2] (9)wrms [m] (10)chi2 (11)rchi2 (12)r_1d [cm/km] (13)r_2d [cm/km^2]\n');
for i=1:size(modinfo,1) 
fprintf(fout1,'%-10.5e  %-16.5e     %d            %-d      %-d    %-12.5e     %-12.5e   %-12.5e       %-12.5e      %-12.5e  %-12.5e  %-12.5e   %-12.5e\n', modinfo(i,1:13));
end 
fclose(fout1);
figure
kappa_temp=num2str(modinfo(:,2)); kappa_txt=cellstr(kappa_temp);
dx=0; dy=0;
plot(modinfo(:,13),modinfo(:,7)*100,'k'); hold on; scatter(modinfo(:,13),modinfo(:,7)*100,'r'); xlabel('Roughness [cm/km^2]'); ...
    ylabel('RMS [cm]'); title('Kappa L-Curve');
text(double(modinfo(:,13)+dx),double((modinfo(:,7)*100)+dy), kappa_txt)

clear Aeq Aineq basenam beq bineq C0 ds fout fout1 foutName foutNameINV  ...
    i id ii ix lb mod_tsu mot_Tgrn Ns num nrealdata num_patch  ...
    options R_ds R_ss R_ts resnorm slip ts ss xunit z
