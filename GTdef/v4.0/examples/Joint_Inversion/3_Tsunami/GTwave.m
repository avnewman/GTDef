function [dart,tsu,modspace,modinfo] = GTwave(filename,flag)

% flag = 1:  window timeseries
 
[dart,tsu,modspace,fault,flt1,addon,kappa,modinfo] = GTwave_open(filename);

for ii=1:length(kappa) 

if flag == 1
%% Window observed data
tic
fprintf('%s\n','......Windowing Observed Tsunami')
[dart]=GTwave_window_observed(dart);
fprintf('%s\n','......Windowing Observed Tsunami COMPLETE')
toc
%% Window Green's Function Data
tic
fprintf('%s\n','......Windowing Greens Functions')
GTwave_window_greens_fun(dart)
fprintf('%s\n','......Windowing Greens Functions COMPLETE')
toc
end
    
%% build data vector
[dart,tsu] = GTwave_build_data_vector(dart,tsu);

%% build GF matrix
[dart,tsu,modspace] = GTwave_build_gf_matrix(dart,tsu,modspace);

[modspace,tsu] = GTwave_add_weights_erros(modspace,tsu);

%% invert waveform
tic
fprintf('%s\n','......Running Inversion')


    beta = kappa(ii)*kappa(ii);

    fprintf(1,'\n............. doing inversion .............\t');
    fprintf(1,'\n............. ............... .............\t');
    fprintf(1,'\n............. kappa = %10d .............\t',kappa(ii));
    [modspace] = GTwave_invert(modspace,fault,kappa(ii),beta);
    
    % in the tsunami only inversion in GTwave, we only are solving for
    % thrust. So for now we will hard-wire that SS and TS both = 0
    xx=modspace.xx;
    ss=zeros(length(xx),1);  ds=xx; ts=zeros(length(xx),1); 
    subfault_slips=[ss, ds, ts];
    modspace.ss=ss;  modspace.ds=ds; modspace.ts=ts;
toc

%% grab stat analysis

    [modspace,modinfo] = GTwave_modinfo(modspace,tsu,beta,kappa(ii),modinfo);
%% Resolution
    [modspace] = GTwave_Resolution(modspace,fault);
    basename='GTwave_res'; 
    foutName = strcat(basename,'_kp',num2str(kappa(ii)),'_patches_R.out');
    fout     = fopen(foutName,'w');


    fprintf(fout,'#(1)ss[m] (2)ds[m] (3)ts[m]  (4)Rss (5)Rds  (6)Rts (7)R_i\n');
    for ii = 1:length(modspace.ri)
        fprintf(fout,'%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\n',modspace.ss(ii), modspace.ds(ii),modspace.ts(ii),full(modspace.R_ss(ii)),full(modspace.R_ds(ii)),full(modspace.R_ts(ii)), modspace.ri(ii)); 
    end
    fclose(fout);
    
%% waveform outfiles

    %GTwave_DART_out(dart,modspace,kappa(ii))

end
%% wrtie modinfo out file
basename='GTwave'; 
foutNameINV = strcat(basename,'_INV.out');
fout1 = fopen(foutNameINV,'w');
fprintf(fout1,'#(1)beta (2)kappa (3)data_num (4)slip_num (5)ndf (6)rss [m^2] (7)rms [m] (8)wrrs [m^2] (9)wrms [m] (10)chi2 (11)rchi2 (12)r_1d [cm/km] (13)r_2d [cm/km^2]\n');
for i=1:size(modinfo,1) 
fprintf(fout1,'%-10.5e  %-16.5e     %d            %-d      %-d    %-12.5e     %-12.5e   %-12.5e       %-12.5e      %-12.5e  %-12.5e  %-12.5e   %-12.5e\n', modinfo(i,1:13));
end 
fclose(fout1);
end

