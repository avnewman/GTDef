function [modspace,modinfo] = GTwave_modinfo(modspace,tsu,beta,kappa,modinfo)
%find slips that are included in the model
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 xx=modspace.xx;  lb=modspace.lb;  ub=modspace.ub;
 Tgrn=modspace.Tgrn;
 
 sm=modspace.sm;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 total = length(xx);
 num = 0;
 for ii = 1:total
      if xx(ii)==0&&isinf(lb(ii))&&isinf(ub(ii))                  % number of slips that are not included
           num = num+1;
      end
 end
 slip_num = total-num;
 
 %%
 
    tsu.out=[];     data_num=0;
    rss = 0;        wrss = 0;   chi2 = 0;
    rms = 0;        wrms = 0;   rchi2 = 0;

 if ~isempty(Tgrn)
     tsu_mod = Tgrn*xx;                                          % model prediction for los
      tsu_mod_err = nan(size(tsu_mod));
      
      tsu.out = [tsu_mod tsu_mod_err  tsu.obs_wgt ];  % reuse los.loc, los.dir and los.wgt
  
      data_ind = find(~isnan(tsu.obs));                           % exclude nan data
      data_num = data_num+length(data_ind);
      tsu_dif2 = (tsu_mod(data_ind)-tsu.obs(data_ind)).^2;        % squared residuals/differences
      tsu_err2 = tsu.obs_err(data_ind).^2;                        % squared errors
      rss  = rss+sum(tsu_dif2);                                   % residual sum of squares
      wrss = wrss+sum(tsu.obs_wgt(data_ind).*tsu_dif2./tsu_err2); % weighted rss 
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
modinfo = [ modinfo ;beta kappa data_num slip_num ndf rss rms wrss wrms chi2 rchi2 r_1d r_2d  ];
end

