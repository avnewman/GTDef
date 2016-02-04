function [ modspace,pnt,los,bsl,nod ] = GTdef_forward(modspace,pnt,los,bsl,nod)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_forward                                     %
% Knowing the slip values, calculate the displacements predicted by 	             %
% the model and report the misfit values too				             %
%                                                                                    %
% INPUT:					  		  	             %
% modspace - model structure                                                         %
% pnt      - point structure	  	                                             %
% bsl      - baseline structure                                                      %
% nod      - node structure                                                          %
%									             %
% OUTPUT:                                                                            %
% add modinfo to modspace                                                            %
% modinfo = [ data_num slip_num ndf rss rms wrss wrms chi2 rchi2 r_1d r_2d strain ]  %
% ---------------------------------------------------------------------------------- %
% add .out to pnt, los, bsl, and nod                                                 %
% pnt.out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]                                %
% los.out  - [lon lat zz ULOS eULOS LOSdirE LOSdirN LOSdirV weight]                  %
% bsl.out  - [lon1 lat1 z1 lon2 lat2 z2 Ue Un Uv Ul eUe eUn eUv eUl wgt]             %
% nod.out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]			             %
%                                                                                    %
% first created by Lujia Feng Tue May  5 17:40:11 EDT 2009		             %
% modified roughness from [cm/km] to [cm^2/km^2] by Lujia Feng		             %
% Tue May 19 01:18:30 EDT 2009			  			             %
% added 1st derivative roughness lfeng Thu Dec  3 20:29:31 EST 2009	             %
% corrected roughness unit from [cm^2/km^2] to [cm/km^2]		             %
% used structure lfeng Wed Feb 22 13:37:09 SGT 2012			             %
% added size check for sm/sm_abs with xx lfeng Oct 21 17:30:35 SGT 2014              %
% added modspace structure lfeng Thu Mar 19 17:32:29 SGT 2015                        %
% changed output to pnt, bsl, nod lfeng Fri Mar 20 16:37:00 SGT 2015                 %
% added InSAR los & Lgrn lfeng Tue Nov  3 19:30:36 SGT 2015                          %
% last modified by Lujia Feng Tue Nov  3 19:37:40 SGT 2015                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xgrn    = modspace.Xgrn;
Lgrn    = modspace.Lgrn;
Bgrn    = modspace.Bgrn;
Ngrn    = modspace.Ngrn;
sm      = modspace.sm;
sm_abs  = modspace.sm_abs;
lb      = modspace.lb;
ub      = modspace.ub;
xx      = modspace.xx;
smooth  = modspace.smooth;

% find slips that are included in the model
total = length(xx);
num = 0;
for ii = 1:total
    if xx(ii)==0&&isinf(lb(ii))&&isinf(ub(ii))			% number of slips that are not included
         num = num+1;
    end
end
slip_num = total-num;

% initialize
data_num = 0; 
rss = 0; wrss = 0; chi2 = 0;
rms = 0; wrms = 0; rchi2 = 0;
pnt.out = []; bsl.out = []; nod.out = [];


%%%%% point %%%%%
if ~isempty(Xgrn)
    Xmod = Xgrn*xx; 						% model prediction for point data
    pnt_mod = reshape(Xmod,[],3);
    pnt_mod_err = nan(size(pnt_mod));
    pnt.out = [ pnt.loc pnt_mod pnt_mod_err pnt.wgt ];		% reuse pnt.loc and pnt.wgt

    data_ind = find(~isnan(pnt.obs));				% exclude nan data
    data_num = data_num+length(data_ind);
    pnt_dif2 = (Xmod(data_ind)-pnt.obs(data_ind)).^2; 		% squared residuals/differences
    pnt_err2 = pnt.obs_err(data_ind).^2; 			% squared errors
    rss = rss+sum(pnt_dif2);					% residual sum of squares
    wrss = wrss+sum(pnt.obs_wgt(data_ind).*pnt_dif2./pnt_err2); % weighted rss 
    chi2 = chi2+sum(pnt_dif2./pnt_err2);			% chi-square
end

%%%%% los displacement %%%%%
if ~isempty(Lgrn)
    los_mod = Lgrn*xx;                                          % model prediction for los
    los_mod_err = nan(size(los_mod));
    los.out = [ los.loc los_mod los_mod_err los.dir los.wgt ];  % reuse los.loc, los.dir and los.wgt

    data_ind = find(~isnan(los.obs));				% exclude nan data
    data_num = data_num+length(data_ind);
    los_dif2 = (los_mod(data_ind)-los.obs(data_ind)).^2;        % squared residuals/differences
    los_err2 = los.obs_err(data_ind).^2; 			% squared errors
    rss  = rss+sum(los_dif2);					% residual sum of squares
    wrss = wrss+sum(los.obs_wgt(data_ind).*los_dif2./los_err2); % weighted rss 
    chi2 = chi2+sum(los_dif2./los_err2);			% chi-square
end

%%%%% baseline %%%%%
if ~isempty(Bgrn)
    Bmod = Bgrn*xx; 
    bsl_mod = reshape(Bmod,[],4);
    bsl_mod_err = nan(size(bsl_mod));
    bsl.out = [ bsl.loc bsl_mod bsl_mod_err bsl.wgt ]; 		% reuse bsl.loc and bsl.wgt

    data_ind = find(~isnan(bsl.obs));
    data_num = data_num+length(data_ind);
    bsl_dif2 = (Bmod(data_ind)-bsl.obs(data_ind)).^2;		% squared residuals/differences
    bsl_err2 = bsl.obs_err(data_ind).^2;			% squared errors
    rss = rss+sum(bsl_dif2);					% residual sum of squares
    wrss = wrss+sum(bsl.obs_wgt(data_ind).*bsl_dif2./bsl_err2);	% weighted rss
    chi2 = chi2+sum(bsl_dif2./bsl_err2);			% chi-square
end

%%%%% node %%%%%
if ~isempty(Ngrn)
    Nmod = Ngrn*xx;
    nod_mod = reshape(Nmod,[],3);
    nod_mod_err = nan(size(nod_mod));
    nod_wgt = zeros(size(nod_mod,1),1);
    nod.out = [ nod.loc nod_mod nod_mod_err nod_wgt ]; 		% reuse nod.loc
end

if ~isempty(sm) && size(sm,2) == size(xx,1)
    sm_all = sum(abs(sm*xx)); 
    sm_num = size(sm,1);
    if strcmp(smooth,'2d')
        r_1d = nan;
    	r_2d = 1e8*0.5*sm_all/sm_num;				% NOTE: roughness [cm/km^2]; 0.5 accounts for both ds and ss
    else
        r_1d = 1e5*0.25*sm_all/sm_num;				% NOTE: r_1d [cm/km]; 0.5 accounts for both ds and ss and twice smoothing (two directions)
	r_2d = nan;
    end
else % lfeng Sat Jun 19 04:55:22 EDT 2010
    r_1d = nan; r_2d = nan;
end

if ~isempty(sm_abs) && size(sm_abs,2) == size(xx,1)
    sm_all = sum(abs(sm_abs*xx)); 
    sm_num = size(sm_abs,1);
    strain = 1e5*sm_all/sm_num;                                 % NOTE: strain [cm/km]; 0.5 accounts for both ds and ss
else
    strain = nan;
end

ndf = data_num-slip_num;					% nominal number of degrees of freedom
rms = sqrt(rss/data_num);					% root mean square
wrms = sqrt(wrss/data_num);
rchi2 = chi2/ndf;						% reduced chi-square
modinfo          = [ data_num slip_num ndf rss rms wrss wrms chi2 rchi2 r_1d r_2d strain ];
modspace.modinfo = [ modspace.modinfo; modinfo ];
