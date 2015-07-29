function [ mod_info,pnt_out,bsl_out,nod_out ] ...
            = GTdef_forward(Xgrn,Bgrn,Ngrn,sm,lb,ub,xx,...
            		    pnt_loc,pnt_obs,pnt_obs_err,pnt_wgt,...
	    		    bsl_loc,bsl_obs,bsl_obs_err,bsl_wgt,nod_loc)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_forward				  %
% Knowing the slip values, calculate the displacements predicted by 	  %
% the model and report the misfit values too				  %
%									  %
% INPUT:					  		  	  %
% Each fault has three (strike, dip, and tensile) components, so          %
% slip_num = flt_num*3                                                    %
%  xx   - final values for ss,ds,ts 	[slip_num*1]                      %
%  lb  - lower bounds for ss,ds,ts 	[slip_num*1]                      %
%  ub  - upper bounds for ss,ds,ts	[slip_num*1]			  %
%  Xgrn - displacements [east;north;vertical] for different sites   	  %
%         from unit slips [(3*nn)*slip_num] 				  %
%         (nn is the  number of sites)  				  %
%  Bgrn - length changes [east;north;vertical;length] for 	  	  %
%         different baselines from unit slips [(4*nn)*slip_num] 	  %
%         (nn is the  number of baselines)  				  %
%  Ngrn - displacements [east;north;vertical] for different nodes   	  %
%         in the profile or grid from unit slips [(3*nn)*slip_num] 	  %
%         (nn is the  number of nodes)  				  %
%  sm  - condensed smoothing matrix with rows of all zeros removed	  %
% Point:    pnt_loc     - [lon lat z]			(nn*3)	  	  % 
%	    pnt_obs     - [east;north;vert]		(3nn*1)		  %
%   	    pnt_obs_err - [east;north;vert]		(3nn*1)    	  %
%   	    pnt_wgt     - [weight]			(nn*1)      	  % 
% Baseline: bsl_loc     - [lon1 lat1 z1 lon2 lat2 z2]	(nn*6)    	  %
%   	    bsl_obs     - [east;north;vert;length]      (4nn*1)    	  %
%   	    bsl_obs_err - [east;north;vert;length]      (4nn*1)    	  %
%   	    bsl_wgt     - [weight]                      (nn*1)      	  %
% Profile   nod_loc	- [lon lat z]			(nn*3)	  	  %
% & Grid:								  %
%									  %
% OUTPUT:                                                                 %
%mod_info = [data_num slip_num ndf rss rms wrss wrms chi2 rchi2 roughness]%
%    slip_num - number of free slips					  %
%    data_num - number of data points (not including nan)                 %
%    ndf      - nominal number of degrees of freedom 			  %
%  Note: since if we introduce smoothing, slips are not independent	  %  
%  we don't really know the real ndf.					  %
%    rss      - residual sum of squares                          	  %
%    rms      - root mean square of rss = sqrt(rss/data_num)              %
%    wrss     - weighted residual sum of squares			  %
%    wrms     = sqrt(wrss/data_num)					  %
%    chi2     - chi-square 		 				  %
%    rchi2    - reduced chi-square = chi2/ndf                             %
%    roughness - the average second-order finite-difference sum of each   % 
%                fault patch (eq. 5 in Jonsson_etal_BSSA_2002) [cm^2/km^2]%
%  pnt_out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]                    %
%  bsl_out  - [lon1 lat1 z1 lon2 lat2 z2 Ue Un Uv Ul eUe eUn eUv eUl wgt] %
%  nod_out  - [lon lat zz Ue Un Uv eUe eUn eUv weight]			  %
%                                                                         %
% first created by Lujia Feng Tue May  5 17:40:11 EDT 2009		  %
% last modified by Lujia Feng Tue May 12 18:57:42 EDT 2009		  %
% modified roughness from [cm/km] to [cm^2/km^2] by Lujia Feng		  %
% Tue May 19 01:18:30 EDT 2009			  			  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
pnt_out = []; bsl_out = []; nod_out = [];

if ~isempty(Xgrn)
    Xmod = Xgrn*xx; 						% model prediction for point data
    pnt_mod = reshape(Xmod,[],3);
    pnt_mod_err = nan(size(pnt_mod));
    pnt_out = [ pnt_loc pnt_mod pnt_mod_err pnt_wgt ];		% reuse pnt_loc and pnt_wgt

    data_ind = find(~isnan(pnt_obs));				% exclude nan data
    data_num = data_num+length(data_ind);
    pnt_dif2 = (Xmod(data_ind)-pnt_obs(data_ind)).^2; 		% squared residuals/differences
    pnt_err2 = pnt_obs_err(data_ind).^2; 			% squared errors
    pnt_wgt3 = [pnt_wgt;pnt_wgt;pnt_wgt]; 			% (3*n)*1 weight vector [east;north;vertical]
    rss = rss+sum(pnt_dif2);					% residual sum of squares
    wrss = wrss+sum(pnt_wgt3(data_ind).*pnt_dif2./pnt_err2); 	% weighted rss 
    chi2 = chi2+sum(pnt_dif2./pnt_err2);			% chi-square
end
if ~isempty(Bgrn)
    Bmod = Bgrn*xx; 
    bsl_mod = reshape(Bmod,[],4);
    bsl_mod_err = nan(size(bsl_mod));
    bsl_out = [ bsl_loc bsl_mod bsl_mod_err bsl_wgt ]; 		% reuse bsl_loc and bsl_wgt

    data_ind = find(~isnan(bsl_obs));
    data_num = data_num+length(data_ind);
    bsl_dif2 = (Bmod(data_ind)-bsl_obs(data_ind)).^2;		% squared residuals/differences
    bsl_err2 = bsl_obs_err(data_ind).^2;			% squared errors
    bsl_wgt4 = [bsl_wgt;bsl_wgt;bsl_wgt;bsl_wgt];		% (4*n)*1 weight vector [east;north;vertical;length]
    rss = rss+sum(bsl_dif2);					% residual sum of squares
    wrss = wrss+sum(bsl_wgt4(data_ind).*bsl_dif2./bsl_err2);	% weighted rss
    chi2 = chi2+sum(bsl_dif2./bsl_err2);			% chi-square
end
if ~isempty(Ngrn)
    Nmod = Ngrn*xx;
    nod_mod = reshape(Nmod,[],3);
    nod_mod_err = nan(size(nod_mod));
    nod_wgt = zeros(size(nod_mod,1),1);
    nod_out = [ nod_loc nod_mod nod_mod_err nod_wgt ]; 		% reuse nod_loc
end
if ~isempty(sm)
    sm_all = sum(abs(sm*xx)); 
    sm_num = size(sm,1);
    roughness = 1e10*0.5*sm_all/sm_num;				% NOTE: roughness [cm^2/km^2]
else
    roughness = nan;
end

ndf = data_num-slip_num;					% nominal number of degrees of freedom
rms = sqrt(rss/data_num);					% root mean square
wrms = sqrt(wrss/data_num);
rchi2 = chi2/ndf;						% reduced chi-square
mod_info = [ data_num slip_num ndf rss rms wrss wrms chi2 rchi2 roughness ];
