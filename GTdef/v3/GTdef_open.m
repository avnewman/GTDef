function [ coord,origin,smooth,surf,beta,grnflag,...
           rigidity,poisson,...
	   earth,edgrn,layer,...
           flt1,flt2,flt3,flt4,flt5,...
	   subflt,dip,...
           pnt,bsl,prf,grd,...
	   sspnt,ssflt1,ssflt2 ] = GTdef_open(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_open.m				                 %
% The format of the input file could be very free  				         %
% Just make sure the first term of each line is an identifiable flag 		         %
% The program will ignore unidentifiable flags 					         %
%										         %
% The function reads in  							         %
% (1) Parameters are not case sensitive						         %
% 		coord         string                  {geo}			         %
%                  geo   - geographic coordiante				         %
% 		   geo_polyconic - geographic using polyconic projection                 %
% 		   local - cartesian coordinate				                 %
%               origin   = [ lon0 lat0 ]  (1x2)                                          %
%               smooth	string		        {2d}				         %
%         	   1d2pc - 1st derivative 2-point central			         %
%         	   1d3pf - 1st derivative 3-point forward                                %
%         	   1d3pb - 1st derivative 3-point backward                               %
%            	      2d - 2nd derivative 3-point central                                %
%               surf          string			{free}		   	         %
%                  fixed - no free-surface					         %
% 		   free  - assume free-surface				                 %
%               kappa 	(1*kappa_num)		{0}				         %
%               beta 	(1*beta_num)		{0}				         %
% Note: beta is used to weight smoothing, usually for 1st derivative		         %
%       kappa is used to weight smoothing, usually for 2nd derivative		         %
%   If not provided by the input file, default values in {} will be used.                %
%               greensfns - output green's functions {off}                               %
%										         %
% (2) Earth Structure:								         %
% Either of the two types of earth structure can be used.			         %
%  earth = homogeneous								         %
%     		rigidity		(scalar)		{30e9 Pa}                %
%		poisson			(scalar)		{0.25}		         %
%  earth = layered		        					         %
%		edgrn.nl        	(scalar)                                         %
% 		edgrn.obsz     		(scalar)                                         %
% 		edgrn.nr	        (scalar)                                         %
% 		edgrn.minr,edgrn.maxr   (scalar)                                         %
% 		edgrn.nz                (scalar)                                         %
% 		edgrn.minz,edgrn.maxz   (scalar)                                         %
%		edgrn.srate		(scalar)				         %
%    		layer - [ id depth vp vs ro ]	(nn*5)				         %
%		                                                                         %
% (3) Faults:									         % 
% read in four types of fault and subfault separately 				         %
% each fault has a structure to store corresponding data			         %
% flt? structure: flt?.name flt?.num flt?.flt					         %
% subflt structure: subflt.name subflt.num & subflt.flt				         %
%										         %
% flt1.flt - [lon1 lat1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]        %
% flt2.flt - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]      %
%      subflt.flt - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		         %
% flt3.flt - [lon1 lat1 z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]  %
% flt4.flt - [lon1 lat1 lon2 lat2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns]%
%      subflt.flt - [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]                 %
% flt5.flt - [ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]                                    %
% dip structure: dip.name dip.num & dip.dip					         %
% dip.dip  - [ dip z1 z2 rows ]	need to be used with dip.name			         %
%										         %
% (4) Data:									         % 
% Point:                                                                                 %
%		pnt.num  - number of point data			(scalar)	         %
% 		pnt.name - names of point data			(cell array)	         %
%		pnt.loc  - [lon lat z]				(nn*3)		         % 
%   	 	pnt.disp - [east north vert]			(nn*3)                   %
%   	 	pnt.err  - [east north vert]			(nn*3)                   %
%   	 	pnt.wgt  - [weight]				(nn*1)                   % 
% Baseline: 									         %
% 		bsl.num  - number of baselines 			(scalar)	         %
% 		bsl.name - names of baselines			(cell array)	         %
% 		bsl.loc  - [lon1 lat1 z1 lon2 lat2 z2]		(nn*6)                   %
%   		bsl.disp - [east north vert absolute]   	(nn*4)                   %
%   		bsl.err  - [east north vert absolute]   	(nn*4)                   %
%   		bsl.wgt  - [weight]                     	(nn*1)                   %
% Profile:									         %
% 		prf.num  - number of profiles 			(scalar)	         %
% 		prf.name - names of profiles			(cell array)	         %
% 		prf.prf  - [lon1 lat1 lon2 lat2 N]		(nn*5)                   %
% Grd:										         %
% 		grd.num  - number of grids 			(scalar)	         %
% 		grd.name - names of grids			(cell array)	         %
%		grd.grd  - [Erot Nrot lon1 lat1 lon2 lat2 Ne Nn](nn*8)                   %
%										         %
% (5) Stress:                                                                            %
% Point:                                                                                 %
%             sspnt.num  - number of points for stress calculation (scalar)              %
%             sspnt.name - names of points                      (cell array)             %
%             sspnt.loc  - [lon lat z]                          (nn*3)                   %
%             sspnt.str  - [strike]                             (nn*1)                   %
%             sspnt.dip  - [dip]                                (nn*1)                   %
%             sspnt.rake - [rake]                               (nn*1)                   %
% 	      sspnt.fric - [friction]                           (nn*1)                   %
% Fault:                                                                                 %
%            ssflt?.fltnum  - number of faults for stress calculation  (scalar)          %    
%            ssflt?.fltname - names of faults                   (cell array)             %
%            ssflt?.flt  - fault parameters                     (flt_num*11)             %
%              									         %
%            ssflt1.flt = [lon1 lat1 z1   z2   len str dip rake fric Nd Ns]              %
%            ssflt2.flt = [lon1 lat1 lon2 lat2 z1  z2  dip rake fric Nd Ns]              %
%                                                                                        %
% Note: 									         %
%   Longitude can be input as [0 360] or [-180 180]. The program will convert 	         %
%   [0 360] to [-180 180] internally.						         %
%   If weight is not assigned, default weight is 1 for all the data.		         %
%										         %
% first created by Lujia Feng Apr 2009					                 %
% added 'coord' flag for coordinate type lfeng Thu Nov  5 20:43:53 EST 2009	         %
% added 'smooth' flag for smoothing lfeng Wed Dec  2 02:26:17 EST 2009 		         %
% added 'beta' flag for smoothing lfeng Wed Dec  2 23:23:21 EST 2009		         %
% added 'dip' flag for bended fault lfeng Mon Dec  7 00:53:06 EST 2009		         %
% added 'freesurface' flag lfeng Wed Dec  9 17:06:48 EST 2009			         %
% added 'fault5' lfeng Fri Dec 11 12:38:41 EST 2009				         %
% changed 'freesurface' to 'surface' flag lfeng Wed Feb 24 12:46:01 EST 2010	         %
% changed 'coord' to string flag lfeng Wed Feb 24 13:40:01 EST 2010		         %
% use cell array of strings for names lfeng Wed Dec  1 14:41:46 EST 2010	         %
% commented out 'fault5' lfeng Wed Dec  1 14:42:53 EST 2010			         %
% added layered earth structure lfeng Mon Feb 20 17:16:32 SGT 2012                       %
% used structures to simplify parameters lfeng Mon Feb 20 17:30:59 SGT 2012              %
% merged fault1 & fault3 and fault2 & fault4 lfeng Tue May  8 16:20:25 SGT 2012          %
% created new fault3 & fault4 for rake lfeng Tue May  8 18:35:44 SGT 2012                %
% added stress lfeng Thu May 17 07:41:29 SGT 2012                                        %
% added polyconic projection lfeng Thu Jun  7 13:43:28 SGT 2012                          %
% changed flt5 to greensfns lfeng Fri Nov 30 14:37:49 SGT 2012                           %
% added saving greensfns lfeng Mon Aug  5 15:51:38 SGT 2013                              %
% added origin lfeng Thu Dec  5 21:40:47 SGT 2013                                        %
% last modified by Lujia FENG Thu Dec  5 21:43:20 SGT 2013                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(filename,'file'), error('GTdef_open ERROR: %s does not exist!',filename); end

fin = fopen(filename,'r');

%%%%%%%%%% set default values %%%%%%%%%%
coord    = 'geo';
origin   = [];
smooth   = '2d';
surf     = 'free';
grnflag  = 'off';
rigidity = 30e9;
poisson  = 0.25;
earth    = 'homogeneous';
edgrn    = [];
layer    = [];

%%%%%%%%%% initialize parameters %%%%%%%%%%
% use CELL ARRAY OF STRINGS for names
flt1.num = 0;  	 flt1.name = {};    flt1.flt = []; 
flt2.num = 0;    flt2.name = {};    flt2.flt = []; 
flt3.num = 0;    flt3.name = {};    flt3.flt = []; 
flt4.num = 0;    flt4.name = {};    flt4.flt = []; 
flt5.num = 0;    flt5.name = {};    flt5.flt = [];   flt5.grname = {};
subflt.num = 0;  subflt.name = {};  subflt.flt = []; 
dip.num = 0;     dip.name = {};     dip.dip = [];
pnt.num = 0;     pnt.name = {};     pnt.loc = [];   pnt.disp = [];  pnt.err = [];   pnt.wgt = []; 
bsl.num = 0;     bsl.name = {};     bsl.loc = [];   bsl.disp = [];  bsl.err = [];   bsl.wgt = []; 
prf.num = 0;     prf.name = {};     prf.prf = []; 
grd.num = 0;     grd.name = {};     grd.grd = [];
sspnt.num  = 0;  sspnt.name  = {};  sspnt.loc  = []; sspnt.str = []; sspnt.dip = []; sspnt.rake = []; sspnt.fric = [];
ssflt1.fltnum = 0;  ssflt1.num = 0;  ssflt1.fltname = {};  ssflt1.flt = []; 
ssflt2.fltnum = 0;  ssflt2.num = 0;  ssflt2.fltname = {};  ssflt2.flt = []; 

kappa_num = 0;
beta_num  = 0;
layer_num = 0;

while(1)   
    % read in one line
    tline = fgetl(fin);
    % test if it is the end of file; exit if yes
    if ischar(tline)~=1, break; end
    % read the 1st term that is a flag for data type; could be some white-spaces before the 1st term
    [flag,remain] = strtok(tline);		% the default delimiter is white-space
    % omit '#' comment lines
    if strncmp(flag,'#',1), continue; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Controlling Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% coordinate %%%%%
    if strcmpi(flag,'coord')
        [coord,remain] = strtok(remain);
        if ~strcmpi(coord,'geo') && ~strcmpi(coord,'geo_polyconic') && ~strcmpi(coord,'local')
            error('GTdef_open ERROR: coordinate system should be either geo, geo_polyconic or local!');
        end
        continue
    end
    %%%%% origin %%%%%
    if strcmpi(flag,'origin')
 	[ lon0,remain ] = GTdef_read1double(remain);
 	[ lat0,remain ] = GTdef_read1double(remain);
	origin = [ lon0 lat0 ];
        continue
    end
    %%%%% smoothing algorithm %%%%%
    if strcmpi(flag,'smooth')
        [smooth,remain] = strtok(remain);
        if ~strcmpi(smooth,'none') && ~strcmpi(smooth,'2d') && ~strcmpi(smooth,'1d2pc') && ~strcmpi(smooth,'1d3pf') && ~strcmpi(smooth,'1d3pb')
            error('GTdef_open ERROR: smooth algorithm should be 2d, 1d2pc, 1d3pf or 1d3pb!');
        end
        continue
    end
    %%%%% surface flag %%%%%
    if strcmpi(flag,'surface')
        [surf,remain] = strtok(remain);
        if ~strcmpi(surf,'none') && ~strcmpi(surf,'free') && ~strcmpi(surf,'fixed')
            error('GTdef_open ERROR: surfce should be either free or fixed!');
        end
        continue
    end
    %%%%% kappa %%%%%
    if strcmpi(flag,'kappa')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
	    while true
                [kk,remain] = strtok(remain);
    		if isempty(kk)||strncmp(kk,'#',1), break; end
	        kappa_num = kappa_num+1;
                kappa(kappa_num) = str2double(kk);             
	    end
	    continue
	end
	%% method 2 %%
	if strcmp(method,'2')
 	    [ k1,remain ] = GTdef_read1double(remain);
 	    [ kn,remain ] = GTdef_read1double(remain);
 	    [ N,remain ] = GTdef_read1double(remain);
	    delta_k = (kn-k1)/(N-1);
	    n0 = kappa_num+1; kappa_num = kappa_num+N; 
	    for ii = n0:kappa_num
	        kappa(ii) = k1+delta_k*(ii-n0);
	    end
	    continue
	end
	continue
    end
    %%%%% beta %%%%%
    if strcmpi(flag,'beta')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
	    while true
                [kk,remain] = strtok(remain);
    		if isempty(kk)||strncmp(kk,'#',1), break; end
	        beta_num = beta_num+1;
                beta(beta_num) = str2double(kk);             
	    end
	    continue
	end
	%% method 2 %%
	if strcmp(method,'2')
 	    [ k1,remain ] = GTdef_read1double(remain);
 	    [ kn,remain ] = GTdef_read1double(remain);
 	    [ N,remain ] = GTdef_read1double(remain);
	    delta_k = (kn-k1)/(N-1);
	    n0 = beta_num+1; beta_num = beta_num+N; 
	    for ii = n0:beta_num
	        beta(ii) = k1+delta_k*(ii-n0);
	    end
	    continue
	end
	continue
    end
    %%%%% green's functions %%%%%
    if strcmpi(flag,'greensfns')
	[grnflag,remain] = strtok(remain);
	if ~strcmpi(grnflag,'on') && ~strcmpi(grnflag,'off')
            error('GTdef_open ERROR: greensfns should be either on or off!');
	end
	continue
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Earth Structures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% earth structure %%%%%
    if strcmpi(flag,'earth')
	[earth,remain] = strtok(remain);
        %%%%% homogeneous OKADA %%%%%
        if strcmpi(earth,'homogeneous')||strcmpi(earth,'homo')
            [ rigidity,remain ] = GTdef_read1double(remain);
	    [ poisson,remain ]  = GTdef_read1double(remain);
	    earth = 'homogeneous';
            continue
        end
        %%%%% layered EDGRN/EDCMP %%%%%
        if strcmpi(earth,'layered');
            [ edgrn.nl,remain ]    = GTdef_read1double(remain);
	    [ edgrn.obsz,remain ]  = GTdef_read1double(remain);
            [ edgrn.nr,remain ]    = GTdef_read1double(remain);
            [ edgrn.minr,remain ]  = GTdef_read1double(remain);
            [ edgrn.maxr,remain ]  = GTdef_read1double(remain);
            [ edgrn.nz,remain ]    = GTdef_read1double(remain);
            [ edgrn.minz,remain ]  = GTdef_read1double(remain);
            [ edgrn.maxz,remain ]  = GTdef_read1double(remain);
            [ edgrn.srate,remain ] = GTdef_read1double(remain);
            continue
        end
        error('GTdef_open ERROR: earth model should be either homogeneous or layered!');
    end
    %%%%% layers %%%%%
    if strcmpi(flag,'layer')
        layer_num = layer_num+1;
	for ii = 1:5
 	    [ layer(layer_num,ii),remain ] = GTdef_read1double(remain);
	end
        continue
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fault Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% fault %%%%%
    if strcmpi(flag,'fault')
        [method,remain] = strtok(remain);
	%% fault 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
	    flt1.num = flt1.num+1; flt1.name = [ flt1.name; name ];
	    for ii = 1:18
 		[ flt1.flt(flt1.num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	%% fault 2 %%
	if strcmp(method,'2')
            [name,remain] = strtok(remain);
	    flt2.num = flt2.num+1; flt2.name = [ flt2.name; name ];
	    for ii = 1:18
 		[ flt2.flt(flt2.num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	%% fault 3 %%
	if strcmp(method,'3')
            [name,remain] = strtok(remain);
	    flt3.num = flt3.num+1; flt3.name = [ flt3.name; name ];
	    for ii = 1:18
 		[ flt3.flt(flt3.num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	%% fault 4 %%
	if strcmp(method,'4')
            [name,remain] = strtok(remain);
	    flt4.num = flt4.num+1; flt4.name = [ flt4.name; name ];
	    for ii = 1:18
 		[ flt4.flt(flt4.num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	%%% fault 5 %%
	if strcmp(method,'5')
            [name,remain] = strtok(remain);
	    flt5.num = flt5.num+1; flt5.name = [ flt5.name; name ];
            [grname,remain] = strtok(remain);
            flt5.grname = [ flt5.grname; grname ];
	    for ii = 1:11
 		[ flt5.flt(flt5.num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	continue
    end
    %%%%% subfault %%%%%
    if strcmpi(flag,'subfault')
        [name,remain] = strtok(remain);
	subflt.num = subflt.num+1; subflt.name = [ subflt.name; name ];
	for ii = 1:11
 	    [ subflt.flt(subflt.num,ii),remain ] = GTdef_read1double(remain);
	end
	continue
    end
    %%%%% dip %%%%%
    if strcmpi(flag,'dip')
        [name,remain] = strtok(remain);
	dip.num = dip.num+1; dip.name = [ dip.name; name ];
	for ii = 1:4
 	    [ dip.dip(dip.num,ii),remain ] = GTdef_read1double(remain);
	end
	continue
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geodetic Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% point %%%%%
    if strcmpi(flag,'point')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
            pnt.num = pnt.num+1; pnt.name = [ pnt.name; name ];
            %% point location [lon lat z] %%
            for ii = 1:3
 	        [ pnt.loc(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements [east north vert] %%
	    pnt.disp(pnt.num,1) = nan; pnt.disp(pnt.num,2) = nan;
 	    [ pnt.disp(pnt.num,3),remain ] = GTdef_read1double(remain);
            %% errors [east north vert] %%
	    pnt.err(pnt.num,1) = nan; pnt.err(pnt.num,2) = nan;
 	    [ pnt.err(pnt.num,3),remain ] = GTdef_read1double(remain);
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        pnt.wgt(pnt.num,1) = 1;
	    else
                pnt.wgt(pnt.num,1) = str2double(str);
	    end
            continue
	end
	%% method 2 %%
	if strcmp(method,'2')
            [name,remain] = strtok(remain);
            pnt.num = pnt.num+1; pnt.name = [ pnt.name; name ];
            %% point location [lon lat z] %%
            for ii = 1:3
 	        [ pnt.loc(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements [east north vert] %%
            for ii = 1:2
 	        [ pnt.disp(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
	    pnt.disp(pnt.num,3) = nan;
            %% errors [east north vert] %%
            for ii = 1:2
 	        [ pnt.err(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
	    pnt.err(pnt.num,3) = nan;
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        pnt.wgt(pnt.num,1) = 1;
	    else
                pnt.wgt(pnt.num,1) = str2double(str);
	    end
            continue
	end
	%% method 3 %%
	if strcmp(method,'3')
            [name,remain] = strtok(remain);
            pnt.num = pnt.num+1; pnt.name = [ pnt.name; name ];
            %% point location [lon lat z] %%
            for ii = 1:3
 	        [ pnt.loc(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements [east north vert] %%
            for ii = 1:3
 	        [ pnt.disp(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% errors [east north vert] %%
            for ii = 1:3
 	        [ pnt.err(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        pnt.wgt(pnt.num,1) = 1;
	    else
                pnt.wgt(pnt.num,1) = str2double(str);
	    end
            continue
	end
	continue
    end
    %%%%% baseline %%%%%
    if strcmpi(flag,'baseline')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
            bsl.num = bsl.num+1; bsl.name = [ bsl.name; name ];
            %% baseline site locations [lon1 lat1 z1 lon2 lat2 z2] %%
            for ii = 1:6
 	        [ bsl.loc(bsl.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements & length change %%
	    bsl.disp(bsl.num,1) = nan; bsl.disp(bsl.num,2) = nan; bsl.disp(bsl.num,3) = nan;
 	    [ bsl.disp(bsl.num,4),remain ] = GTdef_read1double(remain);
            %% errors %%
	    bsl.err(bsl.num,1) = nan; bsl.err(bsl.num,2) = nan; bsl.err(bsl.num,3) = nan;
 	    [ bsl.err(bsl.num,4),remain ] = GTdef_read1double(remain);
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        bsl.wgt(bsl.num,1) = 1;
	    else
                bsl.wgt(bsl.num,1) = str2double(str);
	    end
            continue
	end
	%% method 2 %%
	if strcmp(method,'2')
            [name,remain] = strtok(remain);
            bsl.num = bsl.num+1; bsl.name = [ bsl.name; name ];
            %% baseline site locations [lon1 lat1 z1 lon2 lat2 z2] %%
            for ii = 1:6
 	        [ bsl.loc(bsl.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements & length change %%
            for ii = 1:3
 	        [ bsl.disp(bsl.num,ii),remain ] = GTdef_read1double(remain);
            end
	    bsl.disp(bsl.num,4) = nan;
            %% errors %%
            for ii = 1:3
 	        [ bsl.err(bsl.num,ii),remain ] = GTdef_read1double(remain);
            end
	    bsl.err(bsl.num,4) = nan;
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        bsl.wgt(bsl.num,1) = 1;
	    else
                bsl.wgt(bsl.num,1) = str2double(str);
	    end
            continue
	end
	%% method 3 %%
	if strcmp(method,'3')
            [name,remain] = strtok(remain);
            bsl.num = bsl.num+1; bsl.name = [ bsl.name; name ];
            %% baseline site locations [lon1 lat1 z1 lon2 lat2 z2] %%
            for ii = 1:6
 	        [ bsl.loc(bsl.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements & length change %%
            for ii = 1:4
 	        [ bsl.disp(bsl.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% errors %%
            for ii = 1:4
 	        [ bsl.err(bsl.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        bsl.wgt(bsl.num,1) = 1;
	    else
                bsl.wgt(bsl.num,1) = str2double(str);
	    end
            continue
	end
	continue
    end
    %%%%% profile %%%%%
    if strcmpi(flag,'profile')
        [name,remain] = strtok(remain);
	prf.num = prf.num+1; prf.name = [ prf.name; name ];
	for ii = 1:5
 	    [ prf.prf(prf.num,ii),remain ] = GTdef_read1double(remain);
	end
	continue
    end
    %%%%% grid %%%%%
    if strcmpi(flag,'grid')
        [name,remain] = strtok(remain);
	grd.num = grd.num+1; grd.name = [ grd.name; name ];
	for ii = 1:8
 	    [ grd.grd(grd.num,ii),remain ] = GTdef_read1double(remain);
	end
	continue
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stress Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(flag,'stress')
        [type,remain] = strtok(remain);
        %%%%% point %%%%%
        if strcmpi(type,'point')
            [name,remain] = strtok(remain);
            sspnt.num = sspnt.num+1; sspnt.name = [ sspnt.name; name ];
            %% point location [lon lat z] %%
            for ii = 1:3
                [ sspnt.loc(sspnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% point fault parameters %%
            [ sspnt.str(sspnt.num,1),remain ]  = GTdef_read1double(remain);
            [ sspnt.dip(sspnt.num,1),remain ]  = GTdef_read1double(remain);
            [ sspnt.rake(sspnt.num,1),remain ] = GTdef_read1double(remain);
            [ sspnt.fric(sspnt.num,1),remain ] = GTdef_read1double(remain);
            continue
        end
        %%%%% fault %%%%%
        if strcmpi(type,'fault')
            [method,remain] = strtok(remain);
            %% fault 1 %%
            if strcmp(method,'1')
                [name,remain]  = strtok(remain);
                ssflt1.fltnum  = ssflt1.fltnum+1;
                ssflt1.fltname = [ ssflt1.fltname; name ];
                for ii = 1:11
                    [ ssflt1.flt(ssflt1.fltnum,ii),remain ] = GTdef_read1double(remain);
                end
                continue
            end
            %% fault 2 %%
            if strcmp(method,'2')
                [name,remain]  = strtok(remain);
                ssflt2.fltnum  = ssflt2.fltnum+1;
                ssflt2.fltname = [ ssflt2.fltname; name ];
                for ii = 1:11
                    [ ssflt2.flt(ssflt2.fltnum,ii),remain ] = GTdef_read1double(remain);
                end
                continue
            end
        end
    end
end

% convert kappa^2 to beta
if (kappa_num>0)&&(beta_num>0)
    beta = [ beta kappa.^2 ];
elseif kappa_num>0
    beta = kappa.^2;
elseif (kappa_num==0)&&(beta_num==0)
    beta = 0;
end

% convert longitude [0 360] to [-180 180]
if strcmpi(coord,'geo') || strcmpi(coord,'geo_polyconic')
   if flt1.num~=0, [ flt1.flt(:,1) ] = GTdef_convertlon(flt1.flt(:,1)); end
   if flt2.num~=0
       [ flt2.flt(:,1) ] = GTdef_convertlon(flt2.flt(:,1));
       [ flt2.flt(:,3) ] = GTdef_convertlon(flt2.flt(:,3));
   end
   if flt3.num~=0, [ flt3.flt(:,1) ] = GTdef_convertlon(flt3.flt(:,1)); end
   if flt4.num~=0
       [ flt4.flt(:,1) ] = GTdef_convertlon(flt4.flt(:,1));
       [ flt4.flt(:,3) ] = GTdef_convertlon(flt4.flt(:,3));
   end
   if pnt.num~=0, [ pnt.loc(:,1) ] = GTdef_convertlon(pnt.loc(:,1)); end
   if bsl.num~=0
       [ bsl.loc(:,1) ] = GTdef_convertlon(bsl.loc(:,1));
       [ bsl.loc(:,4) ] = GTdef_convertlon(bsl.loc(:,4));
   end
   if prf.num~=0
       [ prf.prf(:,1) ] = GTdef_convertlon(prf.prf(:,1));
       [ prf.prf(:,3) ] = GTdef_convertlon(prf.prf(:,3));
   end
   if grd.num~=0
       [ grd.grd(:,3) ] = GTdef_convertlon(grd.grd(:,3));
       [ grd.grd(:,5) ] = GTdef_convertlon(grd.grd(:,5));
   end
end

% if layer exists, sort layer according to the ascending id
if ~isempty(layer)
   layer = sortrows(layer,1);
end

fclose(fin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read1double.m				%
% 	      function to read in one double number from a string		%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ value,remain ] = GTdef_read1double(str)

[str1,remain] = strtok(str);
if isempty(str1)||strncmp(str1,'#',1)  
    error('GTdef_open ERROR: The input file is wrong!');
end
value = str2double(str1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_convertlon.m				%
% 	      function to convert longitude [0 360] to [-180 180]		%
% INPUT&OUPUT lon is a row or column vector					%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ lon ] = GTdef_convertlon(lon)

ind = lon>180;
lon(ind) = lon(ind)-360;
