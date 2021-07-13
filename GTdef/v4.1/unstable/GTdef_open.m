function [ modspace,earth,...
           flt1,flt2,flt3,flt4,flt5,flt6,flt7,...
           subflt,addon,...
           pnt,los,bsl,prf,grd,...
           sspnt,ssflt1,ssflt2 ] = GTdef_open(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  GTdef_open.m	                                         %
% The format of the input file could be very free                                        %
% Just make sure the first term of each line is an identifiable flag                     %
% The program will ignore unidentifiable flags                                           %
%                                                                                        %
% Note:                                                                                  %
%   Parameter names are not case sensitive                                               %
%   Longitude can be input as [0 360] or [-180 180]                                      %
%   The program will convert [0 360] to [-180 180] internally                            %
%   If weight is not assigned, default weight is 1 for all the data                      %
%                                                                                        %
% OUTPUT:                                                                                %
% modspace - model structure                                                             %
% earth    - earth structure                                                             %
% flt?     - fault structure                                                             %
% subflt   - subfault structure                                                          %
% pnt      - point structure                                                             %
% los      - los structure                                                               %
% bsl      - baseline structure                                                          %
% prf      - profile structure                                                           %
% grd      - grid structure                                                              %
% sspnt    - to work on                                                                  %
% ssflt1   -                                                                             %
% ssflt2   -                                                                             %
%                                                                                        %
%                                                                                        %
% first created by Lujia Feng Apr 2009                                                   %
% added 'coord' flag for coordinate type lfeng Thu Nov  5 20:43:53 EST 2009              %
% added 'smooth' flag for smoothing lfeng Wed Dec  2 02:26:17 EST 2009                   %
% added 'beta' flag for smoothing lfeng Wed Dec  2 23:23:21 EST 2009                     %
% added 'dip' flag for bended fault lfeng Mon Dec  7 00:53:06 EST 2009                   %
% added 'freesurface' flag lfeng Wed Dec  9 17:06:48 EST 2009                            %
% added 'fault5' lfeng Fri Dec 11 12:38:41 EST 2009                                      %
% changed 'freesurface' to 'surface' flag lfeng Wed Feb 24 12:46:01 EST 2010             %
% changed 'coord' to string flag lfeng Wed Feb 24 13:40:01 EST 2010                      %
% use cell array of strings for names lfeng Wed Dec  1 14:41:46 EST 2010                 %
% commented out 'fault5' lfeng Wed Dec  1 14:42:53 EST 2010                              %
% added layered earth structure lfeng Mon Feb 20 17:16:32 SGT 2012                       %
% used structures to simplify parameters lfeng Mon Feb 20 17:30:59 SGT 2012              %
% merged fault1 & fault3 and fault2 & fault4 lfeng Tue May  8 16:20:25 SGT 2012          %
% created new fault3 & fault4 for rake lfeng Tue May  8 18:35:44 SGT 2012                %
% added stress lfeng Thu May 17 07:41:29 SGT 2012                                        %
% added polyconic projection lfeng Thu Jun  7 13:43:28 SGT 2012                          %
% changed flt5 to greensfns lfeng Fri Nov 30 14:37:49 SGT 2012                           %
% added saving greensfns lfeng Mon Aug  5 15:51:38 SGT 2013                              %
% added origin lfeng Thu Dec  5 21:40:47 SGT 2013                                        %
% added addon to combine dip & strike lfeng Fri Oct 24 14:56:14 SGT 2014                 %
% added sweepAngle lfeng Wed Nov  5 19:32:22 SGT 2014                                    %
% added earth structure lfeng Fri Mar 20 11:11:22 SGT 2015                               %
% added modspace structure lfeng Fri Mar 20 17:24:44 SGT 2015                            %
% added modspace.sdropflag lfeng Thu Mar 26 17:27:41 SGT 2015                            %
% added fault5 for external geometry with Paul Morgan lfeng Wed Jun 17 14:03:00 SGT 2015 %
% added InSAR los lfeng Tue Nov  3 10:50:07 SGT 2015                                     %
% added model resultion subfaults for output anewman Mon Tue May 10 11:55:42 EDT 2016    %
% modified code to remove 'continues'. now mostly works with 'elseif' as well as         %
%   improving handling of unkown types/methods/flags.  Warnings report line# of input.   %
%   Comments can now be UNIX or MATLAB style (# or %). anewman May 11 14:04:07 EDT 2016  %
% added optional .mat file output (see GTdef_input) anewman May 18 17:32:55 UTC 2016     %
% added fault6 for external geometry, changed old fault6 to fault7 lfeng Jun 1 SGT 2016  %
% last modified Lujia Feng Wed Jun  1 17:18:32 SGT 2016                                  %
% added options for reading lsqlin params for fitting                                    %
%  default values are now considerably lower (20, 1e-20)                                 %
%last modified by Andrew Newman Tue May 19 09:58:52 EDT 2020                             %
%last modified by Andrew Newman Tue Jun 16 12:40:00 EDT 2020                             %
%last modified by Andrew Newman Thu Jun 18 18:54:03 EDT 2020                             %
% added project line for input                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(filename,'file'), error('GTdef_open ERROR: %s does not exist!',filename); end

fin = fopen(filename,'r');

%%%%%%%%%% initialize parameters %%%%%%%%%%
% modspace structure defaults
modspace.coord     = 'geo';
modspace.origin    = [];
modspace.smooth    = '2d';
modspace.surf      = 'free';
modspace.grnflag   = 'off';
modspace.sdropflag = 'off';
%modspace.kappa    = 0; will be defined later
%modspace.beta     = 0; will be defined later
% form everything that is needed for x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
modspace.Xgrn = []; modspace.Lgrn = [];
modspace.Bgrn = []; modspace.Ngrn = [];
modspace.Aeq  = []; modspace.beq  = [];
modspace.lb   = []; modspace.ub   = [];
modspace.x0   = []; modspace.xx   = [];
modspace.sm   = []; modspace.sm_abs = []; % sm_abs for calculate absolute 1st derivative (strain)
modspace.modinfo = [];
modspace.res = [];
modspace.Rdiags = [];
modspace.mat = [];
modspace.proj = [];
modspace.lsqlin = [];
%modspace.lsqlin_MaxIter = 2000; modspace.lsqlin_TolFun = 1e-30;  % prior hard-wired values
modspace.lsqlin_MaxIter = 20; modspace.lsqlin_TolFun = 1e-20;  % new lower default values (can be overridden by lsqlin call in input)
% earth structure defaults
earth.type      = 'homogeneous';
earth.rigidity  = 30e9;
earth.poisson   = 0.25;
earth.edgrn     = [];
earth.layer     = [];
earth.edgrnfcts = [];

% fault structure
% - to read in
% use CELL ARRAY OF STRINGS for names
flt1.num = 0;  	 flt1.name = {};    flt1.flt = [];    flt1.sdrop = [];
flt2.num = 0;    flt2.name = {};    flt2.flt = [];    flt2.sdrop = [];
flt3.num = 0;    flt3.name = {};    flt3.flt = [];    flt3.sdrop = [];
flt4.num = 0;    flt4.name = {};    flt4.flt = [];    flt4.sdrop = [];
flt5.num = 0;    flt5.name = {};    flt5.flt = [];    flt5.sdrop = [];    flt5.geoname = {};    flt5.colname = {};
flt6.num = 0;    flt6.name = {};    flt6.flt = [];    flt6.sdrop = [];    flt6.geoname = {};    flt6.colname = {};
flt7.num = 0;    flt7.name = {};    flt7.flt = [];    flt7.sdrop = [];    flt7.grname  = {};
% - to built up later
% ordered along dip first, then along strike
flt1.Min = {}; flt1.xyzflt = {};
flt2.Min = {}; flt2.xyzflt = {};
flt3.Min = {}; flt3.xyzflt = {};
flt4.Min = {}; flt4.xyzflt = {};
flt5.Min = {}; flt5.xyzflt = {};
flt6.Min = {}; flt6.xyzflt = {};
flt7.Min = {}; flt7.xyzflt = {};

% subfault structure
subflt.num   = 0;  subflt.name   = {};  subflt.flt = [];

% addon structure
addon.dipnum = 0;  addon.dipname = {};  addon.dip  = [];
addon.strnum = 0;  addon.strname = {};  addon.str  = [];
addon.crt    = [];

% data structure
% - to read in
pnt.num = 0;     pnt.name = {};     pnt.loc = [];   pnt.disp = [];  pnt.err = [];   pnt.wgt = [];
los.num = 0;     los.name = {};     los.loc = [];   loc.disp = [];  los.err = [];   los.wgt = [];  los.dir = [];
bsl.num = 0;     bsl.name = {};     bsl.loc = [];   bsl.disp = [];  bsl.err = [];   bsl.wgt = [];
prf.num = 0;     prf.name = {};     prf.prf = [];
grd.num = 0;     grd.name = {};     grd.grd = [];
% - to be built later
pnt.crt = []; pnt.obs = []; pnt.obs_err = []; pnt.obs_wgt = []; pnt.coef = [];
los.crt = []; los.obs = []; los.obs_err = []; los.obs_wgt = []; los.coef = [];
bsl.crt = []; bsl.obs = []; bsl.obs_err = []; bsl.obs_wgt = []; bsl.coef = [];
%nod.loc = []; nod.crt = []; nod.lon = [];     nod.lat = [];     nod.name = {};

% stress fault
sspnt.num  = 0;  sspnt.name  = {};  sspnt.loc  = []; sspnt.str = []; sspnt.dip = []; sspnt.rake = []; sspnt.fric = [];
ssflt1.fltnum = 0;  ssflt1.num = 0;  ssflt1.fltname = {};  ssflt1.flt = [];
ssflt2.fltnum = 0;  ssflt2.num = 0;  ssflt2.fltname = {};  ssflt2.flt = [];

kappa_num = 0;
beta_num  = 0;
layer_num = 0;
ln = 0; % line number
while(1)
    % read in one line
    tline = fgetl(fin);
    ln=ln+1;
    % test if it is the end of file; exit if yes
    if ischar(tline)~=1, break; end
    % read the 1st term that is a flag for data type; could be some white-spaces before the 1st term
    [flag,remain] = strtok(tline);		% the default delimiter is white-space
    % omit '#,%, and blank' comment lines
    %if (strncmpi(flag,'#',1)||strncmpi(flag,'%',1)||isempty(flag)), continue; end
    if (GTdef_skip(flag)), continue; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Controlling Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% coordinate %%%%%
    if strcmpi(flag,'coord')
        [modspace.coord,remain] = strtok(remain);
        if ~strcmpi(modspace.coord,'geo') && ~strcmpi(modspace.coord,'geo_polyconic') && ~strcmpi(modspace.coord,'local')
            error('GTdef_open ERROR: coordinate system should be either geo, geo_polyconic or local! Line: %d',ln);
        end
    %%%%% origin %%%%%
    elseif strcmpi(flag,'origin')
 	[ lon0,remain ] = GTdef_read1double(remain);
 	[ lat0,remain ] = GTdef_read1double(remain);
	modspace.origin = [ lon0 lat0 ];
    %%%%% smoothing algorithm %%%%%
    elseif strcmpi(flag,'smooth')
        [modspace.smooth,remain] = strtok(remain);
        if ~strcmpi(modspace.smooth,'none') && ~strcmpi(modspace.smooth,'2d') ...
	&& ~strcmpi(modspace.smooth,'1d2pc') && ~strcmpi(modspace.smooth,'1d3pf') && ~strcmpi(modspace.smooth,'1d3pb')
            error('GTdef_open ERROR: smooth algorithm should be 2d, 1d2pc, 1d3pf or 1d3pb! Line: %d',ln);
        end
    %%%%% surface flag %%%%%
    elseif strcmpi(flag,'surface')
        [modspace.surf,remain] = strtok(remain);
        if ~strcmpi(modspace.surf,'none') && ~strcmpi(modspace.surf,'free') && ~strcmpi(modspace.surf,'fixed')
            error('GTdef_open ERROR: surfce should be either free or fixed!');
        end
    %%%%% lsqlin %%%%%
    elseif strcmpi(flag,'lsqlin')
        %[modspace.lsqlin,remain] = strtok(remain);
        [MaxIter,remain] = strtok(remain);
        [TolFun,remain] = strtok(remain);
        modspace.lsqlin_MaxIter = str2double(MaxIter);
        modspace.lsqlin_TolFun = str2double(TolFun);
    %%%%% kappa %%%%%
    elseif strcmpi(flag,'kappa')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
	    while true
                [kk,remain] = strtok(remain);
    		if GTdef_skip(kk), break; end
	        kappa_num = kappa_num+1;
                kappa(kappa_num) = str2double(kk);
	    end
	%% method 2 %%
        elseif strcmp(method,'2')
 	    [ k1,remain ] = GTdef_read1double(remain);
 	    [ kn,remain ] = GTdef_read1double(remain);
 	    [ N,remain ] = GTdef_read1double(remain);
	    delta_k = (kn-k1)/(N-1);
	    n0 = kappa_num+1; kappa_num = kappa_num+N;
	    for ii = n0:kappa_num
	        kappa(ii) = k1+delta_k*(ii-n0);
	    end
        else
	    warning('Line %d, of %s:  Input kappa method "%s" not recognised. Continuing with only kappa=0.',ln,filename,method)
	    kappa_num=1;
	    kappa=0;
	end
    %%%%% beta %%%%%
    elseif strcmpi(flag,'beta')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
	    while true
                [kk,remain] = strtok(remain);
    		if GTdef_skip(kk), break; end
	        beta_num = beta_num+1;
                beta(beta_num) = str2double(kk);
	    end
	%% method 2 %%
        elseif strcmp(method,'2')
 	    [ k1,remain ] = GTdef_read1double(remain);
 	    [ kn,remain ] = GTdef_read1double(remain);
 	    [ N,remain ] = GTdef_read1double(remain);
	    delta_k = (kn-k1)/(N-1);
	    n0 = beta_num+1; beta_num = beta_num+N;
	    for ii = n0:beta_num
	        beta(ii) = k1+delta_k*(ii-n0);
	    end
        else
	    warning('Line %d, of %s:  Input beta method "%s" not recognised. Continuing with only beta=0.',ln,filename,method)
	    beta_num=1;
	    beta=0;
	end
   %%%%% resolution %%%%%
    elseif strncmpi(flag,'res',3)
        [method,remain] = strtok(remain);
        %% method 1 %%
        if strcmp(method,'1')
          %modspace.res='diags';
          modspace.Rdiags=1;
        %% method 2 %%
        elseif strcmp(method,'2')
           rr = strtok(remain);
           if GTdef_skip(rr)   % process on all data (if empty, or starts with # or %)
             modspace.res='all';
           else
             while (~isempty(remain))
                [rr,remain] = strtok(remain);
                if GTdef_skip(rr) , break ; end
                modspace.res=[modspace.res,str2num(rr)];
             end
           end
        else
          warning('Line %d, of %s: Input resolution method "%s" not recognised. Continuing with resolution method = 1',ln,filename,method)
          modspace.res='diags';
        end
    %%%%% resolution %%%%%

    %%%%% projection %%%%%
    elseif strncmpi(flag,'proj',4)
        [modspace.proj,remain] = strtok(remain);
	if ~strcmpi(modspace.proj,'on') && ~strcmpi(modspace.proj,'off')
            error('GTdef_open ERROR: projection should be either on or off!');
	end
    %%%%% projection %%%%%
    %%%%% matlab file %%%%%
    elseif strncmpi(flag,'mat',3)
	[modspace.mat,remain] = strtok(remain);
	if ~strcmpi(modspace.mat,'on') && ~strcmpi(modspace.mat,'off')
            error('GTdef_open ERROR: matfile should be either on or off!');
	end
    %%%%% green's functions %%%%%
    elseif strncmpi(flag,'green',5)
	[modspace.grnflag,remain] = strtok(remain);
	if ~strcmpi(modspace.grnflag,'on') && ~strcmpi(modspace.grnflag,'off')
            error('GTdef_open ERROR: greensfns should be either on or off!');
	end
    %%%%% stress drop %%%%%
    elseif strcmpi(flag,'stressdrop')
	[modspace.sdropflag,remain] = strtok(remain);
	if ~strcmpi(modspace.grnflag,'on') && ~strcmpi(modspace.grnflag,'off')
            error('GTdef_open ERROR: stressdrop should be either on or off!');
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Earth Structures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% earth structure %%%%%
    elseif strcmpi(flag,'earth')
	[earth.type,remain] = strtok(remain);
        %%%%% homogeneous OKADA %%%%%
        %if strcmpi(earth.type,'homogeneous')||strcmpi(earth.type,'homo')
        if strncmpi(earth.type,'homo',4)  % only 1st 4 letters are needed.
            [ earth.rigidity,remain ] = GTdef_read1double(remain);
	    [ earth.poisson,remain ]  = GTdef_read1double(remain);
	    earth.type = 'homogeneous';
        %%%%% layered EDGRN/EDCMP %%%%%
        elseif strncmpi(earth.type,'lay',3);
            [ edgrn.nl,remain ]    = GTdef_read1double(remain);
	    [ edgrn.obsz,remain ]  = GTdef_read1double(remain);
            [ edgrn.nr,remain ]    = GTdef_read1double(remain);
            [ edgrn.minr,remain ]  = GTdef_read1double(remain);
            [ edgrn.maxr,remain ]  = GTdef_read1double(remain);
            [ edgrn.nz,remain ]    = GTdef_read1double(remain);
            [ edgrn.minz,remain ]  = GTdef_read1double(remain);
            [ edgrn.maxz,remain ]  = GTdef_read1double(remain);
            [ edgrn.srate,remain ] = GTdef_read1double(remain);
	    earth.edgrn = edgrn;
        else
           error('GTdef_open ERROR: earth model should be either homogeneous or layered!');
        end
    %%%%% layers %%%%%
    elseif strncmpi(flag,'lay',3)
        layer_num = layer_num+1;
	for ii = 1:5
 	    [ earth.layer(layer_num,ii),remain ] = GTdef_read1double(remain);
	end
    %%%%% warning.  old usage no longer read %%%%%
    elseif (strcmpi(flag,'rigidity') ||strcmpi(flag,'poisson'))
       warning('Line %d, of %s: Old input call for "%s" is deprecated and currently ignored.  Please see "earth" flag for modern usage.',ln,filename,flag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fault Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% fault %%%%%
    elseif strcmpi(flag,'fault')
        [method,remain] = strtok(remain);
	%% fault 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
	    flt1.num = flt1.num+1; flt1.name = [ flt1.name; name ];
	    for ii = 1:18
 		[ flt1.flt(flt1.num,ii),remain ] = GTdef_read1double(remain);
	    end
	%% fault 2 %%
        elseif strcmp(method,'2')
            [name,remain] = strtok(remain);
	    flt2.num = flt2.num+1; flt2.name = [ flt2.name; name ];
	    for ii = 1:18
 		[ flt2.flt(flt2.num,ii),remain ] = GTdef_read1double(remain);
	    end
	%% fault 3 %%
        elseif strcmp(method,'3')
            [name,remain] = strtok(remain);
	    flt3.num = flt3.num+1; flt3.name = [ flt3.name; name ];
	    for ii = 1:18
 		[ flt3.flt(flt3.num,ii),remain ] = GTdef_read1double(remain);
	    end
	%% fault 4 %%
        elseif strcmp(method,'4')
            [name,remain] = strtok(remain);
	    flt4.num = flt4.num+1; flt4.name = [ flt4.name; name ];
	    for ii = 1:18
 		[ flt4.flt(flt4.num,ii),remain ] = GTdef_read1double(remain);
	    end
	%%% fault 5 %%
        elseif strcmp(method,'5')
            [name,remain] = strtok(remain);
	    flt5.num = flt5.num+1; flt5.name = [ flt5.name; name ];
            [geoname,remain] = strtok(remain);
            flt5.geoname = [ flt5.geoname; geoname ];
            [colname,remain] = strtok(remain);
            flt5.colname = [ flt5.colname; colname ];
	    for ii = 1:11
 		[ flt5.flt(flt5.num,ii),remain ] = GTdef_read1double(remain);
	    end
	%%% fault 6 %%
        elseif strcmp(method,'6')
            [name,remain] = strtok(remain);
	    flt6.num = flt6.num+1; flt6.name = [ flt6.name; name ];
            [geoname,remain] = strtok(remain);
            flt6.geoname = [ flt6.geoname; geoname ];
            [colname,remain] = strtok(remain);
            flt6.colname = [ flt6.colname; colname ];
	    for ii = 1:11
 		[ flt6.flt(flt6.num,ii),remain ] = GTdef_read1double(remain);
	    end
	%%% fault 7 %%
        elseif strcmp(method,'7')
            [name,remain] = strtok(remain);
	    flt7.num = flt7.num+1; flt7.name = [ flt7.name; name ];
            [grname,remain] = strtok(remain);
            flt7.grname = [ flt7.grname; grname ];
	    for ii = 1:11
 		[ flt7.flt(flt7.num,ii),remain ] = GTdef_read1double(remain);
	    end
        else
	    warning('Line %d, of %s: Input fault type "%s" not recognised. Ignoring this fault.',ln,filename,method)
        end
    %%%%% subfault %%%%%
    elseif strcmpi(flag,'subfault')
        [name,remain] = strtok(remain);
	subflt.num = subflt.num+1; subflt.name = [ subflt.name; name ];
	for ii = 1:11
 	    [ subflt.flt(subflt.num,ii),remain ] = GTdef_read1double(remain);
	end
    %%%%% dip %%%%%
    elseif strcmpi(flag,'dip')
        [name,remain] = strtok(remain);
        addon.dipnum  = addon.dipnum+1; addon.dipname = [ addon.dipname; name ];
	% (1)dip (2)z1 (3)z2 (4)rows
        for ii = 1:4
 	    [ addon.dip(addon.dipnum,ii),remain ] = GTdef_read1double(remain);
        end
    %%%%% strike %%%%%
    elseif strcmpi(flag,'strike')
        [name,remain] = strtok(remain);
        addon.strnum  = addon.strnum+1; addon.strname = [ addon.strname; name ];
	% (1)lon1 (2)lat1 (3)lon2 (4)lat2 (5)columns (6)sweepAngle
        for ii = 1:6
 	    [ addon.str(addon.strnum,ii),remain ] = GTdef_read1double(remain);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geodetic Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% point %%%%%
    elseif strcmpi(flag,'point')
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
	    if GTdef_skip(str) % if weight is absent, use default 1
	        pnt.wgt(pnt.num,1) = 1;
	    else
                pnt.wgt(pnt.num,1) = str2double(str);
	    end
	%% method 2 %%
        elseif strcmp(method,'2')
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
	    if GTdef_skip(str)  % if weight is absent, use default 1
	        pnt.wgt(pnt.num,1) = 1;
	    else
                pnt.wgt(pnt.num,1) = str2double(str);
	    end
	%% method 3 %%
        elseif strcmp(method,'3')
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
	    if GTdef_skip(str)  % if weight is absent, use default 1
	        pnt.wgt(pnt.num,1) = 1;
	    else
                pnt.wgt(pnt.num,1) = str2double(str);
	    end
        else
	    warning('Line %d, of %s: Input point type "%s" not recognised. Ignoring this data.',ln,filename,method)
	end
    %%%%% los displacement %%%%%
    elseif strcmpi(flag,'los')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
            los.num = los.num+1; los.name = [ los.name; name ];
            %% InSAR point location [lon lat z] %%
            for ii = 1:3
 	        [ los.loc(los.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% los displacement %%
 	    [ los.disp(los.num,1),remain ] = GTdef_read1double(remain);
            %% los displacement error %%
 	    [ los.err(los.num,1),remain ]  = GTdef_read1double(remain);
            %% unit vector for los direction [east north vert] %%
            for ii = 1:3
 	        [ los.dir(los.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% weight %%
            [str,remain] = strtok(remain);
	    if GTdef_skip(str) % if weight is absent, use default 1
	        los.wgt(los.num,1) = 1;
	    else
                los.wgt(los.num,1) = str2double(str);
	    end
        else
	%% no other methods yet %%
	    warning('Line %d, of %s: Input los type "%s" not recognised. Ignoring this data.',ln,filename,method)
	end
    %%%%% baseline %%%%%
    elseif strncmpi(flag,'base',4)
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
	    if GTdef_skip(str)  % if weight is absent, use default 1
	        bsl.wgt(bsl.num,1) = 1;
	    else
                bsl.wgt(bsl.num,1) = str2double(str);
	    end
	%% method 2 %%
        elseif strcmp(method,'2')
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
	    if GTdef_skip(str)  % if weight is absent, use default 1
	        bsl.wgt(bsl.num,1) = 1;
	    else
                bsl.wgt(bsl.num,1) = str2double(str);
	    end
	%% method 3 %%
        elseif strcmp(method,'3')
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
	    if GTdef_skip(str) % if weight is absent, use default 1
	        bsl.wgt(bsl.num,1) = 1;
	    else
                bsl.wgt(bsl.num,1) = str2double(str);
	    end
        else
	%% no other methods yet %%
	    warning('Line %d, of %s: Input baseline type "%s" not recognised. Ignoring this data.',ln,filename,method)
	end
    %%%%% profile %%%%%
    elseif strncmpi(flag,'prof',4)
        [name,remain] = strtok(remain);
	prf.num = prf.num+1; prf.name = [ prf.name; name ];
	for ii = 1:5
 	    [ prf.prf(prf.num,ii),remain ] = GTdef_read1double(remain);
	end
    %%%%% grid %%%%%
    elseif strcmpi(flag,'grid')
        [name,remain] = strtok(remain);
	grd.num = grd.num+1; grd.name = [ grd.name; name ];
	for ii = 1:8
 	    [ grd.grd(grd.num,ii),remain ] = GTdef_read1double(remain);
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stress Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmpi(flag,'stress')
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
        %%%%% fault %%%%%
        elseif strcmpi(type,'fault')
            [method,remain] = strtok(remain);
            %% fault 1 %%
            if strcmp(method,'1')
                [name,remain]  = strtok(remain);
                ssflt1.fltnum  = ssflt1.fltnum+1;
                ssflt1.fltname = [ ssflt1.fltname; name ];
                for ii = 1:11
                    [ ssflt1.flt(ssflt1.fltnum,ii),remain ] = GTdef_read1double(remain);
                end
            %% fault 2 %%
            elseif strcmp(method,'2')
                [name,remain]  = strtok(remain);
                ssflt2.fltnum  = ssflt2.fltnum+1;
                ssflt2.fltname = [ ssflt2.fltname; name ];
                for ii = 1:11
                    [ ssflt2.flt(ssflt2.fltnum,ii),remain ] = GTdef_read1double(remain);
                end
	    else
	        warning('Line %d, of %s: Input stress fault method "%s" not recognised. Ignoring this line.',ln,filename,method)
            end
        else
	    warning('Line %d, of %s: Input stress type "%s" not recognised. Ignoring this line.',ln,filename,type)
        end
   else
	warning('Line %d, of %s: Input flag "%s" not recognised. Ignoring this line.',ln,filename,flag)
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
modspace.beta  = beta;
modspace.kappa = sqrt(beta);


% convert longitude [0 360] to [-180 180]
if strcmpi(modspace.coord,'geo') || strcmpi(modspace.coord,'geo_polyconic')
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

% set out as the same as input first
flt1.out = flt1.flt;
flt2.out = flt2.flt;
flt3.out = flt3.flt;
flt4.out = flt4.flt;
flt5.out = flt5.flt;
flt6.out = flt6.flt;
flt7.out = flt7.flt;
subflt.out     = subflt.flt;
subflt.outname = subflt.name;

% if layer exists, sort layer according to the ascending id
if ~isempty(earth.layer)
   earth.layer = sortrows(earth.layer,1);
end

fclose(fin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read1double.m				%
% 	      function to read in one double number from a string		%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ value,remain ] = GTdef_read1double(str)

[str1,remain] = strtok(str);
if GTdef_skip(str1)
   error('GTdef_open ERROR: the input file is wrong!');
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GTdef_skip.m                                        %
%        function check whether or not to skip reading input information        %
% INPUT 'string'                                                                %
% OUTPUT logical 0 or 1                                                         %
%   checks if empty, starts with # or %                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ out ] = GTdef_skip(str)
if (strncmpi(str,'#',1)||strncmpi(str,'%',1)||isempty(str))
   out = true();
else
   out = false();
end
