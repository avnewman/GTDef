function [ kappa,rigidity,poisson,...
           flt1_name,flt1_num,flt1,...
	   flt2_name,flt2_num,flt2,...
           flt3_name,flt3_num,flt3,...
	   flt4_name,flt4_num,flt4,...
	   subflt_name,subflt_num,subflt,...
           pnt_name,pnt_num,pnt_loc,pnt_disp,pnt_err,pnt_wgt,...
           bsl_name,bsl_num,bsl_loc,bsl_disp,bsl_err,bsl_wgt,...
           prf_name,prf_num,prf,...
	   grd_name,grd_num,grd ] = GTdef_open(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GTdef_open.m				        %
% 		  function to open the input file for GTdef		        %
%										%
% The format of the input file could be very free. 				%
% Just make sure the first term of each line is an identifiable flag.		%
% The program will ignore unidentifiable flags.					%
%										%
% The function reads in								%
% (1) Parameters: kappa 	(1*kappa_num)		{0}			%
%     		  rigidity	scalar			{30e9 Pa}               %
%		  poisson	scalar			{0.25}			%
%   If not provided by the input file, default values will be used.		%
%										%
% (2) Faults:									% 
% Four types of fault and subfault will be read in separately.			%
%  flt1 - [lon lat z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]     	%
%  flt2 - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]   	%
%  flt3 - [lon1 lat1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]	%
%  flt4 - [lon1 lat1 lon2 lat2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]%
%  subflt - [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  	%
%										%
% (3) Data:									% 
% Point: 	pnt_loc  - [lon lat z]				(nn*3)		% 
%   	 	pnt_disp - [east north vert]			(nn*3)          %
%   	 	pnt_err  - [east north vert]			(nn*3)          %
%   	 	pnt_wgt  - [weight]				(nn*1)          % 
% Baseline: 	bsl_loc  - [lon1 lat1 z1 lon2 lat2 z2]		(nn*6)          %
%   		bsl_disp - [east north vert absolute]   	(nn*4)          %
%   		bsl_err  - [east north vert absolute]   	(nn*4)          %
%   		bsl_wgt  - [weight]                     	(nn*1)          %
% Profile:	prf - [lon1 lat1 lon2 lat2 N]			(nn*5)          %
% Grd:		grd - [Erot Nrot lon1 lat1 lon2 lat2 Ne Nn] 	(nn*8)          %
% Note: 									%
%   Longitude can be input as [0 360] or [-180 180]. The program will convert 	%
%   [0 360] to [-180 180] internally.						%
%   If weight is not assigned, default weight is 1 for all the data.		%
%										%
% first created by Lujia Feng Apr 2009					        %
% last modified by Lujia Feng Fri May  8 01:45:13 EDT 2009			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fin = fopen(filename,'r');

%%%%%%%%%% set default values %%%%%%%%%%
kappa = 0;
rigidity = 30e9;
poisson = 0.25;

%%%%%%%%%% initialize parameters %%%%%%%%%%
flt1_num = 0;  	 flt1_name = '';    flt1 = []; 
flt2_num = 0;    flt2_name = '';    flt2 = []; 
flt3_num = 0;    flt3_name = '';    flt3 = []; 
flt4_num = 0;    flt4_name = '';    flt4 = []; 
subflt_num = 0;  subflt_name = '';  subflt = []; 
pnt_num = 0;     pnt_name = '';     pnt_loc = []; pnt_disp = []; pnt_err = []; pnt_wgt = []; 
bsl_num = 0;     bsl_name = '';     bsl_loc = []; bsl_disp = []; bsl_err = []; bsl_wgt = []; 
prf_num = 0;     prf_name = '';     prf = []; 
grd_num = 0;     grd_name = '';     grd = [];

kappa_num = 0;

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
    %%%%% kappa %%%%%
    if strcmp(flag,'kappa')
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
    %%%%% rigidity %%%%%
    if strcmp(flag,'rigidity')
 	[ rigidity,remain ] = GTdef_read1double(remain);
	continue
    end
    %%%%% poisson %%%%%
    if strcmp(flag,'poisson')
 	[ poisson,remain ] = GTdef_read1double(remain);
        continue
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fault Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% fault %%%%%
    if strcmp(flag,'fault')
        [method,remain] = strtok(remain);
	%% fault 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
	    flt1_num = flt1_num+1; flt1_name = strvcat(flt1_name,name);
	    for ii = 1:16
 		[ flt1(flt1_num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	%% fault 2 %%
	if strcmp(method,'2')
        [name,remain] = strtok(remain);
	    flt2_num = flt2_num+1; flt2_name = strvcat(flt2_name,name);
	    for ii = 1:16
 		[ flt2(flt2_num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	%% fault 3 %%
	if strcmp(method,'3')
        [name,remain] = strtok(remain);
	    flt3_num = flt3_num+1; flt3_name = strvcat(flt3_name,name);
	    for ii = 1:18
 		[ flt3(flt3_num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	%% fault 4 %%
	if strcmp(method,'4')
        [name,remain] = strtok(remain);
	    flt4_num = flt4_num+1; flt4_name = strvcat(flt4_name,name);
	    for ii = 1:18
 		[ flt4(flt4_num,ii),remain ] = GTdef_read1double(remain);
	    end
	    continue
	end
	continue
    end
    %%%%% subfault %%%%%
    if strcmp(flag,'subfault')
        [name,remain] = strtok(remain);
	subflt_num = subflt_num+1; subflt_name = strvcat(subflt_name,name);
	for ii = 1:11
 	    [ subflt(subflt_num,ii),remain ] = GTdef_read1double(remain);
	end
	continue
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Geodetic Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% point %%%%%
    if strcmp(flag,'point')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
            pnt_num = pnt_num+1; pnt_name = strvcat(pnt_name,name);
            %% point location [lon lat z] %%
            for ii = 1:3
 	        [ pnt_loc(pnt_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements [east north vert] %%
	    pnt_disp(pnt_num,1) = nan; pnt_disp(pnt_num,2) = nan;
 	    [ pnt_disp(pnt_num,3),remain ] = GTdef_read1double(remain);
            %% errors [east north vert] %%
	    pnt_err(pnt_num,1) = nan; pnt_err(pnt_num,2) = nan;
 	    [ pnt_err(pnt_num,3),remain ] = GTdef_read1double(remain);
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        pnt_wgt(pnt_num,1) = 1;
	    else
                pnt_wgt(pnt_num,1) = str2double(str);
	    end
            continue
	end
	%% method 2 %%
	if strcmp(method,'2')
            [name,remain] = strtok(remain);
            pnt_num = pnt_num+1; pnt_name = strvcat(pnt_name,name);
            %% point location [lon lat z] %%
            for ii = 1:3
 	        [ pnt_loc(pnt_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements [east north vert] %%
            for ii = 1:2
 	        [ pnt_disp(pnt_num,ii),remain ] = GTdef_read1double(remain);
            end
	    pnt_disp(pnt_num,3) = nan;
            %% errors [east north vert] %%
            for ii = 1:2
 	        [ pnt_err(pnt_num,ii),remain ] = GTdef_read1double(remain);
            end
	    pnt_err(pnt_num,3) = nan;
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        pnt_wgt(pnt_num,1) = 1;
	    else
                pnt_wgt(pnt_num,1) = str2double(str);
	    end
            continue
	end
	%% method 3 %%
	if strcmp(method,'3')
            [name,remain] = strtok(remain);
            pnt_num = pnt_num+1; pnt_name = strvcat(pnt_name,name);
            %% point location [lon lat z] %%
            for ii = 1:3
 	        [ pnt_loc(pnt_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements [east north vert] %%
            for ii = 1:3
 	        [ pnt_disp(pnt_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% errors [east north vert] %%
            for ii = 1:3
 	        [ pnt_err(pnt_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        pnt_wgt(pnt_num,1) = 1;
	    else
                pnt_wgt(pnt_num,1) = str2double(str);
	    end
            continue
	end
	continue
    end
    %%%%% baseline %%%%%
    if strcmp(flag,'baseline')
        [method,remain] = strtok(remain);
	%% method 1 %%
	if strcmp(method,'1')
            [name,remain] = strtok(remain);
            bsl_num = bsl_num+1; bsl_name = strvcat(bsl_name,name);
            %% baseline site locations [lon1 lat1 z1 lon2 lat2 z2] %%
            for ii = 1:6
 	        [ bsl_loc(bsl_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements & length change %%
	    bsl_disp(bsl_num,1) = nan; bsl_disp(bsl_num,2) = nan; bsl_disp(bsl_num,3) = nan;
 	    [ bsl_disp(bsl_num,4),remain ] = GTdef_read1double(remain);
            %% errors %%
	    bsl_err(bsl_num,1) = nan; bsl_err(bsl_num,2) = nan; bsl_err(bsl_num,3) = nan;
 	    [ bsl_err(bsl_num,4),remain ] = GTdef_read1double(remain);
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        bsl_wgt(bsl_num,1) = 1;
	    else
                bsl_wgt(bsl_num,1) = str2double(str);
	    end
            continue
	end
	%% method 2 %%
	if strcmp(method,'2')
            [name,remain] = strtok(remain);
            bsl_num = bsl_num+1; bsl_name = strvcat(bsl_name,name);
            %% baseline site locations [lon1 lat1 z1 lon2 lat2 z2] %%
            for ii = 1:6
 	        [ bsl_loc(bsl_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements & length change %%
            for ii = 1:3
 	        [ bsl_disp(bsl_num,ii),remain ] = GTdef_read1double(remain);
            end
	    bsl_disp(bsl_num,4) = nan;
            %% errors %%
            for ii = 1:3
 	        [ bsl_err(bsl_num,ii),remain ] = GTdef_read1double(remain);
            end
	    bsl_err(bsl_num,4) = nan;
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        bsl_wgt(bsl_num,1) = 1;
	    else
                bsl_wgt(bsl_num,1) = str2double(str);
	    end
            continue
	end
	%% method 3 %%
	if strcmp(method,'3')
            [name,remain] = strtok(remain);
            bsl_num = bsl_num+1; bsl_name = strvcat(bsl_name,name);
            %% baseline site locations [lon1 lat1 z1 lon2 lat2 z2] %%
            for ii = 1:6
 	        [ bsl_loc(bsl_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements & length change %%
            for ii = 1:4
 	        [ bsl_disp(bsl_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% errors %%
            for ii = 1:4
 	        [ bsl_err(bsl_num,ii),remain ] = GTdef_read1double(remain);
            end
            %% weight %%
            [str,remain] = strtok(remain); 
	    if isempty(str)||strncmp(str,'#',1)  % if weight is absent, use default 1
	        bsl_wgt(bsl_num,1) = 1;
	    else
                bsl_wgt(bsl_num,1) = str2double(str);
	    end
            continue
	end
	continue
    end
    %%%%% profile %%%%%
    if strcmp(flag,'profile')
        [name,remain] = strtok(remain);
	prf_num = prf_num+1; prf_name = strvcat(prf_name,name);
	for ii = 1:5
 	    [ prf(prf_num,ii),remain ] = GTdef_read1double(remain);
	end
	continue
    end
    %%%%% grid %%%%%
    if strcmp(flag,'grid')
        [name,remain] = strtok(remain);
	grd_num = grd_num+1; grd_name = strvcat(grd_name,name);
	for ii = 1:8
 	    [ grd(grd_num,ii),remain ] = GTdef_read1double(remain);
	end
	continue
    end
end

% convert longitude [0 360] to [-180 180]
if flt1_num~=0, [ flt1(:,1) ] = GTdef_convertlon(flt1(:,1)); end
if flt2_num~=0
    [ flt2(:,1) ] = GTdef_convertlon(flt2(:,1));
    [ flt2(:,3) ] = GTdef_convertlon(flt2(:,3));
end
if flt3_num~=0, [ flt3(:,1) ] = GTdef_convertlon(flt3(:,1)); end
if flt4_num~=0
    [ flt4(:,1) ] = GTdef_convertlon(flt4(:,1));
    [ flt4(:,3) ] = GTdef_convertlon(flt4(:,3));
end
if pnt_num~=0, [ pnt_loc(:,1) ] = GTdef_convertlon(pnt_loc(:,1)); end
if bsl_num~=0
    [ bsl_loc(:,1) ] = GTdef_convertlon(bsl_loc(:,1));
    [ bsl_loc(:,4) ] = GTdef_convertlon(bsl_loc(:,4));
end
if prf_num~=0
    [ prf(:,1) ] = GTdef_convertlon(prf(:,1));
    [ prf(:,3) ] = GTdef_convertlon(prf(:,3));
end
if grd_num~=0
    [ grd(:,3) ] = GTdef_convertlon(grd(:,3));
    [ grd(:,5) ] = GTdef_convertlon(grd(:,5));
end

fclose(fin);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read1double.m				%
% 	      function to read in one double number from a string		%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ value,remain ] = GTdef_read1double(str)

[str1,remain] = strtok(str);
if isempty(str1)||strncmp(str1,'#',1)  
    disp('The input file is wrong!');
    return; 
end
value = str2double(str1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_convertlon.m				%
% 	      function to convert longitude [0 360] to [-180 180]		%
% INPUT&OUPUT lon is a row or column vector					%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ lon ] = GTdef_convertlon(lon)

ind = find(lon>180); 
lon(ind) = lon(ind)-360;
