function [ modspace,xyzflt,Xgrn,flt ] = GTdef_fault6(modspace,...
           geoname,colname,flt,subflt,Xin,Lin,Bin,Nin,earth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GTdef_fault6				         %
% Fault type 6 uses external geometry with ss+ds                                 %
% Convert type-6 fault to type-3 fault                                           %
% call GTdef_fault3uni to prepare for the inputs to Matlab function 	         %
% x = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0)                                   %
% Additionally, determine the smoothing matrix for the subfaults	         %
%									         %
% INPUT:					  		  	         %
% geoname - geometry file name                                                   %
% colName - name of each column                                                  %
% flt = [ rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns ]                         %
% subflt = [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]                  %
%   dnum - row number for subfaults					         %
%   snum - column number for subfaults				  	         %
%   rake - rake direction                                                        %
%  rs,ts - subfault slips						         %
%  rake0 - lower bounds for rake [rake-90 rake+90]                               %
%  rakeX - upper bounds for rake [rake-90 rake+90]                               %
%  rs0,ts0,rsX,tsX - subfault slip bounds			                 %
% Xin - point site locations in the local cartesian system 	  	         %
%       [n*3] [ xx yy zz ]						         %
% Lin - los point locations in the local cartesian system + los direction        %
%       [n*6] [xx yy zz dirE dirN dirV]                                          %
% Bin - baseline site locations in the local cartesian system 	  	         %
%       [n*6] [ x1 y1 z1 x2 y2 z2 ]					         %
% Nin - grid and profile node locations in the local cartesian system            %
%       [n*3] [ xx yy zz ]						         %
% Note: z1,z2 depth positive downward                                            %
%       zz elevation positive upward                                             %
% ------------------------------------------------------------------------------ %
% Earth Structure:							         %
% Either of the two types of earth structure can be used 		         %
% earth.type = 'homogeneous'                                                     %
%      earth.rigidity		        (scalar)	{30e9 Pa}                %
%      earth.poisson		        (scalar)	{0.25}		         %
% earth.type = 'layered'	        	        		         %
%      earth.edgrn.nl        	        (scalar)                                 %
%      earth.edgrn.obsz     	        (scalar)                                 %
%      earth.edgrn.nr	                (scalar)                                 %
%      earth.edgrn.minr,edgrn.maxr      (scalar)                                 %
%      earth.edgrn.nz                   (scalar)                                 %
%      earth.edgrn.minz,edgrn.maxz      (scalar)                                 %
%      earth.edgrn.srate                (scalar)			         %
%      earth.layer - [ id depth vp vs ro ]	(nn*5)			         %
%      earth.edgrnfcts - green function library                                  %
% ------------------------------------------------------------------------------ %
%                                                                                %
% OUTPUT:                                                                        %
% ------------------------------------------------------------------------------ %
% modspace structure                                                             %
% slip_num = subflt_num*comp_num = Nd*Ns*comp_num                                %
%    comp_num = 2: two (rake, tensile) components                                %
%    comp_num = 3: three (strike, dip, tensile) components                       %
% Xgrn - displacements [east;north;vertical] for different sites   	         %
%        from unit slips [(3*nn)*slip_num] 				         %
%        (nn is the  number of sites)                                            %
% Lgrn - los displacements [los] for different sites   	                         %
%        from unit slips [(1*nn)*slip_num] 				         %
%        (nn is the number of los points)                                        %
% Bgrn - length changes [east;north;vertical;length] for 	  	         %
%        different baselines from unit slips [(4*nn)*slip_num] 	                 %
%        (nn is the  number of baselines)  				         %
% Ngrn - displacements [east;north;vertical] for different nodes   	         %
%        from unit slips [(3*nn)*slip_num] 				         %
%        (nn is the  number of nodes)                                            %
% Aineq  - left-hand  side matrix for linear inequalities  [(flt_num*2)*slip_num]%
% bineq  - right-hand side vector for linear inequalities  [flt_num*2]           %
% Aeq    - left-hand  side matrix for linear equalities    [slip_num*slip_num]   %
% beq    - right-hand side vector for linear equalities    [slip_num*1]          %
% x0     - initial values for ss,ds,ts 	[slip_num*1]                             %
% xx     - final values for ss,ds,ts 	[slip_num*1]                             %
% lb     - lower bounds for ss,ds,ts 	[slip_num*1]                             %
% ub     - upper bounds for ss,ds,ts	[slip_num*1]			         %
% sm     - smoothing matrix for slips   [slip_num*slip_num]                      %
% sm_abs - matrix for calculating the absolute 1st derivative		         %
% smooth - smoothing method						         %
% surf   - surface smoothing setting					         %
% coord  - coordnate system                                                      %
% origin = [ lon0 lat0 ]                                                         %
% ------------------------------------------------------------------------------ %
% xyzflt structure                                                               %
% ------------------------------------------------------------------------------ %
% Xgrn - specifically for this fault only                                        %
%                                                                                %
% first created by Lujia Feng based on GTdef_fault5.m Apr 20  SGT 2016           %
% allowed rake to be not 0 in external geometry lfeng Wed Jun 1 SGT 2016         %
% added Aineq & bineq to modspace lfeng Mon Jun  6 15:46:35 SGT 2016             %
% added xyzflt.compnum for generating sm lfeng Thu Jun 16 15:54:58 SGT 2016      %
% last modified by Lujia Feng Thu Jun 16 17:02:20 SGT 2016                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 11], error('GTdef_fault6 ERROR: need a 1*11 fault vector as input!'); end

smooth = modspace.smooth;
surf   = modspace.surf;

Nd = flt(10);
Ns = flt(11);

% read geometry file
%             1  2  3  4  5   6   7   8    9  10 11    12    13  14  15  16
% newflt3 = [ xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]
if ~isnan(Nd) && ~isnan(Ns)
    [ ~,newflt3,Nd,Ns ] = GTdef_read_geometry(geoname,colname,modspace.origin,modspace.coord,'GTdef_topleft',Nd,Ns);
else
    [ ~,newflt3,Nd,Ns ] = GTdef_read_geometry(geoname,colname,modspace.origin,modspace.coord,'GTdef_topleft');
    flt(10) = Nd;
    flt(11) = Ns;
end

% initialization
% Note: flt is a row vector for the master fault
x1    = newflt3(:,1); y1  = newflt3(:,2); 
z1    = newflt3(:,3); z2  = newflt3(:,4); 
len   = newflt3(:,5); str = newflt3(:,6);  
dip   = newflt3(:,7);
slips = newflt3(:,8:end);                      % slip block
subflt_num = Nd*Ns;                            % subfault num
unit = ones(subflt_num,1);

% if uniform slip == only one patch
if subflt_num==1
    % if slips are not specified in the external file, use GTdef input file
    if ~any(slips(:,2:end)) % excluding rake, so allow rake to be non-zero
        slips  = flt(1:end-2); 
        newflt = [ newflt3(:,1:7) slips ];
    end
    % prjflt=[dnum snum xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX]
    prjflt = [ 1 1 newflt ];
    [ ~,xyzflt ] = GTdef_prjflt3uni(prjflt);

    [ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault3uni(earth,newflt,xyzflt,Xin,Lin,Bin,Nin);
    sm_abs = [];
    sm     = [];
    [ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
    return
end

% only when slips from geometry file are zero (not specified), slips from GTdef input will be used
if ~any(slips(:,2:end)) % excluding rake, so allow rake to be non-zero
    fprintf(1,'slips in GTdef input are used for modeling\n');
    % master fault
    mslips = flt(1:end-2);
    slips  = mslips(unit,:);

    % subflt = [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]
    if ~isempty(subflt)
        num = size(subflt,1); mat = [Nd Ns];
        for ii = 1:num
            dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
            jj = sub2ind(mat,dnum,snum);
            slips(jj,:) = subflt(ii,3:11);
        end
    end
    newflt = [ newflt3(:,1:7) slips ];
else
    fprintf(1,'slips in %s are used for modeling\n',geoname);
    newflt = newflt3;
end

% create dnum & snum
dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns));    smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);        snum = reshape(smat,[],1);

% prjflt=[dnum snum xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX]
prjflt = [ dnum snum x1 y1 z1 z2 len str dip slips ];
[ ~,xyzflt ] = GTdef_prjflt3uni(prjflt);

[ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault3uni(earth,newflt,xyzflt,Xin,Lin,Bin,Nin);

% create smoothing matrices
width = (z2-z1)./abs(sind(dip));
ddip  = sum(width)/subflt_num;
dlen  = sum(len)/subflt_num;
if strcmp(surf,'free')
    if xyzflt.compnum == 2
        [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_2slips(ddip,dlen,Nd,Ns);
    elseif xyzflt.compnum == 3
        [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(ddip,dlen,Nd,Ns);
    end
elseif strcmp(surf,'fixed')
    if xyzflt.compnum == 2
        [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed_2slips(ddip,dlen,Nd,Ns);
    elseif xyzflt.compnum == 3
        [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed_3slips(ddip,dlen,Nd,Ns);
    end
else
    error('GTdef_fault6 ERROR: surface smoothing is wrong!!!');
end

if strcmp(smooth,'1d3pf')
    sm = sm_1d3pf;
elseif strcmp(smooth,'1d3pb')
    sm = sm_1d3pb;
elseif strcmp(smooth,'2d') 
    sm = sm_2d;
else
    error('GTdef_fault6 ERROR: smoothing method is wrong!!!');
end

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them

[ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
