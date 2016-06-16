function [ modspace,xyzflt,Xgrn,flt ] = GTdef_fault5(modspace,...
           geoname,colname,flt,subflt,Xin,Lin,Bin,Nin,earth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                GTdef_fault5				         %
% Fault type 5 uses external geometry with rake+rs                               %
% Convert type-5 fault to type-1 fault                                           %
% call GTdef_fault1uni to prepare for the inputs to Matlab function 	         %
% x = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0)                                   %
% Additionally, determine the smoothing matrix for the subfaults	         %
%									         %
% INPUT:					  		  	         %
% geoname - geometry file name                                                   %
% colName - name of each column                                                  %
% flt    = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns ]                            %
% subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		         %
%   dnum - row number for subfaults					         %
%   snum - column number for subfaults				  	         %
%   ss,ds,ts - subfault slips						         %
%   ss0,ds0,ts0,ssX,dsX,tsX - subfault slip bounds			         %
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
% Each subfault has three (strike, dip, and tensile) components, so              %
% slip_num = subflt_num*3 = Nd*Ns*3                                              %
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
% sm     - smoothing matrix for slips     [slip_num*slip_num]		         %
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
% first created by Lujia Feng & Paul Morgan Wed Jun 17 15:26:59 SGT 2015         %
% output flt for updating Nd & Ns lfeng Tue Jun 23 17:42:15 SGT 2015             %
% added InSAR los Lin & Lgrn lfeng Tue Nov  3 11:46:52 SGT 2015                  %
% added Aineq & bineq to modspace lfeng Mon Jun  6 15:46:35 SGT 2016             %
% last modified by Lujia Feng Mon Jun 13 16:49:07 SGT 2016                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 11], error('GTdef_fault5 ERROR: need a 1*11 fault vector as input!'); end

smooth = modspace.smooth;
surf   = modspace.surf;

Nd = flt(10);
Ns = flt(11);

% read geometry file
%             1  2  3  4  5   6   7   8  9  10 11  12  13  14  15  16
% newflt1 = [ x1 y1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX ]
if ~isnan(Nd) && ~isnan(Ns)
    [ newflt1,~,Nd,Ns ] = GTdef_read_geometry(geoname,colname,modspace.origin,modspace.coord,'GTdef_topleft',Nd,Ns);
else
    [ newflt1,~,Nd,Ns ] = GTdef_read_geometry(geoname,colname,modspace.origin,modspace.coord,'GTdef_topleft');
    flt(10) = Nd;
    flt(11) = Ns;
end

% initialization
% Note: flt is a row vector for the master fault
x1    = newflt1(:,1); y1  = newflt1(:,2); 
z1    = newflt1(:,3); z2  = newflt1(:,4); 
len   = newflt1(:,5); str = newflt1(:,6);  
dip   = newflt1(:,7);
slips = newflt1(:,8:end);                      % slip block
subflt_num = Nd*Ns;                            % subfault num
unit = ones(subflt_num,1);

% if uniform slip == only one patch - the old fault1 case
if subflt_num==1
    % if slips are not specified in the external file, use GTdef input file
    if ~any(slips)
        slips  = flt(1:end-2); 
        newflt = [ newflt1(:,1:7) slips ];
    end
    % prjflt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
    prjflt = [ 1 1 newflt ];
    [ ~,xyzflt ] = GTdef_prjflt1uni(prjflt);

    [ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault1uni(earth,newflt,xyzflt,Xin,Lin,Bin,Nin);
    sm_abs = [];
    sm     = [];
    [ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
    return
end

% only when slips from geometry file are zero (not specified), slips from GTdef input will be used
if ~any(slips)
    fprintf(1,'slips in GTdef input are used for modeling\n');
    % master fault
    mslips = flt(1:end-2);
    slips  = mslips(unit,:);

    % subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]
    if ~isempty(subflt)
        num = size(subflt,1); mat = [Nd Ns];
        for ii = 1:num
            dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
            jj = sub2ind(mat,dnum,snum);
            slips(jj,:) = subflt(ii,3:11);
        end
    end
    newflt = [ newflt1(:,1:7) slips ];
else
    fprintf(1,'slips in %s are used for modeling\n',geoname);
    newflt = newflt1;
end

% create dnum & snum
dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns));    smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);        snum = reshape(smat,[],1);

% prjflt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
prjflt = [ dnum snum x1 y1 z1 z2 len str dip slips ];
[ ~,xyzflt ] = GTdef_prjflt1uni(prjflt);

[ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault1uni(earth,newflt,xyzflt,Xin,Lin,Bin,Nin);

% create smoothing matrices
width = (z2-z1)./abs(sind(dip));
ddip  = sum(width)/subflt_num;
dlen  = sum(len)/subflt_num;
if strcmp(surf,'free')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(ddip,dlen,Nd,Ns);
elseif strcmp(surf,'fixed')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed_3slips(ddip,dlen,Nd,Ns);
else
    error('GTdef_fault5 ERROR: surface smoothing is wrong!!!');
end

if strcmp(smooth,'1d3pf')
    sm = sm_1d3pf;
elseif strcmp(smooth,'1d3pb')
    sm = sm_1d3pb;
elseif strcmp(smooth,'2d') 
    sm = sm_2d;
else
    error('GTdef_fault5 ERROR: smoothing method is wrong!!!');
end

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them

[ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
