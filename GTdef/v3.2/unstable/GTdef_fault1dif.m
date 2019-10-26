function [ modspace,xyzflt,Xgrn ] = GTdef_fault1dif(modspace,...
                       flt,subflt,dipin,Xin,Lin,Bin,Nin,earth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_fault1dif				         %
% Process one distributed-slip type-1 fault, generate its subfaults, and         %
% call GTdef_fault1uni to prepare for the inputs to Matlab function 	         %
% x = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0)                                   %
% Additionally, determine the smoothing matrix for the subfaults	         %
%									         %
% INPUT:					  		  	         %
% flt = [ xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]        %
%   xx,yy - one endpoint among the two endpoints of the master fault             %
%           in the local cartesian coordinate system	  		         %
%   z1  - vertical burial depth (top of fault) >=0                               %  
%   z2  - vertical locking depth (bottom of fault) >=0                           %
%   len - fault length                                                           %
%   str - strike from the endpoint (degree CW from N) [0-360]                    %
%   dip - down from Horiz, right looking from the endpoint [0 180]               %
%   ss  - master-fault strike-slip (left-lateral +)                              %
%   ds  - master-fault dip-slip (thrust +)                                       %
%   ts  - master-fault tensile-slip (opening +)                                  %
%   ss0,ds0,ts0 - lower bounds for master-fault slips			         %
%   ssX,dsX,tsX - upper bounds for master-fault slips			         %
%   Nd  - number of rows defining the subfaults along dip 	                 %
%   Ns  - number of rows defining the subfaults along strike 		         %
% subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		         %
%   dnum - row number for subfaults					         %
%   snum - column number for subfaults				  	         %
%   ss,ds,ts - subfault slips						         %
%   ss0,ds0,ts0,ssX,dsX,tsX - subfault slip bounds			         %
% dipin - dip addon info for the master fault                                    %
%       = [ dip z1 z2 rows ]                                                     %
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
% sm     - smoothing matrix for slips   [slip_num*slip_num]		         %
% sm_abs - matrix for calculating the absolute 1st derivative		         %
% smooth - smoothing method						         %
% surf   - surface smoothing setting					         %
% ------------------------------------------------------------------------------ %
% xyzflt structure                                                               %
% ------------------------------------------------------------------------------ %
% Xgrn - specifically for this fault only                                        %
%                                                                                %
% first created by Lujia Feng Thu Apr 23 21:46:12 EDT 2009		         %
% added 'smooth' by lfeng Wed Dec  2 02:28:33 EST 2009                           %
% added auxiliary sm by lfeng  Wed Dec  2 23:58:49 EST 2009		         %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		         %
% changed 'freesurface' to 'surface' lfeng Wed Feb 24 12:50:51 EST 2010	         %
% added layered earth structure lfeng Tue Feb 28 03:26:05 SGT 2012	         %
% renamed from GTdef_fault3.m to GTdef_fault1dif.m lfeng May 8 SGT 2012          %
% added output Min0 for stress calculation lfeng Thu May 17 SGT 2012             %
% added sm = [] for subflt_num==1 lfeng Tue Oct 21 17:25:33 SGT 2014             %
% added modspace structure lfeng Thu Mar 19 17:32:29 SGT 2015                    %
% added output center points lfeng Thu Mar 19 19:42:23 SGT 2015                  %
% added earth structure lfeng Fri Mar 20 11:39:54 SGT 2015                       %
% corrected GTdef_diffdips error lfeng Mon Mar 23 11:19:57 SGT 2015              %
% added xyzflt lfeng Tue Mar 24 11:21:45 SGT 2015                                %
% changed Xin, Bin, Nin to columnwise lfeng Wed Mar 25 18:49:40 SGT 2015         %
% added Min,SSgrn,DSgrn,TSgrn to xyzflt lfeng Thu Mar 26 15:54:56 SGT 2015       %
% added output Xgrn lfeng Fri Jun 12 12:12:10 SGT 2015                           %
% removed smooth & surf from input lfeng Tue Jun 23 13:02:35 SGT 2015            %
% added InSAR los Lin & Lgrn lfeng Tue Nov  3 11:46:52 SGT 2015                  %
% added Aineq & bineq to modspace lfeng Mon Jun  6 15:46:35 SGT 2016             %
% last modified by Lujia Feng Mon Jun 13 16:18:27 SGT 2016                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], error('GTdef_fault1dif ERROR: need a 1*18 fault vector as input!'); end

smooth = modspace.smooth;
surf   = modspace.surf;

% initialization
% Note: flt is a row vector for the master fault
mx1  = flt(1); my1 = flt(2); 
mz1  = flt(3); mz2 = flt(4); 
mlen = flt(5); mstr = flt(6);  
mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
unit = ones(subflt_num,1);
mx2 = mx1+mlen*sind(mstr);     my2 = my1+mlen*cosd(mstr);	% endpoint 2

% if uniform slip == only one patch - the old fault1 case
if subflt_num==1
    % prjflt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
    prjflt = [ 1 1 mx1 my1 mz1 mz2 mlen mstr mdip mslips ];
    [ ~,xyzflt ] = GTdef_prjflt1uni(prjflt);

    [ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault1uni(earth,[flt(:,1:end-2)],xyzflt,Xin,Lin,Bin,Nin);
    sm_abs = [];
    sm     = [];
    [ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
    return
end

% exclusively specify dips for rows between z1 and z2 depth range
if ~isempty(dipin)
    [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,dipin,Nd,Ns);
else
    xlin  = linspace(mx1,mx2,Ns+1); 
    ylin  = linspace(my1,my2,Ns+1); 
    zlin  = linspace(mz1,mz2,Nd+1)';
    x1mat = xlin(ones(Nd,1),1:end-1);  
    y1mat = ylin(ones(Nd,1),1:end-1);
    z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));
    x1 = reshape(x1mat,[],1); y1 = reshape(y1mat,[],1);
    z1 = reshape(z1mat,[],1); z2 = reshape(z2mat,[],1);
    dip  = mdip*unit; 			% duplicate dips
    dz   = (mz2-mz1)/Nd; 
    ddip = dz/sind(mdip);
end

% duplicate len's, str's,slip's of master-fault
dlen = mlen/Ns;	len = dlen*unit; str = mstr*unit; 
slips = mslips(unit,:);	

% subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]
if ~isempty(subflt)
    num = size(subflt,1); mat = [Nd Ns];
    for ii = 1:num
        dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
        jj = sub2ind(mat,dnum,snum);
	slips(jj,:) = subflt(ii,3:11);
    end
end

% patch order: along dip first, then along strike
newflt = [ x1 y1 z1 z2 len str dip slips ];

dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns));    smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);        snum = reshape(smat,[],1);

% prjflt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
prjflt = [ dnum snum x1 y1 z1 z2 len str dip slips ];
[ ~,xyzflt ] = GTdef_prjflt1uni(prjflt);

[ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault1uni(earth,newflt,xyzflt,Xin,Lin,Bin,Nin);

% create smoothing matrices
if strcmp(surf,'free')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(ddip,dlen,Nd,Ns);
elseif strcmp(surf,'fixed')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed_3slips(ddip,dlen,Nd,Ns);
else
    error('GTdef_fault1dif ERROR: surface smoothing is wrong!!!');
end

if strcmp(smooth,'1d3pf')
    sm = sm_1d3pf;
elseif strcmp(smooth,'1d3pb')
    sm = sm_1d3pb;
elseif strcmp(smooth,'2d') 
    sm = sm_2d;
else
    error('GTdef_fault1dif ERROR: smoothing method is wrong!!!');
end

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them

[ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
