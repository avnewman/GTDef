function [ modspace,xyzflt,Xgrn ] = GTdef_fault2dif(modspace,...
                       flt,subflt,dipin,strin,Xin,Lin,Bin,Nin,earth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_fault2dif				         %
% Process one distributed-slip type-2 fault, generate its subfaults and          %
% call GTdef_fault2uni to prepare for the inputs to Matlab function 	         %
% x = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0)                                   %
% Additionally, determine the smoothing matrix for the subfaults	         %
%									         %
% INPUT:					  		  	         %
% flt = [ x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns ]         %
%   x1,y1 - one endpoint among the two endpoints of the master fault             %
%   x2,y2 - the other endpoint among the two endpoints                           %
%           both in the local cartesian coordinate system	                 %
%   z1  - vertical burial depth (top of fault) >=0                               %  
%   z2  - vertical locking depth (bottom of fault) >=0                           %
%   dip - down from Horiz, right looking from the endpoint 1 [0 180]             %
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
% strin - strike addon info for the master fault                                 %
%       = [ x1 y1 x2 y2 columns sweepAngle ]                                     %
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
% Bgrn - length changes [east;north;vertical;length] for                         %
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
% first created by Lujia Feng Fri Apr 24 10:44:57 EDT 2009		         %
% added "smooth" by lfeng Wed Dec  2 02:39:38 EST 2009                           %
% added auxiliary sm by lfeng Thu Dec  3 00:14:33 EST 2009		         %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		         %
% changed 'freesurface' to 'surface' lfeng Wed Feb 24 13:16:24 EST 2010	         %
% added layered earth structure lfeng Tue Feb 28 03:41:26 SGT 2012	         %
% renamed from GTdef_fault4.m to GTdef_fault2dif.m lfeng May 8 SGT 2012          %
% added output Min0 for stress calculation lfeng Thu May 17 SGT 2012             %
% added sm = [] for subflt_num==1 lfeng Tue Oct 21 17:25:33 SGT 2014             %
% added strin lfeng Fri Oct 24 15:24:11 SGT 2014                                 %
% added sweepAngle lfeng Wed Nov  5 19:35:13 SGT 2014                            %
% added modspace structure lfeng Thu Mar 19 17:32:29 SGT 2015                    %
% added output center points lfeng Thu Mar 19 19:42:23 SGT 2015                  %
% added earth structure lfeng Fri Mar 20 11:39:54 SGT 2015                       %
% added xyzflt lfeng Tue Mar 24 11:21:45 SGT 2015                                %
% changed Xin, Bin, Nin to columnwise lfeng Wed Mar 25 18:49:40 SGT 2015         %
% added Min,SSgrn,DSgrn,TSgrn to xyzflt lfeng Thu Mar 26 15:54:56 SGT 2015       %
% added output Xgrn lfeng Fri Jun 12 12:12:10 SGT 2015                           %
% removed smooth & surf from input lfeng Tue Jun 23 13:02:35 SGT 2015            %
% added InSAR los Lin & Lgrn lfeng Tue Nov  3 11:46:52 SGT 2015                  %
% added Aineq & bineq to modspace lfeng Mon Jun  6 15:46:35 SGT 2016             %
% last modified by Lujia Feng Mon Jun 13 16:23:30 SGT 2016                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], error('GTdef_fault2dif ERROR: need a 1*18 fault vector as input!'); end

smooth = modspace.smooth;
surf   = modspace.surf;

% initialization
% Note: flt is a 1-row vector for the master fault
mx1 = flt(1); my1 = flt(2); mx2 = flt(3); my2 = flt(4); 
mz1 = flt(5); mz2 = flt(6); mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
mstr = GTdef_strike(mx1,my1,mx2,my2);
unit = ones(subflt_num,1);

% if uniform slip == only one patch - the old fault2 case
if subflt_num==1
    % flt = [dnum snum x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
    prjflt = [ 1 1 mx1 my1 mx2 my2 mz1 mz2 mdip mslips ];
    [ ~,xyzflt ] = GTdef_prjflt2uni(prjflt);

    [ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault2uni(earth,[flt(:,1:end-2)],xyzflt,Xin,Lin,Bin,Nin);
    sm_abs = [];
    sm     = [];
    [ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
    return
end

% exclusively specify strikes for columns
if isempty(strin)
    % exclusively specify dips for rows between z1 and z2 depth range
    if ~isempty(dipin)
        [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,dipin,Nd,Ns);
    else
        [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_samedips(mx1,my1,mx2,my2,mz1,mz2,mstr,mdip,Nd,Ns);
    end
    dlen = sqrt((mx1-mx2)^2+(my1-my2)^2)/Ns;
else
    x1=[]; y1=[]; x2=[]; y2=[]; z1=[]; z2=[]; dip=[]; ddip = 0; dlen = 0;
    colsum = sum(strin(:,5)); % can be several mini-segments for each segment
    if colsum ~= Ns, error('GTdef_fault2diff ERROR: strike is not specified correctly!'); end
    % loop through each segment not the mini ones
    strnum = size(strin,1);
    for ii=1:strnum
        sx1 = strin(ii,1); sy1 = strin(ii,2);
        sx2 = strin(ii,3); sy2 = strin(ii,4);
        sNs = strin(ii,5); sweepAngle = strin(ii,6);
        sstr  = GTdef_strike(sx1,sy1,sx2,sy2);
        cdlen = sqrt((sx1-sx2)^2+(sy1-sy2)^2)/sNs;
        if ~isempty(dipin)
            [ cx1,cy1,cx2,cy2,cz1,cz2,cdip,cddip ] = GTdef_diffdips(sx1,sy1,sx2,sy2,mz1,mz2,sstr,dipin,Nd,sNs,sweepAngle);
        else
            [ cx1,cy1,cx2,cy2,cz1,cz2,cdip,cddip ] = GTdef_samedips(sx1,sy1,sx2,sy2,mz1,mz2,sstr,mdip,Nd,sNs,sweepAngle);
        end
        x1 = [x1;cx1]; y1 = [y1;cy1]; 
        x2 = [x2;cx2]; y2 = [y2;cy2]; 
        z1 = [z1;cz1]; z2 = [z2;cz2]; 
        dip  = [dip;cdip];
        ddip = ddip+cddip;
        dlen = dlen+cdlen;
    end
    ddip = ddip/strnum;
    dlen = dlen/strnum;
end

% form subfaults
slips = mslips(unit,:);	% duplicate slips of master-fault

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
newflt = [ x1 y1 x2 y2 z1 z2 dip slips ];

dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns));    smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);        snum = reshape(smat,[],1);

% flt = [dnum snum x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
prjflt = [ dnum snum x1 y1 x2 y2 z1 z2 dip slips ];
[ ~,xyzflt ] = GTdef_prjflt2uni(prjflt);

[ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault2uni(earth,newflt,xyzflt,Xin,Lin,Bin,Nin);

% create smoothing matrices
if strcmp(surf,'free')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(ddip,dlen,Nd,Ns);
elseif strcmp(surf,'fixed')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed_3slips(ddip,dlen,Nd,Ns);
else
    error('GTdef_fault2dif ERROR: surface smoothing is wrong!!!');
end

if strcmp(smooth,'1d3pf')
    sm = sm_1d3pf;
elseif strcmp(smooth,'1d3pb')
    sm = sm_1d3pb;
elseif strcmp(smooth,'2d') 
    sm = sm_2d;
else
    error('GTdef_fault2dif ERROR: smoothing method is wrong!!!');
end

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them

[ modspace ] = GTdef_addall(modspace,Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0);
