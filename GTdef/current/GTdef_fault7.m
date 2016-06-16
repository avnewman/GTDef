function [ Xgrn,Lgrn,Bgrn,Ngrn,sm,sm_abs,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = GTdef_fault7(flt,subflt,vertices,grnfns,smooth,surf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 GTdef_fault7				         %
% Fault type 7 takes external Greens functions from FEM tools e.g. PyLith        %
% for irregular surface discretized into small patches                           %
% with different local strike and dip                                            %
% x = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0)                                   %
% Additionally, determine the smoothing matrix for the subfaults	         %
%									         %
% INPUT:					  		  	         %
%  flt = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns ]		                 %
%    ss  - master-fault strike-slip (left-lateral +)                             %
%    ds  - master-fault dip-slip (thrust +)                                      %
%    ts  - master-fault tensile-slip (opening +)                                 %
%    ss0,ds0,ts0 - lower bounds for master-fault slips			         %
%    ssX,dsX,tsX - upper bounds for master-fault slips			         %
%    Nd  - number of rows defining the subfaults along dip 	                 %
%    Ns  - number of rows defining the subfaults along strike 		         %
%    ddip - average length along dip					         %
%    dlen - average length laong strike					         %
%  subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		         %
%    dnum - row number for subfaults					         %
%    snum - column number for subfaults				  	         %
%    ss,ds,ts - subfault slips						         %
%    ss0,ds0,ts0,ssX,dsX,tsX - subfault slip bounds			         %
%vertices - subfault/vertices location                                           %
%         = [ id dnum snum xx yy zz ]                                            %
% grnfns  - array of size                        (vertexNum*pntNum*9)            %
%         for each patch-site pair                                               % 
%         = [ ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]            %
%  smooth - smoothing method						         %
%  surf   - surface smoothing setting					         %
%                                                                                %
% OUTPUT:                                                                        %
% Each subfault has three (strike, dip, and tensile) components, so              %
%  slip_num = subflt_num*3 = Nd*Ns*3                                             %
%  Xgrn - displacements [east;north;vertical] for different sites   	         %
%         from unit slips [(3*nn)*slip_num] 				         %
%         (nn is the  number of sites)  				         %
%  Aineq  - left-hand  side matrix for linear inequalities [(flt_num*2)*slip_num]%
%  bineq  - right-hand side vector for linear inequalities [flt_num*2]           %
%  Aeq    - left-hand side matrix for linear equalities    [slip_num*slip_num]   %
%  beq    - right-hand side vector for linear equalities   [slip_num*1]          %
%  x0     - initial values for ss,ds,ts    [slip_num*1]                          %
%  lb     - lower bounds for ss,ds,ts      [slip_num*1]                          %
%  ub     - upper bounds for ss,ds,ts      [slip_num*1]			         %
%  sm     - smoothing matrix for slips     [slip_num*slip_num]		         %
%  sm_abs - matrix for calculating the absolute 1st derivative		         %
%                                                                                %
% first created by Lujia Feng Fri Dec 11 11:47:47 EST 2009		         %
% changed from GTdef_fault5 to GTdef_fault6 lfeng Wed Jun 17 SGT 2015            %
% last modified by Lujia Feng Wed Jun 17 11:04:57 SGT 2015                       %
% renamed from GTdef_fault6 to GTdef_fault7 lfeng Apr 20 11:59:26 SGT 2016       %
% added Lgrn = [] lfeng Fri Jun 10 01:01:28 SGT 2016                             %
% added Aineq & bineq to modspace lfeng Mon Jun 13 17:05:53 SGT 2016             %
% need to code for los Lgrn lfeng Tue Nov  3 14:32:13 SGT 2015                   %
% last modified by Lujia Feng Thu Jun 16 00:57:37 SGT 2016                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 11], error('GTdef_fault7 ERROR: need a 1*11 fault vector as input!'); end

% initialization
Xgrn = []; Lgrn = []; Bgrn = []; Ngrn = [];

% Note: flt is a 1-row vector for the master fault
mslips = flt(1:9);				% slip block
Nd = round(flt(10)); Ns = round(flt(11));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
unit = ones(subflt_num,1);
slips = mslips(unit,:);				% duplicate slips of master-fault

% form subfaults
if ~isempty(subflt)
    num = size(subflt,1); mat = [Nd Ns];
    for ii = 1:num
        dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
        jj = sub2ind(mat,dnum,snum);
	slips(jj,:) = subflt(ii,3:11);
    end
end

% check if GTdef input Nd & Ns == Greens functions Nd & Ns
vertNd = max(vertices(:,2));
vertNs = max(vertices(:,3));
if vertNd ~= Nd
    error('GTdef_fault7 ERROR: GTdef input Nd is inconsistent with greens function database Nd!');
end
if vertNs ~= Ns
    error('GTdef_fault7 ERROR: GTdef input Ns is inconsistent with greens function database Ns!');
end

% initialization
x0 = [ slips(:,1); slips(:,2); slips(:,3) ];	% [ss;ds;ts]
lb = [ slips(:,4); slips(:,6); slips(:,8) ];	% [ss0;ds0;ts0]
ub = [ slips(:,5); slips(:,7); slips(:,9) ];	% [ssX;dsX;tsX]

comp_num = 3;                                   % component number = 3 (strike, dip, tensile)
slip_num = subflt_num*comp_num;			% slip num
Aineq = zeros(subflt_num*2,slip_num); bineq = zeros(subflt_num*2,1);  % useful only for fault type 3 & 4
Aeq   = zeros(slip_num,slip_num);     beq   = zeros(slip_num,1);

% fix slips that are not free
for ii = 1:slip_num
    if lb(ii)==ub(ii)			        % lb==ub means the slip is fixed
    	lb(ii) = -Inf; ub(ii) = Inf;	        % relax lb==ub; lsqlin() has a probl with lb==ub
	Aeq(ii,ii) = 1; 		        % use equalities instead of lb==ub
	beq(ii) = x0(ii);
    end
end

% form Xgrn
if ~isempty(grnfns)
    % strike-slip
    cgrn = grnfns(:,:,1:3);
    cgrn = permute(cgrn,[2 3 1]);
    cgrn = reshape(cgrn,[],subflt_num);
    Xgrn = [ Xgrn cgrn ];
    % dip-slip
    cgrn = grnfns(:,:,4:6);
    cgrn = permute(cgrn,[2 3 1]);
    cgrn = reshape(cgrn,[],subflt_num);
    Xgrn = [ Xgrn cgrn ];
    % opening
    cgrn = grnfns(:,:,7:9);
    cgrn = permute(cgrn,[2 3 1]);
    cgrn = reshape(cgrn,[],subflt_num);
    Xgrn = [ Xgrn cgrn ];
end

% create smoothing matrices
[ ddip,dlen ] = PyLith_dist_greensfns(vertices);
if strcmp(surf,'free')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_3slips(ddip,dlen,Nd,Ns);
elseif strcmp(surf,'fixed')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed_3slips(ddip,dlen,Nd,Ns);
else
    error('GTdef_fault7 ERROR: surface smoothing is wrong!!!');
end

if strcmp(smooth,'1d3pf')
    sm = sm_1d3pf;
elseif strcmp(smooth,'1d3pb')
    sm = sm_1d3pb;
elseif strcmp(smooth,'2d') 
    sm = sm_2d;
else
    error('GTdef_fault7 ERROR: smoothing method is wrong!!!');
end

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them
