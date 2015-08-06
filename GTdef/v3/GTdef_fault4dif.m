function [ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0,Min0 ] = GTdef_fault4dif(flt,subflt,dipin,Xin,Bin,Nin,...
								      earth,mu,nu,edgrn,edgrnfcts,smooth,surf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_fault4dif				       %
% Process one distributed-slip type-4 fault, generate its subfaults and        %
% call GTdef_fault4uni to prepare for the inputs to Matlab function 	       %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					       %
% Here we have no inequalities, so set A=[];b=[]			       %
% Additionally, determine the smoothing matrix for the subfaults	       %
%									       %
% INPUT:					  		  	       %
% flt = [ x1 y1 x2 y2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX Nd Ns ] %
%    x1,y1 - one endpoint among the two endpoints of the master fault          %
%    x2,y2 - the other endpoint among the two endpoints 		       %
%            both in the local cartesian coordinate system	               %
%    z1  - vertical burial depth (top of fault) >=0                            %  
%    z2  - vertical locking depth (bottom of fault) >=0                        %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]          %
%   rake - Aki-Richards convention                                             %
%    rs  - rake-slip (rake direction +)                                        %
%    ts  - tensile-slip (opening +)                                            %
%  rake0,rakeX - rake is usually fixed, currently dummy parameters             %
%      rs0,ts0 - lower bounds for slips				               %
%      rsX,tsX - upper bounds for slips				               %
%    Nd  - number of rows defining the subfaults along dip 	               %
%    Ns  - number of rows defining the subfaults along strike 		       %
% subflt = [ dnum snum rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]                %
%   dnum - row number for subfaults					       %
%   snum - column number for subfaults				  	       %
%   rake - rake direction                                                      %
%  rs,ts - subfault slips						       %
%      rake0,rakeX - rake is usually fixed, currently dummy parameters         %
%  rs0,ts0,rsX,tsX - subfault slip bounds			               %
%  Xin - point site locations in the local cartesian system 	  	       %
%        [3*n] [ xx;yy;zz ]						       %
%  Bin - baseline site locations in the local cartesian system 	  	       %
%        [6*n] [ x1;y1;z1;x2;y2;z2 ]					       %
%  Nin - grid and profile node locations in the local cartesian system         %
%        [3*n] [ xx;yy;zz ]						       %
%  (earth = homogeneous)						       %
%  mu  - shear modulus							       %
%  nu  - poisson's ratio                                                       %
%  (earth = layered)		        				       %
%  edgrn, edgrnfcts                                                            %
%  smooth - smoothing method						       %
%  surf   - surface smoothing setting					       %
%                                                                              %
% OUTPUT:                                                                      %
% Each fault has two (rake and tensile) components, so                         %
%  slip_num = subflt_num*2 = Nd*Ns*2                                           %
%  Xgrn - displacements [east;north;vertical] for different sites   	       %
%         from unit slips [(3*nn)*slip_num] 				       %
%         (nn is the  number of sites)  				       %
%  Bgrn - length changes [east;north;vertical;length] for 	  	       %
%         different baselines from unit slips [(4*nn)*slip_num] 	       %
%         (nn is the  number of baselines)  				       %
%  Ngrn - displacements [east;north;vertical] for different nodes   	       %
%         from unit slips [(3*nn)*slip_num] 				       %
%         (nn is the  number of nodes)  				       %
%  Aeq - left-hand side matrix for linear equalities  [slip_num*slip_num]      %
%  beq - right-hand side vector for linear equalities [slip_num*1]             %
%  x0  - initial values for ss,ds,ts 	[slip_num*1]                           %
%  lb  - lower bounds for ss,ds,ts 	[slip_num*1]                           %
%  ub  - upper bounds for ss,ds,ts	[slip_num*1]			       %
%  sm  - smoothing matrix for slips     [slip_num*slip_num]		       %
%  sm_abs - matrix for calculating the absolute 1st derivative		       %
%                                                                              %
% first created by Lujia Feng Thu May 10 15:04:55 SGT 2012                     %
% added output Min0 for stress calculation lfeng Thu May 17 SGT 2012           %
% added sm = [] for subflt_num==1 lfeng Tue Oct 21 17:25:33 SGT 2014           %
% last modified by Lujia Feng Tue Oct 21 17:34:57 SGT 2014                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], error('GTdef_fault4dif ERROR: need a 1*18 fault vector as input!'); end

% initialization
% Note: flt is a 1-row vector for the master fault
mx1 = flt(1); my1 = flt(2); mx2 = flt(3); my2 = flt(4); 
mz1 = flt(5); mz2 = flt(6); mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
mstr = GTdef_strike(mx1,my1,mx2,my2);
unit = ones(subflt_num,1);

% if uniform slip == only one patch
if subflt_num==1
    [ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0,Min0 ] = GTdef_fault4uni([flt(:,1:end-2)],Xin,Bin,Nin,earth,mu,nu,edgrn,edgrnfcts);
    sm_abs = [];
    sm     = [];
    return
end

% exclusively specify dips for rows between z1 and z2 depth range
if ~isempty(dipin)
    [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,mdip,Nd,Ns,dipin);
else
    dip = mdip*unit; 			% duplicate dips
    zlin = linspace(mz1,mz2,Nd+1)';
    xlin = linspace(mx1,mx2,Ns+1); ylin = linspace(my1,my2,Ns+1); 
    x1mat = xlin(ones(Nd,1),1:end-1); y1mat = ylin(ones(Nd,1),1:end-1);	% endpoint 1
    x2mat = xlin(ones(Nd,1),2:end);   y2mat = ylin(ones(Nd,1),2:end);	% endpoint 2
    z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));	% depths
    x1 = reshape(x1mat,[],1);  y1 = reshape(y1mat,[],1);
    x2 = reshape(x2mat,[],1);  y2 = reshape(y2mat,[],1);
    z1 = reshape(z1mat,[],1);  z2 = reshape(z2mat,[],1);
    dz = (mz2-mz1)/Nd; ddip = dz/sind(mdip);
end

% form subfaults
slips = mslips(unit,:);				% duplicate slips of master-fault
dlen = sqrt((mx1-mx2)^2+(my1-my2)^2)/Ns;

%  subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]
if ~isempty(subflt)
    num = size(subflt,1); mat = [Nd Ns];
    for ii = 1:num
        dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
        jj = sub2ind(mat,dnum,snum);
	slips(jj,:) = subflt(ii,3:11);
    end
end
newflt = [ x1 y1 x2 y2 z1 z2 dip slips ];

% treat subfault sequence as fault type 2
[ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0,Min0 ] = GTdef_fault4uni(newflt,Xin,Bin,Nin,earth,mu,nu,edgrn,edgrnfcts);


% create smoothing matrices
if strcmp(surf,'free')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free_2slips(ddip,dlen,Nd,Ns);
elseif strcmp(surf,'fixed')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed_2slips(ddip,dlen,Nd,Ns);
else
    error('GTdef_fault4dif ERROR: Surface smoothing is wrong!!!');
end

if strcmp(smooth,'1d3pf')
    sm = sm_1d3pf;
elseif strcmp(smooth,'1d3pb')
    sm = sm_1d3pb;
elseif strcmp(smooth,'2d') 
    sm = sm_2d;
else
    error('GTdef_fault4dif ERROR: Smoothing method is wrong!!!');
end

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them
