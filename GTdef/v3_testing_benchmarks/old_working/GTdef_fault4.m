function [ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_fault4(flt,subflt,dipin,Xin,Bin,Nin,mu,nu,smooth,surf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_fault4				  %
% Process one type-4 fault, generate its subfaults and, call GTdef_fault2 %
% to prepare for the inputs to Matlab function 				  %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					  %
% Here we have no inequalities, so set A=[];b=[]			  %
% Additionally, determine the smoothing matrix for the subfaults	  %
%									  %
% INPUT:					  		  	  %
%  flt = [ x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns ] %
%    x1,y1 - one endpoint among the two endpoints of the master fault     %
%    x2,y2 - the other endpoint among the two endpoints 		  %
%            both in the local cartesian coordinate system	          %
%    z1  - vertical burial depth (top of fault)                           %  
%    z2  - vertical locking depth (bottom of fault)                       %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]     %
%    ss  - master-fault strike-slip (left-lateral +)                      %
%    ds  - master-fault dip-slip (thrust +)                               %
%    ts  - master-fault tensile-slip (opening +)                          %
%    ss0,ds0,ts0 - lower bounds for master-fault slips			  %
%    ssX,dsX,tsX - upper bounds for master-fault slips			  %
%    Nd  - number of rows defining the subfaults along dip 	          %
%    Ns  - number of rows defining the subfaults along strike 		  %
%  subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  %
%    dnum - row number for subfaults					  %
%    snum - column number for subfaults				  	  %
%    ss,ds,ts - subfault slips						  %
%    ss0,ds0,ts0,ssX,dsX,tsX - subfault slip bounds			  %
%  Xin - point site locations in the local cartesian system 	  	  %
%        [3*n] [ xx;yy;zz ]						  %
%  Bin - baseline site locations in the local cartesian system 	  	  %
%        [6*n] [ x1;y1;z1;x2;y2;z2 ]					  %
%  Nin - grid and profile node locations in the local cartesian system    %
%        [3*n] [ xx;yy;zz ]						  %
%  mu  - shear modulus							  %
%  nu  - poisson's ratio                                                  %
%  smooth - smoothing method						  %
%  surf   - surface smoothing setting					  %
%                                                                         %
% OUTPUT:                                                                 %
% Each subfault has three (strike, dip, and tensile) components, so       %
%  slip_num = subflt_num*3 = Nd*Ns*3                                      %
%  Xgrn - displacements [east;north;vertical] for different sites   	  %
%         from unit slips [(3*nn)*slip_num] 				  %
%         (nn is the  number of sites)  				  %
%  Bgrn - length changes [east;north;vertical;length] for 	  	  %
%         different baselines from unit slips [(4*nn)*slip_num] 	  %
%         (nn is the  number of baselines)  				  %
%  Ngrn - displacements [east;north;vertical] for different nodes   	  %
%         from unit slips [(3*nn)*slip_num] 				  %
%         (nn is the  number of nodes)  				  %
%  Aeq - left-hand side matrix for linear equalities  [slip_num*slip_num] %
%  beq - right-hand side vector for linear equalities [slip_num*1]        %
%  x0  - initial values for ss,ds,ts 	[slip_num*1]                      %
%  lb  - lower bounds for ss,ds,ts 	[slip_num*1]                      %
%  ub  - upper bounds for ss,ds,ts	[slip_num*1]			  %
%  sm  - smoothing matrix for slips     [slip_num*slip_num]		  %
%  sm_abs - matrix for calculating the absolute 1st derivative		  %
%                                                                         %
% first created by Lujia Feng Fri Apr 24 10:44:57 EDT 2009		  %
% added "smooth" by Lujia Feng Wed Dec  2 02:39:38 EST 2009		  %
% added auxiliary sm by lfeng Thu Dec  3 00:14:33 EST 2009		  %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		  %
% changed 'freesurface' to 'surface' lfeng Wed Feb 24 13:16:24 EST 2010	  %
% last modified by Lujia Feng Wed Dec  1 18:57:58 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], error('need a 1*18 fault vector for GTdef_fault4'); end

% initialization
% Note: flt is a 1-row vector for the master fault
mx1 = flt(1); my1 = flt(2); mx2 = flt(3); my2 = flt(4); 
mz1 = flt(5); mz2 = flt(6); mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
mstr = GTdef_strike(mx1,my1,mx2,my2);
unit = ones(subflt_num,1);

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
newflt = [ x1 y1 x2 y2 z1 z2 dip slips];

% treat subfault sequence as fault type 2
[ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_fault2(newflt,Xin,Bin,Nin,mu,nu);


% create smoothing matrices
if strcmp(surf,'free')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_free(ddip,dlen,Nd,Ns);
elseif strcmp(surf,'fixed')
    [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed(ddip,dlen,Nd,Ns);
else
    error('Surface smoothing is wrong!!!');
end

if strcmp(smooth,'1d3pf')
    sm = sm_1d3pf;
elseif strcmp(smooth,'1d3pb')
    sm = sm_1d3pb;
elseif strcmp(smooth,'2d') 
    sm = sm_2d;
else
    error('Smoothing method is wrong!!!');
end

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them
