function [ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_fault3(flt,subflt,Xin,Bin,Nin,mu,nu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_fault3				  %
% Process one type-3 fault, generate its subfaults, and call GTdef_fault1 %
% to prepare for the inputs to Matlab function 				  %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					  %
% Here we have no inequalities, so set A=[];b=[]			  %
% Additionally, determine the smoothing matrix for the subfaults	  %
%									  %
% INPUT:					  		  	  %
%  flt = [ x1 y1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns]%
%    x1,y1 - one endpoint among the two endpoints of the master fault     %
%            in the local cartesian coordinate system	  		  %
%    z1  - vertical burial depth (top of fault)                           %  
%    z2  - vertical locking depth (bottom of fault)                       %
%    len - fault length                                                   %
%    str - strike from the endpoint (degree CW from N) [0-360]            %
%    dip - down from Horiz, right looking from the endpoint [0 180]       %
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
%                                                                         %
% related function GTdef_fault4.m					  %
% first created by Lujia Feng Thu Apr 23 21:46:12 EDT 2009		  %
% last modified by Lujia Feng Mon May 11 13:39:54 EDT 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], disp('need a 1*18 fault vector for GTdef_fault3'); return; end

% initialization
% Note: flt is a row vector for the master fault
x1 = flt(1); y1 = flt(2); mz1 = flt(3); mz2 = flt(4); 
mlen = flt(5);  mstr = flt(6);  mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num

% form subfaults
x2 = x1+mlen*sind(mstr);     y2 = y1+mlen*cosd(mstr);	% endpoint 2
xlin = linspace(x1,x2,Ns+1); ylin = linspace(y1,y2,Ns+1); 
zlin = linspace(mz1,mz2,Nd+1)';
xmat = xlin(ones(Nd,1),1:end-1);  ymat = ylin(ones(Nd,1),1:end-1);
z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));
xx = reshape(xmat,[],1);  yy = reshape(ymat,[],1);
z1 = reshape(z1mat,[],1); z2 = reshape(z2mat,[],1);
dlen = mlen/Ns;	dz = (mz2-mz1)/Nd; ddip = dz/sind(mdip);
unit = ones(subflt_num,1);
len = dlen*unit; str = mstr*unit; dip = mdip*unit;
slips = mslips(unit,:);				% duplicate slips of master-fault

%  subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]
if ~isempty(subflt)
    num = size(subflt,1); mat = [Nd Ns];
    for ii = 1:num
        dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
        jj = sub2ind(mat,dnum,snum);
	slips(jj,:) = subflt(ii,3:11);
    end
end
newflt = [ xx yy z1 z2 len str dip slips];

% treat subfault sequence as fault type 1
[ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_fault1(newflt,Xin,Bin,Nin,mu,nu);

% create smoothing matrix for one type of slip
[ sm0 ] = GTdef_smooth2d(ddip,dlen,Nd,Ns);	       % [ sm0  0  0  ]
% create smoothing matrix for three types of slips  sm = [  0  sm0 0  ]
[ sm ] = GTdef_add_diagonal(sm0,sm0);		       % [  0  0  sm0 ]
[ sm ] = GTdef_add_diagonal(sm,sm0);

ind_fixed = find(lb==-Inf);	% index for fixed slips 
sm(ind_fixed) = 0;		% don't do smoothing for them
