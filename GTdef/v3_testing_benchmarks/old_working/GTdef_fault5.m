function [ Xgrn,Bgrn,Ngrn,sm,sm_abs,Aeq,beq,lb,ub,x0 ] = GTdef_fault5(flt,subflt,bndry,Xin,Bin,Nin,mu,nu,smooth,surf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_fault5				  %
% Fault type 5 is discretized into small patches with different local     %
% strike and dip from its master faul ,but they are still touching	  %
% Process one type-5 fault, combine its existing irregular subpatches	  %
% to prepare for the inputs to Matlab function 				  %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					  %
% Here we have no inequalities, so set A=[];b=[]			  %
% Additionally, determine the smoothing matrix for the subfaults	  %
%									  %
% INPUT:					  		  	  %
%  flt = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns dlen slen ]		  %
%    ss  - master-fault strike-slip (left-lateral +)                      %
%    ds  - master-fault dip-slip (thrust +)                               %
%    ts  - master-fault tensile-slip (opening +)                          %
%    ss0,ds0,ts0 - lower bounds for master-fault slips			  %
%    ssX,dsX,tsX - upper bounds for master-fault slips			  %
%    Nd  - number of rows defining the subfaults along dip 	          %
%    Ns  - number of rows defining the subfaults along strike 		  %
%    ddip - average length along dip					  %
%    dlen - average length laong strike					  %
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
% first created by Lujia Feng Fri Dec 11 11:47:47 EST 2009		  %
% changed 'freesurface' to 'surface' lfeng Wed Feb 24 12:50:51 EST 2010	  %
% last modified by Lujia Feng Wed Dec  1 18:58:31 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 13], error('need a 1*13 fault vector for GTdef_fault5'); end

% initialization
% Note: flt is a 1-row vector for the master fault
mslips = flt(1:9);				% slip block
Nd = round(flt(10)); Ns = round(flt(11));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
ddip = flt(1:12); dlen = flt(1:13);
unit = ones(subflt_num,1);
slips = mslips(unit,:);				% duplicate slips of master-fault

% form subfaults
slips = mslips(unit,:);				% duplicate slips of master-fault
if ~isempty(subflt)
    num = size(subflt,1); mat = [Nd Ns];
    for ii = 1:num
        dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
        jj = sub2ind(mat,dnum,snum);
	slips(jj,:) = subflt(ii,3:11);
    end
end

if ~isempty(bndry)
    num = size(bndry,1); mat = [Nd Ns];
    for ii = 1:num
        dd = round(bndry(ii,1)); ss = round(bndry(ii,2));
        jj = sub2ind(mat,dd,ss);
	x1(jj,1) = bndry(ii,3);  y1(jj,1) = bndry(ii,4);
	x2(jj,1) = bndry(ii,12); y2(jj,1) = bndry(ii,13);
	z1(jj,1) = bndry(ii,5);  z2(jj,1) = bndry(ii,8);
	% center points for dip calculation
	xtopc(jj,1) = 0.5*(x1(jj,1)+x2(jj,1)); 
	ytopc(jj,1) = 0.5*(y1(jj,1)+y2(jj,1));
	xbotc(jj,1) = 0.5*(bndry(ii,6)+bndry(ii,9)); 
	ybotc(jj,1) = 0.5*(bndry(ii,7)+bndry(ii,10));
    end
end

dist3d = sqrt((xtopc-xbotc).^2+(ytopc-ybotc).^2+(z1-z2).^2);
heigh = z2-z1;
dip = asind(heigh./dist3d);

num = size(dip,1);
for ii = 1:num
     if isnan(dip(ii))
         x1(ii)=0; y1(ii)=0; x2(ii)=0; y2(ii)=0; dip(ii)=0; slips(ii,:)=0;
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
