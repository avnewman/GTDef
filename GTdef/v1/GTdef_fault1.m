function [ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ] = GTdef_fault1(flt,Xin,Bin,Nin,mu,nu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_fault1				  %
% Process the type-1 faults and prepare for the inputs to Matlab function %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					  %
% Here we have no inequalities, so set A=[];b=[]			  %
%									  %
% INPUT:					  		  	  %
%  flt = [ xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX ]     %
%    xx,yy - one endpoint among the two endpoints of the faults           %
%            in the local cartesian coordinate system	  		  %
%    z1  - vertical burial depth (top of fault)                           %  
%    z2  - vertical locking depth (bottom of fault)                       %
%    len - fault length                                                   %
%    str - strike from the endpoint (degree CW from N) [0-360]       	  %
%    dip - down from Horiz, right looking from the endpoint [0 180]       %
%    ss  - strike-slip (left-lateral +)                                   %
%    ds  - dip-slip (thrust +)                                            %
%    ts  - tensile-slip (opening +)                                       %
%    ss0,ds0,ts0 - lower bounds for slips				  %
%    ssX,dsX,tsX - upper bounds for slips				  %
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
% Each fault has three (strike, dip, and tensile) components, so          %
%  slip_num = flt_num*3                                                   %
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
%                                                                         %
% related function GTdef_fault2.m					  %
% first created by Lujia Feng Wed Apr 22 16:52:43 EDT 2009	   	  %
% last modified by Lujia Feng Fri May  8 02:46:05 EDT 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=16, disp('need a n*16 fault vector for GTdef_fault1'); return; end

% initialization
[ flt_num ] = size(flt,1);			% fault num
slip_num = flt_num*3;				% slip num
xx = flt(:,1); yy = flt(:,2); z1 = flt(:,3); z2 = flt(:,4); 
len = flt(:,5); str = flt(:,6); dip = flt(:,7);
x0  = [ flt(:,8); flt(:,9); flt(:,10) ];	% [ss;ds;ts]
lb  = [ flt(:,11); flt(:,13); flt(:,15) ];	% [ss0;ds0;ts0]
ub  = [ flt(:,12); flt(:,14); flt(:,16) ];	% [ssX;dsX;tsX]

xunit = x0; xunit(xunit~=0) = 1;		% set unit slips for nonzero initial slips
Aeq = zeros(slip_num,slip_num); beq = zeros(slip_num,1);

% fix slips that are not free
for ii = 1:slip_num
    if lb(ii)==ub(ii)			% lb==ub means the slip is fixed
    	lb(ii) = -Inf; ub(ii) = Inf;	% relax lb==ub; lsqlin() has a probl with lb==ub
	Aeq(ii,ii) = 1; 		% use equalities instead of lb==ub
	beq(ii) = x0(ii);
    elseif x0(ii)==0			% initial value is 0, but not fixed
	xunit(ii) = 1;
    end
end

% convert fault info to M matrix taken by disloc3d_mod2.m
xunit3 = reshape(xunit,[],3);
ss = xunit3(:,1); ds = xunit3(:,2); ts = xunit3(:,3);
% M = [len width depth dip str east north ss ds ts];	[flt_num*10]
[ M ] = GTdef_fault1p_to_okada(xx,yy,z1,z2,len,str,dip,ss,ds,ts);
% unit slip for ss only; zero slips for ds,ts
Mss = M; Mss(:,9:10) = 0;
% unit slip for ds only; zero slips for ss,ts
Mds = M; Mds(:,[8 10]) = 0;
% unit slip for ss only; zero slips for ds,ts
Mts = M; Mts(:,8:9) = 0;
% Min = [len;width;depth;dip;str;east;north;ss;ds;ts];	[10*slip_num]
Min = [ Mss' Mds' Mts' ];

% form green matrices
Xgrn = []; Bgrn = []; Ngrn = [];
% point data
if ~isempty(Xin)
    zero_col = zeros(3*size(Xin,2),1);
    for ii = 1:slip_num
        if xunit(ii)~=0
            [U,D,S,FLAG,FLAG2] = disloc3d_mod2(Min(:,ii),Xin,mu,nu);
    	    Xgrn(:,ii) = reshape(U',[],1);
        else
            Xgrn(:,ii) = zero_col;		% no need to calculate, if we don't consider this slip
        end
    end
end
% baseline data
if ~isempty(Bin)
    zero_col = zeros(4*size(Bin,2),1);
    x1 = Bin(1,:);  y1 = Bin(2,:);  z1 = Bin(3,:);			% 1th sites of the baselines
    x2 = Bin(4,:);  y2 = Bin(5,:);  z2 = Bin(6,:);                     	% 2nd sites of the baselines
    z0 = zeros(1,length(x1));
    Bin1 = [x1;y1;z0];  Bin2 = [x2;y2;z0];				% set depth=0 for inputting to disloc3d_mod2 
    dis0 = sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);			% base distance between two sites
    for ii = 1:slip_num
        if xunit(ii)~=0
            [U1,D,S,FLAG,FLAG2] = disloc3d_mod2(Min(:,ii),Bin1,mu,nu);	
            [U2,D,S,FLAG,FLAG2] = disloc3d_mod2(Min(:,ii),Bin2,mu,nu);	
	    dU = U2-U1; dx = dU(1,:); dy = dU(2,:); dz = dU(3,:);
	    dis = sqrt((x2-x1+dx).^2+(y2-y1+dy).^2+(z2-z1+dz).^2);
	    ddis = dis-dis0;
    	    Bgrn(:,ii) = [ reshape(dU',[],1);ddis' ];
        else
            Bgrn(:,ii) = zero_col;		% no need to calculate, if we don't consider this slip
        end
    end
end
% profile and grid node data
if ~isempty(Nin)
    zero_col = zeros(3*size(Nin,2),1);
    for ii = 1:slip_num
        if xunit(ii)~=0
            [U,D,S,FLAG,FLAG2] = disloc3d_mod2(Min(:,ii),Nin,mu,nu);
    	    Ngrn(:,ii) = reshape(U',[],1);
        else
            Ngrn(:,ii) = zero_col;		% no need to calculate, if we don't consider this slip
        end
    end
end

% a zero row-vector for smoothing
sm = zeros(1,slip_num);
