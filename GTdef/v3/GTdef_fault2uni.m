function [ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0,Min0 ] = GTdef_fault2uni(flt,Xin,Bin,Nin,earth,mu,nu,edgrn,edgrnfcts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_fault2uni				  %
% Process uniform-slip type-2 faults and                                  %
% prepare for the inputs to Matlab function                               %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)					  %
% Here we have no inequalities, so set A=[];b=[]			  %
%									  %
% INPUT:					  		  	  %
%  flt = [ x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX ]       %
%    x1,y1 - one endpoint among the two endpoints of the faults           %
%            in the local cartesian coordinate system	  		  %
%    x2,y2 - the other endpoint among the two endpoints of the faults     %
%    z1  - vertical burial depth (top of fault) >=0                       %  
%    z2  - vertical locking depth (bottom of fault) >=0                   %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]     %
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
%  (earth = homogeneous)						  %
%  mu  - shear modulus							  %
%  nu  - poisson's ratio                                                  %
%  (earth = layered)		        				  %
%  edgrn, edgrnfcts                                                       %
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
% related function GTdef_fault1uni.m					  %
% first created by Lujia Feng Thu Apr 23 16:14:59 EDT 2009		  %
% preallocate and do parfor lfeng Thu Nov 11 11:56:20 EST 2010	  	  %
% added layered earth model lfeng Tue Feb 28 03:29:36 SGT 2012		  %
% renamed from GTdef_fault2.m to GTdef_fault2uni.m lfeng May 8 SGT 2012   %
% added output Min0 for stress calculation lfeng Thu May 17 SGT 2012      %
% last modified by Lujia Feng Thu May 17 08:44:13 SGT 2012                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=16, error('GTdef_fault2uni ERROR: need a n*16 fault vector as input!'); end

% initialization
[ flt_num ] = size(flt,1);			% fault num
comp_num = 3;                                   % component number = 3 (strike, dip, tensile)
slip_num = flt_num*comp_num;			% slip num
x1 = flt(:,1); y1 = flt(:,2); x2  = flt(:,3); y2 = flt(:,4); 
z1 = flt(:,5); z2 = flt(:,6); dip = flt(:,7);
x0  = [ flt(:,8);  flt(:,9);  flt(:,10) ];	% [ss;ds;ts]
lb  = [ flt(:,11); flt(:,13); flt(:,15) ];	% [ss0;ds0;ts0]
ub  = [ flt(:,12); flt(:,14); flt(:,16) ];	% [ssX;dsX;tsX]

xunit = x0; xunit(xunit~=0) = 1;		% set unit slips for nonzero initial slips
Aeq = zeros(slip_num,slip_num); beq = zeros(slip_num,1);

% fix slips that are not free
for ii = 1:slip_num
   if lb(ii)==ub(ii)			% lb==ub means the slip is fixed
       lb(ii) = -Inf; ub(ii) = Inf;	% relax lb==ub; lsqlin() has a probl with lb==ub
       Aeq(ii,ii) = 1; 		        % use equalities instead of lb==ub
       beq(ii) = x0(ii);
   elseif x0(ii)==0			% initial value is 0, but not fixed
       xunit(ii) = 1;
   end
end

% convert fault info to M matrix taken by disloc3d_mod2.m
xunit_mat = reshape(xunit,[],comp_num);
ss = xunit_mat(:,1); ds = xunit_mat(:,2); ts = xunit_mat(:,3);
if strcmpi(earth,'homogeneous')
   % M = [len width depth dip str east north ss ds ts];	[flt_num*10]
   [ Min0 ] = GTdef_fault2p_to_okada(x1,y1,x2,y2,z1,z2,dip,flt(:,8),flt(:,9),flt(:,10)); 
   % unit slips
   [ M ]    = GTdef_fault2p_to_okada(x1,y1,x2,y2,z1,z2,dip,ss,ds,ts);
   % unit slip for ss only; zero slips for ds,ts
   Mss = M; Mss(:,9:10) = 0;
   % unit slip for ds only; zero slips for ss,ts
   Mds = M; Mds(:,[8 10]) = 0;
   % unit slip for ts only; zero slips for ss,ds
   Mts = M; Mts(:,8:9) = 0;
   % Min = [len;width;depth;dip;str;east;north;ss;ds;ts];	[10*slip_num]
   Min = [ Mss' Mds' Mts' ];
else
   % no opening for EDGRN, so ignore any input for ts
   ts(:,:) = 0;
   xunit = [ ss; ds; ts ];
   % M = [slip north east depth length width str dip rake]
   [ Min0 ] = GTdef_fault2p_to_edcmp(x1,y1,x2,y2,z1,z2,dip,flt(:,8),flt(:,9),ts); 
   % unit slip for ss only; zero slips for ds,ts
   [ Mss ] = GTdef_fault2p_to_edcmp(x1,y1,x2,y2,z1,z2,dip,ss,ts,ts);
   % unit slip for ds only; zero slips for ss,ts
   [ Mds ] = GTdef_fault2p_to_edcmp(x1,y1,x2,y2,z1,z2,dip,ts,ds,ts);
   % zero slips for ss,ds,ts
   Mts = Mds; Mts(:,[1 9]) = 0;
   % Min = [slip;north;east;depth;length;width;str;dip;rake] 	[9*slip_num]
   Min = [ Mss' Mds' Mts' ];
end

% form green matrix
Xgrn = []; Bgrn = []; Ngrn = [];
% point data
if ~isempty(Xin)
    Xgrn = zeros(3*size(Xin,2),slip_num);	% preallocate (fast) to avoid dynamically allocating (slow)
    parfor ii = 1:slip_num			% parfor for parallel computing on multiple workers
        if xunit(ii)~=0				% no need to calculate, if we don't consider this slip
	    if strcmpi(earth,'homogeneous')
            [ U,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Xin,mu,nu);
        else
            % pntsrc = [ pxs pys pzs pmoment ]
            [ pntsrc ]  = GTdef_edcmp_discretize(Min(:,ii),edgrn);
            [ U,~,~,~ ] = GTdef_edcmp(edgrn,[],edgrnfcts,pntsrc,Xin);
	    end
    	    Xgrn(:,ii) = reshape(U',[],1);
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
    Bgrn = zeros(4*size(Bin,2),slip_num);				% preallocate (fast) to avoid dynamically allocating (slow)
    parfor ii = 1:slip_num						% parfor for parallel computing on multiple workers
        if xunit(ii)~=0							% no need to calculate, if we don't consider this slip
	    if strcmpi(earth,'homogeneous')
                [ U1,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Bin1,mu,nu);	
                [ U2,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Bin2,mu,nu);	
	    else
   		[ pntsrc ] = GTdef_edcmp_discretize(Min(:,ii),edgrn);
		[ U1,~,~,~ ] = GTdef_edcmp(edgrn,edgrnfcts,pntsrc,Bin1);
		[ U2,~,~,~ ] = GTdef_edcmp(edgrn,edgrnfcts,pntsrc,Bin2);
	    end
	    dU = U2-U1; dx = dU(1,:); dy = dU(2,:); dz = dU(3,:);
	    dis = sqrt((x2-x1+dx).^2+(y2-y1+dy).^2+(z2-z1+dz).^2);
	    ddis = dis-dis0;
    	    Bgrn(:,ii) = [ reshape(dU',[],1);ddis' ];
        end
    end
end

% profile and grid node data
if ~isempty(Nin)
    Ngrn = zeros(3*size(Nin,2),slip_num);	% preallocate (fast) to avoid dynamically allocating (slow)
    parfor ii = 1:slip_num                      % parfor for parallel computing on multiple workers
        if xunit(ii)~=0				% no need to calculate, if we don't consider this slip
	    if strcmpi(earth,'homogeneous')
                [ U,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Nin,mu,nu);
	    else
   		[ pntsrc ] = GTdef_edcmp_discretize(Min(:,ii),edgrn);
		[ U,~,~,~ ] = GTdef_edcmp(edgrn,edgrnfcts,pntsrc,Nin);
	    end
    	    Ngrn(:,ii) = reshape(U',[],1);
        end
    end
end

% a zero row-vector for smoothing
sm = zeros(1,slip_num);
