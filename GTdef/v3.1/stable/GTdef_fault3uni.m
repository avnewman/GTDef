function [ xyzflt,Xgrn,Lgrn,Bgrn,Ngrn,sm,Aineq,bineq,Aeq,beq,lb,ub,x0 ] = ...
           GTdef_fault3uni(earth,flt,xyzflt,Xin,Lin,Bin,Nin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                GTdef_fault3uni                                 %
% Process uniform-slip type-3 faults and                                         %
% prepare for inputs to Matlab function                                          %
% x = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0)                                   %
%									         %
% INPUT:					  		  	         %
% earth  structure                                                               %
% xyzflt structure                                                               %
% ------------------------------------------------------------------------------ %
% flt = [ xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]       %
%   xx,yy - one endpoint of the faults in the local cartesian coordinate         %
%   z1  - vertical burial depth (top of fault) >=0                               %
%   z2  - vertical locking depth (bottom of fault) >=0                           %
%   len - fault length                                                           %
%   str - strike from the endpoint (degree CW from N) [0-360]       	         %
%   dip - down from Horiz, right looking from the endpoint [0 180]               %
%   rake- Aki-Richards convention                                                %
%   rs  - rake-slip (rake direction +)                                           %
%   ts  - tensile-slip (opening +)                                               %
% rake0 - lower bounds for rake [rake-90 rake+90]                                %
% rakeX - upper bounds for rake [rake-90 rake+90]                                %
% rs0,ts0 - lower bounds for slips				                 %
% rsX,tsX - upper bounds for slips				                 %
% Xin - point site locations in the local cartesian system 	  	         %
%       [3*n] [ xx yy zz ]						         %
% Lin - los point locations in the local cartesian system + los direction        %
%       [n*6] [xx yy zz dirE dirN dirV]                                          %
% Bin - baseline site locations in the local cartesian system 	  	         %
%       [6*n] [ x1 y1 z1 x2 y2 z2 ]					         %
% Nin - grid and profile node locations in the local cartesian system            %
%       [3*n] [ xx yy zz ]						         %
% Note: z1,z2 depth positive downward                                            %
%       zz elevation positive upward                                             %
%                                                                                %
% OUTPUT:                                                                        %
% slip_num = flt_num*comp_num                                                    %
%    comp_num = 2: two (rake, tensile) components                                %
%    comp_num = 3: three (strike, dip, tensile) components                       %
% ------------------------------------------------------------------------------ %
% Add fields to xyzflt structure                                                 %
%                                          Cin     Min                           %
% SSgrn - strike-slip stress kernels  [6*flt_num*flt_num]                        %
% DSgrn - dip-slip stress kernels     [6*flt_num*flt_num]                        %
% TSgrn - tensile-slip stress kernels [6*flt_num*flt_num]                        %
% Min   - subflt info in cartesian coordinate                                    %
%   for  Okada                                                                   %
%       = [len width depth dip str east north ss ds ts]     [flt_num*10]         %
%   for  layered model                                                           %
%       = [slip north east depth length width str dip rake] [flt_num*9]          %
% compnum - component number                                                     %
% ------------------------------------------------------------------------------ %
% Xgrn - displacements [east;north;vertical] for different sites   	         %
%        from unit slips [(3*nn)*slip_num] 				         %
%        (nn is the  number of sites)  	                                         %
% Lgrn - los displacements [los] for different sites   	                         %
%        from unit slips [(1*nn)*slip_num] 				         %
%        (nn is the number of los points)                                        %
% Bgrn - length changes [east;north;vertical;length] for 	  	         %
%        different baselines from unit slips [(4*nn)*slip_num] 	                 %
%        (nn is the  number of baselines)  				         %
% Ngrn - displacements [east;north;vertical] for different nodes   	         %
%        from unit slips [(3*nn)*slip_num] 				         %
%        (nn is the  number of nodes)                                            %
% Aineq - left-hand  side matrix for linear inequalities  [(flt_num*2)*slip_num] %
% bineq - right-hand side vector for linear inequalities  [flt_num*2]            %
% Aeq   - left-hand  side matrix for linear equalities    [slip_num*slip_num]    %
% beq   - right-hand side vector for linear equalities    [slip_num*1]           %
% x0    - initial values for ss,ds,ts 	[slip_num*1]                             %
% lb    - lower bounds for ss,ds,ts 	[slip_num*1]                             %
% ub    - upper bounds for ss,ds,ts	[slip_num*1]			         %
%                                                                                %
% first created by Lujia Feng Wed May  9 10:45:34 SGT 2012                       %
% added output Min0 for stress calculation lfeng Thu May 17 SGT 2012             %
% added earth structure lfeng Fri Mar 20 11:39:54 SGT 2015                       %
% added Min,SSgrn,DSgrn,TSgrn to xyzflt lfeng Thu Mar 26 15:54:56 SGT 2015       %
% more work for SSgrn,DSgrn,TSgrn lfeng Thu Mar 26 16:51:49 SGT 2015             %
% added los Lin & Lgrn lfeng Tue Nov  3 11:46:52 SGT 2015                        %
% added Aineq & bineq to modspace lfeng Mon Jun 13 14:01:37 SGT 2016             %
% allowed rake0 & rakeX to be 90 deg within rake lfeng Mon Jun 13 SGT 2016       %
% last modified by Lujia Feng Thu Jun 16 23:18:20 SGT 2016                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=16, error('GTdef_fault3uni ERROR: need a n*16 fault vector as input!'); end

etype     = earth.type;
mu        = earth.rigidity;
nu        = earth.poisson;
edgrn     = earth.edgrn;
edgrnfcts = earth.edgrnfcts;

% Cin - center point of patches in the local cartesian system
%     = [n*3] [ xx yy zz ]
Cin = xyzflt.xyzctr;

% initialization
[ flt_num ] = size(flt,1);		   % fault num
xx  = flt(:,1); yy  = flt(:,2); z1  = flt(:,3); z2   = flt(:,4); 
len = flt(:,5); str = flt(:,6); dip = flt(:,7); 
rake  = flt(:,8);
rake0 = flt(:,11);                         % rake lower bounds
rakeX = flt(:,12);                         % rake uppder bounds

% do not invert for rake
if isequal(rake0,rakeX)                    % rake0==rakeX means rake is fixed
    comp_num = 2;                          % component number = 2 (rake, tensile)
    slip_num = flt_num*comp_num;           % slip num

    x0  = [ flt(:,9); flt(:,10) ];         % [rs;ts]
    lb  = [ flt(:,13); flt(:,15) ];        % [rs0;ts0]
    ub  = [ flt(:,14); flt(:,16) ];        % [rsX;tsX]

    xunit = x0; xunit(xunit~=0) = 1;       % set unit slips for nonzero initial slips
    Aineq = zeros(flt_num*2,slip_num); bineq = zeros(flt_num*2,1);
    Aeq   = zeros(slip_num,slip_num);  beq   = zeros(slip_num,1);

    % fix slips that are not free
    for ii = 1:slip_num
        if lb(ii)==ub(ii)                  % lb==ub means the slip is fixed
            lb(ii) = -Inf; ub(ii) = Inf;   % relax lb==ub; lsqlin() has a probl with lb==ub
            Aeq(ii,ii) = 1; 		   % use equalities instead of lb==ub
            beq(ii) = x0(ii);
        elseif x0(ii)==0                   % initial value is 0, but not fixed
            xunit(ii) = 1;
        end
    end
    
    % convert fault info to M matrix taken by disloc3d_mod2.m or GTdef_edcmp.m
    xunit_mat = reshape(xunit,[],comp_num);
    rs = xunit_mat(:,1); ts = xunit_mat(:,2);
    ss = rs.*cosd(rake); ds = rs.*sind(rake);
    ss0 = flt(:,9).*cosd(rake); ds0 = flt(:,9).*sind(rake);
    if strcmpi(etype,'homogeneous')
        % M = [len width depth dip str east north ss ds ts]          [flt_num*10]
        [ Min0 ] = GTdef_fault1p_to_okada(xx,yy,z1,z2,len,str,dip,ss0,ds0,flt(:,10)); 
        % unit slips
        [ M ] = GTdef_fault1p_to_okada(xx,yy,z1,z2,len,str,dip,ss,ds,ts);
        % unit slip for rs only; zero slips for ts
        Mrs = M; Mrs(:,10)  = 0;
        % unit slip for ts only; zero slips for rs
        Mts = M; Mts(:,8:9) = 0;
        % Min = [len;width;depth;dip;str;east;north;ss;ds;ts] 	     [10*slip_num]
        Min = [ Mrs' Mts' ];
    else
        % no opening for EDGRN, so ignore any input for ts
        ts(:,:) = 0;
        xunit = [ ss; ds; ts ];
        % M = [slip north east depth length width str dip rake]
        [ Min0 ] = GTdef_fault1p_to_edcmp(xx,yy,z1,z2,len,str,dip,ss0,ds0,ts);
        % unit slip for rs only; zero slips for ts
        [ Mrs ] = GTdef_fault1p_to_edcmp(xx,yy,z1,z2,len,str,dip,ss,ds,ts);
        % zero slips for ss,ds,ts
        Mts = Mrs; Mts(:,[1 9]) = 0;
        % Min = [slip;north;east;depth;length;width;str;dip;rake]    [9*slip_num]
        Min = [ Mrs' Mts' ];
    end
else
% invert for rake, rotate ss-ds coordinate to rs-rs90 coordinate
%
% rs90      ds    rs
%  .         ^    /           rs90 - 90 anticlockwise relative to rs
%    .       |   /
%      .     |  /
%        .   | /
%          . |/
%  -------------------->  ss
%            |
%            |
%            |
%            |
%
    comp_num = 3;                          % component number = 3 (strike, dip, tensile)
    slip_num = flt_num*comp_num;           % slip num

    % set rs90 the same as rs for now to get correct x0,lb,ub
    x0 = [ flt(:,9);  zeros(flt_num,1);  flt(:,10) ];	% [rs;rs90;ts]     ~ [ss;ds;ts]
    lb = [ flt(:,13); flt(:,13); flt(:,15) ];	% [rs0;rs900;ts0] ~ [ss0;ds0;ts0]
    ub = [ flt(:,14); flt(:,14); flt(:,16) ];	% [rsX;rs90X;tsX] ~ [ssX;dsX;tsX]

    xunit = x0; xunit(xunit~=0) = 1;       % set unit slips for nonzero initial slips
    Aeq   = zeros(slip_num,slip_num); beq   = zeros(slip_num,1);
    % fix slips that are not free
    for ii = 1:slip_num
        if lb(ii)==ub(ii)                  % lb==ub means the slip is fixed
            lb(ii) = -Inf; ub(ii) = Inf;   % relax lb==ub; lsqlin() has a probl with lb==ub
            Aeq(ii,ii) = 1; 		   % use equalities instead of lb==ub
            beq(ii)    = x0(ii);
        elseif x0(ii)==0                   % initial value is 0, but not fixed
            xunit(ii) = 1;
        end
    end
    
    % relax bounds for rs90
    lb(flt_num+1:flt_num*2) =  -Inf(flt_num,1); % set lb for rs90 to be -Inf
    ub(flt_num+1:flt_num*2) =   Inf(flt_num,1); % set ub for rs90 to be  Inf

    % two inequality constraints
    Aineq = zeros(flt_num*2,slip_num); bineq = zeros(flt_num*2,1);
    for ii = 1:flt_num
        % tangent ineqaulity 1 
        % rs90/rs >= tand(rakevary0) -> rs*tand(rakevary0) - rs90 <= 0
        rakevary0 = rake0(ii) - rake(ii);
        if rakevary0 < 0
            if rakevary0 <= -90  
                fprintf(1,'GTdef_fault3uni WARNING: rake0-rake<=-90 is not allowed, reset rake0-rake to be -89.99\n');
                rakevary0 = -89.99;
            end
            Aineq(ii,ii)         = tand(rakevary0);
            Aineq(ii,ii+flt_num) = -1;
            % bineq(ii) = 0;
        else
            fprintf(1,'GTdef_fault3uni WARNING: rake0>=rake is not allowed, reset rake0 to be rake\n');
            Aineq(ii,ii)         =  0; 
            Aineq(ii,ii+flt_num) = -1;
        end
        % tangent ineqaulity 2 
        % rs90/rs <= tand(rakevaryX) -> -rs*tand(rakevaryX) + rs90 <= 0
        rakevaryX = rakeX(ii) - rake(ii);
        if rakevaryX > 0
            if rakevaryX >= 90  
                fprintf(1,'GTdef_fault3uni WARNING: rakeX-rake>=90 is not allowed, reset rakeX-rake to be 89.99\n');
                rakevaryX = 89.99;
            end
            Aineq(ii+flt_num,ii)         = -tand(rakevaryX);
            Aineq(ii+flt_num,ii+flt_num) = 1;
            % bineq(ii+flt_num) = 0;
        else
            fprintf(1,'GTdef_fault3uni WARNING: rakeX<=rake is not allowed, reset rakeX to be rake\n');
            Aineq(ii+flt_num,ii)         = 0;
            Aineq(ii+flt_num,ii+flt_num) = 1;
        end
    end

    % convert fault info to M matrix taken by disloc3d_mod2.m or GTdef_edcmp.m
    xunit_mat = reshape(xunit,[],comp_num);
    rs   = xunit_mat(:,1); 
    rs90 = xunit_mat(:,2); 
    ts   = xunit_mat(:,3);

    ss0 = flt(:,9).*cosd(rake); ds0 = flt(:,9).*sind(rake); % initial values
    if strcmpi(etype,'homogeneous')
        % M = [len width depth dip str east north ss ds ts]     [flt_num*10]
        [ Min0 ] = GTdef_fault1p_to_okada(xx,yy,z1,z2,len,str,dip,ss0,ds0,flt(:,10)); 

        % unit slip for rs,ts
        ss = rs.*cosd(rake); ds = rs.*sind(rake);
        [ M ] = GTdef_fault1p_to_okada(xx,yy,z1,z2,len,str,dip,ss,ds,ts);
        % unit slip for rs only; zero slips for rs90,ts
        Mrs   = M; Mrs(:,10)  = 0;
        % unit slip for ts only; zero slips for rs,rs90
        Mts   = M; Mts(:,8:9) = 0;

        % unit slip for rs90 only; zero slips for rs,ts
        ss = rs90.*cosd(rake+90); ds = rs90.*sind(rake+90); ts(:,:) = 0;
        [ Mrs90 ] = GTdef_fault1p_to_okada(xx,yy,z1,z2,len,str,dip,ss,ds,ts);

        % Min = [len;width;depth;dip;str;east;north;ss;ds;ts] 	[10*slip_num]
        Min = [ Mrs' Mrs90' Mts' ];
    else
        % no opening for EDGRN, so ignore any input for ts
        ts(:,:) = 0; xunit(flt_num*2+1:end) = 0;
        % M = [slip north east depth length width str dip rake]
        [ Min0 ] = GTdef_fault1p_to_edcmp(xx,yy,z1,z2,len,str,dip,ss0,ds0,ts);

        % unit slip for rs only; zero slips for rs90,ts
        ss = rs.*cosd(rake); ds = rs.*sind(rake);
        [ Mrs ] = GTdef_fault1p_to_edcmp(xx,yy,z1,z2,len,str,dip,ss,ds,ts);
        Mts = Mrs; Mts(:,[1 9]) = 0;

        % unit slip for rs90 only; zero slips for rs,ts
        ss = rs90.*cosd(rake+90); ds = rs90.*sind(rake+90);
        [ Mrs90 ] = GTdef_fault1p_to_edcmp(xx,yy,z1,z2,len,str,dip,ss,ds,ts);

        % Min = [slip;north;east;depth;length;width;str;dip;rake]    [9*slip_num]
        Min = [ Mrs' Mrs90' Mts' ];
    end
end



% form stress kernel
SSgrn = [];
DSgrn = [];
TSgrn = [];

% form green matrices
Xgrn = []; Lgrn = []; Bgrn = []; Ngrn = [];

% point data
if ~isempty(Xin)
    Xin  = Xin';
    Xgrn = zeros(3*size(Xin,2),slip_num);	% preallocate (fast) to avoid dynamically allocating (slow)
    parfor ii = 1:slip_num			% parfor for parallel computing on multiple workers
        if xunit(ii)~=0				% no need to calculate, if we don't consider this slip
	    if strcmpi(etype,'homogeneous')
                [ U,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Xin,mu,nu);	
            else
                % pntsrc = [ pxs pys pzs pmoment ]
                [ pntsrc ] = GTdef_edcmp_discretize(Min(:,ii),edgrn);
                [ U,~,~,~ ] = GTdef_edcmp(edgrn,[],edgrnfcts,pntsrc,Xin);
	    end
    	    Xgrn(:,ii) = reshape(U',[],1);
        end
    end
end

% InSAR los data
if ~isempty(Lin)
    Ldir = Lin(:,4:end);                        % LOS direction from ground to satellite [east north vert]
    Lin  = Lin(:,1:3)';                         % LOS point locations
    Lgrn = zeros(1*size(Lin,2),slip_num);	% preallocate (fast) to avoid dynamically allocating (slow)
    parfor ii = 1:slip_num			% parfor for parallel computing on multiple workers
        if xunit(ii)~=0				% no need to calculate, if we don't consider this slip
	    if strcmpi(etype,'homogeneous')
                [ U,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Lin,mu,nu);	
            else
                % pntsrc = [ pxs pys pzs pmoment ]
                [ pntsrc ] = GTdef_edcmp_discretize(Min(:,ii),edgrn);
                [ U,~,~,~ ] = GTdef_edcmp(edgrn,[],edgrnfcts,pntsrc,Lin);
	    end
            % project deformation in 3D to LOS using dot product
            %        [ dirE dirN dirV ]              
            % Ldir = |  :    :    :   |      [lospnt_num*3]
            %        [  :    :    :   ]              
            % ------------------------------------------------
            %     [ Ue ... ]    [ Ux ... ]
            % U = | Un ... | or | Uy ... |   [3*lospnt_num]
            %     [ Uu ... ]    [ Uz ... ]
            % ------------------------------------------------
            ULOS       = dot(Ldir,U',2);        % 2 means along row dim
    	    Lgrn(:,ii) = ULOS;
        end
    end
end

% baseline data
if ~isempty(Bin)
    Bin = Bin';
    x1 = Bin(1,:);  y1 = Bin(2,:);  z1 = Bin(3,:);			% 1th sites of the baselines
    x2 = Bin(4,:);  y2 = Bin(5,:);  z2 = Bin(6,:);                     	% 2nd sites of the baselines
    z0 = zeros(1,length(x1));
    Bin1 = [x1;y1;z0];  Bin2 = [x2;y2;z0];				% set depth=0 for inputting to disloc3d_mod2 
    dis0 = sqrt((x2-x1).^2+(y2-y1).^2+(z2-z1).^2);			% base distance between two sites
    Bgrn = zeros(4*size(Bin,2),slip_num);				% preallocate (fast) to avoid dynamically allocating (slow)
    parfor ii = 1:slip_num						% parfor for parallel computing on multiple workers
        if xunit(ii)~=0							% no need to calculate, if we don't consider this slip
	    if strcmpi(etype,'homogeneous')
                [ U1,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Bin1,mu,nu);	
                [ U2,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Bin2,mu,nu);	
	    else
   		[ pntsrc ] = GTdef_edcmp_discretize(Min(:,ii),edgrn);
		[ U1,~,~,~ ] = GTdef_edcmp(edgrn,[],edgrnfcts,pntsrc,Bin1);
		[ U2,~,~,~ ] = GTdef_edcmp(edgrn,[],edgrnfcts,pntsrc,Bin2);
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
    Nin  = Nin';
    Ngrn = zeros(3*size(Nin,2),slip_num);	% preallocate (fast) to avoid dynamically allocating (slow)
    parfor ii = 1:slip_num                      % parfor for parallel computing on multiple workers
        if xunit(ii)~=0				% no need to calculate, if we don't consider this slip
	    if strcmpi(etype,'homogeneous')
                [ U,~,~,~,~,~ ] = disloc3d_mod2(Min(:,ii),Nin,mu,nu);
	    else
   		[ pntsrc ] = GTdef_edcmp_discretize(Min(:,ii),edgrn);
		[ U,~,~,~ ] = GTdef_edcmp(edgrn,[],edgrnfcts,pntsrc,Nin);
	    end
    	    Ngrn(:,ii) = reshape(U',[],1);
        end
    end
end

% a zero row-vector for smoothing
sm = zeros(1,slip_num);

xyzflt.Min     = Min0;
xyzflt.SSgrn   = SSgrn; 
xyzflt.DSgrn   = DSgrn;
xyzflt.TSgrn   = TSgrn;
xyzflt.compnum = comp_num;
