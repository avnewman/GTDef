function [ newflt,prjflt,xyzflt ] = GTdef_prjflt2dif(flt,subflt,dipin,strin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_prjflt2dif				  %
%									  %
% Project distributed-slip fault2 geometry and slip information           %
% onto surface geographic coordinate				          %
%									  %
% INPUT:					  		  	  %
% flt = [ x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns ]  %
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
% Note: z1,z2 depth positive downward                                     %
% subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  %
%    dnum - row number for subfaults					  %
%    snum - column number for subfaults				  	  %
%    ss,ds,ts - subfault slips						  %
%    ss0,ds0,ts0,ssX,dsX,tsX - subfault slip bounds			  %
% dipin - dip addon info for the master fault                             %
%       = [ dip z1 z2 rows ]                                              %
% strin - strike addon info for the master fault                          %
%       = [ x1 y1 x2 y2 columns sweepAngle ]                              %
%                                                                         %
% OUTPUT:                                                                 %
% newflt - all fault patches                                              %
%        = [ dnum snum x1 y1 x2 y2 z1 z2 dip slips ]                      %
% prjflt = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1                %
%            xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr           %
%            ss ds ts rake rs es ns ]                                     %
%   [xtop1 ytop1],[xbot1 ybot1],[xbot2 ybot2],and [xtop2 ytop2]	are the   %
%   surface projections of four points confining the fault interface      %
%   They are in a counterclockwise sense looking from the RHS of endpoint %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	  %
%   [xx yy] is surface projection of one fault endpoint 	  	  %
% ----------------------------------------------------------------------- %
% flt?.xyzflt stucture for fault patch info in cartesian coordinate       %
% xyzflt.xyzctr  - location of center points in cartesian coordinate      %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.xyztop1 - location of top upper corners in cartesian coordinate  %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.suv     - unit vector in the strike direction                    %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.duv     - unit vector in the dip direction                       %
%                = [ xx yy zz ]                    (flt_num*3)            %
% xyzflt.nuv     - normal unit vector                                     %
%                = [ xx yy zz ]                    (flt_num*3)            %
% Note: zz elevation positive upward                                      %
%									  %
% first created by Lujia Feng Fri Sat Dec  5 00:59:25 EST 2009		  %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		  %
% renamed from GTdef_prjfault4.m to GTdef_prjflt2dif.m lfeng May 14 2012  %
% output xyzctr & xyztop1 lfeng Tue Jun 12 18:19:51 SGT 2012              %
% output newflt lfeng Wed Jun 13 15:55:40 SGT 2012                        %
% added strin lfeng Fri Oct 24 18:42:43 SGT 2014                          %
% added sweepAngle lfeng Wed Nov 12 10:57:14 SGT 2014                     %
% added xyzflt lfeng Tue Mar 24 11:21:45 SGT 2015                         %
% last modified by Lujia Feng Tue Mar 24 13:04:58 SGT 2015                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], error('GTdef_prjflt2dif ERROR: need a 1*18 fault vector as input!');  end

% initialization
% Note: flt is a 1-row vector for the master fault
mx1 = flt(1); my1 = flt(2); mx2 = flt(3); my2 = flt(4); 
mz1 = flt(5); mz2 = flt(6); mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); 
Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
mstr = GTdef_strike(mx1,my1,mx2,my2);
unit = ones(subflt_num,1);

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
    x1=[]; y1=[]; x2=[]; y2=[]; z1=[]; z2=[]; dip=[];
    colsum = sum(strin(:,5)); % can be several mini-segments for each segment
    if colsum ~= Ns, error('GTdef_prjflt2dif ERROR: strike is not specified correctly!'); end
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
        dip = [dip;cdip];
    end
end

% form subfaults
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

dlin = round(linspace(1,Nd,Nd)'); 
slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns)); 
smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);  
snum = reshape(smat,[],1);
newflt = [ dnum snum x1 y1 x2 y2 z1 z2 dip slips ];

[ prjflt,xyzflt ] = GTdef_prjflt2uni(newflt);
