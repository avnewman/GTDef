function [ newflt,prjflt,xyzflt ] = GTdef_prjflt1dif(flt,subflt,dipin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GTdef_prjflt1dif				  %
%									  %
% Project distributed-slip fault1 geometry and slip information           %
% onto surface geographic coordinate                                      %
%									  %
% INPUT:					  		  	  %
% flt = [ x1 y1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX Nd Ns ]%
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
% Note: z1,z2 depth positive downward                                     %
%  subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]		  %
%    dnum - row number for subfaults					  %
%    snum - column number for subfaults				  	  %
%    ss,ds,ts - subfault slips						  %
%    ss0,ds0,ts0,ssX,dsX,tsX - subfault slip bounds			  %
% dipin - dip addon info for the master fault                             %
%       = [ dip z1 z2 rows ]                                              %
%                                                                         %
% OUTPUT:                                                                 %
% newflt - all fault patches                                              %
%        = [ dnum snum x1 y1 z1 z2 len str dip slips ]                    %
% prjflt = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1                %
%            xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr           %
%            ss ds ts rake ts es ns ]                                     %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	  %
%   are the surface projection of four points that confine 		  %
%   the fault interface 						  %
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
% first created by Lujia Feng Fri Dec  4 23:54:23 EST 2009		  %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		  %
% renamed from GTdef_prjfault3.m to GTdef_prjflt1dif.m lfeng May 13 2012  %
% output xyzctr & xyztop1 lfeng Tue Jun 12 18:19:51 SGT 2012              %
% output newflt lfeng Wed Jun 13 15:53:43 SGT 2012                        %
% corrected GTdef_diffdips error lfeng Mon Mar 23 11:19:57 SGT 2015       %
% added xyzflt lfeng Tue Mar 24 11:21:45 SGT 2015                         %
% last modified by Lujia Feng Tue Mar 24 13:05:32 SGT 2015                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt)~=[1 18], error('GTdef_prjflt1dif ERROR: need a 1*18 fault vector as input!'); end

% initialization
% Note: flt is a row vector for the master fault
mx1 = flt(1); my1 = flt(2); mz1 = flt(3); mz2 = flt(4); 
mlen = flt(5);  mstr = flt(6);  mdip = flt(7);
mslips = flt(8:16);				% slip block
Nd = round(flt(17)); Ns = round(flt(18));	% number of rows and columns
subflt_num = Nd*Ns;				% subfault num
unit = ones(subflt_num,1);
mx2 = mx1+mlen*sind(mstr);     my2 = my1+mlen*cosd(mstr);	% endpoint 2

% exclusively specify dips for rows between z1 and z2 depth range
if ~isempty(dipin)
    [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,dipin,Nd,Ns);
else
    xlin = linspace(mx1,mx2,Ns+1); 
    ylin = linspace(my1,my2,Ns+1); 
    zlin = linspace(mz1,mz2,Nd+1)';
    x1mat = xlin(ones(Nd,1),1:end-1);  y1mat = ylin(ones(Nd,1),1:end-1);
    z1mat = zlin(1:end-1,ones(1,Ns));  z2mat = zlin(2:end,ones(1,Ns));
    z1 = reshape(z1mat,[],1); z2 = reshape(z2mat,[],1);
    x1 = reshape(x1mat,[],1); y1 = reshape(y1mat,[],1);
    dip = mdip*unit; 			% duplicate dips
end

dlen = mlen/Ns; len = dlen*unit; str = mstr*unit;
slips = mslips(unit,:);				% duplicate slips of master-fault

% subflt = [ dnum snum ss ds ts ss0 ssX ds0 dsX ts0 tsX ]
if ~isempty(subflt)
    num = size(subflt,1); mat = [Nd Ns];
    for ii = 1:num
        dnum = round(subflt(ii,1)); snum = round(subflt(ii,2));
        jj = sub2ind(mat,dnum,snum);
	slips(jj,:) = subflt(ii,3:11);
    end
end

% create dnum & snum
dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
dmat = dlin(1:end,ones(1,Ns)); smat = slin(ones(Nd,1),1:end);
dnum = reshape(dmat,[],1);  snum = reshape(smat,[],1);

% newflt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
newflt = [ dnum snum x1 y1 z1 z2 len str dip slips ];
[ prjflt,xyzflt ] = GTdef_prjflt1uni(newflt);
