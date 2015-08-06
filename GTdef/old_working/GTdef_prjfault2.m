function [ prjflt ] = GTdef_prjfault2(flt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_prjfault2				  %
%									  %
% Project fault2 geometry and slip information onto surface geographic	  %
% coordinate								  %
%									  %
% INPUT:					  		  	  %
% flt=[dnum snum x1 y1 x2 y2 z1 z2 dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]  %
%    dnum - row number for subfaults					  %
%           dnum is 0 for single fault					  %
%    snum - column number for subfaults				  	  %
%           snum is 0 for sinlge fualt					  %
%    x1,y1 - one endpoint among the two endpoints of the faults           %
%            in the local cartesian coordinate system	  		  %
%    x2,y2 - the other endpoint among the two endpoints of the faults     %
%    z1  - vertical burial depth (top of fault)                           %  
%    z2  - vertical locking depth (bottom of fault)                       %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]     %
%    ss  - strike-slip (left-lateral +)                                   %
%    ds  - dip-slip (thrust +)                                            %
%    ts  - tensile-slip (opening +)                                       %
%    ss0,ds0,ts0 - lower bounds for slips				  %
%    ssX,dsX,tsX - upper bounds for slips				  %
%                                                                         %
% OUTPUT:                                                                 %
%  prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2   %
%             xctr  yctr ss ds ts ]  					  %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	  %
%   are the surface projection of four points that confine 		  %
%   the fault interface 						  %
%   They are in a counterclockwise sense looking from the RHS of endpoint %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	  %
%   [x1 y1] and [x2 y2] are surface projection of fault endpoints	  %
%									  %
% first created by Lujia Feng Fri Dec  4 22:45:34 EST 2009		  %
% removed rad2deg & find, use any lfeng Wed Dec  1 15:40:19 EST 2010	  %
% last modified by Lujia Feng Wed Dec  1 16:01:00 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=18, error('need a n*18 fault vector for GTdef_prjfault2'); end

[ flt_num ] = size(flt,1);			% fault num
dnum = flt(:,1); snum = flt(:,2);
x1 = flt(:,3); y1 = flt(:,4); x2 = flt(:,5); y2 = flt(:,6); 
z1 = flt(:,7); z2 = flt(:,8); dip = flt(:,9);
x0  = [ flt(:,10) flt(:,11) flt(:,12) ];		% [ss ds ts]

if z1<z2	% burial depth must be smaller than locking depth
    if min(dip>=0)&&min(dip<=180)
        str = 90-180*atan2(y2-y1,x2-x1)/pi;   % degree CW from N [0-360]
        %%%%% vertical dipping %%%%%
        ind_v = dip==90;
	if any(ind_v)
            xtop1(ind_v,1) = x1(ind_v); ytop1(ind_v,1) = y1(ind_v);
            xbot1(ind_v,1) = x1(ind_v); ybot1(ind_v,1) = y1(ind_v);
            xtop2(ind_v,1) = x2(ind_v); ytop2(ind_v,1) = y2(ind_v);
            xbot2(ind_v,1) = x2(ind_v); ybot2(ind_v,1) = y2(ind_v);
	end
	    
        %%%%% other dipping (looking from the endpoint) %%%%%
        ind = dip~=90;
	if any(ind)
            % points correspond to [x1 y1]
            xtop1(ind,1) = x1(ind)+z1(ind)./tand(dip(ind)).*cosd(str(ind));
            ytop1(ind,1) = y1(ind)-z1(ind)./tand(dip(ind)).*sind(str(ind));
            xbot1(ind,1) = x1(ind)+z2(ind)./tand(dip(ind)).*cosd(str(ind));
            ybot1(ind,1) = y1(ind)-z2(ind)./tand(dip(ind)).*sind(str(ind));
            % points correspond to [x2 y2]
            xtop2(ind,1) = x2(ind)+z1(ind)./tand(dip(ind)).*cosd(str(ind));
            ytop2(ind,1) = y2(ind)-z1(ind)./tand(dip(ind)).*sind(str(ind));
            xbot2(ind,1) = x2(ind)+z2(ind)./tand(dip(ind)).*cosd(str(ind));
            ybot2(ind,1) = y2(ind)-z2(ind)./tand(dip(ind)).*sind(str(ind));
	end

	%  the center point
	xctr(:,1) = 0.5*(xtop1(:,1)+xbot2(:,1));
	yctr(:,1) = 0.5*(ytop1(:,1)+ybot2(:,1));
	%xctr2(:,1) = 0.5*(xtop2(:,1)+xbot1(:,1));
	%yctr2(:,1) = 0.5*(ytop2(:,1)+ybot1(:,1));
	prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr x0 ]; 
    else
        error('Dip must be within [0 180]!');
    end
else
   error('Burial depth can not be greater than locking depth!');
end
