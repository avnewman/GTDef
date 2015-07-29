function [ prjflt ] = GTdef_prjfault1_EW(flt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GTdef_prjfault1_EW				  %
%									  %
% Project fault1 geometry and slip information onto surface geographic	  %
% coordinate ASSUMING FAULT STRIKE EW					  %
% Note: this is used for GMT grdrotater, not real fault orientation	  %
%                                                             		  %
% INPUT:					  		  	  %
% flt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]%
%    dnum - row number for subfaults					  %
%           dnum is 0 for single fault					  %
%    snum - column number for subfaults				  	  %
%           snum is 0 for sinlge fualt					  %
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
%                                                                         %
% OUTPUT:                                                                 %
%  prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2   %
%             xctr  yctr ss ds ts ]  					  %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	  %
%   are the surface projection of four points that confine 		  %
%   the fault interface 						  %
%   They are in a counterclockwise sense looking from the RHS of endpoint %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	  %
%   [xx yy] is surface projection of one fault endpoint 	  	  %
%                                                                         %
% based on GTdef_prjfault1.m lfeng Tue Jul 27 02:21:25 EDT 2010		  %
% removed disp find, use any lfeng Wed Dec  1 16:07:09 EST 2010		  %
% last modified by lfeng Wed Dec  1 16:07:17 EST 2010			  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=18, error('need a n*18 fault vector for GTdef_prjfault1_EW');  end

[ flt_num ] = size(flt,1);			% fault num
dnum = flt(:,1); snum = flt(:,2);
xx = flt(:,3); yy = flt(:,4); z1 = flt(:,5); z2 = flt(:,6); 
len = flt(:,7); str = flt(:,8); dip = flt(:,9);

str(:) = 90; 

x0 = [ flt(:,10) flt(:,11) flt(:,12) ];		% [ss ds ts]

if z1<z2	% burial depth must be smaller than locking depth
    if min(dip>=0)&&min(dip<=180)
        %%%%% vertical dipping %%%%%
        ind_v = dip==90;
	if any(ind_v)
            xtop1(ind_v,1) = xx(ind_v); ytop1(ind_v,1) = yy(ind_v);
            xbot1(ind_v,1) = xx(ind_v); ybot1(ind_v,1) = yy(ind_v);
            xtop2(ind_v,1) = xx(ind_v)+len(ind_v).*sind(str(ind_v)); 
	    ytop2(ind_v,1) = yy(ind_v)+len(ind_v).*cosd(str(ind_v));
            xbot2(ind_v,1) = xtop2(ind_v,1); ybot2(ind_v,1) = xtop2(ind_v,1);
	end
	    
        %%%%% other dipping (looking from the endpoint) %%%%%
        ind = dip~=90;
	if any(ind)
            % points correspond to [xx yy]
            xtop1(ind,1) = xx(ind)+z1(ind)./tand(dip(ind)).*cosd(str(ind));
            ytop1(ind,1) = yy(ind)-z1(ind)./tand(dip(ind)).*sind(str(ind));
            xbot1(ind,1) = xx(ind)+z2(ind)./tand(dip(ind)).*cosd(str(ind));
            ybot1(ind,1) = yy(ind)-z2(ind)./tand(dip(ind)).*sind(str(ind));
	    % [xx2 yy2] is the other endpoint of fault projection at the surface
            xx2(ind,1) = xx(ind)+len(ind).*sind(str(ind)); 
	    yy2(ind,1) = yy(ind)+len(ind).*cosd(str(ind));
            xtop2(ind,1) = xx2(ind)+z1(ind)./tand(dip(ind)).*cosd(str(ind));
            ytop2(ind,1) = yy2(ind)-z1(ind)./tand(dip(ind)).*sind(str(ind));
            xbot2(ind,1) = xx2(ind)+z2(ind)./tand(dip(ind)).*cosd(str(ind));
            ybot2(ind,1) = yy2(ind)-z2(ind)./tand(dip(ind)).*sind(str(ind));
	end

	%  the center point
	xctr(:,1) = 0.5*(xtop1(:,1)+xbot2(:,1));
	yctr(:,1) = 0.5*(ytop1(:,1)+ybot2(:,1));
	%xctr2(:,1) = 0.5*(xtop2(:,1)+xbot1(:,1));
	%yctr2(:,1) = 0.5*(ytop2(:,1)+ybot1(:,1));
	prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr x0 ]; 
    else
        error('dip must be within [0 180]!');
    end
else
   error('Burial depth is greater than locking depth!');
end
