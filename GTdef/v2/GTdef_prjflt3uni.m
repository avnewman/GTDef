function [ prjflt,xyzctr,xyztop1 ] = GTdef_prjflt3uni(flt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GTdef_prjflt3uni				   %
%									           %
% Project uniform-slip fault3 geometry and slip information                        %
% onto surface geographic coordinate				                   %
%                                                             		           %
% INPUT:					  		  	           %
% flt = [dnum snum xx yy z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX] %
%    dnum - row number for subfaults					           %
%    snum - column number for subfaults				  	           %
%    xx,yy - one endpoint among the two endpoints of the faults                    %
%            in the local cartesian coordinate system	  		           %
%    z1  - vertical burial depth (top of fault)                                    %  
%    z2  - vertical locking depth (bottom of fault)                                %
%    len - fault length                                                            %
%    str - strike from the endpoint (degree CW from N) [0-360]       	           %
%    dip - down from Horiz, right looking from the endpoint [0 180]                %
%   rake - Aki-Richards convention                                                 %
%    rs  - rake-slip (rake direction +)                                            %
%    ts  - tensile-slip (opening +)                                                %
%  rake0,rakeX - rake is usually fixed, currently dummy parameters                 %
%      rs0,ts0 - lower bounds for slips				                   %
%      rsX,tsX - upper bounds for slips				                   %
%                                                                                  %
% OUTPUT:                                                                          %
% prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2             %
%             xctr yctr rake rs ts ]  					           %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	           %
%   are the surface projection of four points that confine 		           %
%   the fault interface 						           %
%   They are in a counterclockwise sense looking from the RHS of endpoint          %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	           %
%   [xx yy] is surface projection of one fault endpoint 	  	           %
% xyzctr  = [ xx;yy;zz ] (3*nn) for center points                                  %
% xyztop1 = [ xx;yy;zz ] (3*nn) for top upper corners                              %
%   zz elevation positive upward                                                   %
%                                                                                  %
% first created by Lujia Feng Mon May 14 01:11:31 SGT 2012                         %
% output xyzctr & xyztop1 lfeng Tue Jun 12 18:06:23 SGT 2012                       %
% last modified by Lujia Feng Tue Jun 12 18:10:26 SGT 2012                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=18  
    error('GTdef_prjflt3uni ERROR: need a n*18 fault vector as input!'); 
end

[ flt_num ] = size(flt,1);			% fault num
dnum = flt(:,1); snum = flt(:,2);
xx = flt(:,3); yy = flt(:,4); z1 = flt(:,5); z2 = flt(:,6); 
len = flt(:,7); str = flt(:,8); dip = flt(:,9);
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

	% center points
	xctr(:,1) = 0.5*(xtop1(:,1)+xbot2(:,1));
	yctr(:,1) = 0.5*(ytop1(:,1)+ybot2(:,1));
        zctr(:,1) = -0.5*(z1(:,1)+z2(:,1));
	xyzctr  = [xctr';yctr';zctr'];
	% top upper corners
	xyztop1 = [xtop1';ytop1';-z1'];
	%xctr2(:,1) = 0.5*(xtop2(:,1)+xbot1(:,1));
	%yctr2(:,1) = 0.5*(ytop2(:,1)+ybot1(:,1));

	prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2 xctr yctr x0 ]; 
    else
        error('GTdef_prjflt3uni ERROR: dip must be within [0 180]!');
    end
else
    error('GTdef_prjflt3uni ERROR: burial depth is greater than locking depth!');
end
