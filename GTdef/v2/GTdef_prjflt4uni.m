function [ prjflt,xyzctr,xyztop1 ] = GTdef_prjflt4uni(flt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_prjflt4uni				       %
%									       %
% Project uniform-slip fault4 geometry and slip information                    %
% onto surface geographic coordinate                                           %	
%									       %
% INPUT:					  		  	       %
% flt=[dnum snum x1 y1 x2 y2 z1 z2 dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX] %
%    dnum - row number for subfaults					       %
%    snum - column number for subfaults				  	       %
%    x1,y1 - one endpoint among the two endpoints of the faults                %
%            in the local cartesian coordinate system	  		       %
%    x2,y2 - the other endpoint among the two endpoints of the faults          %
%    z1  - vertical burial depth (top of fault)                                %  
%    z2  - vertical locking depth (bottom of fault)                            %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]          %
%    rake- Aki-Richards convention                                             %
%    rs  - rake-slip (rake direction +)                                        %
%    ts  - tensile-slip (opening +)                                            %
%  rake0,rakeX - rake is usually fixed, currently dummy parameters             %
%      rs0,ts0 - lower bounds for slips				               %
%      rsX,tsX - upper bounds for slips				               %
%                                                                              %
% OUTPUT:                                                                      %
% prjflt = [ dnum snum xtop1 ytop1 xbot1 ybot1 xbot2 ybot2 xtop2 ytop2         %
%             xctr yctr rake rs ts ]  					       %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	       %
%   are the surface projection of four points that confine 		       %
%   the fault interface 						       %
%   They are in a counterclockwise sense looking from the RHS of endpoint      %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	       %
%   [x1 y1] and [x2 y2] are surface projection of fault endpoints	       %
% xyzctr  = [ xx;yy;zz ] (3*nn) for center points                              %
% xyztop1 = [ xx;yy;zz ] (3*nn) for top upper corners                          %
%   zz elevation positive upward                                               %
%									       %
% first created by Lujia Feng Mon May 14 01:27:29 SGT 2012                     %
% output xyzctr & xyztop1 lfeng Tue Jun 12 18:14:04 SGT 2012                   %
% last modified by Lujia Feng Tue Jun 12 18:15:42 SGT 2012                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=18, error('GTdef_prjflt4uni ERROR: need a n*18 fault vector as input!'); end

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

	%  center points
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
        error('GTdef_prjflt4uni ERROR: dip must be within [0 180]!');
    end
else
   error('GTdef_prjflt4uni ERROR: burial depth can not be greater than locking depth!');
end
