function [ prjflt,xyzflt ] = GTdef_prjflt3uni(flt)

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
% Note: z1,z2 depth positive downward                                              %
%                                                                                  %
% OUTPUT:                                                                          %
% prjflt = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1                         %
%            xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr                    %
%            ss ds ts rake rs ]                                                    %
%   [xtop1 ytop1], [xbot1 ybot1], [xbot2 ybot2], and [xtop2 ytop2]	           %
%   are the surface projection of four points that confine 		           %
%   the fault interface 						           %
%   They are in a counterclockwise sense looking from the RHS of endpoint          %
%   [xtop1 ytop1] and [xbot1 ybot1] correspond [xx yy] at z1 and z2	           %
%   [xx yy] is surface projection of one fault endpoint 	  	           %
% -------------------------------------------------------------------------------- %
% flt?.xyzflt stucture for fault patch info in cartesian coordinate                %
% xyzflt.xyzctr  - location of center points in cartesian coordinate               %
%                = [ xx yy zz ]                    (flt_num*3)                     %
% xyzflt.xyztop1 - location of top upper corners in cartesian coordinate           %
%                = [ xx yy zz ]                    (flt_num*3)                     %
% xyzflt.suv     - unit vector in the strike direction                             %
%                = [ xx yy zz ]                    (flt_num*3)                     %
% xyzflt.duv     - unit vector in the dip direction                                %
%                = [ xx yy zz ]                    (flt_num*3)                     %
% xyzflt.nuv     - normal unit vector                                              %
%                = [ xx yy zz ]                    (flt_num*3)                     %
% Note: zz elevation positive upward                                               %
%                                                                                  %
% first created by Lujia Feng Mon May 14 01:11:31 SGT 2012                         %
% output xyzctr & xyztop1 lfeng Tue Jun 12 18:06:23 SGT 2012                       %
% output depths lfeng Wed Nov 12 19:39:45 SGT 2014                                 %
% added xyzflt lfeng Tue Mar 24 11:21:45 SGT 2015                                  %
% added ss & ds to prjflt lfeng Wed Apr 27 18:00:37 SGT 2016                       %
% added es & ns to prjflt lfeng Thu Jun  2 16:38:40 SGT 2016                       %
% last modified by Lujia Feng Wed Apr 27 18:16:57 SGT 2016                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(flt,2)~=18  
    error('GTdef_prjflt3uni ERROR: need a n*18 fault vector as input!'); 
end

[ fltnum ] = size(flt,1);			% fault num
dnum = flt(:,1); 
snum = flt(:,2);
xx   = flt(:,3); 
yy   = flt(:,4); 
z1   = flt(:,5); 
z2   = flt(:,6); 
len  = flt(:,7); 
str  = flt(:,8); 
dip  = flt(:,9);
% rake + rs (rake slip) + ts
rake = flt(:,10);
rs   = flt(:,11);
ts   = flt(:,12);
% ss + ds
ss = rs.*cosd(rake);
ds = rs.*sind(rake);
% slip in east + north directions
[ es,ns ] = rotate_xy(ss,ds,str-360-90);

x0 = [ ss ds ts rake rs es ns ];

if z1<z2	% burial depth must be smaller than locking depth
    if min(dip>=0)&&min(dip<=180)
        %%%%% vertical dipping %%%%%
        ind_v = dip==90;
	if any(ind_v)
            xtop1(ind_v,1) = xx(ind_v); ytop1(ind_v,1) = yy(ind_v);   ztop1(ind_v,1) = z1(ind_v);
            xbot1(ind_v,1) = xx(ind_v); ybot1(ind_v,1) = yy(ind_v);   zbot1(ind_v,1) = z2(ind_v);
            xtop2(ind_v,1) = xx(ind_v)+len(ind_v).*sind(str(ind_v));  ztop2(ind_v,1) = z1(ind_v);
	    ytop2(ind_v,1) = yy(ind_v)+len(ind_v).*cosd(str(ind_v));
            xbot2(ind_v,1) = xtop2(ind_v,1); ybot2(ind_v,1) = xtop2(ind_v,1); zbot2(ind_v,1) = z2(ind_v);
	end
	    
        %%%%% other dipping (looking from the endpoint) %%%%%
        ind = dip~=90;
	if any(ind)
            % points correspond to [xx yy]
            xtop1(ind,1) = xx(ind)+z1(ind)./tand(dip(ind)).*cosd(str(ind));
            ytop1(ind,1) = yy(ind)-z1(ind)./tand(dip(ind)).*sind(str(ind));
            ztop1(ind,1) = z1(ind);
            xbot1(ind,1) = xx(ind)+z2(ind)./tand(dip(ind)).*cosd(str(ind));
            ybot1(ind,1) = yy(ind)-z2(ind)./tand(dip(ind)).*sind(str(ind));
            zbot1(ind,1) = z2(ind);
	    % [xx2 yy2] is the other endpoint of fault projection at the surface
            xx2(ind,1) = xx(ind)+len(ind).*sind(str(ind)); 
	    yy2(ind,1) = yy(ind)+len(ind).*cosd(str(ind));
            xtop2(ind,1) = xx2(ind)+z1(ind)./tand(dip(ind)).*cosd(str(ind));
            ytop2(ind,1) = yy2(ind)-z1(ind)./tand(dip(ind)).*sind(str(ind));
            ztop2(ind,1) = z1(ind);
            xbot2(ind,1) = xx2(ind)+z2(ind)./tand(dip(ind)).*cosd(str(ind));
            ybot2(ind,1) = yy2(ind)-z2(ind)./tand(dip(ind)).*sind(str(ind));
            zbot2(ind,1) = z2(ind);
	end

        % center points
        xctr(:,1) = 0.5*(xtop1(:,1)+xbot2(:,1));
        yctr(:,1) = 0.5*(ytop1(:,1)+ybot2(:,1));
        zctr(:,1) = 0.5*(z1(:,1)+z2(:,1));
        xyzctr  = [xctr yctr -zctr];

        % top upper corners
        xyztop1 = [xtop1 ytop1 -z1];
        %xctr2(:,1) = 0.5*(xtop2(:,1)+xbot1(:,1));
        %yctr2(:,1) = 0.5*(ytop2(:,1)+ybot1(:,1));

        % unit vector in the strike direction
        suv = [ sind(str) cosd(str) zeros(fltnum,1) ];

        % unit vector in the dip direction (thrust positive)
        duv = [ -cosd(dip).*cosd(str) cosd(dip).*sind(str) sind(dip) ];

        % unit vector in the normal direction
        nuv = [ sind(dip).*cosd(str) -sind(dip).*sind(str) cosd(dip) ];
        
        prjflt = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1 xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr x0 ];
    else
        error('GTdef_prjflt3uni ERROR: dip must be within [0 180]!');
    end
else
    error('GTdef_prjflt3uni ERROR: burial depth is greater than locking depth!');
end

xyzflt.xyzctr  = xyzctr;
xyzflt.xyztop1 = xyztop1;
xyzflt.suv     = suv;
xyzflt.duv     = suv;
xyzflt.nuv     = nuv;
