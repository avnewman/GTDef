function [ M ] = GTdef_fault1p_to_edcmp(xx,yy,z1,z2,len,str,dip,ss,ds,ts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       GTdef_fault1p_to_edcmp                            		%
% Convert the fault info described by one endpoint to the format taken    		%
% by Fortran code EDCMP							  		%
% Note: coordinates are different							%
% GTdef: X = east,  Y = north, Z = downward						%
% EDCMP: X = north, Y = east,  Z = downward	                                        %
% Note: north & east are different							%
%                                                                                       %
% INPUT											%
% GTdef 1-point fault input:                                                		%
%      xx,yy - surface projection of one endpoint of the fault        			%
%              but EDCMP needs xx yy depth of the reference point			%
%      z1  - vertical burial depth (top of fault) >=0                         		%  
%      z2  - vertical locking depth (bottom of fault) >=0                     		%
%      len - fault length                                                 		%
%      str - strike from the endpoint (degree CW from N) [0-360]          		%
%      dip - down from Horiz, right looking from the endpoint [0 180]     		%
%      ss  - strike slip (left-lateral +)                                 		%
%      ds  - dip-slip (thrust +)                                          		%
%      ts  - tensile-slip (opening +)                                     		%
% Note: they can be scalars or column vectors; the endpoint can be either 		%
%       the left one or the right one					  		%
%---------------------------------------------------------------------------------------%
% OUTPUT										%
%        M = [slip north east depth length width str dip rake] [slipnum*9]		%
% EDCMP format output is                                      				%
%        M = [slip;north;east;depth;length;width;str;dip;rake] [9*slipnum]		%
%    [xs ys zs] is the upper reference point for strike					%
%   length - strike direction								%
%    width - dip direction								%
%     slip - dislocation value [m]							%
%   strike - CW from North [deg]							%
%      dip - downward from horizontal [deg] [0 90]					%
%     rake - CCW from strike [deg] [-180 180]						%
%                   N                                                                   %
%                  /                                                                    %
%                 /| strike                                                             %
%         Ref:-> @------------------------                                              %
%                |\        p .            \ W                                           %
%                :-\      i .              \ i                                          %
%                |  \    l .                \ d                                         %
%                :90 \  S .                  \ t                                        %
%                |-dip\  .                    \ h					%
%                :     \. | rake               \                                        %
%                Z      -------------------------                                       %
%                              L e n g t h                                              %
%    Note that if one of the parameters length and width = 0, then a line source        %
%    will be considered and the dislocation parameter Slip has the unit m^2; if         %
%    both length and width = 0, then a point source will be considered and the          %
%    Slip has the unit m^3.                                                             %
% Note: opening is not considered in EDCMP						%
%---------------------------------------------------------------------------------------%
%											%
% REFERENCE  										%
% Wang, R., Martin, F. L., & Roth, F. (2003)						%
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		%
%                                                                         		%
% related function: GTdef_fault2p_to_edcmp()                              		%
% first created by Mon Feb 27 11:32:02 SGT 2012				  		%
% last modified by Wed Feb 29 12:39:12 SGT 2012						%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if z1<z2	% burial depth must be smaller than locking depth
    if min(dip>=0)&&min(dip<=180)
        width = (z2-z1)./abs(sind(dip));
        depth = z1;				% z>=0
	slip  = sqrt(ss.^2+ds.^2);		% ts is ignored
	rake  = atan2(ds,ss).*180/pi;		% P = atan2(Y,X)
        %%%%% right dipping (looking from the endpoint) %%%%%
        ind_r = dip<=90;
        if any(ind_r);
	    north(ind_r,1) = yy(ind_r); 
	    east(ind_r,1)  = xx(ind_r);		% switch x y
	end
        %%%%% left dipping (looking from the endpoint) %%%%%
        ind_l = dip>90;
        if any(ind_l);
            east(ind_l,1)  = xx(ind_l)+len(ind_l).*sind(str(ind_l));		
            north(ind_l,1) = yy(ind_l)+len(ind_l).*cosd(str(ind_l));
	    % change dip to [0 90]
	    dip(ind_l,1) = 180 - dip(ind_l,1);
	    % reverse strike
	    str(ind_l,1) = str(ind_l)+180;
	    ind = str>=360;
	    str(ind,1)   = str(ind)-360;
	end
	% convert xx yy at surface to at burial depth
	ind_d = depth~=0;
	if any(ind_d)
            east(ind_d,1)  = east(ind_d)+depth(ind_d)./tand(dip(ind_d)).*cosd(str(ind_d));
            north(ind_d,1) = north(ind_d)-depth(ind_d)./tand(dip(ind_d)).*sind(str(ind_d));
	end
	% in EDGRN: slip - dislocation; north - xs; east - ys; depth - zs
        M = [ slip north east depth len width str dip rake ];
    else
        error('GTdef_fault1p_to_edcmp ERROR: dip must be within [0 180]!');
    end
else
   error('GTdef_fault1p_to_edcmp ERROR: burial depth can not be greater than locking depth!');
end
