function [ M ] = GTdef_fault2p_to_okada(x1,y1,x2,y2,z1,z2,dip,ss,ds,ts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       GTdef_fault2p_to_okada                            %
% Convert the fault info described by two endpoints to the format         %
% taken by Okada's Fortran code disloc3d_mod2.m that Ting Chen @ GT       %
% converted to Matlab                                                     %
% (1) 2-point fault input:                                                %
%      x1,y1 - endpoint 1; x2,y2 - endpoint 2 		                  %
%      z1  - vertical burial depth (top of fault)                         %  
%      z2  - vertical locking depth (bottom of fault)                     %
%      dip - down from Horiz, right looking from the endpoint 1 [0 180]   %
%      ss  - strike slip (left-lateral +)                                 %
%      ds  - dip-slip (thrust +)                                          %
%      ts  - tensile-slip (opening +)                                     %
% Note: they can be scalars or column vectors; x1,y1 can be either        %
%       the left endpoint or the right endpoint				  %
% (2) disloc3d_mod2 format output in M:                                   %
%      M = [length;width;depth;dip;str;east;north;ss;ds;ts]               %
%      [east north depth] is the bottom center for dip > 0,               %
%      and is the top center for dip < 0                                  %
% Note: the dip Okada's code takes is [-90 90]                            %
%                                                                         %
% based on Peter Cervelli's Fortran code (Ting Chen converted to Matlab)  %
% related function: GTdef_fault2p_to_okada()                              %
% first created by Lujia Feng Apr 21, 2009                                %
% added more comments about the endpoints Mon May 18 14:54:29 EDT 2009	  %
% removed rad2deg find disp, use any lfeng Wed Dec 1 15:48:34 EST 2010	  %
% last modified by Lujia Feng Wed Dec  1 15:58:17 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if z1<z2	% burial depth must be smaller than locking depth
    if min(dip>=0)&&min(dip<=180)
        width = (z2-z1)./abs(sind(dip));
        len = sqrt((x2-x1).^2+(y2-y1).^2);
        str = 90-180*atan2(y2-y1,x2-x1)/pi;   % degree CW from N [0-360]
        %%%%% right dipping (looking from endpoint 1) %%%%%
        ind_r = dip<90;
        % use the bottom center
	% (cx,cy) is the anchor point at the locking depth corresponding to endpoint 1
        if any(ind_r);
            cx(ind_r,1) = x1(ind_r)+z2(ind_r)./tand(dip(ind_r)).*cosd(str(ind_r)); 	% "." can be ignored here actually
            cy(ind_r,1) = y1(ind_r)-z2(ind_r)./tand(dip(ind_r)).*sind(str(ind_r));
            depth(ind_r,1) = z2(ind_r);
	end
        
        %%%%% vertical dipping %%%%%
        ind_v = dip==90;
        % use the bottom center
	% (cx,cy) is the anchor point at the locking depth corresponding to endpoint 1
        if any(ind_v);
            cx(ind_v,1) = x1(ind_v);
            cy(ind_v,1) = y1(ind_v);
            depth(ind_v,1) = z2(ind_v);
	end
	    
        %%%%% left dipping (looking from endpoint 1) %%%%%
        ind_l = dip>90;
        dip(ind_l) = dip(ind_l) - 180;      	%  dip [degree] conversion
        % use the top center
	% (cx,cy) is the anchor point at the burial depth corresponding to endpoint 1
        if any(ind_l);
            cx(ind_l,1) = x1(ind_l)+z1(ind_l)./tand(dip(ind_l)).*cosd(str(ind_l));
            cy(ind_l,1) = y1(ind_l)-z1(ind_l)./tand(dip(ind_l)).*sind(str(ind_l));
            depth(ind_l,1) = z1(ind_l);
            ss(ind_l) = -1*ss(ind_l);           % opening is still opening
            ds(ind_l) = -1*ds(ind_l);
	end
	    
        east  = cx+0.5.*len.*sind(str);		
        north = cy+0.5.*len.*cosd(str);
        M = [ len width depth dip str east north ss ds ts ];
    else
        error('Dip must be within [0 180]!');
    end
else
   error('Burial depth can not be greater than locking depth!');
end
