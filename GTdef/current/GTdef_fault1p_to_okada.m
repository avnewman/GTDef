function [ M ] = GTdef_fault1p_to_okada(xx,yy,z1,z2,len,str,dip,ss,ds,ts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       GTdef_fault1p_to_okada                            %
% Convert the fault info described by one endpoint to the format taken    %
% by Okada's Fortran code disloc3d_mod2.m that Ting Chen @ GT       	  %
% converted to Matlab                                                     %
% (1) 1-point fault input:                                                %
%      xx,yy - one endpoint of the fault        			  %
%      z1  - vertical burial depth (top of fault)                         %  
%      z2  - vertical locking depth (bottom of fault)                     %
%      len - fault length                                                 %
%      str - strike from the endpoint (degree CW from N) [0-360]          %
%      dip - down from Horiz, right looking from the endpoint [0 180]     %
%      ss  - strike slip (left-lateral +)                                 %
%      ds  - dip-slip (thrust +)                                          %
%      ts  - tensile-slip (opening +)                                     %
% Note: they can be scalars or column vectors; the endpoint can be either %
%       the left one or the right one					  %
% (2) disloc3d_mod2 format output in M:                                   %
%      M = [length;width;depth;dip;str;east;north;ss;ds;ts]               %
%      [east north depth] is the bottom center for dip > 0,               %
%      and is the top center for dip < 0                                  %
% Note: the dip Okada's code takes is [-90 90]                            %
%                                                                         %
% based on Peter Cervelli's Fortran code (Ting Chen converted to Matlab)  %
% related function: GTdef_fault2p_to_okada()                              %
% first created by Lujia Feng Apr 21, 2009                                %
% added more comments about the endpoint Mon May 18 14:17:44 EDT 2009	  %
% removed find and use any lfeng Wed Dec  1 15:55:19 EST 2010		  %
% last modified by Lujia Feng Wed Dec  1 15:58:47 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if z1<z2	% burial depth must be smaller than locking depth
    if min(dip>=0)&&min(dip<=180)
        width = (z2-z1)./abs(sind(dip));
        %%%%% right dipping (looking from the endpoint) %%%%%
        ind_r = dip<90;
        % use the bottom center
	% (cx,cy) is the anchor point at the locking depth corresponding to the endpoint
        if any(ind_r);
            cx(ind_r,1) = xx(ind_r)+z2(ind_r)./tand(dip(ind_r)).*cosd(str(ind_r));	% "." can be ignored here actually
            cy(ind_r,1) = yy(ind_r)-z2(ind_r)./tand(dip(ind_r)).*sind(str(ind_r));
            depth(ind_r,1) = z2(ind_r);
	end
        
        %%%%% vertical dipping %%%%%
        ind_v = dip==90;
        % use the bottom center
	% (cx,cy) is the anchor point at the locking depth corresponding to the endpoint
        if any(ind_v);
            cx(ind_v,1) = xx(ind_v);
            cy(ind_v,1) = yy(ind_v);
            depth(ind_v,1) = z2(ind_v);
	end
	    
        %%%%% left dipping (looking from the endpoint) %%%%%
        ind_l = dip>90;
        dip(ind_l) = dip(ind_l) - 180;      %  dip [degree] conversion
        % use the top center
	% (cx,cy) is the anchor point at the burial depth corresponding to the endpoint
        if any(ind_l);
            cx(ind_l,1) = xx(ind_l)+z1(ind_l)./tand(dip(ind_l)).*cosd(str(ind_l));
            cy(ind_l,1) = yy(ind_l)-z1(ind_l)./tand(dip(ind_l)).*sind(str(ind_l));
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
