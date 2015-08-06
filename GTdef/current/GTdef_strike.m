function [ str ] = GTdef_strike(x1,y1,x2,y2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_strike				  %
% Calculate the trike from [x1,y1] to [x2,y2]				  %
% INPUT:					  		  	  %
% [x1,y1] and [x2,y2] are both in cartesian coordiante		          %
% they can be scalars and column vectors				  %
%									  %
% OUTPUT:								  %
% strke is CW from N [0 360]                                              %
%									  %
% first created by Lujia Feng Mon Dec  7 06:45:49 EST 2009		  %
% added range control lfeng Wed Aug  5 12:51:39 SGT 2015                  %
% last modified by Lujia Feng Wed Aug  5 12:53:08 SGT 2015                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%str   = 90-180*atan2(y2-y1,x2-x1)/pi;   % degree CW from N [0-360]

[ row,col ] = size(x1);
if col~=1&&row==1, error('GTdef_strike ERROR: inputs have to be vectors!'); end

for ii =1:row
   if y1(ii,1)==y2(ii,1)&&x1(ii,1)>x2(ii,1)
      str(ii,1) = 270;
   elseif y1(ii,1)==y2(ii,1)&&x1(ii,1)<x2(ii,1)
      str(ii,1) = 90;
   elseif y2(ii,1)-y1(ii,1)>0
      str(ii,1) = atand((x2(ii,1)-x1(ii,1))/(y2(ii,1)-y1(ii,1))); 
   else
      str(ii,1) = atand((x2(ii,1)-x1(ii,1))/(y2(ii,1)-y1(ii,1)))+180; 
   end
end

% make str value within [0-360]
ind = str<0;
str(ind) = str(ind) + 360;
ind = str>360;
str(ind) = str(ind) - 360;
