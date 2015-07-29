function [ str ] = GTdef_strike(x1,y1,x2,y2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_strike				  %
% Calculate the trike from [x1,y1] to [x2,y2]				  %
% INPUT:					  		  	  %
% [x1,y1] and [x2,y2] are both in cartesian coordiante		          %
% they can be scalars and column vectors				  %
%									  %
% OUTPUT:								  %
% strke is CW from N							  %
%									  %
% first created by Lujia Feng Mon Dec  7 06:45:49 EST 2009		  %
% last modified by Lujia Feng Tue Dec  8 18:12:53 EST 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ row,col ] = size(x1);
if col~=1&&row==1, error('Inputs have to be vectors for GTdef_strike'); end

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
