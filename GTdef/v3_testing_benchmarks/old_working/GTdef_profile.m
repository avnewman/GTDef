function [ x_prf,y_prf,nod_name ] = GTdef_profile(prf,prf_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_profile				  %
% Determine all the profile nodes from a profile definition		  %
%									  %
% INPUT:					  		  	  %
%  prf = [x1 y1 x2 y2 num] (1*5)					  %
%    x1,y1 - the lower-left corner of the profile         		  %
%    x2,y2 - the upper-right corner of the profile         		  %
%    num   - number of profile nodes				 	  % 
%  prf_name - the name of the profile					  %
%                                                                         %
% OUTPUT:                                                                 %
%  x_prf - row vector for x values of the profile nodes			  %
%  y_prf - row vector for y values of the profile nodes			  %
%                                                                         %
% first created by Lujia Feng Fri May  1 15:51:45 EDT 2009		  %
% use cell array of strings for names lfeng Wed Dec  1 17:21:46 EST 2010  %
% last modified by Lujia Feng Tue May  5 15:06:51 EDT 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = prf(1); y1 = prf(2); x2 = prf(3); y2 = prf(4); num = prf(5);

x_prf = linspace(x1,x2,num); 
y_prf = linspace(y1,y2,num);
nod_name = {};
for ii=1:num
    name = strcat(prf_name,'_',int2str(ii));
    nod_name = [ nod_name; name ];   
end
