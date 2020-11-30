function [ x_grd,y_grd,nod_name ] = GTdef_grid(grd,grd_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GTdef_grid				  %
% Determine all the grid nodes from a grid definition			  %
%									  %
% INPUT:					  		  	  %
%  grd = [ xrot,yrot,x1,y1,x2,y2,xnum,ynum ] (1*8)			  %
%    xrot  - rotation angle for x axis (east) of grid [degree] (-90 90)	  %
%    yrot  - rotation angle for y axis (north) of grid [degree] (-90 90)  %
%             + CW; - CCW; xrot and yrot can not be 90 or -90		  %
%    x1,y1 - the lower-left corner of the grid         		  	  %
%    x2,y2 - the upper-right corner of the grid         		  %
%    x3,y3 - the lower-right corner of the grid         		  %
%                      (x2-x1)-tand(yrot)*(y2-y1)			  %
%          x3 = x1 + --------------------------------			  %
%                        1+tand(xrot)*tand(yrot)			  %
%									  %
%                      (y2-y1)+tand(xrot)*(x2-x1)			  %
%          y3 = y2 - --------------------------------			  %
%                        1+tand(xrot)*tand(yrot)			  %
%									  %
%    x4,y4 - the upper-left corner of the grid         		  	  %
%                      (x2-x1)-tand(yrot)*(y2-y1)			  %
%          x4 = x2 - --------------------------------			  %
%                        1+tand(xrot)*tand(yrot)			  %
%									  %
%                      (y2-y1)+tand(xrot)*(x2-x1)			  %
%          y4 = y1 + --------------------------------			  %
%                        1+tand(xrot)*tand(yrot)			  %
%									  %
%    xnum  - number of grid nodes along x axis				  %
%    ynum  - number of grid nodes along y axis				  %
%                                                                         %
% OUTPUT:                                                                 %
%  x_grd - row vector for x values of the grid nodes			  %
%  y_grd - row vector for y values of the grid nodes			  %
%                                                                         %
% first created by Lujia Feng Fri May  1 14:15:53 EDT 2009		  %
% use cell array of strings for names lfeng Wed Dec  1 17:20:21 EST 2010  %
% last modified by Lujia Feng Wed Dec  1 17:20:29 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xrot = grd(1); yrot = grd(2); 
x1 = grd(3);   y1 = grd(4);  x2 = grd(5); y2 = grd(6);
xnum = grd(7); ynum = grd(8);

% calculate the lower-right and upper-left corners of the grid
bot = 1+tand(xrot)*tand(yrot);
dx = (x2-x1-tand(yrot)*(y2-y1))/bot;
dy = (y2-y1+tand(xrot)*(x2-x1))/bot;

x3 = x1+dx; x4 = x2-dx;
y3 = y2-dy; y4 = y1+dy;

% nodes on the bottom side of the grid
xbot = linspace(x1,x3,xnum); ybot = linspace(y1,y3,xnum);
% nodes on the top side of the grid
xtop = linspace(x4,x2,xnum); ytop = linspace(y4,y2,xnum);

for ii = 1:xnum
   x_col = linspace(xtop(ii),xbot(ii),ynum);
   y_col = linspace(ytop(ii),ybot(ii),ynum);
   xmat(:,ii) = x_col;
   ymat(:,ii) = y_col;
end
x_grd = reshape(xmat,1,[]);
y_grd = reshape(ymat,1,[]);

nod_name = {};
for ii=1:xnum
   for jj=1:ynum
       name = strcat(grd_name,'_',int2str(jj),'_',int2str(ii));
       nod_name = [ nod_name; name ];   
   end
end
