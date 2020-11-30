function [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_samedips(mx1,my1,mx2,my2,mz1,mz2,mstr,mdip,Nd,Ns,sweepAngle)
%                                                        1   2   3   4   5   6   7    8    9  10 11

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GTdef_samedips.m				  %
% Calculate x1,y1,x2,y2,z1,z2,dip,ddip of subfaults for planar faults	  %
%								          %
% INPUT:					  		  	  %
%    master fault info							  %
%    mx1,my1,mx2,my2,mz1,mz2,mstr,mdip,Nd,Ns                              %
%    sweepAngle - angle between normal direction of fault trace and east  %
%                E=0 N=90 W=180 S=270 [deg]                               %
%									  %
% OUTPUT: all column vectors					          %
%    x1,y1 - one endpoint among the two endpoints                         %
%    x2,y2 - the other endpoint among the two endpoints 		  %
%            both in the local cartesian coordinate system	          %
%    z1  - vertical burial depth (top of fault)                           %
%    z2  - vertical locking depth (bottom of fault)                       %
%    dip - down from Horiz, right looking from the endpoint 1 [0 180]     %
%   ddip - average patch size along-dip                                   %
%									  %
% first created by Lujia Feng Fri Oct 24 15:50:59 SGT 2014                %
% added sweepAngle lfeng Wed Nov  5 19:44:17 SGT 2014                     %
% last modified by Lujia Feng Wed Nov 12 18:04:12 SGT 2014                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 10
   sweepAngle = 0;
end

xlin  = linspace(mx1,mx2,Ns+1); 
ylin  = linspace(my1,my2,Ns+1); 
zlin  = linspace(mz1,mz2,Nd+1)';

x1mat = xlin(ones(Nd,1),1:end-1); y1mat = ylin(ones(Nd,1),1:end-1);	% endpoint 1
x2mat = xlin(ones(Nd,1),2:end);   y2mat = ylin(ones(Nd,1),2:end);	% endpoint 2
z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));	% depths

% adjust endpoints according to sweepAngle except for 1st row
if sweepAngle~=0
   for ii =2:Nd
      % get the top left corner of each patch
      dwidth    = z1mat(ii-1,:)/tand(mdip);
      x1TopLeft = x1mat(ii-1,:)+dwidth.*cosd(mstr);
      y1TopLeft = y1mat(ii-1,:)-dwidth.*sind(mstr);
      x2TopRigt = x2mat(ii-1,:)+dwidth.*cosd(mstr);
      y2TopRigt = y2mat(ii-1,:)-dwidth.*sind(mstr);
      % apply sweep angle
      dwidth    = (z2mat(ii-1,:)-z1mat(ii-1,:))/tand(mdip);
      x3TopLeft = x1TopLeft+dwidth.*cosd(sweepAngle);
      y3TopLeft = y1TopLeft+dwidth.*sind(sweepAngle);
      x4TopRigt = x2TopRigt+dwidth.*cosd(sweepAngle);
      y4TopRigt = y2TopRigt+dwidth.*sind(sweepAngle);
      % get back the interception point at surface
      dwidth      = z1mat(ii,:)/tand(mdip);
      x1mat(ii,:) = x3TopLeft-dwidth.*cosd(mstr);
      y1mat(ii,:) = y3TopLeft+dwidth.*sind(mstr);
      x2mat(ii,:) = x4TopRigt-dwidth.*cosd(mstr);
      y2mat(ii,:) = y4TopRigt+dwidth.*sind(mstr);
   end
end

x1 = reshape(x1mat,[],1);  y1 = reshape(y1mat,[],1);
x2 = reshape(x2mat,[],1);  y2 = reshape(y2mat,[],1);
z1 = reshape(z1mat,[],1);  z2 = reshape(z2mat,[],1);

dip  = mdip*ones(Nd*Ns,1);      % duplicate dips
dz   = (mz2-mz1)/Nd; 
ddip = dz/sind(mdip);
