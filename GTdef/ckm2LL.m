function [lon,lat] = ckm2LL(xx,yy,lon0,lat0,rot)

%  This matlab program 
%  (1) rotates the Cartesion coordinate by rot [degree]
%      to the coordinate with horizontal x-axis 
%  (2) transforms local Cartesian coordinate to lon,lat
%      using (lon0,lat0) as origin
%  Input:
%  (1) xx,yy [m] can be scalars, vectors or matrices,
%      but they must have the same size
%  (2) (lon0,lat0) [degree] is used as the origin reference point
%  (3) rot [radian] is the rotation angle from the old coordinate
%      system (east is x+ ;north is y+) to the new one.
%      Counterclockwise rotation is positive
%      Clockwise rotation is negative
%  Output:
%  (1) lon,lat [degree] may be scalars, vectors or matrices
%      depending on the input xx,yy [m]
%  Note: usually used with LL2ckm()  
%  (1) use the same (lon0,lat0) for both
%  (2) use rot for LL2ckm(); use -rot for ckm2LL()
%  Other related functions: ckm2LLd(), LL2ckmd() and rotateckmd()
%  first written by Amenda Thomas 2006
%  modified and commented by Lujia Feng Thu Aug 28 17:15:21 EDT 2008
%  last modified by lfeng Wed Nov 17 16:14:06 EST 2010

R = 6378137;                % R - Earth's radius at equator [m]
% assume a perfect spheric Earth
%r = R;
% assume an oblate ellipsoid Earth [WGS 84 (World Geodetic System)]
ff = 1/298.257;             % ff - flattening factor
r = R*(1-ff*sin(lat0*pi/180)^2);  % r - radius at lat [m]

mpd = r*pi/180;             % mpd - meters per degree
cos_rot = cos(rot);         % cos_rot - cos of rotation angle
sin_rot = sin(rot);         % sin_rot - sin of rotation angle

if size(xx)==size(yy)
   % rotate the coordinate system by rot [degree]
   % back to the state of horizontal x-axis
   x_rot = xx*cos_rot + yy*sin_rot;
   y_rot =-xx*sin_rot + yy*cos_rot;
   % transform from xx,yy to lon,lat using (lon0,lat0) as origin
   lat = lat0 + y_rot/mpd;
   lon = lon0 + x_rot/mpd./cos(lat0*pi/180);
else
   error('xx and yy are not consistent!');
end
