function [x_rot,y_rot] = LL2ckm(lon,lat,lon0,lat0,rot)

%  This matlab program 
%  (1) transforms lon,lat to local Cartesian coordinate
%      using (lon0,lat0) as origin
%  (2) rotates the Cartesion coordinate by rot [degree]
%  Input:
%  (1) lon,lat [degree] can be scalars, vectors or matrices,
%      but they must have the same size
%  (2) (lon0,lat0) [degree] is used as the origin reference point
%  (3) rot [radian] is the rotation angle from the old coord
%      system (east is x+ ;north is y+) to the new one.
%      Counterclockwise rotation is positive
%      Clockwise rotation is negative
%  Output:
%  (1) x_rot,y_rot [m] may be scalars, vectors or matrices
%      depending on the input lon,lat [degree]
%  Note: usually used with ckm2LL()
%  (1) use the same (lon0,lat0) for both
%  (2) use rot for LL2ckm(); use -rot for ckm2LL()
%  Other related functions: ckm2LLd(), LL2ckmd() and rotateckmd()
%
%  first written by Amenda Thomas 2006
%  modified and commented by Lujia Feng Wed Aug 27 14:42:26 EDT 2008
%  last modified by lfeng Wed Nov 17 16:31:28 EST 2010

R = 6378137;                % R - Earth's radius at equator [m]
% assume a perfect spheric Earth
%r = R;
% assume an oblate ellipsoid Earth [WGS 84 (World Geodetic System)]
ff = 1/298.257;             % ff - flattening factor
r = R*(1-ff*sin(lat0*pi/180)^2);  % r - radius at lat [m]

mpd = r*pi/180;             % mpd - meters per degree
cos_rot = cos(rot);         % cos_rot - cos of rotation angle
sin_rot = sin(rot);         % sin_rot - sin of rotation angle

if size(lon)==size(lat)
   % transform from lon,lat to xx,yy using (lon0,lat0) as origin
   yy = (lat-lat0)*mpd;
   xx = (lon-lon0)*mpd.*cos(lat0*pi/180);
   % rotate the coordinate system by rot [degree]
   x_rot = xx*cos_rot + yy*sin_rot;
   y_rot =-xx*sin_rot + yy*cos_rot;
else
   error('lon and lat are not consistent!'); 
end
