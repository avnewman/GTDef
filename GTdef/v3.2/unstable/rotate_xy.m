function [x_rot,y_rot] = rotate_xy(xx,yy,rot)

% rotates the cartesian coordinate by rot [degree]
%  Input:
%  (1) xx,yy can be scalars, vectors or matrices,
%      but they must have the same size
%  (2) rot [degree] is the rotation angle from the old coordinate
%      system (east is x+ ;north is y+) to the new one.
%      Counterclockwise rotation is positive
%      Clockwise rotation is negative
%  Output:
%  (1) x_rot,y_rot may be scalars, vectors or matrices
%      depending on the input xx,yy

cos_rot = cosd(rot);        % cos_rot - cos of rotation angle
sin_rot = sind(rot);        % sin_rot - sin of rotation angle
% rotate the coordinate system by rot [degree]
x_rot = xx.*cos_rot + yy.*sin_rot;
y_rot =-xx.*sin_rot + yy.*cos_rot;
