%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_convertlon.m                                %
%             function to convert longitude [0 360] to [-180 180]               %
% INPUT&OUPUT lon is a row or column vector                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ lon ] = GTdef_convertlon(lon)

ind = lon>180;
lon(ind) = lon(ind)-360;
end