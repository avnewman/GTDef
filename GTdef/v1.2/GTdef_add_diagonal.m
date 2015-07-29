function [ AB ] = GTdef_add_diagonal(A,B)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_add_diagonal				  %
% Add two matrices diagonally						  %
% AB =	[ A  0 ]							  %
%	[ 0  B ]							  %
%									  %
% first created by Lujia Feng Fri Apr 24 22:52:07 EDT 2009		  %
% last modified by Lujia Feng Fri Apr 24 22:55:47 EDT 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m0 n0] = size(A);    [m1 n1] = size(B);
upright = zeros(m0,n1); botleft = zeros(m1,n0);
AB = [ A upright; botleft B ];
