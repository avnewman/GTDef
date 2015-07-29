function [ sm ] = GTdef_sm1d_2pctr_free(dd,ds,Nd,Ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        GTdef_sm1d_2pctr_free			  	  %
% Calculate the first-order derivatives of the slips of the subfaults 	  %
% over the fault surface using 2-point central finite-difference	  %
% approximation 			  				  %
% Note:									  %
%   Use zero slips for virtual boundary patches; 			  %
%   Use the same slips for virtual free-surface patches 		  %
%   as the most top patches.   					  	  %
%   Free-surface assumption is good when coseismic rupture reaches surface%
%   Only consider one type of slip here.			  	  %
%									  %
% INPUT:								  %
%   dd - distance along dip between two vertically adjacent patches	  %
%   ds - distance along strike between two horizontally adjacent patches  %
%   Nd - number of patches along dip				          %
%   Ns - number of patches along strike					  %
% 									  %
% OUTPUT:								  %
%   sm - smoothing matrix [nn*nn]  (nn = Ns*Nd)				  %
% Note: order slips columnwise similar to matlab reshape function	  %
%    _____________							  %
%  + |	 |   |   |                                                        %
%    | 1 | 3 | 5 |                                                        %
%    |___|___|___|                                                        %
%    |	 |   |   |                                                        %
%    | 2 | 4 | 6 |                                                        %
%  - |___|___|___| +                                                      %
%								          %
% REFERENCE:								  %
% Finite-difference 2-point central approximation of 1st derivative	  %
% (ith row and jth colunm)						  %
%      Si+1,j - Si-1,j	           Si,j+1 - Si,j-1		  	  %
% -------------------------  +  -------------------------		  %
%          2*dd				 2*ds				  %
% The coordinate has cartesian convection 				  %
% (right is x+; up is y+)					  	  %
%								 	  %
% first created by Lujia Feng Tue Dec  1 14:46:13 EST 2009		  %
% last modified by Lujia Feng Tue Dec  1 23:04:11 EST 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = Ns*Nd;			% total patch/slip number
mx = 0.5/ds; 			% the patches on left and right
my = 0.5/dd;			% the patches on top and bottom

if nn==1, error('Only one patch. No need to smooth!'); end

% only one horizontal layer
if Nd==1			
    d0 = my*ones(nn,1);		% diagonals = 0
    d1 = mx*ones(nn,1); 	% diagonal +/- 1
    B = [ -d1 d0 d1 ];		% diagonal columns
    ind = [ -1 0 1 ];		% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% only one vertical layer
if Ns==1			
    d0 = [my;0*ones(nn-1,1)];	% diagonals = 0
    d1 = my*ones(nn,1); 	% diagonal +/- 1
    B = [ d1 d0 -d1 ];		% diagonal columns
    ind = [ -1 0 1 ];		% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% more layers
d0 = zeros(Nd,Ns);		% diagonals = 0
d0(1,:) = my;
d0 = reshape(d0,[],1);
dn = mx*ones(nn,1);		% diagonal +/- Nd
% sub-diagonal = diagonal - 1
t1 = my*ones(Nd-1,Ns); 	
t2 = zeros(1,Ns);
t3 = [ t1;t2 ];
d10 = reshape(t3,[],1);
% super-diagonal = diagonal + 1
d01 = -d10(end:-1:1);
B = [ -dn d10 d0 d01 dn ]; 	% diagonal columns
ind = [ -Nd -1 0 1 Nd ]; 	% index for diagonal columns
sm = spdiags(B,ind,nn,nn);
