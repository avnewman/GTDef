function [ sm ] = GTdef_smooth2d(dd,ds,Nd,Ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_smooth2d				  %
% Calculate the second-order derivatives of the slips of the subfaults 	  %
% over the fault surface. 						  %
% Note:									  %
%   Use zero slips for virtual boundary patches; 			  %
%   Use the same slips for virtual free-surface patches 		  %
%     as the most top patches.   					  %
%   Only consider one type of slip here.			  	  %
%									  %
% INPUT:								  %
%   dd - distance along dip between two vertically adjacent patches	  %
%   ds - distance along strike between two horizontally adjacent patches  %
%   Nd - number of patches along dip				          %
%   Ns - number of patches along strike					  %
% 									  %
% OUTPUT:								  %
%   sm - smoothing matrix						  %
%								          %
% REFERENCE:								  %
% Jonsson, et al. (2002), BSSA, 92(4), 1377-1389 eq (A1)		  %
%  Si-1,j - 2Si,j + Si+1,j	 Si,j-1 - 2Si,j + Si,j+1		  %
% -------------------------  +  -------------------------		  %
%          dx*dx			 dy*dy				  %
%								 	  %
% related function: GTdef_smooth1d()					  %
% first created by Lujia Feng Fri Apr 24 18:40:18 EDT 2009		  %
% last modified by Lujia Feng Fri May  8 17:35:31 EDT 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = Ns*Nd;			% total patch/slip number
mx = 1/ds^2; 			% the patches on left and right
my = 1/dd^2;			% the patches on top and bottom
mxy = -2*mx-2*my;		% the patch in the center
mxy_free = -2*mx-my;		% the free-surface patch in the center

if nn==1, disp('Only one patch. No need to smooth!'); return; end

% only one horizontal layer
if Nd==1			
    d0 = mxy_free*ones(nn,1);	% diagonals = 0
    d1 = mx*ones(nn,1); 	% diagonal +/- 1
    B = [ d1 d0 d1 ];		% diagonal columns
    ind = [ -1 0 1 ];		% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% only one vertical layer
if Ns==1			
    d0 = [mxy_free;mxy*ones(nn-1,1)];	% diagonals = 0
    d1 = my*ones(nn,1); 	% diagonal +/- 1
    B = [ d1 d0 d1 ];		% diagonal columns
    ind = [ -1 0 1 ];		% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% more layers
d0 = mxy*ones(Nd,Ns);		% diagonals = 0
d0(1,:) = mxy_free;
d0 = reshape(d0,[],1);
dn = mx*ones(nn,1);		% diagonal +/- Nd
% diagonal +/- 1
t1 = my*ones(Nd-1,Ns); 	
t2 = zeros(1,Ns);
t3 = [ t1;t2 ];
d10 = reshape(t3,[],1);
d01 = d10(end:-1:1);
B = [ dn d10 d0 d01 dn ]; 	% diagonal columns
ind = [ -Nd -1 0 1 Nd ]; 	% index for diagonal columns
sm = spdiags(B,ind,nn,nn);
