function [ sm ] = GTdef_sm1d_3pfwd_rtdw(dd,ds,Nd,Ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     GTdef_sm1d_3pfwd_rtdw            			  %
% This file is actually the same as GTdef_sm1d_3pfwd_rtdw_free.m	  %
% Calculate the first-order derivatives of the slips of the subfaults 	  %
% over the fault surface using 3-point forward finite-difference	  %
% approximation 			  				  %
% Note:									  %
%   Use zero slips for virtual boundary patches including the top boundary%
%   Do not assume the free-surface assumption				  %
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
%    ______________________						  %
%    |	    |      |      |                                               %
%    |  1   | 4(-) |  7   |                                               %
%    |______|______|______|                                               %
%    |	    |      |      |                                               %
%    | 2(-) |  5   | 8(+) |                                               %
%    |______|______|______|                                               %
%    |	    |      |      |                                               %
%    |  3   | 6(+) |  9   |                                               %
%    |______|______|______|                                               %
%								          %
% Finite-difference 3-point forward approximation of 1st derivative	  %
% (ith row and jth colunm)						  %
%   4*Si+1,j - 3*Si,j - Si+2,j	     4*Si,j+1 - 3*Si,j - Si,j+2	  	  %
% ------------------------------ +  -----------------------------	  %
%             2*dd			       2*ds			  %
% The coordinate has cartesian convection 				  %
% (right is x+; down is y+)					  	  %
%								 	  %
% first created by Lujia Feng Wed Dec  9 19:44:08 EST 2009		  %
% last modified by Lujia Feng Wed Dec  9 19:45:44 EST 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = Ns*Nd;			% total patch/slip number
mx = 0.5/ds; 			% the patches on left and right
my = 0.5/dd;			% the patches on top and bottom
mx0 = -3*mx;  my0 = -3*my;	% coeff for the point where the derivative is calculated
mx1 = 4*mx;   my1 = 4*my;       % coeff for 1st point forward from point 0
mx2 = -mx;    my2 = -my;	% coeff for 2nd point forward from point 0

if nn==1, error('Only one patch. No need to smooth!');  end

% only one horizontal layer
if Nd==1			
    d0 = (mx0+my0)*ones(nn,1);	% diagonals = 0
    d1 = mx1*ones(nn,1); 	% diagonal + 1
    d2 = mx2*ones(nn,1);	% diagonal + 2
    B = [ d0 d1 d2 ];		% diagonal columns
    ind = [ 0 1 2 ];		% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% only one vertical layer
if Ns==1			
    d0 = (mx0+my0)*ones(nn,1);	% diagonals = 0
    d1 = my1*ones(nn,1); 	% diagonal - 1
    d2 = my2*ones(nn,1);		% diagonal - 2
    B = [ d0 d1 d2 ];			% diagonal columns
    ind = [ 0 1 2 ];			% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% two horizontal layers
if Nd==2			% includes 2x2
    d0 = (mx0+my0)*ones(Nd,Ns);		% diagonals = 0
    d0 = reshape(d0,[],1);
    dn  = mx1*ones(nn,1);		% diagonal + Nd
    dn2 = mx2*ones(nn,1);		% diagonal + 2*Nd
    d01 = my1*ones(Nd,Ns);
    d01(1,:) = 0;
    d01 = reshape(d01,[],1);
    B = [ d0 d01 dn dn2 ]; 	% diagonal columns
    ind = [ 0 1 Nd 2*Nd ]; 	% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% two vertical layers
if Ns==2		
    d0 = (mx0+my0)*ones(Nd,Ns);		% diagonals = 0
    d0 = reshape(d0,[],1);
    dn  = mx1*ones(nn,1);		% diagonal + Nd
    d01 = my1*ones(Nd,Ns);
    d01(1,:) = 0;
    d01 = reshape(d01,[],1);
    d02 = my2*ones(Nd,Ns);
    d02(1,:) = 0; d02(2,:) = 0;
    d02 = reshape(d02,[],1);
    B = [ d0 d01 d02 dn ]; 	% diagonal columns
    ind = [ 0 1 2 Nd ]; 	% index for diagonal columns
    sm = spdiags(B,ind,nn,nn);
    return;
end

% more layers
d0 = (mx0+my0)*ones(Nd,Ns);		% diagonals = 0
d0 = reshape(d0,[],1);
dn  = mx1*ones(nn,1);		% diagonal + Nd
dn2 = mx2*ones(nn,1);		% diagonal + 2*Nd

d01 = my1*ones(Nd,Ns);
d01(1,:) = 0;
d01 = reshape(d01,[],1);
d02 = my2*ones(Nd,Ns);
d02(1,:) = 0; d02(2,:) = 0;
d02 = reshape(d02,[],1);

B = [ d0 d01 d02 dn dn2 ]; 	% diagonal columns
ind = [ 0 1 2 Nd 2*Nd ]; 	% index for diagonal columns
sm = spdiags(B,ind,nn,nn);
