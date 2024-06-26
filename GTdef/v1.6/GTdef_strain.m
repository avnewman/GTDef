function [ sm ] = GTdef_strain(dd,ds,Nd,Ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_strain				  %
% Calculate the first-order derivatives of the slips of the subfaults 	  %
% over the fault surface using 2-point center approximation		  %
% Separate x and y directions						  %
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
%    _____________							  %
%  + |	 |   |   |                                                        %
%    | 1 | 3 | 5 |                                                        %
%    |___|___|___|                                                        %
%    |	 |   |   |                                                        %
%    | 2 | 4 | 6 |                                                        %
%  - |___|___|___| +                                                      %
%									  %
%       Si,j - Si-1,j							  %
% (1)  ------------------						  %
%            dx                                                           %
%                                                                         %
%       Si,j - Si,j-1                                                     %
% (2)  ------------------                                                 %
%            dx                                                           %
% Take absolute value later in GTdef_forward.m				  %
%								 	  %
% first created by Lujia Feng Wed Dec  9 17:20:47 EST 2009		  %
% last modified by Lujia Feng Wed Dec  9 18:11:22 EST 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = Ns*Nd;			% total patch/slip number
mx = 0.5/ds; 			% the patches on left and right
my = 0.5/dd;			% the patches on top and bottom

if nn==1, error('Only one patch. No need to smooth!'); end

% x direction (from left to right)
dn = mx*ones(nn,1);		% diagonal +/- Nd
B = [ -dn dn ]; 		% diagonal columns
ind = [ -Nd Nd ]; 		% index for diagonal columns
smx = spdiags(B,ind,nn,nn);

% y direction (from down to up)
t1 = my*ones(Nd-1,Ns); 		% sub-diagonal = diagonal - 1
t2 = zeros(1,Ns);
t3 = [ t1;t2 ];
d10 = reshape(t3,[],1);
d01 = -d10(end:-1:1);		% super-diagonal = diagonal + 1
B = [ d10 d01 ]; 	% diagonal columns
ind = [ -1 1 ]; 	% index for diagonal columns
smy = spdiags(B,ind,nn,nn);

sm = [ smx;smy ];  
