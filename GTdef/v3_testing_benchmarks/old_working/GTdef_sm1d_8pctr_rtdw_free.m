function [ sm ] = GTdef_sm1d_8pctr_rtdw_free(dd,ds,Nd,Ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       GTdef_sm1d_8pctr_rtdw_free			  %
% Calculate the first-order derivatives of the slips of the subfaults 	  %
% over the fault surface using 8-point central finite-difference	  %
% approximation 			  				  %
% Note:									  %
%   Use zero slips for virtual boundary patches; 			  %
%   Use the same slips for virtual free-surface patches 		  %
%     as the most top patches.   					  %
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
%    ______________________						  %
%    |	    |      |      |                                               %
%    | 1(-) | 4(-) | 7(+) |                                               %
%    |______|______|______|                                               %
%    |	    |      |      |                                               %
%    | 2(-) |  5   | 8(+) |                                               %
%    |______|______|______|                                               %
%    |	    |      |      |                                               %
%    | 3(-) | 6(+) | 9(+) |                                               %
%    |______|______|______|                                               %
%								          %
% Finite-difference 8-point central approximation of 1st derivative	  %
%  (8) - (2)     (6) - (4)     (7) - (3)      (9) - (1)                   %         
% ----------- + ----------- + ------------ + -----------                  %
%     2*dx	    2*dy          2*dxy	         2*dxy                    %
% The coordinate has cartesian convection 				  %
% (right is x+; up is y+)					  	  %
%								 	  %
% first created by Lujia Feng Thu Dec  3 17:41:09 EST 2009		  %
% last modified by Fri Dec  4 00:59:35 EST 2009				  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = Ns*Nd;			% total patch/slip number
mx = 0.5/ds; 			% the patches on left and right
my = 0.5/dd;			% the patches on top and bottom
mxy = 0.5/sqrt(ds^2+dd^2);

if nn==1, error('Only one patch. No need to smooth!'); end

% more layers
%%%%%%%%% diagonals = 0 %%%%%%%%%
d0 = zeros(Nd,Ns);		
d0(1,:) = -my;
d0 = reshape(d0,[],1);
%%%%%%%%% sub-diagonal = - 1 %%%%%%%%%
t1 = -my*ones(Nd-1,Ns); 	
t2 = zeros(1,Ns);
t3 = [ t1;t2 ];
d10 = reshape(t3,[],1);
%%%%%%%%% super-diagonal = + 1 %%%%%%%%%
d01 = -d10(end:-1:1);

%%%%%%%%% diagonal + Nd %%%%%%%%%
d0n = mx*ones(Nd,Ns);		
d0n(1,:) = mx+mxy;
d0n = reshape(d0n,[],1);
%%%%%%%%% Nd + 1 %%%%%%%%%
d0n1 = mxy*ones(Nd,Ns);		
d0n1(1,:) = 0;
d0n1 = reshape(d0n1,[],1);
%%%%%%%%% Nd - 1 %%%%%%%%%
d0n_1 = mxy*ones(Nd,Ns);		
d0n_1(Nd,:) = 0;
d0n_1 = reshape(d0n_1,[],1);

%%%%%%%%% diagonal - Nd %%%%%%%%%
dn0 = -mx*ones(Nd,Ns);		
dn0(1,:) = mxy-mx;
dn0 = reshape(dn0,[],1);
%%%%%%%%% -Nd + 1 %%%%%%%%%
dn10 = -mxy*ones(Nd,Ns);		
dn10(1,:) = 0;
dn10 = reshape(dn10,[],1);
%%%%%%%%% -Nd - 1 %%%%%%%%%
dn_10 = -mxy*ones(Nd,Ns);		
dn_10(Nd,:) = 0;
dn_10 = reshape(dn_10,[],1);

ind = [ -Nd-1 -Nd -Nd+1 -1 0 1 Nd-1 Nd Nd+1 ]; 	% index for diagonal columns
B = [ dn_10 dn0 dn10 d10 d0 d01 d0n_1 d0n d0n1 ];
sm = spdiags(B,ind,nn,nn);
