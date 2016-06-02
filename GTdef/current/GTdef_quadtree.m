function [Dnew,Dbox]=GTdef_quadtree(D,X,Y,ndown,dDmax,mpd,method,dimwgt,LOSd,outfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_quadtree				                 %
%  Given a structured matrix, D, with values described by the position vectors, X, Y     %
%  data can be downsampled by ndown, and then split into equidimensional boxes when      %
%  data within an initial box changes by more than dDmax.  Boxes will not be smaller     %
%  than mpd.  Optionally, you can write out files with:                                  %
%     .in)  [format for GTdef_input] consisting of:                                      %
%         los type name lon lat z  ULOS  eULOS LOSdirE  LOSdirN  LOSdirV  weight         %
%     .txt) the value, position defined by the center of mass of the data, the std.      %
%        deviation of the data dimension of the box, and number of data within box.      %
%     .box) a 'gmt-style' file that contains the description of the boxes used for       %
%        plotting.									 %
%     *prefix defined by 'outfile'							 %
% INPUT:					  		  	                 %
%   D  - matrix of data                                                                  %
%   X,Y - Vectors defining positions within strucutre matrix D                           %
%   ndown - integer defining subsampling {0}                                             %
%   dDmax - float defining  maximum allowable change parameter {1}                       %
%   mpd - integer defining minimum pixel dimension (2^mpd) {2}                           %
%   method - integer defining method for defining dDmax {1}                              %
%       1 = difference                                                                   %
%       2 = percent change (10=10% across entire range of D)                             %
%   dimwgt - integer defining whether weight                                             %
%	0 = set all to 1                                                                 %
%       1 = use dimension of one side, normalized by min dimension  (2x2 = dim 2)        %
%       2 = use power of dimension, normed by min power (2^pwr; 2x2 = pwr 1; 8x8 = pwr 3;%
%   LOSd - (optional) unit vector giving look direction to satelite ([E N U]) {[0 0 1.]} %
%   outfile - (optional) string definining outfile prefixes                              %
% OUTPUT: (same info as in outfiles above)	  		  	                 %
%   Dnew structure                                                                       %
%   Dbox structure                                                                       %
% first created by Andrew Newman Tue Apr 19 17:23:30 EDT 2016 			         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   weighting method
     if (exist('method','var')==1)
	     if(method>2)
	             method =1;
		     disp('GTdef_quadtree WARNING: method only has option 1 or 2. setting to 1');
	     elseif(method==2)
		     Drange=max(max(D)')-min(min(D)');
		     dDmax=dDmax/100.*Drange;   % rescale by range
		     disp(['GTdef_quadtree MESSAGE: new dDmax is ',num2str(dDmax),'. Scaled by range.']);
	     end
     else
	     method=1;
     end
%   outputfile
     % create outfiles
     if (exist('outfile','var')==1)
       FIN=fopen(strcat(outfile,'.in'),'w');
       FTXT=fopen(strcat(outfile,'.txt'),'w');
       FBOX=fopen(strcat(outfile,'.box'),'w');
     end

% initialize parameters
Dnew.D    = []; Dnew.err  = [];
Dnew.x    = []; Dnew.y    = [];
Dnew.Dpos = []; Dnew.dim  = [];
Dnew.Ndat = []; 
Dbox.x1   = []; Dbox.y1   = [];
Dbox.x2   = []; Dbox.y2   = [];

%downsample test
  Ddown=downsample(downsample(D,ndown)',ndown)';
  Xdown=downsample(X,ndown);Ydown=downsample(Y,ndown);
 
% determine size of array
[Nx,Ny]=size(Ddown);
% pad array to next power of 2^n
  p=nextpow2(max([Nx Ny]));
  Dpad=padarray(Ddown,[2^p-Nx 2^p-Ny],NaN,'post');

% quadtree decomposition for boxes with more than umax displacement, and within mmin and mmax array size
mmin=2^mpd;
mmax=2^(p-1);
SD=qtdecomp(Dpad,dDmax,[mmin mmax]);

% recipe for plotting from http://www.mathworks.com/help/images/ref/qtdecomp.html
dim=0; n=0;
j=1;
while dim <= mmax
   n=n+1 ;
   dim=2^n ;
 % EAST 
   numblocks=length(find(SD==dim));
   if (numblocks > 0)
      values = repmat (uint8(1), [dim dim numblocks]);
      values(2:dim,2:dim,:)=0;
      [vals,r,c] = qtgetblk(Dpad,SD,dim);  
       for i=1:length(r)
       % extract
          % n NaN
          nnan=sum(sum(isnan(vals(:,:,i))));
          if nnan<dim^2  % only do other stuff if data exists
              
              % number of data in quadtree box
              ndata=dim^2-nnan;

              % center of mass position based on existance of data
              com=centerOfMass(int16(not(isnan(vals(:,:,i))))); % row and column of center of mass of resolved pixels.
              rcom=round(com(1));
              ccom=round(com(2));
              % value at the center of mass position
              vcom=vals(rcom,ccom,i);
              
              % mean
              val=nanmean(reshape(vals(:,:,i),1,numel(vals(:,:,i))));
              % stddev
              stderr=nanstd(reshape(vals(:,:,i),1,numel(vals(:,:,i))));
              
              xpos=Xdown(c(i)+ccom); % b/c of padding, the mid-point of a box can be outside of the original dimensions
              ypos=Ydown(r(i)+rcom); 
              
              x1=Xdown(c(i)); y1=Ydown(r(i));
              x2=Xdown(min(length(Xdown),c(i)+dim)); y2=Ydown(min(length(Ydown),r(i)+dim));
              
	      % disp(j)
              Dnew.D(j)    = val;   Dnew.err(j)  = stderr;
              Dnew.x(j)    = xpos;  Dnew.y(j)    = ypos;
              Dnew.Dpos(j) = vcom;  Dnew.dim(j)  = dim;
              Dnew.Ndat(j) = ndata; 
              Dbox.x1(j)   = x1;    Dbox.y1(j)   = y1;
              Dbox.x2(j)   = x2;    Dbox.y2(j)   = y2;

	      % weight output based on number of data, or relative box size?
	       wgt=1.;
	       if (dimwgt==1)
	          wgt=dim/mmin;
	       elseif (dimwgt==2)
	          wgt=n/mpd;
	       end

              if (exist('outfile','var')==1)
                fprintf(FIN,'los  1 %s  %13.8f %13.8f  0.0000   %9.5f %9.5f   %6.4f %6.4f %6.4f  %6.2f\n',outfile, xpos,ypos, val,stderr, LOSd(1), LOSd(2), LOSd(3), wgt);
                fprintf(FTXT,'%9.5f %9.5f %9.5f %13.8f %13.8f    %6d %7d\n',val,stderr,vcom, xpos,ypos,dim,ndata);
                fprintf(FBOX,'%10.5f %10.5f\n%10.5f %10.5f\n%10.5f %10.5f\n%10.5f %10.5f\n%10.5f %10.5f\n>\n',x1,y1,x2,y1,x2,y2,x1,y2,x1,y1);
	      end
	      j=j+1;
          end
       end
   end
end

 if (exist('outfile','var')==1)
    fclose(FIN);
    fclose(FTXT);
    fclose(FBOX);
 end

end % MAIN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = centerOfMass(A,varargin) 
% CENTEROFMASS finds the center of mass of the N-dimensional input array 
% 
%   CENTEROFMASS(A) finds the gray-level-weighted center of mass of the 
%   N-dimensional numerical array A. A must be real and finite. A warning 
%   is issued if A contains any negative values. Any NaN elements of A will 
%   automatically be ignored. CENTEROFMASS produces center of mass 
%   coordinates in units of pixels. An empty array is returned if the 
%   center of mass is undefined. 
% 
%   The center of mass is reported under the assumption that the first 
%   pixel in each array dimension is centered at 1. 
% 
%   Also note that numerical arrays other than DOUBLE and SINGLE are 
%   converted to SINGLE in order to prevent numerical roundoff error. 
% 
%   Examples: 
%       A = rgb2gray(imread('saturn.png')); 
%       C = centerOfMass(A); 
% 
%       figure; imagesc(A); colormap gray; axis image 
%       hold on; plot(C(2),C(1),'rx') 
% 
%   See also:  
% 
% 
% 
%   Jered R Wells 
%   2013/05/07 
%   jered [dot] wells [at] gmail [dot] com 
% 
%   v1.0 
% 
%   UPDATES 
%       YYYY/MM/DD - jrw - v1.1 
% 
%   Copyright (c) 2013, Jered  Wells 
%   All rights reserved. 
%    
%   Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are 
%   met: 
%    
%       * Redistributions of source code must retain the above copyright 
%         notice, this list of conditions and the following disclaimer. 
%       * Redistributions in binary form must reproduce the above copyright 
%         notice, this list of conditions and the following disclaimer in 
%         the documentation and/or other materials provided with the distribution 
%    
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%   POSSIBILITY OF SUCH DAMAGE. 
%    
% 
 
%% INPUT CHECK 
narginchk(0,1); 
nargoutchk(0,1); 
fname = 'centerOfMass'; 
 
% Checked required inputs 
validateattributes(A,{'numeric'},{'real','finite'},fname,'A',1); 
 
%% INITIALIZE VARIABLES 
A(isnan(A)) = 0; 
if ~(strcmpi(class(A),'double') || strcmpi(class(A),'single')) 
    A = single(A); 
end 
if any(A(:)<0) 
    warning('MATLAB:centerOfMass:neg','Array A contains negative values.'); 
end 
 
%% PROCESS 
sz = size(A); 
nd = ndims(A); 
M = sum(A(:)); 
C = zeros(1,nd); 
if M==0 
    C = []; 
else 
    for ii = 1:nd 
        shp = ones(1,nd); 
        shp(ii) = sz(ii); 
        rep = sz; 
        rep(ii) = 1; 
        ind = repmat(reshape(1:sz(ii),shp),rep); 
        C(ii) = sum(ind(:).*A(:))./M; 
    end 
end 
 
% Assemble the VARARGOUT cell array 
varargout = {C}; 
 
end % MAIN 
