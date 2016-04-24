%   datafile type
        [D, RD]=geotiffread('plot_quadtree/ChileQuake_SOLUTION_Vertical_Component.tiff');
        clear RD;
%   geometric constraints
       	X=linspace(-71.7485,-69.9865,1762);
       	Y=linspace(-29.8725,-33.6675,3795);  % qtree decomp starts at top left (pixel 1 1)

	grdwrite2(X,Y,flipud(D),'test.grd')
%   downsample
    ndown=4 ;
%   branching method
     method=1 ;  % absolute difference
     method=2 ;  % percent difference
%   branching threshold
     dDmax=30;    % max difference value
     mpd=2;     % smallest size by power (2^mpd box)
%   dimensional weight
     dimwgt=1 ;  %weight by dimension  n/mmin=weight (smallest weight = 1)
     dimwgt=0 ;  %all = 1
     dimwgt=2 ;  %weight by power 2^n where n/mpd=weight  (smallest weight = 1)
%   LOS vector
     LOSd=[0 0 1.];
     
%   outputfile
       OFpre='test';

[Dnew,Dbox]=GTdef_quadtree(D,X,Y,ndown,dDmax,mpd,method,dimwgt,LOSd,OFpre);
