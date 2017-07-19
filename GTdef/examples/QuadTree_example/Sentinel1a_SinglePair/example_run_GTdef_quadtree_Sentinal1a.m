%   datafile type
      %   [D, RD]=geotiffread('plot_quadtree/ChileQuake_SOLUTION_Vertical_Component.tiff');
       	    %X=linspace(-71.7485,-69.9865,1762);
       	    %Y=linspace(-29.8725,-33.6675,3795);  % qtree decomp starts at top left (pixel 1 1)
        [D, RD]=geotiffread('S1_ChileQuake_Unw_CM_20150824_20150917.tiff');
	    X=linspace(RD.LongitudeLimits(1),RD.LongitudeLimits(2),RD.RasterSize(2)); 
	    Y=linspace(RD.LatitudeLimits(2),RD.LatitudeLimits(1),RD.RasterSize(1)); 
        %clear RD;
%   geometric constraints

%   downsample
    ndown=2 ;
%   branching method
     method=1 ;  % absolute difference
     method=2 ;  % percent difference
%   branching threshold
     dDmax=20;    % max difference value
     mpd=4;     % smallest size by power (2^mpd box)
%   dimensional weight
     dimwgt=1 ;  %weight by dimension  n/mmin=weight (smallest weight = 1)
     dimwgt=0 ;  %all = 1
     dimwgt=2 ;  %weight by power 2^n where n/mpd=weight  (smallest weight = 1)

%   LOS vector
     LOSd=[-0.634 0.170 -0.755]; % Desecending orbit at 15° of long and nominally 41° incidence. (transition between IW2 and IW3 for Sentinal)
%   outputfile
       OFpre='20150824_20150917';
	grdwrite2(X,Y,flipud(D),strcat(OFpre,'.grd'))

%[Dnew,Dbox]=GTdef_quadtree(D,X,Y,ndown,dDmax,mpd,method,dimwgt,LOSd,OFpre);
LOStops.UR1p= [-69.09465 -29.26375]; % lon lat  (Digitized)
LOStops.UR3p= [-70.65689 -28.81895]; % lon lat  (Digitized)
LOStops.azim= [15.0]; % for descending orbit (at least in Chile)
%LOStops.alt= [693]; % altitude in km (at epogee)
LOStops.IW1inc=[30.86]; LOStops.IW3inc=[41.75];  % values used are for lowest altitude
[Dnew,Dbox]=GTdef_quadtree_vLOS(D,X,Y,ndown,dDmax,mpd,method,dimwgt,LOStops,OFpre);


% should write a wrapper to recalculate the LOS that corrects for SAT position
