function [ vertices,cells,gpsites,siteList,grnfns ] = PyLith_greensfns(fltName,pntName,siteName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               PyLith_greensfns                                %
% read PyLith output hdf5 files for fault and gps points                        %
% to compile and reorder a green function database                              %
% order patches along strike first then along dip                               %
% order sites according to file order in siteName                               %
%                                                                               %
% INPUT:                                                                        %
% PyLith fault output hierarchy convention                                      %
% groups are highest; datasets are underneath groups; only datasets readable    %
% output by PyLith utility                                                      %
% h5dump -n filename: dump the hierarchy of an hdf5 file                        %
% h5dump -H filename: dump the hierarchy with dataset dimensions and attributes %
% ----------------------------------------------------------------------------- %
% fltName  - hdf5 file for fault patches                                        %
% pntName  - hdf5 file for discrete points                                      %
% siteName - site location file that corresponds to spatialdb                   %
% 1      2          3           4             5                6                %
% Site	 XX         YY 	        Height	      Lon              Lat              %
% ABGS   131.54919  209.68823   -0.00100      99.387520914     0.220824642      %
% BITI    42.68039  388.48297   -0.00100      97.811362862     1.078602706      %
%                                                                               %
% STRUCTURES:                                                                   %
% vertices.num  - number of vertices               (scalar)                     %
% vertices.step - number of time steps             (scalar)                     %
% vertices.id   - 1-based id number for vertices   (vertexNum*1)                %
% vertices.loc  - location of vertices in global coordinate system              %
%               = [ x y z ]                        (vertexNum*3)                %
% vertices.slip - slip of vertices in fault coordinate                          %
%               = [ leftlateral reverse opening ]  (timeStep*vertexNum*3)       %
% vertices.dT   - traction change                  (timeStep*vertexNum*3)       %
%               = [ shear-leftlateral shear-updip normal ]                      %
% cells.num     - number of cells                  (scalar)                     %
% cells.id      - 1-based id number for cells      (vertexNum*1)                %
% cells.topo    - topology of cells (4 vertices)   (cellNum*4)                  %
% gpsites.num   - number of gps sites              (scalar)                     %
% gpsites.step  - number of time steps             (scalar)                     %
% gpsites.id    - 1-based id number for gps sites (PyLith randomizes site order)%
% gpsites.loc   - location of gps sites in global coordinate system             %
%               = [ x y z ]                        (pntNum*3)                   %
% gpsites.disp  - greens function for gps sites                                 %
%               = [ dx dy dz ]                     (timeStep*pntNum*3)          %
% siteList      - site names stored as a cell      (pntNum*1)                   %
% ----------------------------------------------------------------------------- %
%                                                                               %
% OUTPUT:                                                                       %
% grnfns        - array of size                    (vertexNum*pntNum*9)         %
%               for each patch-site pair                                        % 
%               = [ ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]     %
%                                                                               %
% first created by lfeng Thu Nov 29 10:11:17 SGT 2012                           %
% last modified by lfeng Fri Nov 30 12:35:35 SGT 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ vertices,cells ] = PyLith_readhdf5_fault(fltName);
[ gpsites ]        = PyLith_readhdf5_pnts(pntName);
[ siteList,loc ]   = GPS_readsites(siteName);

% find Nd & Ns using number of vertices and cells
% solve equation ax^2+bx+c=0 pp = [ a b c ]
% Nd^2+Nd(1+cellNum-vertNum)+cellNum=0
pp   = [ 1 1+cells.num-vertices.num cells.num ];
sols = roots(pp);
% assume strike vertices > dip vertices
Nd   = min(sols)+1;   % without 1 is number of elements
Ns   = max(sols)+1;
vertices.Nd = Nd;
vertices.Ns = Ns;

% reorder patches along strike (y) increasingly along dip (x) increasingly first then 
newloc       = vertices.loc*1e-3;                            % from m to km
newloc(:,2)  = round(0.1*(newloc(:,2)))*10;                  % change y to closest 10,20,30,...,90
[ ~,ind ]    = sortrows(newloc,[2 1]);
vertices.loc = vertices.loc(ind,:);
vertices.id  = vertices.id(ind);

% loop through all vertices to find their corresponding time step
slip = round(vertices.slip);                                 % very important to round to closest integer!!!!!
vertices.time = zeros(vertices.num,3);                       % 0 means no corresponding time step [ ss ds ts ]
for ii=1:vertices.num
    id = vertices.id(ii);
    % strike-slip
    timeInd = find(slip(:,id,1)==1);
    if isscalar(timeInd)      % should be only one value = 1
        vertices.time(ii,1) = timeInd;
    elseif ~isempty(timeInd)
        num = length(timeInd);
        error('PyLith_greensfns ERROR: vertex %d have %f strike-slip Greens functions!',id,num);
    end
    % dip-slip
    timeInd = find(slip(:,id,2)==1);
    if isscalar(timeInd)      % should be only one value = 1
        vertices.time(ii,2) = timeInd;
    elseif ~isempty(timeInd)
        num = length(timeInd);
        error('PyLith_greensfns ERROR: vertex %d have %f dip-slip Greens functions!',id,num);
    end
    % opening
    timeInd = find(slip(:,id,3)==1);
    if isscalar(timeInd)      % should be only one value = 1
        vertices.time(ii,3) = timeInd;
    elseif ~isempty(timeInd)
        num = length(timeInd);
        error('PyLith_greensfns ERROR: vertex %d have %f opening Greens functions!',id,num);
    end
end

% check greens functions
vertInd = find(vertices.time(:,1)>0 & vertices.time(:,2)>0);
vertNum = length(vertInd);
if vertNum ~= vertices.num
    error('PyLith_greensfns ERROR: Greens function vertices number %d mismatches vertices number %d!',...
           vertNum,vertices.num);
end

%%----------------------- visual check -----------------------
%-----xx = vertices.loc(vertInd,1);
%-----yy = vertices.loc(vertInd,2);
%-----zz = vertices.loc(vertInd,3);
%-----plot3(xx,yy,zz);

% reorder gps sites according to PyLith input because PyLith randomizes them 
newInd = zeros(gpsites.num,1);
for ii=1:gpsites.num
    xx = loc(ii,1); yy = loc(ii,2);
    dist = sqrt((gpsites.loc(:,1)-xx).^2+(gpsites.loc(:,2)-yy).^2);
    [ ~,ind ] = sort(dist);
    newInd(ii) = ind(1); 
end
gpsites.loc = gpsites.loc(newInd,:);
gpsites.id  = gpsites.id(newInd);

% reorder greens function field [ dx dy dz ] (vertNum*pntNum*(strike-slip dip-slip opening)
grnfns = zeros(vertices.num,gpsites.num,9);
for ii=1:vertices.num
    % strike-slip
    step = vertices.time(ii,1);
    if step>0
       for jj=1:gpsites.num
	   siteId = gpsites.id(jj);
           grnfns(ii,jj,1) = gpsites.disp(step,siteId,1);  % xx
           grnfns(ii,jj,2) = gpsites.disp(step,siteId,2);  % yy
           grnfns(ii,jj,3) = gpsites.disp(step,siteId,3);  % zz
       end
    end
    % dip-slip
    step = vertices.time(ii,2);
    if step>0
       for jj=1:gpsites.num
	   siteId = gpsites.id(jj);
           grnfns(ii,jj,4) = gpsites.disp(step,siteId,1);  % xx
           grnfns(ii,jj,5) = gpsites.disp(step,siteId,2);  % yy
           grnfns(ii,jj,6) = gpsites.disp(step,siteId,3);  % zz
       end
    end
    % dip-slip
    step = vertices.time(ii,3);
    if step>0
       for jj=1:gpsites.num
	   siteId = gpsites.id(jj);
           grnfns(ii,jj,7) = gpsites.disp(step,siteId,1);  % xx
           grnfns(ii,jj,8) = gpsites.disp(step,siteId,2);  % yy
           grnfns(ii,jj,9) = gpsites.disp(step,siteId,3);  % zz
       end
    end
end

%%----------------------- visual check -----------------------
%siteId = 20;
%xx = vertices.loc(:,1);
%yy = vertices.loc(:,2);
%% dip-slip
%%zz = sqrt(grnfns(:,siteId,1).^2 + grnfns(:,siteId,2).^2); % horizontal
%%zz = grnfns(:,siteId,3);                                  % vertical
%% dip-slip
%zz = sqrt(grnfns(:,siteId,4).^2 + grnfns(:,siteId,5).^2);  % horizontal
%%zz = grnfns(:,siteId,6);                                  % vertical
%plot3(xx,yy,zz,'ob'); hold on
%pntx = gpsites.loc(siteId,1);
%pnty = gpsites.loc(siteId,2);
%plot3(pntx,pnty,max(zz),'or');
