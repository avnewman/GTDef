function [ gpsites ] = PyLith_readhdf5_pnts(fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PyLith_readhdf5_pnts                                 %
% read PyLith output hdf5 files for gps points                                  %
%                                                                               %
% PyLith fault output hierarchy convention                                      %
% groups are highest; datasets are underneath groups; only datasets readable    %
% output by PyLith utility                                                      %
% ----------------------------------------------------------------------------- %
% h5dump -n MTW_greensfns-gpspoints.h5                                          %
% FILE_CONTENTS {                                                               %
%  group      /                                                                 %
%  group      /geometry                                                         %
%  dataset    /geometry/vertices                (vertices here for gps sites)   %
%  dataset    /time                             (time step = vertices id)       %
%  group      /topology                                                         %
%  dataset    /topology/cells                   (gps sites id)                  %
%  group      /vertex_fields                                                    %
%  dataset    /vertex_fields/displacement       (Greens' functions)             %
% }                                                                             %
% ----------------------------------------------------------------------------- %
% INPUT:                                                                        %
% fileName - hdf5 file for discrete points                                      %
%                                                                               %
% OUTPUT: structures                                                            %
% gpsites.num   - number of gps sites                                           %
% gpsites.step  - number of time steps                                          %
% gpsites.id    - 1-based id number for gps sites (PyLith randomizes site order)%
% gpsites.loc   - location of gps sites in global coordinate system             %
%               = [ x y z ]                                                     %
% gpsites.disp  - greens function for gps sites                                 %
%               = [ dx dy dz ]                                                  %
%                                                                               %
% first created by Lujia Feng Wed Nov 28 22:38:36 SGT 2012                      %
% last modified by Lujia Feng Thu Nov 29 12:08:59 SGT 2012                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fileName,'file'), error('PyLith_readhdf5_pnts ERROR: %s does not exist!',fileName); end
fin = fopen(fileName,'r');

% display the file structure
%h5disp(fileName);

% use Matlab high-level access functions
% h5read reads dataset
gpsites.id   = h5read(fileName,'/topology/cells');                      % gps sites id
gpsites.id   = gpsites.id'+1;                                           % +1 from h5 index to matlab index
gpsites.loc  = h5read(fileName,'/geometry/vertices');                   % [3*vertexNum]
gpsites.loc  = gpsites.loc';				                % convert to columnwise [vertexNum*3]
gpsites.disp = h5read(fileName,'/vertex_fields/displacement');	        % displacement [3*vertexNum*timeStep]
gpsites.disp = permute(gpsites.disp,[3 2 1]);                           % [timeStep*vertexNum*3]
gpsites.num  = size(gpsites.loc,1);
gpsites.step = length(h5read(fileName,'/time'));                        % number of time steps

fclose(fin);
