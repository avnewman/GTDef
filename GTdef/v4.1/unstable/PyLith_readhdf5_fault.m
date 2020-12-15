function [ vertices,cells ] = PyLith_readhdf5_fault(fileName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PyLith_readhdf5_fault                                %
% read PyLith output hdf5 files for fault subsets                               %
% Note: HDF5 writes in row major order, Matlab reads in column major order      %
%                                                                               %
% INPUT:                                                                        %
% PyLith fault output hierarchy convention                                      %
% groups are highest; datasets are underneath groups; only datasets readable    %
% output by PyLith utility                                                      %
% ----------------------------------------------------------------------------- %
% h5dump -n MTW_greensfns-fault.h5                                              %
% FILE_CONTENTS {                                                               %
%    group      /                                                               %
%    group      /geometry                                                       %
%    dataset    /geometry/vertices              (vertices on fault)             %
%    dataset    /time                           (time step)                     %
%    group      /topology                       (cell topology with 4 vertices) %
%    dataset    /topology/cells                                                 %
%    group      /vertex_fields                                                  %
%    dataset    /vertex_fields/slip             (slip for vertices)             %
%    dataset    /vertex_fields/traction_change  (traction change for vertices)  %
%    or                                                                         %
%    dataset    /vertex_fields/traction         (traction for vertices)         %
% }                                                                             %
% ----------------------------------------------------------------------------- %
% fileName - hdf5 file for fault patches                                        %
%                                                                               %
% OUTPUT: structures                                                            %
% vertices.num  - number of vertices                                            %
% vertices.step - number of time steps                                          %
% vertices.id   - 1-based id number for vertices                                %
% vertices.loc  - location of vertices in global coordinate system              %
%               = [ x y z ]                                                     %
% vertices.slip - slip of vertices in fault coordinate (stepnum*vertnum*3)      %
%               = [ leftlateral reverse opening ]                               %
% vertices.dT   - traction change (stepnum*vertnum*3)                           %
%               = [ shear-leftlateral shear-updip normal ]                      %
% cells.num     - number of cells                                               %
% cells.id      - 1-based id number for cells                                   %
% cells.topo    - topology of cells (4 vertices)                                %
%                                                                               %
% first created by Lujia Feng Fri Nov 25 10:18:14 SGT 2011                      %
% added field existence check lfeng Sat Feb  2 17:01:43 SGT 2013                %
% last modified by Lujia Feng Sat Feb  2 17:02:16 SGT 2013                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fileName,'file'), error('PyLith_readhdf5_fault ERROR: %s does not exist!',fileName); end
fin = fopen(fileName,'r');
 
% display the file structure
%h5disp(fileName);
%info = h5info(fileName);

% use Matlab high-level access functions
% h5read reads dataset
vertices.loc  = h5read(fileName,'/geometry/vertices');                  % [3*vertexNum]
vertices.loc  = vertices.loc';						% convert to columnwise [vertexNum*3]
vertices.slip = h5read(fileName,'/vertex_fields/slip');	                % prescribed slip [3*vertexNum*timeStep]
vertices.slip = permute(vertices.slip,[3 2 1]);                         % [timeStep*vertexNum*3]
vertices.num  = size(vertices.loc,1);
vertices.id   = [ 1:vertices.num ]';                                    % id starts from 0 as C convention
vertices.step = length(h5read(fileName,'/time'));                       % number of time steps
if hdf5_exists(fileName,'/vertex_fields/traction_change')
   vertices.dT   = h5read(fileName,'/vertex_fields/traction_change');   % traction change [3*vertexNum*timeStep]
   vertices.dT   = permute(vertices.dT,[3 2 1]);                        % [timeStep*vertexNum*3]
end
if hdf5_exists(fileName,'/vertex_fields/traction')
   vertices.T   = h5read(fileName,'/vertex_fields/traction');           % traction [3*vertexNum*timeStep]
   vertices.T   = permute(vertices.T,[3 2 1]);                          % [timeStep*vertexNum*3]
end

cells.topo    = h5read(fileName,'/topology/cells');                     % topology for cells
cells.topo    = cells.topo'+1;                                          % +1 from h5 index to matlab index
cells.num     = size(cells.topo,1);
cells.id      = [ 1:cells.num ]';                                       % id starts from 0 as C convention in hdf5 files
                                                                        % id is 1-based in Matlab
fclose(fin);
