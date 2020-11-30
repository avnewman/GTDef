function [ ddip,dlen ] = PyLith_dist_greensfns(vertices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PyLith_dist_greensfns                             %
% determine average distances along dip and strike                              %
% for greens functions database                                                 %
%                                                                               %
% INPUT:                                                                        %
% vertices - [ id dnum snum xx yy zz ]            (vertNum*6)                   %
%                                                                               %
% OUTPUT:                                                                       %
% ddip - average distance along dip between two vertices	                %
% dlen - distance along strike between two vertices                             %
%                                                                               %
% first created by lfeng Fri Nov 30 19:19:14 SGT 2012                           %
% last modified by lfeng Fri Nov 30 19:56:07 SGT 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns = max(vertices(:,3));
Nd = max(vertices(:,2));
% m->km
%vertices(:,4:6) = vertices(:,4:6)*1e-3;

% calculate along dip average distance
dipMat = [];
for ii=1:Ns
    ind  = vertices(:,3)==ii;
    vert = vertices(ind,:);
    xx1  = vert(1:end-1,4);
    yy1  = vert(1:end-1,5);
    zz1  = vert(1:end-1,6);
    xx2  = vert(2:end,4);
    yy2  = vert(2:end,5);
    zz2  = vert(2:end,6);
    dist = sqrt((xx1-xx2).^2+(yy1-yy2).^2+(zz1-zz2).^2);
    dipMat = [ dipMat dist ];
end
ddip = mean(reshape(dipMat,[],1));

%%--------------- visual check ---------------
%figure;
%[X,Y] = meshgrid(1:Ns,1:Nd-1);
%surf(X,Y,dipMat);

% calculate along strike average distance
strMat = [];
for ii=1:Nd
    ind  = vertices(:,2)==ii;
    vert = vertices(ind,:);
    xx1  = vert(1:end-1,4);
    yy1  = vert(1:end-1,5);
    zz1  = vert(1:end-1,6);
    xx2  = vert(2:end,4);
    yy2  = vert(2:end,5);
    zz2  = vert(2:end,6);
    dist = sqrt((xx1-xx2).^2+(yy1-yy2).^2+(zz1-zz2).^2);
    strMat = [ strMat; dist' ];
end
dlen = mean(reshape(strMat,[],1));

%%--------------- visual check ---------------
%figure;
%[X,Y] = meshgrid(1:Ns-1,1:Nd);
%surf(X,Y,strMat);
