function [ siteloc,grnList,grnfns ] = PyLith_trim_greensfns(usedList,siteList,siteloc,grnList,grnfns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PyLith_trim_greensfns                             %
% trim PyLith greensfns related data using usedList                             %
% issue warnings if any sites can not be found in greensfns                     %
%                                                                               %
% INPUT:                                                                        %
% usedList - used site names stored as a cell     (pntNum*1)                    %
% siteList - site names stored as a cell          (pntNum*1)                    %
% siteloc  - [ xx yy zz ]                         (pntNum*3)                    %
% vertices - [ id dnum snum xx yy zz ]            (vertNum*6)                   %
% grnList  - site names corresponding to grnfns   (vertexNum*pntNum*1)          %
% grnfns   - array of size                        (vertexNum*pntNum*9)          %
%            for each patch-site pair                                           % 
%          = [ ss_dx ss_dy ss_dz ds_dx ds_dy ds_dz ts_dx ts_dy ts_dz ]          %
%                                                                               %
% OUTPUT:                                                                       %
%                                                                               %
% first created by lfeng Fri Nov 30 17:15:19 SGT 2012                           %
% last modified by lfeng Fri Nov 30 18:35:27 SGT 2012                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usedNum = length(usedList);
usedInd = zeros(usedNum,1);
for ii=1:usedNum
    site = usedList{ii};
    ind  = find(strcmpi(site,siteList)==true);
    if isscalar(ind)      % should be only one value = 1
        usedInd(ii) = ind;
    elseif isempty(ind)
        error('PyLith_trim_greensfns ERROR: site %s is not in the Greens function database!',site);
    else
        error('PyLith_trim_greensfns ERROR: site %s exists more than once in the Greens function database!',site);
    end
end

siteNum = size(siteloc,1);
siteloc = siteloc(usedInd,:);
grnList = reshape(grnList,siteNum,[]);
grnList = grnList(usedInd,:);
grnList = reshape(grnList,[],1);
grnfns  = grnfns(:,usedInd,:);
