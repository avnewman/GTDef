function [ matrix1,ind ] = GTdef_condense(matrix0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_condense				  %
% Condense the smoothing matrix by removing rows of all zeros		  %
% INPUT:					  		  	  %
%   matrix0 - original large matrix					  %
%									  %
% OUTPUT:								  %
%   matrix1 - new small condensed matrix				  %
%   ind     - index of rows that are kept                                 %
%									  %
% first created by Lujia Feng Thu Dec  3 01:13:59 EST 2009		  %
% output ind lfeng Thu Jun 16 16:44:47 SGT 2016                           %
% last modified by Lujia Feng Thu Dec  3 01:17:03 EST 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ row col ] = size(matrix0);
zero_row = zeros(1,col);
ind = [];
for ii = 1:row
    if ~isequal(matrix0(ii,:),zero_row)
        ind = [ ind ii ];
    end
end
matrix1 = matrix0(ind,:);
