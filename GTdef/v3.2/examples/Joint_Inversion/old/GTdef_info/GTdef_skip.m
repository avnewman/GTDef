%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GTdef_skip.m                                        %
%        function check whether or not to skip reading input information        %
% INPUT 'string'                                                                %
% OUTPUT logical 0 or 1                                                         %
%   checks if empty, starts with # or %                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ out ] = GTdef_skip(str)
if (strncmpi(str,'#',1)||strncmpi(str,'%',1)||isempty(str))
   out = true();
else
   out = false();
end