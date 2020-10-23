    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read1double.m                               %
%             function to read in one double number from a string               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ value,remain ] = GTdef_read1double(str)

[str1,remain] = strtok(str);
if GTdef_skip(str1)
   error('GTdef_open ERROR: the input file is wrong!');
end
value = str2double(str1);
end