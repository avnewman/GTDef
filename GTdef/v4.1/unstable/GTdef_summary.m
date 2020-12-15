function [] = GTdef_summary(filename,modspace)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_summary.m				        %
% 		  function to output the summary of models			%
%										%
% INPUT:									%
% modspace - model structure                                                    %
%										%
% OUTPUT:                                                                       %
% an output file                                                                %
%										%
% first created by Lujia Feng Mon May 11 15:25:46 EDT 2009		        %
% added 1st derivative roughness lfeng Thu Dec  3 01:33:11 EST 2009		%
% added units in the output lfeng Wed Jul 21 17:07:56 EDT 2010			%
% added modspace structure lfeng Fri Mar 20 20:34:08 SGT 2015                   %
% added output filename lfeng Tue Jun 23 18:53:50 SGT 2015                      %
% last modified by Lujia Feng Tue Jun 23 18:54:30 SGT 2015                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen(filename,'w');

modSummary = [ modspace.beta; modspace.kappa; modspace.modinfo' ];
fprintf(fout,'#(1)beta (2)kappa (3)data_num (4)slip_num (5)ndf (6)rss [m^2] (7)rms [m] (8)wrrs [m^2] (9)wrms [m] (10)chi2 (11)rchi2 (12)r_1d [cm/km] (13)r_2d [cm/km^2] (14)strain [cm/km]\n');
fprintf(fout,'%-10.5e  %-16.5e %d %-d %-d %-12.5e %-12.5e  %-12.5e  %-12.5e  %-12.5e  %-12.5e  %-12.5e  %-12.5e  %-12.5e\n',modSummary);
fclose(fout);

fprintf(1,'GTdef_summary output %s\n',filename);
