function [ ] = GTdef_write_edgrn_input(finp_name,edgrn,layer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           GTdef_write_edgrn_input.m					%
% prepare input file for fortran code EDGRN by Wang et al. 2003				%
% need to first compile fortran code EDGRN as executable edgrn2.0			%
%											%
% INPUT											%
% Layered earth structure:								%
%   edgrn  structure                                                    		%	
%   layer - [ id depth vp vs ro ]	(nn*5)						%
%											%
% OUTPUT										%
% Point source library should cover the whole model region 				%
% Note: no while lines allowed in the input file!!					%
%       Fortran uses 0d+3 'd' for double presion 'e' for single precision		%
%       Using either 'd' or 'e' produced the same results				%
%											%
% REFERENCE  										%
% Wang, R., Martin, F. L., & Roth, F. (2003)						%
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		%
%		                                                                	%
% first created by lfeng Thu Feb 23 16:10:39 SGT 2012					%
% modified based Emma Hill's make_edgrn_input_file.m lfeng Thu Feb 23 16:45:14 SGT 2012 %
% last modified by lfeng Mon Feb 27 02:23:03 SGT 2012					%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen(finp_name,'w');

fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'#     PARAMETERS FOR THE OBSERVATION PROFILE\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'# obs_depth - uniform depth of observation points [m]\n');
fprintf(fout,'# nr        - number of equidistant radial distances\n');
fprintf(fout,'# r1,r2     - minimum and maximum distances [m]\n');
fprintf(fout,'# nzs       - number of equidistant source depths\n');
fprintf(fout,'# zs1,zs2   - minimum and maximum source depths [m]\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'%-14.5e \t\t\t\t\t|dble: obs_depth;\n',edgrn.obsz);
fprintf(fout,'%-8.0f   %8.2e   %14.5e \t\t|int: nr; dble: r1, r2;\n',   edgrn.nr,edgrn.minr,edgrn.maxr);
fprintf(fout,'%-8.0f   %8.2e   %14.5e \t\t|int: nzs; dble: zs1, zs2;\n',edgrn.nz,edgrn.minz,edgrn.maxz);
fprintf(fout,'#\n');

fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'#     WAVENUMBER INTEGRATION PARAMETERS\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'# srate     - sampling rate for wavenumber integration [10-128]\n');
fprintf(fout,'# The larger the value is, the more accurate the results will be\n');
fprintf(fout,'# More computation time will be required)\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'%.0f \t\t\t\t\t\t|dble: srate;\n',edgrn.srate);
fprintf(fout,'#\n');

fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'#     OUTPUT FILES\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'# outputs - output directory\n');
fprintf(fout,'# grnfile - three files for fundamental Greens functions\n');
fprintf(fout,'# length of names < 80 characters\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'''./edgrnfcts/''  ''izmhs.ss''  ''izmhs.ds''  ''izmhs.cl'' \t|char: outputs,grnfile(3);\n');
fprintf(fout,'#\n');

fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'#     MULTILAYERED MODEL PARAMETERS\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'# no_model_lines - total number of the data lines\n');
fprintf(fout,'# then a table for the layered model\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'%.0f \t\t\t\t\t\t|int: no_model_lines;\n',edgrn.nl);
fprintf(fout,'#\n');

fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'# no   depth[m]       vp[m/s]         vs[m/s]        ro[kg/m^3]\n');
fprintf(fout,'#------------------------------------------------------------------------------\n');
fprintf(fout,'%.0f   %12.4e  %-12.4e  %-12.4e  %-12.4e\n',layer');
