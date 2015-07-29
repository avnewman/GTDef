function [] = GTdef_summary(filename,kappa,mod_info)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_summary.m				        %
% 		  function to output the summary of models			%
%										%
% INPUT:									%
% kappa    - a row vector of kappa values                                       %
% mod_info - the inverted model corresponding to each kappa			%
%          = [ data_num slip_num ndf rss rms wrrs wrms chi2 rchi2 roughness ]	%
%										%
% OUTPUT: an output file                                                        %
%										%
% first created by Lujia Feng Mon May 11 15:25:46 EDT 2009		        %
% last modified by Lujia Feng Mon May 11 15:30:48 EDT 2009			%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen(filename,'w');

mod_summary = [ kappa;mod_info' ];
fprintf(fout,'#kappa data_num slip_num ndf rss rms wrrs wrms chi2 rchi2 roughness\n');
fprintf(fout,'%-10.5f  %d %-d %-d     %-12.5e  %-12.5e  %-12.5e  %-12.5e  %-12.5e  %-12.5e  %-12.5e\n',mod_summary);

fclose(fout);
