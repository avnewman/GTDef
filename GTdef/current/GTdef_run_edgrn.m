function [ ] = GTdef_run_edgrn(finName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             	  GTdef_run_edgrn.m                                     %
% (1) read in input file for GTdef							%
% (2) run GTdef_write_edgrn_input.m to output input file for edgrn                      % 
%                                                                                       %
% INPUT											%
% GTdef input file									%
%                                                                                       %
% OUTPUT										%
% "edgrnfcts" folder contains point source library that should cover the whole region 	%
% Note: no while lines allowed in the input file!!					%
%                                                                                       %
% REFERENCE  										%
% Wang, R., Martin, F. L., & Roth, F. (2003)						%
% Computers & Geosciences, 29(2), 195-207. doi:10.1016/S0098-3004(02)00111-5		%
%                                                                                       %
% first created by Lujia Feng Fri Feb 24 18:10:17 SGT 2012                              %
% modified GTdef_open.m lfeng Tue Jun 14 17:01:14 SGT 2016                              %
% last modified by Lujia Feng Tue Jun 14 17:06:05 SGT 2016                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
[ modspace,earth,...
  flt1,flt2,flt3,flt4,flt5,flt6,flt7,...
  subflt,addon,...
  pnt,los,bsl,prf,grd,...
  sspnt,ssflt1,ssflt2 ] = GTdef_open(finName);
toc

[ ~,basename,~ ] = fileparts(finName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% layered earth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if layered model is used, green function library is built up here
fprintf(1,'\n..... calculating point source library ......\n');
if strcmp(earth.type,'layered')
    fedgrnName = [ basename '_edgrn.inp' ];
    GTdef_write_edgrn_input(fedgrnName,earth.edgrn,earth.layer);
    folderName = 'edgrnfcts';
    % create green function folder if it does not exist
    if ~exist(folderName,'dir'), mkdir(folderName); end
    system(['echo ' fedgrnName ' | ~/code/fortran/edgrn_edcmp/edgrn2.0']);
    %system(['echo ' fedgrnName ' | ./edgrn2.0']);
end
