function [ fout ] = GTdef_GTdef2subfaults_head(outFlag,basename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_GTdef2subfaults_head  				%
%                        write output header information                        %
% added inverse2 lfeng Tue Aug  4 17:08:20 SGT 2015                             %
% last modified by Lujia Feng Wed Aug  5 01:13:28 SGT 2015                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(outFlag,'relax')
   foutName = [ basename '_subfaults.flt' ];
   fprintf(1,'The output file is %s\n',foutName);
   fout = fopen(foutName,'w');
   fprintf(fout,'#------------------------------------------------------------------------------------------------\n');
   fprintf(fout,'# (1)no (2)slip (3)north (4)east (5)depth (6)length (7)width (8)strike (9)dip (10)rake\n');
   fprintf(fout,'#          [m]     [m]      [m]     [m]       [m]      [m]      [deg]   [deg]     [deg]\n');
   fprintf(fout,'#------------------------------------------------------------------------------------------------\n');
elseif strcmp(outFlag,'inverse2')
   foutName = [ basename '_subfaults.input' ];
   fprintf(1,'The output file is %s\n',foutName);
   fout = fopen(foutName,'w');
   fprintf(fout,'#------------------------------------------------------------------------------------------------\n');
   fprintf(fout,'# % (1)midLon (2)midLat (3)midDepth (4)length (5)width (6)rake (7)strike (8)dip (9)slip\n');
   fprintf(fout,'#      [deg]     [deg]       [m]        [m]       [m]    [deg]    [deg]    [deg]   [m]\n');
   fprintf(fout,'#------------------------------------------------------------------------------------------------\n');
else
   foutName = [ basename '_subfaults.out' ];
   fprintf(1,'Output file is %s\n',foutName);
   fout = fopen(foutName,'w');
   fprintf(fout,'#------------------------------------------------------------------------------------------------\n');
   fprintf(fout,'# (1)no (2)dnum (3)snum (4)lon (5)lat (6)z (7)length (8)width (9)strike (10)dip (11)rake (12)slip\n');
   fprintf(fout,'#                         [deg]  [deg]  [m]    [m]       [m]     [deg]     [deg]   [deg]     [m]\n');
   fprintf(fout,'#------------------------------------------------------------------------------------------------\n');
end
