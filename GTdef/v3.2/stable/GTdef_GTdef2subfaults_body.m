function [ ] = GTdef_GTdef2subfaults_body(outFlag,fout,out)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GTdef_GTdef2subfaults_body  				%
%                        write output main information                          %
% added inverse2 lfeng Tue Aug  4 17:08:20 SGT 2015                             %
% last modified by Lujia Feng Tue Aug  4 17:56:32 SGT 2015                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(outFlag,'relax')
   % (1)no (2)slip (3)north (4)east (5)depth (6)length (7)width (8)strike (9)dip (10)rake
   fprintf(fout,'%-8d %14.6f %14.6e %14.6e %14.6e %14.6e %14.6e %12.6f %12.6f %12.6f\n',out');
elseif strcmp(outFlag,'inverse2')
   % (1)midLon (2)midLat (3)midDepth (4)length (5)width (6)rake (7)strike (8)dip (9)slip
   fprintf(fout,'%14.8f %14.8f %14.4f %14.5e %14.4e %12.6f %8.2f %8.2f %14.6f\n',out');
else
   % (1)no (2)dnum (3)snum (4)lon (5)lat (6)z (7)length (8)width (9)strike (10)dip (11)rake (12)slip
   fprintf(fout,'%3d %3d %3d %14.8f %12.8f  %14.6e   %11.5e %12.5e  %10.6f %10.6f %10.6f %12.8f\n',out');
end
fclose(fout);
