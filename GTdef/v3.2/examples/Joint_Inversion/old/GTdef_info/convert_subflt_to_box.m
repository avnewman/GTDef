function [  ] = convert_subflt_to_box( fault_param,slip_file,Outfile)
% written by AL Williamson
%
% convert roberto's subfaults to a gmt supported multisegement file
% initial lat and lon are the center of a subfault
%clc; clear all; close all;
%[Subflts]=load('fault_param184_DS.txt');
[Subflts]=load(fault_param);

V=slip_file;


FID=fopen(Outfile,'w');

length(Subflts)
for ii = 1:length(Subflts)
    FLT=Subflts(ii,:);
    
    Lon_Cen=FLT(2); Lat_Cen=FLT(1); Z_Cen=FLT(3)*1000; len=FLT(4)*1000; width=FLT(5)*1000; dip=FLT(6); strike=FLT(7); ds=FLT(9);
    
    [lon_t, lat_t] = Haversine( Lon_Cen, Lat_Cen, strike - 90, (width/2)/1000);
    [lon0,lat0]=Haversine(lon_t, lat_t, strike-180,  (len/2)/1000);   
   
   fprintf(FID,'%s%.4f \n', '> -Z',V(ii));
   %fprintf(FID,'%s \n','>');
   fprintf(FID,'%.3f %.3f \n',lon0,lat0);
    
   [lon1,lat1]=Haversine(lon0,lat0,strike,len/1000);
   fprintf(FID,'%.3f %.3f \n',lon1,lat1);
   
   [lon2,lat2]=Haversine(lon1,lat1,strike+90,width/1000);
   fprintf(FID,'%.3f %.3f \n',lon2,lat2);
   
   [lon3,lat3]=Haversine(lon2,lat2,strike-180,len/1000);
   fprintf(FID,'%.3f %.3f \n',lon3,lat3);
   fprintf(FID,'%.3f %.3f \n',lon0,lat0);

end
fclose(FID);
end