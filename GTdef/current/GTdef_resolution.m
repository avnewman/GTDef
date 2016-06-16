function [] = GTdef_resolution(finName,modspace,flt1,flt2,flt3,flt4,flt5,flt6,subflt,addon,pnt,los,bsl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               GTdef_resolution.m                              % 
%                                                                               %
% (1) Replacates the functionality of GTdef_project, but includes the model     %
%     resolution matrix, by component as the last thee fields.                  %
%     model interdependence is also output (using 2)                            %
%                                                                               %
% INPUT                                                                         %
%   modspace - model structure                                                  %
%   finName  - its preferred to use the same prefix as kappa-specific ones      %
%   flt?     - all fault types                                                  %
%   subflt   - special subfault characterization                                %
% OUTPUT								        %
% (1)  create an output file that contains the model resolution matrix (diags)  %
%       format is the same as GTdef_project, but adds the resolution at end     %
%      files are named '*_kp???_patches_R.out                                   %
% (2)  create an output file that contains the data resolution matrix (diags)   %
%      files are named '*_kp???_patches_D.out                                   %
% (3)  create a similar output but describes the interdependence of model       %
%      parameters for a given model node                                        %
%      files are named '*_kp???_patches_R_??.out                                %
%                                                                               %
% PARAMETERS                                                                    %
%   same as GTdef_project (but also looks for info in modspace.resolflag        %
%                                                                               %
% first created by Andrew Newman Thu May  5 09:21:55 AST 2016                   %
% added fault6 lfeng Thu Jun  2 10:50:44 SGT 2016                               %
% last modified by Andrew Newman Thu May  5 09:21:55 AST 2016                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coord  = modspace.coord;
origin = modspace.origin;
R      = modspace.R;
Rdiag  = diag(full(R));
N      = modspace.N;
Ndiag  = diag(full(N));
res    = modspace.resolflag;
[ ~,basename,~ ] = fileparts(finName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% set origin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the middle point of the region as origin for the local cartesian coordinate
if strcmpi(coord,'geo') || strcmpi(coord,'geo_polyconic')
   if isempty(origin)
      lonlist = []; latlist = [];
      if flt1.num~=0
          lonlist = [ lonlist;flt1.flt(:,1) ]; latlist = [ latlist;flt1.flt(:,2) ];
      end
      if flt2.num~=0
          lonlist = [ lonlist;flt2.flt(:,1) ]; latlist = [ latlist;flt2.flt(:,2) ];
      end
      if flt3.num~=0
          lonlist = [ lonlist;flt3.flt(:,1) ]; latlist = [ latlist;flt3.flt(:,2) ];
      end
      if flt4.num~=0
          lonlist = [ lonlist;flt4.flt(:,1) ]; latlist = [ latlist;flt4.flt(:,2) ];
      end
      lon0 = 0.5*(min(lonlist)+max(lonlist));
      lat0 = 0.5*(min(latlist)+max(latlist));
   else
      lon0 = origin(1); lat0 = origin(2);
   end
elseif strcmpi(coord,'local')~=1
   error('GTdef_project ERROR: coordinate input is wrong!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% strike variations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addon.crt = [];
if addon.strnum~=0
fprintf(1,'\n....... processing strike variations .......\t');
tic
    % convert strike controlling points from geographic to local cartesian coordinate
    switch coord
       case 'geo'
          [sx1,sy1] = LL2ckmd(addon.str(:,1),addon.str(:,2),lon0,lat0,0);
          [sx2,sy2] = LL2ckmd(addon.str(:,3),addon.str(:,4),lon0,lat0,0);
       case 'geo_polyconic'
          [sx1,sy1] = latlon_to_xy(addon.str(:,1),addon.str(:,2),lon0,lat0,0);
          [sx2,sy2] = latlon_to_xy(addon.str(:,3),addon.str(:,4),lon0,lat0,0);
       case 'local'
          sx1 = addon.str(:,1); sy1 = addon.str(:,2);
          sx2 = addon.str(:,3); sy2 = addon.str(:,4);
    end
    addon.crt = [ sx1 sy1 sx2 sy2 addon.str(:,[5 6]) ];
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% point data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnt.crt = [];
if pnt.num~=0
fprintf(1,'\n......... processing the point data .........\t');
tic
    % convert point data from geographic to local cartesian coordinate
    switch coord
       case 'geo'
          [pxx,pyy] = LL2ckmd(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [pxx,pyy] = latlon_to_xy(pnt.loc(:,1),pnt.loc(:,2),lon0,lat0,0);
       case 'local'
          pxx = pnt.loc(:,1); pyy = pnt.loc(:,2);
    end
    pzz = zeros(length(pxx),1);
    pnt.crt = [pxx pyy pzz];                   	% cartesian - 3*n matrix [xx;yy;zz]; it is just Xin

    % output point file
    fpntName = strcat(basename,'_point.out');
    fpnt = fopen(fpntName,'w');
    fprintf(fpnt,'#(1)point (2)3 (3)name (4)lon (5)lat (6)z (7)Ue (8)Un (9)Uv (10)eUe (11)eUn (12)eUv (13)weight (14)fault name (15)Dstr1  (16)Dstr2 (17)Ddip (18)Dvert\n'); 
toc
end


xsect = []; xsectName = {};
prjfltAll = []; 
fltAllName = {}; 
prjflt1  = []; prjflt2 = []; prjflt3 = []; prjflt4 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    switch coord
       case 'geo'
          [x1,y1] = LL2ckmd(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [x1,y1] = latlon_to_xy(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
       case 'local'
          x1 = flt1.flt(:,1); y1 = flt1.flt(:,2);
    end
    newflt1 = [ x1 y1 flt1.flt(:,3:end) ];
    for ii = 1:flt1.num
       Nd = flt1.flt(ii,17); Ns = flt1.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt1.name{ii};
       cflt   = newflt1(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,1,cfname,cflt);
       end
       % cross-section projection
       [ xsect1Name,xsect1 ] = GTdef_xsection(1,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect1 ];
       xsectName = [ xsectName; xsect1Name ];   
       % surface projection
       [ ~,prjflt1,~ ] = GTdef_prjflt1dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:));
       prjfltAll = [ prjfltAll; prjflt1 ];
       name1    = cell(fltNum,1);
       for ii = 1:fltNum, name1{ii} = cfname; end
       fltAllName = [ fltAllName; name1 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2.num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    switch coord
       case 'geo'
          [x2_1,y2_1] = LL2ckmd(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
          [x2_2,y2_2] = LL2ckmd(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
       case 'geo_polyconic'
          [x2_1,y2_1] = latlon_to_xy(flt2.flt(:,1),flt2.flt(:,2),lon0,lat0,0);
          [x2_2,y2_2] = latlon_to_xy(flt2.flt(:,3),flt2.flt(:,4),lon0,lat0,0);
       case 'local'
          x2_1 = flt2.flt(:,1); y2_1 = flt2.flt(:,2);
          x2_2 = flt2.flt(:,3); y2_2 = flt2.flt(:,4);
    end
    newflt2 = [ x2_1 y2_1 x2_2 y2_2 flt2.flt(:,5:end) ];
    for ii = 1:flt2.num
       Nd = flt2.flt(ii,17); Ns = flt2.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt2.name{ii};
       cflt   = newflt2(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);
       % find strikes for the master fault
       strInd = strcmpi(cfname,addon.strname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,2,cfname,cflt,addon.crt(strInd,:));
       end
       % cross-section projection
       [ xsect2Name,xsect2 ] = GTdef_xsection(2,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect2 ];
       xsectName = [ xsectName; xsect2Name ];   
       % surface projection
       [ ~,prjflt2,~ ] = GTdef_prjflt2dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
       prjfltAll = [ prjfltAll; prjflt2 ];
       name2    = cell(fltNum,1);
       for ii = 1:fltNum, name2{ii} = cfname; end
       fltAllName = [ fltAllName; name2 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault3 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt3.num~=0
fprintf(1,'\n.......... processing fault type-3 ..........\t');
tic
    % convert to cartesian coordinate
    switch coord
       case 'geo'
          [x3,y3] = LL2ckmd(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
       case 'geo_polyconic'
          [x3,y3] = latlon_to_xy(flt3.flt(:,1),flt3.flt(:,2),lon0,lat0,0);
       case 'local'
          x3 = flt3.flt(:,1); y3 = flt3.flt(:,2);
    end
    newflt3 = [x3 y3 flt3.flt(:,3:end)];
    for ii = 1:flt3.num
       Nd = flt3.flt(ii,17); Ns = flt3.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt3.name{ii};
       cflt   = newflt3(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,3,cfname,cflt);
       end
       % cross-section projection
       [ xsect3Name,xsect3 ] = GTdef_xsection(3,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect3 ];
       xsectName = [ xsectName; xsect3Name ];   
       % surface projection
       [ ~,prjflt3,~ ] = GTdef_prjflt3dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:));
       prjfltAll = [ prjfltAll; prjflt3 ];
       name3 = cell(fltNum,1);
       for ii = 1:fltNum, name3{ii} = cfname; end
       fltAllName = [ fltAllName; name3 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4.num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    switch coord
       case 'geo'
          [x4_1,y4_1] = LL2ckmd(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
          [x4_2,y4_2] = LL2ckmd(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
       case 'geo_polyconic'
          [x4_1,y4_1] = latlon_to_xy(flt4.flt(:,1),flt4.flt(:,2),lon0,lat0,0);
          [x4_2,y4_2] = latlon_to_xy(flt4.flt(:,3),flt4.flt(:,4),lon0,lat0,0);
       case 'local'
          x4_1 = flt4.flt(:,1); y4_1 = flt4.flt(:,2);
          x4_2 = flt4.flt(:,3); y4_2 = flt4.flt(:,4);
    end
    newflt4 = [x4_1 y4_1 x4_2 y4_2 flt4.flt(:,5:end)];
    for ii = 1:flt4.num
       Nd = flt4.flt(ii,17); Ns = flt4.flt(ii,18); fltNum = Nd*Ns;
       cfname = flt4.name{ii};
       cflt   = newflt4(ii,:);
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       % find dips for the master fault
       dipInd = strcmpi(cfname,addon.dipname);
       % find strikes for the master fault
       strInd = strcmpi(cfname,addon.strname);

       % point projection
       if ~isempty(pnt.crt)
          GTdef_prjpnt(fpnt,pnt,4,cfname,cflt,addon.crt(strInd,:));
       end
       % cross-section projection
       [ xsect4Name,xsect4 ] = GTdef_xsection(4,cfname,cflt,addon.dip(dipInd,:));
       xsect     = [ xsect; xsect4 ];
       xsectName = [ xsectName; xsect4Name ];   
       % surface projection
       [ ~,prjflt4,~ ] = GTdef_prjflt4dif(cflt,subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
       prjfltAll = [ prjfltAll; prjflt4 ];
       name4 = cell(fltNum,1);
       for ii = 1:fltNum, name4{ii} = cfname; end
       fltAllName = [ fltAllName; name4 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault5 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt5.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    for ii = 1:flt5.num
       cfname  = flt5.name{ii};
       cflt    = flt5.flt(ii,:);
       % find geometry file for the fault
       geoname = flt5.geoname{ii};
       % find column names for the fault
       colname = flt5.colname{ii};
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       
       % point & cross-section projections do ont apply to fault5

       % surface projection
       [ ~,prjflt5,~ ] = GTdef_prjflt5(modspace,geoname,colname,cflt,subflt.flt(subInd,:));
       prjfltAll = [ prjfltAll; prjflt5 ];
       fltNum   = size(prjflt5,1);
       name5    = cell(fltNum,1);
       for ii = 1:fltNum, name5{ii} = cfname; end
       fltAllName = [ fltAllName; name5 ];   
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault6 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt6.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    for ii = 1:flt6.num
       cfname  = flt6.name{ii};
       cflt    = flt6.flt(ii,:);
       % find geometry file for the fault
       geoname = flt6.geoname{ii};
       % find column names for the fault
       colname = flt6.colname{ii};
       % find subfaults for the master fault
       subInd = strcmpi(cfname,subflt.name);
       
       % point & cross-section projections do ont apply to fault6

       % surface projection
       [ ~,prjflt6,~ ] = GTdef_prjflt6(modspace,geoname,colname,cflt,subflt.flt(subInd,:));
       prjfltAll = [ prjfltAll; prjflt6 ];
       fltNum   = size(prjflt6,1);
       name6    = cell(fltNum,1);
       for ii = 1:fltNum, name6{ii} = cfname; end
       fltAllName = [ fltAllName; name6 ];   
    end
toc
end
if ~isempty(pnt.crt)  
    fclose(fpnt); 
    fprintf(1,'GTdef_project output %s\n',fpntName); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Surface Model Resolution File
foutName = strcat(basename,'_patches_R.out');
fout     = fopen(foutName,'w');
[rd3rd,~]=size(Rdiag);
rd3rd=rd3rd/3;
if ~isempty(prjfltAll)
    %              (1)  (2)  (3)   (4)   (5)   (6)   (7)   (8)   (9)   (10)  (11)  (12)  (13)  (14)  (15) (16) (17) (18) (19) (20) (21)  (22)
    % prjfltAll = [ dnum snum xtop1 ytop1 ztop1 xbot1 ybot1 zbot1 xbot2 ybot2 zbot2 xtop2 ytop2 ztop2 xctr yctr zctr ss   ds   ts  rake   rs]
    if strcmpi(coord,'geo')
        xx = prjfltAll(:,[3 6 9 12 15]);  yy = prjfltAll(:,[4 7 10 13 16]);
        [lon,lat] = ckm2LLd(xx,yy,lon0,lat0,0);
        newprjfltAll = [ prjfltAll(:,1:2) lon(:,1) lat(:,1) prjfltAll(:,5) lon(:,2) lat(:,2) prjfltAll(:,8) lon(:,3) lat(:,3) prjfltAll(:,11) ...
	                lon(:,4) lat(:,4) prjfltAll(:,14) lon(:,5) lat(:,5) prjfltAll(:,17) prjfltAll(:,18:22) ];
    end
    if strcmpi(coord,'geo_polyconic')
        xx = prjfltAll(:,[3 6 9 12 15]);  yy = prjfltAll(:,[4 7 10 13 16]);
        lon = zeros(size(xx));         lat = zeros(size(yy));
        for ii=1:5
            [lon(:,ii),lat(:,ii)] = xy_to_latlon(xx(:,ii),yy(:,ii),lon0,lat0,0);
        end        
        newprjfltAll = [ prjfltAll(:,1:2) lon(:,1) lat(:,1) prjfltAll(:,5) lon(:,2) lat(:,2) prjfltAll(:,8) lon(:,3) lat(:,3) prjfltAll(:,11) ...
	                lon(:,4) lat(:,4) prjfltAll(:,14) lon(:,5) lat(:,5) prjfltAll(:,17) prjfltAll(:,18:22) ];
    end
    if strcmpi(coord,'local')
        newprjfltAll = prjfltAll;
    end
    fprintf(fout,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)ztop1 (7)xbot1 (8)ybot1 (9)zbot1 (10)xbot2 (11)ybot2 (12)zbot2 (13)xtop2 (14)ytop2 (15)ztop2 (16)xctr (17)yctr (18)zcrt (19)ss[m] (20)ds[m] (21)ts[m] (22)rake[deg] (23)rs[m] (24)Rss (25)Rds  (26)Rts\n');
    [ row,~ ] = size(newprjfltAll);
    for ii =1:row
        name = fltAllName{ii};
        flt  = newprjfltAll(ii,:);
        Rd   = [Rdiag(ii),Rdiag(ii+rd3rd),Rdiag(ii+2*rd3rd)];
        %             1     2   3   4      5      6      7      8      9      10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26	
        fprintf(fout,'%-10s %4d %4d %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %10.5f %10.5f %10.5f %10.5f %10.5f %6.4f %6.4f %6.4f\n',name,flt,Rd);
    end
end

fclose(fout);
fprintf(1,'GTdef_resolution Model Resolution output %s\n',foutName); 

%   Data Resolution File 
foutName = strcat(basename,'_data_R.out');
fout     = fopen(foutName,'w');
prow=0; lrow=0; brow=0;

if ( ~ pnt.num == 0 )
    prow = pnt.num;
    fprintf(fout,'#(1)pnt (2)lon (3)lat (4)elev (5)obs_e[m] (6)obs_n[m] (7)obs_u[m] (8)err_e[m] (9)err_n[m] (10)err_u[m] (11)wgt (12)pred_e[m] (13)pred_n[m] (14)pred_u[m] (15)Rdata\n');
     for ii =1:prow
         nr=Ndiag(ii); 
                %      1     2        3      4      5      6      7        8      9     10       11      12     13     14       15  
        fprintf(fout,'pnt   %12.5f %11.5f %12.2f   %12.5f %12.5f %12.5f   %11.5f %11.5f %11.5f   %6.3f   %12.5f %12.5f %12.5f   %6.4f\n',pnt.loc(ii,:),pnt.obs(ii,:), pnt.err(ii,:), pnt.wgt(ii), pnt.out(ii,[4:6]), nr);
    end
end
if ( ~ los.num == 0 )
    lrow = los.num;
    fprintf(fout,'#(1)los  (2)lon (3)lat (4)elev (5)obs[m] (6)err[m] (7)wgt (8)look_e (9)look_n  (10)look_u (11)pred[m] (12)Rdata\n');
     for ii =1:lrow
         nr=Ndiag(ii+prow);
                %      1     2        3      4      5      6      7      8      9      10     11     12  
        fprintf(fout,'los   %12.5f %11.5f %12.2f   %12.5f %11.5f %6.3f   %7.5f %7.5f %7.3f   %12.5f %6.4f\n',los.loc(ii,:),los.obs(ii), los.err(ii), los.wgt(ii), los.dir(ii,:),los.out(ii,4), nr);
    end
end
% needs to be checked 
if ( ~ bsl.num == 0 )
    brow = bsl.num;
    fprintf(fout,'#(1)bsl  (2)lon1 (3)lat1 (4)elev1   (5)lon2 (6)lat2 (7)elev2  (8)obs[m] (9)err[m] (10)wgt  (11)pred[m] (12)Rdata\n');
     for ii =1:brow
         nr=Ndiag(ii+prow+lrow);
                %      1     2      3      4        5      6      7        8      9      10      11    12  
        fprintf(fout,'bsl   %12.5f %11.5f %12.2f   %12.5f %11.5f %12.2f   %12.5f %11.5f %6.3f   %12.5f %6.4f\n',bsl.loc(ii,:),bsl.obs(ii), bsl.err(ii), bsl.wgt(ii),bsl.out(ii,4), nr);
    end
end

fprintf(1,'GTdef_resolution Data Resolution output %s\n',foutName); 

% Create a series of Resolution Spread files, showing the model-dependency
% (row parameters of the Resolution Matrix, for a given model. Currently,
% this is identified by a specific subfault number as described by a value
% within the Ms vector 

R=full(R); nres=0;
if (isstr(res))
    if strcmp(res,'all')
        nres=rd3rd;
        useres=0;
    end
else
    nres=length(res);
    useres=1;
end
if (nres>0)
	for jj=1:nres
           if(useres==1)    % get nodes from list
  	         pt=res(jj);    % get all nodes
           else
             pt=jj;
           end
           foutName = strcat(basename,'_patches_R_',num2str(pt),'.out');
           fout     = fopen(foutName,'w');
           fprintf(fout,'#(1)name (2)dnum (3)snum (4)xtop1 (5)ytop1 (6)ztop1 (7)xbot1 (8)ybot1 (9)zbot1 (10)xbot2 (11)ybot2 (12)zbot2 (13)xtop2 (14)ytop2 (15)ztop2 (16)xctr (17)yctr (18)zcrt (19)ss[m] (20)ds[m] (21)ts[m] (22)rake[deg] (23)rs[m] (24)Rss (25)Rds  (26)Rts\n');
           [ row,~ ] = size(newprjfltAll);
           for ii =1:row
               name = fltAllName{ii};
               flt  = newprjfltAll(ii,:);
               Rd   = [R(pt,ii),R(pt+rd3rd,ii+rd3rd),R(pt+2*rd3rd,ii+2*rd3rd)]; % version outputs only the spread of the resolution across like slip-direction (Rss only reports Rss, etc..)
               Rd_ss   = [R(pt,ii),R(pt,ii+rd3rd),R(pt,ii+2*rd3rd)];  % full matrix spread relative to ss patch component.
               Rd_ds   = [R(pt+rd3rd,ii),R(pt+rd3rd,ii+rd3rd),R(pt+rd3rd,ii+2*rd3rd)];  % full matrix spread relative to ds patch component.
               Rd_ts   = [R(pt+2*rd3rd,ii),R(pt+2*rd3rd,ii+rd3rd),R(pt+2*rd3rd,ii+2*rd3rd)];  % full matrix spread relative to ts patch component.
               %             1     2   3   4      5      6      7      8      9      10     11     12     13     14     15     16     17     18     19     20     21     22     23     24     25     26	
               fprintf(fout,'%-10s %4d %4d %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %12.5f %11.5f %12.3e %10.5f %10.5f %10.5f %10.5f %10.5f %6.4f %6.4f %6.4f\n',name,flt,full(Rd));
           end
           fclose(fout);
           fprintf(1,'GTdef_resolution Model Resolution output %s\n',foutName); 
	end
end
