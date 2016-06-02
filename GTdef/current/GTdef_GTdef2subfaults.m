function [ ] = GTdef_GTdef2subfaults(finName,outFlag)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         GTdef_GTdef2subfaults                        %
%                                                                      %
% Convert GTdef format to lists of subfaults                           %
% only fault2 has been verified!                                       %
%                                                                      %
% INPUT:                                                               %
% finName - GTdef input/output file                                    %
% outFlag - output format                                              %
% -------------------------------------------------------------------- %
% (1) 'GTdef_center' or 'GTdef_topleft' uses geographic coordinates    %
% 1  2    3    4   5   6 7      8     9      10  11   12               %
% no dnum snum lon lat z length width strike dip rake slip             %
% 'GTdef_center' uses center point                                     %
% 'GTdef_topleft' uses top left point                                  %
%                                                                      %
% (2) 'relax' uses local coordinate                                    %
% 1  2    3     4    5     6      7     8      9   10                  %
% no slip north east depth length width strike dip rake                %
% 'relax' always uses top left point                                   %
%                                                                      %
% (3) 'inverse2' uses geographic coordinates and center points         % 
% 1      2      3        4      5     6    7      8   9                %
% midLon midLat midDepth length width rake strike dip slip(cm)         %    
%                                                                      %
% This is for geometry not for slip                                    %
% 1      2      3        4      5   6      7                           %
% midLon midLat midDepth strike dip length width                       %
% -------------------------------------------------------------------- %
%                                                                      %
% OUTPUT:                                                              %
% a file ended with '_subfaults.out' or '_subfaults.flt'               %
%                                                                      %
% first created by Lujia Feng Mon Jul  9 11:26:06 SGT 2012             %
% added origin Lujia Feng Thu Dec  5 21:55:50 SGT 2013                 %
% added strike lfeng Wed Nov 12 18:45:05 SGT 2014                      %
% added modspace & xyzflt lfeng Tue Mar 24 11:21:45 SGT 2015           %
% changed pntpos to outFlag lfeng Thu Mar 26 18:30:32 SGT 2015         %
% modified for multiple faults lfeng Thu Mar 26 19:00:01 SGT 2015      %
% added inverse2 output lfeng Tue Aug  4 17:31:04 SGT 2015             %
% last modified by Lujia Feng Tue Aug  4 17:53:43 SGT 2015             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read in %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'.......... reading the input file ...........\t');
tic
[ modspace,~,...
  flt1,flt2,flt3,flt4,flt5,flt6,flt7,...
  subflt,addon,...
  ~,~,~,~,~,...
  ~,~,~ ] = GTdef_open(finName);
toc

coord  = modspace.coord;
origin = modspace.origin;

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
    error('GTdef_GTdef2subfaults ERROR: coordinate input is wrong!!!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% strike variations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if addon.strnum~=0
fprintf(1,'\n....... processing strike variations .......\t');
tic
    % convert strike controlling points from geographic to local cartesian coordinate
    switch modspace.coord
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
else
    addon.crt = addon.str;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault1 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt1.num~=0
fprintf(1,'\n.......... processing fault type-1 ..........\t');
tic
    % convert to cartesian coordinate
    switch coord
        case 'geo'
           [x1,y1] = LL2ckmd(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
        case 'geo_polyconic'
           [x1,y1] = latlon_to_xy(flt1.flt(:,1),flt1.flt(:,2),lon0,lat0,0);
        case 'local'
           x1 = flt1.flt(:,1); y1 = flt1.flt(:,2);
    end
    tmpflt = [x1 y1 flt1.flt(:,3:end)];
    for ii = 1:flt1.num
        Nd      = flt1.flt(ii,17); 
	Ns      = flt1.flt(ii,18); 
	fltnum  = Nd*Ns;
        fltname = flt1.name{ii};
        % find the subfaults for the master fault
        subInd = strcmpi(fltname,subflt.name);
        % find dips for the master fault
        dipInd = strcmpi(fltname,addon.dipname);
	% project faults
        [ newflt,~,xyzflt ] = GTdef_prjflt1dif(tmpflt(ii,:),subflt.flt(subInd,:),addon.dip(dipInd,:));
        %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
        % newflt = [ dnum snum x1  y1  z1  z2  len str dip slips ];
        % slips  = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX ];
        dnum  = newflt(:,1);  snum = newflt(:,2);
        x1    = newflt(:,3);  y1   = newflt(:,4);
        z1    = newflt(:,5);  z2   = newflt(:,6);
        len   = newflt(:,7);  str  = newflt(:,8);
        dp    = newflt(:,9);
        ss    = newflt(:,10); ds   = newflt(:,11);
        width = (z2-z1)./abs(sind(dp));
        slip  = sqrt(ss.^2+ds.^2);
        rake  = atan2(ds,ss).*180/pi;
	subfltnum = size(newflt,1);
	inum  = [ 1:subfltnum ]';
        % different output formats
        if strfind(outFlag,'GTdef')
	   if strfind(outFlag,'center') 
              xx    =  xyzflt.xyzctr(:,1);  
              yy    =  xyzflt.xyzctr(:,2); 
              depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           end
           if strfind(outFlag,'topleft')
              xx    =  xyzflt.xyztop1(:,1);  
              yy    =  xyzflt.xyztop1(:,2); 
              depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           end
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ inum dnum snum lon lat depth len width str dp rake slip ];
        elseif strcmp(outFlag,'relax')
           east  =  xyzflt.xyztop1(:,1);  
           north =  xyzflt.xyztop1(:,2); 
           depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           out = [ inum slip north east depth len width str dp rake ];
        elseif strcmp(outFlag,'inverse2')
           xx    =  xyzflt.xyzctr(:,1);  
           yy    =  xyzflt.xyzctr(:,2); 
           depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ lon lat depth len width rake str dp slip ];
        else
           error('GTdef_GTdef2subfaults ERROR: output format is wrong!!!');
        end

        % output to a file
        [ fout ] = GTdef_GTdef2subfaults_head(outFlag,fltname);
        GTdef_GTdef2subfaults_body(outFlag,fout,out)
    end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault2 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt2.num~=0
fprintf(1,'\n.......... processing fault type-2 ..........\t');
tic
    % convert to cartesian coordinate
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
    tmpflt = [ x2_1 y2_1 x2_2 y2_2 flt2.flt(:,5:end) ];
    for ii = 1:flt2.num
        Nd      = flt2.flt(ii,17); 
	Ns      = flt2.flt(ii,18); 
	fltnum  = Nd*Ns;
        fltname = flt2.name{ii};
        % find the subfaults for the master fault
        subInd = strcmpi(fltname,subflt.name);
        % find dips for the master fault
        dipInd = strcmpi(fltname,addon.dipname);
        % find strikes for the master fault
        strInd = strcmpi(fltname,addon.strname);
	% project faults
        [ newflt,~,xyzflt ] = GTdef_prjflt2dif(tmpflt(ii,:),subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
        %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
        % newflt = [ dnum snum x1  y1  x2  y2  z1  z2  dip slips ];
        % slips  = [ ss ds ts ss0 ssX ds0 dsX ts0 tsX ];
        dnum  = newflt(:,1);  snum = newflt(:,2);
        x1    = newflt(:,3);  y1   = newflt(:,4);
        x2    = newflt(:,5);  y2   = newflt(:,6);
        z1    = newflt(:,7);  z2   = newflt(:,8);
        dp    = newflt(:,9);
        ss    = newflt(:,10); ds   = newflt(:,11);
        str   = GTdef_strike(x1,y1,x2,y2);
        len   = sqrt((x2-x1).^2+(y2-y1).^2);
        width = (z2-z1)./abs(sind(dp));
        slip  = sqrt(ss.^2+ds.^2);
        rake  = atan2(ds,ss).*180/pi;
	subfltnum = size(newflt,1);
	inum  = [ 1:subfltnum ]';
        % different output formats
        if strfind(outFlag,'GTdef')
	   if strfind(outFlag,'center') 
              xx    =  xyzflt.xyzctr(:,1);  
              yy    =  xyzflt.xyzctr(:,2); 
              depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           end
           if strfind(outFlag,'topleft')
              xx    =  xyzflt.xyztop1(:,1);  
              yy    =  xyzflt.xyztop1(:,2); 
              depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           end
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ inum dnum snum lon lat depth len width str dp rake slip ];
        elseif strcmp(outFlag,'relax')
           east  =  xyzflt.xyztop1(:,1);  
           north =  xyzflt.xyztop1(:,2); 
           depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           out = [ inum slip north east depth len width str dp rake ];
        elseif strcmp(outFlag,'inverse2')
           xx    =  xyzflt.xyzctr(:,1);  
           yy    =  xyzflt.xyzctr(:,2); 
           depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ lon lat depth len width rake str dp slip ];
        else
           error('GTdef_GTdef2subfaults ERROR: output format is wrong!!!');
        end

        % output to a file
        [ fout ] = GTdef_GTdef2subfaults_head(outFlag,fltname);
        GTdef_GTdef2subfaults_body(outFlag,fout,out)
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
    tmpflt = [x3 y3 flt3.flt(:,3:end)];
    for ii = 1:flt3.num
        Nd      = flt3.flt(ii,17); 
	Ns      = flt3.flt(ii,18); 
	fltnum  = Nd*Ns;
        fltname = flt3.name{ii};
        % find subfaults for the master fault
        subInd = strcmpi(fltname,subflt.name);
        % find dips for the master fault
        dipInd = strcmpi(fltname,addon.dipname);
        [ newflt,~,xyzflt ] = GTdef_prjflt3dif(tmpflt(ii,:),subflt.flt(subInd,:),addon.dip(dipInd,:));
        %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
        % newflt = [ dnum snum x1  y1  z1  z2  len str dip slips];
        % slips  = [ rake rs ts rake0 rakeX rs0 rsX ts0 tsX ];
        dnum  = newflt(:,1);  snum = newflt(:,2);
        x1    = newflt(:,3);  y1   = newflt(:,4);
        z1    = newflt(:,5);  z2   = newflt(:,6);
        len   = newflt(:,7);  str  = newflt(:,8);
        dp    = newflt(:,9);
        rake  = newflt(:,10); slip = newflt(:,11);
        width = (z2-z1)./abs(sind(dp));
	subfltnum = size(newflt,1);
	inum  = [ 1:subfltnum ]';
        % different output formats
        if strfind(outFlag,'GTdef')
	   if strfind(outFlag,'center') 
              xx    =  xyzflt.xyzctr(:,1);  
              yy    =  xyzflt.xyzctr(:,2); 
              depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           end
           if strfind(outFlag,'topleft')
              xx    =  xyzflt.xyztop1(:,1);  
              yy    =  xyzflt.xyztop1(:,2); 
              depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           end
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ inum dnum snum lon lat depth len width str dp rake slip ];
        elseif strcmp(outFlag,'relax')
           east  =  xyzflt.xyztop1(:,1);  
           north =  xyzflt.xyztop1(:,2); 
           depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           out = [ inum slip north east depth len width str dp rake ];
        elseif strcmp(outFlag,'inverse2')
           xx    =  xyzflt.xyzctr(:,1);  
           yy    =  xyzflt.xyzctr(:,2); 
           depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ lon lat depth len width rake str dp slip ];
        else
           error('GTdef_GTdef2subfaults ERROR: output format is wrong!!!');
        end

        % output to a file
        [ fout ] = GTdef_GTdef2subfaults_head(outFlag,fltname);
        GTdef_GTdef2subfaults_body(outFlag,fout,out)
    end
toc
end

newflt = []; xyztop = []; xyzctr = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fault4 data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flt4.num~=0
fprintf(1,'\n.......... processing fault type-4 ..........\t');
tic
    % convert to cartesian coordinate
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
    tmpflt = [ x4_1 y4_1 x4_2 y4_2 flt4.flt(:,5:end) ];
    for ii = 1:flt4.num
        Nd      = flt4.flt(ii,17); 
        Ns      = flt4.flt(ii,18); 
        fltnum  = Nd*Ns;
        fltname = flt4.name{ii};
        % find the subfaults for the master fault
        subInd = strcmpi(fltname,subflt.name);
        % find dips for the master fault
        dipInd = strcmpi(fltname,addon.dipname);
        % find strikes for the master fault
        strInd = strcmpi(fltname,addon.strname);

        [ newflt,~,xyzflt ] = GTdef_prjflt4dif(tmpflt(ii,:),subflt.flt(subInd,:),addon.dip(dipInd,:),addon.crt(strInd,:));
        %            (1)  (2)  (3) (4) (5) (6) (7) (8) (9) (10-18)
        % newflt = [ dnum snum x1  y1  x2  y2  z1  z2  dip slips ];
        % slips  = [ rake rs ts rake0 rakeX rs0 rsX ts0 tsX ];
        dnum  = newflt(:,1);  snum = newflt(:,2);
        x1    = newflt(:,3);  y1   = newflt(:,4);
        x2    = newflt(:,5);  y2   = newflt(:,6);
        z1    = newflt(:,7);  z2   = newflt(:,8);
        dp    = newflt(:,9);
        rake  = newflt(:,10); slip = newflt(:,11);
        str   = GTdef_strike(x1,y1,x2,y2);
        len   = sqrt((x2-x1).^2+(y2-y1).^2);
        width = (z2-z1)./abs(sind(dp));
	subfltnum = size(newflt,1);
	inum  = [ 1:subfltnum ]';
        % different output formats
        if strfind(outFlag,'GTdef')
	   if strfind(outFlag,'center') 
              xx    =  xyzflt.xyzctr(:,1);  
              yy    =  xyzflt.xyzctr(:,2); 
              depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           end
           if strfind(outFlag,'topleft')
              xx    =  xyzflt.xyztop1(:,1);  
              yy    =  xyzflt.xyztop1(:,2); 
              depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           end
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ inum dnum snum lon lat depth len width str dp rake slip ];
        elseif strcmp(outFlag,'relax')
           east  =  xyzflt.xyztop1(:,1);  
           north =  xyzflt.xyztop1(:,2); 
           depth = -xyzflt.xyztop1(:,3); % from positive up to positive down (depth)
           out = [ inum slip north east depth len width str dp rake ];
        elseif strcmp(outFlag,'inverse2')
           xx    =  xyzflt.xyzctr(:,1);  
           yy    =  xyzflt.xyzctr(:,2); 
           depth = -xyzflt.xyzctr(:,3); % from positive up to positive down (depth)
           % convert to geographic coordinate if necessary
           switch coord
               case 'geo'
               [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
               case 'geo_polyconic'
               [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
           end
           out = [ lon lat depth len width rake str dp slip ];
        else
           error('GTdef_GTdef2subfaults ERROR: output format is wrong!!!');
        end

        % output to a file
        [ fout ] = GTdef_GTdef2subfaults_head(outFlag,fltname);
        GTdef_GTdef2subfaults_body(outFlag,fout,out)
    end
toc
end
