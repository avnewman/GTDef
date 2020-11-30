function [ newflt1,newflt3,Nd,Ns ] = GTdef_read_geometry(finName,colName,origin,coord,outFlag,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             GTdef_read_geometry                               %
% convert any fault geometry or model to GTdef fault type 1 and/or fault type 3 %
%                                                                               %
% INPUT:                                                                        %
% (1) finName - external fault geometry file                                    %
% (2) colName  - name of each column                                            %
%                                                                               %
%      l ------------ r  ----> strike                                           %
%        |          |    |                                                      %
%        |    c     |    |                                                      %
%        |          |    v                                                      %
%        ------------    dip                                                    %
%                                                                               %
%  depth must be positive                                                       %
%  specify unit for dep & slip by attaching 'm' 'km' 'cm' 'mm'                  %
%  use '=' to set up a constant value for len [m], wid [m], str, & dip          %
%  Columns with unrecognized flags won't be read in                             %
% dnum,snum  - patch number along dip and strike                                %
% clon,clat  - geographic coordinates for center points                         % 
% llon,llat  - geographic coordinates for topleft points                        %  
% depm,depkm - depth in m or km                                                 % 
% str        - strike in degree                                                 %
% dip        - dip in degree                                                    %
% len        - length in meter                                                  %
% wid        - width in meter                                                   %
% dsm,dscm,smm  - dip slip in m, cm, or mm                                      %
% ssm,sscm,ssmm - strike slip in m, cm, or mm                                   %
% tsm,tscm,tsmm - tensile slip in m, cm, or mm                                  %
% rake       - rake in degree                                                   %
% slipm,slipcm,slipmm - rake slip in m, cm, or mm                               %
%                                                                               %
% Examples:                                                                     %
% USGS:                                                                         %
% colName = 'len=12000,wid=12000,clat,clon,depkm,slipcm,rake,str,dip'           %
% inverse2:                                                                     %
% colName = 'clon,clat,depm,str,dip,len,wid'                                    %
%                                                                               %
% (2) origin = [lon0 lat0]                                                      %
%            = [ 99.9695 -4.2290 ]                                              %
% (3) coord  = 'geo', 'geo_polyconic', or 'local'                               %
%                                                                               %
% (4) outFlag - output format                                                   %
%     if outFlag = [], do not output any files                                  %
% ----------------------------------------------------------------------------- %
%     a 'GTdef_center' or 'GTdef_topleft' uses geographic coordinates           %
%     1  2    3    4   5   6 7      8     9      10  11   12                    %
%     no dnum snum lon lat z length width strike dip rake slip                  %
%     'GTdef_center' uses center point                                          %
%     'GTdef_topleft' uses top left point                                       %
%                                                                               %
%     b 'relax' uses local coordinate                                           %
%     1  2    3     4    5     6      7     8      9   10                       %
%     no slip north east depth length width strike dip rake                     %
%     'relax' always uses top left point                                        %
% ----------------------------------------------------------------------------- %
%                                                                               %
% (5,6) Nd,Ns - optional                                                        %
% If it is too difficult for the code to figure out Nd & Ns automatically,      %
% need to mannually specify Nd & Ns                                             %
%                                                                               %
% OUTPUT:                                                                       %
% (1) geometry array                                                            %
% patch order: along dip first, then along strike                               %
%                                                 slips                         %
%                                    _______________|_______________            %
%                                    |                             |            %
% newflt1 = [ x1 y1 z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX ]        %
% newflt3 = [ x1 y1 z1 z2 len str dip rake rs ts rake0 rakeX rs0 rsX ts0 tsX ]  %
%                                    |_____________________________________|    %
%                                                   |                           %
%                                                 slips                         %
%                                                                               %
% Note: x1 & y1 for newflt1 need to be projected back to surface!!!!            %
%   but x1 & y1 for out should not be projected!!!!!                            %
%                                                                               %
% (2) save to subfault file in GTdef or relax format                            %
%                                                                               %
% first created by Lujia Feng & Paul Morgan Wed Jun 17 14:36:12 SGT 2015        %
% added units to dep & slip values lfeng & Paul Fri Jun 19 15:21:25 SGT 2015    %
% added '=' for constant varialbes lfeng Fri Jun 19 17:36:29 SGT 2015           %
% added output fault geometry lfeng Tue Jun 23 11:58:17 SGT 2015                %
% corrected x1 & y1 error lfeng with pmorgan lfeng Wed Jun 24 14:53:13 SGT 2015 %
% corrected negative depth error with pmorgan lfeng Thu Jul  2 11:06:52 SGT 2015%
% added dnum & snum to read in GTdef own geometry output lfeng Aug  5 SGT 2015  %
% added converting to fault type 3 lfeng Wed Jun  1 12:12:43 SGT 2016           %
% last modified by Lujia Feng Wed Jun  1 12:28:44 SGT 2016                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(varargin)
   Nd = varargin{1};
   Ns = varargin{2};
end 

defltCell = {'dnum' 'snum' 'clon' 'clat' 'llon' 'llat' 'depm' 'depkm' 'str' ...
             'dip' 'len' 'wid' ...
             'dsm' 'dscm' 'dsmm' 'ssm' 'sscm' 'ssmm'  'tsm' 'tscm' 'tsmm' ...
             'rake' 'slipm' 'slipcm' 'slipmm'};

if ~exist(finName,'file'), error('GTdef_read_geometry ERROR: %s does not exist!',finName); end
[ ~,colsCell ] = regexp(colName,',','match','split');
colnum = length(colsCell);

% check constant values '='
indCell = strfind(colsCell,'=');
newcolsCell = {};
for ii=1:colnum
   ind = indCell{ii};
   % '=' used
   if ~isempty(ind)
      [ ~,constCell ] = regexp(colsCell{ii},'=','match','split');
      par   = constCell{1}; 
      value = str2num(constCell{2});
      ind  = find(strcmp(par,{'str' 'dip' 'len' 'wid'}));
      fprintf(1,'\n');
      switch ind
         case 1
            str = value;
            fprintf(1,'GTdef_read_geometry WARNING: constant value %f is used for strike\n',str);
         case 2
            dip = value;
            fprintf(1,'GTdef_read_geometry WARNING: constant value %f is used for dip\n',dip);
         case 3
            len = value;
            fprintf(1,'GTdef_read_geometry WARNING: constant value %f is used for length\n',len);
         case 4
            wid = value;
            fprintf(1,'GTdef_read_geometry WARNING: constant value %f is used for width\n',wid);
	 otherwise
            fprintf(1,'GTdef_read_geometry WARNING: constant value can be used only for str, dip len, and wid\n');
      end
   % '=' not used
   else
      newcolsCell = [ newcolsCell colsCell{ii} ];
   end
end

% select only the columns we need
colsCell = newcolsCell;
colnum   = length(colsCell);
newcolsCell = {};
colPattern = '';
for ii=1:colnum
   if any(strcmp(colsCell{ii},defltCell))
      colPattern = [ colPattern ' %f' ];
      newcolsCell = [ newcolsCell colsCell{ii} ];
   else
      colPattern = [ colPattern ' %*s' ];
   end
end

% read in geometry
fin = fopen(finName,'r');
geoCell = textscan(fin,colPattern,'CommentStyle','#');
fclose(fin);
geoData = cell2mat(geoCell);
fltnum  = size(geoData,1);

% assign columns to correct variable
% read in lat & lon
ind  = strcmp('llon',newcolsCell);
llon = geoData(:,ind);
ind  = strcmp('llat',newcolsCell);
llat = geoData(:,ind);
ind  = strcmp('clon',newcolsCell);
clon = geoData(:,ind);
ind  = strcmp('clat',newcolsCell);
clat = geoData(:,ind);

% read in dnum & snum
ind = strcmp('dnum',newcolsCell);
if any(ind)
   dnum = geoData(:,ind);
end
ind = strcmp('snum',newcolsCell);
if any(ind)
   snum = geoData(:,ind);
end


% read in dep
ind = strcmp('depm',newcolsCell);
if any(ind)
   dep = geoData(:,ind);
else
   ind = strcmp('depkm',newcolsCell);
   if any(ind)
      dep = geoData(:,ind)*1e3;   % convert km to m
   else
      error('GTdef_read_geometry ERROR: need to specify depm or depkm!');
   end
end

% read in str
ind = strcmp('str',newcolsCell);
if any(ind)
   str = geoData(:,ind);
else
   if isscalar(str)
      str = str.*ones(fltnum,1);
   else
      error('GTdef_read_geometry ERROR: need to specify str!');
   end
end

% read in dip
ind = strcmp('dip',newcolsCell);
if any(ind)
   dip = geoData(:,ind);
else
   if isscalar(dip)
      dip = dip.*ones(fltnum,1);
   else
      error('GTdef_read_geometry ERROR: need to specify dip!');
   end
end

% read in len
ind = strcmp('len',newcolsCell);
if any(ind)
   len = geoData(:,ind);
else
   if isscalar(len)
      len = len.*ones(fltnum,1);
   else
      error('GTdef_read_geometry ERROR: need to specify len!');
   end
end

% read in wid
ind = strcmp('wid',newcolsCell);
if any(ind)
   wid = geoData(:,ind);
else
   if isscalar(wid)
      wid = wid.*ones(fltnum,1);
   else
      error('GTdef_read_geometry ERROR: need to specify wid!');
   end
end

% read in ss
ind  = strcmp('ssm',newcolsCell);
if any(ind)
   ss = geoData(:,ind);
else
   ind  = strcmp('sscm',newcolsCell);
   if any(ind)
      ss = geoData(:,ind)*1e-2;
   else
      ind  = strcmp('ssmm',newcolsCell);
      if any(ind)
         ss = geoData(:,ind)*1e-3;
      else
         ss = zeros(fltnum,1);
      end
   end
end

% read in ds
ind  = strcmp('dsm',newcolsCell);
if any(ind)
   ds = geoData(:,ind);
else
   ind  = strcmp('dscm',newcolsCell);
   if any(ind)
      ds = geoData(:,ind)*1e-2;
   else
      ind  = strcmp('dsmm',newcolsCell);
      if any(ind)
         ds = geoData(:,ind)*1e-3;
      else
         ds = zeros(fltnum,1);
      end
   end
end

% read in ts
ind  = strcmp('tsm',newcolsCell);
if any(ind)
   ts = geoData(:,ind);
else
   ind  = strcmp('tscm',newcolsCell);
   if any(ind)
      ts = geoData(:,ind)*1e-2;
   else
      ind  = strcmp('tsmm',newcolsCell);
      if any(ind)
         ts = geoData(:,ind)*1e-3;
      else
         ts = zeros(fltnum,1);
      end
   end
end

% read in rake
ind  = strcmp('rake',newcolsCell);
rake = geoData(:,ind);

% read in slip
ind  = strcmp('slipm',newcolsCell);
if any(ind)
   slip = geoData(:,ind);
else
   ind  = strcmp('slipcm',newcolsCell);
   if any(ind)
      slip = geoData(:,ind)*1e-2;
   else
      ind  = strcmp('slipmm',newcolsCell);
      if any(ind)
         slip = geoData(:,ind)*1e-3;
      else
         slip = [];
      end
   end
end

% convert rake+slip to ds+ss
if ~isempty(rake) && ~isempty(slip)
   ss = slip.*cosd(rake);
   ds = slip.*sind(rake);
% convert ds+ss to rake+slip
elseif ~isempty(ds) && ~isempty(ss)
   slip = sqrt(ss.^2+ds.^2);
   rake = atan2(ds,ss).*180/pi;
end

% from geographic coord to cartesian
lon0 = origin(1);
lat0 = origin(2);
switch coord
   case 'geo'
      if ~isempty(llon) && ~isempty(llat)
         [x1,y1] = LL2ckmd(llon,llat,lon0,lat0,0);
      elseif ~isempty(clon) && ~isempty(clat)
         [cx,cy] = LL2ckmd(clon,clat,lon0,lat0,0);
      else
         error('GTdef_read_geometry ERROR: need to specify llon,llat or clon,clat!');
      end
   case 'geo_polyconic'
      if ~isempty(llon) && ~isempty(llat)
         [x1,y1] = latlon_to_xy(llon,llat,lon0,lat0,0);
      elseif ~isempty(clon) && ~isempty(clat)
         [cx,cy] = latlon_to_xy(clon,clat,lon0,lat0,0);
      else
         error('GTdef_read_geometry ERROR: need to specify llon,llat or clon,clat!');
      end
   case 'local'
      if ~isempty(llon) && ~isempty(llat)
         x1 = llon; y1 = llat;
      elseif ~isempty(clon) && ~isempty(clat)
         cx = clon; cy = clat;
      else
         error('GTdef_read_geometry ERROR: need to specify llon,llat or clon,clat!');
      end
end

% convert center point to topleft point
if ~isempty(clon) && ~isempty(clat)
   ll  = 0.5*len;
   ww  = 0.5.*wid.*cosd(dip);
   dx  = ll.*sind(str) + ww.*cosd(str);
   dy  = ll.*cosd(str) - ww.*sind(str); 
   x1  = cx - dx;
   y1  = cy - dy;
   cz  = dep;
   z1  = cz - 0.5.*wid.*sind(dip);
   z2  = z1 + wid.*sind(dip);
end

% convert topleft point to center point
if ~isempty(llon) && ~isempty(llat)
   ll  = 0.5.*len;
   ww  = 0.5.*wid.*cosd(dip);
   dx  = ll.*sind(str) + ww.*cosd(str);
   dy  = ll.*cosd(str) - ww.*sind(str); 
   cx  = x1 + dx;
   cy  = y1 + dy;
   z1  = dep;
   z2  = z1 + wid.*sind(dip);
   cz  = z1 + 0.5.*wid.*sind(dip);
end

% check if depth z1 is negative
ind1 = z1<0;
if any(ind1)
   ind2 = z1<-1;
   if any(ind2)
      error('GTdef_read_geometry ERROR: depth z1 can not be negative!');
   else
      z1(ind1) = 0;
   end
end

% project topleft point to surface
if dip==90
   surf_x1 = x1; 
   surf_y1 = y1;
else
   dd  =  z1./tand(dip);
   dx  =  dd.*cosd(str);
   dy  = -dd.*sind(str); 
   surf_x1 = x1 - dx;
   surf_y1 = y1 - dy;
end

slips0X = zeros(fltnum,6); % [ss0 ssX ds0 dsX ts0 tsX] or [rake0 rakeX rs0 rsX ts0 tsX]
newflt1 = [ surf_x1 surf_y1 z1 z2 len str dip ss    ds  ts slips0X ];
newflt3 = [ surf_x1 surf_y1 z1 z2 len str dip rake slip ts slips0X ];

%plot3(cx,cy,cz,'ko'); hold on;

% find out Nd & Ns using depth, assuming depths all the same!!!!!!
% if not the case, use your brain to figure out Nd & Ns
uniquedep = unique(dep);
if exist('Nd')==0 && exist('Ns')==0
   Nd = length(uniquedep);
   Ns = fltnum/Nd;
   [ depCounts,~ ] = histcounts(dep,[-1e3; uniquedep+10]);
   if any(depCounts-Ns)
      error('GTdef_read_geometry ERROR: the number of patches does not match!');
   end
end
fprintf(1,'The number of patches is %d: %d along-dip; %d along-strike\n',fltnum,Nd,Ns);

% output different formats
if isempty(outFlag)
   out = [];
else
   if strfind(outFlag,'GTdef')
      if strfind(outFlag,'center') 
         xx = cx;  
         yy = cy; 
         zz = cz; % depth is positive down
      end
      if strfind(outFlag,'topleft')
         xx = x1;
         yy = y1; 
         zz = z1; % depth is positive down 
      end
      % convert cartesian to geographic coordinate if necessary
      switch coord
          case 'geo'
             [ lon,lat ] = ckm2LLd(xx,yy,lon0,lat0,0);
          case 'geo_polyconic'
             [ lon,lat ] = xy_to_latlon(xx,yy,lon0,lat0,0);
      end
      out = [ lon lat zz len wid str dip rake slip ]; % not including inum dnum snum 
   elseif strcmp(outFlag,'relax')
      east  = x1;  
      north = y1; 
      depth = z1; % depth is positive down
      out = [ slip north east dep len wid str dip rake ]; % not including inum
   else
      error('GTdef_read_geometry ERROR: output format is wrong!!!');
   end
end

% if dnum & snum not given, then the program need to figure out dnum & snum
if exist('dnum')==0 && exist('snum')==0
   % check along-dip or along-strike
   % if along-dip first
   if dep(1:Nd)==uniquedep 
      % check right-hand or not                                   
      xx12 = x1(1+Nd) - x1(1);                                    
      yy12 = y1(1+Nd) - y1(1);                                    
      dist = sqrt(xx12.^2+yy12.^2);                               
      if yy12>=0                                                  
         azim = atand(xx12./yy12);                                
         if azim<0, azim = 360+azim; end                          
      else                                                        
         azim = 180+atand(xx12./yy12);
      end
      azimdif = abs(azim-str(1));
   
      % if right-handed, need no further action
      if azimdif < 90
         fprintf(1,'\nThe input fault is along-dip & right-handed\n');
         fprintf(1,'                                     \n');
         fprintf(1,'   l ------------ r     ----> strike \n');
         fprintf(1,'    | | |        |      |            \n');
         fprintf(1,'    | | | .....> |      |            \n');
         fprintf(1,'    | V V        |      v            \n');
         fprintf(1,'    |            |     dip           \n');
         fprintf(1,'    --------------                   \n');
         if azimdif>45, fprintf(1,'GTdef_read_geometry WARNING: the difference between azimuth and strike is >45\n'); end
      % left-handed, need to swap the columns
      else                                                        
         fprintf(1,'\nThe input fault is along-dip & left-handed\n');
         fprintf(1,'                                     \n');
         fprintf(1,'   l ------------ r     ----> strike \n');
         fprintf(1,'    |        | | |      |            \n');
         fprintf(1,'    | <..... | | |      |            \n');
         fprintf(1,'    |        V V |      v            \n');
         fprintf(1,'    |            |     dip           \n');
         fprintf(1,'    --------------                   \n');
         ind = 1:fltnum;                                          
         indmat = reshape(ind,Nd,Ns);                             
         indmat = indmat(:,end:-1:1);                             
         ind = reshape(indmat,[],1);                              
         newflt1 = newflt1(ind,:);                                
         newflt3 = newflt3(ind,:);                                
         if ~isempty(out)                                         
             out = out(ind,:);                                    
         end
      end                                                                   
   % if along-strike first
   else 
      % check right-hand or not
      xx12 = x1(2) - x1(1);
      yy12 = y1(2) - y1(1);                                       
      dist = sqrt(xx12.^2+yy12.^2);                               
      if yy12>=0                                                  
         azim = atand(xx12./yy12);                                
         if azim<0, azim = 360+azim; end                          
      else                                                        
         azim = 180+atand(xx12./yy12);                            
      end                                                         
      azimdif = abs(azim-str(1));
   
      % if right-handed, transpose to along-dip first
      if azimdif < 90
         fprintf(1,'\nThe input fault is along-strike & right-handed\n');
         fprintf(1,'                                     \n');
         fprintf(1,'   l ------------ r     ----> strike \n');
         fprintf(1,'    | -------->  |      |            \n');
         fprintf(1,'    | -------->  |      |            \n');
         fprintf(1,'    |     .      |      v            \n');
         fprintf(1,'    |     .      |     dip           \n');
         fprintf(1,'    --------------                   \n');
         if azimdif>45, fprintf(1,'GTdef_read_geometry WARNING: the difference between azimuth and strike is >45\n'); end
         ind     = 1:fltnum;
         indmat  = reshape(ind,Ns,Nd);
         indmat  = indmat';
         ind     = reshape(indmat,[],1);
         newflt1 = newflt1(ind,:);
         newflt3 = newflt3(ind,:);
         if ~isempty(out)
            out = out(ind,:);
         end
      % left-handed, need to swap the columns and transpose       
      else                                                        
         fprintf(1,'\nThe input fault is along-strike & left-handed\n');
         fprintf(1,'                                     \n');
         fprintf(1,'   l ------------ r     ----> strike \n');
         fprintf(1,'    | <--------- |      |            \n');
         fprintf(1,'    | <--------- |      |            \n');
         fprintf(1,'    |     .      |      v            \n');
         fprintf(1,'    |     .      |     dip           \n');
         fprintf(1,'    --------------                   \n');
         ind     = 1:fltnum;                                      
         indmat  = reshape(ind,Ns,Nd);                            
         indmat  = indmat';                                       
         indmat  = indmat(:,end:-1:1);                            
         ind     = reshape(indmat,[],1);                          
         newflt1 = newflt1(ind,:);                                
         newflt3 = newflt3(ind,:);                                
         if ~isempty(out)
            out = out(ind,:);
         end
      end
   end
   inum  = [ 1:fltnum ]';
   dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
   dmat = dlin(1:end,ones(1,Ns));    smat = slin(ones(Nd,1),1:end);
   dnum = reshape(dmat,[],1);        snum = reshape(smat,[],1);
% reorder patches according to dnum & snum
else
   numMat = [ dnum snum ];
   [ numMat,ind ] = sortrows(numMat,[2 1]);
   if ~isempty(out)
      out = out(ind,:);
   end
   inum = [ 1:fltnum ]';
   dnum = dnum(ind);
   snum = snum(ind);
end

%plot3(newflt1(:,1),newflt1(:,2),newflt1(:,3),'r');

%---------------------------------- output to a file ----------------------------------
if ~isempty(out)
   if strfind(outFlag,'GTdef')
      % GTdef out = [ inum dnum snum lon lat zz len wid str dip rake slip ];
      out = [ inum dnum snum out ];
   elseif strcmp(outFlag,'relax')
      % relax out = [ inum slip north east dep len wid str dip rake ];
      out = [ inum out ];
   end
   [ ~,basename,~ ] = fileparts(finName);
   [ fout ] = GTdef_GTdef2subfaults_head(outFlag,basename);
   GTdef_GTdef2subfaults_body(outFlag,fout,out);
end
