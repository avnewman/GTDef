function [dart,tsu,modspace,fault,flt1,addon,kappa,modinfo] = GTwave_open(filename)

fin = fopen(filename,'r');
%%%%%%%%%%%%%%%%%%%%%%%
kappa_num=0;

dart=[]; dart.id=[]; dart.loc=[]; dart.window=[]; dart.num=0;
tsu.obs=[];  tsu.wgt=[];  tsu.err=[]; tsu.obs_wgt=[]; tsu.obs_err=[]; tsu.coef=[];

modspace.C=[]; modspace.d=[];
modspace.Tgrn=[]; modspace.R=[];  modspace.R_ss=[]; modspace.R_ds=[]; modspace.R_ts=[];
modspace.ri=[];

modinfo=[];
fault=[]; 
flt1=[]; flt1.flt=[]; flt1.name=[]; flt1.num=0;

addon=[];
addon.dip=[]; addon.dipname=[]; addon.dipnum=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ln = 0; % line number
while(1)   
     % read in one line
     tline = fgetl(fin);
     ln=ln+1;
     % test if it is the end of file; exit if yes
     if ischar(tline)~=1, break; end
     
     % read the 1st term that is a flag for data type; could be some white-spaces before the 1st term
     [flag,remain] = strtok(tline);             % the default delimiter is white-space
     if (strncmpi(flag,'#',1)||strncmpi(flag,'%',1)||isempty(flag)), continue; end
     
     if strcmpi(flag,'kappa')
         [method,remain] = strtok(remain);
         %% method 1 %%
         if strcmp(method,'1')
             while true
                 [kk,remain] = strtok(remain);
                 %if GTdef_skip(kk), break; end
                 
                 if (strncmpi(kk,'#',1)||strncmpi(kk,'%',1)||isempty(kk)),break;end  
                 
                 kappa_num = kappa_num+1;
                 kappa(kappa_num) = str2double(kk);             
             end
         else
          fprintf('%s\n','kappa method 2 not yet supported');
         end
     
     %% FAULT DATA
         elseif strcmpi(flag,'fault')
         [method,remain] = strtok(remain);
        %% fault 1 %%
         if strcmp(method,'1')
             [name,remain] = strtok(remain);
             flt1.num = flt1.num+1; flt1.name = [ flt1.name; name ];
             for ii = 1:18
              [tmp,remain] = strtok(remain);
              flt1.flt(flt1.num,ii) = str2double(tmp);
             end
         else
             fprintf('%s\n','Other Fault types not yet supported. Use Fault 1');

         end
     %% DIP DATA
     elseif strcmpi(flag,'dip')
     [name,remain] = strtok(remain);
        addon.dipnum  = addon.dipnum+1; addon.dipname = [ addon.dipname; name ];
        % (1)dip (2)z1 (3)z2 (4)rows
         for ii = 1:4
             [tmp,remain] = strtok(remain);
              addon.dip(addon.dipnum,ii) = str2double(tmp);
             
         end
     
     %% DART DATA
       elseif strcmpi(flag,'DART')
         [name,remain] = strtok(remain);
         dart.id=[dart.id ; str2double(name)];
        
         
         dart.num=dart.num+1;
         for ii=1:2
         % dart location
            [tmp,remain] = strtok(remain);
            dart.loc(dart.num,ii) = str2double(tmp);
         end
         for ii=1:2
         % dart windowing
            [tmp,remain] = strtok(remain);
            dart.window(dart.num,ii) = str2double(tmp);
         end
         % error
            [tmp,remain] = strtok(remain);
            tsu.err(dart.num) = str2double(tmp);
         %weight
            [tmp,remain] = strtok(remain);
            tsu.wgt(dart.num) = str2double(tmp);
     end
end


fault.lon=flt1.flt(1); fault.Nd=flt1.flt(17);
fault.lat=flt1.flt(2); fault.Ns=flt1.flt(18);
fault.len=flt1.flt(5);
end

