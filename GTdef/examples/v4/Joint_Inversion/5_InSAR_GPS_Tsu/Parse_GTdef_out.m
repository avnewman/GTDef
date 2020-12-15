function [] = Parse_GTdef_out(filename)
%filename='GT1_kp5000.out'

fin = fopen(filename,'r');

% use CELL ARRAY OF STRINGS for names
flt1.num = 0;    flt1.name = {};    flt1.flt = [];    flt1.sdrop = [];
% data structure
% - to read in
pnt.num = 0;     pnt.name = {};     pnt.loc = [];   pnt.disp = [];  pnt.err = [];   pnt.wgt = [];
los.num = 0;     los.name = {};     los.loc = [];   loc.disp = [];  los.err = [];   los.wgt = [];  los.dir = [];
% - to be built later
pnt.crt = []; pnt.obs = []; pnt.obs_err = []; pnt.obs_wgt = []; pnt.coef = [];
los.crt = []; los.obs = []; los.obs_err = []; los.obs_wgt = []; los.coef = [];

ln = 0; % line number
faultinfo.coord=[];
    while(1)
    % read in one line
    tline = fgetl(fin);
    ln=ln+1;
    % test if it is the end of file; exit if yes
    if ischar(tline)~=1, break; end
    % read the 1st term that is a flag for data type; could be some white-spaces before the 1st term
    [flag,remain] = strtok(tline);              % the default delimiter is white-space
    % omit '#,%, and blank' comment lines
    if (GTdef_skip(flag)), continue; end
    if strcmpi(flag,'coord')
        faultinfo.coord=tline;
    
    elseif strcmpi(flag,'fault')
        [method,remain] = strtok(remain);
        %% fault 1 %%
        if strcmp(method,'1')
            [name,remain] = strtok(remain);
            flt1.num = flt1.num+1; flt1.name = [ flt1.name; name ];
            for ii = 1:18
                [ flt1.flt(flt1.num,ii),remain ] = GTdef_read1double(remain);
            end
        end
        
        
        
         %%%%% point %%%%%
    elseif strcmpi(flag,'point')
        [method,remain] = strtok(remain);
 %% method 3 %%
        if strcmp(method,'1')
            [name,remain] = strtok(remain);
        elseif strcmp(method,'2')
            [name,remain] = strtok(remain);

        elseif strcmp(method,'3')
            [name,remain] = strtok(remain);
            pnt.num = pnt.num+1; pnt.name = [ pnt.name; name ];
            %% point location [lon lat z] %%
            for ii = 1:3
                [ pnt.loc(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% displacements [east north vert] %%
            for ii = 1:3
                [ pnt.disp(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% errors [east north vert] %%
            for ii = 1:3
                [ pnt.err(pnt.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% weight %%
            [str,remain] = strtok(remain);
            if GTdef_skip(str)  % if weight is absent, use default 1
                pnt.wgt(pnt.num,1) = 1;
            else
                pnt.wgt(pnt.num,1) = str2double(str);
            end
        end
        
        
        %%%%% los displacement %%%%%
    elseif strcmpi(flag,'los')
        [method,remain] = strtok(remain);
        %% method 1 %%
        if strcmp(method,'1')
            [name,remain] = strtok(remain);
            los.num = los.num+1; los.name = [ los.name; name ];
            %% InSAR point location [lon lat z] %%
            for ii = 1:3
                [ los.loc(los.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% los displacement %%
            [ los.disp(los.num,1),remain ] = GTdef_read1double(remain);
            %% los displacement error %%
            [ los.err(los.num,1),remain ]  = GTdef_read1double(remain);
            %% unit vector for los direction [east north vert] %%
            for ii = 1:3
                [ los.dir(los.num,ii),remain ] = GTdef_read1double(remain);
            end
            %% weight %%
            [str,remain] = strtok(remain);
            if GTdef_skip(str) % if weight is absent, use default 1
                los.wgt(los.num,1) = 1;
            else
                los.wgt(los.num,1) = str2double(str);
            end
        else
        %% no other methods yet %%
            warning('Line %d, of %s: Input los type "%s" not recognised. Ignoring this data.',ln,filename,method)
        end
        
    end
    
    end
    
 

    %% move lon and lat to the corner of each subfautl
   for ii = 1:length(flt1.flt(:,1))
        d(ii)=flt1.flt(ii,3)/tand( flt1.flt(ii,7));
       [sflon(ii),sflat(ii)]=Haversine(flt1.flt(ii,1),flt1.flt(ii,2),flt1.flt(ii,6)+90,d(ii)/1000);
   end
   
   %% build subfault multisegment file for GMT
   Outfile=[filename, '_GMTsubflt'];
   convert_subflt_to_box( 'fault_param184.txt',flt1.flt(:,9),Outfile);
   
end