function [] = GTdef_output(filename,earth,modspace,beta,...
           		   flt1,flt2,flt3,flt4,flt5,flt6,flt7,...
                           subflt,addon,...
          		   pnt,los,bsl,prf,grd,nod)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   GTdef_output.m                                       %
% output GTdef model results				                                 %
%										         %
% INPUT:									         %
% beta     - current beta value                                                          %
% earth    - earth structure                                                             %
% modspace - model structure                                                             %
% flt?     - fault structure                                                             %
% subflt   - subfault structure                                                          %
% addon    - addon structure                                                             %
% pnt      - point structure	  	                                                 %
% los      - los structure                                                               %
% bsl      - baseline structure                                                          %
% prf      - profile structure	  	                                                 %
% grd      - grid structure                                                              %
% nod      - node structure                                                              %
%										         %
% OUTPUT: an output file                                                                 %
%                                                                                        %
% first created by Lujia Feng Wed May  6 20:58:39 EDT 2009			         %
% added beta lfeng Wed Dec  2 23:42:51 EST 2009					         %
% added 1st derivative r_1d lfeng Thu Dec  3 01:27:48 EST 2009			         %
% added flag 'dip' by lfeng Mon Dec  7 01:05:28 EST 2009		  	         %
% added fault5 by lfeng Fri Dec 11 13:00:11 EST 2009				         %
% corrected the wrong 'freesurface' flag and changed it to 'surface' flag	         %
%    lfeng Wed Feb 24 13:26:18 EST 2010					                 %
% added the units in the output	lfeng Wed Jul 21 17:06:07 EDT 2010		         %
% use cell array of strings for names lfeng Wed Dec  1 17:36:42 EST 2010	         %
% use 4 digits after . for slip constraints lfeng Thu Dec  9 03:52:06 EST 2010           %
% test existence before output lfeng Thu Apr 14 12:58:16 EDT 2011		         %
% used structure & added layered output lfeng Wed Feb 22 14:31:31 SGT 2012	         %
% merged fault1 & fault3 and fault2 & fault4 lfeng Thu May 10 17:03:16 SGT 2012          %
% modified fault5 lfeng Sat Dec  1 00:32:44 SGT 2012                                     %
% added addon to combine dip & strike lfeng Fri Oct 24 14:56:14 SGT 2014                 %
% added earth structure lfeng Fri Mar 20 11:39:54 SGT 2015                               %
% output stress drop lfeng Wed Mar 25 17:37:14 SGT 2015                                  %
% added sdropflag lfeng Thu Mar 26 17:35:06 SGT 2015                                     %
% added fault5 external geometry lfeng Tue Jun 23 13:33:47 SGT 2015                      %
% added output origin lfeng with Paul M.  Wed Jun 24 12:34:52 SGT 2015                   %
% added InSAR los & Lgrn lfeng Tue Nov  3 23:17:37 SGT 2015                              %
% added fault6, changed old fault6 to fault7 lfeng Thu Jun  2 10:38:24 SGT 2016          %
% last modified by Lujia Feng Thu Jun  2 10:47:44 SGT 2016                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

etype     = earth.type;
rigidity  = earth.rigidity;
poisson   = earth.poisson;
edgrn     = earth.edgrn;
layer     = earth.layer;

coord     = modspace.coord;
origin    = modspace.origin;
smooth    = modspace.smooth;
surf      = modspace.surf;
modinfo   = modspace.modinfo(end,:);
sdropflag = modspace.sdropflag;

fout = fopen(filename,'w');

%%%%%%%%%% coordinates %%%%%%%%%%
if ~isempty(coord)
    fprintf(fout,'coord   \t%s\n',coord);
end
if ~isempty(origin)
    fprintf(fout,'origin   %f   %f\n',origin);
end

%%%%%%%%%% model info %%%%%%%%%%
if ~isempty(modinfo)
    fprintf(fout,'#data_num\t%d\n#slip_num\t%d\n#ndf     \t%d\n#rss      \t%-12.5e [m^2]\n#rms      \t%-12.5e [m]\n#wrss   \t%-12.5e [m^2]\n#wrms   \t%-12.5e [m]\n#chi2    \t%-12.5e\n#rchi2   \t%-12.5e\n#r_1d   \t%-12.5e [cm/km]\n#r_2d   \t%-12.5e [cm/km^2]\n#strain   \t%-12.5e [cm/km]\n',modinfo);
    fprintf(fout,'\n');
end

%%%%%%%%%% earth model %%%%%%%%%%
if strcmpi(etype,'homogeneous')||strcmpi(etype,'homo')
    fprintf(fout,'earth\thomogeneous\t%-10.2e\t%-6.4f\n',rigidity,poisson);
    fprintf(fout,'\n');
elseif strcmpi(etype,'layered')
    fprintf(fout,'earth\tlayered   %.0f %-8.4e  %.0f %-8.4e %-8.4e  %.0f %-8.4e %-8.4e   %.4f\n',...
            edgrn.nl,edgrn.obsz,edgrn.nr,edgrn.minr,edgrn.maxr,edgrn.nz,edgrn.minz,edgrn.maxz,edgrn.srate);
    fprintf(fout,'\n');
end

%%%%%%%%%% layers %%%%%%%%%%
if ~isempty(layer)
    fprintf(fout,'layer   %.0f   %8.4e  %-8.4e  %-8.4e  %-8.4e\n',layer');
    fprintf(fout,'\n');
end

%%%%%%%%%% smooth parameters %%%%%%%%%%
if ~isempty(beta)
    kappa = sqrt(beta);
    fprintf(fout,'kappa   \t 1  \t%-12.2f\nbeta     \t  1   \t%-12.2f\n\n',kappa,beta);
end

if ~isempty(smooth)
    fprintf(fout,'smooth   \t%s\n',smooth);
end
if ~isempty(surf)
    fprintf(fout,'surface  \t%s\n',surf);
end
fprintf(fout,'\n');

%%%%%%%%%%  fault 1  %%%%%%%%%%
% note: output flt1.out
for ii =1:flt1.num
    fltname = flt1.name{ii};
    if strcmp(sdropflag,'on')
       fprintf(fout,'# stress drop %s %10.4f MPa\n',fltname,flt1.sdrop(ii)*1e-6);
    end
    fprintf(fout,'fault 1 %s  %-14.8f %-12.8f %-6.4e %12.4e %-12.4e %-5.2f %-5.2f  %10.5f %-8.5f %-8.5f  %-5.4f %-5.4f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',fltname,flt1.out(ii,:));
    ind = strcmpi(fltname,subflt.outname);
    num = sum(ind);
    subflt1 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.5f %12.5f %12.5f  %8.4f %8.4f  %8.4f %8.4f  %8.4f %8.4f\n',fltname,subflt1(jj,:));
    end
end

%%%%%%%%%%  fault 2  %%%%%%%%%%
% note: output flt2.out
for ii =1:flt2.num
    fltname = flt2.name{ii};
    if strcmp(sdropflag,'on')
       fprintf(fout,'# stress drop %s %10.4f MPa\n',fltname,flt2.sdrop(ii)*1e-6);
    end
    fprintf(fout,'fault 2 %s  %-14.8f %-12.8f %14.8f %-12.8f %12.4e %-12.4e %-5.2f  %10.5f %-8.5f %-8.5f  %-5.4f %-5.4f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',fltname,flt2.out(ii,:));
    ind = strcmpi(fltname,subflt.outname);
    num = sum(ind);
    subflt2 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.5f %12.5f %12.5f  %8.4f %8.4f  %8.4f %8.4f  %8.4f %8.4f\n',fltname,subflt2(jj,:));
    end
end

%%%%%%%%%%  fault 3  %%%%%%%%%%
for ii =1:flt3.num
    fltname = flt3.name{ii};
    if strcmp(sdropflag,'on')
       fprintf(fout,'# stress drop %s %10.4f MPa\n',fltname,flt3.sdrop(ii)*1e-6);
    end
    fprintf(fout,'fault 3 %s  %-14.8f %-12.8f %-6.4e %12.4e %-12.4e %-5.2f %-5.2f  %10.2f %-8.5f %-8.5f  %-7.2f %-7.2f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',fltname,flt3.out(ii,:));
    ind = strcmpi(fltname,subflt.outname);
    num = sum(ind);
    subflt3 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.2f %12.5f %12.5f  %7.2f %7.2f  %8.4f %8.4f  %8.4f %8.4f\n',fltname,subflt3(jj,:));
    end
end

%%%%%%%%%%  fault 4  %%%%%%%%%%
for ii =1:flt4.num
    fltname = flt4.name{ii};
    if strcmp(sdropflag,'on')
       fprintf(fout,'# stress drop %s %10.4f MPa\n',fltname,flt4.sdrop(ii)*1e-6);
    end
    fprintf(fout,'fault 4 %s  %-14.8f %-12.8f %14.8f %-12.8f %12.4e %-12.4e %-5.2f  %10.2f %-8.5f %-8.5f  %-7.2f %-7.2f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d %d\n',fltname,flt4.out(ii,:));
    ind = strcmpi(fltname,subflt.outname);
    num = sum(ind);
    subflt4 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.2f %12.5f %12.5f  %7.2f %7.2f  %8.4f %8.4f  %8.4f %8.4f\n',fltname,subflt4(jj,:));
    end
end

%%%%%%%%%%  fault 5  %%%%%%%%%%
for ii =1:flt5.num
    fltname = flt5.name{ii};
    geoname = flt5.geoname{ii};
    colname = flt5.colname{ii};
    if strcmp(sdropflag,'on')
       fprintf(fout,'# stress drop %s %10.4f MPa\n',fltname,flt5.sdrop(ii)*1e-6);
    end
    fprintf(fout,'fault 5   %s   %s   %s   %10.5f %-8.5f %-8.5f  %-5.4f %-5.4f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d  %d \n',fltname,geoname,colname,flt5.out(ii,:));
    ind = strcmpi(fltname,subflt.outname);
    num = sum(ind);
    subflt5 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.5f %12.5f %12.5f  %8.4f %8.4f  %8.4f %8.4f  %8.4f %8.4f\n',fltname,subflt5(jj,:));
    end
end

%%%%%%%%%%  fault 6  %%%%%%%%%%
for ii =1:flt6.num
    fltname = flt6.name{ii};
    geoname = flt6.geoname{ii};
    colname = flt6.colname{ii};
    if strcmp(sdropflag,'on')
       fprintf(fout,'# stress drop %s %10.4f MPa\n',fltname,flt6.sdrop(ii)*1e-6);
    end
    fprintf(fout,'fault 6   %s   %s   %s   %10.2f %-8.5f %-8.5f  %-7.2f %-7.2f  %-5.4f %-5.4f  %-5.4f %-5.4f    %d  %d \n',fltname,geoname,colname,flt6.out(ii,:));
    ind = strcmpi(fltname,subflt.outname);
    num = sum(ind);
    subflt6 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.5f %12.5f %12.5f  %7.2f %7.2f  %8.4f %8.4f  %8.4f %8.4f\n',fltname,subflt6(jj,:));
    end
end

%%%%%%%%%%  fault 7  %%%%%%%%%%
for ii =1:flt7.num
    fltname = flt7.name{ii};
    grnname = flt7.grname{ii};
    fprintf(fout,'fault 6   %s   %s   %10.5f %10.5f %10.5f  %5.4f %5.4f  %5.4f %5.4f  %5.4f %5.4f    %d  %d \n',fltname,grnname,flt7.out(ii,:));
    ind = strcmpi(fltname,subflt.outname);
    num = sum(ind);
    subflt7 = subflt.out(ind,:);
    for jj = 1:num
        fprintf(fout,'     subfault %s  %5d %5d  %12.5f %12.5f %12.5f  %8.4f %8.4f  %8.4f %8.4f  %8.4f %8.4f\n',fltname,subflt7(jj,:));
    end
end

%%%%%%%%%%  dip  %%%%%%%%%%
for ii =1:addon.dipnum
    fprintf(fout,'dip   %s  %12.4f  %14.4e  %14.4e  %4d\n',addon.dipname{ii},addon.dip(ii,:));
end

%%%%%%%%%%  strike  %%%%%%%%%%
for ii =1:addon.strnum
    fprintf(fout,'strike   %s  %14.8f %12.8f %14.8f %12.8f %4d %14.3f\n',addon.strname{ii},addon.str(ii,:));
end

%%%%%%%%%%  point  %%%%%%%%%%
if strcmp(coord,'local')
    for ii =1:pnt.num
        fprintf(fout,'point 3 %s\t%14.5e %14.5e %14.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n', pnt.name{ii},pnt.out(ii,:));
    end
else
    for ii =1:pnt.num
        fprintf(fout,'point 3 %s\t%14.8f %12.8f %12.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n', pnt.name{ii},pnt.out(ii,:));
    end
end

%%%%%%%%%%  los point  %%%%%%%%%%
if strcmp(coord,'local')
    for ii =1:los.num
        fprintf(fout,'los 1 %s\t%14.5e %14.5e %14.5e  %10.5f %10.5f   %10.5f %10.5f %10.5f  %-5.2f\n', los.name{ii},los.out(ii,:));
    end
else
    for ii =1:los.num
        fprintf(fout,'los 1 %s\t%14.8f %12.8f %12.5e  %10.5f %10.5f   %10.5f %10.5f %10.5f  %-5.2f\n', los.name{ii},los.out(ii,:));
    end
end

%%%%%%%%%%  baseline  %%%%%%%%%%
if strcmp(coord,'local')
    for ii =1:bsl.num
        fprintf(fout,'baseline 3 %s\t%14.5e %14.5e %14.5e  %14.8f %12.8f %6.4e  %10.5f %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f %10.5f  %5.2f\n', bsl.name{ii},bsl.out(ii,:));
    end
else
    for ii =1:bsl.num
        fprintf(fout,'baseline 3 %s\t%14.8f %12.8f %12.5e  %14.8f %12.8f %6.4e  %10.5f %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f %10.5f  %5.2f\n', bsl.name{ii},bsl.out(ii,:));

    end
end

%%%%%%%%%%  profile  %%%%%%%%%%
for ii =1:prf.num
    name = prf.name{ii};
    name_len = length(name);
    fprintf(fout,'#profile %s  %-14.8f %-12.8f  %14.8f %-12.8f    %d\n',name,prf.prf(ii,:));
    ind = strncmpi(name,nod.name,name_len);
    num = sum(ind);
    cnod_name = nod.name(ind);
    cnod_out  = nod.out(ind,:);
    if strcmp(coord,'local')
        for jj = 1:num
            fprintf(fout,'point 3 %s\t%14.5e %14.5e %14.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));
        end
    else
        for jj = 1:num
            fprintf(fout,'point 3 %s\t%14.8f %12.8f %12.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));
        end
    end
end

%%%%%%%%%%  grid  %%%%%%%%%%
for ii =1:grd.num
    name = grd.name{ii};
    name_len = length(name);
    fprintf(fout,'#grid %s %-5.2f %-5.2f    %12.8f %-12.8f  %10.4f %-8.4f    %d  %d\n',name,grd.grd(ii,:));
    ind = strncmpi(name,nod.name,name_len);
    num = sum(ind);
    cnod_name = nod.name(ind);
    cnod_out  = nod.out(ind,:);
    if strcmp(coord,'local')
        for jj = 1:num
        	fprintf(fout,'point 3 %s\t%14.5e %14.5e %14.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));
        end
    else
        for jj = 1:num
        	fprintf(fout,'point 3 %s\t%14.8f %12.8f %12.5e  %10.5f %10.5f %10.5f  %10.5f %10.5f %10.5f  %-5.2f\n',cnod_name{jj},cnod_out(jj,:));

        end
    end
end
fclose(fout);

fprintf(1,'GTdef_output output %s\n',filename);
