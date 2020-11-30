function [] = GTdef_output_stress(filename,sspnt,ssflt1,ssflt2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_output_stress.m				 %
% 		   function to output stress calculation results                 %
% INPUT:									 %
% sspnt  - stress point structure                                                %
%    sspnt.num, sspnt.name, sspnt.loc                                            %
%    sspnt.crt = [ xx; yy; zz ]                                                  %
%    sspnt.str,sspnt.dip,sspnt.rake,sspnt.fric                                   %
%    sspnt.shear,sspnt.normal,sspnt.coulomb					 %
% ssflt1 - stress fault1 structure                                               %
%    ssflt1.num, ssflt1.name                                                     %
%    ssflt1.flt = [lon1 lat1 z1 z2 len str dip rake fric Nd Ns]                  %
%    ssflt1.str,ssflt1.dip,ssflt1.rake,ssflt1.fric                               %
%    ssflt1.shear,ssflt1.normal,ssflt1.coulomb					 %
% ssflt2 - stress fault2 structure                                               %
%    ssflt2.num, ssflt2.name                                                     %
%    ssflt2.flt = [lon1 lat1 lon2 lat2 z1 z2 dip rake fric Nd Ns]                %
%    ssflt2.str,ssflt2.dip,ssflt2.rake,ssflt2.fric                               %
%    ssflt2.shear,ssflt2.normal,ssflt2.coulomb					 %
%                                                                                %
% OUTPUT: an output file                                                         %
%                                                                                %
% first created by Lujia Feng Fri May 18 10:46:52 SGT 2012                       %
% output fault lfeng Mon Jun 11 11:46:05 SGT 2012                                %
% last modified by Lujia Feng Mon Jun 11 15:52:17 SGT 2012                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fout = fopen(filename,'w');

fprintf(fout,'#stress calculation\n');
fprintf(fout,'#(1)stress (2)point/fault (3)name (4)lon (5)lat (6)z (7)strike[deg] (8)dip[deg] (9)rake[deg] (10)friction (11)shear[Pa] (12)normal[Pa] (13)coulomb[Pa] (14)Nd (15)Ns\n');

%%%%%%%%%%  stress calculation at points  %%%%%%%%%%
% Nd=0 Ns=0 for points
% stress point
if sspnt.num~=0
    for ii =1:sspnt.num
        fprintf(fout,'stress point %s    %14.8f  %-12.8f %-6.4e  %6.2f %-6.2f %-6.2f %4.2f  %12.5e %12.5e %12.5e %-3d %-3d\n',...
                sspnt.name{ii},sspnt.loc(ii,:),sspnt.str(ii,:),sspnt.dip(ii,:),sspnt.rake(ii,:),sspnt.fric(ii,:),...
	            sspnt.shear(ii,:),sspnt.normal(ii,:),sspnt.coulomb(ii,:),0,0);
    end
end

%%%%%%%%%%  stress calculation for fault1 %%%%%%%%%%
if ssflt1.num~=0
    fprintf(fout,'\n');
    for ii =1:ssflt1.num
        fprintf(fout,'stress fault %s    %14.8f  %-12.8f %-6.4e  %6.2f %-6.2f %-6.2f %4.2f  %12.5e %12.5e %12.5e %-3d %-3d\n',...
                ssflt1.name{ii},ssflt1.loc(ii,:),ssflt1.str(ii,:),ssflt1.dip(ii,:),ssflt1.rake(ii,:),ssflt1.fric(ii,:),...
	            ssflt1.shear(ii,:),ssflt1.normal(ii,:),ssflt1.coulomb(ii,:),ssflt1.dnum(ii,:),ssflt1.snum(ii,:));
    end
end

%%%%%%%%%%  stress calculation for fault2 %%%%%%%%%%
if ssflt2.num~=0
    fprintf(fout,'\n');
    for ii =1:ssflt2.num
        fprintf(fout,'stress fault %s    %14.8f  %-12.8f %-6.4e  %6.2f %-6.2f %-6.2f %4.2f  %12.5e %12.5e %12.5e %-3d %-3d\n',...
                ssflt2.name{ii},ssflt2.loc(ii,:),ssflt2.str(ii,:),ssflt2.dip(ii,:),ssflt2.rake(ii,:),ssflt2.fric(ii,:),...
                ssflt2.shear(ii,:),ssflt2.normal(ii,:),ssflt2.coulomb(ii,:),ssflt2.dnum(ii,:),ssflt2.snum(ii,:));
    end
end

fclose(fout);
