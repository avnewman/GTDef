function [ ssflt1 ] = GTdef_stressfault1(ssflt1,addon) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_stressfault1                             %
% Process type-1 fault and generate its subfaults                         %
% use the center points for stress calculation                            %
%                                                                         %
% INPUT:                                                                  %
% ssflt1 structure                                                        %
%   ssflt1.fltnum ssflt1.fltname ssflt1.flt                               %
% addon structure                                                         %
% addon.dipname addon.dipnum & addon.dip                                  %
% addon.dip - [ dip z1 z2 rows ]                                          %
% addon.strname addon.strnum & addon.str                                  %
% addon.str - [ lon1 lat1 lon2 lat2 columns sweepAngle ]                  %
%                                                                         % 
% OUTPUT:                                                                 %
% ssflt1.crt  ssflt1.str  ssflt1.dip  ssflt1.rake ssflt1.fric             %
% ssflt1.dnum ssflt1.snum                                                 %
%                                                                         %
% first created by Lujia Feng Mon Jun 11 12:52:51 SGT 2012                %
% changed dip_struct to addon lfeng Fri Oct 24 15:10:01 SGT 2014          %
% modified outputs of GTdef_prjflt1uni.m lfeng Wed Apr 27 SGT 2016        %
% last modified by Lujia Feng Wed Apr 27 23:20:49 SGT 2016                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization
ssflt1.num  = 0;  ssflt1.name = {};
ssflt1.loc  = []; ssflt1.crt  = []; 
ssflt1.str  = []; ssflt1.dip  = [];   
ssflt1.rake = []; ssflt1.fric = [];  
ssflt1.dnum = []; ssflt1.snum = [];

% prepare the parameters
for ii = 1:ssflt1.fltnum
    % ssflt1.flt = [x1 y1 z1 z2 len str dip rake fric Nd Ns]         
    flt = ssflt1.flt(ii,:);
    % Note: flt is a row vector for the master fault
    mx1  = flt(1); my1  = flt(2); 
    mz1  = flt(3); mz2  = flt(4); 
    mlen = flt(5); mstr = flt(6);  mdip = flt(7);
    mrak = flt(8); mfrc = flt(9);
    Nd = round(flt(10)); Ns = round(flt(11));
    subflt_num = Nd*Ns;
    mx2 = mx1+mlen*sind(mstr);     my2 = my1+mlen*cosd(mstr);	% endpoint 2

    % fault parameters
    unit = ones(subflt_num,1);
    dlen = mlen/Ns; len = dlen*unit; str = mstr*unit;
    slips = zeros(subflt_num,9);
    % master fault name
    cflt_name = ssflt1.fltname{ii};
    for kk = 1:subflt_num, cname{kk} = cflt_name; end

    % exclusively specify dips for rows between z1 and z2 depth range
    if addon.dipnum~=0
        % find dips for the master fault
        dipInd = strcmpi(cflt_name,addon.dipname);
    else
        dipInd = [];
    end
    if ~isempty(dipInd)
        [ x1,y1,x2,y2,z1,z2,dip,ddip ] = GTdef_diffdips(mx1,my1,mx2,my2,mz1,mz2,mstr,mdip,Nd,Ns,addon.dip(dipInd,:));
    else
        dip = mdip*unit; 			% duplicate dips
        zlin = linspace(mz1,mz2,Nd+1)';
        z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));
        z1 = reshape(z1mat,[],1); z2 = reshape(z2mat,[],1);
        xlin = linspace(mx1,mx2,Ns+1); ylin = linspace(my1,my2,Ns+1); 
        x1mat = xlin(ones(Nd,1),1:end-1);  y1mat = ylin(ones(Nd,1),1:end-1);
        x1 = reshape(x1mat,[],1);  y1 = reshape(y1mat,[],1);
        dz = (mz2-mz1)/Nd; ddip = dz/sind(mdip);
    end

    dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
    dmat = dlin(1:end,ones(1,Ns)); smat = slin(ones(Nd,1),1:end);
    dnum = reshape(dmat,[],1);  snum = reshape(smat,[],1);

    % calculate the center point
    % flt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
    newflt = [ dnum snum x1 y1 z1 z2 len str dip slips ];
    [ ~,xyzflt ] = GTdef_prjflt1uni(newflt);
    xyzctr = xyzflt.xyzctr;

    % master fault name
    ssflt1.name = [ ssflt1.name; cname ];   
    ssflt1.num  = ssflt1.num + subflt_num;
    ssflt1.crt  = [ ssflt1.crt xyzctr ]; 
    ssflt1.dnum = [ ssflt1.dnum; dnum ];
    ssflt1.snum = [ ssflt1.snum; snum ];
    ssflt1.str  = [ ssflt1.str; str ]; 
    ssflt1.dip  = [ ssflt1.dip; dip ];   
    ssflt1.rake = [ ssflt1.rake; unit*mrak ]; 
    ssflt1.fric = [ ssflt1.fric; unit*mfrc ];  
end
