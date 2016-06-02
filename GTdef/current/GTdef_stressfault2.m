function [ ssflt2 ] = GTdef_stressfault2(ssflt2,addon) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          GTdef_stressfault2				  %
% Process type-2 fault and generate its subfaults                         %
% use the center points for stress calculation                            %
%                                                                         %
% INPUT:                                                                  %
% ssflt2 structure                                                        %
%   ssflt2.fltnum ssflt2.fltname ssflt2.flt                               %
% addon structure                                                         %
% addon.dipname addon.dipnum & addon.dip                                  %
% addon.dip - [ dip z1 z2 rows ]                                          %
% addon.strname addon.strnum & addon.str                                  %
% addon.str - [ lon1 lat1 lon2 lat2 columns sweepAngle ]                  %
%                                                                         % 
% OUTPUT:                                                                 %
% ssflt2.crt  ssflt2.str  ssflt2.dip  ssflt2.rake ssflt2.fric             %
% ssflt2.dnum ssflt2.snum                                                 %
%                                                                         %
% first created by Lujia Feng Mon Jun 11 14:10:31 SGT 2012                %
% changed dip_struct to addon lfeng Fri Oct 24 15:10:01 SGT 2014          %
% modified outputs of GTdef_prjflt2uni.m lfeng Wed Apr 27 SGT 2016        %
% last modified by Lujia Feng Wed Apr 27 23:21:43 SGT 2016                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialization
ssflt2.num  = 0;  ssflt2.name = {};
ssflt2.loc  = []; ssflt2.crt  = []; 
ssflt2.str  = []; ssflt2.dip = [];   
ssflt2.rake = []; ssflt2.fric = [];  
ssflt2.dnum = []; ssflt2.snum = [];

% prepare the parameters
for ii = 1:ssflt2.fltnum
    % ssflt2.flt = [x1 y1 x2 y2 z1 z2 dip rake fric Nd Ns]       
    flt = ssflt2.flt(ii,:);
    % Note: flt is a row vector for the master fault
    mx1 = flt(1); my1 = flt(2); mx2 = flt(3); my2 = flt(4); 
    mz1 = flt(5); mz2 = flt(6); mdip = flt(7);
    mrak = flt(8); mfrc = flt(9);
    Nd = round(flt(10)); Ns = round(flt(11));
    subflt_num = Nd*Ns; 
    mstr = GTdef_strike(mx1,my1,mx2,my2);

    unit  = ones(subflt_num,1);
    slips = zeros(subflt_num,9);
    % master fault name
    cflt_name = ssflt2.fltname{ii};
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
        z1mat = zlin(1:end-1,ones(1,Ns)); z2mat = zlin(2:end,ones(1,Ns));	% depths
        xlin = linspace(mx1,mx2,Ns+1); ylin = linspace(my1,my2,Ns+1); 
        x1mat = xlin(ones(Nd,1),1:end-1); y1mat = ylin(ones(Nd,1),1:end-1);	% endpoint 1
        x2mat = xlin(ones(Nd,1),2:end);   y2mat = ylin(ones(Nd,1),2:end);	% endpoint 2
        x1 = reshape(x1mat,[],1);  y1 = reshape(y1mat,[],1);
        x2 = reshape(x2mat,[],1);  y2 = reshape(y2mat,[],1);
        z1 = reshape(z1mat,[],1);  z2 = reshape(z2mat,[],1);
    end

    dlin = round(linspace(1,Nd,Nd)'); slin = round(linspace(1,Ns,Ns));
    dmat = dlin(1:end,ones(1,Ns)); smat = slin(ones(Nd,1),1:end);
    dnum = reshape(dmat,[],1);  snum = reshape(smat,[],1);

    % calculate the center point
    % flt=[dnum snum xx yy z1 z2 len str dip ss ds ts ss0 ssX ds0 dsX ts0 tsX]
    newflt = [ dnum snum x1 y1 x2 y2 z1 z2 dip slips ];
    [ ~,xyzflt ] = GTdef_prjflt2uni(newflt);
    xyzctr = xyzflt.xyzctr;
    
    ssflt2.name = [ ssflt2.name; cname ];   
    ssflt2.num  = ssflt2.num + subflt_num;
    ssflt2.crt  = [ ssflt2.crt xyzctr ]; 
    ssflt2.dnum = [ ssflt2.dnum; dnum ];
    ssflt2.snum = [ ssflt2.snum; snum ];
    ssflt2.str  = [ ssflt2.str; mstr*unit ]; 
    ssflt2.dip  = [ ssflt2.dip; dip ];   
    ssflt2.rake = [ ssflt2.rake; unit*mrak ]; 
    ssflt2.fric = [ ssflt2.fric; unit*mfrc ];  
end
