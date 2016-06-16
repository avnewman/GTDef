function [ modspace ] = GTdef_invert(modspace,pnt,los,bsl,beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              GTdef_invert                               %
% Invert the slip values using                                            %
% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)                                    %
% Here we have no inequalities, so set A=[];b=[]                          %
%                                                                         %
% INPUT:                                                                  %
% modspace - model structure                                              %
% pnt      - point structure                                              %
% los      - los structure                                                %
% bsl      - baseline structure                                           %
% beta     - current beta value                                           %
%                                                                         %
% OUTPUT:                                                                 %
% add .xx to modspace                                                     %
% modspace.xx - final values for ss,ds,ts from inversion     [slip_num*1] %
%                                                                         %
% first created by Lujia Feng Tue May  5 20:21:47 EDT 2009		  %
% used beta = beta^2 instead of beta lfeng Wed Dec  2 23:35:12 EST 2009   %
% used structure lfeng Wed Feb 22 19:36:43 SGT 2012			  %
% added modspace structure lfeng Thu Mar 19 17:32:29 SGT 2015             %
% added condensing sm & sm_abs lfeng Fri Mar 20 20:28:15 SGT 2015         %
% added InSAR los & Lgrn lfeng Tue Nov  3 23:17:37 SGT 2015               %
% added Model and Data Resolution Matrices, output of weighted Green's    %
% function matrix (and other data into modspace) anewman May 10 EDT 2016  %
% added Aineq & bineq to modspace lfeng Mon Jun  6 15:46:35 SGT 2016      %
% added condensing Aineq & Aeq lfeng Thu Jun 16 16:48:37 SGT 2016         %
% last modified by Lujia Feng Thu Jun 16 16:49:07 SGT 2016                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xgrn  = modspace.Xgrn;  
Lgrn  = modspace.Lgrn;  
Bgrn  = modspace.Bgrn;
lb    = modspace.lb;
ub    = modspace.ub;
x0    = modspace.x0;

% condense inequality equations
[ Aineq,ind ] = GTdef_condense(modspace.Aineq);
bineq = modspace.bineq(ind);
modspace.Aineq = Aineq;
modspace.bineq = bineq;

% condense equality equations
[ Aeq,ind ] = GTdef_condense(modspace.Aeq);
beq = modspace.beq(ind);
modspace.Aeq = Aeq;
modspace.beq = beq;

% condense the smoothing matrix by removing rows of all zeros
sm     = GTdef_condense(modspace.sm);
sm_abs = GTdef_condense(modspace.sm_abs);
modspace.sm     = sm;
modspace.sm_abs = sm_abs;

C = []; d = [];
% add error & weight to the equations
if ~isempty(Xgrn)
    eqNum = size(Xgrn,1); % equation number
    for ii = 1:eqNum
        Xgrn(ii,:)  = Xgrn(ii,:).*pnt.coef(ii);
        pnt.obs(ii) = pnt.obs(ii).*pnt.coef(ii);
    end
    ind = find(~isnan(pnt.obs));		% exclude nan values
    C = [ C;Xgrn(ind,:) ]; d = [ d;pnt.obs(ind) ];
end

if ~isempty(Lgrn)
    eqNum = size(Lgrn,1);
    for ii = 1:eqNum
        Lgrn(ii,:)  = Lgrn(ii,:).*los.coef(ii);
        los.obs(ii) = los.obs(ii).*los.coef(ii);
    end
    ind = find(~isnan(los.obs));		% exclude nan values
    C = [ C;Lgrn(ind,:) ]; d = [ d;los.obs(ind) ];
end

if ~isempty(Bgrn)
    eqNum = size(Bgrn,1);
    for ii = 1:eqNum
        Bgrn(ii,:)    = Bgrn(ii,:).*bsl.coef(ii,1);
	bsl.obs(ii,1) = bsl.obs(ii,1).*bsl.coef(ii,1);
    end
    ind = find(~isnan(bsl.obs));		% exclude nan values
    C = [ C;Bgrn(ind,:) ]; d = [ d;bsl.obs(ind) ];
end

nrealdata=size(C,1); % identify the size of the real-data part of the data vector (not zeros for the damped part.).
% add beta to the smoothing matrix
sm = sm.*beta;

% combine the green matrices and the smoothing matrix
d_sm = zeros(size(sm,1),1);
C = [ C;sm ]; d = [ d;d_sm ];

options = optimset('MaxIter',2000,'TolFun',1e-30);
[xx,resnorm] = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0,options);
%[xx,resnorm,res,exitflag,output] = lsqlin(C,d,Aineq,bineq,Aeq,beq,lb,ub,x0);
fprintf(1,'resnorm = %-12.5e\n',resnorm);

modspace.xx = xx;

%% New for model resolution matrix

ix=find(isfinite(lb));
id=[[1:nrealdata]';ix+nrealdata];
%C0=C(id,ix);  % operator matrix removing non-inverted data (includes smoothing)
C0=C(id,:);    % operator matrix (includes smoothing)
d=d(1:nrealdata);
x0=x0(ix);

% Develop weighted and damped over-determined generalized inverse matrix, Gg
% for the model resolution matrix  (See Menke, 2012 -- Chapter 4)
G=C0(1:nrealdata,:);     % strip off the smoothing part of the matrix
%e2I=sm(ix,ix).^2; (version that ignored non-inverted data)      
e2I=sm.^2;               % weighted smoothing matrix
Gg=inv(G'*G+e2I)*G';     % overdetermined.
%Gg = G'*inv(G*G'+e2I);  % same thing, but for underdetermined 
R=Gg*G ;            % Model Resolution Matrix 
N=G*Gg ;            % Data Resolution Matrix


modspace.resnorm = resnorm;
modspace.C       = C;
modspace.d       = d;
modspace.Gg      = Gg;
modspace.G       = G;
modspace.e2I     = e2I;
modspace.N       = N;
modspace.R       = R;
