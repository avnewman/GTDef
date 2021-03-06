function [ modspace ] = GTdef_addall(modspace,Xgrn1,Lgrn1,Bgrn1,Ngrn1,sm1,sm_abs1,Aeq1,beq1,lb1,ub1,x01)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_addall				                 %
% Add new green1,Aeq1,beq1,lb1,ub1,x01 to pre-existing green0,Aeq0,beq0,                 %
% lb0,ub0,x00 for x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)			                 %
% Here we have no inequalities, so set A=[];b=[]			                 %
%									                 %
% INPUT:					  		  	                 %
% modspace structure                                                                     %
%                                                                                        %
% OUTPUT:					  		  	                 %
% modspace structure                                                                     %
%                                                                                        %
% first created by Lujia Feng Fri Apr 24 14:38:50 EDT 2009		                 %
% added auxiliary sm lfeng Thu Dec  3 00:58:41 EST 2009			                 %
% replaced auxiliary sm with sm_abs lfeng Wed Dec  9 16:47:55 EST 2009	                 %
% added modspace structure lfeng Thu Mar 19 17:32:29 SGT 2015                            %
% added InSAR los green functions Lgrn lfeng Tue Nov  3 14:14:35 SGT 2015                %
% last modified by Lujia Feng Tue Nov  3 14:16:57 SGT 2015                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xgrn0   = modspace.Xgrn;
Lgrn0   = modspace.Lgrn;
Bgrn0   = modspace.Bgrn;
Ngrn0   = modspace.Ngrn;
sm0     = modspace.sm;
sm_abs0 = modspace.sm_abs;
Aeq0    = modspace.Aeq;
beq0    = modspace.beq;
lb0     = modspace.lb;
ub0     = modspace.ub;
x00     = modspace.x0;

%%%%%%%%%%%%%%%%%%%% green %%%%%%%%%%%%%%%%%%%%
Xgrn = [ Xgrn0 Xgrn1 ];
Lgrn = [ Lgrn0 Lgrn1 ];
Bgrn = [ Bgrn0 Bgrn1 ];
Ngrn = [ Ngrn0 Ngrn1 ];

%%%%%%%%%%%%%%%%%%%%% sm  %%%%%%%%%%%%%%%%%%%%%
% sm =	[ sm0  0  ]
%	[  0  sm1 ]
[ sm ] = GTdef_add_diagonal(sm0,sm1);

%%%%%%%%%%%%%%%%%%%%% sm_abs  %%%%%%%%%%%%%%%%%%%%%
% sm =	[ sm0_abs  0  ]
%	[  0  sm1_abs ]
[ sm_abs ] = GTdef_add_diagonal(sm_abs0,sm_abs1);

%%%%%%%%%%%%%%%%%%%%% Aeq %%%%%%%%%%%%%%%%%%%%%
% Aeq =	[ Aeq0  0  ]
%	[  0  Aeq1 ]
[ Aeq ] = GTdef_add_diagonal(Aeq0,Aeq1);

%%%%%%%%%%%%%%%%%%%%% beq %%%%%%%%%%%%%%%%%%%%%
beq = [ beq0;beq1 ];

%%%%%%%%%%%%%%%%%%%%% lb  %%%%%%%%%%%%%%%%%%%%%
lb = [ lb0;lb1 ];

%%%%%%%%%%%%%%%%%%%%% ub  %%%%%%%%%%%%%%%%%%%%%
ub = [ ub0;ub1 ];

%%%%%%%%%%%%%%%%%%%%% x0  %%%%%%%%%%%%%%%%%%%%%
x0 = [ x00;x01 ];

modspace.Xgrn   = Xgrn;
modspace.Lgrn   = Lgrn;
modspace.Bgrn   = Bgrn;
modspace.Ngrn   = Ngrn;
modspace.sm     = sm;
modspace.sm_abs = sm_abs;
modspace.Aeq    = Aeq;
modspace.beq    = beq;
modspace.lb     = lb; 
modspace.ub     = ub; 
modspace.x0     = x0; 
