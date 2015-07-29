function [ Xgrn,Bgrn,Ngrn,sm,Aeq,beq,lb,ub,x0 ]... 
         = GTdef_addall(Xgrn0,Bgrn0,Ngrn0,sm0,Aeq0,beq0,lb0,ub0,x00,...
	                Xgrn1,Bgrn1,Ngrn1,sm1,Aeq1,beq1,lb1,ub1,x01)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_addall				  %
% Add new green1,Aeq1,beq1,lb1,ub1,x01 to pre-existing green0,Aeq0,beq0,  %
% lb0,ub0,x00 for x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0)			  %
% Here we have no inequalities, so set A=[];b=[]			  %
%									  %
% PARAMETERS:					  		  	  %
%  Xgrn - displacements [east;north;vertical] for different sites   	  %
%         from unit slips [(3*nn)*slip_num] 				  %
%         (nn is the  number of sites)  				  %
%  Bgrn - length changes [east;north;vertical;absolute] for 	  	  %
%         different baselines from unit slips [(4*nn)*slip_num] 	  %
%         (nn is the  number of baselines)  				  %
%  Ngrn - displacements [east;north;vertical] for different nodes   	  %
%         from unit slips [(3*nn)*slip_num] 				  %
%         (nn is the  number of nodes)  				  %
%  sm  - smoothing matrix						  %
%  Aeq - left-hand side matrix for linear equalities  [slip_num*slip_num] %
%  beq - right-hand side vector for linear equalities [slip_num*1]        %
%  x0  - initial values for ss,ds,ts 	[slip_num*1]                      %
%  lb  - lower bounds for ss,ds,ts 	[slip_num*1]                      %
%  ub  - upper bounds for ss,ds,ts	[slip_num*1]			  %
%                                                                         %
% first created by Lujia Feng Fri Apr 24 14:38:50 EDT 2009		  %
% last modified by Lujia Feng Fri May  1 19:28:53 EDT 2009		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% green %%%%%%%%%%%%%%%%%%%%
Xgrn = [ Xgrn0 Xgrn1 ];
Bgrn = [ Bgrn0 Bgrn1 ];
Ngrn = [ Ngrn0 Ngrn1 ];

%%%%%%%%%%%%%%%%%%%%% sm  %%%%%%%%%%%%%%%%%%%%%
% sm =	[ sm0  0  ]
%	[  0  sm1 ]
[m0 n0] = size(sm0);    [m1 n1] = size(sm1);
upright = zeros(m0,n1); botleft = zeros(m1,n0);
sm = [ sm0 upright; botleft sm1 ];

%%%%%%%%%%%%%%%%%%%%% Aeq %%%%%%%%%%%%%%%%%%%%%
% Aeq =	[ Aeq0  0  ]
%	[  0  Aeq1 ]
[m0 n0] = size(Aeq0);   [m1 n1] = size(Aeq1);
upright = zeros(m0,n1); botleft = zeros(m1,n0);
Aeq = [ Aeq0 upright; botleft Aeq1 ];

%%%%%%%%%%%%%%%%%%%%% beq %%%%%%%%%%%%%%%%%%%%%
beq = [ beq0;beq1 ];

%%%%%%%%%%%%%%%%%%%%% lb  %%%%%%%%%%%%%%%%%%%%%
lb = [ lb0;lb1 ];

%%%%%%%%%%%%%%%%%%%%% ub  %%%%%%%%%%%%%%%%%%%%%
ub = [ ub0;ub1 ];

%%%%%%%%%%%%%%%%%%%%% x0  %%%%%%%%%%%%%%%%%%%%%
x0 = [ x00;x01 ];
