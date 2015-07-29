function [ sm_1d3pf,sm_1d3pb,sm_2d,sm_abs ] = GTdef_sm_fixed(dd,ds,Nd,Ns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            GTdef_sm_fixed				  %
% Create the smoothing matrices for the four types of smoothing methods   %
% Motion on surface is fixed						  %
% INPUT:					  		  	  %
%   dd - distance along dip between two vertically adjacent patches	  %
%   ds - distance along strike between two horizontally adjacent patches  %
%   Nd - number of patches along dip				          %
%   Ns - number of patches along strike					  %
%									  %
% OUTPUT:								  %
%   sm_1d2pc,sm_1d3pf,sm_1d3pb,sm_2d 		[3nn*3nn]  (nn = Ns*Nd)	  %
%   take into account 3 slips						  %
%									  %
% first created by Lujia Feng Thu Dec  3 00:23:51 EST 2009		  %
% changed from GTdef_sm.m to GTdef_sm_fixed.m lfeng Wed Feb 24  2010	  %
% last modified by Lujia Feng Wed Feb 24 12:54:08 EST 2010		  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% First: create smoothing matrix for one type of slip
% Then:  create smoothing matrix for three types of slips  		    
%     [ sm0  0  0  ]
%sm = [  0  sm0 0  ]
%     [  0  0  sm0 ]

%%%%%% 2-point center %%%%%%
%[ sm0_1d2pc ] = GTdef_sm1d_2pctr(dd,ds,Nd,Ns);  
%[ sm_1d2pc ] = GTdef_add_diagonal(sm0_1d2pc,sm0_1d2pc);		       
%[ sm_1d2pc ] = GTdef_add_diagonal(sm_1d2pc,sm0_1d2pc);

%%%%%%% 8-point center %%%%%%
%[ sm0_8_1 ] = GTdef_sm1d_8pctr_uprt(dd,ds,Nd,Ns);  
%[ sm0_8_2 ] = GTdef_sm1d_8pctr_rtdw(dd,ds,Nd,Ns);  
%[ sm_8_1 ] = GTdef_add_diagonal(sm0_8_1,sm0_8_1);
%[ sm_8_1 ] = GTdef_add_diagonal(sm_8_1,sm0_8_1);
%[ sm_8_2 ] = GTdef_add_diagonal(sm0_8_2,sm0_8_2);
%[ sm_8_2 ] = GTdef_add_diagonal(sm_8_2,sm0_8_2);
%sm_1d8pc = [ sm_8_1;sm_8_2 ];

%%%%%% 3-point forward %%%%%%
[ sm0_3f_1 ] = GTdef_sm1d_3pfwd_uprt(dd,ds,Nd,Ns);  
[ sm0_3f_2 ] = GTdef_sm1d_3pfwd_rtdw(dd,ds,Nd,Ns);  

[ sm_3f_1 ] = GTdef_add_diagonal(sm0_3f_1,sm0_3f_1);		     
[ sm_3f_1 ] = GTdef_add_diagonal(sm_3f_1,sm0_3f_1);
[ sm_3f_2 ] = GTdef_add_diagonal(sm0_3f_2,sm0_3f_2);		     
[ sm_3f_2 ] = GTdef_add_diagonal(sm_3f_2,sm0_3f_2);
sm_1d3pf = [ sm_3f_1;sm_3f_2 ];

%%%%%% 3-point backward %%%%%%
[ sm0_3b_1 ] = GTdef_sm1d_3pbwd_uprt(dd,ds,Nd,Ns);  
[ sm0_3b_2 ] = GTdef_sm1d_3pbwd_rtdw(dd,ds,Nd,Ns);  

[ sm_3b_1 ] = GTdef_add_diagonal(sm0_3b_1,sm0_3b_1);		     
[ sm_3b_1 ] = GTdef_add_diagonal(sm_3b_1,sm0_3b_1);
[ sm_3b_2 ] = GTdef_add_diagonal(sm0_3b_2,sm0_3b_2);		     
[ sm_3b_2 ] = GTdef_add_diagonal(sm_3b_2,sm0_3b_2);
sm_1d3pb = [ sm_3b_1;sm_3b_2 ];

%%%%%% 2-point center %%%%%%
[ sm0_2d ] = GTdef_sm2d(dd,ds,Nd,Ns);  
[ sm_2d ] = GTdef_add_diagonal(sm0_2d,sm0_2d);		     
[ sm_2d ] = GTdef_add_diagonal(sm_2d,sm0_2d);

%%%%%% strain sm %%%%%%
[ sm0_abs ] = GTdef_strain(dd,ds,Nd,Ns);  
[ sm_abs ] = GTdef_add_diagonal(sm0_abs,sm0_abs);		     
[ sm_abs ] = GTdef_add_diagonal(sm_abs,sm0_abs);
