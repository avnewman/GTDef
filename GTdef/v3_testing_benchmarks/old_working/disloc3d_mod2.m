function [U,D,S,FLAG,FLAG2] = disloc3d_mod2(M,Xin,mu,nu)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UXT=0; UYT=0; UZT=0;

    UXXT=0; UXYT=0; UXZT=0;
    UYXT=0; UYYT=0; UYZT=0;
    UZXT=0; UZYT=0; UZZT=0;

%% CHECK THE ARGUMENTS
    IJ=size(M);
    I=IJ(1);
    NMOD=IJ(2);
    if(I~=10)
      disp('m must be 10x1 model vector');
    end

    IJ=size(Xin);
    I=IJ(1);
    NSTAT=IJ(2);
    if(I~=3)
      disp('x must be 3xn');
    end

    IJ=size(mu);
    I=IJ(1);
    J=IJ(2);
    if((I~=1)||(J~=1))
      disp('mu must be a scalar.');
    end

    IJ=size(nu);
    I=IJ(1);
    J=IJ(2);
    if((I~=1)||(J~=1))
      disp('nu must be a scalar.');
    end

%% CREATE A MATRIX FOR RETURN ARGUMENT
    for I=1:NSTAT
        for J=1:3
          U(J,I)=0;
        end
        for J=1:9
          D(J,I)=0;
        end
        for J=1:6
          S(J,I)=0;
        end
        FLAG(I)=0;
        for J=1:NMOD
          FLAG2(J,I)=0;
        end
    end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MU=mu;
    NU=nu;
    LAMBDA=2*MU*NU/(1-2*NU);
    DEG2RAD=3.14159265358979/180;
    ALPHA = (LAMBDA+MU)/(LAMBDA+2*MU);


%% LOOP OVER STATIONS
    for I=1:NSTAT
        for J=1:3
          STAT(J)=Xin(J,I);
        end
        if(STAT(3)>0)
          disp('Positive depth given.');
        else
            for J=1:NMOD
                for K=1:10
                  MODEL(K)=M(K,J);
                end
            
                STRIKE = (MODEL(5)-90)*DEG2RAD;
    		CS = cos(STRIKE);
    		SS = sin(STRIKE);
                DIP = MODEL(4);
                CD = cos(DIP*DEG2RAD);
                SD = sin(DIP*DEG2RAD);
                DISL1=MODEL(8);
                DISL2=MODEL(9);
                DISL3=MODEL(10);
                DEPTH=MODEL(3)-0.5*MODEL(2)*SD;
                AL1=MODEL(1)/2;
                AL2=AL1;
                AW1=MODEL(2)/2;
                AW2=AW1;
                X=CS*(-MODEL(6) + STAT(1)) - SS*(-MODEL(7) + STAT(2));
                Y=-0.5*CD*MODEL(2) + SS*(-MODEL(6) + STAT(1)) + CS*(-MODEL(7)+STAT(2));
                Z=STAT(3);
        
                if((MODEL(3)-SD*MODEL(2)<0)|(MODEL(1)<=0)|(MODEL(2)<=0)|(MODEL(3)<0))
                  disp('Unphysical model.');
                else
                [UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,FLAGOUT] =...
                    DC3D(ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3);
                
                FLAG2(J,I)=FLAGOUT+5+5*FLAGOUT;
                
                A = CS*UX +  SS*UY;
                B = -SS*UX + CS*UY;
                C = UZ;
                AAA = CS^2*UXX + CS*SS*(UXY + UYX) + SS^2*UYY;
                BBB = CS^2*UXY - SS^2*UYX + CS*SS*(-UXX + UYY);
                CCC = CS*UXZ + SS*UYZ;
                DDD = -(SS*(CS*UXX + SS*UXY)) + CS*(CS*UYX + SS*UYY);
                EEE = SS^2*UXX - CS*SS*(UXY + UYX) + CS^2*UYY;
                FFF = -(SS*UXZ) + CS*UYZ;
                GGG = CS*UZX + SS*UZY;
                HHH = -(SS*UZX) + CS*UZY;
                III = UZZ;
                
                UXT=UXT + A;
                UYT=UYT + B;
                UZT=UZT + C;
               
                UXXT=UXXT + AAA;
                UXYT=UXYT + BBB;
                UXZT=UXZT + CCC;
                UYXT=UYXT + DDD;
                UYYT=UYYT + EEE;
                UYZT=UYZT + FFF;
                UZXT=UZXT + GGG;
                UZYT=UZYT + HHH;
                UZZT=UZZT + III;
                end
            end
        
        UOUT=[UXT UYT UZT];
        DOUT=[UXXT UXYT UXZT UYXT UYYT UYZT UZXT UZYT UZZT];
        THETA=DOUT(1)+DOUT(5)+DOUT(9);
                           
        SOUT(1)=LAMBDA*THETA+2*MU*DOUT(1);
        SOUT(2)=MU*(DOUT(2)+DOUT(4));
        SOUT(3)=MU*(DOUT(3)+DOUT(7));
        SOUT(4)=LAMBDA*THETA+2*MU*DOUT(5);
        SOUT(5)=MU*(DOUT(6)+DOUT(8));
        SOUT(6)=LAMBDA*THETA+2*MU*DOUT(9);
        
        for J=1:3
            U(J,I)=UOUT(J);
        end
        for J=1:9
            D(J,I)=DOUT(J);
        end
        for J=1:6
            S(J,I)=SOUT(J);
        end

        FLAG(I)=FLAGOUT;
        
        UXT=0; UYT=0; UZT=0;

        UXXT=0; UXYT=0; UXZT=0;
        UYXT=0; UYYT=0; UYZT=0;
        UZXT=0; UZYT=0; UZZT=0;
        end
    end

  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET] = ...
    DC3D(ALPHA,X,Y,Z,DEPTH,DIP,AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3)

%********************************************************************
%*****                                                          *****
%*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****
%*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****    
%*****                         CODED BY  Y.OKADA ... SEP 1991   *****    
%*****                         REVISED   Y.OKADA ... NOV 1991   *****    
%*****                                                          *****    
%********************************************************************    
%                                                                        
%***** INPUT                                                             
%*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)            
%*****   X,Y,Z : COORDINATE OF OBSERVING POINT                           
%*****   DEPTH : SOURCE DEPTH                                            
%*****   DIP   : DIP-ANGLE (DEGREE)                                      
%*****   AL1,AL2   : FAULT LENGTH (-STRIKE,+STRIKE)                      
%*****   AW1,AW2   : FAULT WIDTH  ( DOWNDIP, UPDIP)                      
%*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS               
%                                                                        
%***** OUTPUT                                                            
%*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)                
%*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /              
%*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) ) 
%*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     
%*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )

    global SD CD XI2 ET2 Q2 R;
    F0=0.0;

    if(Z>0.0)
        disp('  POSITIVE Z WAS GIVEN IN SUB-DC3D  ');
    end
    
    for I=1:12
      U(I)=F0;
      DUA(I)=F0;
      DUB(I)=F0;
      DUC(I)=F0;
    end

    AALPHA=ALPHA;
    DDIP=DIP;
    DCCON0(AALPHA,DDIP);

%======================================                                  
%=====  REAL-SOURCE CONTRIBUTION  =====                                  
%======================================                                      
    
    D=DEPTH+Z;
    P=Y*CD+D*SD;
    Q=Y*SD-D*CD;
    JXI=0;
    JET=0;
    if((X+AL1)*(X-AL2)<=0) 
      JXI=1;
    end
    if((P+AW1)*(P-AW2)<=0)
      JET=1;
    end
    DD1=DISL1;
    DD2=DISL2;
    DD3=DISL3;

    for K=1:2
      if(K==1)
        ET=P+AW1;
      end
      if(K==2)
        ET=P-AW2;
      end
        for J=1:2
          if(J==1)
            XI=X+AL1;
          end
          if(J==2)
            XI=X-AL2;
          end

          DCCON2(XI,ET,Q,SD,CD);

	%=======================================                                 
	%=====  IN CASE OF SINGULAR (R=0)  =====                                 
	%=======================================                                 
          if((JXI==1) & (Q==F0) & (ET==F0))
            UX=F0;
            UY=F0;
            UZ=F0;
            UXX=F0;
            UYX=F0;
            UZX=F0;
            UXY=F0;
            UYY=F0;
            UZY=F0;
            UXZ=F0;
            UYZ=F0;
            UZZ=F0;
            IRET=1;
            return;
          end
          if((JET==1) & (Q==F0) & (XI==F0))
            UX=F0;
            UY=F0;
            UZ=F0;
            UXX=F0;
            UYX=F0;
            UZX=F0;
            UXY=F0;
            UYY=F0;
            UZY=F0;
            UXZ=F0;
            UYZ=F0;
            UZZ=F0;
            IRET=1;
            return;
          end
	%=======================================

          DUA=UA(XI,ET,Q,DD1,DD2,DD3);
          
          for I=1:3:10
            DU(I)=-DUA(I);
            DU(I+1)=-DUA(I+1)*CD+DUA(I+2)*SD;
            DU(I+2)=-DUA(I+1)*SD-DUA(I+2)*CD;
            if(I<10)
              continue
            else
               DU(I)=-DU(I);
               DU(I+1)=-DU(I+1);
               DU(I+2)=-DU(I+2);
             end
          end
            for I=1:12
              if((J+K)~=3)
                U(I)=U(I)+DU(I);
              else
                U(I)=U(I)-DU(I);
              end
           end
        end
    end      
    
%=======================================                                 
%=====  IMAGE-SOURCE CONTRIBUTION  =====                                 
%=======================================                                 

    ZZ=Z;
    D=DEPTH-Z;
    P=Y*CD+D*SD;
    Q=Y*SD-D*CD;
    JET=0;
    if((P+AW1)*(P-AW2)<=0)
      JET=1;
    end
    for K=1:2
      if(K==1)
        ET=P+AW1;
      else
        ET=P-AW2;
      end
        for J=1:2
          if(J==1)
            XI=X+AL1;
        else
            XI=X-AL2;
        end
        DCCON2(XI,ET,Q,SD,CD);
        DUA=UA(XI,ET,Q,DD1,DD2,DD3);
        DUB=UB(XI,ET,Q,DD1,DD2,DD3);
        DUC=UC(XI,ET,Q,ZZ,DD1,DD2,DD3);
        
            for I=1:3:10
              DU(I)=DUA(I)+DUB(I)+Z*DUC(I);
              DU(I+1)=(DUA(I+1)+DUB(I+1)+Z*DUC(I+1))*CD...
                  -(DUA(I+2)+DUB(I+2)+Z*DUC(I+2))*SD;
              DU(I+2)=(DUA(I+1)+DUB(I+1)-Z*DUC(I+1))*SD...
                  +(DUA(I+2)+DUB(I+2)-Z*DUC(I+2))*CD;
              if(I==10)
                DU(10)=DU(10)+DUC(1);
                DU(11)=DU(11)+DUC(2)*CD-DUC(3)*SD;
                DU(12)=DU(12)-DUC(2)*SD-DUC(3)*CD;
              end
            end
            for I=1:12
              if((J+K)~=3)
                U(I)=U(I)+DU(I);
              end
              if((J+K)==3)
                U(I)=U(I)-DU(I);
              end
            end
        end
    end

    UX=U(1);
    UY=U(2);
    UZ=U(3);
    UXX=U(4);
    UYX=U(5);
    UZX=U(6);
    UXY=U(7);
    UYY=U(8);
    UZY=U(9);
    UXZ=U(10);
    UYZ=U(11);
    UZZ=U(12);
    IRET=0;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function DCCON0(ALPHA,DIP)

%*******************************************************************     
%*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****     
%*******************************************************************     
%                                                                        
%***** INPUT                                                             
%*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)            
%*****   DIP   : DIP-ANGLE (DEGREE)                                      
%### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO    
%                                                                        

    global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D;
    F0=0.0;
    F1=1.0;
    F2=2.0;
    PI2=6.283185307179586;
    EPS=1.0E-6;

    ALP1=(F1-ALPHA)/F2;
    ALP2=ALPHA/F2;
    ALP3=(F1-ALPHA)/ALPHA;
    ALP4=F1-ALPHA;
    ALP5=ALPHA;

    P18=PI2/360.0;
    SD=sin(DIP*P18);
    CD=cos(DIP*P18);
    if(abs(CD)<EPS)
      CD=F0;
        if(SD>F0)
          SD=F1;
        end
        if(SD<F0)
          SD=-F1;
        end
    end

    SDSD=SD*SD;
    CDCD=CD*CD;
    SDCD=SD*CD;
    S2D=F2*SDCD;
    C2D=CDCD-SDSD;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = UA(XI,ET,Q,DISL1,DISL2,DISL3)

%********************************************************************   
%*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   
%*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   
%********************************************************************   
%                                                                       
%***** INPUT                                                            
%*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  
%*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              
%***** OUTPUT                                                           
%*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     
%                                                                       

    global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D;
    global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32 EY EZ FY FZ GY GZ HY HZ;
    F0=0.0;
    F2=2.0;
    PI2=6.283185307179586;
    for I=1:12
      U(I)=F0;
    end
    XY=XI*Y11;
    QX=Q*X11;
    QY=Q*Y11;

%======================================                                  
%=====  STRIKE-SLIP CONTRIBUTION  =====                                  
%======================================                                  

    if(DISL1~=F0) 
        DU( 1)=    TT/F2 +ALP2*XI*QY;
        DU( 2)=           ALP2*Q/R;
        DU( 3)= ALP1*ALE -ALP2*Q*QY;
        DU( 4)=-ALP1*QY  -ALP2*XI2*Q*Y32;
        DU( 5)=          -ALP2*XI*Q/R3;
        DU( 6)= ALP1*XY  +ALP2*XI*Q2*Y32;
        DU( 7)= ALP1*XY*SD        +ALP2*XI*FY+D/F2*X11;
        DU( 8)=                    ALP2*EY;
        DU( 9)= ALP1*(CD/R+QY*SD) -ALP2*Q*FY;
        DU(10)= ALP1*XY*CD        +ALP2*XI*FZ+Y/F2*X11;
        U(11)=                    ALP2*EZ;
        DU(12)=-ALP1*(SD/R-QY*CD) -ALP2*Q*FZ;
        for I=1:12
          U(I)=U(I)+DISL1/PI2*DU(I);
        end
    end

%======================================                                  
%=====    DIP-SLIP CONTRIBUTION   =====                                  
%======================================                                  

if(DISL2~=F0)
    DU( 1)=           ALP2*Q/R;
    DU( 2)=    TT/F2 +ALP2*ET*QX;
    DU( 3)= ALP1*ALX -ALP2*Q*QX;
    DU( 4)=        -ALP2*XI*Q/R3;
    DU( 5)= -QY/F2 -ALP2*ET*Q/R3;
    DU( 6)= ALP1/R +ALP2*Q2/R3;
    DU( 7)=                      ALP2*EY;
    DU( 8)= ALP1*D*X11+XY/F2*SD +ALP2*ET*GY;
    DU( 9)= ALP1*Y*X11          -ALP2*Q*GY;
    DU(10)=                      ALP2*EZ;
    DU(11)= ALP1*Y*X11+XY/F2*CD +ALP2*ET*GZ;
    DU(12)=-ALP1*D*X11          -ALP2*Q*GZ;
    for I=1:12
        U(I)=U(I)+DISL2/PI2*DU(I);
    end
end

%========================================                                
%=====  TENSILE-FAULT CONTRIBUTION  =====                                
%========================================                                

    if(DISL3~=F0)
        DU( 1)=-ALP1*ALE -ALP2*Q*QY;
        DU( 2)=-ALP1*ALX -ALP2*Q*QX;
        DU( 3)=    TT/F2 -ALP2*(ET*QX+XI*QY);
        DU( 4)=-ALP1*XY  +ALP2*XI*Q2*Y32;
        DU( 5)=-ALP1/R   +ALP2*Q2/R3;
        DU( 6)=-ALP1*QY  -ALP2*Q*Q2*Y32;
        DU( 7)=-ALP1*(CD/R+QY*SD)  -ALP2*Q*FY;
        DU( 8)=-ALP1*Y*X11         -ALP2*Q*GY;
        DU( 9)= ALP1*(D*X11+XY*SD) +ALP2*Q*HY;
        DU(10)= ALP1*(SD/R-QY*CD)  -ALP2*Q*FZ;
        DU(11)= ALP1*D*X11         -ALP2*Q*GZ;
        DU(12)= ALP1*(Y*X11+XY*CD) +ALP2*Q*HZ;
        for I=1:12
          U(I)=U(I)+DISL3/PI2*DU(I);
        end
    end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = UB(XI,ET,Q,DISL1,DISL2,DISL3)

%********************************************************************    
%*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****    
%*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****    
%********************************************************************     
%                                                                        
%***** INPUT                                                             
%*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                   
%*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS                
%***** OUTPUT                                                            
%*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                      
%

    global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D;
    global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32 EY EZ FY FZ GY GZ HY HZ;
    F0=0.0;
    F1=1.0;
    F2=2.0;
    PI2=6.283185307179586;

    RD=R+D;
    D11=F1/(R*RD);
    AJ2=XI*Y/RD*D11;
    AJ5=-(D+Y*Y/RD)*D11;
    if(CD~=F0)
        if(XI==F0)
          AI4=F0;
        else
        X=sqrt(XI2+Q2);
        AI4=F1/CDCD*(XI/RD*SDCD+F2*atan((ET*(X+Q*CD)+X*(R+X)*SD)/(XI*(R+X)*CD)));
        end
    AI3=(Y*CD/RD-ALE+SD*log(RD))/CDCD;
    AK1=XI*(D11-Y11*SD)/CD;
    AK3=(Q*Y11-Y*D11)/CD;
    AJ3=(AK1-AJ2*SD)/CD;
    AJ6=(AK3-AJ5*SD)/CD;
    else
    RD2=RD*RD;
    AI3=(ET/RD+Y*Q/RD2-ALE)/F2;
    AI4=XI*Y/RD2/F2;
    AK1=XI*Q/RD*D11;
    AK3=SD/RD*(XI2*D11-F1);
    AJ3=-XI/RD2*(Q2*D11-F1/F2);
    AJ6=-Y/RD2*(XI2*D11-F1/F2);
    end

    XY=XI*Y11;
    AI1=-XI/RD*CD-AI4*SD;
    AI2= log(RD)+AI3*SD;
    AK2= F1/R+AK3*SD;
    AK4= XY*CD-AK1*SD;
    AJ1= AJ5*CD-AJ6*SD;
    AJ4=-XY-AJ2*CD+AJ3*SD;

    for I=1:12
      U(I)=F0;
    end

    QX=Q*X11;
    QY=Q*Y11;

%======================================                                  
%=====  STRIKE-SLIP CONTRIBUTION  =====                                  
%======================================                                  
    if(DISL1~=F0)
        DU( 1)=-XI*QY-TT -ALP3*AI1*SD;
        DU( 2)=-Q/R      +ALP3*Y/RD*SD;
        DU( 3)= Q*QY     -ALP3*AI2*SD;
        DU( 4)= XI2*Q*Y32 -ALP3*AJ1*SD;
        DU( 5)= XI*Q/R3   -ALP3*AJ2*SD;
        DU( 6)=-XI*Q2*Y32 -ALP3*AJ3*SD;
        DU( 7)=-XI*FY-D*X11 +ALP3*(XY+AJ4)*SD;
        DU( 8)=-EY          +ALP3*(F1/R+AJ5)*SD;
        DU( 9)= Q*FY        -ALP3*(QY-AJ6)*SD;
        DU(10)=-XI*FZ-Y*X11 +ALP3*AK1*SD;
        DU(11)=-EZ          +ALP3*Y*D11*SD;
        DU(12)= Q*FZ        +ALP3*AK2*SD;
        for I=1:12
          U(I)=U(I)+DISL1/PI2*DU(I);
        end
    end

%======================================                                  
%=====    DIP-SLIP CONTRIBUTION   =====                                  
%======================================                                  

    if(DISL2~=F0)
       DU( 1)=-Q/R      +ALP3*AI3*SDCD;
        DU( 2)=-ET*QX-TT -ALP3*XI/RD*SDCD;
        DU( 3)= Q*QX     +ALP3*AI4*SDCD;
        DU( 4)= XI*Q/R3     +ALP3*AJ4*SDCD;
        DU( 5)= ET*Q/R3+QY  +ALP3*AJ5*SDCD;
        DU( 6)=-Q2/R3       +ALP3*AJ6*SDCD;
        DU( 7)=-EY          +ALP3*AJ1*SDCD;
        DU( 8)=-ET*GY-XY*SD +ALP3*AJ2*SDCD;
        DU( 9)= Q*GY        +ALP3*AJ3*SDCD;
        DU(10)=-EZ          -ALP3*AK3*SDCD;
        DU(11)=-ET*GZ-XY*CD -ALP3*XI*D11*SDCD;
        DU(12)= Q*GZ        -ALP3*AK4*SDCD;
        for I=1:12
          U(I)=U(I)+DISL2/PI2*DU(I);
        end
    end

%========================================                                
%=====  TENSILE-FAULT CONTRIBUTION  =====                                
%========================================                                

    if(DISL3~=F0)
        DU( 1)= Q*QY           -ALP3*AI3*SDSD;
        DU( 2)= Q*QX           +ALP3*XI/RD*SDSD;
        DU( 3)= ET*QX+XI*QY-TT -ALP3*AI4*SDSD;
        DU( 4)=-XI*Q2*Y32 -ALP3*AJ4*SDSD;
        DU( 5)=-Q2/R3     -ALP3*AJ5*SDSD;
        DU( 6)= Q*Q2*Y32  -ALP3*AJ6*SDSD;
        DU( 7)= Q*FY -ALP3*AJ1*SDSD;
        DU( 8)= Q*GY -ALP3*AJ2*SDSD;
        DU( 9)=-Q*HY -ALP3*AJ3*SDSD;
        DU(10)= Q*FZ +ALP3*AK3*SDSD;
        DU(11)= Q*GZ +ALP3*XI*D11*SDSD;
        DU(12)=-Q*HZ +ALP3*AK4*SDSD;
        for I=1:12
          U(I)=U(I)+DISL3/PI2*DU(I);
        end
    end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function U = UC(XI,ET,Q,Z,DISL1,DISL2,DISL3)

%********************************************************************    
%*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   
%*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****    
%********************************************************************    
%                                                                        
%***** INPUT                                                             
%*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM               
%*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS               
%***** OUTPUT                                                            
%*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                      

    global ALP1 ALP2 ALP3 ALP4 ALP5 SD CD SDSD CDCD SDCD S2D C2D;
    global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32 EY EZ FY FZ GY GZ HY HZ;
    F0=0.0;
    F1=1.0;
    F2=2.0;
    F3=3.0;
    PI2=6.283185307179586;

    C=D+Z;                                        
    X53=(8.0*R2+9.0*R*XI+F3*XI2)*X11*X11*X11/R2;
    Y53=(8.0*R2+9.0*R*ET+F3*ET2)*Y11*Y11*Y11/R2;
    H=Q*CD-Z;
    Z32=SD/R3-H*Y32;
    Z53=F3*SD/R5-H*Y53;
    Y0=Y11-XI2*Y32;
    Z0=Z32-XI2*Z53;
    PPY=CD/R3+Q*Y32*SD;
    PPZ=SD/R3-Q*Y32*CD;
    QQ=Z*Y32+Z32+Z0;
    QQY=F3*C*D/R5-QQ*SD;
    QQZ=F3*C*Y/R5-QQ*CD+Q*Y32;
    XY=XI*Y11;
    QX=Q*X11;
    QY=Q*Y11;    
    QR=F3*Q/R5;
    CQX=C*Q*X53;
    CDR=(C+D)/R3;
    YY0=Y/R3-Y0*CD;
    for I=1:12
      U(I)=F0;
    end

%======================================                                  
%=====  STRIKE-SLIP CONTRIBUTION  =====                                  
%======================================                                  

    if(DISL1~=F0)
        DU( 1)= ALP4*XY*CD           -ALP5*XI*Q*Z32;
        DU( 2)= ALP4*(CD/R+F2*QY*SD) -ALP5*C*Q/R3;
        DU( 3)= ALP4*QY*CD           -ALP5*(C*ET/R3-Z*Y11+XI2*Z32);
        DU( 4)= ALP4*Y0*CD                  -ALP5*Q*Z0;
        DU( 5)=-ALP4*XI*(CD/R3+F2*Q*Y32*SD) +ALP5*C*XI*QR;
        DU( 6)=-ALP4*XI*Q*Y32*CD            +ALP5*XI*(F3*C*ET/R5-QQ);
        DU( 7)=-ALP4*XI*PPY*CD    -ALP5*XI*QQY;
        DU( 8)= ALP4*F2*(D/R3-Y0*SD)*SD-Y/R3*CD-ALP5*(CDR*SD-ET/R3-C*Y*QR);
        DU( 9)=-ALP4*Q/R3+YY0*SD  +ALP5*(CDR*CD+C*D*QR-(Y0*CD+Q*Z0)*SD);
        DU(10)= ALP4*XI*PPZ*CD    -ALP5*XI*QQZ;
        DU(11)= ALP4*F2*(Y/R3-Y0*CD)*SD+D/R3*CD -ALP5*(CDR*CD+C*D*QR);
        DU(12)=         YY0*CD    -ALP5*(CDR*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD);
        for I=1:12
          U(I)=U(I)+DISL1/PI2*DU(I);
        end
    end

%======================================                                  
%=====    DIP-SLIP CONTRIBUTION   =====                                  
%======================================                                  

    if(DISL2~=F0)
        DU( 1)= ALP4*CD/R -QY*SD -ALP5*C*Q/R3;
        DU( 2)= ALP4*Y*X11       -ALP5*C*ET*Q*X32;
        DU( 3)=     -D*X11-XY*SD -ALP5*C*(X11-Q2*X32);
        DU( 4)=-ALP4*XI/R3*CD +ALP5*C*XI*QR +XI*Q*Y32*SD;
        DU( 5)=-ALP4*Y/R3     +ALP5*C*ET*QR;
        DU( 6)=    D/R3-Y0*SD +ALP5*C/R3*(F1-F3*Q2/R2);
        DU( 7)=-ALP4*ET/R3+Y0*SDSD -ALP5*(CDR*SD-C*Y*QR);
        DU( 8)= ALP4*(X11-Y*Y*X32) -ALP5*C*((D+F2*Q*CD)*X32-Y*ET*Q*X53);
        DU( 9)=  XI*PPY*SD+Y*D*X32 +ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53);
        DU(10)=      -Q/R3+Y0*SDCD -ALP5*(CDR*CD+C*D*QR);
        DU(11)= ALP4*Y*D*X32       -ALP5*C*((Y-F2*Q*SD)*X32+D*ET*Q*X53);
        DU(12)=-XI*PPZ*SD+X11-D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53);
        for I=1:12
          U(I)=U(I)+DISL2/PI2*DU(I);
        end
    end

%========================================                                
%=====  TENSILE-FAULT CONTRIBUTION  =====                                
%========================================                                

    if(DISL3~=F0)
        DU( 1)=-ALP4*(SD/R+QY*CD)   -ALP5*(Z*Y11-Q2*Z32);
        DU( 2)= ALP4*F2*XY*SD+D*X11 -ALP5*C*(X11-Q2*X32);
        DU( 3)= ALP4*(Y*X11+XY*CD)  +ALP5*Q*(C*ET*X32+XI*Z32);
        DU( 4)= ALP4*XI/R3*SD+XI*Q*Y32*CD+ALP5*XI*(F3*C*ET/R5-F2*Z32-Z0);
        DU( 5)= ALP4*F2*Y0*SD-D/R3 +ALP5*C/R3*(F1-F3*Q2/R2);
        DU( 6)=-ALP4*YY0           -ALP5*(C*ET*QR-Q*Z0);
        DU( 7)= ALP4*(Q/R3+Y0*SDCD)   +ALP5*(Z/R3*CD+C*D*QR-Q*Z0*SD);
        DU( 8)=-ALP4*F2*XI*PPY*SD-Y*D*X32+ALP5*C*((Y+F2*Q*SD)*X32-Y*Q2*X53);
        DU( 9)=-ALP4*(XI*PPY*CD-X11+Y*Y*X32)+ALP5*(C*((D+F2*Q*CD)*X32-Y*ET*Q*X53)+XI*QQY);
        DU(10)=  -ET/R3+Y0*CDCD -ALP5*(Z/R3*SD-C*Y*QR-Y0*SDSD+Q*Z0*CD);
        DU(11)= ALP4*F2*XI*PPZ*SD-X11+D*D*X32-ALP5*C*((D-F2*Q*CD)*X32-D*Q2*X53);
        DU(12)= ALP4*(XI*PPZ*CD+Y*D*X32)+ALP5*(C*((Y-F2*Q*SD)*X32+D*ET*Q*X53)+XI*QQZ);
        for I=1:12
          U(I)=U(I)+DISL3/PI2*DU(I);
        end
    end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function DCCON2(XI,ET,Q,SD,CD)

%**********************************************************************  
%*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   *****  
%**********************************************************************  
%                                                                        
%***** INPUT                                                             
%*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                   
%*****   SD,CD   : SIN, COS OF DIP-ANGLE                                 
%                                                                       
%### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZER0 
%

    global XI2 ET2 Q2 R R2 R3 R5 Y D TT ALX ALE X11 Y11 X32 Y32 EY EZ FY FZ GY GZ HY HZ;
    F0=0.0;
    F1=1.0;
    F2=2.0;
    EPS=1.0E-6;

    if(abs(XI)<EPS)
      XI=F0;
    end
    if(abs(ET)<EPS)
      ET=F0;
    end
    if(abs(Q)<EPS)
      Q=F0;
    end
    XI2=XI*XI;
    ET2=ET*ET;
    Q2=Q*Q;
    R2=XI2+ET2+Q2;
    R=sqrt(R2);
    if(R==F0)
      return;
    end
    R3=R*R2;
    R5=R3*R2;
    Y=ET*CD+Q*SD;
    D=ET*SD-Q*CD;
        
    if(Q==F0)
      TT=F0;
    else
      TT=atan(XI*ET/(Q*R));
    end

    if((XI<F0) & (Q==F0) & (ET==F0))
      ALX=-log(R-XI);
      X11=F0;
      X32=F0;
    else
      RXI=R+XI;
      ALX=log(RXI);
      X11=F1/(R*RXI);
      X32=(R+RXI)*X11*X11/R;
    end

    if((ET<F0) & (Q==F0) & (XI==F0))
      ALE=-log(R-ET);
      Y11=F0;
      Y32=F0;
    else
      RET=R+ET;
      ALE=log(RET);
      Y11=F1/(R*RET);
      Y32=(R+RET)*Y11*Y11/R;
    end

    EY=SD/R-Y*Q/R3;
    EZ=CD/R+D*Q/R3;
    FY=D/R3+XI2*Y32*SD;
    FZ=Y/R3+XI2*Y32*CD;
    GY=F2*X11*SD-Y*Q*X32;
    GZ=F2*X11*CD+D*Q*X32;
    HY=D*Q*X32+XI*Q*Y32*SD;
    HZ=Y*Q*X32+XI*Q*Y32*CD;
