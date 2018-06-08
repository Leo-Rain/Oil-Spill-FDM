      PROGRAM Hydrodynamic
	  !    Hydrodynamic & OIL Dynamic 
      IMPLICIT NONE 

    REAL  Depth,DEN,NUM,CP,Z0,vol,Cd,ZX,VolumOIL,HWave,Period,Slope,Rip
	REAL  Dx,Dy,Dt,Velx,Vely,T,ADVx,ADVy,Gravity,Omega
   	REAL  B,W,PP,G,Pw,Pa,Po,Pi
	REAL  Crx,Cry,Thickness,Cfilm,MIN,CENTERx,LENGTH,CENTERy,L,SIGMA
	REAL  AA,BB,CC,DD,EE,SSin,SSout,Beta,K,TETA,Ce
	REAL,ALLOCATABLE:: A1(:,:),B1(:,:),C1(:,:),D1(:,:),E1(:,:),F1(:,:),G1(:,:),H1(:,:)
	REAL,ALLOCATABLE:: XX(:),BX(:),AX(:,:),CX1(:,:),MX1(:,:),MX2(:,:)
	REAL,ALLOCATABLE:: Vw(:,:),Twx(:,:),Twy(:,:),WINDteta(:,:),LANDA(:,:),EVAP(:,:)
	REAL,ALLOCATABLE:: Z(:,:),H(:,:),D(:,:),U(:,:),V(:,:),WINDx(:,:),WINDy(:,:),Friction(:,:)
	REAL,ALLOCATABLE:: ADVv(:,:),DIFx(:,:),DIFy(:,:),Vp(:,:),Up(:,:),Coriolis(:,:)
	REAL,ALLOCATABLE:: AAA(:),BBB(:),CCC(:),DDD(:),SI(:),CX(:),C(:),CXX(:)
	REAL,ALLOCATABLE:: UU(:,:),VV(:,:),ZZ(:,:),HH(:,:),S(:,:),SX(:,:),WET(:,:),Ah(:,:),Cchezy(:,:)
	INTEGER X,Y,Q,M,N,I,J,P,R
	INTEGER Li,Lj,Ri,Rj,NUMCOORDINATE
     
	OPEN (10, FILE="INPUT.TXT")
	OPEN (20, FILE="TOTAL.DAT")
	OPEN (30, FILE="Boundary.TXT",status="unknown")	
	OPEN (40, FILE="WET.TXT",status="unknown")
	OPEN (50, FILE="Island.TXT",status="unknown")
	OPEN (60, FILE="DRY.TXT",status="unknown")
	OPEN (70, FILE="Bathymetry.TXT",status="unknown")	
	READ (10,*) B,K,Ce,G,Pi,T,X,Y,Dx,Dy,Dt,Velx,Vely,Depth,Z0,Pa,Pw,Po,Thickness
	READ (10,*) VolumOIL,Cfilm,MIN,CENTERx,CENTERy,LENGTH,SSin,SSout,HWave,Period,Slope,Rip

    ALLOCATE (XX(MAX(2*X,2*Y)),BX(MAX(2*X,2*Y)),AX(MAX(2*X,2*Y),MAX(2*X,2*Y)))
	ALLOCATE (CX1(X,Y),MX1(X,Y),MX2(X,Y),ADVv(X,Y),DIFx(X,Y),DIFy(X,Y),CXX(MAX(X,Y)))
	ALLOCATE (Vw(X,Y),Twx(X,Y),Twy(X,Y),WINDteta(X,Y),LANDA(X,Y),EVAP(X,Y))
	ALLOCATE (Vp(X,Y),Up(X,Y),Coriolis(X,Y),WINDx(X,Y),WINDy(X,Y),Friction(X,Y))
	ALLOCATE (Z(X,Y),H(X,Y),U(X,Y),V(X,Y),D(X,Y),ZZ(X,Y),HH(X,Y),UU(X,Y),VV(X,Y),S(X,Y),SX(X,Y)) 
	ALLOCATE (AAA(MAX(X,Y)),BBB(MAX(X,Y)),CCC(MAX(X,Y)),DDD(MAX(X,Y)),SI(MAX(X,Y)),CX(MAX(X,Y)),C(MAX(X,Y)))
	ALLOCATE (A1(X,Y),B1(X,Y),C1(X,Y),D1(X,Y),E1(X,Y),F1(X,Y),G1(X,Y),H1(X,Y),WET(X,Y),Ah(X,Y),Cchezy(X,Y))

      WRITE(20,*), 'Title="Hydrodynamic"'
      WRITE(20,*), 'variables="i","j","z","s","h","d","u","v","w"'

	  !******** Initial & Boundary Condition ***********

!     WET=0
	  
!	  DO I=1,X

!	  READ (30,*) Li,Lj,Ri,Rj
	  
!	  DO J=Lj,Rj

	  WET=1
	  
!	  END DO
!	  END DO

!	  READ (50,*) NUMCOORDINATE

!	  DO I=1,NUMCOORDINATE

!	  READ (50,*) Li,Lj,Ri,Rj
	  
!	  DO J=Lj,Rj

!	  WET(Li,J)=0
	  
!	  END DO
!	  END DO

!	  DO I=1,X 
	  	  
!	  WRITE (40,'(50(I1))'), INT(WET(I,1:Y))
 
!	  END DO

!	  CLOSE (40)

      DO P=1,X*Y 

	  READ (70,*) I,J,D(I,J)

	  END DO
 	  	  
	  DO I=1,X 
	  DO J=1,Y

	  U(I,J)=0.0
	  V(I,J)=0.0
	  UU(I,J)=0.0
	  VV(I,J)=0.0
      Z(I,J)=Z0
	  H(I,J)=D(I,J)+Z(I,J)
	  IF (H(I,J)<0.00) H(I,J)=0
	  IF (H(I,J)==0)   WET(I,J)=0 
	  IF (H(I,J)>0)    WET(I,J)=1
!	  IF (WET(I,J)==0) H(I,J)=0 
	  IF (WET(I,J)==0) S(I,J)=0 
	  IF (WET(I,J)==0) U(I,J)=0  
	  IF (WET(I,J)==0) V(I,J)=0
	  L=INT(LENGTH/2)
      SIGMA=SQRT(-(L*Dx)**2/(2*LOG(MIN)))
!	  LANDA=0.2055e-3*Dx**1.15
	  S(I,J)=0 
	  SX(I,J)=0
	  EVAP(I,J)=0

	  END DO 
	  END DO

	  PRINT*, H(1,1)

	  DO J=25,31
	  DO I=116,123 

	  S(I,J)=0.01
	  
	  END DO
	  END DO

	  W=2.*Pi/(24.*3600.) 
	  Omega=2.*W*SIN(Pi*B)
	  
	  DO I=1,X
	  DO J=1,Y

	  Cd=(0.63+0.066*Vw(I,J))*1.E-3
	  Twx(I,J)=Cd*Pa*ABS(Vw(I,J))*Vw(I,J)*COS(Pi*WINDteta(I,J))
	  Twy(I,J)=Cd*Pa*ABS(Vw(I,J))*Vw(I,J)*SIN(Pi*WINDteta(I,J))

	  END DO 
	  END DO
     
	  !******** Main Loop ***********

	  DO P=1,(T/dT)

	  Q=MOD(P,100)
	  IF (Q.eq.0.or.P.eq.1) THEN ;	  PP=PP+1
      WRITE(20,*), 'zone t="zone number-' ,PP,'",i=',X,',j=',Y,''
      WRITE(20,*), 'f=point'

	  DO J=1,Y  ;	DO I=1,X

	      WRITE(20,*),I,J,Z(I,J),S(I,J),H(I,J),D(I,J),U(I,J),V(I,J),0

      END DO;	  END DO;	  END IF

	  
	  !************** Hydrodynamic *************
	  
	  !********* Alternating Directin Implicit *******

	  !******** Initial & Boundary Condition ***********

		DO J=1,Y
	    DO I=1,X

		Z(X,J)=HWave/2.*COS(-2*Pi*P*dT/Period)
		H(I,J)=D(I,J)+Z(I,J)

		Z(I,1)=0.1*Z(I,2)
		Z(I,Y)=0.1*Z(I,Y-1)		
		U(I,1)=0.1*U(I,2)
		U(I,Y)=0.1*U(I,Y-1)
		V(I,1)=0.1*V(I,2)
		V(I,Y)=0.1*V(I,Y-1)
!		V(I,Y-1)=0 !.1*V(I,Y-2)

		IF (H(I,J)<0.00) H(I,J)=0
        IF (H(I,J)==0)   WET(I,J)=0 
	    IF (H(I,J)>0)    WET(I,J)=1
!		IF (WET(I,J)==0) H(I,J)=0
		IF (WET(I,J)==0) S(I,J)=0
		IF (WET(I,J)==0) U(I,J)=0
		IF (WET(I-1,J)==0.and.WET(I,J)==1) U(I-1,J)=U(I,J)
!		IF (WET(I,J)==0) V(I,J)=0
		IF (H(I,J)==0)   Ah(I,J)=0  

		IF (WET(I,J)==1) THEN

		Cchezy(I,J)=-18*LOG10(K/(H(I,J)*12))
	    Ah(I,J)=0.0376 !Ce*H(I,J)*SQRT(G*(U(I,J)**2+V(I,J)**2))/Cchezy(I,J)

		END IF
 

		END DO
		END DO

		advv=0
		vp=0
		up=0
		WINDx=0
		WINDy=0
		Coriolis=0
		DIFx=0
		DIFy=0
		Friction=0

	  !********* 1st Half Time Step *******

	  DO J=1,Y
	  DO I=1,X

	  IF (WET(I,J)==1) THEN

	  IF (I<X .and. J>1 .and. J<Y) THEN 
	  ADVv(I,J)=((V(I,J))*U(I,J)-(V(I,J-1))*U(I,J-1))/(Dy)
	  Vp(I,J)=(V(I,J)+V(I+1,J)+V(I,J-1)+V(I+1,J-1))/4.
	  ELSE IF (I==X .and. J>1.and. J<Y)     THEN
	  ADVv(I,J)=((V(I,J))*U(I,J)-(V(I,J-1))*U(I,J-1))/(Dy)
	  Vp(I,J)=(V(I,J)+V(I,J-1))/2.
	  ELSE IF (I<X .and. J==Y)    THEN
	  ADVv(I,J)=((V(I,J))*U(I,J)-(V(I,J-1))*U(I,J-1))/(Dy)
	  Vp(I,J)=(V(I,J)+V(I+1,J)+V(I,J-1)+V(I+1,J-1))/4.
	  ELSE IF (I==X .and. J==Y)     THEN 
	  ADVv(I,J)=((V(I,J))*U(I,J)-(V(I,J-1))*U(I,J-1))/(Dy)
	  Vp(I,J)=(V(I,J)+V(I,J-1))/2.
	  ELSE IF (I<X .and. J==1)     THEN 
	  ADVv(I,J)=((V(I,J))*U(I,J))/(Dy)
	  Vp(I,J)=(V(I,J)+V(I+1,J))/2.
	  ELSE
	  ADVv(I,J)=((V(I,J))*U(I,J))/(Dy)
	  Vp(I,J)=V(I,J)
	  END IF

	  Coriolis(I,J)=Omega*Vp(I,J)

	  END IF
	  
	  END DO
	  END DO

	  DO J=1,Y
	  DO I=1,X

	  IF (WET(I,J)==1) THEN

	  IF (I<X) THEN 
	  WINDx(I,J)=Twx(I,J)/(Pw*.5*(H(I,J)+H(I+1,J)))
	  Friction(I,J)=G*SQRT(U(I,J)**2.+Vp(I,J)**2.)/(Cchezy(I,J)**2*.5*(H(I,J)+H(I+1,J)))
	  ELSE 
	  WINDx(I,J)=Twx(I,J)/(Pw*.5*H(I,J))
	  Friction(I,J)=G*SQRT(U(I,J)**2.+Vp(I,J)**2.)/(Cchezy(I,J)**2*.5*H(I,J))
	  END IF

	  IF (J>1 .and. J<Y) THEN
	  DIFy(I,J)=Ah(I,J)*(U(I,J+1)+U(I,J-1)-2.*U(I,J))/(Dy**2.)
	  ELSE IF (J==1) THEN
	  DIFy(I,J)=Ah(I,J)*(U(I,J+1)-U(I,J))/(Dy**2.)
	  ELSE 
	  DIFy(I,J)=Ah(I,J)*(U(I,J-1)-U(I,J))/(Dy**2.)
	  END IF

	  END IF

	  END DO
	  END DO

	  !************ Coefficients Set Up ************

	  DO J=1,Y
	  DO I=1,X

	  IF (I>1) THEN
	  A1(I,J)=-(H(I,J)+H(I-1,J))/(2.*Dx)
	  ELSE
	  A1(I,J)=-H(I,J)/(2.*Dx)
	  END IF

	  B1(I,J)=2./Dt

	  IF (I<X) THEN
	  C1(I,J)=+(H(I,J)+H(I+1,J))/(2.*Dx)
	  ELSE
	  C1(I,J)=+H(I,J)/(2.*Dx)
	  END IF

	  D1(I,J)=(-U(I,J)/(Dx)-Ah(I,J)/(Dx**2.))
	  E1(I,J)=-G/Dx
	  F1(I,J)=+U(I,J)/(Dx)+2./Dt+Friction(I,J)+2.*Ah(I,J)/(Dx**2.)
	  G1(I,J)=+G/Dx
	  H1(I,J)=(-Ah(I,J)/(Dx**2.))

	  END DO
	  END DO

	  !************ Matrix Set Up ***************

	  DO J=1,Y
	  DO I=1,X

	  MX1(I,J)=2.*U(I,J)/Dt-ADVv(I,J)+Coriolis(I,J)+WINDx(I,J)+DIFy(I,J)

	  END DO 
	  END DO

	  DO J=1,Y
	  DO I=1,X

	  IF (J>1 .and. J<Y) THEN
	  CX1(I,J)=2.*Z(I,J)/Dt-((H(I,J)+H(I,J+1))*V(I,J)-(H(I,J)+H(I,J-1))*V(I,J-1))/(2.*Dy)
	  ELSE IF (J==1) THEN 
	  CX1(I,J)=2.*Z(I,J)/Dt-((H(I,J)+H(I,J+1))*V(I,J))/(2.*Dy)
	  ELSE 
	  CX1(I,J)=2.*Z(I,J)/Dt-(H(I,J)*V(I,J)-(H(I,J)+H(I,J-1))*V(I,J-1))/(2.*Dy)
	  END IF

	  END DO 
	  END DO

	  DO J=1,Y	 ! **** X Direction Sweep  START
	     
  	  DO I=1,2*X

	  IF (I==2) THEN
	  BX(I)=MX1(INT(I/2.),J)-D1(1,J)*UU(1,J)
	  ELSE IF (I==1) THEN
	  BX(I)=CX1(INT(I/2.+1/2.),J)-A1(1,J)*UU(1,J)
	  ELSE

	  IF (MOD(I,2)==0) THEN
	  IF (I==2*X) THEN
	  BX(I)=MX1(INT(I/2.),J)-G1(X,J)*ZZ(X,J)-H1(X,J)*UU(X,J)
	  ELSE
	  BX(I)=MX1(INT(I/2.),J)
	  END IF
	  ELSE 
	  BX(I)=CX1(INT(I/2.+1/2.),J)
	  END IF
	  END IF
	  END DO

	  AX=0
	  DO I=1,2*X
	  IF (MOD(I,2)/=0) THEN
	  AX(I,I)  =B1(INT(I/2.+1/2.),J)
	  AX(I,I+1)=C1(INT(I/2.+1/2.),J)
	  IF (I>1) THEN
	  AX(I,I-1)=A1(INT(I/2.+1/2.),J)
	  END IF
	  ELSE 
	  AX(I,I)  =F1(INT(I/2.),J)
	  IF (I<2*X-1) THEN
	  AX(I,I+1)=G1(INT(I/2.),J)
	  END IF
	  IF (I<2*X) THEN
	  AX(I,I+2)=H1(INT(I/2.),J)
	  END IF
	  IF (I>1) THEN
	  AX(I,I-1)=E1(INT(I/2.),J)
	  END IF
	  IF (I>2) THEN
	  AX(I,I-2)=D1(INT(I/2.),J)
	  END IF
	  END IF
	  END DO	  

      !********* PentaDiagonal *******

	  CALL PentaDiagonal (2*X,AX,BX,XX)

      CONTINUE     
	                           
	  !**********************************

	  DO M=1,2*X
	  IF (MOD(M,2)==0) THEN
	  UU(INT(M/2.),J)=XX(M)
	  ELSE 
	  ZZ(INT(M/2.+1/2.),J)=XX(M)
	  END IF
	  END DO

	  END DO   ! **** X Direction Sweep  START

	  U=UU
	  Z=ZZ

	  !************* Oil Dynamic *************

	  ADVv=0
	  DIFy=0
	  CXX=0
	  SX=0

	  DO I=1,X
	  DO J=1,Y

	  LANDA(I,J)=(G*S(I,J)**2)*(Pw-Po)*Po/(Pw*Cfilm)

	  IF (WET(I,J)==1) THEN

	  IF (J>1 .and. J<Y) THEN 
	  ADVv(I,J)=((S(I,J)+S(I,J+1))*(V(I,J)+0.03*Vw(I,J)*SIN(Pi*WINDteta(I,J)))-(S(I,J)+S(I,J-1))*(V(I,J-1)+0.03*Vw(I,J-1)*SIN(Pi*WINDteta(I,J-1))))/(2.*Dy)
	  DIFy(I,J)=LANDA(I,J)*(S(I,J+1)+S(I,J-1)-2.*S(I,J))/(Dy**2.)
	  ELSE IF (J==1)     THEN
	  ADVv(I,J)=((S(I,J)+S(I,J+1))*(V(I,J)+0.03*Vw(I,J)*SIN(Pi*WINDteta(I,J))))/(2.*Dy)
	  DIFy(I,J)=LANDA(I,J)*(S(I,J+1)-2.*S(I,J))/(Dy**2.)
	  ELSE IF (J==Y)     THEN
	  ADVv(I,J)=(S(I,J)*(V(I,J)+0.03*Vw(I,J)*SIN(Pi*WINDteta(I,J)))-(S(I,J)+S(I,J-1))*(V(I,J-1)+0.03*Vw(I,J-1)*SIN(Pi*WINDteta(I,J-1))))/(2.*Dy)
	  DIFy(I,J)=LANDA(I,J)*(S(I,J-1)-2.*S(I,J))/(Dy**2.)
	  END IF

	  END IF

	  END DO
	  END DO

	  !********** Coefficients Set Up ************

	  DO I=1,X
	  DO J=1,Y

	  IF (I>1) THEN
	  A1(I,J)=-(U(I-1,J)+0.03*Vw(I-1,J)*COS(Pi*WINDteta(I-1,J)))/(2.*Dx)-LANDA(I,J)/(Dx**2.)
	  ELSE
	  A1(I,J)=-LANDA(I,J)/(Dx**2.)
	  END IF

	  IF (I>1) THEN
	  B1(I,J)=2./Dt+(U(I,J)+0.03*Vw(I,J)*COS(Pi*WINDteta(I,J)))/(2.*Dx)-(U(I-1,J)+0.03*Vw(I-1,J)*COS(Pi*WINDteta(I-1,J)))/(2.*Dx)+2*LANDA(I,J)/(Dx**2.)
	  ELSE
	  B1(I,J)=2./Dt+(U(I,J)+0.03*Vw(I,J)*COS(Pi*WINDteta(I,J)))/(2.*Dx)+2*LANDA(I,J)/(Dx**2.)
	  END IF

	  C1(I,J)=(U(I,J)+0.03*Vw(I,J)*COS(Pi*WINDteta(I,J)))/(2.*Dx)-LANDA(I,J)/(Dx**2.)

	  END DO
	  END DO

	  !************ Matrix Set Up ***************

!	  DO I=1,X
!	  DO J=1,Y

!	  EVAP(I,J)=0

!	  END DO
!	  END DO

	  DO I=1,X
	  DO J=1,Y

	  CX1(I,J)=2.*S(I,J)/Dt+DIFy(I,J)-ADVv(I,J)+2.*EVAP(I,J)/Dt

	  END DO
	  END DO

	  DO J=1,Y
	  	    
  	  DO I=1,X

	  IF (I==1) THEN
	  DDD(I)=CX1(I,J)-A1(1,J)*SSin
	  ELSE IF (I==X) THEN
	  DDD(I)=CX1(I,J)-C1(X,J)*SSout
	  ELSE
	  DDD(I)=CX1(I,J)
	  END IF

	  AAA(I)=A1(I,J)
	  BBB(I)=B1(I,J)
	  CCC(I)=C1(I,J)

	  END DO

	  !************* Tridiagonal Matrix **************

	  CALL Thomas (X,AAA,BBB,CCC,DDD,CXX)

      CONTINUE  
	  
	  !***********************************************
	  
	  DO M=1,X
	  
	  IF (CXX(M)<0) CXX(M)=0

	  SX(M,J)=CXX(M) 
	  
	  END DO

	  END DO
	  
	  S=SX	  

	  !******** Initial & Boundary Condition ***********

		DO J=1,Y
	    DO I=1,X

		Z(X,J)=HWave/2.*COS(-2*Pi*P*dT/Period)
		H(I,J)=D(I,J)+Z(I,J)

		IF (H(I,J)<0.00) H(I,J)=0
        IF (H(I,J)==0)   WET(I,J)=0 
	    IF (H(I,J)>0)    WET(I,J)=1
!  	    IF (WET(I,J)==0) H(I,J)=0
		IF (WET(I,J)==0) S(I,J)=0 
		IF (WET(I,J)==0) U(I,J)=0
!		IF (WET(I-1,J)==0.and.WET(I,J)==1) U(I-1,J)=U(I,J)  
		IF (WET(I,J)==0) V(I,J)=0
		IF (H(I,J)==0)   Ah(I,J)=0

		IF (WET(I,J)==1) THEN

		Cchezy(I,J)=-18*LOG10(K/(H(I,J)*12))
	    Ah(I,J)=0.0376 !Ce*H(I,J)*SQRT(G*(U(I,J)**2+V(I,J)**2))/Cchezy(I,J)

		END IF

		END DO
		END DO

		advv=0
		vp=0
		up=0
		WINDx=0
		WINDy=0
		Coriolis=0
		DIFx=0
		DIFy=0
		Friction=0


      !********* 2nd Half Time Step *******

	  DO I=1,X
	  DO J=1,Y

	  IF (WET(I,J)==1) THEN

	  IF (I>1 .and. I<X .and. J<Y)        THEN 
	  ADVv(I,J)=((U(I,J))*V(I,J)-(U(I-1,J))*V(I-1,J))/(Dx)
	  Up(I,J)=(U(I,J)+U(I,J+1)+U(I-1,J)+U(I-1,J+1))/4.
	  ELSE IF (I>1 .and. I<X .and. J==Y)  THEN
	  ADVv(I,J)=((U(I,J))*V(I,J)-(U(I-1,J))*V(I-1,J))/(Dx)
	  Up(I,J)=(U(I,J)+U(I-1,J))/2.
	  ELSE IF (I==X .and. J<Y)    THEN
	  ADVv(I,J)=((U(I,J))*V(I,J)-(U(I-1,J))*V(I-1,J))/(Dx)
	  Up(I,J)=(U(I,J)+U(I,J+1)+U(I-1,J)+U(I-1,J+1))/4.
	  ELSE IF (I==X .and. J==Y)     THEN 
	  ADVv(I,J)=((U(I,J))*V(I,J)-(U(I-1,J))*V(I-1,J))/(Dx)
	  Up(I,J)=(U(I,J)+U(I-1,J))/2.
	  ELSE IF (I==1 .and. J<Y)     THEN
	  ADVv(I,J)=((U(I,J))*V(I,J))/(Dx)
	  Up(I,J)=(U(I,J)+U(I,J+1))/2.
	  ELSE
	  ADVv(I,J)=((U(I,J))*V(I,J))/(Dx)
	  Up(I,J)=U(I,J)
	  END IF

	  Coriolis(I,J)=-Omega*Up(I,J)

	  END IF

	  END DO
	  END DO

	  DO I=1,X
	  DO J=1,Y

	  IF (WET(I,J)==1) THEN

	  IF (J<Y) THEN 
	  WINDy(I,J)=Twy(I,J)/(Pw*.5*(H(I,J)+H(I,J+1)))
	  Friction(I,J)=G*SQRT(V(I,J)**2.+Up(I,J)**2.)/(Cchezy(I,J)**2*.5*(H(I,J)+H(I,J+1)))
	  ELSE 
	  WINDy(I,J)=Twy(I,J)/(Pw*.5*H(I,J))
	  Friction(I,J)=G*SQRT(V(I,J)**2.+Up(I,J)**2.)/(Cchezy(I,J)**2*.5*H(I,J))
	  END IF

	  IF (I>1 .and. I<X) THEN
	  DIFx(I,J)=Ah(I,J)*(V(I+1,J)+V(I-1,J)-2.*V(I,J))/(Dx**2.)
	  ELSE IF (I==1) THEN
	  DIFx(I,J)=Ah(I,J)*(V(I+1,J)-V(I,J))/(Dx**2.)
	  ELSE 
	  DIFx(I,J)=Ah(I,J)*(V(I-1,J)-V(I,J))/(Dx**2.)
	  END IF

	  END IF

	  END DO
	  END DO

	  !************ Coefficients Set Up ************

	  DO I=1,X
	  DO J=1,Y

	  IF (J>1) THEN
	  A1(I,J)=-(H(I,J)+H(I,J-1))/(2.*Dy)
	  ELSE
	  A1(I,J)=-H(I,J)/(2.*Dy)
	  END IF

	  B1(I,J)=2./Dt

	  IF (J<Y) THEN
	  C1(I,J)=+(H(I,J)+H(I,J+1))/(2.*Dy)
	  ELSE
	  C1(I,J)=+H(I,J)/(2.*Dy)
	  END IF

	  D1(I,J)=(-V(I,J)/(Dy)-Ah(I,J)/(Dy**2.))
	  E1(I,J)=-G/(Dy)
	  F1(I,J)=+V(I,J)/(Dy)+2./Dt+Friction(I,J)+2.*Ah(I,J)/(Dy**2.)
	  G1(I,J)=+G/(Dy)
	  H1(I,J)=(-Ah(I,J)/(Dy**2.))

	  END DO
	  END DO

	  !************ Matrix Set Up ***************

	  DO I=1,X
	  DO J=1,Y

	  MX2(I,J)=2.*V(I,J)/Dt-ADVv(I,J)+Coriolis(I,J)+WINDy(I,J)+DIFx(I,J)

	  IF (I>1 .and. I<X) THEN
	  CX1(I,J)=2.*Z(I,J)/Dt-((H(I,J)+H(I+1,J))*U(I,J)-(H(I,J)+H(I-1,J))*U(I-1,J))/(2.*Dx)
	  ELSE IF (I==1) THEN 
	  CX1(I,J)=2.*Z(I,J)/Dt-((H(I,J)+H(I+1,J))*U(I,J))/(2.*Dx)
	  ELSE 
	  CX1(I,J)=2.*Z(I,J)/Dt-(H(I,J)*U(I,J)-(H(I,J)+H(I-1,J))*U(I-1,J))/(2.*Dx)
	  END IF  

	  END DO 
	  END DO

	  DO I=1,X   ! **** Y Direction Sweep  START

  	  DO J=1,2*Y

	  IF (J==2) THEN
	  BX(J)=MX2(I,INT(J/2.))-D1(I,1)*VV(I,1)
	  ELSE IF (J==1) THEN
	  BX(J)=CX1(I,INT(J/2.+1/2.))-A1(I,1)*VV(I,1)
	  ELSE
	  
	  IF (MOD(J,2)==0) THEN
	  IF (J==2*Y) THEN
	  BX(J)=MX2(I,INT(J/2.))-G1(I,Y)*ZZ(I,Y)-H1(I,Y)*VV(I,Y)
	  ELSE
	  BX(J)=MX2(I,INT(J/2.))
	  END IF
	  ELSE 
	  BX(J)=CX1(I,INT(J/2.+1/2.))
	  END IF
	  END IF
	  END DO

	  AX=0
	  DO J=1,2*Y
	  IF (MOD(J,2)/=0) THEN
	  AX(J,J)  =B1(I,INT(J/2.+1/2.))
	  AX(J,J+1)=C1(I,INT(J/2.+1/2.))
	  IF (J>1) THEN
	  AX(J,J-1)=A1(I,INT(J/2.+1/2.))
	  END IF
	  ELSE 
	  AX(J,J)  =F1(I,INT(J/2.))
	  IF (J<2*Y-1) THEN
	  AX(J,J+1)=G1(I,INT(J/2.))
	  END IF
	  IF (J<2*Y) THEN
	  AX(J,J+2)=H1(I,INT(J/2.))
	  END IF
	  IF (J>1) THEN
	  AX(J,J-1)=E1(I,INT(J/2.))
	  END IF
	  IF (J>2) THEN
	  AX(J,J-2)=D1(I,INT(J/2.))
	  END IF
	  END IF
	  END DO

      !********* PentaDiagonal *******

	  CALL PentaDiagonal (2*Y,AX,BX,XX)

	  CONTINUE                              


	  !**********************************

	  DO N=1,2*Y

	  IF (MOD(N,2)==0) THEN
	  VV(I,INT(N/2.))=XX(N)
	  ELSE 
	  ZZ(I,INT(N/2.+1/2.))=XX(N)
	  END IF
	  END DO

	  END DO   ! **** Y Direction Sweep  END

	  V=VV
	  Z=ZZ

	  !************* Oil Dynamic *************

	  ADVv=0
	  DIFy=0
	  CXX=0
	  SX=0

	  DO I=1,X
	  DO J=1,Y

	  LANDA(I,J)=(G*S(I,J)**2)*(Pw-Po)*Po/(Pw*Cfilm)

	  IF (WET(I,J)==1) THEN

	  IF (I>1 .and. I<X) THEN 
	  ADVv(I,J)=((S(I,J)+S(I+1,J))*(U(I,J)+0.03*Vw(I,J)*COS(Pi*WINDteta(I,J)))-(S(I,J)+S(I-1,J))*(U(I-1,J)+0.03*Vw(I-1,J)*COS(Pi*WINDteta(I-1,J))))/(2.*Dx)
	  DIFy(I,J)=LANDA(I,J)*(S(I+1,J)+S(I-1,J)-2.*S(I,J))/(Dx**2.)
	  ELSE IF (I==1)     THEN
	  ADVv(I,J)=((S(I,J)+S(I+1,J))*(U(I,J)+0.03*Vw(I,J)*COS(Pi*WINDteta(I,J))))/(2.*Dx)
	  DIFy(I,J)=LANDA(I,J)*(S(I+1,J)-2.*S(I,J))/(Dx**2.)
	  ELSE IF (I==X)     THEN
	  ADVv(I,J)=(S(I,J)*(U(I,J)+0.03*Vw(I,J)*COS(Pi*WINDteta(I,J)))-(S(I,J)+S(I-1,J))*(U(I-1,J)+0.03*Vw(I-1,J)*COS(Pi*WINDteta(I-1,J))))/(2.*Dx)
	  DIFy(I,J)=LANDA(I,J)*(S(I-1,J)-2.*S(I,J))/(Dx**2.)
	  END IF

	  END IF

	  END DO
	  END DO

	  !********** Coefficients Set Up ************

	  DO I=1,X
	  DO J=1,Y

	  IF (J>1) THEN
	  A1(I,J)=-(V(I,J-1)+0.03*Vw(I,J-1)*SIN(Pi*WINDteta(I,J-1)))/(2.*Dy)-LANDA(I,J)/(Dy**2.)
	  ELSE
	  A1(I,J)=-LANDA(I,J)/(Dy**2.)
	  END IF

	  IF (J>1) THEN
	  B1(I,J)=2./Dt+(V(I,J)+0.03*Vw(I,J)*SIN(Pi*WINDteta(I,J)))/(2.*Dy)-(V(I,J-1)+0.03*Vw(I,J-1)*SIN(Pi*WINDteta(I,J-1)))/(2.*Dy)+2*LANDA(I,J)/(Dy**2.)
	  ELSE
	  B1(I,J)=2./Dt+(V(I,J)+0.03*Vw(I,J)*SIN(Pi*WINDteta(I,J)))/(2.*Dy)+2*LANDA(I,J)/(Dy**2.)
	  END IF

	  C1(I,J)=(V(I,J)+0.03*Vw(I,J)*SIN(Pi*WINDteta(I,J)))/(2.*Dy)-LANDA(I,J)/(Dy**2.)

	  END DO
	  END DO

	  !************ Matrix Set Up ***************

!	  DO I=1,X
!	  DO J=1,Y

!	  EVAP(I,J)=0

!	  END DO
!	  END DO

	  DO I=1,X
	  DO J=1,Y

	  CX1(I,J)=2.*S(I,J)/Dt+DIFy(I,J)-ADVv(I,J)+2.*EVAP(I,J)/Dt

	  END DO
	  END DO

	  DO I=1,X
	  	    
  	  DO J=1,Y

	  IF (J==1) THEN
	  DDD(J)=CX1(I,J)-A1(I,1)*SSin
	  ELSE IF (J==Y) THEN
	  DDD(J)=CX1(I,J)-C1(I,Y)*SSout
	  ELSE
	  DDD(J)=CX1(I,J)
	  END IF

	  AAA(J)=A1(I,J)
	  BBB(J)=B1(I,J)
	  CCC(J)=C1(I,J)

	  END DO

	  !************* Tridiagonal Matrix **************

	  CALL Thomas (Y,AAA,BBB,CCC,DDD,CXX)

      CONTINUE  
	  
	  !***********************************************
	  
	  DO M=1,Y

	  IF (CXX(M)<0) CXX(M)=0
	  
	  SX(I,M)=CXX(M) 
	  
	  END DO

	  END DO
	  
	  S=SX

	  !********* Output *******
	    	  
	  PRINT *, P
	  END DO

	  DO I=1,X
	  
	  WRITE (60,'(50(I1))'), INT(WET(I,1:Y))

	  END DO

	  !********* SUBROUTINE Penta-Diagonal *******

	  ! A Computational Algorithm For Solving Nearly Penta-Diagonal Matrix 

	  CONTAINS

	  SUBROUTINE PentaDiagonal (Arrays,A,BX,XX)

	  INTEGER M,R,Arrays
      REAL  S,T
      REAL  BX(MAX(2*X,2*Y)),A(MAX(2*X,2*Y),MAX(2*X,2*Y))
	  REAL  C(MAX(2*X,2*Y)),E(MAX(2*X,2*Y)),F(MAX(2*X,2*Y))
	  REAL  H(MAX(2*X,2*Y)),V(MAX(2*X,2*Y)),Z(MAX(2*X,2*Y)),SIGMA(MAX(2*X,2*Y))
	  REAL, INTENT(OUT) :: XX(MAX(2*X,2*Y))

	  S=0.
	  T=0.
	  C(1)=A(1,1)
	  F(1)=0.
	  E(1)=A(1,2)
	  V(1)=T
	  H(1)=S/C(1)
	  F(2)=A(2,1)/C(1)
	  C(2)=A(2,2)-F(2)*E(1)
	  V(2)=-F(2)*V(1)
	  H(2)=-E(1)*H(1)/C(2)
	  E(2)=A(2,3)-F(2)*A(1,3)

	  DO M=3,Arrays-1
	  F(M)=A(M,M-1)/C(M-1)-A(M,M-2)*E(M-2)/(C(M-2)*C(M-1))
	  IF (M<Arrays-1) THEN
	  E(M)=A(M,M+1)-F(M)*A(M-1,M+1)
	  END IF
	  C(M)=A(M,M)-F(M)*E(M-1)-A(M,M-2)*A(M-2,M)/C(M-2)
	  END DO

	  DO M=3,Arrays-1
	  IF (M==Arrays-2) THEN
	  V(M)=A(M,M+2)-A(M,M-2)*V(M-2)/C(M-2)-F(M)*V(M-1)
	  H(M)=A(M+2,M)/C(M)-(A(M-2,M)*H(M-2)/C(M)+E(M-1)*H(M-1)/C(M))
	  ELSE IF (M==Arrays-1) THEN
	  V(M)=A(M,M+1)-A(M,M-2)*V(M-2)/C(M-2)-F(M)*V(M-1)
	  H(M)=A(M+1,M)/C(M)-(A(M-2,M)*H(M-2)/C(M)+E(M-1)*H(M-1)/C(M))
	  ELSE 
	  V(M)=A(M,M-2)*V(M-2)/C(M-2)-F(M)*V(M-1)
	  H(M)=-(A(M-2,M)*H(M-2)/C(M)+E(M-1)*H(M-1)/C(M))
	  END IF
	  END DO

	  SIGMA=H*V
	  C(Arrays)=A(Arrays,Arrays)-(SUM(SIGMA)-H(Arrays)*V(Arrays))

	  Z(1)=BX(1)
	  Z(2)=BX(2)-F(2)*Z(1)

	  DO M=3,Arrays-1
	  Z(M)=BX(M)-F(M)*Z(M-1)-A(M,M-2)*Z(M-2)/C(M-2)
	  END DO

	  SIGMA=H*Z
	  Z(Arrays)=BX(Arrays)-(SUM(SIGMA)-H(Arrays)*Z(Arrays))

	  XX(Arrays)=Z(Arrays)/C(Arrays)
	  XX(Arrays-1)=(Z(Arrays-1)-V(Arrays-1)*XX(Arrays))/C(Arrays-1)
	  XX(Arrays-2)=(Z(Arrays-2)-E(Arrays-2)*XX(Arrays-1)-V(Arrays-2)*XX(Arrays))/C(Arrays-2)

	  DO M=Arrays-3,1,-1
	  XX(M)=(Z(M)-E(M)*XX(M+1)-A(M,M+2)*XX(M+2)-V(M)*XX(Arrays))/C(M)
	  END DO

	  RETURN

	  END SUBROUTINE PentaDiagonal

	  !********* SUBROUTINE Tri-Diagonal *******

	  ! A Computational Algorithm For Solving Tri-Diagonal Matrix

	  SUBROUTINE Thomas (Arrays,AAA,BBB,CCC,DDD,CXX)

	  INTEGER M,R,Arrays
	  REAL  AAA(MAX(X,Y)),BBB(MAX(X,Y)),CCC(MAX(X,Y)),DDD(MAX(X,Y))
      REAL, INTENT(OUT) :: CXX(MAX(X,Y))

	  IF (BBB(1)==0.0) STOP 666
	  ZX=BBB(1)
      C(1)=DDD(1)/ZX

      DO  M=2,Arrays
         CX(M)=CCC(M-1)/ZX
         ZX=BBB(M)-AAA(M)*CX(M)
		 IF (ZX==0.0) EXIT
         C(M)=(DDD(M)-AAA(M)*C(M-1))/ZX
	  END DO

	  DO  M=Arrays-1,1,-1
          CXX(M)=C(M)-CX(M+1)*C(M+1)
      END DO	  

	  RETURN

	  END SUBROUTINE Thomas

     !********* SUBROUTINE Guass Seidel*******

	  SUBROUTINE GuassSeidel (Arrays,A,BX,XX)

	  INTEGER M,R,Iteration,Arrays
      REAL  DiffMax,Diff,Temp,SUM,Test
      REAL  BX(MAX(2*X,2*Y)),A(MAX(2*X,2*Y),MAX(2*X,2*Y))
	  REAL, INTENT(OUT) :: XX(MAX(2*X,2*Y))

	  TEST=0.0001
	  XX=0
	  Iteration=0
	  DO
      DiffMax=0
      Iteration=Iteration+1
	  IF (Iteration > 100) STOP 1005
      DO M=1,Arrays
	  Temp=XX(M)
	  SUM=0
	  DO R=1,Arrays
	  IF (R /= M) THEN
	  SUM=SUM+XX(R)*AX(M,R)
	  END IF
	  END DO
      XX(M)=(BX(M)-SUM)/AX(M,M)
	  Diff=ABS(Temp-XX(M))
	  IF (Diff > DiffMax) DiffMax=Diff
      END DO
	  IF (DiffMax < Test) EXIT
	  END DO

	  RETURN

	  END SUBROUTINE GuassSeidel

	  !********* END PROGRAM  *******

	  END PROGRAM Hydrodynamic

