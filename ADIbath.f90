      PROGRAM Bathymetry
	  !    Bathymetry 
      IMPLICIT NONE 

      REAL  S,L,Pi,dX,dY,B,X,Y
      INTEGER I,J,P,K
     
	  OPEN (10, FILE="INbath.TXT",status="unknown")
	  OPEN (20, FILE="Bathymetry.DAT",status="unknown")
	  OPEN (30, FILE="Bathymetry.TXT",status="unknown")
   	  READ (10,*) P,K,S,L,dX,dY
	  
	  Pi=3.1415926535897932
	   
      WRITE(20,*), 'Title="Hydrodynamic"'
      WRITE(20,*), 'variables="i","j","b"'
      WRITE(20,*), 'zone t="zone number-' ,1,'",i=',P,',j=',K,''
      WRITE(20,*), 'f=point'	  
	  
	  DO J=1,K
	  DO I=1,P

	  IF (I<=37) THEN 
	  B=+0.225
      ELSE IF (I>99) THEN
	  B=-((I-100)*0.3/45) 	  
	  ELSE 
	  B=-S*((I-37)*dX-0.875*SIN(Pi*(I-37)*dX/(3.5))*SIN(2.*Pi*(J)*dY/(L)))+0.225
	  END IF

	  WRITE(20,*),146-I,J,B
	  WRITE(30,*),146-I,J,B

      END DO
	  END DO

	  END PROGRAM Bathymetry

