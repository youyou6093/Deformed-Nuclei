
C------------------------------------------------------------------------
c     ********************** RACAH ROUTINE **********************
c
C 	CALCULATE VARIOUS RACAH COEFFICIENTS (i.e., 3J,6J,9J)
C           EVERYTHING IN DOUBLE PRECISION (MAY-16-1990)
c-----------------------------------------------------------------------
        REAL*8 FUNCTION TJ(A1,A2,A3,B1,B2,B3)
        IMPLICIT REAL*8 (A-H,O-Z)
        IF(TR(A1,A2,A3).LT.0.)GO TO 8
        IF(B1+B2+B3.NE.0.)GO TO 8
        IF(PAR(2.*(A1+A2+A3)).EQ.-1.)GO TO 8
        N1=ABS(MIN(A3-A2+B1,A3-A1-B2,0.0D0))
        N2=MIN(A1+A2-A3,A1-B1,A2+B2)
        IF(N2.LT.N1)GO TO 8
        N1=N1+1
        N2=N2+1
        S=0.
        DO 10 I=N1,N2
        Y=I-1
        T=BINOM(A3+A1-A2+Y,A1-B1)*BINOM(A3+A2-A1+Y,A2+B2)/BINOM(2.*Y,Y)
        T=T*BINOM(A1-B1,Y)*BINOM(A2+B2,Y)/BINOM(2.*A3+2.*Y,2.*Y)
10      S=S+PAR(Y)*T*BINOM(A1+A2-A3,Y)*BINOM(2.*A3+2.*Y,A3+A1-A2+Y)
        T=BINOM(A1+A2+A3,2.*A3)/(BINOM(A1+A2+A3,2.*A1)*BINOM(A1+A2+A3,2.
     1   *A2)*(A1+A2+A3+1.))
        T=SQRT(T)/SQRT(BINOM(2.*A1,A1+B1)*BINOM(2.*A2,A2+B2)*BINOM(2.*A3
     1   ,A3+B3))
        TJ=PAR(A1-A2-B3)*T*S
        RETURN
8       TJ=0.
        RETURN
        END

        FUNCTION CG(A1,B1,A2,B2,A3,B3)
        IMPLICIT REAL*8 (A-H,O-Z)
        IF(TR(A1,A2,A3).LT.0.)GO TO 8
        IF(B1+B2.NE.B3)GO TO 8
        IF(PAR(2.*(A1+A2+A3)).EQ.-1.)GO TO 8
        N1=ABS(MIN(A3-A2+B1,A3-A1-B2,0.0D0))
        N2=MIN(A1+A2-A3,A1-B1,A2+B2)
        IF(N2.LT.N1)GO TO 8
        N1=N1+1
        N2=N2+1
        S=0.
        DO 10 I=N1,N2
        Y=I-1
        T=BINOM(A3+A1-A2+Y,A1-B1)*BINOM(A3+A2-A1+Y,A2+B2)/BINOM(2.*Y,Y)
        T=T*BINOM(A1-B1,Y)*BINOM(A2+B2,Y)/BINOM(2.*A3+2.*Y,2.*Y)
10      S=S+PAR(Y)*T*BINOM(A1+A2-A3,Y)*BINOM(2.*A3+2.*Y,A3+A1-A2+Y)
        T=BINOM(A1+A2+A3,2.*A3)/(BINOM(A1+A2+A3,2.*A1)*BINOM(A1+A2+A3,2.
     1   *A2)*(A1+A2+A3+1.))
        T=SQRT(T)/SQRT(BINOM(2.*A1,A1+B1)*BINOM(2.*A2,A2+B2)*BINOM(2.*A3
     1   ,A3+B3))
        CG=SQRT(2.*A3+1.)*T*S
        RETURN
8       CG=0.
        RETURN
        END

      FUNCTION ODD(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      P(A) = PAR(A)
      ODD = (1.0 - P(X))/2.0
      RETURN
      END

      FUNCTION EVEN(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      P(A) = PAR(A)
      EVEN = (1.0+P(X))/2.0
      RETURN
      END

      FUNCTION DELTA(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X) 10,20,10
20    DELTA = 1.0
      RETURN
10    DELTA = 0.0
      RETURN
      END

        FUNCTION SJ(A1,A2,A3,B1,B2,B3)
        IMPLICIT REAL*8 (A-H,O-Z)
        D(A,B,C)=1./SQRT(BINOM(2.*A,A+B-C)*BINOM(A+B+C,2.*A)*(A+B+C+1.))
        IF(TR(A1,A2,A3).LT.0..OR.TR(B1,B2,A3).LT.0..OR.TR(A1,B2,B3).LT.0
     1   ..OR.TR(B1,A2,B3).LT.0.)GO TO 9
        N1=ABS(MIN(-A1-B1+A3+B3,-A2-B2+A3+B3,0.0D0))
        N2=MIN(A1+A2+B1+B2+1.,A1+A2-A3,B1+B2-A3,A1+B2-B3,B1+A2-B3)
        IF(N2.LT.N1)GO TO 9
        N1=N1+1
        N2=N2+1
        S=0.
        DO 10 I=N1,N2
        Y=I-1
        T=(BINOM(A1+A2-A3,Y)*BINOM(B1+B2-A3,Y)/BINOM(A1+B1+A2+B2+2.*
     1   Y+1.,3.*Y))
        T=T*(BINOM(A1+B2-B3,Y)*BINOM(B1+A2-B3,Y)/BINOM(3.*Y,Y))
         T=T*(BINOM(A2-B1+B3+Y,A1+A2-A3)*BINOM(B1-A2+B3+Y,B1+B2-A3)/BINO
     1   M(2.*Y,Y))
        T=T*BINOM(2.*B3+2.*Y,A2-B1+B3+Y)*BINOM(A1+B1+A2+B2-2.*B3,A1+B2-
     1   B3)*BINOM(A1+B1+A2+B2+2.*Y,2.*B3+2.*Y)*(A1+B1+A2+B2+2.*Y+1.)
10      S=S+PAR(Y)*T
        T=D(A1,A2,A3)*S
        T=T*D(A1,B2,B3)
        T=T*D(B1,A2,B3)
        SJ=PAR(A1+A2+B1+B2)*T*D(B1,B2,A3)
        RETURN
9       SJ=0.
        RETURN
        END

      FUNCTION DFACN(A)
      IMPLICIT REAL*8 (A-H,O-Z)
      DFACN = DFAC(A)
      RETURN
      END

      FUNCTION FAC(X)
      IMPLICIT REAL*8 (A-H,O-Z)
      IF(X) 3,4,4
3     FAC = 0.0
      RETURN
4     S = 1.0
      DO 1 I = 1,1000
      Y = I
      IF(X-Y) 2,1,1
1     S = S*Y
2      FAC = S
      RETURN
      END

        FUNCTION BINOM(A,B)
        IMPLICIT REAL*8 (A-H,O-Z)
        IF(A.EQ.B)GO TO 12
        IF(B.EQ.0.)GO TO 12
        IF(A-B.GT.B)GO TO 9
        M=A-B
        GO TO 11
9       M=B
11      BINOM=1.
        X=A-FLOAT(M)
        DO 10 I=1,M
        AI=I
10      BINOM=BINOM*(X+AI)/AI
        RETURN
12      BINOM=1.
        RETURN
        END

        FUNCTION DFAC(A)
        IMPLICIT REAL*8 (A-H,O-Z)
        IF(A)6,7,8
6       DFAC=0.
        RETURN
7       DFAC=1.
        RETURN
8       N=A
        DFAC=1.
        DO 9 I=1,N,2
        AI=I
9       DFAC=AI*DFAC
        RETURN
        END

	FUNCTION par(x)
        IMPLICIT REAL*8 (A-H,O-Z)
	eps = 0.01
	xx = ABS( x )
	int1 = xx / 2.0 + eps
	int1 = 2 * int1
	int2 = xx + eps
	IF ( int1 .EQ. int2 ) THEN
	  par = + 1.0
	ELSE 
	  par = - 1.0
	END IF
	RETURN
	END

	FUNCTION TR(J1,J2,J3)
        IMPLICIT REAL*8 (A-H,O-Z)
	REAL*8  J1, J2, J3, S1, S2, S3
	S1 = J1 + J2 - J3 + 0.1
	S2 = J2 + J3 - J1 + 0.1
	S3 = J3 + J1 - J2 + 0.1
	IF((S1.LT.0.0).OR.(S2.LT.0.0).OR.(S3.LT.0.0))  GOTO 10
	TR = + 1.0
	RETURN
   10	TR = - 1.0
	RETURN
	END

	FUNCTION CLEBSCH(J1,J2,J3)
        IMPLICIT REAL*8 (A-H,O-Z)
CMRZ	CALCULATES THE VALUE OF SQRT((2J1+1)(2J2+1))*TJ(J1,0.5,J2,-0.5,J3,0.)

	REAL*8 J1,J2,J3
	IF (TR(J1,J2,J3))  5,5,1
    1   IF (PAR(J1+J2+J3))  2,2,3
    2   Z = J3 + 1.0
	GO TO 4
    3   Z = J3
    4   Y=(FAC(J1+J2-J3)*FAC(J1+J3-J2)*FAC(J2+J3-J1))/FAC(J1+J2+J3+1.)
	Y=SQRT(Y)
	Y=2.0*Y*FAC((Z+J1+J2)/2.0)/
     #  (FAC((J1+J2-Z)/2.)*FAC((J1+Z-J2-1.)/2.)*FAC((J2+Z-J1-1.)/2.))
	CLEBSCH = Y*PAR((J1+J2+Z-2.0)/2.0)
	RETURN
    5   CLEBSCH = 0.0
	RETURN  
	END

	FUNCTION wsj(a1,a2,b2,b1,a3,b3)
        IMPLICIT REAL*8 (A-H,O-Z)
	wsj = par(a1+a2+b1+b2) * sj(a1,a2,a3,b1,b2,b3)
	RETURN
	END

	FUNCTION ssj(a1,a2,a0,b2,b1)
        IMPLICIT REAL*8 (A-H,O-Z)
	     IF ( (b1 .GT. a1) .AND. (b2 .GT. a2) ) THEN
	  GOTO 10
	ELSE IF ( (b1 .GT. a1) .AND. (b2 .LT. a2) ) THEN
	  GOTO 20
	ELSE IF ( (b1 .LT. a1) .AND. (b2 .GT. a2) ) THEN
	  GOTO 30
	ELSE IF ( (b1 .LT. a1) .AND. (b2 .LT. a2) ) THEN
	  GOTO 40
	END IF

   10	ssj = ( a0 + b1 + b2 + 1.0 ) * ( b1 + b2 - a0 )
	ssj = ssj / ( (2.0*b1) * (2.0*b1+1.0) * (2.0*b2) * (2.0*b2+1.0) )
	ssj = par( a0 + b1 + b2 ) * SQRT( ssj )
	RETURN

   20	ssj = ( a0 + a2 - a1 ) * ( a0 + a1 - a2 + 1.0 )
	ssj = ssj / ( (2.0*a1+1.0)*(2.0*a1+2.0)*(2.0*a2)*(2.0*a2+1.0) )
	ssj = par( a0 + a1 + a2 ) * SQRT( ssj )
	RETURN

   30	ssj = ( a0 + b2 - b1 ) * ( a0 + b1 - b2 + 1.0 )
	ssj = ssj / ( (2.0*b1+1.0)*(2.0*b1+2.0)*(2.0*b2)*(2.0*b2+1.0) )
	ssj = par( a0 + b1 + b2 ) * SQRT( ssj )
	RETURN

   40	ssj = ( a0 + a1 + a2 + 1.0 ) * ( a1 + a2 - a0 )
	ssj = ssj / ( (2.0*a1) * (2.0*a1+1.0) * (2.0*a2) * (2.0*a2+1.0) )
	ssj = par( a0 + a1 + a2 ) * SQRT( ssj )
	RETURN
	END	
C
C
C  TO CALC 9J		1	2	3
C			4	5	6
C			7	8	9
C
C  SUM OVER RACAH'S
C
C
	FUNCTION coef9(a1,a2,a3,a4,a5,a6,a7,a8,a9)
        IMPLICIT REAL*8 (A-H,O-Z)
	coef9 = 0.0

	bmin = MAX( ABS(a1-a9) , ABS(a4-a8) , ABS(a2-a6) )
	bmax = MIN( a1+a9 , a4+a8 , a2+a6 )
	dif = bmax - bmin + 0.01
	IF ( dif .LT. 0.0 )  RETURN
	nj = dif

	DO 100  kj = 0, nj
	b = bmin + kj
  100	coef9 = coef9 + (2.0*b+1.0)*wsj(a1,a9,a4,a8,b,a7)
     &		*wsj(a2,a6,a8,a4,b,a5)*wsj(a1,a9,a2,a6,b,a3)

	RETURN
	END
c-----------------------------------------------------------------------










