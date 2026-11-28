C     Number editing check for prustpro: reads lines "kind hex" or
C     "S 11-character-field" and prints the Fortran results.
      PROGRAM NUMCHK
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 IB
      CHARACTER*1 F(11),G(11),H(11),P(11),C
      CHARACTER*11 S
      CHARACTER*80 LINE
   10 READ(*,'(A80)',END=90) LINE
      C = LINE(1:1)
      IF(C.EQ.'S') THEN
      S = LINE(3:13)
      DO I=1,11
      F(I)=S(I:I)
      ENDDO
      CALL IN9(X,F)
      IB = TRANSFER(X,IB)
      WRITE(*,'(Z16.16)') IB
      GO TO 10
      ENDIF
      READ(LINE(3:18),'(Z16)') IB
      X = TRANSFER(IB,X)
      DO I=1,11
      F(I)='#'
      G(I)='#'
      H(I)='#'
      P(I)='#'
      ENDDO
      CALL OUT9(X,F)
      CALL OUT9G(X,G)
      CALL OUT10(X,H)
      CALL OUT7(X,P)
      X1=X
      X2=X
      X3=X
      CALL INCORE10(X1)
      CALL INCORE9(X2)
      CALL INCORE7(X3)
      WRITE(*,20) F,G,H,P,TRANSFER(X1,IB),TRANSFER(X2,IB),
     1 TRANSFER(X3,IB)
   20 FORMAT(4(11A1,'|'),3(Z16.16,1X))
      GO TO 10
   90 STOP
      END
