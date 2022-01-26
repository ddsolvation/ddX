Module Complex_Bessel

!      REMARK ON ALGORITHM 644, COLLECTED ALGORITHMS FROM ACM.
!      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
!      VOL. 21, NO. 4, December, 1995, P.  388--393.

! Code converted using TO_F90 by Alan Miller
! Date: 2002-02-08  Time: 17:53:05
! Latest revision - 16 April 2002

IMPLICIT NONE
INTEGER, PARAMETER, PUBLIC  :: dp = SELECTED_REAL_KIND(12, 60)

PRIVATE
PUBLIC  :: cbesh, cbesi, cbesj, cbesk, cbesy, cairy, cbiry, gamln


CONTAINS


SUBROUTINE cbesh(z, fnu, kode, m, n, cy, nz, ierr)
!***BEGIN PROLOGUE  CBESH
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  890801, 930101   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  H-BESSEL FUNCTIONS,BESSEL FUNCTIONS OF COMPLEX ARGUMENT,
!             BESSEL FUNCTIONS OF THIRD KIND,HANKEL FUNCTIONS
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE THE H-BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!***DESCRIPTION

!   ON KODE=1, CBESH COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
!   HANKEL (BESSEL) FUNCTIONS CY(J)=H(M,FNU+J-1,Z) FOR KINDS M=1
!   OR 2, REAL, NONNEGATIVE ORDERS FNU+J-1, J=1,...,N, AND COMPLEX
!   Z.NE.CMPLX(0.0E0,0.0E0) IN THE CUT PLANE -PI < ARG(Z) <= PI.
!   ON KODE=2, CBESH COMPUTES THE SCALED HANKEL FUNCTIONS

!   CY(I)=H(M,FNU+J-1,Z)*EXP(-MM*Z*I)       MM=3-2M,      I**2=-1.

!   WHICH REMOVES THE EXPONENTIAL BEHAVIOR IN BOTH THE UPPER
!   AND LOWER HALF PLANES.  DEFINITIONS AND NOTATION ARE FOUND IN
!   THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).

!   INPUT
!     Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI < ARG(Z) <= PI
!     FNU    - ORDER OF INITIAL H FUNCTION, FNU >= 0.0E0
!     KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!              KODE= 1  RETURNS
!                       CY(J)=H(M,FNU+J-1,Z),      J=1,...,N
!                  = 2  RETURNS
!                       CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
!                            J=1,...,N  ,  I**2=-1
!     M      - KIND OF HANKEL FUNCTION, M=1 OR 2
!     N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1

!   OUTPUT
!     CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!              VALUES FOR THE SEQUENCE
!              CY(J)=H(M,FNU+J-1,Z)  OR
!              CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))  J=1,...,N
!              DEPENDING ON KODE, I**2=-1.
!     NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
!              NZ= 0   , NORMAL RETURN
!              NZ > 0 , FIRST NZ COMPONENTS OF CY SET TO ZERO DUE TO UNDERFLOW,
!                        CY(J)=CMPLX(0.0,0.0) J=1,...,NZ WHEN Y > 0.0 AND M=1
!                        OR Y < 0.0 AND M=2. FOR THE COMPLEMENTARY HALF PLANES,
!                        NZ STATES ONLY THE NUMBER OF UNDERFLOWS.
!     IERR    -ERROR FLAG
!              IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!              IERR=1, INPUT ERROR   - NO COMPUTATION
!              IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 TOO
!                      LARGE OR ABS(Z) TOO SMALL OR BOTH
!              IERR=3, ABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
!                      BUT LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!                      PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!              IERR=4, ABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF
!                      COMPLETE LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!              IERR=5, ERROR              - NO COMPUTATION,
!                      ALGORITHM TERMINATION CONDITION NOT MET

!***LONG DESCRIPTION

!    THE COMPUTATION IS CARRIED OUT BY THE RELATION

!    H(M,FNU,Z)=(1/MP)*EXP(-MP*FNU)*K(FNU,Z*EXP(-MP))
!        MP=MM*HPI*I,  MM=3-2*M,  HPI=PI/2,  I**2=-1

!    FOR M=1 OR 2 WHERE THE K BESSEL FUNCTION IS COMPUTED FOR THE
!    RIGHT HALF PLANE RE(Z) >= 0.0. THE K FUNCTION IS CONTINUED
!    TO THE LEFT HALF PLANE BY THE RELATION

!    K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
!    MP=MR*PI*I, MR=+1 OR -1, RE(Z) > 0, I**2=-1

!    WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.

!    EXPONENTIAL DECAY OF H(M,FNU,Z) OCCURS IN THE UPPER HALF Z PLANE FOR
!    M=1 AND THE LOWER HALF Z PLANE FOR M=2.  EXPONENTIAL GROWTH OCCURS IN THE
!    COMPLEMENTARY HALF PLANES.  SCALING BY EXP(-MM*Z*I) REMOVES THE
!    EXPONENTIAL BEHAVIOR IN THE WHOLE Z PLANE FOR Z TO INFINITY.

!    FOR NEGATIVE ORDERS,THE FORMULAE

!          H(1,-FNU,Z) = H(1,FNU,Z)*EXP( PI*FNU*I)
!          H(2,-FNU,Z) = H(2,FNU,Z)*EXP(-PI*FNU*I)
!                    I**2=-1

!    CAN BE USED.

!    IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELEMENTARY
!    FUNCTIONS. WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS LARGE, LOSSES OF
!    SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
!    CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN LOSSES
!    EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG IERR=3 IS
!    TRIGGERED WHERE UR=R1MACH(4)=UNIT ROUNDOFF.  ALSO IF EITHER IS LARGER
!    THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS LOST AND IERR=4.
!    IN ORDER TO USE THE INT FUNCTION, ARGUMENTS MUST BE FURTHER RESTRICTED
!    NOT TO EXCEED THE LARGEST MACHINE INTEGER, U3=I1MACH(9).
!    THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS RESTRICTED BY MIN(U2,U3).
!    ON 32 BIT MACHINES, U1,U2, AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6,
!    2.1E+9 IN SINGLE PRECISION ARITHMETIC AND 1.3E+8, 1.8D+16, 2.1E+9 IN
!    DOUBLE PRECISION ARITHMETIC RESPECTIVELY.  THIS MAKES U2 AND U3
!    LIMITING IN THEIR RESPECTIVE ARITHMETICS.  THIS MEANS THAT ONE CAN
!    EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
!    IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!    SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.

!    THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX BESSEL
!    FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT ROUNDOFF,1.0E-18)
!    IS THE NOMINAL PRECISION AND 10**S REPRESENTS THE INCREASE IN ERROR DUE
!    TO ARGUMENT REDUCTION IN THE ELEMENTARY FUNCTIONS.
!    HERE, S=MAX(1, ABS(LOG10(ABS(Z))), ABS(LOG10(FNU))) APPROXIMATELY
!    (I.E. S=MAX(1,ABS(EXPONENT OF ABS(Z),ABS(EXPONENT OF FNU)) ).
!    HOWEVER, THE PHASE ANGLE MAY HAVE ONLY ABSOLUTE ACCURACY.
!    THIS IS MOST LIKELY TO OCCUR WHEN ONE COMPONENT (IN ABSOLUTE VALUE)
!    IS LARGER THAN THE OTHER BY SEVERAL ORDERS OF MAGNITUDE.
!    IF ONE COMPONENT IS 10**K LARGER THAN THE OTHER, THEN ONE CAN EXPECT
!    ONLY MAX(ABS(LOG10(P))-K, 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER
!    WAY, WHEN K EXCEEDS THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN
!    THE SMALLER COMPONENT.  HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE
!    ACCURACY BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!    COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE MAGNITUDE
!    OF THE LARGER COMPONENT.  IN THESE EXTREME CASES, THE PRINCIPAL PHASE
!    ANGLE IS ON THE ORDER OF +P, -P, PI/2-P, OR -PI/2+P.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!        I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!      COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!        BY D. E. AMOS, SAND83-0083, MAY 1983.

!      COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!        AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY 1983

!      A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!        AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-1018, MAY, 1985

!      A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!        AND NONNEGATIVE ORDER BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE,
!        VOL. 12, NO. 3, SEPTEMBER 1986, PP 265-273.

!***ROUTINES CALLED  CACON,CBKNU,CBUNK,CUOIK,I1MACH,R1MACH
!***END PROLOGUE  CBESH

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: cy(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: ierr

COMPLEX (dp)  :: zn, zt, csgn
REAL (dp)     :: aa, alim, aln, arg, az, cpn, dig, elim, fmm, fn, fnul,  &
                 rhpi, rl, r1m5, sgn, spn, tol, ufl, xn, xx, yn, yy,  &
                 bb, ascle, rtol, atol
INTEGER       :: i, inu, inuh, ir, k, k1, k2, mm, mr, nn, nuf, nw

REAL (dp), PARAMETER  :: hpi = 1.57079632679489662_dp

!***FIRST EXECUTABLE STATEMENT  CBESH
nz = 0
xx = REAL(z, KIND=dp)
yy = AIMAG(z)
ierr = 0
IF (xx == 0.0_dp .AND. yy == 0.0_dp) ierr = 1
IF (fnu < 0.0_dp) ierr = 1
IF (m < 1 .OR. m > 2) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (n < 1) ierr = 1
IF (ierr /= 0) RETURN
nn = n
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
!-----------------------------------------------------------------------
tol = MAX(EPSILON(0.0_dp), 1.0D-18)
k1 = MINEXPONENT(0.0_dp)
k2 = MAXEXPONENT(0.0_dp)
r1m5 = LOG10( REAL( RADIX(0.0_dp), KIND=dp) )
k = MIN(ABS(k1), ABS(k2))
elim = 2.303_dp * (k*r1m5 - 3.0_dp)
k1 = DIGITS(0.0_dp) - 1
aa = r1m5 * k1
dig = MIN(aa, 18.0_dp)
aa = aa * 2.303_dp
alim = elim + MAX(-aa, -41.45_dp)
fnul = 10.0_dp + 6.0_dp * (dig - 3.0_dp)
rl = 1.2_dp * dig + 3.0_dp
fn = fnu + (nn-1)
mm = 3 - m - m
fmm = mm
zn = z * CMPLX(0.0_dp, -fmm, KIND=dp)
xn = REAL(zn, KIND=dp)
yn = AIMAG(zn)
az = ABS(z)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
aa = 0.5_dp / tol
bb = HUGE(0) * 0.5_dp
aa = MIN(aa,bb)
IF (az <= aa) THEN
  IF (fn <= aa) THEN
    aa = SQRT(aa)
    IF (az > aa) ierr = 3
    IF (fn > aa) ierr = 3
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
!-----------------------------------------------------------------------
    ufl = TINY(0.0_dp) * 1.0D+3
    IF (az >= ufl) THEN
      IF (fnu <= fnul) THEN
        IF (fn > 1.0_dp) THEN
          IF (fn <= 2.0_dp) THEN
            IF (az > tol) GO TO 10
            arg = 0.5_dp * az
            aln = -fn * LOG(arg)
            IF (aln > elim) GO TO 50
          ELSE
            CALL cuoik(zn, fnu, kode, 2, nn, cy, nuf, tol, elim, alim)
            IF (nuf < 0) GO TO 50
            nz = nz + nuf
            nn = nn - nuf
!-----------------------------------------------------------------------
!     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
!     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
!-----------------------------------------------------------------------
            IF (nn == 0) GO TO 40
          END IF
        END IF

        10 IF (.NOT.(xn < 0.0_dp .OR. (xn == 0.0_dp .AND. yn < 0.0_dp  &
               .AND. m == 2))) THEN
!-----------------------------------------------------------------------
!     RIGHT HALF PLANE COMPUTATION, XN >= 0.  .AND.  (XN.NE.0.  .OR.
!     YN >= 0.  .OR.  M=1)
!-----------------------------------------------------------------------
          CALL cbknu(zn, fnu, kode, nn, cy, nz, tol, elim, alim)
          GO TO 20
        END IF
!-----------------------------------------------------------------------
!     LEFT HALF PLANE COMPUTATION
!-----------------------------------------------------------------------
        mr = -mm
        CALL cacon(zn, fnu, kode, mr, nn, cy, nw, rl, fnul, tol, elim, alim)
        IF (nw < 0) GO TO 60
        nz = nw
      ELSE
!-----------------------------------------------------------------------
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
!-----------------------------------------------------------------------
        mr = 0
        IF (.NOT.(xn >= 0.0_dp .AND. (xn /= 0.0_dp .OR. yn >= 0.0_dp  &
               .OR. m /= 2))) THEN
          mr = -mm
          IF (xn == 0.0_dp .AND. yn < 0.0_dp) zn = -zn
        END IF
        CALL cbunk(zn, fnu, kode, mr, nn, cy, nw, tol, elim, alim)
        IF (nw < 0) GO TO 60
        nz = nz + nw
      END IF
!-----------------------------------------------------------------------
!     H(M,FNU,Z) = -FMM*(I/HPI)*(ZT**FNU)*K(FNU,-Z*ZT)

!     ZT=EXP(-FMM*HPI*I) = CMPLX(0.0,-FMM), FMM=3-2*M, M=1,2
!-----------------------------------------------------------------------
      20 sgn = SIGN(hpi,-fmm)
!-----------------------------------------------------------------------
!     CALCULATE EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
      inu = fnu
      inuh = inu / 2
      ir = inu - 2 * inuh
      arg = (fnu - (inu-ir)) * sgn
      rhpi = 1.0_dp / sgn
      cpn = rhpi * COS(arg)
      spn = rhpi * SIN(arg)
!     ZN = CMPLX(-SPN,CPN)
      csgn = CMPLX(-spn, cpn, KIND=dp)
!     IF (MOD(INUH,2).EQ.1) ZN = -ZN
      IF (MOD(inuh,2) == 1) csgn = -csgn
      zt = CMPLX(0.0_dp, -fmm, KIND=dp)
      rtol = 1.0_dp / tol
      ascle = ufl * rtol
      DO  i = 1, nn
!       CY(I) = CY(I)*ZN
!       ZN = ZN*ZT
        zn = cy(i)
        aa = REAL(zn, KIND=dp)
        bb = AIMAG(zn)
        atol = 1.0_dp
        IF (MAX(ABS(aa),ABS(bb)) <= ascle) THEN
          zn = zn * rtol
          atol = tol
        END IF
        zn = zn * csgn
        cy(i) = zn * atol
        csgn = csgn * zt
      END DO
      RETURN

      40 IF (xn >= 0.0_dp) RETURN
    END IF

    50 ierr = 2
    nz = 0
    RETURN

    60 IF (nw == -1) GO TO 50
    nz = 0
    ierr = 5
    RETURN
  END IF
END IF
nz = 0
ierr = 4
RETURN
END SUBROUTINE cbesh



SUBROUTINE cbesi(z, fnu, kode, n, cy, nz, ierr)
!***BEGIN PROLOGUE  CBESI
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  890801, 930101   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  I-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION OF THE FIRST KIND
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE I-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!***DESCRIPTION

!   ON KODE=1, CBESI COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
!   BESSEL FUNCTIONS CY(J)=I(FNU+J-1,Z) FOR REAL (dp), NONNEGATIVE
!   ORDERS FNU+J-1, J=1,...,N AND COMPLEX Z IN THE CUT PLANE
!   -PI < ARG(Z) <= PI. ON KODE=2, CBESI RETURNS THE SCALED FUNCTIONS

!   CY(J)=EXP(-ABS(X))*I(FNU+J-1,Z)   J = 1,...,N , X=REAL(Z)

!   WITH THE EXPONENTIAL GROWTH REMOVED IN BOTH THE LEFT AND
!   RIGHT HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND
!   NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
!   FUNCTIONS (REF.1)

!   INPUT
!     Z      - Z=CMPLX(X,Y),  -PI < ARG(Z) <= PI
!     FNU    - ORDER OF INITIAL I FUNCTION, FNU >= 0.0_dp
!     KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!              KODE= 1  RETURNS
!                       CY(J)=I(FNU+J-1,Z), J=1,...,N
!                  = 2  RETURNS
!                       CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X)), J=1,...,N
!     N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1

!   OUTPUT
!     CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!              VALUES FOR THE SEQUENCE
!              CY(J)=I(FNU+J-1,Z)  OR
!              CY(J)=I(FNU+J-1,Z)*EXP(-ABS(X))  J=1,...,N
!              DEPENDING ON KODE, X=REAL(Z)
!     NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
!              NZ= 0   , NORMAL RETURN
!              NZ > 0 , LAST NZ COMPONENTS OF CY SET TO ZERO
!                        DUE TO UNDERFLOW, CY(J)=CMPLX(0.0,0.0),
!                        J = N-NZ+1,...,N
!     IERR   - ERROR FLAG
!              IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!              IERR=1, INPUT ERROR   - NO COMPUTATION
!              IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z) TOO
!                      LARGE ON KODE=1
!              IERR=3, ABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
!                      BUT LOSSES OF SIGNIFICANCE BY ARGUMENT
!                      REDUCTION PRODUCE LESS THAN HALF OF MACHINE
!                      ACCURACY
!              IERR=4, ABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTA-
!                      TION BECAUSE OF COMPLETE LOSSES OF SIGNIFI-
!                      CANCE BY ARGUMENT REDUCTION
!              IERR=5, ERROR              - NO COMPUTATION,
!                      ALGORITHM TERMINATION CONDITION NOT MET

!***LONG DESCRIPTION

!   THE COMPUTATION IS CARRIED OUT BY THE POWER SERIES FOR
!   SMALL ABS(Z), THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z),
!   THE MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN AND A
!   NEUMANN SERIES FOR IMTERMEDIATE MAGNITUDES, AND THE
!   UNIFORM ASYMPTOTIC EXPANSIONS FOR I(FNU,Z) AND J(FNU,Z)
!   FOR LARGE ORDERS. BACKWARD RECURRENCE IS USED TO GENERATE
!   SEQUENCES OR REDUCE ORDERS WHEN NECESSARY.

!   THE CALCULATIONS ABOVE ARE DONE IN THE RIGHT HALF PLANE AND
!   CONTINUED INTO THE LEFT HALF PLANE BY THE FORMULA

!   I(FNU,Z*EXP(M*PI)) = EXP(M*PI*FNU)*I(FNU,Z)  REAL(Z) > 0.0
!                 M = +I OR -I,  I**2=-1

!   FOR NEGATIVE ORDERS,THE FORMULA

!        I(-FNU,Z) = I(FNU,Z) + (2/PI)*SIN(PI*FNU)*K(FNU,Z)

!   CAN BE USED.  HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE FUNCTION
!   CHANGES RADICALLY.  WHEN FNU IS A LARGE POSITIVE INTEGER,THE MAGNITUDE OF
!   I(-FNU,Z) = I(FNU,Z) IS A LARGE NEGATIVE POWER OF TEN.  BUT WHEN FNU IS NOT
!   AN INTEGER, K(FNU,Z) DOMINATES IN MAGNITUDE WITH A LARGE POSITIVE POWER OF
!   TEN AND THE MOST THAT THE SECOND TERM CAN BE REDUCED IS BY UNIT ROUNDOFF
!   FROM THE COEFFICIENT. THUS, WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF
!   A LARGE INTEGER FOR FNU.  HERE, LARGE MEANS FNU > ABS(Z).

!   IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELEMENTARY
!   FUNCTIONS.  WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS
!   LARGE, LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
!   CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN
!   LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR FLAG
!   IERR=3 IS TRIGGERED WHERE UR=EPSILON(0.0_dp)=UNIT ROUNDOFF.  ALSO
!   IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS
!   LOST AND IERR=4. IN ORDER TO USE THE INT FUNCTION, ARGUMENTS
!   MUST BE FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE
!   INTEGER, U3=HUGE(0). THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS
!   RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2, AND U3
!   ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION
!   ARITHMETIC AND 1.3E+8, 1.8D+16, 2.1E+9 IN DOUBLE PRECISION
!   ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
!   THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT ONE CAN EXPECT
!   TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS
!   IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!   SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.

!   THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX BESSEL
!   FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P = MAX(UNIT ROUNDOFF,1.0E-18)
!   IS THE NOMINAL PRECISION AND 10**S REPRESENTS THE INCREASE IN ERROR DUE TO
!   ARGUMENT REDUCTION IN THE ELEMENTARY FUNCTIONS.  HERE, S =
!   MAX(1,ABS(LOG10(ABS(Z))), ABS(LOG10(FNU))) APPROXIMATELY
!   (I.E. S = MAX(1,ABS(EXPONENT OF ABS(Z), ABS(EXPONENT OF FNU)) ).
!   HOWEVER, THE PHASE ANGLE MAY HAVE ONLY ABSOLUTE ACCURACY.
!   THIS IS MOST LIKELY TO OCCUR WHEN ONE COMPONENT (IN ABSOLUTE VALUE) IS
!   LARGER THAN THE OTHER BY SEVERAL ORDERS OF MAGNITUDE.
!   IF ONE COMPONENT IS 10**K LARGER THAN THE OTHER, THEN ONE CAN EXPECT ONLY
!   MAX(ABS(LOG10(P))-K, 0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K
!   EXCEEDS THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!   COMPONENT.  HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY BECAUSE,
!   IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER COMPONENT WILL NOT
!   (AS A RULE) DECREASE BELOW P TIMES THE MAGNITUDE OF THE LARGER COMPONENT.
!   IN THESE EXTREME CASES, THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF
!   +P, -P, PI/2-P, OR -PI/2+P.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!         I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!       COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!         BY D. E. AMOS, SAND83-0083, MAY 1983.

!       COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!         AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY 1983

!       A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!         AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-1018, MAY 1985

!       A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!         AND NONNEGATIVE ORDER BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE,
!         VOL. 12, NO. 3, SEPTEMBER 1986, PP 265-273.

!***ROUTINES CALLED  CBINU,I1MACH,R1MACH
!***END PROLOGUE  CBESI

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: cy(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: ierr

COMPLEX (dp)  :: csgn, zn
REAL (dp)     :: aa, alim, arg, dig, elim, fnul, rl, r1m5, s1, s2,  &
                 tol, xx, yy, az, fn, bb, ascle, rtol, atol
INTEGER       :: i, inu, k, k1, k2, nn

REAL (dp), PARAMETER     :: pi = 3.14159265358979324_dp
COMPLEX (dp), PARAMETER  :: cone = (1.0_dp, 0.0_dp)

!***FIRST EXECUTABLE STATEMENT  CBESI
ierr = 0
nz = 0
IF (fnu < 0.0_dp) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (n < 1) ierr = 1
IF (ierr /= 0) RETURN
xx = REAL(z, KIND=dp)
yy = AIMAG(z)
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
!-----------------------------------------------------------------------
tol = MAX(EPSILON(0.0_dp), 1.0D-18)
k1 = MINEXPONENT(0.0_dp)
k2 = MAXEXPONENT(0.0_dp)
r1m5 = LOG10( REAL( RADIX(0.0_dp), KIND=dp) )
k = MIN(ABS(k1),ABS(k2))
elim = 2.303_dp * (k*r1m5 - 3.0_dp)
k1 = DIGITS(0.0_dp) - 1
aa = r1m5 * k1
dig = MIN(aa, 18.0_dp)
aa = aa * 2.303_dp
alim = elim + MAX(-aa, -41.45_dp)
rl = 1.2_dp * dig + 3.0_dp
fnul = 10.0_dp + 6.0_dp * (dig - 3.0_dp)
az = ABS(z)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
aa = 0.5_dp / tol
bb = HUGE(0) * 0.5_dp
aa = MIN(aa,bb)
IF (az <= aa) THEN
  fn = fnu + (n-1)
  IF (fn <= aa) THEN
    aa = SQRT(aa)
    IF (az > aa) ierr = 3
    IF (fn > aa) ierr = 3
    zn = z
    csgn = cone
    IF (xx < 0.0_dp) THEN
      zn = -z
!-----------------------------------------------------------------------
!     CALCULATE CSGN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
      inu = fnu
      arg = (fnu - inu) * pi
      IF (yy < 0.0_dp) arg = -arg
      s1 = COS(arg)
      s2 = SIN(arg)
      csgn = CMPLX(s1, s2, KIND=dp)
      IF (MOD(inu,2) == 1) csgn = -csgn
    END IF
    CALL cbinu(zn, fnu, kode, n, cy, nz, rl, fnul, tol, elim, alim)
    IF (nz >= 0) THEN
      IF (xx >= 0.0_dp) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE
!-----------------------------------------------------------------------
      nn = n - nz
      IF (nn == 0) RETURN
      rtol = 1.0_dp / tol
      ascle = TINY(0.0_dp) * rtol * 1.0E+3
      DO  i = 1, nn
!       CY(I) = CY(I)*CSGN
        zn = cy(i)
        aa = REAL(zn, KIND=dp)
        bb = AIMAG(zn)
        atol = 1.0_dp
        IF (MAX(ABS(aa),ABS(bb)) <= ascle) THEN
          zn = zn * rtol
          atol = tol
        END IF
        zn = zn * csgn
        cy(i) = zn * atol
        csgn = -csgn
      END DO
      RETURN
    END IF
    IF (nz /= -2) THEN
      nz = 0
      ierr = 2
      RETURN
    END IF
    nz = 0
    ierr = 5
    RETURN
  END IF
END IF
nz = 0
ierr = 4
RETURN
END SUBROUTINE cbesi



SUBROUTINE cbesj(z, fnu, kode, n, cy, nz, ierr)
!***BEGIN PROLOGUE  CBESJ
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  890801, 930101   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  J-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
!             BESSEL FUNCTION OF FIRST KIND
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE THE J-BESSEL FUNCTION OF A COMPLEX ARGUMENT
!***DESCRIPTION

!    ON KODE=1, CBESJ COMPUTES AN N MEMBER  SEQUENCE OF COMPLEX
!    BESSEL FUNCTIONS CY(I) = J(FNU+I-1,Z) FOR REAL (dp), NONNEGATIVE
!    ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
!    -PI < ARG(Z) <= PI.  ON KODE=2, CBESJ RETURNS THE SCALED FUNCTIONS

!    CY(I) = EXP(-ABS(Y))*J(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)

!    WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
!    LOWER HALF PLANES FOR Z TO INFINITY.  DEFINITIONS AND NOTATION
!    ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).

!    INPUT
!      Z      - Z=CMPLX(X,Y),  -PI < ARG(Z) <= PI
!      FNU    - ORDER OF INITIAL J FUNCTION, FNU >= 0.0_dp
!      KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!               KODE= 1  RETURNS
!                        CY(I)=J(FNU+I-1,Z), I=1,...,N
!                   = 2  RETURNS
!                        CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...
!      N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1

!    OUTPUT
!      CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!               VALUES FOR THE SEQUENCE
!               CY(I)=J(FNU+I-1,Z)  OR
!               CY(I)=J(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
!               DEPENDING ON KODE, Y=AIMAG(Z).
!      NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW,
!               NZ= 0   , NORMAL RETURN
!               NZ > 0 , LAST NZ COMPONENTS OF CY SET TO ZERO
!                         DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
!                         I = N-NZ+1,...,N
!      IERR   - ERROR FLAG
!               IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!               IERR=1, INPUT ERROR   - NO COMPUTATION
!               IERR=2, OVERFLOW      - NO COMPUTATION, AIMAG(Z)
!                       TOO LARGE ON KODE=1
!               IERR=3, ABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
!                       BUT LOSSES OF SIGNIFICANCE BY ARGUMENT
!                       REDUCTION PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!               IERR=4, ABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE
!                       OF COMPLETE LOSSES OF SIGNIFICANCE BY ARGUMENT
!                       REDUCTION
!               IERR=5, ERROR              - NO COMPUTATION,
!                       ALGORITHM TERMINATION CONDITION NOT MET

!***LONG DESCRIPTION

!    THE COMPUTATION IS CARRIED OUT BY THE FORMULA

!    J(FNU,Z)=EXP( FNU*PI*I/2)*I(FNU,-I*Z)    AIMAG(Z) >= 0.0

!    J(FNU,Z)=EXP(-FNU*PI*I/2)*I(FNU, I*Z)    AIMAG(Z) < 0.0

!    WHERE I**2 = -1 AND I(FNU,Z) IS THE I BESSEL FUNCTION.

!    FOR NEGATIVE ORDERS,THE FORMULA

!         J(-FNU,Z) = J(FNU,Z)*COS(PI*FNU) - Y(FNU,Z)*SIN(PI*FNU)

!    CAN BE USED.  HOWEVER,FOR LARGE ORDERS CLOSE TO INTEGERS, THE FUNCTION
!    CHANGES RADICALLY.  WHEN FNU IS A LARGE POSITIVE INTEGER, THE MAGNITUDE
!    OF J(-FNU,Z)=J(FNU,Z)*COS(PI*FNU) IS A LARGE NEGATIVE POWER OF TEN.
!    BUT WHEN FNU IS NOT AN INTEGER, Y(FNU,Z) DOMINATES IN MAGNITUDE WITH A
!    LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE SECOND TERM CAN BE
!    REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT.  THUS, WIDE CHANGES CAN
!    OCCUR WITHIN UNIT ROUNDOFF OF A LARGE INTEGER FOR FNU.  HERE, LARGE MEANS
!    FNU > ABS(Z).

!    IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELEMENTARY
!    FUNCTIONS.  WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS LARGE, LOSSES OF
!    SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.  CONSEQUENTLY, IF EITHER ONE
!    EXCEEDS U1=SQRT(0.5/UR), THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY
!    AND AN ERROR FLAG IERR=3 IS TRIGGERED WHERE UR = EPSILON(0.0_dp) = UNIT
!    ROUNDOFF.  ALSO IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE
!    IS LOST AND IERR=4.  IN ORDER TO USE THE INT FUNCTION, ARGUMENTS MUST BE
!    FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE INTEGER, U3 = HUGE(0).
!    THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS RESTRICTED BY MIN(U2,U3).
!    ON 32 BIT MACHINES, U1,U2, AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9
!    IN SINGLE PRECISION ARITHMETIC AND 1.3E+8, 1.8D+16, 2.1E+9 IN DOUBLE
!    PRECISION ARITHMETIC RESPECTIVELY.  THIS MAKES U2 AND U3 LIMITING IN
!    THEIR RESPECTIVE ARITHMETICS.  THIS MEANS THAT ONE CAN EXPECT TO RETAIN,
!    IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS IN SINGLE AND ONLY 7
!    DIGITS IN DOUBLE PRECISION ARITHMETIC.
!    SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.

!    THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX BESSEL
!    FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P = MAX(UNIT ROUNDOFF, 1.0E-18)
!    IS THE NOMINAL PRECISION AND 10**S REPRESENTS THE INCREASE IN ERROR DUE
!    TO ARGUMENT REDUCTION IN THE ELEMENTARY FUNCTIONS.  HERE,
!    S = MAX(1,ABS(LOG10(ABS(Z))), ABS(LOG10(FNU))) APPROXIMATELY
!    (I.E. S = MAX(1, ABS(EXPONENT OF ABS(Z), ABS(EXPONENT OF FNU)) ).
!    HOWEVER, THE PHASE ANGLE MAY HAVE ONLY ABSOLUTE ACCURACY.  THIS IS MOST
!    LIKELY TO OCCUR WHEN ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN
!    THE OTHER BY SEVERAL ORDERS OF MAGNITUDE.  IF ONE COMPONENT IS 10**K
!    LARGER THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K, 0)
!    SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS THE EXPONENT
!    OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER COMPONENT.
!    HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY BECAUSE, IN COMPLEX
!    ARITHMETIC WITH PRECISION P, THE SMALLER COMPONENT WILL NOT (AS A RULE)
!    DECREASE BELOW P TIMES THE MAGNITUDE OF THE LARGER COMPONENT.
!    IN THESE EXTREME CASES, THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P,
!    -P, PI/2-P, OR -PI/2+P.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!         I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!       COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!         BY D. E. AMOS, SAND83-0083, MAY 1983.

!       COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!         AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY 1983

!       A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!         AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-1018, MAY 1985

!       A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!         AND NONNEGATIVE ORDER BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE,
!         VOL. 12, NO. 3, SEPTEMBER 1986, PP 265-273.

!***ROUTINES CALLED  CBINU,I1MACH,R1MACH
!***END PROLOGUE  CBESJ

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: cy(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: ierr

COMPLEX (dp)  :: ci, csgn, zn
REAL (dp)     :: aa, alim, arg, dig, elim, fnul, rl, r1, r1m5, r2,  &
                 tol, yy, az, fn, bb, ascle, rtol, atol
INTEGER       :: i, inu, inuh, ir, k1, k2, nl, k

REAL (dp), PARAMETER  :: hpi = 1.570796326794896619_dp

!***FIRST EXECUTABLE STATEMENT  CBESJ
ierr = 0
nz = 0
IF (fnu < 0.0_dp) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (n < 1) ierr = 1
IF (ierr /= 0) RETURN
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
!-----------------------------------------------------------------------
tol = MAX(EPSILON(0.0_dp), 1.0D-18)
k1 = MINEXPONENT(0.0_dp)
k2 = MAXEXPONENT(0.0_dp)
r1m5 = LOG10( REAL( RADIX(0.0_dp), KIND=dp) )
k = MIN(ABS(k1),ABS(k2))
elim = 2.303_dp * (k*r1m5 - 3.0_dp)
k1 = DIGITS(0.0_dp) - 1
aa = r1m5 * k1
dig = MIN(aa, 18.0_dp)
aa = aa * 2.303_dp
alim = elim + MAX(-aa, -41.45_dp)
rl = 1.2_dp * dig + 3.0_dp
fnul = 10.0_dp + 6.0_dp * (dig - 3.0_dp)
ci = CMPLX(0.0_dp, 1.0_dp, KIND=dp)
yy = AIMAG(z)
az = ABS(z)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
aa = 0.5_dp / tol
bb = HUGE(0) * 0.5_dp
aa = MIN(aa,bb)
fn = fnu + (n-1)
IF (az <= aa) THEN
  IF (fn <= aa) THEN
    aa = SQRT(aa)
    IF (az > aa) ierr = 3
    IF (fn > aa) ierr = 3
!-----------------------------------------------------------------------
!     CALCULATE CSGN=EXP(FNU*HPI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
    inu = fnu
    inuh = inu / 2
    ir = inu - 2 * inuh
    arg = (fnu - (inu-ir)) * hpi
    r1 = COS(arg)
    r2 = SIN(arg)
    csgn = CMPLX(r1, r2, KIND=dp)
    IF (MOD(inuh,2) == 1) csgn = -csgn
!-----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE
!-----------------------------------------------------------------------
    zn = -z * ci
    IF (yy < 0.0_dp) THEN
      zn = -zn
      csgn = CONJG(csgn)
      ci = CONJG(ci)
    END IF
    CALL cbinu(zn, fnu, kode, n, cy, nz, rl, fnul, tol, elim, alim)
    IF (nz >= 0) THEN
      nl = n - nz
      IF (nl == 0) RETURN
      rtol = 1.0_dp / tol
      ascle = TINY(0.0_dp) * rtol * 1.0E+3
      DO  i = 1, nl
!       CY(I)=CY(I)*CSGN
        zn = cy(i)
        aa = REAL(zn, KIND=dp)
        bb = AIMAG(zn)
        atol = 1.0_dp
        IF (MAX(ABS(aa),ABS(bb)) <= ascle) THEN
          zn = zn * rtol
          atol = tol
        END IF
        zn = zn * csgn
        cy(i) = zn * atol
        csgn = csgn * ci
      END DO
      RETURN
    END IF
    IF (nz /= -2) THEN
      nz = 0
      ierr = 2
      RETURN
    END IF
    nz = 0
    ierr = 5
    RETURN
  END IF
END IF
nz = 0
ierr = 4
RETURN
END SUBROUTINE cbesj



SUBROUTINE cbesk(z, fnu, kode, n, cy, nz, ierr)
!***BEGIN PROLOGUE  CBESK
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  890801, 930101   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  K-BESSEL FUNCTION,COMPLEX BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION OF THE SECOND KIND,
!             BESSEL FUNCTION OF THE THIRD KIND
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE K-BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!***DESCRIPTION

!   ON KODE=1, CBESK COMPUTES AN N MEMBER SEQUENCE OF COMPLEX BESSEL FUNCTIONS
!   CY(J)=K(FNU+J-1,Z) FOR REAL (dp), NONNEGATIVE ORDERS FNU+J-1, J=1,...,N
!   AND COMPLEX Z.NE.CMPLX(0.0,0.0) IN THE CUT PLANE -PI < ARG(Z) <= PI.
!   ON KODE=2, CBESK RETURNS THE SCALED K FUNCTIONS,

!   CY(J)=EXP(Z)*K(FNU+J-1,Z) , J=1,...,N,

!   WHICH REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND RIGHT HALF
!   PLANES FOR Z TO INFINITY.  DEFINITIONS AND NOTATION ARE FOUND IN THE NBS
!   HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).

!   INPUT
!     Z      - Z=CMPLX(X,Y),Z.NE.CMPLX(0.,0.),-PI < ARG(Z) <= PI
!     FNU    - ORDER OF INITIAL K FUNCTION, FNU >= 0.0_dp
!     N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
!     KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!              KODE= 1  RETURNS
!                       CY(I)=K(FNU+I-1,Z), I=1,...,N
!                  = 2  RETURNS
!                       CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N

!   OUTPUT
!     CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!              VALUES FOR THE SEQUENCE
!              CY(I)=K(FNU+I-1,Z), I=1,...,N OR
!              CY(I)=K(FNU+I-1,Z)*EXP(Z), I=1,...,N
!              DEPENDING ON KODE
!     NZ     - NUMBER OF COMPONENTS SET TO ZERO DUE TO UNDERFLOW.
!              NZ= 0   , NORMAL RETURN
!              NZ > 0 , FIRST NZ COMPONENTS OF CY SET TO ZERO
!                        DUE TO UNDERFLOW, CY(I)=CMPLX(0.0,0.0),
!                        I=1,...,N WHEN X >= 0.0.  WHEN X < 0.0, NZ STATES
!                        ONLY THE NUMBER OF UNDERFLOWS IN THE SEQUENCE.
!     IERR   - ERROR FLAG
!              IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!              IERR=1, INPUT ERROR   - NO COMPUTATION
!              IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
!                      TOO LARGE OR ABS(Z) IS TOO SMALL OR BOTH
!              IERR=3, ABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE, BUT
!                      LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION PRODUCE
!                      LESS THAN HALF OF MACHINE ACCURACY
!              IERR=4, ABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION BECAUSE OF
!                      COMPLETE LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!              IERR=5, ERROR              - NO COMPUTATION,
!                      ALGORITHM TERMINATION CONDITION NOT MET

!***LONG DESCRIPTION

!   EQUATIONS OF THE REFERENCE ARE IMPLEMENTED FOR SMALL ORDERS DNU AND
!   DNU+1.0 IN THE RIGHT HALF PLANE X >= 0.0.  FORWARD RECURRENCE GENERATES
!   HIGHER ORDERS.  K IS CONTINUED TO THE LEFT HALF PLANE BY THE RELATION

!   K(FNU,Z*EXP(MP)) = EXP(-MP*FNU)*K(FNU,Z)-MP*I(FNU,Z)
!   MP=MR*PI*I, MR=+1 OR -1, RE(Z) > 0, I**2=-1

!   WHERE I(FNU,Z) IS THE I BESSEL FUNCTION.

!   FOR LARGE ORDERS, FNU > FNUL, THE K FUNCTION IS COMPUTED
!   BY MEANS OF ITS UNIFORM ASYMPTOTIC EXPANSIONS.

!   FOR NEGATIVE ORDERS, THE FORMULA

!                 K(-FNU,Z) = K(FNU,Z)

!   CAN BE USED.

!   CBESK ASSUMES THAT A SIGNIFICANT DIGIT SINH(X) FUNCTION IS AVAILABLE.

!   IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELEMENTARY
!   FUNCTIONS.  WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS LARGE, LOSSES OF
!   SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.
!   CONSEQUENTLY, IF EITHER ONE EXCEEDS U1=SQRT(0.5/UR), THEN LOSSES EXCEEDING
!   HALF PRECISION ARE LIKELY AND AN ERROR FLAG IERR=3 IS TRIGGERED WHERE
!   UR = EPSILON(0.0_dp) = UNIT ROUNDOFF.  ALSO IF EITHER IS LARGER THAN
!   U2 = 0.5/UR, THEN ALL SIGNIFICANCE IS LOST AND IERR=4.  IN ORDER TO USE
!   THE INT FUNCTION, ARGUMENTS MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
!   LARGEST MACHINE INTEGER, U3=HUGE(0).  THUS, THE MAGNITUDE OF Z AND FNU+N-1
!   IS RESTRICTED BY MIN(U2,U3).  ON 32 BIT MACHINES, U1,U2, AND U3 ARE
!   APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION ARITHMETIC AND
!   1.3E+8, 1.8D+16, 2.1E+9 IN DOUBLE PRECISION ARITHMETIC RESPECTIVELY.
!   THIS MAKES U2 AND U3 LIMITING IN THEIR RESPECTIVE ARITHMETICS.  THIS MEANS
!   THAT ONE CAN EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO
!   DIGITS IN SINGLE AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!   SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.

!   THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX BESSEL
!   FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P = MAX(UNIT ROUNDOFF,1.0E-18)
!   IS THE NOMINAL PRECISION AND 10**S REPRESENTS THE INCREASE IN ERROR DUE TO
!   ARGUMENT REDUCTION IN THE ELEMENTARY FUNCTIONS.  HERE, S =
!   MAX(1,ABS(LOG10(ABS(Z))), ABS(LOG10(FNU))) APPROXIMATELY (I.E. S =
!   MAX(1,ABS(EXPONENT OF ABS(Z),ABS(EXPONENT OF FNU)) ).
!   HOWEVER, THE PHASE ANGLE MAY HAVE ONLY ABSOLUTE ACCURACY.  THIS IS MOST
!   LIKELY TO OCCUR WHEN ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE
!   OTHER BY SEVERAL ORDERS OF MAGNITUDE.  IF ONE COMPONENT IS 10**K LARGER
!   THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K, 0)
!   SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS THE EXPONENT
!   OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER COMPONENT.  HOWEVER, THE
!   PHASE ANGLE RETAINS ABSOLUTE ACCURACY BECAUSE, IN COMPLEX ARITHMETIC WITH
!   PRECISION P, THE SMALLER COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P
!   TIMES THE MAGNITUDE OF THE LARGER COMPONENT.  IN THESE EXTREME CASES,
!   THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P, OR -PI/2+P.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!           I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!         COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!           BY D. E. AMOS, SAND83-0083, MAY 1983.

!         COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!           AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY 1983.

!         A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!           AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-1018, MAY 1985

!         A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!           AND NONNEGATIVE ORDER BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE,
!           VOL. 12, NO. 3, SEPTEMBER 1986, PP 265-273.

!***ROUTINES CALLED  CACON,CBKNU,CBUNK,CUOIK,I1MACH,R1MACH
!***END PROLOGUE  CBESK

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: cy(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: ierr

REAL (dp)  :: aa, alim, aln, arg, az, dig, elim, fn, fnul, rl, r1m5,  &
              tol, ufl, xx, yy, bb
INTEGER    :: k, k1, k2, mr, nn, nuf, nw

!***FIRST EXECUTABLE STATEMENT  CBESK
ierr = 0
nz = 0
xx = REAL(z, KIND=dp)
yy = AIMAG(z)
IF (yy == 0.0_dp .AND. xx == 0.0_dp) ierr = 1
IF (fnu < 0.0_dp) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (n < 1) ierr = 1
IF (ierr /= 0) RETURN
nn = n
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU
!-----------------------------------------------------------------------
tol = MAX(EPSILON(0.0_dp), 1.0D-18)
k1 = MINEXPONENT(0.0_dp)
k2 = MAXEXPONENT(0.0_dp)
r1m5 = LOG10( REAL( RADIX(0.0_dp), KIND=dp) )
k = MIN(ABS(k1),ABS(k2))
elim = 2.303_dp * (k*r1m5 - 3.0_dp)
k1 = DIGITS(0.0_dp) - 1
aa = r1m5 * k1
dig = MIN(aa, 18.0_dp)
aa = aa * 2.303_dp
alim = elim + MAX(-aa, -41.45_dp)
fnul = 10.0_dp + 6.0_dp * (dig - 3.0_dp)
rl = 1.2_dp * dig + 3.0_dp
az = ABS(z)
fn = fnu + (nn-1)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
aa = 0.5_dp / tol
bb = HUGE(0) * 0.5_dp
aa = MIN(aa,bb)
IF (az <= aa) THEN
  IF (fn <= aa) THEN
    aa = SQRT(aa)
    IF (az > aa) ierr = 3
    IF (fn > aa) ierr = 3
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON THE LAST MEMBER OF THE SEQUENCE
!-----------------------------------------------------------------------
!     UFL = EXP(-ELIM)
    ufl = TINY(0.0_dp) * 1.0E+3
    IF (az >= ufl) THEN
      IF (fnu <= fnul) THEN
        IF (fn > 1.0_dp) THEN
          IF (fn <= 2.0_dp) THEN
            IF (az > tol) GO TO 10
            arg = 0.5_dp * az
            aln = -fn * LOG(arg)
            IF (aln > elim) GO TO 30
          ELSE
            CALL cuoik(z, fnu, kode, 2, nn, cy, nuf, tol, elim, alim)
            IF (nuf < 0) GO TO 30
            nz = nz + nuf
            nn = nn - nuf
!-----------------------------------------------------------------------
!     HERE NN=N OR NN=0 SINCE NUF=0,NN, OR -1 ON RETURN FROM CUOIK
!     IF NUF=NN, THEN CY(I)=CZERO FOR ALL I
!-----------------------------------------------------------------------
            IF (nn == 0) GO TO 20
          END IF
        END IF

        10 IF (xx >= 0.0_dp) THEN
!-----------------------------------------------------------------------
!     RIGHT HALF PLANE COMPUTATION, REAL(Z) >= 0.
!-----------------------------------------------------------------------
          CALL cbknu(z, fnu, kode, nn, cy, nw, tol, elim, alim)
          IF (nw < 0) GO TO 40
          nz = nw
          RETURN
        END IF
!-----------------------------------------------------------------------
!     LEFT HALF PLANE COMPUTATION
!     PI/2 < ARG(Z) <= PI AND -PI < ARG(Z) < -PI/2.
!-----------------------------------------------------------------------
        IF (nz /= 0) GO TO 30
        mr = 1
        IF (yy < 0.0_dp) mr = -1
        CALL cacon(z, fnu, kode, mr, nn, cy, nw, rl, fnul, tol, elim, alim)
        IF (nw < 0) GO TO 40
        nz = nw
        RETURN
      END IF
!-----------------------------------------------------------------------
!     UNIFORM ASYMPTOTIC EXPANSIONS FOR FNU > FNUL
!-----------------------------------------------------------------------
      mr = 0
      IF (xx < 0.0_dp) THEN
        mr = 1
        IF (yy < 0.0_dp) mr = -1
      END IF
      CALL cbunk(z, fnu, kode, mr, nn, cy, nw, tol, elim, alim)
      IF (nw < 0) GO TO 40
      nz = nz + nw
      RETURN

      20 IF (xx >= 0.0_dp) RETURN
    END IF

    30 nz = 0
    ierr = 2
    RETURN

    40 IF (nw == -1) GO TO 30
    nz = 0
    ierr = 5
    RETURN
  END IF
END IF
nz = 0
ierr = 4
RETURN
END SUBROUTINE cbesk



SUBROUTINE cbesy(z, fnu, kode, n, cy, nz, ierr)

! N.B. Argument CWRK has been removed.

!***BEGIN PROLOGUE  CBESY
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  890801, 930101  (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  Y-BESSEL FUNCTION,BESSEL FUNCTION OF COMPLEX ARGUMENT,
!             BESSEL FUNCTION OF SECOND KIND
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE THE Y-BESSEL FUNCTION OF A COMPLEX ARGUMENT
!***DESCRIPTION

!   ON KODE=1, CBESY COMPUTES AN N MEMBER SEQUENCE OF COMPLEX
!   BESSEL FUNCTIONS CY(I)=Y(FNU+I-1,Z) FOR REAL (dp), NONNEGATIVE
!   ORDERS FNU+I-1, I=1,...,N AND COMPLEX Z IN THE CUT PLANE
!   -PI < ARG(Z) <= PI.
!   ON KODE=2, CBESY RETURNS THE SCALED FUNCTIONS

!   CY(I)=EXP(-ABS(Y))*Y(FNU+I-1,Z)   I = 1,...,N , Y=AIMAG(Z)

!   WHICH REMOVE THE EXPONENTIAL GROWTH IN BOTH THE UPPER AND
!   LOWER HALF PLANES FOR Z TO INFINITY. DEFINITIONS AND NOTATION
!   ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL FUNCTIONS (REF. 1).

!   INPUT
!     Z      - Z=CMPLX(X,Y), Z.NE.CMPLX(0.,0.),-PI < ARG(Z) <= PI
!     FNU    - ORDER OF INITIAL Y FUNCTION, FNU >= 0.0_dp
!     KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!              KODE= 1  RETURNS
!                       CY(I)=Y(FNU+I-1,Z), I=1,...,N
!                  = 2  RETURNS
!                       CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y)), I=1,...,N
!                       WHERE Y=AIMAG(Z)
!     N      - NUMBER OF MEMBERS OF THE SEQUENCE, N >= 1
!     CWRK   - A COMPLEX WORK VECTOR OF DIMENSION AT LEAST N

!   OUTPUT
!     CY     - A COMPLEX VECTOR WHOSE FIRST N COMPONENTS CONTAIN
!              VALUES FOR THE SEQUENCE
!              CY(I)=Y(FNU+I-1,Z)  OR
!              CY(I)=Y(FNU+I-1,Z)*EXP(-ABS(Y))  I=1,...,N
!              DEPENDING ON KODE.
!     NZ     - NZ=0 , A NORMAL RETURN
!              NZ > 0 , NZ COMPONENTS OF CY SET TO ZERO DUE TO
!              UNDERFLOW (GENERALLY ON KODE=2)
!     IERR   - ERROR FLAG
!              IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!              IERR=1, INPUT ERROR   - NO COMPUTATION
!              IERR=2, OVERFLOW      - NO COMPUTATION, FNU+N-1 IS
!                      TOO LARGE OR ABS(Z) IS TOO SMALL OR BOTH
!              IERR=3, ABS(Z) OR FNU+N-1 LARGE - COMPUTATION DONE
!                      BUT LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!                      PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!              IERR=4, ABS(Z) OR FNU+N-1 TOO LARGE - NO COMPUTATION
!                      BECAUSE OF COMPLETE LOSSES OF SIGNIFICANCE
!                      BY ARGUMENT REDUCTION
!              IERR=5, ERROR              - NO COMPUTATION,
!                      ALGORITHM TERMINATION CONDITION NOT MET

!***LONG DESCRIPTION

!   THE COMPUTATION IS CARRIED OUT IN TERMS OF THE I(FNU,Z) AND
!   K(FNU,Z) BESSEL FUNCTIONS IN THE RIGHT HALF PLANE BY

!       Y(FNU,Z) = I*CC*I(FNU,ARG) - (2/PI)*CONJG(CC)*K(FNU,ARG)

!       Y(FNU,Z) = CONJG(Y(FNU,CONJG(Z)))

!   FOR AIMAG(Z) >= 0 AND AIMAG(Z) < 0 RESPECTIVELY, WHERE
!   CC=EXP(I*PI*FNU/2), ARG=Z*EXP(-I*PI/2) AND I**2=-1.

!   FOR NEGATIVE ORDERS,THE FORMULA

!       Y(-FNU,Z) = Y(FNU,Z)*COS(PI*FNU) + J(FNU,Z)*SIN(PI*FNU)

!   CAN BE USED.  HOWEVER,FOR LARGE ORDERS CLOSE TO HALF ODD INTEGERS THE
!   FUNCTION CHANGES RADICALLY.  WHEN FNU IS A LARGE POSITIVE HALF ODD INTEGER,
!   THE MAGNITUDE OF Y(-FNU,Z) = J(FNU,Z)*SIN(PI*FNU) IS A LARGE NEGATIVE
!   POWER OF TEN.  BUT WHEN FNU IS NOT A HALF ODD INTEGER, Y(FNU,Z) DOMINATES
!   IN MAGNITUDE WITH A LARGE POSITIVE POWER OF TEN AND THE MOST THAT THE
!   SECOND TERM CAN BE REDUCED IS BY UNIT ROUNDOFF FROM THE COEFFICIENT.
!   THUS, WIDE CHANGES CAN OCCUR WITHIN UNIT ROUNDOFF OF A LARGE HALF
!   ODD INTEGER.  HERE, LARGE MEANS FNU > ABS(Z).

!   IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELEMENTARY
!   FUNCTIONS.  WHEN THE MAGNITUDE OF Z OR FNU+N-1 IS LARGE, LOSSES OF
!   SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR.  CONSEQUENTLY, IF EITHER ONE
!   EXCEEDS U1=SQRT(0.5/UR), THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY
!   AND AN ERROR FLAG IERR=3 IS TRIGGERED WHERE UR = EPSILON(0.0_dp) = UNIT
!   ROUNDOFF.  ALSO IF EITHER IS LARGER THAN U2=0.5/UR, THEN ALL SIGNIFICANCE
!   IS LOST AND IERR=4.  IN ORDER TO USE THE INT FUNCTION, ARGUMENTS MUST BE
!   FURTHER RESTRICTED NOT TO EXCEED THE LARGEST MACHINE INTEGER, U3 = HUGE(0).
!   THUS, THE MAGNITUDE OF Z AND FNU+N-1 IS RESTRICTED BY MIN(U2,U3).
!   ON 32 BIT MACHINES, U1,U2, AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9
!   IN SINGLE PRECISION ARITHMETIC AND 1.3E+8, 1.8D+16, 2.1E+9 IN DOUBLE
!   PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMITING IN
!   THEIR RESPECTIVE ARITHMETICS.  THIS MEANS THAT ONE CAN EXPECT TO RETAIN,
!   IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS IN SINGLE AND ONLY
!   7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!   SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.

!   THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX BESSEL
!   FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P = MAX(UNIT ROUNDOFF,1.0E-18)
!   IS THE NOMINAL PRECISION AND 10**S REPRESENTS THE INCREASE IN ERROR DUE TO
!   ARGUMENT REDUCTION IN THE ELEMENTARY FUNCTIONS.  HERE, S =
!   MAX(1,ABS(LOG10(ABS(Z))), ABS(LOG10(FNU))) APPROXIMATELY (I.E. S =
!   MAX(1,ABS(EXPONENT OF ABS(Z),ABS(EXPONENT OF FNU)) ).
!   HOWEVER, THE PHASE ANGLE MAY HAVE ONLY ABSOLUTE ACCURACY.  THIS IS MOST
!   LIKELY TO OCCUR WHEN ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE
!   OTHER BY SEVERAL ORDERS OF MAGNITUDE.  IF ONE COMPONENT IS 10**K LARGER
!   THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K, 0)
!   SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS THE EXPONENT
!   OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER COMPONENT.
!   HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY BECAUSE, IN COMPLEX
!   ARITHMETIC WITH PRECISION P, THE SMALLER COMPONENT WILL NOT (AS A RULE)
!   DECREASE BELOW P TIMES THE MAGNITUDE OF THE LARGER COMPONENT.  IN THESE
!   EXTREME CASES, THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P,
!   PI/2-P, OR -PI/2+P.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!        I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!      COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!        BY D. E. AMOS, SAND83-0083, MAY 1983.

!      COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!        AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY 1983

!      A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!        AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-1018, MAY 1985

!      A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!        AND NONNEGATIVE ORDER BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE,
!        VOL. 12, NO. 3, SEPTEMBER 1986, PP 265-273.

!***ROUTINES CALLED  CBESI,CBESK,I1MACH,R1MACH
!***END PROLOGUE  CBESY

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: cy(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: ierr

COMPLEX (dp)  :: ci, csgn, cspn, cwrk(n), ex, zu, zv, zz, zn
REAL (dp)     :: arg, elim, ey, r1, r2, tay, xx, yy, ascle, rtol,  &
                 atol, tol, aa, bb, ffnu, rhpi, r1m5
INTEGER       :: i, ifnu, k, k1, k2, nz1, nz2, i4
COMPLEX (dp), PARAMETER  :: cip(4) = (/ (1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp),  &
                                        (-1.0_dp, 0.0_dp), (0.0_dp, -1.0_dp) /)
REAL (dp), PARAMETER     :: hpi = 1.57079632679489662_dp

!***FIRST EXECUTABLE STATEMENT  CBESY
xx = REAL(z, KIND=dp)
yy = AIMAG(z)
ierr = 0
nz = 0
IF (xx == 0.0_dp .AND. yy == 0.0_dp) ierr = 1
IF (fnu < 0.0_dp) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (n < 1) ierr = 1
IF (ierr /= 0) RETURN
ci = CMPLX(0.0_dp, 1.0_dp, KIND=dp)
zz = z
IF (yy < 0.0_dp) zz = CONJG(z)
zn = -ci * zz
CALL cbesi(zn, fnu, kode, n, cy, nz1, ierr)
IF (ierr == 0 .OR. ierr == 3) THEN
  CALL cbesk(zn, fnu, kode, n, cwrk, nz2, ierr)
  IF (ierr == 0 .OR. ierr == 3) THEN
    nz = MIN(nz1, nz2)
    ifnu = fnu
    ffnu = fnu - ifnu
    arg = hpi * ffnu
    csgn = CMPLX(COS(arg), SIN(arg), KIND=dp)
    i4 = MOD(ifnu, 4) + 1
    csgn = csgn * cip(i4)
    rhpi = 1.0_dp / hpi
    cspn = CONJG(csgn) * rhpi
    csgn = csgn * ci
    IF (kode /= 2) THEN
      DO  i = 1, n
        cy(i) = csgn * cy(i) - cspn * cwrk(i)
        csgn = ci * csgn
        cspn = -ci * cspn
      END DO
      IF (yy < 0.0_dp) cy(1:n) = CONJG(cy(1:n))
      RETURN
    END IF

    r1 = COS(xx)
    r2 = SIN(xx)
    ex = CMPLX(r1, r2, KIND=dp)
    tol = MAX(EPSILON(0.0_dp), 1.0D-18)
    k1 = MINEXPONENT(0.0_dp)
    k2 = MAXEXPONENT(0.0_dp)
    k = MIN(ABS(k1),ABS(k2))
    r1m5 = LOG10( REAL( RADIX(0.0_dp), KIND=dp) )
!-----------------------------------------------------------------------
!     ELIM IS THE APPROXIMATE EXPONENTIAL UNDER- AND OVERFLOW LIMIT
!-----------------------------------------------------------------------
    elim = 2.303_dp * (k*r1m5 - 3.0_dp)
    ey = 0.0_dp
    tay = ABS(yy+yy)
    IF (tay < elim) ey = EXP(-tay)
    cspn = ex * ey * cspn
    nz = 0
    rtol = 1.0_dp / tol
    ascle = TINY(0.0_dp) * rtol * 1.0E+3
    DO  i = 1, n
!----------------------------------------------------------------------
!       CY(I) = CSGN*CY(I)-CSPN*CWRK(I): PRODUCTS ARE COMPUTED IN
!       SCALED MODE IF CY(I) OR CWRK(I) ARE CLOSE TO UNDERFLOW TO
!       PREVENT UNDERFLOW IN AN INTERMEDIATE COMPUTATION.
!----------------------------------------------------------------------
      zv = cwrk(i)
      aa = REAL(zv, KIND=dp)
      bb = AIMAG(zv)
      atol = 1.0_dp
      IF (MAX(ABS(aa),ABS(bb)) <= ascle) THEN
        zv = zv * rtol
        atol = tol
      END IF
      zv = zv * cspn
      zv = zv * atol
      zu = cy(i)
      aa = REAL(zu, KIND=dp)
      bb = AIMAG(zu)
      atol = 1.0_dp
      IF (MAX(ABS(aa),ABS(bb)) <= ascle) THEN
        zu = zu * rtol
        atol = tol
      END IF
      zu = zu * csgn
      zu = zu * atol
      cy(i) = zu - zv
      IF (yy < 0.0_dp) cy(i) = CONJG(cy(i))
      IF (cy(i) == CMPLX(0.0_dp, 0.0_dp, KIND=dp) .AND. ey == 0.0_dp) nz = nz + 1
      csgn = ci * csgn
      cspn = -ci * cspn
    END DO
    RETURN
  END IF
END IF
nz = 0
RETURN
END SUBROUTINE cbesy



SUBROUTINE cairy(z, id, kode, ai, nz, ierr)
!***BEGIN PROLOGUE  CAIRY
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  890801, 930101   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE AIRY FUNCTIONS AI(Z) AND DAI(Z) FOR COMPLEX Z
!***DESCRIPTION

!   ON KODE=1, CAIRY COMPUTES THE COMPLEX AIRY FUNCTION AI(Z) OR
!   ITS DERIVATIVE DAI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY.  ON
!   KODE=2, A SCALING OPTION EXP(ZTA)*AI(Z) OR EXP(ZTA)*
!   DAI(Z)/DZ IS PROVIDED TO REMOVE THE EXPONENTIAL DECAY IN
!   -PI/3 < ARG(Z) < PI/3 AND THE EXPONENTIAL GROWTH IN
!   PI/3 < ABS(ARG(Z)) < PI WHERE ZTA=(2/3)*Z*SQRT(Z)

!   WHILE THE AIRY FUNCTIONS AI(Z) AND DAI(Z)/DZ ARE ANALYTIC IN
!   THE WHOLE Z PLANE, THE CORRESPONDING SCALED FUNCTIONS DEFINED
!   FOR KODE=2 HAVE A CUT ALONG THE NEGATIVE REAL AXIS.
!   DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF
!   MATHEMATICAL FUNCTIONS (REF. 1).

!   INPUT
!     Z      - Z=CMPLX(X,Y)
!     ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
!     KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!              KODE= 1  RETURNS
!                       AI=AI(Z)                ON ID=0 OR
!                       AI=DAI(Z)/DZ            ON ID=1
!                  = 2  RETURNS
!                       AI=EXP(ZTA)*AI(Z)       ON ID=0 OR
!                       AI=EXP(ZTA)*DAI(Z)/DZ   ON ID=1 WHERE
!                       ZTA=(2/3)*Z*SQRT(Z)

!   OUTPUT
!     AI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND KODE
!     NZ     - UNDERFLOW INDICATOR
!              NZ= 0   , NORMAL RETURN
!              NZ= 1   , AI=CMPLX(0.0,0.0) DUE TO UNDERFLOW IN
!                        -PI/3 < ARG(Z) < PI/3 ON KODE=1
!     IERR   - ERROR FLAG
!              IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!              IERR=1, INPUT ERROR   - NO COMPUTATION
!              IERR=2, OVERFLOW      - NO COMPUTATION, REAL(ZTA)
!                      TOO LARGE WITH KODE=1.
!              IERR=3, ABS(Z) LARGE      - COMPUTATION COMPLETED
!                      LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!                      PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!              IERR=4, ABS(Z) TOO LARGE  - NO COMPUTATION
!                      COMPLETE LOSS OF ACCURACY BY ARGUMENT REDUCTION
!              IERR=5, ERROR              - NO COMPUTATION,
!                      ALGORITHM TERMINATION CONDITION NOT MET


!***LONG DESCRIPTION

!   AI AND DAI ARE COMPUTED FOR ABS(Z) > 1.0 FROM THE K BESSEL FUNCTIONS BY
!      AI(Z) = C*SQRT(Z)*K(1/3,ZTA) , DAI(Z) = -C*Z*K(2/3,ZTA)
!                     C = 1.0/(PI*SQRT(3.0))
!                   ZTA = (2/3)*Z**(3/2)

!   WITH THE POWER SERIES FOR ABS(Z) <= 1.0.

!   IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELEMENTARY
!   FUNCTIONS.  WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES OF SIGNIFICANCE BY
!   ARGUMENT REDUCTION OCCUR.  CONSEQUENTLY, IF THE MAGNITUDE OF ZETA =
!   (2/3)*Z**1.5 EXCEEDS U1 = SQRT(0.5/UR), THEN LOSSES EXCEEDING HALF
!   PRECISION ARE LIKELY AND AN ERROR FLAG IERR=3 IS TRIGGERED WHERE UR =
!   EPSILON(0.0_dp) = UNIT ROUNDOFF.  ALSO, IF THE MAGNITUDE OF ZETA IS LARGER
!   THAN U2=0.5/UR, THEN ALL SIGNIFICANCE IS LOST AND IERR=4.  IN ORDER TO USE
!   THE INT FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
!   LARGEST INTEGER, U3 = HUGE(0).  THUS, THE MAGNITUDE OF ZETA MUST BE
!   RESTRICTED BY MIN(U2,U3).  ON 32 BIT MACHINES, U1, U2, AND U3 ARE
!   APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE PRECISION ARITHMETIC AND
!   1.3E+8, 1.8D+16, 2.1E+9 IN DOUBLE PRECISION ARITHMETIC RESPECTIVELY.
!   THIS MAKES U2 AND U3 LIMITING IN THEIR RESPECTIVE ARITHMETICS.  THIS MEANS
!   THAT THE MAGNITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
!   DOUBLE PRECISION ARITHMETIC.  THIS ALSO MEANS THAT ONE CAN EXPECT TO
!   RETAIN, IN THE WORST CASES ON 32 BIT MACHINES, NO DIGITS IN SINGLE
!   PRECISION AND ONLY 7 DIGITS IN DOUBLE PRECISION ARITHMETIC.
!   SIMILAR CONSIDERATIONS HOLD FOR OTHER MACHINES.

!   THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!   BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P = MAX(UNIT
!   ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRESENTS
!   THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!   ELEMENTARY FUNCTIONS.  HERE, S = MAX(1,ABS(LOG10(ABS(Z))),
!   ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
!   ABS(Z),ABS(EXPONENT OF FNU)) ).  HOWEVER, THE PHASE ANGLE MAY
!   HAVE ONLY ABSOLUTE ACCURACY.  THIS IS MOST LIKELY TO OCCUR WHEN
!   ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!   SEVERAL ORDERS OF MAGNITUDE.  IF ONE COMPONENT IS 10**K LARGER
!   THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
!   0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
!   THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!   COMPONENT.  HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
!   BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!   COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
!   MAGNITUDE OF THE LARGER COMPONENT.  IN THESE EXTREME CASES,
!   THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P, OR -PI/2+P.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!      I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!    COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!      AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY 1983

!    A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!      AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-1018, MAY 1985

!    A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!      AND NONNEGATIVE ORDER BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE,
!      VOL. 12, NO. 3, SEPTEMBER 1986, PP 265-273.

!***ROUTINES CALLED  CACAI,CBKNU,I1MACH,R1MACH
!***END PROLOGUE  CAIRY

COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(IN)        :: id
INTEGER, INTENT(IN)        :: kode
COMPLEX (dp), INTENT(OUT)  :: ai
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: ierr

COMPLEX (dp)  :: csq, cy(1), s1, s2, trm1, trm2, zta, z3
REAL (dp)     :: aa, ad, ak, alim, atrm, az, az3, bk, ck, dig, dk, d1, d2, &
                 elim, fid, fnu, rl, r1m5, sfac, tol, zi, zr, z3i, z3r, bb, &
                 alaz
INTEGER       :: iflag, k, k1, k2, mr, nn
REAL (dp), PARAMETER  :: tth = 6.66666666666666667D-01,  &
           c1 = 3.55028053887817240D-01, c2 = 2.58819403792806799D-01,  &
           coef = 1.83776298473930683D-01
COMPLEX (dp), PARAMETER  :: cone = (1.0_dp, 0.0_dp)

!***FIRST EXECUTABLE STATEMENT  CAIRY
ierr = 0
nz = 0
IF (id < 0 .OR. id > 1) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (ierr /= 0) RETURN
az = ABS(z)
tol = MAX(EPSILON(0.0_dp), 1.0D-18)
fid = id
IF (az <= 1.0_dp) THEN
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(Z) <= 1.
!-----------------------------------------------------------------------
  s1 = cone
  s2 = cone
  IF (az < tol) GO TO 30
  aa = az * az
  IF (aa >= tol/az) THEN
    trm1 = cone
    trm2 = cone
    atrm = 1.0_dp
    z3 = z * z * z
    az3 = az * aa
    ak = 2.0_dp + fid
    bk = 3.0_dp - fid - fid
    ck = 4.0_dp - fid
    dk = 3.0_dp + fid + fid
    d1 = ak * dk
    d2 = bk * ck
    ad = MIN(d1,d2)
    ak = 24.0_dp + 9.0_dp * fid
    bk = 30.0_dp - 9.0_dp * fid
    z3r = REAL(z3, KIND=dp)
    z3i = AIMAG(z3)
    DO  k = 1, 25
      trm1 = trm1 * CMPLX(z3r/d1, z3i/d1, KIND=dp)
      s1 = s1 + trm1
      trm2 = trm2 * CMPLX(z3r/d2, z3i/d2, KIND=dp)
      s2 = s2 + trm2
      atrm = atrm * az3 / ad
      d1 = d1 + ak
      d2 = d2 + bk
      ad = MIN(d1,d2)
      IF (atrm < tol*ad) EXIT
      ak = ak + 18.0_dp
      bk = bk + 18.0_dp
    END DO
  END IF

  IF (id /= 1) THEN
    ai = s1 * c1 - z * s2 * c2
    IF (kode == 1) RETURN
    zta = z * SQRT(z) * tth
    ai = ai * EXP(zta)
    RETURN
  END IF
  ai = -s2 * c2
  IF (az > tol) ai = ai + z * z * s1 * c1/(1.0_dp + fid)
  IF (kode == 1) RETURN
  zta = z * SQRT(z) * tth
  ai = ai * EXP(zta)
  RETURN
END IF
!-----------------------------------------------------------------------
!     CASE FOR ABS(Z) > 1.0
!-----------------------------------------------------------------------
fnu = (1.0_dp + fid) / 3.0_dp
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!-----------------------------------------------------------------------
k1 = MINEXPONENT(0.0_dp)
k2 = MAXEXPONENT(0.0_dp)
r1m5 = LOG10( REAL( RADIX(0.0_dp), KIND=dp) )
k = MIN(ABS(k1),ABS(k2))
elim = 2.303_dp * (k*r1m5 - 3.0_dp)
k1 = DIGITS(0.0_dp) - 1
aa = r1m5 * k1
dig = MIN(aa,18.0_dp)
aa = aa * 2.303_dp
alim = elim + MAX(-aa,-41.45_dp)
rl = 1.2_dp * dig + 3.0_dp
alaz = LOG(az)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
aa = 0.5_dp / tol
bb = HUGE(0) * 0.5_dp
aa = MIN(aa,bb)
aa = aa ** tth
IF (az > aa) GO TO 70
aa = SQRT(aa)
IF (az > aa) ierr = 3
csq = SQRT(z)
zta = z * csq * tth
!-----------------------------------------------------------------------
!     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
!-----------------------------------------------------------------------
iflag = 0
sfac = 1.0_dp
zi = AIMAG(z)
zr = REAL(z, KIND=dp)
ak = AIMAG(zta)
IF (zr < 0.0_dp) THEN
  bk = REAL(zta, KIND=dp)
  ck = -ABS(bk)
  zta = CMPLX(ck, ak, KIND=dp)
END IF
IF (zi == 0.0_dp) THEN
  IF (zr <= 0.0_dp) THEN
    zta = CMPLX(0.0_dp, ak, KIND=dp)
  END IF
END IF
aa = REAL(zta, KIND=dp)
IF (aa < 0.0_dp .OR. zr <= 0.0_dp) THEN
  IF (kode /= 2) THEN
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
    IF (aa <= -alim) THEN
      aa = -aa + 0.25_dp * alaz
      iflag = 1
      sfac = tol
      IF (aa > elim) GO TO 50
    END IF
  END IF
!-----------------------------------------------------------------------
!     CBKNU AND CACAI RETURN EXP(ZTA)*K(FNU,ZTA) ON KODE=2
!-----------------------------------------------------------------------
  mr = 1
  IF (zi < 0.0_dp) mr = -1
  CALL cacai(zta, fnu, kode, mr, 1, cy, nn, rl, tol, elim, alim)
  IF (nn < 0) GO TO 60
  nz = nz + nn
ELSE
  IF (kode /= 2) THEN
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
    IF (aa >= alim) THEN
      aa = -aa - 0.25_dp * alaz
      iflag = 2
      sfac = 1.0_dp / tol
      IF (aa < -elim) GO TO 40
    END IF
  END IF
  CALL cbknu(zta, fnu, kode, 1, cy, nz, tol, elim, alim)
END IF
s1 = cy(1) * coef
IF (iflag == 0) THEN
  IF (id /= 1) THEN
    ai = csq * s1
    RETURN
  END IF
  ai = -z * s1
  RETURN
END IF
s1 = s1 * sfac
IF (id /= 1) THEN
  s1 = s1 * csq
  ai = s1 / sfac
  RETURN
END IF
s1 = -s1 * z
ai = s1 / sfac
RETURN

30 aa = 1.0E+3 * TINY(0.0_dp)
s1 = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
IF (id /= 1) THEN
  IF (az > aa) s1 = c2 * z
  ai = c1 - s1
  RETURN
END IF
ai = -c2
aa = SQRT(aa)
IF (az > aa) s1 = z * z * 0.5_dp
ai = ai + s1 * c1
RETURN

40 nz = 1
ai = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
RETURN

50 nz = 0
ierr = 2
RETURN

60 IF (nn == -1) GO TO 50
nz = 0
ierr = 5
RETURN

70 ierr = 4
nz = 0
RETURN
END SUBROUTINE cairy



SUBROUTINE cbiry(z, id, kode, bi, ierr)
!***BEGIN PROLOGUE  CBIRY
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  890801, 930101   (YYMMDD)
!***CATEGORY NO.  B5K
!***KEYWORDS  AIRY FUNCTION,BESSEL FUNCTIONS OF ORDER ONE THIRD
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE AIRY FUNCTIONS BI(Z) AND DBI(Z) FOR COMPLEX Z
!***DESCRIPTION

!   ON KODE=1, CBIRY COMPUTES THE COMPLEX AIRY FUNCTION BI(Z) OR ITS
!   DERIVATIVE DBI(Z)/DZ ON ID=0 OR ID=1 RESPECTIVELY.  ON KODE=2,
!   A SCALING OPTION EXP(-AXZTA)*BI(Z) OR EXP(-AXZTA)*DBI(Z)/DZ
!   IS PROVIDED TO REMOVE THE EXPONENTIAL BEHAVIOR IN BOTH THE LEFT AND
!   RIGHT HALF PLANES WHERE ZTA = (2/3)*Z*SQRT(Z) = CMPLX(XZTA,YZTA)
!   AND AXZTA=ABS(XZTA).
!   DEFINITIONS AND NOTATION ARE FOUND IN THE NBS HANDBOOK OF MATHEMATICAL
!   FUNCTIONS (REF. 1).

!   INPUT
!     Z      - Z=CMPLX(X,Y)
!     ID     - ORDER OF DERIVATIVE, ID=0 OR ID=1
!     KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
!              KODE= 1  RETURNS
!                       BI=BI(Z)                 ON ID=0 OR
!                       BI=DBI(Z)/DZ             ON ID=1
!                  = 2  RETURNS
!                       BI=EXP(-AXZTA)*BI(Z)     ON ID=0 OR
!                       BI=EXP(-AXZTA)*DBI(Z)/DZ ON ID=1 WHERE
!                       ZTA=(2/3)*Z*SQRT(Z)=CMPLX(XZTA,YZTA)
!                       AND AXZTA=ABS(XZTA)

!   OUTPUT
!     BI     - COMPLEX ANSWER DEPENDING ON THE CHOICES FOR ID AND KODE
!     IERR   - ERROR FLAG
!              IERR=0, NORMAL RETURN - COMPUTATION COMPLETED
!              IERR=1, INPUT ERROR   - NO COMPUTATION
!              IERR=2, OVERFLOW      - NO COMPUTATION, REAL(Z)
!                      TOO LARGE WITH KODE=1
!              IERR=3, ABS(Z) LARGE      - COMPUTATION COMPLETED
!                      LOSSES OF SIGNIFICANCE BY ARGUMENT REDUCTION
!                      PRODUCE LESS THAN HALF OF MACHINE ACCURACY
!              IERR=4, ABS(Z) TOO LARGE  - NO COMPUTATION
!                      COMPLETE LOSS OF ACCURACY BY ARGUMENT
!                      REDUCTION
!              IERR=5, ERROR              - NO COMPUTATION,
!                      ALGORITHM TERMINATION CONDITION NOT MET

!***LONG DESCRIPTION

!      BI AND DBI ARE COMPUTED FOR ABS(Z) > 1.0 FROM THE I BESSEL
!      FUNCTIONS BY

!             BI(Z)=C*SQRT(Z)*( I(-1/3,ZTA) + I(1/3,ZTA) )
!            DBI(Z)=C *  Z  * ( I(-2/3,ZTA) + I(2/3,ZTA) )
!                            C=1.0/SQRT(3.0)
!                            ZTA=(2/3)*Z**(3/2)

!      WITH THE POWER SERIES FOR ABS(Z) <= 1.0.

!      IN MOST COMPLEX VARIABLE COMPUTATION, ONE MUST EVALUATE ELE-
!      MENTARY FUNCTIONS. WHEN THE MAGNITUDE OF Z IS LARGE, LOSSES
!      OF SIGNIFICANCE BY ARGUMENT REDUCTION OCCUR. CONSEQUENTLY, IF
!      THE MAGNITUDE OF ZETA=(2/3)*Z**1.5 EXCEEDS U1=SQRT(0.5/UR),
!      THEN LOSSES EXCEEDING HALF PRECISION ARE LIKELY AND AN ERROR
!      FLAG IERR=3 IS TRIGGERED WHERE UR=EPSILON(0.0_dp)=UNIT ROUNDOFF.
!      ALSO, IF THE MAGNITUDE OF ZETA IS LARGER THAN U2=0.5/UR, THEN
!      ALL SIGNIFICANCE IS LOST AND IERR=4. IN ORDER TO USE THE INT
!      FUNCTION, ZETA MUST BE FURTHER RESTRICTED NOT TO EXCEED THE
!      LARGEST INTEGER, U3=HUGE(0). THUS, THE MAGNITUDE OF ZETA
!      MUST BE RESTRICTED BY MIN(U2,U3). ON 32 BIT MACHINES, U1,U2,
!      AND U3 ARE APPROXIMATELY 2.0E+3, 4.2E+6, 2.1E+9 IN SINGLE
!      PRECISION ARITHMETIC AND 1.3E+8, 1.8D+16, 2.1E+9 IN DOUBLE
!      PRECISION ARITHMETIC RESPECTIVELY. THIS MAKES U2 AND U3 LIMIT-
!      ING IN THEIR RESPECTIVE ARITHMETICS. THIS MEANS THAT THE MAG-
!      NITUDE OF Z CANNOT EXCEED 3.1E+4 IN SINGLE AND 2.1E+6 IN
!      DOUBLE PRECISION ARITHMETIC. THIS ALSO MEANS THAT ONE CAN
!      EXPECT TO RETAIN, IN THE WORST CASES ON 32 BIT MACHINES,
!      NO DIGITS IN SINGLE PRECISION AND ONLY 7 DIGITS IN DOUBLE
!      PRECISION ARITHMETIC.

!      THE APPROXIMATE RELATIVE ERROR IN THE MAGNITUDE OF A COMPLEX
!      BESSEL FUNCTION CAN BE EXPRESSED BY P*10**S WHERE P=MAX(UNIT
!      ROUNDOFF,1.0E-18) IS THE NOMINAL PRECISION AND 10**S REPRE-
!      SENTS THE INCREASE IN ERROR DUE TO ARGUMENT REDUCTION IN THE
!      ELEMENTARY FUNCTIONS. HERE, S=MAX(1,ABS(LOG10(ABS(Z))),
!      ABS(LOG10(FNU))) APPROXIMATELY (I.E. S=MAX(1,ABS(EXPONENT OF
!      ABS(Z),ABS(EXPONENT OF FNU)) ). HOWEVER, THE PHASE ANGLE MAY
!      HAVE ONLY ABSOLUTE ACCURACY. THIS IS MOST LIKELY TO OCCUR WHEN
!      ONE COMPONENT (IN ABSOLUTE VALUE) IS LARGER THAN THE OTHER BY
!      SEVERAL ORDERS OF MAGNITUDE. IF ONE COMPONENT IS 10**K LARGER
!      THAN THE OTHER, THEN ONE CAN EXPECT ONLY MAX(ABS(LOG10(P))-K,
!      0) SIGNIFICANT DIGITS; OR, STATED ANOTHER WAY, WHEN K EXCEEDS
!      THE EXPONENT OF P, NO SIGNIFICANT DIGITS REMAIN IN THE SMALLER
!      COMPONENT. HOWEVER, THE PHASE ANGLE RETAINS ABSOLUTE ACCURACY
!      BECAUSE, IN COMPLEX ARITHMETIC WITH PRECISION P, THE SMALLER
!      COMPONENT WILL NOT (AS A RULE) DECREASE BELOW P TIMES THE
!      MAGNITUDE OF THE LARGER COMPONENT. IN THESE EXTREME CASES,
!      THE PRINCIPAL PHASE ANGLE IS ON THE ORDER OF +P, -P, PI/2-P,
!      OR -PI/2+P.

!***REFERENCES  HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND
!           I. A. STEGUN, NBS AMS SERIES 55, U.S. DEPT. OF COMMERCE, 1955.

!         COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!           AND LARGE ORDER BY D. E. AMOS, SAND83-0643, MAY 1983

!         A SUBROUTINE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!           AND NONNEGATIVE ORDER BY D. E. AMOS, SAND85-1018, MAY 1985

!         A PORTABLE PACKAGE FOR BESSEL FUNCTIONS OF A COMPLEX ARGUMENT
!           AND NONNEGATIVE ORDER BY D. E. AMOS, ACM TRANS. MATH. SOFTWARE,
!           VOL. 12, NO. 3, SEPTEMBER 1986, PP 265-273.

!***ROUTINES CALLED  CBINU,I1MACH,R1MACH
!***END PROLOGUE  CBIRY

COMPLEX (dp), INTENT(IN)   :: z
INTEGER, INTENT(IN)        :: id
INTEGER, INTENT(IN)        :: kode
COMPLEX (dp), INTENT(OUT)  :: bi
INTEGER, INTENT(OUT)       :: ierr

COMPLEX (dp)  :: csq, cy(2), s1, s2, trm1, trm2, zta, z3
REAL (dp)     :: aa, ad, ak, alim, atrm, az, az3, bb, bk, ck, dig, dk,  &
                 d1, d2, elim, fid, fmr, fnu, fnul, rl, r1m5,  &
                 sfac, tol, zi, zr, z3i, z3r
INTEGER       :: k, k1, k2, nz
REAL (dp), PARAMETER  :: tth = 6.66666666666666667D-01,  &
    c1 = 6.14926627446000736D-01, c2 = 4.48288357353826359D-01,  &
    coef = 5.77350269189625765D-01, pi = 3.141592653589793238_dp
REAL (dp), PARAMETER  :: cone = (1.0_dp,0.0_dp)

!***FIRST EXECUTABLE STATEMENT  CBIRY
ierr = 0
nz = 0
IF (id < 0 .OR. id > 1) ierr = 1
IF (kode < 1 .OR. kode > 2) ierr = 1
IF (ierr /= 0) RETURN
az = ABS(z)
tol = MAX(EPSILON(0.0_dp), 1.0D-18)
fid = id
IF (az <= 1.0_dp) THEN
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(Z) <= 1.
!-----------------------------------------------------------------------
  s1 = cone
  s2 = cone
  IF (az < tol) GO TO 30
  aa = az * az
  IF (aa >= tol/az) THEN
    trm1 = cone
    trm2 = cone
    atrm = 1.0_dp
    z3 = z * z * z
    az3 = az * aa
    ak = 2.0_dp + fid
    bk = 3.0_dp - fid - fid
    ck = 4.0_dp - fid
    dk = 3.0_dp + fid + fid
    d1 = ak * dk
    d2 = bk * ck
    ad = MIN(d1,d2)
    ak = 24.0_dp + 9.0_dp * fid
    bk = 30.0_dp - 9.0_dp * fid
    z3r = REAL(z3, KIND=dp)
    z3i = AIMAG(z3)
    DO  k = 1, 25
      trm1 = trm1 * CMPLX(z3r/d1, z3i/d1, KIND=dp)
      s1 = s1 + trm1
      trm2 = trm2 * CMPLX(z3r/d2, z3i/d2, KIND=dp)
      s2 = s2 + trm2
      atrm = atrm * az3 / ad
      d1 = d1 + ak
      d2 = d2 + bk
      ad = MIN(d1,d2)
      IF (atrm < tol*ad) EXIT
      ak = ak + 18.0_dp
      bk = bk + 18.0_dp
    END DO
  END IF

  IF (id /= 1) THEN
    bi = s1 * c1 + z * s2 * c2
    IF (kode == 1) RETURN
    zta = z * SQRT(z) * tth
    aa = REAL(zta, KIND=dp)
    aa = -ABS(aa)
    bi = bi * EXP(aa)
    RETURN
  END IF
  bi = s2 * c2
  IF (az > tol) bi = bi + z * z * s1 * c1/(1.0_dp+fid)
  IF (kode == 1) RETURN
  zta = z * SQRT(z) * tth
  aa = REAL(zta, KIND=dp)
  aa = -ABS(aa)
  bi = bi * EXP(aa)
  RETURN
END IF
!-----------------------------------------------------------------------
!     CASE FOR ABS(Z) > 1.0
!-----------------------------------------------------------------------
fnu = (1.0_dp+fid) / 3.0_dp
!-----------------------------------------------------------------------
!     SET PARAMETERS RELATED TO MACHINE CONSTANTS.
!     TOL IS THE APPROXIMATE UNIT ROUNDOFF LIMITED TO 1.0E-18.
!     ELIM IS THE APPROXIMATE EXPONENTIAL OVER- AND UNDERFLOW LIMIT.
!     EXP(-ELIM) < EXP(-ALIM)=EXP(-ELIM)/TOL    AND
!     EXP(ELIM) > EXP(ALIM)=EXP(ELIM)*TOL       ARE INTERVALS NEAR
!     UNDERFLOW AND OVERFLOW LIMITS WHERE SCALED ARITHMETIC IS DONE.
!     RL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC EXPANSION FOR LARGE Z.
!     DIG = NUMBER OF BASE 10 DIGITS IN TOL = 10**(-DIG).
!     FNUL IS THE LOWER BOUNDARY OF THE ASYMPTOTIC SERIES FOR LARGE FNU.
!-----------------------------------------------------------------------
k1 = MINEXPONENT(0.0_dp)
k2 = MAXEXPONENT(0.0_dp)
r1m5 = LOG10( REAL( RADIX(0.0_dp), KIND=dp) )
k = MIN(ABS(k1),ABS(k2))
elim = 2.303_dp * (k*r1m5 - 3.0_dp)
k1 = DIGITS(0.0_dp) - 1
aa = r1m5 * k1
dig = MIN(aa,18.0_dp)
aa = aa * 2.303_dp
alim = elim + MAX(-aa,-41.45_dp)
rl = 1.2_dp * dig + 3.0_dp
fnul = 10.0_dp + 6.0_dp * (dig - 3.0_dp)
!-----------------------------------------------------------------------
!     TEST FOR RANGE
!-----------------------------------------------------------------------
aa = 0.5_dp / tol
bb = HUGE(0) * 0.5_dp
aa = MIN(aa,bb)
aa = aa ** tth
IF (az > aa) GO TO 60
aa = SQRT(aa)
IF (az > aa) ierr = 3
csq = SQRT(z)
zta = z * csq * tth
!-----------------------------------------------------------------------
!     RE(ZTA) <= 0 WHEN RE(Z) < 0, ESPECIALLY WHEN IM(Z) IS SMALL
!-----------------------------------------------------------------------
sfac = 1.0_dp
zi = AIMAG(z)
zr = REAL(z, KIND=dp)
ak = AIMAG(zta)
IF (zr < 0.0_dp) THEN
  bk = REAL(zta, KIND=dp)
  ck = -ABS(bk)
  zta = CMPLX(ck, ak, KIND=dp)
END IF
IF (zi == 0.0_dp .AND. zr <= 0.0_dp) zta = CMPLX(0.0_dp, ak, KIND=dp)
aa = REAL(zta, KIND=dp)
IF (kode /= 2) THEN
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
  bb = ABS(aa)
  IF (bb >= alim) THEN
    bb = bb + 0.25_dp * LOG(az)
    sfac = tol
    IF (bb > elim) GO TO 40
  END IF
END IF
fmr = 0.0_dp
IF (aa < 0.0_dp .OR. zr <= 0.0_dp) THEN
  fmr = pi
  IF (zi < 0.0_dp) fmr = -pi
  zta = -zta
END IF
!-----------------------------------------------------------------------
!     AA=FACTOR FOR ANALYTIC CONTINUATION OF I(FNU,ZTA)
!     KODE=2 RETURNS EXP(-ABS(XZTA))*I(FNU,ZTA) FROM CBINU
!-----------------------------------------------------------------------
CALL cbinu(zta,fnu,kode,1,cy,nz,rl,fnul,tol,elim,alim)
IF (nz < 0) GO TO 50
aa = fmr * fnu
z3 = CMPLX(sfac, 0.0_dp, KIND=dp)
s1 = cy(1) * CMPLX(COS(aa), SIN(aa), KIND=dp) * z3
fnu = (2.0_dp - fid) / 3.0_dp
CALL cbinu(zta, fnu, kode, 2, cy, nz, rl, fnul, tol, elim, alim)
cy(1) = cy(1) * z3
cy(2) = cy(2) * z3
!-----------------------------------------------------------------------
!     BACKWARD RECUR ONE STEP FOR ORDERS -1/3 OR -2/3
!-----------------------------------------------------------------------
s2 = cy(1) * CMPLX(fnu+fnu, 0.0_dp, KIND=dp) / zta + cy(2)
aa = fmr * (fnu-1.0_dp)
s1 = (s1 + s2*CMPLX(COS(aa), SIN(aa), KIND=dp)) * coef
IF (id /= 1) THEN
  s1 = csq * s1
  bi = s1 / sfac
  RETURN
END IF
s1 = z * s1
bi = s1 / sfac
RETURN

30 aa = c1 * (1.0_dp-fid) + fid * c2
bi = CMPLX(aa, 0.0_dp, KIND=dp)
RETURN

40 nz = 0
ierr = 2
RETURN

50 IF (nz == -1) GO TO 40
nz = 0
ierr = 5
RETURN

60 ierr = 4
nz = 0
RETURN
END SUBROUTINE cbiry



SUBROUTINE cunik(zr, fnu, ikflg, ipmtr, tol, init, phi, zeta1, zeta2, total, &
                 cwrk)

!***BEGIN PROLOGUE  CUNIK
!***REFER TO  CBESI,CBESK

!  CUNIK COMPUTES PARAMETERS FOR THE UNIFORM ASYMPTOTIC EXPANSIONS OF
!  THE I AND K FUNCTIONS ON IKFLG= 1 OR 2 RESPECTIVELY BY

!  W(FNU,ZR) = PHI*EXP(ZETA)*SUM

!  WHERE     ZETA = -ZETA1 + ZETA2       OR
!                    ZETA1 - ZETA2

!  THE FIRST CALL MUST HAVE INIT=0.  SUBSEQUENT CALLS WITH THE SAME ZR
!  AND FNU WILL RETURN THE I OR K FUNCTION ON IKFLG= 1 OR 2 WITH NO CHANGE
!  IN INIT.  CWRK IS A COMPLEX WORK ARRAY.  IPMTR=0 COMPUTES ALL PARAMETERS.
!  IPMTR=1 COMPUTES PHI, ZETA1, ZETA2.

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  CUNIK

COMPLEX (dp), INTENT(IN)      :: zr
REAL (dp), INTENT(IN)         :: fnu
INTEGER, INTENT(IN)           :: ikflg
INTEGER, INTENT(IN)           :: ipmtr
REAL (dp), INTENT(IN)         :: tol
INTEGER, INTENT(IN OUT)       :: init
COMPLEX (dp), INTENT(OUT)     :: phi
COMPLEX (dp), INTENT(IN OUT)  :: zeta1
COMPLEX (dp), INTENT(IN OUT)  :: zeta2
COMPLEX (dp), INTENT(IN OUT)  :: total
COMPLEX (dp), INTENT(IN OUT)  :: cwrk(16)

COMPLEX (dp)  :: cfn, crfn, s, sr, t, t2,  zn
REAL (dp)     :: ac, rfn, test, tstr, tsti
INTEGER       :: i, j, k, l
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp, 0.0_dp)
REAL (dp), PARAMETER     :: con(2) = (/ 3.98942280401432678D-01, &
                                        1.25331413731550025_dp /)
REAL (dp), PARAMETER  :: c(120) = (/  &
  1.00000000000000000_dp, -2.08333333333333333D-01, 1.25000000000000000D-01, 3.34201388888888889D-01, &
 -4.01041666666666667D-01, 7.03125000000000000D-02, -1.02581259645061728_dp, 1.84646267361111111_dp,  &
 -8.91210937500000000D-01, 7.32421875000000000D-02, 4.66958442342624743_dp, -1.12070026162229938D+01,  &
  8.78912353515625000_dp, -2.36408691406250000_dp, 1.12152099609375000D-01, -2.82120725582002449D+01,  &
  8.46362176746007346D+01, -9.18182415432400174D+01, 4.25349987453884549D+01, -7.36879435947963170_dp,  &
  2.27108001708984375D-01, 2.12570130039217123D+02, -7.65252468141181642D+02, 1.05999045252799988D+03, &
 -6.99579627376132541D+02, 2.18190511744211590D+02, -2.64914304869515555D+01, 5.72501420974731445D-01, &
 -1.91945766231840700D+03, 8.06172218173730938D+03, -1.35865500064341374D+04, 1.16553933368645332D+04, &
 -5.30564697861340311D+03, 1.20090291321635246D+03, -1.08090919788394656D+02, 1.72772750258445740_dp,  &
  2.02042913309661486D+04, -9.69805983886375135D+04, 1.92547001232531532D+05, -2.03400177280415534D+05,  &
  1.22200464983017460D+05, -4.11926549688975513D+04, 7.10951430248936372D+03, -4.93915304773088012D+02,  &
  6.07404200127348304_dp, -2.42919187900551333D+05, 1.31176361466297720D+06, -2.99801591853810675D+06,  &
  3.76327129765640400D+06, -2.81356322658653411D+06, 1.26836527332162478D+06, -3.31645172484563578D+05,  &
  4.52187689813627263D+04, -2.49983048181120962D+03, 2.43805296995560639D+01, 3.28446985307203782D+06,  &
 -1.97068191184322269D+07, 5.09526024926646422D+07, -7.41051482115326577D+07, 6.63445122747290267D+07,  &
 -3.75671766607633513D+07, 1.32887671664218183D+07, -2.78561812808645469D+06, 3.08186404612662398D+05,  &
 -1.38860897537170405D+04, 1.10017140269246738D+02, -4.93292536645099620D+07, 3.25573074185765749D+08,  &
 -9.39462359681578403D+08, 1.55359689957058006D+09, -1.62108055210833708D+09, 1.10684281682301447D+09,  &
 -4.95889784275030309D+08, 1.42062907797533095D+08, -2.44740627257387285D+07, 2.24376817792244943D+06, &
 -8.40054336030240853D+04, 5.51335896122020586D+02, 8.14789096118312115D+08, -5.86648149205184723D+09,  &
  1.86882075092958249D+10, -3.46320433881587779D+10, 4.12801855797539740D+10, -3.30265997498007231D+10,  &
  1.79542137311556001D+10, -6.56329379261928433D+09, 1.55927986487925751D+09, -2.25105661889415278D+08,  &
  1.73951075539781645D+07, -5.49842327572288687D+05, 3.03809051092238427D+03, -1.46792612476956167D+10,  &
  1.14498237732025810D+11, -3.99096175224466498D+11, 8.19218669548577329D+11, -1.09837515608122331D+12,  &
  1.00815810686538209D+12, -6.45364869245376503D+11, 2.87900649906150589D+11, -8.78670721780232657D+10,  &
  1.76347306068349694D+10, -2.16716498322379509D+09, 1.43157876718888981D+08, -3.87183344257261262D+06,  &
  1.82577554742931747D+04, 2.86464035717679043D+11, -2.40629790002850396D+12, 9.10934118523989896D+12,  &
 -2.05168994109344374D+13, 3.05651255199353206D+13, -3.16670885847851584D+13, 2.33483640445818409D+13,  &
 -1.23204913055982872D+13, 4.61272578084913197D+12, -1.19655288019618160D+12, 2.05914503232410016D+11,  &
 -2.18229277575292237D+10, 1.24700929351271032D+09, -2.91883881222208134D+07, 1.18838426256783253D+05 /)

IF (init == 0) THEN
!-----------------------------------------------------------------------
!     INITIALIZE ALL VARIABLES
!-----------------------------------------------------------------------
  rfn = 1.0_dp / fnu
  crfn = rfn
  cwrk = czero
!     T = ZR*CRFN
!-----------------------------------------------------------------------
!     OVERFLOW TEST (ZR/FNU TOO SMALL)
!-----------------------------------------------------------------------
  tstr = REAL(zr, KIND=dp)
  tsti = AIMAG(zr)
  test = TINY(0.0_dp) * 1.0E+3
  ac = fnu * test
  IF (ABS(tstr) <= ac .AND. ABS(tsti) <= ac) THEN
    ac = 2.0_dp * ABS(LOG(test)) + fnu
    zeta1 = ac
    zeta2 = fnu
    phi = cone
    RETURN
  END IF
  t = zr * crfn
  s = cone + t * t
  sr = SQRT(s)
  cfn = fnu
  zn = (cone+sr) / t
  zeta1 = cfn * LOG(zn)
  zeta2 = cfn * sr
  t = cone / sr
  sr = t * crfn
  cwrk(16) = SQRT(sr)
  phi = cwrk(16) * con(ikflg)
  IF (ipmtr /= 0) RETURN
  t2 = cone / s
  cwrk(1) = cone
  crfn = cone
  ac = 1.0_dp
  l = 1
  DO  k = 2, 15
    s = czero
    DO  j = 1, k
      l = l + 1
      s = s * t2 + c(l)
    END DO
    crfn = crfn * sr
    cwrk(k) = crfn * s
    ac = ac * rfn
    tstr = REAL(cwrk(k), KIND=dp)
    tsti = AIMAG(cwrk(k))
    test = ABS(tstr) + ABS(tsti)
    IF (ac < tol .AND. test < tol) GO TO 30
  END DO
  k = 15

  30 init = k
END IF

IF (ikflg /= 2) THEN
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE I FUNCTION
!-----------------------------------------------------------------------
  total = SUM( cwrk(1:init) )
  phi = cwrk(16) * con(1)
  RETURN
END IF
!-----------------------------------------------------------------------
!     COMPUTE SUM FOR THE K FUNCTION
!-----------------------------------------------------------------------
s = czero
t = cone
DO  i = 1, init
  s = s + t * cwrk(i)
  t = -t
END DO
total = s
phi = cwrk(16) * con(2)
RETURN
END SUBROUTINE cunik



SUBROUTINE cuoik(z, fnu, kode, ikflg, n, y, nuf, tol, elim, alim)
!***BEGIN PROLOGUE  CUOIK
!***REFER TO  CBESI,CBESK,CBESH

!   CUOIK COMPUTES THE LEADING TERMS OF THE UNIFORM ASYMPTOTIC EXPANSIONS
!   FOR THE I AND K FUNCTIONS AND COMPARES THEM (IN LOGARITHMIC FORM)
!   TO ALIM AND ELIM FOR OVER AND UNDERFLOW, WHERE ALIM < ELIM.
!   IF THE MAGNITUDE, BASED ON THE LEADING EXPONENTIAL, IS LESS THAN ALIM OR
!   GREATER THAN -ALIM, THEN THE RESULT IS ON SCALE.
!   IF NOT, THEN A REFINED TEST USING OTHER MULTIPLIERS (IN LOGARITHMIC FORM)
!   IS MADE BASED ON ELIM.  HERE EXP(-ELIM) = SMALLEST MACHINE NUMBER*1000
!   AND EXP(-ALIM) = EXP(-ELIM)/TOL

!   IKFLG=1 MEANS THE I SEQUENCE IS TESTED
!        =2 MEANS THE K SEQUENCE IS TESTED
!   NUF = 0 MEANS THE LAST MEMBER OF THE SEQUENCE IS ON SCALE
!       =-1 MEANS AN OVERFLOW WOULD OCCUR
!   IKFLG=1 AND NUF > 0 MEANS THE LAST NUF Y VALUES WERE SET TO ZERO
!           THE FIRST N-NUF VALUES MUST BE SET BY ANOTHER ROUTINE
!   IKFLG=2 AND NUF = N MEANS ALL Y VALUES WERE SET TO ZERO
!   IKFLG=2 AND 0 < NUF < N NOT CONSIDERED.  Y MUST BE SET BY ANOTHER ROUTINE

!***ROUTINES CALLED  CUCHK,CUNHJ,CUNIK,R1MACH
!***END PROLOGUE  CUOIK

COMPLEX (dp), INTENT(IN)      :: z
REAL (dp), INTENT(IN)         :: fnu
INTEGER, INTENT(IN)           :: kode
INTEGER, INTENT(IN)           :: ikflg
INTEGER, INTENT(IN)           :: n
COMPLEX (dp), INTENT(IN OUT)  :: y(n)
INTEGER, INTENT(OUT)          :: nuf
REAL (dp), INTENT(IN)         :: tol
REAL (dp), INTENT(IN)         :: elim
REAL (dp), INTENT(IN)         :: alim

COMPLEX (dp)  :: arg, asum, bsum, cwrk(16), cz, phi, sum, zb, zeta1, zeta2, &
                 zn, zr
REAL (dp)     :: aarg, aphi, ascle, ax, ay, fnn, gnn, gnu, rcz, x, yy
INTEGER       :: iform, init, nn, nw
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp, 0.0_dp)
REAL (dp), PARAMETER     :: aic = 1.265512123484645396_dp

nuf = 0
nn = n
x = REAL(z, KIND=dp)
zr = z
IF (x < 0.0_dp) zr = -z
zb = zr
yy = AIMAG(zr)
ax = ABS(x) * 1.73205080756887_dp
ay = ABS(yy)
iform = 1
IF (ay > ax) iform = 2
gnu = MAX(fnu, 1.0_dp)
IF (ikflg /= 1) THEN
  fnn = nn
  gnn = fnu + fnn - 1.0_dp
  gnu = MAX(gnn, fnn)
END IF
!-----------------------------------------------------------------------
!     ONLY THE MAGNITUDE OF ARG AND PHI ARE NEEDED ALONG WITH THE
!     REAL PARTS OF ZETA1, ZETA2 AND ZB.  NO ATTEMPT IS MADE TO GET
!     THE SIGN OF THE IMAGINARY PART CORRECT.
!-----------------------------------------------------------------------
IF (iform /= 2) THEN
  init = 0
  CALL cunik(zr, gnu, ikflg, 1, tol, init, phi, zeta1, zeta2, sum, cwrk)
  cz = -zeta1 + zeta2
ELSE
  zn = -zr * CMPLX(0.0_dp, 1.0_dp, KIND=dp)
  IF (yy <= 0.0_dp) THEN
    zn = CONJG(-zn)
  END IF
  CALL cunhj(zn, gnu, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
  cz = -zeta1 + zeta2
  aarg = ABS(arg)
END IF
IF (kode == 2) cz = cz - zb
IF (ikflg == 2) cz = -cz
aphi = ABS(phi)
rcz = REAL(cz, KIND=dp)
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
IF (rcz <= elim) THEN
  IF (rcz >= alim) THEN
    rcz = rcz + LOG(aphi)
    IF (iform == 2) rcz = rcz - 0.25_dp * LOG(aarg) - aic
    IF (rcz > elim) GO TO 80
  ELSE
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
    IF (rcz >= -elim) THEN
      IF (rcz > -alim) GO TO 40
      rcz = rcz + LOG(aphi)
      IF (iform == 2) rcz = rcz - 0.25_dp * LOG(aarg) - aic
      IF (rcz > -elim) GO TO 30
    END IF

    10 y(1:nn) = czero
    nuf = nn
    RETURN

    30 ascle = 1.0E+3 * TINY(0.0_dp) / tol
    cz = cz + LOG(phi)
    IF (iform /= 1) THEN
      cz = cz - 0.25_dp * LOG(arg) - aic
    END IF
    ax = EXP(rcz) / tol
    ay = AIMAG(cz)
    cz = ax * CMPLX(COS(ay), SIN(ay), KIND=dp)
    CALL cuchk(cz, nw, ascle, tol)
    IF (nw == 1) GO TO 10
  END IF

  40 IF (ikflg == 2) RETURN
  IF (n == 1) RETURN
!-----------------------------------------------------------------------
!     SET UNDERFLOWS ON I SEQUENCE
!-----------------------------------------------------------------------
  50 gnu = fnu + (nn-1)
  IF (iform /= 2) THEN
    init = 0
    CALL cunik(zr, gnu, ikflg, 1, tol, init, phi, zeta1, zeta2, sum, cwrk)
    cz = -zeta1 + zeta2
  ELSE
    CALL cunhj(zn, gnu, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
    cz = -zeta1 + zeta2
    aarg = ABS(arg)
  END IF
  IF (kode == 2) cz = cz - zb
  aphi = ABS(phi)
  rcz = REAL(cz, KIND=dp)
  IF (rcz >= -elim) THEN
    IF (rcz > -alim) RETURN
    rcz = rcz + LOG(aphi)
    IF (iform == 2) rcz = rcz - 0.25_dp * LOG(aarg) - aic
    IF (rcz > -elim) GO TO 70
  END IF

  60 y(nn) = czero
  nn = nn - 1
  nuf = nuf + 1
  IF (nn == 0) RETURN
  GO TO 50

  70 ascle = 1.0E+3 * TINY(0.0_dp) / tol
  cz = cz + LOG(phi)
  IF (iform /= 1) THEN
    cz = cz - 0.25_dp * LOG(arg) - aic
  END IF
  ax = EXP(rcz) / tol
  ay = AIMAG(cz)
  cz = ax * CMPLX(COS(ay), SIN(ay), KIND=dp)
  CALL cuchk(cz, nw, ascle, tol)
  IF (nw == 1) GO TO 60
  RETURN
END IF

80 nuf = -1
RETURN
END SUBROUTINE cuoik



SUBROUTINE cwrsk(zr, fnu, kode, n, y, nz, cw, tol, elim, alim)
!***BEGIN PROLOGUE  CWRSK
!***REFER TO  CBESI,CBESK

!     CWRSK COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY
!     NORMALIZING THE I FUNCTION RATIOS FROM CRATI BY THE WRONSKIAN

!***ROUTINES CALLED  CBKNU,CRATI,R1MACH
!***END PROLOGUE  CWRSK

COMPLEX (dp), INTENT(IN)   :: zr
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
COMPLEX (dp), INTENT(OUT)  :: cw(2)
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: cinu, cscl, ct, c1, c2, rct, st
REAL (dp)     :: act, acw, ascle, s1, s2, yy
INTEGER       :: i, nw

!-----------------------------------------------------------------------
!     I(FNU+I-1,Z) BY BACKWARD RECURRENCE FOR RATIOS
!     Y(I)=I(FNU+I,Z)/I(FNU+I-1,Z) FROM CRATI NORMALIZED BY THE
!     WRONSKIAN WITH K(FNU,Z) AND K(FNU+1,Z) FROM CBKNU.
!-----------------------------------------------------------------------
nz = 0
CALL cbknu(zr, fnu, kode, 2, cw, nw, tol, elim, alim)
IF (nw == 0) THEN
  CALL crati(zr, fnu, n, y, tol)
!-----------------------------------------------------------------------
!     RECUR FORWARD ON I(FNU+1,Z) = R(FNU,Z)*I(FNU,Z),
!     R(FNU+J-1,Z)=Y(J),  J=1,...,N
!-----------------------------------------------------------------------
  cinu = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
  IF (kode /= 1) THEN
    yy = AIMAG(zr)
    s1 = COS(yy)
    s2 = SIN(yy)
    cinu = CMPLX(s1, s2, KIND=dp)
  END IF
!-----------------------------------------------------------------------
!     ON LOW EXPONENT MACHINES THE K FUNCTIONS CAN BE CLOSE TO BOTH THE
!     UNDER AND OVERFLOW LIMITS AND THE NORMALIZATION MUST BE SCALED TO
!     PREVENT OVER OR UNDERFLOW.  CUOIK HAS DETERMINED THAT THE RESULT
!     IS ON SCALE.
!-----------------------------------------------------------------------
  acw = ABS(cw(2))
  ascle = 1.0E+3 * TINY(0.0_dp) / tol
  cscl = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
  IF (acw <= ascle) THEN
    cscl = CMPLX(1.0_dp/tol, 0.0_dp, KIND=dp)
  ELSE
    ascle = 1.0_dp / ascle
    IF (acw >= ascle) THEN
      cscl = CMPLX(tol, 0.0_dp, KIND=dp)
    END IF
  END IF
  c1 = cw(1) * cscl
  c2 = cw(2) * cscl
  st = y(1)
!-----------------------------------------------------------------------
!     CINU=CINU*(CONJG(CT)/ABS(CT))*(1.0_dp/ABS(CT) PREVENTS
!     UNDER- OR OVERFLOW PREMATURELY BY SQUARING ABS(CT)
!-----------------------------------------------------------------------
  ct = zr * (c2 + st*c1)
  act = ABS(ct)
  rct = CMPLX(1.0_dp/act, 0.0_dp, KIND=dp)
  ct = CONJG(ct) * rct
  cinu = cinu * rct * ct
  y(1) = cinu * cscl
  IF (n == 1) RETURN
  DO  i = 2, n
    cinu = st * cinu
    st = y(i)
    y(i) = cinu * cscl
  END DO
  RETURN
END IF
nz = -1
IF (nw == -2) nz = -2
RETURN
END SUBROUTINE cwrsk



SUBROUTINE cmlri(z, fnu, kode, n, y, nz, tol)
!***BEGIN PROLOGUE  CMLRI
!***REFER TO  CBESI,CBESK

!     CMLRI COMPUTES THE I BESSEL FUNCTION FOR RE(Z) >= 0.0 BY THE
!     MILLER ALGORITHM NORMALIZED BY A NEUMANN SERIES.

!***ROUTINES CALLED  GAMLN,R1MACH
!***END PROLOGUE  CMLRI

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: tol

COMPLEX (dp)  :: ck, cnorm, pt, p1, p2, rz, sum
REAL (dp)     :: ack, ak, ap, at, az, bk, fkap, fkk, flam, fnf, rho,  &
                 rho2, scle, tfnf, tst, x
INTEGER       :: i, iaz, ifnu, inu, itime, k, kk, km, m

COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp), &
                            ctwo = (2.0_dp,0.0_dp)

scle = 1.0E+3 * TINY(0.0_dp) / tol
nz = 0
az = ABS(z)
x = REAL(z, KIND=dp)
iaz = az
ifnu = fnu
inu = ifnu + n - 1
at = iaz + 1
ck = CMPLX(at, 0.0_dp, KIND=dp) / z
rz = ctwo / z
p1 = czero
p2 = cone
ack = (at + 1.0_dp) / az
rho = ack + SQRT(ack*ack - 1.0_dp)
rho2 = rho * rho
tst = (rho2+rho2) / ((rho2-1.0_dp)*(rho-1.0_dp))
tst = tst / tol
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR INDEX FOR SERIES
!-----------------------------------------------------------------------
ak = at
DO  i = 1, 80
  pt = p2
  p2 = p1 - ck * p2
  p1 = pt
  ck = ck + rz
  ap = ABS(p2)
  IF (ap > tst*ak*ak) GO TO 20
  ak = ak + 1.0_dp
END DO
GO TO 90

20 i = i + 1
k = 0
IF (inu >= iaz) THEN
!-----------------------------------------------------------------------
!     COMPUTE RELATIVE TRUNCATION ERROR FOR RATIOS
!-----------------------------------------------------------------------
  p1 = czero
  p2 = cone
  at = inu + 1
  ck = CMPLX(at, 0.0_dp, KIND=dp) / z
  ack = at / az
  tst = SQRT(ack/tol)
  itime = 1
  DO  k = 1, 80
    pt = p2
    p2 = p1 - ck * p2
    p1 = pt
    ck = ck + rz
    ap = ABS(p2)
    IF (ap >= tst) THEN
      IF (itime == 2) GO TO 40
      ack = ABS(ck)
      flam = ack + SQRT(ack*ack - 1.0_dp)
      fkap = ap / ABS(p1)
      rho = MIN(flam,fkap)
      tst = tst * SQRT(rho/(rho*rho - 1.0_dp))
      itime = 2
    END IF
  END DO
  GO TO 90
END IF
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE AND SUM NORMALIZING RELATION
!-----------------------------------------------------------------------
40 k = k + 1
kk = MAX(i+iaz, k+inu)
fkk = kk
p1 = czero
!-----------------------------------------------------------------------
!     SCALE P2 AND SUM BY SCLE
!-----------------------------------------------------------------------
p2 = CMPLX(scle, 0.0_dp, KIND=dp)
fnf = fnu - ifnu
tfnf = fnf + fnf
bk = gamln(fkk+tfnf+1.0_dp) - gamln(fkk+1.0_dp) - gamln(tfnf+1.0_dp)
bk = EXP(bk)
sum = czero
km = kk - inu
DO  i = 1, km
  pt = p2
  p2 = p1 + CMPLX(fkk+fnf, 0.0_dp, KIND=dp) * rz * p2
  p1 = pt
  ak = 1.0_dp - tfnf / (fkk+tfnf)
  ack = bk * ak
  sum = sum + CMPLX(ack+bk, 0.0_dp, KIND=dp) * p1
  bk = ack
  fkk = fkk - 1.0_dp
END DO
y(n) = p2
IF (n /= 1) THEN
  DO  i = 2, n
    pt = p2
    p2 = p1 + CMPLX(fkk+fnf, 0.0_dp, KIND=dp) * rz * p2
    p1 = pt
    ak = 1.0_dp - tfnf / (fkk+tfnf)
    ack = bk * ak
    sum = sum + CMPLX(ack+bk, 0.0_dp, KIND=dp) * p1
    bk = ack
    fkk = fkk - 1.0_dp
    m = n - i + 1
    y(m) = p2
  END DO
END IF
IF (ifnu > 0) THEN
  DO  i = 1, ifnu
    pt = p2
    p2 = p1 + CMPLX(fkk+fnf, 0.0_dp, KIND=dp) * rz * p2
    p1 = pt
    ak = 1.0_dp - tfnf / (fkk+tfnf)
    ack = bk * ak
    sum = sum + CMPLX(ack+bk, 0.0_dp, KIND=dp) * p1
    bk = ack
    fkk = fkk - 1.0_dp
  END DO
END IF
pt = z
IF (kode == 2) pt = pt - x
p1 = -fnf * LOG(rz) + pt
ap = gamln(1.0_dp+fnf)
pt = p1 - ap
!-----------------------------------------------------------------------
!     THE DIVISION EXP(PT)/(SUM+P2) IS ALTERED TO AVOID OVERFLOW
!     IN THE DENOMINATOR BY SQUARING LARGE QUANTITIES
!-----------------------------------------------------------------------
p2 = p2 + sum
ap = ABS(p2)
p1 = CMPLX(1.0_dp/ap, 0.0_dp, KIND=dp)
ck = EXP(pt) * p1
pt = CONJG(p2) * p1
cnorm = ck * pt
y(1:n) = y(1:n) * cnorm
RETURN

90 nz = -2
RETURN
END SUBROUTINE cmlri



SUBROUTINE cunhj(z, fnu, ipmtr, tol, phi, arg, zeta1, zeta2, asum, bsum)
!***BEGIN PROLOGUE  CUNHJ
!***REFER TO  CBESI,CBESK

!  REFERENCES
!      HANDBOOK OF MATHEMATICAL FUNCTIONS BY M. ABRAMOWITZ AND I.A.
!      STEGUN, AMS55, NATIONAL BUREAU OF STANDARDS, 1965, CHAPTER 9.

!      ASYMPTOTICS AND SPECIAL FUNCTIONS BY F.W.J. OLVER, ACADEMIC
!      PRESS, N.Y., 1974, PAGE 420

!  ABSTRACT
!      CUNHJ COMPUTES PARAMETERS FOR BESSEL FUNCTIONS C(FNU,Z) =
!      J(FNU,Z), Y(FNU,Z) OR H(I,FNU,Z) I=1,2 FOR LARGE ORDERS FNU
!      BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSION

!      C(FNU,Z) = C1*PHI*( ASUM*AIRY(ARG) + C2*BSUM*DAIRY(ARG) )

!      FOR PROPER CHOICES OF C1, C2, AIRY AND DAIRY WHERE AIRY IS
!      AN AIRY FUNCTION AND DAIRY IS ITS DERIVATIVE.

!            (2/3)*FNU*ZETA**1.5 = ZETA1-ZETA2,

!      ZETA1 = 0.5*FNU*LOG((1+W)/(1-W)), ZETA2 = FNU*W FOR SCALING
!      PURPOSES IN AIRY FUNCTIONS FROM CAIRY OR CBIRY.

!      MCONJ = SIGN OF AIMAG(Z), BUT IS AMBIGUOUS WHEN Z IS REAL AND
!      MUST BE SPECIFIED.  IPMTR=0 RETURNS ALL PARAMETERS.  IPMTR =
!      1 COMPUTES ALL EXCEPT ASUM AND BSUM.

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  CUNHJ

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: ipmtr
REAL (dp), INTENT(IN)      :: tol
COMPLEX (dp), INTENT(OUT)  :: phi
COMPLEX (dp), INTENT(OUT)  :: arg
COMPLEX (dp), INTENT(OUT)  :: zeta1
COMPLEX (dp), INTENT(OUT)  :: zeta2
COMPLEX (dp), INTENT(OUT)  :: asum
COMPLEX (dp), INTENT(OUT)  :: bsum

COMPLEX (dp)  :: cfnu, cr(14), dr(14), p(30), przth, ptfn, rtzta, rzth,  &
                 suma, sumb, tfn, t2, up(14), w, w2, za, zb, zc, zeta, zth
REAL (dp)     :: ang, ap(30), atol, aw2, azth, btol, fn13, fn23, pp, rfn13, &
                 rfnu, rfnu2, wi, wr, zci, zcr, zetai, zetar, zthi, zthr,  &
                 asumr, asumi, bsumr, bsumi, test, tstr, tsti, ac
INTEGER       :: ias, ibs, is, j, jr, ju, k, kmax, kp1, ks, l, lr, lrp1,  &
                 l1, l2, m
REAL (dp), PARAMETER  :: ar(14) = (/  &
    1.00000000000000000_dp, 1.04166666666666667D-01,  &
    8.35503472222222222D-02, 1.28226574556327160D-01,  &
    2.91849026464140464D-01, 8.81627267443757652D-01,  &
    3.32140828186276754_dp, 1.49957629868625547D+01,  &
    7.89230130115865181D+01, 4.74451538868264323D+02,  &
    3.20749009089066193D+03, 2.40865496408740049D+04,  &
    1.98923119169509794D+05, 1.79190200777534383D+06 /)
REAL (dp), PARAMETER  :: br(14) = (/  &
    1.00000000000000000_dp, -1.45833333333333333D-01,  &
   -9.87413194444444444D-02, -1.43312053915895062D-01,  &
   -3.17227202678413548D-01, -9.42429147957120249D-01,  &
   -3.51120304082635426_dp, -1.57272636203680451D+01,  &
   -8.22814390971859444D+01, -4.92355370523670524D+02,  &
   -3.31621856854797251D+03, -2.48276742452085896D+04,  &
   -2.04526587315129788D+05, -1.83844491706820990D+06 /)
REAL (dp), PARAMETER  :: c(105) = (/  &
    1.00000000000000000_dp, -2.08333333333333333D-01, 1.25000000000000000D-01, &
    3.34201388888888889D-01, -4.01041666666666667D-01, 7.03125000000000000D-02, &
   -1.02581259645061728_dp, 1.84646267361111111_dp, -8.91210937500000000D-01,  &
    7.32421875000000000D-02, 4.66958442342624743_dp, -1.12070026162229938D+01,  &
    8.78912353515625000_dp, -2.36408691406250000_dp, 1.12152099609375000D-01,  &
   -2.82120725582002449D+01, 8.46362176746007346D+01, -9.18182415432400174D+01,  &
    4.25349987453884549D+01, -7.36879435947963170_dp, 2.27108001708984375D-01,  &
    2.12570130039217123D+02, -7.65252468141181642D+02, 1.05999045252799988D+03,  &
   -6.99579627376132541D+02, 2.18190511744211590D+02, -2.64914304869515555D+01, &
    5.72501420974731445D-01, -1.91945766231840700D+03, 8.06172218173730938D+03, &
   -1.35865500064341374D+04, 1.16553933368645332D+04, -5.30564697861340311D+03, &
    1.20090291321635246D+03, -1.08090919788394656D+02, 1.72772750258445740_dp,  &
    2.02042913309661486D+04, -9.69805983886375135D+04, 1.92547001232531532D+05, &
   -2.03400177280415534D+05, 1.22200464983017460D+05, -4.11926549688975513D+04,  &
    7.10951430248936372D+03, -4.93915304773088012D+02, 6.07404200127348304_dp,  &
   -2.42919187900551333D+05, 1.31176361466297720D+06, -2.99801591853810675D+06, &
    3.76327129765640400D+06, -2.81356322658653411D+06, 1.26836527332162478D+06, &
   -3.31645172484563578D+05, 4.52187689813627263D+04, -2.49983048181120962D+03,  &
    2.43805296995560639D+01, 3.28446985307203782D+06, -1.97068191184322269D+07, &
    5.09526024926646422D+07, -7.41051482115326577D+07, 6.63445122747290267D+07, &
   -3.75671766607633513D+07, 1.32887671664218183D+07, -2.78561812808645469D+06, &
    3.08186404612662398D+05, -1.38860897537170405D+04, 1.10017140269246738D+02, &
   -4.93292536645099620D+07, 3.25573074185765749D+08, -9.39462359681578403D+08, &
    1.55359689957058006D+09, -1.62108055210833708D+09, 1.10684281682301447D+09, &
   -4.95889784275030309D+08, 1.42062907797533095D+08, -2.44740627257387285D+07, &
    2.24376817792244943D+06, -8.40054336030240853D+04, 5.51335896122020586D+02,  &
    8.14789096118312115D+08, -5.86648149205184723D+09, 1.86882075092958249D+10, &
   -3.46320433881587779D+10, 4.12801855797539740D+10, -3.30265997498007231D+10,  &
    1.79542137311556001D+10, -6.56329379261928433D+09, 1.55927986487925751D+09, &
   -2.25105661889415278D+08, 1.73951075539781645D+07, -5.49842327572288687D+05,  &
    3.03809051092238427D+03, -1.46792612476956167D+10, 1.14498237732025810D+11, &
   -3.99096175224466498D+11, 8.19218669548577329D+11, -1.09837515608122331D+12, &
    1.00815810686538209D+12, -6.45364869245376503D+11, 2.87900649906150589D+11, &
   -8.78670721780232657D+10, 1.76347306068349694D+10, -2.16716498322379509D+09,  &
    1.43157876718888981D+08, -3.87183344257261262D+06, 1.82577554742931747D+04 /)
REAL (dp), PARAMETER  :: alfa1(30) = (/  &
 -4.44444444444444444D-03, -9.22077922077922078D-04, -8.84892884892884893D-05, &
  1.65927687832449737D-04, 2.46691372741792910D-04, 2.65995589346254780D-04,  &
  2.61824297061500945D-04, 2.48730437344655609D-04, 2.32721040083232098D-04,  &
  2.16362485712365082D-04, 2.00738858762752355D-04, 1.86267636637545172D-04,  &
  1.73060775917876493D-04, 1.61091705929015752D-04, 1.50274774160908134D-04,  &
  1.40503497391269794D-04, 1.31668816545922806D-04, 1.23667445598253261D-04,  &
  1.16405271474737902D-04, 1.09798298372713369D-04, 1.03772410422992823D-04,  &
  9.82626078369363448D-05, 9.32120517249503256D-05, 8.85710852478711718D-05,  &
  8.42963105715700223D-05, 8.03497548407791151D-05, 7.66981345359207388D-05,  &
  7.33122157481777809D-05, 7.01662625163141333D-05, 6.72375633790160292D-05 /)
REAL (dp), PARAMETER  :: alfa2(30) = (/  &
  6.93735541354588974D-04, 2.32241745182921654D-04, -1.41986273556691197D-05, &
 -1.16444931672048640D-04, -1.50803558053048762D-04, -1.55121924918096223D-04, &
 -1.46809756646465549D-04, -1.33815503867491367D-04, -1.19744975684254051D-04, &
 -1.06184319207974020D-04, -9.37699549891194492D-05, -8.26923045588193274D-05, &
 -7.29374348155221211D-05, -6.44042357721016283D-05, -5.69611566009369048D-05, &
 -5.04731044303561628D-05, -4.48134868008882786D-05, -3.98688727717598864D-05, &
 -3.55400532972042498D-05, -3.17414256609022480D-05, -2.83996793904174811D-05, &
 -2.54522720634870566D-05, -2.28459297164724555D-05, -2.05352753106480604D-05, &
 -1.84816217627666085D-05, -1.66519330021393806D-05, -1.50179412980119482D-05, &
 -1.35554031379040526D-05, -1.22434746473858131D-05, -1.10641884811308169D-05 /)
REAL (dp), PARAMETER  :: alfa3(30) = (/  &
 -3.54211971457743841D-04, -1.56161263945159416D-04, 3.04465503594936410D-05, &
  1.30198655773242693D-04, 1.67471106699712269D-04, 1.70222587683592569D-04,  &
  1.56501427608594704D-04, 1.36339170977445120D-04, 1.14886692029825128D-04,  &
  9.45869093034688111D-05, 7.64498419250898258D-05, 6.07570334965197354D-05,  &
  4.74394299290508799D-05, 3.62757512005344297D-05, 2.69939714979224901D-05,  &
  1.93210938247939253D-05, 1.30056674793963203D-05, 7.82620866744496661D-06,  &
  3.59257485819351583D-06, 1.44040049814251817D-07, -2.65396769697939116D-06, &
 -4.91346867098485910D-06, -6.72739296091248287D-06, -8.17269379678657923D-06, &
 -9.31304715093561232D-06, -1.02011418798016441D-05, -1.08805962510592880D-05, &
 -1.13875481509603555D-05, -1.17519675674556414D-05, -1.19987364870944141D-05 /)
REAL (dp), PARAMETER  :: alfa4(30) = (/  &
  3.78194199201772914D-04, 2.02471952761816167D-04, -6.37938506318862408D-05, &
 -2.38598230603005903D-04, -3.10916256027361568D-04, -3.13680115247576316D-04, &
 -2.78950273791323387D-04, -2.28564082619141374D-04, -1.75245280340846749D-04, &
 -1.25544063060690348D-04, -8.22982872820208365D-05, -4.62860730588116458D-05, &
 -1.72334302366962267D-05, 5.60690482304602267D-06, 2.31395443148286800D-05,  &
  3.62642745856793957D-05, 4.58006124490188752D-05, 5.24595294959114050D-05,  &
  5.68396208545815266D-05, 5.94349820393104052D-05, 6.06478527578421742D-05,  &
  6.08023907788436497D-05, 6.01577894539460388D-05, 5.89199657344698500D-05,  &
  5.72515823777593053D-05, 5.52804375585852577D-05, 5.31063773802880170D-05,  &
  5.08069302012325706D-05, 4.84418647620094842D-05, 4.60568581607475370D-05 /)
REAL (dp), PARAMETER  :: alfa5(30) = (/  &
 -6.91141397288294174D-04, -4.29976633058871912D-04, 1.83067735980039018D-04, &
  6.60088147542014144D-04, 8.75964969951185931D-04, 8.77335235958235514D-04,  &
  7.49369585378990637D-04, 5.63832329756980918D-04, 3.68059319971443156D-04,  &
  1.88464535514455599D-04, 3.70663057664904149D-05, -8.28520220232137023D-05, &
 -1.72751952869172998D-04, -2.36314873605872983D-04, -2.77966150694906658D-04, &
 -3.02079514155456919D-04, -3.12594712643820127D-04, -3.12872558758067163D-04, &
 -3.05678038466324377D-04, -2.93226470614557331D-04, -2.77255655582934777D-04, &
 -2.59103928467031709D-04, -2.39784014396480342D-04, -2.20048260045422848D-04, &
 -2.00443911094971498D-04, -1.81358692210970687D-04, -1.63057674478657464D-04, &
 -1.45712672175205844D-04, -1.29425421983924587D-04, -1.14245691942445952D-04 /)
REAL (dp), PARAMETER  :: alfa6(30) = (/  &
  1.92821964248775885D-03, 1.35592576302022234D-03, -7.17858090421302995D-04, &
 -2.58084802575270346D-03, -3.49271130826168475D-03, -3.46986299340960628D-03, &
 -2.82285233351310182D-03, -1.88103076404891354D-03, -8.89531718383947600D-04, &
  3.87912102631035228D-06, 7.28688540119691412D-04, 1.26566373053457758D-03,  &
  1.62518158372674427D-03, 1.83203153216373172D-03, 1.91588388990527909D-03,  &
  1.90588846755546138D-03, 1.82798982421825727D-03, 1.70389506421121530D-03,  &
  1.55097127171097686D-03, 1.38261421852276159D-03, 1.20881424230064774D-03,  &
  1.03676532638344962D-03, 8.71437918068619115D-04, 7.16080155297701002D-04,  &
  5.72637002558129372D-04, 4.42089819465802277D-04, 3.24724948503090564D-04,  &
  2.20342042730246599D-04, 1.28412898401353882D-04, 4.82005924552095464D-05 /)
REAL (dp)  :: alfa(180)
REAL (dp), PARAMETER  :: beta1(30) = (/  &
  1.79988721413553309D-02, 5.59964911064388073D-03, 2.88501402231132779D-03,  &
  1.80096606761053941D-03, 1.24753110589199202D-03, 9.22878876572938311D-04,  &
  7.14430421727287357D-04, 5.71787281789704872D-04, 4.69431007606481533D-04,  &
  3.93232835462916638D-04, 3.34818889318297664D-04, 2.88952148495751517D-04,  &
  2.52211615549573284D-04, 2.22280580798883327D-04, 1.97541838033062524D-04,  &
  1.76836855019718004D-04, 1.59316899661821081D-04, 1.44347930197333986D-04,  &
  1.31448068119965379D-04, 1.20245444949302884D-04, 1.10449144504599392D-04,  &
  1.01828770740567258D-04, 9.41998224204237509D-05, 8.74130545753834437D-05,  &
  8.13466262162801467D-05, 7.59002269646219339D-05, 7.09906300634153481D-05,  &
  6.65482874842468183D-05, 6.25146958969275078D-05, 5.88403394426251749D-05 /)
REAL (dp), PARAMETER  :: beta2(30) = (/  &
 -1.49282953213429172D-03, -8.78204709546389328D-04, -5.02916549572034614D-04, &
 -2.94822138512746025D-04, -1.75463996970782828D-04, -1.04008550460816434D-04, &
 -5.96141953046457895D-05, -3.12038929076098340D-05, -1.26089735980230047D-05, &
 -2.42892608575730389D-07, 8.05996165414273571D-06, 1.36507009262147391D-05,  &
  1.73964125472926261D-05, 1.98672978842133780D-05, 2.14463263790822639D-05,  &
  2.23954659232456514D-05, 2.28967783814712629D-05, 2.30785389811177817D-05,  &
  2.30321976080909144D-05, 2.28236073720348722D-05, 2.25005881105292418D-05,  &
  2.20981015361991429D-05, 2.16418427448103905D-05, 2.11507649256220843D-05,  &
  2.06388749782170737D-05, 2.01165241997081666D-05, 1.95913450141179244D-05,  &
  1.90689367910436740D-05, 1.85533719641636667D-05, 1.80475722259674218D-05 /)
REAL (dp), PARAMETER  :: beta3(30) = (/  &
  5.52213076721292790D-04, 4.47932581552384646D-04, 2.79520653992020589D-04,  &
  1.52468156198446602D-04, 6.93271105657043598D-05, 1.76258683069991397D-05,  &
 -1.35744996343269136D-05, -3.17972413350427135D-05, -4.18861861696693365D-05, &
 -4.69004889379141029D-05, -4.87665447413787352D-05, -4.87010031186735069D-05, &
 -4.74755620890086638D-05, -4.55813058138628452D-05, -4.33309644511266036D-05, &
 -4.09230193157750364D-05, -3.84822638603221274D-05, -3.60857167535410501D-05, &
 -3.37793306123367417D-05, -3.15888560772109621D-05, -2.95269561750807315D-05, &
 -2.75978914828335759D-05, -2.58006174666883713D-05, -2.41308356761280200D-05, &
 -2.25823509518346033D-05, -2.11479656768912971D-05, -1.98200638885294927D-05, &
 -1.85909870801065077D-05, -1.74532699844210224D-05, -1.63997823854497997D-05 /)
REAL (dp), PARAMETER  :: beta4(30) = (/  &
 -4.74617796559959808D-04, -4.77864567147321487D-04, -3.20390228067037603D-04, &
 -1.61105016119962282D-04, -4.25778101285435204D-05, 3.44571294294967503D-05, &
  7.97092684075674924D-05, 1.03138236708272200D-04, 1.12466775262204158D-04,  &
  1.13103642108481389D-04, 1.08651634848774268D-04, 1.01437951597661973D-04,  &
  9.29298396593363896D-05, 8.40293133016089978D-05, 7.52727991349134062D-05,  &
  6.69632521975730872D-05, 5.92564547323194704D-05, 5.22169308826975567D-05,  &
  4.58539485165360646D-05, 4.01445513891486808D-05, 3.50481730031328081D-05,  &
  3.05157995034346659D-05, 2.64956119950516039D-05, 2.29363633690998152D-05,  &
  1.97893056664021636D-05, 1.70091984636412623D-05, 1.45547428261524004D-05,  &
  1.23886640995878413D-05, 1.04775876076583236D-05, 8.79179954978479373D-06 /)
REAL (dp), PARAMETER  :: beta5(30) = (/  &
  7.36465810572578444D-04, 8.72790805146193976D-04, 6.22614862573135066D-04,  &
  2.85998154194304147D-04, 3.84737672879366102D-06, -1.87906003636971558D-04, &
 -2.97603646594554535D-04, -3.45998126832656348D-04, -3.53382470916037712D-04, &
 -3.35715635775048757D-04, -3.04321124789039809D-04, -2.66722723047612821D-04, &
 -2.27654214122819527D-04, -1.89922611854562356D-04, -1.55058918599093870D-04, &
 -1.23778240761873630D-04, -9.62926147717644187D-05, -7.25178327714425337D-05, &
 -5.22070028895633801D-05, -3.50347750511900522D-05, -2.06489761035551757D-05, &
 -8.70106096849767054D-06, 1.13698686675100290D-06, 9.16426474122778849D-06,  &
  1.56477785428872620D-05, 2.08223629482466847D-05, 2.48923381004595156D-05,  &
  2.80340509574146325D-05, 3.03987774629861915D-05, 3.21156731406700616D-05 /)
REAL (dp), PARAMETER  :: beta6(30) = (/  &
 -1.80182191963885708D-03, -2.43402962938042533D-03, -1.83422663549856802D-03, &
 -7.62204596354009765D-04, 2.39079475256927218D-04, 9.49266117176881141D-04,  &
  1.34467449701540359D-03, 1.48457495259449178D-03, 1.44732339830617591D-03,  &
  1.30268261285657186D-03, 1.10351597375642682D-03, 8.86047440419791759D-04,  &
  6.73073208165665473D-04, 4.77603872856582378D-04, 3.05991926358789362D-04,  &
  1.60315694594721630D-04, 4.00749555270613286D-05, -5.66607461635251611D-05, &
 -1.32506186772982638D-04, -1.90296187989614057D-04, -2.32811450376937408D-04, &
 -2.62628811464668841D-04, -2.82050469867598672D-04, -2.93081563192861167D-04, &
 -2.97435962176316616D-04, -2.96557334239348078D-04, -2.91647363312090861D-04, &
 -2.83696203837734166D-04, -2.73512317095673346D-04, -2.61750155806768580D-04 /)
REAL (dp), PARAMETER  :: beta7(30) = (/  &
  6.38585891212050914D-03, 9.62374215806377941D-03, 7.61878061207001043D-03,  &
  2.83219055545628054D-03, -2.09841352012720090D-03, -5.73826764216626498D-03, &
 -7.70804244495414620D-03, -8.21011692264844401D-03, -7.65824520346905413D-03, &
 -6.47209729391045177D-03, -4.99132412004966473D-03, -3.45612289713133280D-03, &
 -2.01785580014170775D-03, -7.59430686781961401D-04, 2.84173631523859138D-04, &
  1.10891667586337403D-03, 1.72901493872728771D-03, 2.16812590802684701D-03,  &
  2.45357710494539735D-03, 2.61281821058334862D-03, 2.67141039656276912D-03,  &
  2.65203073395980430D-03, 2.57411652877287315D-03, 2.45389126236094427D-03,  &
  2.30460058071795494D-03, 2.13684837686712662D-03, 1.95896528478870911D-03,  &
  1.77737008679454412D-03, 1.59690280765839059D-03, 1.42111975664438546D-03 /)
REAL (dp)  :: beta(210)
REAL (dp), PARAMETER  :: gama(30) = (/  &
  6.29960524947436582D-01, 2.51984209978974633D-01, 1.54790300415655846D-01, &
  1.10713062416159013D-01, 8.57309395527394825D-02, 6.97161316958684292D-02, &
  5.86085671893713576D-02, 5.04698873536310685D-02, 4.42600580689154809D-02, &
  3.93720661543509966D-02, 3.54283195924455368D-02, 3.21818857502098231D-02, &
  2.94646240791157679D-02, 2.71581677112934479D-02, 2.51768272973861779D-02, &
  2.34570755306078891D-02, 2.19508390134907203D-02, 2.06210828235646240D-02, &
  1.94388240897880846D-02, 1.83810633800683158D-02, 1.74293213231963172D-02, &
  1.65685837786612353D-02, 1.57865285987918445D-02, 1.50729501494095594D-02, &
  1.44193250839954639D-02, 1.38184805735341786D-02, 1.32643378994276568D-02, &
  1.27517121970498651D-02, 1.22761545318762767D-02, 1.18338262398482403D-02 /)
REAL (dp), PARAMETER  :: ex1 = 3.33333333333333333D-01,  &
    ex2 = 6.66666666666666667D-01, hpi = 1.57079632679489662_dp,  &
    pi = 3.14159265358979324_dp, thpi = 4.71238898038468986_dp
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp)

! Associate arrays alfa & beta

alfa(  1: 30) = alfa1
alfa( 31: 60) = alfa2
alfa( 61: 90) = alfa3
alfa( 91:120) = alfa4
alfa(121:150) = alfa5
alfa(151:180) = alfa6
beta(  1: 30) = beta1
beta( 31: 60) = beta2
beta( 61: 90) = beta3
beta( 91:120) = beta4
beta(121:150) = beta5
beta(151:180) = beta6
beta(181:210) = beta7

rfnu = 1.0_dp / fnu
!     ZB = Z*CMPLX(RFNU,0.0_dp)
!-----------------------------------------------------------------------
!     OVERFLOW TEST (Z/FNU TOO SMALL)
!-----------------------------------------------------------------------
tstr = REAL(z, KIND=dp)
tsti = AIMAG(z)
test = TINY(0.0_dp) * 1.0E+3
ac = fnu * test
IF (ABS(tstr) <= ac .AND. ABS(tsti) <= ac) THEN
  ac = 2.0_dp * ABS(LOG(test)) + fnu
  zeta1 = ac
  zeta2 = fnu
  phi = cone
  arg = cone
  RETURN
END IF
zb = z * rfnu
rfnu2 = rfnu * rfnu
!-----------------------------------------------------------------------
!     COMPUTE IN THE FOURTH QUADRANT
!-----------------------------------------------------------------------
fn13 = fnu ** ex1
fn23 = fn13 * fn13
rfn13 = 1.0_dp/fn13
w2 = cone - zb * zb
aw2 = ABS(w2)
IF (aw2 > 0.25_dp) GO TO 110
!-----------------------------------------------------------------------
!     POWER SERIES FOR ABS(W2) <= 0.25_dp
!-----------------------------------------------------------------------
k = 1

p(1) = cone
suma = gama(1)
ap(1) = 1.0_dp
IF (aw2 >= tol) THEN
  DO  k = 2, 30
    p(k) = p(k-1) * w2
    suma = suma + p(k) * gama(k)
    ap(k) = ap(k-1) * aw2
    IF (ap(k) < tol) GO TO 20
  END DO
  k = 30
END IF

20 kmax = k
zeta = w2 * suma
arg = zeta * fn23
za = SQRT(suma)
zeta2 = SQRT(w2) * fnu
zeta1 = zeta2 * (cone + zeta*za*ex2)
za = za + za
phi = SQRT(za) * rfn13
IF (ipmtr /= 1) THEN
!-----------------------------------------------------------------------
!     SUM SERIES FOR ASUM AND BSUM
!-----------------------------------------------------------------------
  sumb = czero
  DO k = 1, kmax
    sumb = sumb + p(k)*beta(k)
  END DO
  asum = czero
  bsum = sumb
  l1 = 0
  l2 = 30
  btol = tol * (ABS(REAL(bsum)) + ABS(AIMAG(bsum)))
  atol = tol
  pp = 1.0_dp
  ias = 0
  ibs = 0
  IF (rfnu2 >= tol) THEN
    DO  is = 2, 7
      atol = atol / rfnu2
      pp = pp * rfnu2
      IF (ias /= 1) THEN
        suma = czero
        DO  k = 1, kmax
          m = l1 + k
          suma = suma + p(k) * alfa(m)
          IF (ap(k) < atol) EXIT
        END DO
        asum = asum + suma * pp
        IF (pp < tol) ias = 1
      END IF
      IF (ibs /= 1) THEN
        sumb = czero
        DO  k = 1, kmax
          m = l2 + k
          sumb = sumb + p(k) * beta(m)
          IF (ap(k) < atol) EXIT
        END DO
        bsum = bsum + sumb * pp
        IF (pp < btol) ibs = 1
      END IF
      IF (ias == 1 .AND. ibs == 1) EXIT
      l1 = l1 + 30
      l2 = l2 + 30
    END DO
  END IF

  asum = asum + cone
  pp = rfnu * rfn13
  bsum = bsum * pp
END IF

100 RETURN
!-----------------------------------------------------------------------
!     ABS(W2) > 0.25_dp
!-----------------------------------------------------------------------
110 w = SQRT(w2)
wr = REAL(w, KIND=dp)
wi = AIMAG(w)
IF (wr < 0.0_dp) wr = 0.0_dp
IF (wi < 0.0_dp) wi = 0.0_dp
w = CMPLX(wr, wi, KIND=dp)
za = (cone+w) / zb
zc = LOG(za)
zcr = REAL(zc, KIND=dp)
zci = AIMAG(zc)
IF (zci < 0.0_dp) zci = 0.0_dp
IF (zci > hpi) zci = hpi
IF (zcr < 0.0_dp) zcr = 0.0_dp
zc = CMPLX(zcr, zci, KIND=dp)
zth = (zc-w) * 1.5_dp
cfnu = CMPLX(fnu, 0.0_dp, KIND=dp)
zeta1 = zc * cfnu
zeta2 = w * cfnu
azth = ABS(zth)
zthr = REAL(zth, KIND=dp)
zthi = AIMAG(zth)
ang = thpi
IF (zthr < 0.0_dp .OR. zthi >= 0.0_dp) THEN
  ang = hpi
  IF (zthr /= 0.0_dp) THEN
    ang = ATAN(zthi/zthr)
    IF (zthr < 0.0_dp) ang = ang + pi
  END IF
END IF
pp = azth ** ex2
ang = ang * ex2
zetar = pp * COS(ang)
zetai = pp * SIN(ang)
IF (zetai < 0.0_dp) zetai = 0.0_dp
zeta = CMPLX(zetar, zetai, KIND=dp)
arg = zeta * fn23
rtzta = zth / zeta
za = rtzta / w
phi = SQRT(za+za) * rfn13
IF (ipmtr == 1) GO TO 100
tfn = CMPLX(rfnu, 0.0_dp, KIND=dp) / w
rzth = CMPLX(rfnu, 0.0_dp, KIND=dp) / zth
zc = rzth * ar(2)
t2 = cone / w2
up(2) = (t2*c(2) + c(3)) * tfn
bsum = up(2) + zc
asum = czero
IF (rfnu >= tol) THEN
  przth = rzth
  ptfn = tfn
  up(1) = cone
  pp = 1.0_dp
  bsumr = REAL(bsum, KIND=dp)
  bsumi = AIMAG(bsum)
  btol = tol * (ABS(bsumr) + ABS(bsumi))
  ks = 0
  kp1 = 2
  l = 3
  ias = 0
  ibs = 0
  DO  lr = 2, 12, 2
    lrp1 = lr + 1
!-----------------------------------------------------------------------
!     COMPUTE TWO ADDITIONAL CR, DR, AND UP FOR TWO MORE TERMS IN
!     NEXT SUMA AND SUMB
!-----------------------------------------------------------------------
    DO  k = lr, lrp1
      ks = ks + 1
      kp1 = kp1 + 1
      l = l + 1
      za = CMPLX(c(l), 0.0_dp, KIND=dp)
      DO  j = 2, kp1
        l = l + 1
        za = za * t2 + c(l)
      END DO
      ptfn = ptfn * tfn
      up(kp1) = ptfn * za
      cr(ks) = przth * br(ks+1)
      przth = przth * rzth
      dr(ks) = przth * ar(ks+2)
    END DO
    pp = pp * rfnu2
    IF (ias /= 1) THEN
      suma = up(lrp1)
      ju = lrp1
      DO  jr = 1, lr
        ju = ju - 1
        suma = suma + cr(jr) * up(ju)
      END DO
      asum = asum + suma
      asumr = REAL(asum, KIND=dp)
      asumi = AIMAG(asum)
      test = ABS(asumr) + ABS(asumi)
      IF (pp < tol .AND. test < tol) ias = 1
    END IF
    IF (ibs /= 1) THEN
      sumb = up(lr+2) + up(lrp1) * zc
      ju = lrp1
      DO  jr = 1, lr
        ju = ju - 1
        sumb = sumb + dr(jr) * up(ju)
      END DO
      bsum = bsum + sumb
      bsumr = REAL(bsum, KIND=dp)
      bsumi = AIMAG(bsum)
      test = ABS(bsumr) + ABS(bsumi)
      IF (pp < btol .AND. test < tol) ibs = 1
    END IF
    IF (ias == 1 .AND. ibs == 1) GO TO 170
  END DO
END IF

170 asum = asum + cone
bsum = -bsum * rfn13 / rtzta
GO TO 100
END SUBROUTINE cunhj



SUBROUTINE cseri(z, fnu, kode, n, y, nz, tol, elim, alim)
!***BEGIN PROLOGUE  CSERI
!***REFER TO  CBESI,CBESK

!     CSERI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE POWER SERIES FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) <= 2*SQRT(FNU+1). NZ=0 IS A NORMAL RETURN.
!     NZ > 0 MEANS THAT THE LAST NZ COMPONENTS WERE SET TO ZERO
!     DUE TO UNDERFLOW. NZ < 0 MEANS UNDERFLOW OCCURRED, BUT THE
!     CONDITION ABS(Z) <= 2*SQRT(FNU+1) WAS VIOLATED AND THE
!     COMPUTATION MUST BE COMPLETED IN ANOTHER ROUTINE WITH N=N-ABS(NZ).

!***ROUTINES CALLED  CUCHK,GAMLN,R1MACH
!***END PROLOGUE  CSERI

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: ak1, ck, coef, crsc, cz, hz, rz, s1, s2, w(2)
REAL (dp)     :: aa, acz, ak, arm, ascle, atol, az, dfnu,  &
                 fnup, rak1, rs, rtr1, s, ss, x
INTEGER       :: i, ib, iflag, il, k, l, m, nn, nw
REAL (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp)

nz = 0
az = ABS(z)
IF (az /= 0.0_dp) THEN
  x = REAL(z, KIND=dp)
  arm = 1.0D+3 * TINY(0.0_dp)
  rtr1 = SQRT(arm)
  crsc = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
  iflag = 0
  IF (az >= arm) THEN
    hz = z * 0.5_dp
    cz = czero
    IF (az > rtr1) cz = hz * hz
    acz = ABS(cz)
    nn = n
    ck = LOG(hz)

    10 dfnu = fnu + (nn-1)
    fnup = dfnu + 1.0_dp
!-----------------------------------------------------------------------
!     UNDERFLOW TEST
!-----------------------------------------------------------------------
    ak1 = ck * dfnu
    ak = gamln(fnup)
    ak1 = ak1 - ak
    IF (kode == 2) ak1 = ak1 - x
    rak1 = REAL(ak1, KIND=dp)
    IF (rak1 > -elim) GO TO 30

    20 nz = nz + 1
    y(nn) = czero
    IF (acz > dfnu) GO TO 120
    nn = nn - 1
    IF (nn == 0) RETURN
    GO TO 10

    30 IF (rak1 <= -alim) THEN
      iflag = 1
      ss = 1.0_dp / tol
      crsc = CMPLX(tol, 0.0_dp, KIND=dp)
      ascle = arm * ss
    END IF
    ak = AIMAG(ak1)
    aa = EXP(rak1)
    IF (iflag == 1) aa = aa * ss
    coef = aa * CMPLX(COS(ak), SIN(ak), KIND=dp)
    atol = tol * acz / fnup
    il = MIN(2,nn)
    DO  i = 1, il
      dfnu = fnu + (nn-i)
      fnup = dfnu + 1.0_dp
      s1 = cone
      IF (acz >= tol*fnup) THEN
        ak1 = cone
        ak = fnup + 2.0_dp
        s = fnup
        aa = 2.0_dp

        40 rs = 1.0_dp / s
        ak1 = ak1 * cz * rs
        s1 = s1 + ak1
        s = s + ak
        ak = ak + 2.0_dp
        aa = aa * acz * rs
        IF (aa > atol) GO TO 40
      END IF
      m = nn - i + 1
      s2 = s1 * coef
      w(i) = s2
      IF (iflag /= 0) THEN
        CALL cuchk(s2, nw, ascle, tol)
        IF (nw /= 0) GO TO 20
      END IF
      y(m) = s2 * crsc
      IF (i /= il) coef = coef * dfnu / hz
    END DO
    IF (nn <= 2) RETURN
    k = nn - 2
    ak = k
    rz = (cone+cone) / z
    IF (iflag == 1) GO TO 80
    ib = 3

    60 DO  i = ib, nn
      y(k) = CMPLX(ak+fnu, 0.0_dp, KIND=dp) * rz * y(k+1) + y(k+2)
      ak = ak - 1.0_dp
      k = k - 1
    END DO
    RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD WITH SCALED VALUES
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION ABOVE THE
!     UNDERFLOW LIMIT = ASCLE = TINY(0.0_dp)*CSCL*1.0E+3
!-----------------------------------------------------------------------
    80 s1 = w(1)
    s2 = w(2)
    DO  l = 3, nn
      ck = s2
      s2 = s1 + CMPLX(ak+fnu, 0.0_dp, KIND=dp) * rz * s2
      s1 = ck
      ck = s2 * crsc
      y(k) = ck
      ak = ak - 1.0_dp
      k = k - 1
      IF (ABS(ck) > ascle) GO TO 100
    END DO
    RETURN

    100 ib = l + 1
    IF (ib > nn) RETURN
    GO TO 60
  END IF
  nz = n
  IF (fnu == 0.0_dp) nz = nz - 1
END IF
y(1) = czero
IF (fnu == 0.0_dp) y(1) = cone
IF (n == 1) RETURN
y(2:n) = czero
RETURN
!-----------------------------------------------------------------------
!     RETURN WITH NZ < 0 IF ABS(Z*Z/4) > FNU+N-NZ-1 COMPLETE
!     THE CALCULATION IN CBINU WITH N=N-ABS(NZ)
!-----------------------------------------------------------------------
120 nz = -nz
RETURN
END SUBROUTINE cseri



SUBROUTINE casyi(z, fnu, kode, n, y, nz, rl, tol, elim, alim)
!***BEGIN PROLOGUE  CASYI
!***REFER TO  CBESI,CBESK

!     CASYI COMPUTES THE I BESSEL FUNCTION FOR REAL(Z) >= 0.0 BY
!     MEANS OF THE ASYMPTOTIC EXPANSION FOR LARGE ABS(Z) IN THE
!     REGION ABS(Z) > MAX(RL,FNU*FNU/2).  NZ=0 IS A NORMAL RETURN.
!     NZ < 0 INDICATES AN OVERFLOW ON KODE=1.

!***ROUTINES CALLED  R1MACH
!***END PROLOGUE  CASYI

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: rl
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: ak1, ck, cs1, cs2, cz, dk, ez, p1, rz, s2
REAL (dp)     :: aa, acz, aez, ak, arg, arm, atol, az, bb, bk, dfnu,  &
                 dnu2, fdn, rtr1, s, sgn, sqk, x, yy
INTEGER       :: i, ib, il, inu, j, jl, k, koded, m, nn

REAL (dp), PARAMETER     :: pi = 3.14159265358979324_dp, rtpi = 0.159154943091895336_dp
! rtpi = reciprocal of 2.pi
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp)

nz = 0
az = ABS(z)
x = REAL(z, KIND=dp)
arm = 1.0D+3 * TINY(0.0_dp)
rtr1 = SQRT(arm)
il = MIN(2,n)
dfnu = fnu + (n-il)
!-----------------------------------------------------------------------
!     OVERFLOW TEST
!-----------------------------------------------------------------------
ak1 = rtpi / z
ak1 = SQRT(ak1)
cz = z
IF (kode == 2) cz = z - x
acz = REAL(cz, KIND=dp)
IF (ABS(acz) <= elim) THEN
  dnu2 = dfnu + dfnu
  koded = 1
  IF (.NOT.(ABS(acz) > alim .AND. n > 2)) THEN
    koded = 0
    ak1 = ak1 * EXP(cz)
  END IF
  fdn = 0.0_dp
  IF (dnu2 > rtr1) fdn = dnu2 * dnu2
  ez = z * 8.0_dp
!-----------------------------------------------------------------------
!     WHEN Z IS IMAGINARY, THE ERROR TEST MUST BE MADE RELATIVE TO THE
!     FIRST RECIPROCAL POWER SINCE THIS IS THE LEADING TERM OF THE
!     EXPANSION FOR THE IMAGINARY PART.
!-----------------------------------------------------------------------
  aez = 8.0_dp * az
  s = tol / aez
  jl = rl + rl + 2
  yy = AIMAG(z)
  p1 = czero
  IF (yy /= 0.0_dp) THEN
!-----------------------------------------------------------------------
!     CALCULATE EXP(PI*(0.5+FNU+N-IL)*I) TO MINIMIZE LOSSES OF
!     SIGNIFICANCE WHEN FNU OR N IS LARGE
!-----------------------------------------------------------------------
    inu = fnu
    arg = (fnu - inu) * pi
    inu = inu + n - il
    ak = -SIN(arg)
    bk = COS(arg)
    IF (yy < 0.0_dp) bk = -bk
    p1 = CMPLX(ak, bk, KIND=dp)
    IF (MOD(inu,2) == 1) p1 = -p1
  END IF
  DO  k = 1, il
    sqk = fdn - 1.0_dp
    atol = s * ABS(sqk)
    sgn = 1.0_dp
    cs1 = cone
    cs2 = cone
    ck = cone
    ak = 0.0_dp
    aa = 1.0_dp
    bb = aez
    dk = ez
    DO  j = 1, jl
      ck = ck * sqk / dk
      cs2 = cs2 + ck
      sgn = -sgn
      cs1 = cs1 + ck * sgn
      dk = dk + ez
      aa = aa * ABS(sqk) / bb
      bb = bb + aez
      ak = ak + 8.0_dp
      sqk = sqk - ak
      IF (aa <= atol) GO TO 20
    END DO
    GO TO 60

    20 s2 = cs1
    IF (x+x < elim) s2 = s2 + p1 * cs2 * EXP(-z-z)
    fdn = fdn + 8.0_dp * dfnu + 4.0_dp
    p1 = -p1
    m = n - il + k
    y(m) = s2 * ak1
  END DO
  IF (n <= 2) RETURN
  nn = n
  k = nn - 2
  ak = k
  rz = (cone+cone) / z
  ib = 3
  DO  i = ib, nn
    y(k) = CMPLX(ak+fnu, 0.0_dp, KIND=dp) * rz * y(k+1) + y(k+2)
    ak = ak - 1.0_dp
    k = k - 1
  END DO
  IF (koded == 0) RETURN
  ck = EXP(cz)
  y(1:nn) = y(1:nn) * ck
  RETURN
END IF
nz = -1
RETURN

60 nz = -2
RETURN
END SUBROUTINE casyi



SUBROUTINE cbunk(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
!***BEGIN PROLOGUE  CBUNK
!***REFER TO  CBESK,CBESH

!     CBUNK COMPUTES THE K BESSEL FUNCTION FOR FNU > FNUL.
!     ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION FOR K(FNU,Z)
!     IN CUNK1 AND THE EXPANSION FOR H(2,FNU,Z) IN CUNK2

!***ROUTINES CALLED  CUNK1,CUNK2
!***END PROLOGUE  CBUNK

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: mr
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN OUT)  :: tol
REAL (dp), INTENT(IN OUT)  :: elim
REAL (dp), INTENT(IN OUT)  :: alim

REAL (dp)  :: ax, ay, xx, yy

nz = 0
xx = REAL(z, KIND=dp)
yy = AIMAG(z)
ax = ABS(xx) * 1.7321_dp
ay = ABS(yy)
IF (ay <= ax) THEN
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR K(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  CALL cunk1(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
ELSE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR H(2,FNU,Z*EXP(M*HPI)) FOR LARGE FNU
!     APPLIED IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I AND HPI=PI/2
!-----------------------------------------------------------------------
  CALL cunk2(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
END IF
RETURN
END SUBROUTINE cbunk



SUBROUTINE cunk1(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
!***BEGIN PROLOGUE  CUNK1
!***REFER TO  CBESK

!     CUNK1 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE
!     RIGHT HALF PLANE TO THE LEFT HALF PLANE BY MEANS OF THE
!     UNIFORM ASYMPTOTIC EXPANSION.
!     MR INDICATES THE DIRECTION OF ROTATION FOR ANALYTIC CONTINUATION.
!     NZ=-1 MEANS AN OVERFLOW WILL OCCUR

!***ROUTINES CALLED  CS1S2,CUCHK,CUNIK,R1MACH
!***END PROLOGUE  CUNK1

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: mr
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: cfn, ck, crsc, cs, cscl, csgn, cspn, csr(3), css(3),  &
                 cwrk(16,3), cy(2), c1, c2, phi(2), rz, sum(2), s1, s2, &
                 zeta1(2), zeta2(2), zr, phid, zeta1d, zeta2d, sumd
REAL (dp)     :: ang, aphi, asc, ascle, bry(3), cpn, c2i, c2m, c2r,  &
                 fmr, fn, fnf, rs1, sgn, spn, x
INTEGER       :: i, ib, iflag, ifn, il, init(2), inu, iuf, k, kdflg, kflag, &
                 kk, m, nw, j, ipard, initd, ic
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp)
REAL (dp), PARAMETER     :: pi = 3.14159265358979324_dp

kdflg = 1
nz = 0
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!-----------------------------------------------------------------------
cscl = 1.0_dp/tol
crsc = tol
css(1) = cscl
css(2) = cone
css(3) = crsc
csr(1) = crsc
csr(2) = cone
csr(3) = cscl
bry(1) = 1.0E+3 * TINY(0.0_dp) / tol
bry(2) = 1.0_dp / bry(1)
bry(3) = HUGE(0.0_dp)
x = REAL(z, KIND=dp)
zr = z
IF (x < 0.0_dp) zr = -z
j = 2
DO  i = 1, n
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
  j = 3 - j
  fn = fnu + (i-1)
  init(j) = 0
  CALL cunik(zr, fn, 2, 0, tol, init(j), phi(j), zeta1(j), zeta2(j), sum(j), &
             cwrk(1:,j))
  IF (kode /= 1) THEN
    cfn = fn
    s1 = zeta1(j) - cfn * (cfn/(zr + zeta2(j)))
  ELSE
    s1 = zeta1(j) - zeta2(j)
  END IF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) <= elim) THEN
    IF (kdflg == 1) kflag = 2
    IF (ABS(rs1) >= alim) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
      aphi = ABS(phi(j))
      rs1 = rs1 + LOG(aphi)
      IF (ABS(rs1) > elim) GO TO 10
      IF (kdflg == 1) kflag = 1
      IF (rs1 >= 0.0_dp) THEN
        IF (kdflg == 1) kflag = 3
      END IF
    END IF
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    s2 = phi(j) * sum(j)
    c2r = REAL(s1, KIND=dp)
    c2i = AIMAG(s1)
    c2m = EXP(c2r) * REAL(css(kflag), KIND=dp)
    s1 = c2m * CMPLX(COS(c2i), SIN(c2i), KIND=dp)
    s2 = s2 * s1
    IF (kflag == 1) THEN
      CALL cuchk(s2, nw, bry(1), tol)
      IF (nw /= 0) GO TO 10
    END IF
    cy(kdflg) = s2
    y(i) = s2 * csr(kflag)
    IF (kdflg == 2) GO TO 30
    kdflg = 2
    CYCLE
  END IF

  10 IF (rs1 > 0.0_dp) GO TO 150
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  IF (x < 0.0_dp) GO TO 150
  kdflg = 1
  y(i) = czero
  nz = nz + 1
  IF (i /= 1) THEN
    IF (y(i-1) /= czero) THEN
      y(i-1) = czero
      nz = nz + 1
    END IF
  END IF
END DO
i = n

30 rz = 2.0_dp / zr
ck = fn * rz
ib = i + 1
IF (n >= ib) THEN
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
!     ON UNDERFLOW
!-----------------------------------------------------------------------
  fn = fnu + (n-1)
  ipard = 1
  IF (mr /= 0) ipard = 0
  initd = 0
  CALL cunik(zr, fn, 2, ipard, tol, initd, phid, zeta1d, zeta2d, sumd,  &
             cwrk(1:,3))
  IF (kode /= 1) THEN
    cfn = fn
    s1 = zeta1d - cfn * (cfn/(zr + zeta2d))
  ELSE
    s1 = zeta1d - zeta2d
  END IF
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) <= elim) THEN
    IF (ABS(rs1) < alim) GO TO 50
!-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!-----------------------------------------------------------------------
    aphi = ABS(phid)
    rs1 = rs1 + LOG(aphi)
    IF (ABS(rs1) < elim) GO TO 50
  END IF
  IF (rs1 > 0.0_dp) GO TO 150
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  IF (x < 0.0_dp) GO TO 150
  nz = n
  y(1:n) = czero
  RETURN
!-----------------------------------------------------------------------
!     RECUR FORWARD FOR REMAINDER OF THE SEQUENCE
!-----------------------------------------------------------------------
  50 s1 = cy(1)
  s2 = cy(2)
  c1 = csr(kflag)
  ascle = bry(kflag)
  DO  i = ib, n
    c2 = s2
    s2 = ck * s2 + s1
    s1 = c2
    ck = ck + rz
    c2 = s2 * c1
    y(i) = c2
    IF (kflag < 3) THEN
      c2r = REAL(c2, KIND=dp)
      c2i = AIMAG(c2)
      c2r = ABS(c2r)
      c2i = ABS(c2i)
      c2m = MAX(c2r,c2i)
      IF (c2m > ascle) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1 * c1
        s2 = c2
        s1 = s1 * css(kflag)
        s2 = s2 * css(kflag)
        c1 = csr(kflag)
      END IF
    END IF
  END DO
END IF
IF (mr == 0) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z) < 0.0_dp
!-----------------------------------------------------------------------
nz = 0
fmr = mr
sgn = -SIGN(pi, fmr)
!-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCIONS RESP.
!-----------------------------------------------------------------------
csgn = CMPLX(0.0_dp, sgn, KIND=dp)
inu = fnu
fnf = fnu - inu
ifn = inu + n - 1
ang = fnf * sgn
cpn = COS(ang)
spn = SIN(ang)
cspn = CMPLX(cpn, spn, KIND=dp)
IF (MOD(ifn,2) == 1) cspn = -cspn
asc = bry(1)
kk = n
iuf = 0
kdflg = 1
ib = ib - 1
ic = ib - 1
DO  k = 1, n
  fn = fnu + (kk-1)
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
  m = 3
  IF (n > 2) GO TO 80

  70 initd = init(j)
  phid = phi(j)
  zeta1d = zeta1(j)
  zeta2d = zeta2(j)
  sumd = sum(j)
  m = j
  j = 3 - j
  GO TO 90

  80 IF (.NOT.(kk == n .AND. ib < n)) THEN
    IF (kk == ib .OR. kk == ic) GO TO 70
    initd = 0
  END IF

  90 CALL cunik(zr, fn, 1, 0, tol, initd, phid, zeta1d, zeta2d, sumd,  &
                cwrk(1:,m))
  IF (kode /= 1) THEN
    cfn = fn
    s1 = -zeta1d + cfn * (cfn/(zr + zeta2d))
  ELSE
    s1 = -zeta1d + zeta2d
  END IF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) > elim) GO TO 110
  IF (kdflg == 1) iflag = 2
  IF (ABS(rs1) >= alim) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    aphi = ABS(phid)
    rs1 = rs1 + LOG(aphi)
    IF (ABS(rs1) > elim) GO TO 110
    IF (kdflg == 1) iflag = 1
    IF (rs1 >= 0.0_dp) THEN
      IF (kdflg == 1) iflag = 3
    END IF
  END IF
  s2 = csgn * phid * sumd
  c2r = REAL(s1, KIND=dp)
  c2i = AIMAG(s1)
  c2m = EXP(c2r) * REAL(css(iflag), KIND=dp)
  s1 = c2m * CMPLX(COS(c2i), SIN(c2i), KIND=dp)
  s2 = s2 * s1
  IF (iflag == 1) THEN
    CALL cuchk(s2, nw, bry(1), tol)
    IF (nw /= 0) s2 = 0.0_dp
  END IF

  100 cy(kdflg) = s2
  c2 = s2
  s2 = s2 * csr(iflag)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
  s1 = y(kk)
  IF (kode /= 1) THEN
    CALL cs1s2(zr, s1, s2, nw, asc, alim, iuf)
    nz = nz + nw
  END IF
  y(kk) = s1 * cspn + s2
  kk = kk - 1
  cspn = -cspn
  IF (c2 == czero) THEN
    kdflg = 1
    CYCLE
  END IF
  IF (kdflg == 2) GO TO 130
  kdflg = 2
  CYCLE

  110 IF (rs1 > 0.0_dp) GO TO 150
  s2 = czero
  GO TO 100
END DO
k = n

130 il = n - k
IF (il == 0) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!-----------------------------------------------------------------------
s1 = cy(1)
s2 = cy(2)
cs = csr(iflag)
ascle = bry(iflag)
fn = inu + il
DO  i = 1, il
  c2 = s2
  s2 = s1 + (fn + fnf) * rz * s2
  s1 = c2
  fn = fn - 1.0_dp
  c2 = s2 * cs
  ck = c2
  c1 = y(kk)
  IF (kode /= 1) THEN
    CALL cs1s2(zr, c1, c2, nw, asc, alim, iuf)
    nz = nz + nw
  END IF
  y(kk) = c1 * cspn + c2
  kk = kk - 1
  cspn = -cspn
  IF (iflag < 3) THEN
    c2r = REAL(ck, KIND=dp)
    c2i = AIMAG(ck)
    c2r = ABS(c2r)
    c2i = ABS(c2i)
    c2m = MAX(c2r, c2i)
    IF (c2m > ascle) THEN
      iflag = iflag + 1
      ascle = bry(iflag)
      s1 = s1 * cs
      s2 = ck
      s1 = s1 * css(iflag)
      s2 = s2 * css(iflag)
      cs = csr(iflag)
    END IF
  END IF
END DO
RETURN

150 nz = -1
RETURN
END SUBROUTINE cunk1



SUBROUTINE cunk2(z, fnu, kode, mr, n, y, nz, tol, elim, alim)
!***BEGIN PROLOGUE  CUNK2
!***REFER TO  CBESK

!  CUNK2 COMPUTES K(FNU,Z) AND ITS ANALYTIC CONTINUATION FROM THE RIGHT HALF
!  PLANE TO THE LEFT HALF PLANE BY MEANS OF THE UNIFORM ASYMPTOTIC EXPANSIONS
!  FOR H(KIND,FNU,ZN) AND J(FNU,ZN) WHERE ZN IS IN THE RIGHT HALF PLANE,
!  KIND=(3-MR)/2, MR=+1 OR -1.
!  HERE ZN=ZR*I OR -ZR*I WHERE ZR=Z IF Z IS IN THE RIGHT HALF PLANE OR ZR=-Z
!  IF Z IS IN THE LEFT HALF PLANE.  MR INDICATES THE DIRECTION OF ROTATION FOR
!  ANALYTIC CONTINUATION.
!  NZ=-1 MEANS AN OVERFLOW WILL OCCUR

!***ROUTINES CALLED  CAIRY,CS1S2,CUCHK,CUNHJ,R1MACH
!***END PROLOGUE  CUNK2

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: mr
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: ai, arg(2), asum(2), bsum(2), cfn, ck, cs, csgn, cspn,  &
                 csr(3), css(3), cy(2), c1, c2, dai, phi(2), rz, s1, s2, &
                 zb, zeta1(2), zeta2(2), zn, zr, phid, argd, zeta1d, zeta2d, &
                 asumd, bsumd
REAL (dp)     :: aarg, ang, aphi, asc, ascle, bry(3), car, cpn, c2i,  &
                 c2m, c2r, crsc, cscl, fmr, fn, fnf, rs1, sar, sgn, spn, x, yy
INTEGER       :: i, ib, iflag, ifn, il, in, inu, iuf, k, kdflg, kflag, kk,  &
                 nai, ndai, nw, idum, j, ipard, ic
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp),  &
                            ci = (0.0_dp,1.0_dp),   &
                            cr1 = (1.0_dp, 1.73205080756887729_dp),  &
                            cr2 = (-0.5_dp, -8.66025403784438647D-01)
REAL (dp), PARAMETER     :: hpi = 1.57079632679489662_dp,  &
                            pi = 3.14159265358979324_dp,  &
                            aic = 1.26551212348464539_dp
COMPLEX (dp), PARAMETER  :: cip(4) = (/ (1.0_dp,0.0_dp), (0.0_dp,-1.0_dp),  &
                                        (-1.0_dp,0.0_dp), (0.0_dp,1.0_dp) /)

kdflg = 1
nz = 0
!-----------------------------------------------------------------------
!     EXP(-ALIM)=EXP(-ELIM)/TOL=APPROX. ONE PRECISION GREATER THAN
!     THE UNDERFLOW LIMIT
!-----------------------------------------------------------------------
cscl = 1.0_dp/tol
crsc = tol
css(1) = cscl
css(2) = cone
css(3) = crsc
csr(1) = crsc
csr(2) = cone
csr(3) = cscl
bry(1) = 1.0E+3 * TINY(0.0_dp) / tol
bry(2) = 1.0_dp / bry(1)
bry(3) = HUGE(0.0_dp)
x = REAL(z, KIND=dp)
zr = z
IF (x < 0.0_dp) zr = -z
yy = AIMAG(zr)
zn = -zr * ci
zb = zr
inu = fnu
fnf = fnu - inu
ang = -hpi * fnf
car = COS(ang)
sar = SIN(ang)
cpn = -hpi * car
spn = -hpi * sar
c2 = CMPLX(-spn, cpn, KIND=dp)
kk = MOD(inu,4) + 1
cs = cr1 * c2 * cip(kk)
IF (yy <= 0.0_dp) THEN
  zn = CONJG(-zn)
  zb = CONJG(zb)
END IF
!-----------------------------------------------------------------------
!     K(FNU,Z) IS COMPUTED FROM H(2,FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT.  FOURTH QUADRANT VALUES (YY <= 0.0_dp) ARE COMPUTED BY
!     CONJUGATION SINCE THE K FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!-----------------------------------------------------------------------
j = 2
DO  i = 1, n
!-----------------------------------------------------------------------
!     J FLIP FLOPS BETWEEN 1 AND 2 IN J = 3 - J
!-----------------------------------------------------------------------
  j = 3 - j
  fn = fnu + (i-1)
  CALL cunhj(zn, fn, 0, tol, phi(j), arg(j), zeta1(j), zeta2(j), asum(j), &
             bsum(j))
  IF (kode /= 1) THEN
    cfn = CMPLX(fn, 0.0_dp, KIND=dp)
    s1 = zeta1(j) - cfn * (cfn/(zb + zeta2(j)))
  ELSE
    s1 = zeta1(j) - zeta2(j)
  END IF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) <= elim) THEN
    IF (kdflg == 1) kflag = 2
    IF (ABS(rs1) >= alim) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
      aphi = ABS(phi(j))
      aarg = ABS(arg(j))
      rs1 = rs1 + LOG(aphi) - 0.25_dp * LOG(aarg) - aic
      IF (ABS(rs1) > elim) GO TO 10
      IF (kdflg == 1) kflag = 1
      IF (rs1 >= 0.0_dp) THEN
        IF (kdflg == 1) kflag = 3
      END IF
    END IF
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
    c2 = arg(j) * cr2
    CALL cairy(c2, 0, 2, ai, nai, idum)
    CALL cairy(c2, 1, 2, dai, ndai, idum)
    s2 = cs * phi(j) * (ai*asum(j) + cr2*dai*bsum(j))
    c2r = REAL(s1, KIND=dp)
    c2i = AIMAG(s1)
    c2m = EXP(c2r) * REAL(css(kflag), KIND=dp)
    s1 = c2m * CMPLX(COS(c2i), SIN(c2i), KIND=dp)
    s2 = s2 * s1
    IF (kflag == 1) THEN
      CALL cuchk(s2, nw, bry(1), tol)
      IF (nw /= 0) GO TO 10
    END IF
    IF (yy <= 0.0_dp) s2 = CONJG(s2)
    cy(kdflg) = s2
    y(i) = s2 * csr(kflag)
    cs = -ci * cs
    IF (kdflg == 2) GO TO 30
    kdflg = 2
    CYCLE
  END IF

  10 IF (rs1 > 0.0_dp) GO TO 150
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  IF (x < 0.0_dp) GO TO 150
  kdflg = 1
  y(i) = czero
  cs = -ci * cs
  nz = nz + 1
  IF (i /= 1) THEN
    IF (y(i-1) /= czero) THEN
      y(i-1) = czero
      nz = nz + 1
    END IF
  END IF
END DO
i = n

30 rz = 2.0_dp / zr
ck = fn * rz
ib = i + 1
IF (n >= ib) THEN
!-----------------------------------------------------------------------
!     TEST LAST MEMBER FOR UNDERFLOW AND OVERFLOW, SET SEQUENCE TO ZERO
!     ON UNDERFLOW
!-----------------------------------------------------------------------
  fn = fnu + (n-1)
  ipard = 1
  IF (mr /= 0) ipard = 0
  CALL cunhj(zn, fn, ipard, tol, phid, argd, zeta1d, zeta2d, asumd, bsumd)
  IF (kode /= 1) THEN
    cfn = fn
    s1 = zeta1d - cfn * (cfn/(zb + zeta2d))
  ELSE
    s1 = zeta1d - zeta2d
  END IF
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) <= elim) THEN
    IF (ABS(rs1) < alim) GO TO 50
!-----------------------------------------------------------------------
!     REFINE ESTIMATE AND TEST
!-----------------------------------------------------------------------
    aphi = ABS(phid)
    aarg = ABS(argd)
    rs1 = rs1 + LOG(aphi) - 0.25_dp * LOG(aarg) - aic
    IF (ABS(rs1) < elim) GO TO 50
  END IF
  IF (rs1 > 0.0_dp) GO TO 150
!-----------------------------------------------------------------------
!     FOR X < 0.0, THE I FUNCTION TO BE ADDED WILL OVERFLOW
!-----------------------------------------------------------------------
  IF (x < 0.0_dp) GO TO 150
  nz = n
  y(1:n) = czero
  RETURN
!-----------------------------------------------------------------------
!     SCALED FORWARD RECURRENCE FOR REMAINDER OF THE SEQUENCE
!-----------------------------------------------------------------------
  50 s1 = cy(1)
  s2 = cy(2)
  c1 = csr(kflag)
  ascle = bry(kflag)
  DO  i = ib, n
    c2 = s2
    s2 = ck * s2 + s1
    s1 = c2
    ck = ck + rz
    c2 = s2 * c1
    y(i) = c2
    IF (kflag < 3) THEN
      c2r = REAL(c2, KIND=dp)
      c2i = AIMAG(c2)
      c2r = ABS(c2r)
      c2i = ABS(c2i)
      c2m = MAX(c2r,c2i)
      IF (c2m > ascle) THEN
        kflag = kflag + 1
        ascle = bry(kflag)
        s1 = s1 * c1
        s2 = c2
        s1 = s1 * css(kflag)
        s2 = s2 * css(kflag)
        c1 = csr(kflag)
      END IF
    END IF
  END DO
END IF
IF (mr == 0) RETURN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION FOR RE(Z) < 0.0_dp
!-----------------------------------------------------------------------
nz = 0
fmr = mr
sgn = -SIGN(pi, fmr)
!-----------------------------------------------------------------------
!     CSPN AND CSGN ARE COEFF OF K AND I FUNCTIONS RESP.
!-----------------------------------------------------------------------
csgn = CMPLX(0.0_dp, sgn, KIND=dp)
IF (yy <= 0.0_dp) csgn = CONJG(csgn)
ifn = inu + n - 1
ang = fnf * sgn
cpn = COS(ang)
spn = SIN(ang)
cspn = CMPLX(cpn, spn, KIND=dp)
IF (MOD(ifn,2) == 1) cspn = -cspn
!-----------------------------------------------------------------------
!     CS=COEFF OF THE J FUNCTION TO GET THE I FUNCTION.  I(FNU,Z) IS
!     COMPUTED FROM EXP(I*FNU*HPI)*J(FNU,-I*Z) WHERE Z IS IN THE FIRST
!     QUADRANT.  FOURTH QUADRANT VALUES (YY <= 0.0_dp) ARE COMPUTED BY
!     CONJUGATION SINCE THE I FUNCTION IS REAL ON THE POSITIVE REAL AXIS
!-----------------------------------------------------------------------
cs = CMPLX(car, -sar, KIND=dp) * csgn
in = MOD(ifn,4) + 1
c2 = cip(in)
cs = cs * CONJG(c2)
asc = bry(1)
kk = n
kdflg = 1
ib = ib - 1
ic = ib - 1
iuf = 0
DO  k = 1, n
!-----------------------------------------------------------------------
!     LOGIC TO SORT OUT CASES WHOSE PARAMETERS WERE SET FOR THE K
!     FUNCTION ABOVE
!-----------------------------------------------------------------------
  fn = fnu + (kk-1)
  IF (n > 2) GO TO 80

  70 phid = phi(j)
  argd = arg(j)
  zeta1d = zeta1(j)
  zeta2d = zeta2(j)
  asumd = asum(j)
  bsumd = bsum(j)
  j = 3 - j
  GO TO 90

  80 IF (.NOT.(kk == n .AND. ib < n)) THEN
    IF (kk == ib .OR. kk == ic) GO TO 70
    CALL cunhj(zn, fn, 0, tol, phid, argd, zeta1d, zeta2d, asumd, bsumd)
  END IF

  90 IF (kode /= 1) THEN
    cfn = fn
    s1 = -zeta1d + cfn * (cfn/(zb + zeta2d))
  ELSE
    s1 = -zeta1d + zeta2d
  END IF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) > elim) GO TO 110
  IF (kdflg == 1) iflag = 2
  IF (ABS(rs1) >= alim) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    aphi = ABS(phid)
    aarg = ABS(argd)
    rs1 = rs1 + LOG(aphi) - 0.25_dp * LOG(aarg) - aic
    IF (ABS(rs1) > elim) GO TO 110
    IF (kdflg == 1) iflag = 1
    IF (rs1 >= 0.0_dp) THEN
      IF (kdflg == 1) iflag = 3
    END IF
  END IF
  CALL cairy(argd, 0, 2, ai, nai, idum)
  CALL cairy(argd, 1, 2, dai, ndai, idum)
  s2 = cs * phid * (ai*asumd + dai*bsumd)
  c2r = REAL(s1, KIND=dp)
  c2i = AIMAG(s1)
  c2m = EXP(c2r) * REAL(css(iflag), KIND=dp)
  s1 = c2m * CMPLX(COS(c2i), SIN(c2i), KIND=dp)
  s2 = s2 * s1
  IF (iflag == 1) THEN
    CALL cuchk(s2, nw, bry(1), tol)
    IF (nw /= 0) s2 = 0.0_dp
  END IF

  100 IF (yy <= 0.0_dp) s2 = CONJG(s2)
  cy(kdflg) = s2
  c2 = s2
  s2 = s2 * csr(iflag)
!-----------------------------------------------------------------------
!     ADD I AND K FUNCTIONS, K SEQUENCE IN Y(I), I=1,N
!-----------------------------------------------------------------------
  s1 = y(kk)
  IF (kode /= 1) THEN
    CALL cs1s2(zr, s1, s2, nw, asc, alim, iuf)
    nz = nz + nw
  END IF
  y(kk) = s1 * cspn + s2
  kk = kk - 1
  cspn = -cspn
  cs = -cs * ci
  IF (c2 == czero) THEN
    kdflg = 1
    CYCLE
  END IF
  IF (kdflg == 2) GO TO 130
  kdflg = 2
  CYCLE

  110 IF (rs1 > 0.0_dp) GO TO 150
  s2 = czero
  GO TO 100
END DO
k = n

130 il = n - k
IF (il == 0) RETURN
!-----------------------------------------------------------------------
!     RECUR BACKWARD FOR REMAINDER OF I SEQUENCE AND ADD IN THE
!     K FUNCTIONS, SCALING THE I SEQUENCE DURING RECURRENCE TO KEEP
!     INTERMEDIATE ARITHMETIC ON SCALE NEAR EXPONENT EXTREMES.
!-----------------------------------------------------------------------
s1 = cy(1)
s2 = cy(2)
cs = csr(iflag)
ascle = bry(iflag)
fn = inu + il
DO  i = 1, il
  c2 = s2
  s2 = s1 + CMPLX(fn+fnf, 0.0_dp, KIND=dp) * rz * s2
  s1 = c2
  fn = fn - 1.0_dp
  c2 = s2 * cs
  ck = c2
  c1 = y(kk)
  IF (kode /= 1) THEN
    CALL cs1s2(zr, c1, c2, nw, asc, alim, iuf)
    nz = nz + nw
  END IF
  y(kk) = c1 * cspn + c2
  kk = kk - 1
  cspn = -cspn
  IF (iflag < 3) THEN
    c2r = REAL(ck, KIND=dp)
    c2i = AIMAG(ck)
    c2r = ABS(c2r)
    c2i = ABS(c2i)
    c2m = MAX(c2r, c2i)
    IF (c2m > ascle) THEN
      iflag = iflag + 1
      ascle = bry(iflag)
      s1 = s1 * cs
      s2 = ck
      s1 = s1 * css(iflag)
      s2 = s2 * css(iflag)
      cs = csr(iflag)
    END IF
  END IF
END DO
RETURN

150 nz = -1
RETURN
END SUBROUTINE cunk2



SUBROUTINE cbuni(z, fnu, kode, n, y, nz, nui, nlast, fnul, tol, elim, alim)
!***BEGIN PROLOGUE  CBUNI
!***REFER TO  CBESI,CBESK

!   CBUNI COMPUTES THE I BESSEL FUNCTION FOR LARGE ABS(Z) > FNUL AND
!   FNU+N-1 < FNUL.  THE ORDER IS INCREASED FROM FNU+N-1 GREATER THAN FNUL
!   BY ADDING NUI AND COMPUTING ACCORDING TO THE UNIFORM ASYMPTOTIC EXPANSION
!   FOR I(FNU,Z) ON IFORM=1 AND THE EXPANSION FOR J(FNU,Z) ON IFORM=2

!***ROUTINES CALLED  CUNI1,CUNI2,R1MACH
!***END PROLOGUE  CBUNI

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(IN)        :: nui
INTEGER, INTENT(OUT)       :: nlast
REAL (dp), INTENT(IN)      :: fnul
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: cscl, cscr, cy(2), rz, st, s1, s2
REAL (dp)     :: ax, ay, dfnu, fnui, gnu, xx, yy, ascle, bry(3), str, sti, &
                 stm
INTEGER       :: i, iflag, iform, k, nl, nw

nz = 0
xx = REAL(z, KIND=dp)
yy = AIMAG(z)
ax = ABS(xx) * 1.73205080756887_dp
ay = ABS(yy)
iform = 1
IF (ay > ax) iform = 2
IF (nui == 0) GO TO 40
fnui = nui
dfnu = fnu + (n-1)
gnu = dfnu + fnui
IF (iform /= 2) THEN
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  CALL cuni1(z, gnu, kode, 2, cy, nw, nlast, fnul, tol, elim, alim)
ELSE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU APPLIED
!     IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I AND HPI=PI/2
!-----------------------------------------------------------------------
  CALL cuni2(z, gnu, kode, 2, cy, nw, nlast, fnul, tol, elim, alim)
END IF
IF (nw >= 0) THEN
  IF (nw /= 0) GO TO 50
  ay = ABS(cy(1))
!----------------------------------------------------------------------
!     SCALE BACKWARD RECURRENCE, BRY(3) IS DEFINED BUT NEVER USED
!----------------------------------------------------------------------
  bry(1) = 1.0E+3 * TINY(0.0_dp) / tol
  bry(2) = 1.0_dp / bry(1)
  bry(3) = bry(2)
  iflag = 2
  ascle = bry(2)
  ax = 1.0_dp
  cscl = ax
  IF (ay <= bry(1)) THEN
    iflag = 1
    ascle = bry(1)
    ax = 1.0_dp / tol
    cscl = ax
  ELSE
    IF (ay >= bry(2)) THEN
      iflag = 3
      ascle = bry(3)
      ax = tol
      cscl = ax
    END IF
  END IF
  ay = 1.0_dp / ax
  cscr = ay
  s1 = cy(2) * cscl
  s2 = cy(1) * cscl
  rz = 2.0_dp / z
  DO  i = 1, nui
    st = s2
    s2 = CMPLX(dfnu+fnui, 0.0_dp, KIND=dp) * rz * s2 + s1
    s1 = st
    fnui = fnui - 1.0_dp
    IF (iflag < 3) THEN
      st = s2 * cscr
      str = REAL(st, KIND=dp)
      sti = AIMAG(st)
      str = ABS(str)
      sti = ABS(sti)
      stm = MAX(str,sti)
      IF (stm > ascle) THEN
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1 * cscr
        s2 = st
        ax = ax * tol
        ay = 1.0_dp / ax
        cscl = ax
        cscr = ay
        s1 = s1 * cscl
        s2 = s2 * cscl
      END IF
    END IF
  END DO
  y(n) = s2 * cscr
  IF (n == 1) RETURN
  nl = n - 1
  fnui = nl
  k = nl
  DO  i = 1, nl
    st = s2
    s2 = CMPLX(fnu+fnui, 0.0_dp, KIND=dp) * rz * s2 + s1
    s1 = st
    st = s2 * cscr
    y(k) = st
    fnui = fnui - 1.0_dp
    k = k - 1
    IF (iflag < 3) THEN
      str = REAL(st, KIND=dp)
      sti = AIMAG(st)
      str = ABS(str)
      sti = ABS(sti)
      stm = MAX(str,sti)
      IF (stm > ascle) THEN
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1 * cscr
        s2 = st
        ax = ax * tol
        ay = 1.0_dp / ax
        cscl = ax
        cscr = ay
        s1 = s1 * cscl
        s2 = s2 * cscl
      END IF
    END IF
  END DO
  RETURN
END IF

30 nz = -1
IF (nw == -2) nz = -2
RETURN

40 IF (iform /= 2) THEN
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR I(FNU,Z) FOR LARGE FNU APPLIED IN
!     -PI/3 <= ARG(Z) <= PI/3
!-----------------------------------------------------------------------
  CALL cuni1(z, fnu, kode, n, y, nw, nlast, fnul, tol, elim, alim)
ELSE
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR J(FNU,Z*EXP(M*HPI)) FOR LARGE FNU APPLIED
!     IN PI/3 < ABS(ARG(Z)) <= PI/2 WHERE M=+I OR -I AND HPI=PI/2
!-----------------------------------------------------------------------
  CALL cuni2(z, fnu, kode, n, y, nw, nlast, fnul, tol, elim, alim)
END IF
IF (nw < 0) GO TO 30
nz = nw
RETURN

50 nlast = n
RETURN
END SUBROUTINE cbuni



SUBROUTINE cuni1(z, fnu, kode, n, y, nz, nlast, fnul, tol, elim, alim)
!***BEGIN PROLOGUE  CUNI1
!***REFER TO  CBESI,CBESK

!     CUNI1 COMPUTES I(FNU,Z)  BY MEANS OF THE UNIFORM ASYMPTOTIC
!     EXPANSION FOR I(FNU,Z) IN -PI/3 <= ARG Z <= PI/3.

!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION. NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I)=CZERO FOR I=NLAST+1,N

!***ROUTINES CALLED  CUCHK,CUNIK,CUOIK,R1MACH
!***END PROLOGUE  CUNI1

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: nlast
REAL (dp), INTENT(IN)      :: fnul
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: cfn, crsc, cscl, csr(3), css(3), c1, c2, cwrk(16), phi,  &
                 rz, sum, s1, s2, zeta1, zeta2, cy(2)
REAL (dp)     :: aphi, ascle, bry(3), c2i, c2m, c2r, fn, rs1, yy
INTEGER       :: i, iflag, init, k, m, nd, nn, nuf, nw
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp)

nz = 0
nd = n
nlast = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAGNITUDE
!     ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM) = EXP(ELIM)*TOL
!-----------------------------------------------------------------------
cscl = CMPLX(1.0_dp/tol, 0.0_dp, KIND=dp)
crsc = tol
css(1) = cscl
css(2) = cone
css(3) = crsc
csr(1) = crsc
csr(2) = cone
csr(3) = cscl
bry(1) = 1.0E+3 * TINY(0.0_dp) / tol
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
fn = MAX(fnu, 1.0_dp)
init = 0
CALL cunik(z, fn, 1, 1, tol, init, phi, zeta1, zeta2, sum, cwrk)
IF (kode /= 1) THEN
  cfn = fn
  s1 = -zeta1 + cfn * (cfn/(z + zeta2))
ELSE
  s1 = -zeta1 + zeta2
END IF
rs1 = REAL(s1, KIND=dp)
IF (ABS(rs1) > elim) GO TO 70

10 nn = MIN(2,nd)
DO  i = 1, nn
  fn = fnu + (nd-i)
  init = 0
  CALL cunik(z, fn, 1, 0, tol, init, phi, zeta1, zeta2, sum, cwrk)
  IF (kode /= 1) THEN
    cfn = fn
    yy = AIMAG(z)
    s1 = -zeta1 + cfn * (cfn/(z+zeta2)) + CMPLX(0.0_dp, yy, KIND=dp)
  ELSE
    s1 = -zeta1 + zeta2
  END IF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) > elim) GO TO 50
  IF (i == 1) iflag = 2
  IF (ABS(rs1) >= alim) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    aphi = ABS(phi)
    rs1 = rs1 + LOG(aphi)
    IF (ABS(rs1) > elim) GO TO 50
    IF (i == 1) iflag = 1
    IF (rs1 >= 0.0_dp) THEN
      IF (i == 1) iflag = 3
    END IF
  END IF
!-----------------------------------------------------------------------
!     SCALE S1 IF ABS(S1) < ASCLE
!-----------------------------------------------------------------------
  s2 = phi * sum
  c2r = REAL(s1, KIND=dp)
  c2i = AIMAG(s1)
  c2m = EXP(c2r) * REAL(css(iflag), KIND=dp)
  s1 = c2m * CMPLX(COS(c2i), SIN(c2i), KIND=dp)
  s2 = s2 * s1
  IF (iflag == 1) THEN
    CALL cuchk(s2, nw, bry(1), tol)
    IF (nw /= 0) GO TO 50
  END IF
  m = nd - i + 1
  cy(i) = s2
  y(m) = s2 * csr(iflag)
END DO
IF (nd > 2) THEN
  rz = 2.0_dp / z
  bry(2) = 1.0_dp / bry(1)
  bry(3) = HUGE(0.0_dp)
  s1 = cy(1)
  s2 = cy(2)
  c1 = csr(iflag)
  ascle = bry(iflag)
  k = nd - 2
  fn = k
  DO  i = 3, nd
    c2 = s2
    s2 = s1 + CMPLX(fnu+fn, 0.0_dp, KIND=dp) * rz * s2
    s1 = c2
    c2 = s2 * c1
    y(k) = c2
    k = k - 1
    fn = fn - 1.0_dp
    IF (iflag < 3) THEN
      c2r = REAL(c2, KIND=dp)
      c2i = AIMAG(c2)
      c2r = ABS(c2r)
      c2i = ABS(c2i)
      c2m = MAX(c2r,c2i)
      IF (c2m > ascle) THEN
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1 * c1
        s2 = c2
        s1 = s1 * css(iflag)
        s2 = s2 * css(iflag)
        c1 = csr(iflag)
      END IF
    END IF
  END DO
END IF

40 RETURN
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
50 IF (rs1 <= 0.0_dp) THEN
  y(nd) = czero
  nz = nz + 1
  nd = nd - 1
  IF (nd == 0) GO TO 40
  CALL cuoik(z, fnu, kode, 1, nd, y, nuf, tol, elim, alim)
  IF (nuf >= 0) THEN
    nd = nd - nuf
    nz = nz + nuf
    IF (nd == 0) GO TO 40
    fn = fnu + (nd-1)
    IF (fn >= fnul) GO TO 10
    nlast = nd
    RETURN
  END IF
END IF

60 nz = -1
RETURN

70 IF (rs1 > 0.0_dp) GO TO 60
nz = n
y(1:n) = czero
RETURN
END SUBROUTINE cuni1



SUBROUTINE cuni2(z, fnu, kode, n, y, nz, nlast, fnul, tol, elim, alim)
!***BEGIN PROLOGUE  CUNI2
!***REFER TO  CBESI,CBESK

!     CUNI2 COMPUTES I(FNU,Z) IN THE RIGHT HALF PLANE BY MEANS OF
!     UNIFORM ASYMPTOTIC EXPANSION FOR J(FNU,ZN) WHERE ZN IS Z*I
!     OR -Z*I AND ZN IS IN THE RIGHT HALF PLANE ALSO.

!     FNUL IS THE SMALLEST ORDER PERMITTED FOR THE ASYMPTOTIC
!     EXPANSION.  NLAST=0 MEANS ALL OF THE Y VALUES WERE SET.
!     NLAST.NE.0 IS THE NUMBER LEFT TO BE COMPUTED BY ANOTHER
!     FORMULA FOR ORDERS FNU TO FNU+NLAST-1 BECAUSE FNU+NLAST-1 < FNUL.
!     Y(I) = CZERO FOR I=NLAST+1,N

!***ROUTINES CALLED  CAIRY,CUCHK,CUNHJ,CUOIK,R1MACH
!***END PROLOGUE  CUNI2

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
INTEGER, INTENT(OUT)       :: nlast
REAL (dp), INTENT(IN)      :: fnul
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: ai, arg, asum, bsum, cfn, cid, crsc, cscl, csr(3), css(3), &
                 cy(2), c1, c2, dai, phi, rz, s1, s2, zb, zeta1, zeta2,  &
                 zn, zar
REAL (dp)     :: aarg, ang, aphi, ascle, ay, bry(3), car, c2i, c2m,  &
                 c2r, fn, rs1, sar, yy
INTEGER       :: i, iflag, in, inu, j, k, nai, nd, ndai, nn, nuf, nw, idum
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp), &
                            ci = (0.0_dp,1.0_dp)
COMPLEX (dp), PARAMETER  :: cip(4) = (/ (1.0_dp,0.0_dp), (0.0_dp,1.0_dp),  &
                                        (-1.0_dp,0.0_dp), (0.0_dp,-1.0_dp) /)
REAL (dp), PARAMETER     :: hpi = 1.57079632679489662_dp, aic = 1.265512123484645396_dp

nz = 0
nd = n
nlast = 0
!-----------------------------------------------------------------------
!     COMPUTED VALUES WITH EXPONENTS BETWEEN ALIM AND ELIM IN MAGNITUDE
!     ARE SCALED TO KEEP INTERMEDIATE ARITHMETIC ON SCALE,
!     EXP(ALIM) = EXP(ELIM)*TOL
!-----------------------------------------------------------------------
cscl = CMPLX(1.0_dp/tol, 0.0_dp, KIND=dp)
crsc = tol
css(1) = cscl
css(2) = cone
css(3) = crsc
csr(1) = crsc
csr(2) = cone
csr(3) = cscl
bry(1) = 1.0E+3 * TINY(0.0_dp) / tol
yy = AIMAG(z)
!-----------------------------------------------------------------------
!     ZN IS IN THE RIGHT HALF PLANE AFTER ROTATION BY CI OR -CI
!-----------------------------------------------------------------------
zn = -z * ci
zb = z
cid = -ci
inu = fnu
ang = hpi * (fnu - inu)
car = COS(ang)
sar = SIN(ang)
c2 = CMPLX(car, sar, KIND=dp)
zar = c2
in = inu + n - 1
in = MOD(in,4)
c2 = c2 * cip(in+1)
IF (yy <= 0.0_dp) THEN
  zn = CONJG(-zn)
  zb = CONJG(zb)
  cid = -cid
  c2 = CONJG(c2)
END IF
!-----------------------------------------------------------------------
!     CHECK FOR UNDERFLOW AND OVERFLOW ON FIRST MEMBER
!-----------------------------------------------------------------------
fn = MAX(fnu,1.0_dp)
CALL cunhj(zn, fn, 1, tol, phi, arg, zeta1, zeta2, asum, bsum)
IF (kode /= 1) THEN
  cfn = fnu
  s1 = -zeta1 + cfn * (cfn/(zb + zeta2))
ELSE
  s1 = -zeta1 + zeta2
END IF
rs1 = REAL(s1, KIND=dp)
IF (ABS(rs1) > elim) GO TO 70

10 nn = MIN(2,nd)
DO  i = 1, nn
  fn = fnu + (nd-i)
  CALL cunhj(zn, fn, 0, tol, phi, arg, zeta1, zeta2, asum, bsum)
  IF (kode /= 1) THEN
    cfn = fn
    ay = ABS(yy)
    s1 = -zeta1 + cfn * (cfn/(zb+zeta2)) + CMPLX(0.0_dp, ay, KIND=dp)
  ELSE
    s1 = -zeta1 + zeta2
  END IF
!-----------------------------------------------------------------------
!     TEST FOR UNDERFLOW AND OVERFLOW
!-----------------------------------------------------------------------
  rs1 = REAL(s1, KIND=dp)
  IF (ABS(rs1) > elim) GO TO 50
  IF (i == 1) iflag = 2
  IF (ABS(rs1) >= alim) THEN
!-----------------------------------------------------------------------
!     REFINE  TEST AND SCALE
!-----------------------------------------------------------------------
    aphi = ABS(phi)
    aarg = ABS(arg)
    rs1 = rs1 + LOG(aphi) - 0.25_dp * LOG(aarg) - aic
    IF (ABS(rs1) > elim) GO TO 50
    IF (i == 1) iflag = 1
    IF (rs1 >= 0.0_dp) THEN
      IF (i == 1) iflag = 3
    END IF
  END IF
!-----------------------------------------------------------------------
!     SCALE S1 TO KEEP INTERMEDIATE ARITHMETIC ON SCALE NEAR
!     EXPONENT EXTREMES
!-----------------------------------------------------------------------
  CALL cairy(arg, 0, 2, ai, nai, idum)
  CALL cairy(arg, 1, 2, dai, ndai, idum)
  s2 = phi * (ai*asum + dai*bsum)
  c2r = REAL(s1, KIND=dp)
  c2i = AIMAG(s1)
  c2m = EXP(c2r) * REAL(css(iflag), KIND=dp)
  s1 = c2m * CMPLX(COS(c2i), SIN(c2i), KIND=dp)
  s2 = s2 * s1
  IF (iflag == 1) THEN
    CALL cuchk(s2, nw, bry(1), tol)
    IF (nw /= 0) GO TO 50
  END IF
  IF (yy <= 0.0_dp) s2 = CONJG(s2)
  j = nd - i + 1
  s2 = s2 * c2
  cy(i) = s2
  y(j) = s2 * csr(iflag)
  c2 = c2 * cid
END DO
IF (nd > 2) THEN
  rz = 2.0_dp / z
  bry(2) = 1.0_dp / bry(1)
  bry(3) = HUGE(0.0_dp)
  s1 = cy(1)
  s2 = cy(2)
  c1 = csr(iflag)
  ascle = bry(iflag)
  k = nd - 2
  fn = k
  DO  i = 3, nd
    c2 = s2
    s2 = s1 + (fnu + fn) * rz * s2
    s1 = c2
    c2 = s2 * c1
    y(k) = c2
    k = k - 1
    fn = fn - 1.0_dp
    IF (iflag < 3) THEN
      c2r = REAL(c2, KIND=dp)
      c2i = AIMAG(c2)
      c2r = ABS(c2r)
      c2i = ABS(c2i)
      c2m = MAX(c2r,c2i)
      IF (c2m > ascle) THEN
        iflag = iflag + 1
        ascle = bry(iflag)
        s1 = s1 * c1
        s2 = c2
        s1 = s1 * css(iflag)
        s2 = s2 * css(iflag)
        c1 = csr(iflag)
      END IF
    END IF
  END DO
END IF

40 RETURN

50 IF (rs1 <= 0.0_dp) THEN
!-----------------------------------------------------------------------
!     SET UNDERFLOW AND UPDATE PARAMETERS
!-----------------------------------------------------------------------
  y(nd) = czero
  nz = nz + 1
  nd = nd - 1
  IF (nd == 0) GO TO 40
  CALL cuoik(z, fnu, kode, 1, nd, y, nuf, tol, elim, alim)
  IF (nuf >= 0) THEN
    nd = nd - nuf
    nz = nz + nuf
    IF (nd == 0) GO TO 40
    fn = fnu + (nd-1)
    IF (fn >= fnul) THEN
!      FN = AIMAG(CID)
!      J = NUF + 1
!      K = MOD(J,4) + 1
!      S1 = CIP(K)
!      IF (FN < 0.0_dp) S1 = CONJG(S1)
!      C2 = C2*S1
      in = inu + nd - 1
      in = MOD(in,4) + 1
      c2 = zar * cip(in)
      IF (yy <= 0.0_dp) c2 = CONJG(c2)
      GO TO 10
    END IF
    nlast = nd
    RETURN
  END IF
END IF

60 nz = -1
RETURN

70 IF (rs1 > 0.0_dp) GO TO 60
nz = n
y(1:n) = czero
RETURN
END SUBROUTINE cuni2



SUBROUTINE cs1s2(zr, s1, s2, nz, ascle, alim, iuf)
!***BEGIN PROLOGUE  CS1S2
!***REFER TO  CBESK,CAIRY

!     CS1S2 TESTS FOR A POSSIBLE UNDERFLOW RESULTING FROM THE ADDITION OF
!     THE I AND K FUNCTIONS IN THE ANALYTIC CONTINUATION FORMULA WHERE S1=K
!     FUNCTION AND S2=I FUNCTION.
!     ON KODE=1 THE I AND K FUNCTIONS ARE DIFFERENT ORDERS OF MAGNITUDE,
!     BUT FOR KODE=2 THEY CAN BE OF THE SAME ORDER OF MAGNITUDE AND THE
!     MAXIMUM MUST BE AT LEAST ONE PRECISION ABOVE THE UNDERFLOW LIMIT.

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  CS1S2

COMPLEX (dp), INTENT(IN)      :: zr
COMPLEX (dp), INTENT(IN OUT)  :: s1
COMPLEX (dp), INTENT(IN OUT)  :: s2
INTEGER, INTENT(OUT)          :: nz
REAL (dp), INTENT(IN)         :: ascle
REAL (dp), INTENT(IN)         :: alim
INTEGER, INTENT(IN OUT)       :: iuf

COMPLEX (dp)  :: c1, s1d
REAL (dp)     :: aa, aln, as1, as2, xx

COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp)

nz = 0
as1 = ABS(s1)
as2 = ABS(s2)
aa = REAL(s1, KIND=dp)
aln = AIMAG(s1)
IF (aa /= 0.0_dp .OR. aln /= 0.0_dp) THEN
  IF (as1 /= 0.0_dp) THEN
    xx = REAL(zr, KIND=dp)
    aln = -xx - xx + LOG(as1)
    s1d = s1
    s1 = czero
    as1 = 0.0_dp
    IF (aln >= -alim) THEN
      c1 = LOG(s1d) - zr - zr
      s1 = EXP(c1)
      as1 = ABS(s1)
      iuf = iuf + 1
    END IF
  END IF
END IF
aa = MAX(as1,as2)
IF (aa > ascle) RETURN
s1 = czero
s2 = czero
nz = 1
iuf = 0
RETURN
END SUBROUTINE cs1s2



SUBROUTINE cshch(z, csh, cch)
!***BEGIN PROLOGUE  CSHCH
!***REFER TO  CBESK,CBESH

!     CSHCH COMPUTES THE COMPLEX HYPERBOLIC FUNCTIONS CSH=SINH(X+I*Y)
!     AND CCH=COSH(X+I*Y), WHERE I**2=-1.

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  CSHCH

COMPLEX (dp), INTENT(IN OUT)  :: z
COMPLEX (dp), INTENT(OUT)     :: csh
COMPLEX (dp), INTENT(OUT)     :: cch

REAL (dp)  :: cchi, cchr, ch, cn, cshi, cshr, sh, sn, x, y

x = REAL(z, KIND=dp)
y = AIMAG(z)
sh = SINH(x)
ch = COSH(x)
sn = SIN(y)
cn = COS(y)
cshr = sh * cn
cshi = ch * sn
csh = CMPLX(cshr, cshi, KIND=dp)
cchr = ch * cn
cchi = sh * sn
cch = CMPLX(cchr, cchi, KIND=dp)
RETURN
END SUBROUTINE cshch



SUBROUTINE crati(z, fnu, n, cy, tol)
!***BEGIN PROLOGUE  CRATI
!***REFER TO  CBESI,CBESK,CBESH

!   CRATI COMPUTES RATIOS OF I BESSEL FUNCTIONS BY BACKWARD RECURRENCE.
!   THE STARTING INDEX IS DETERMINED BY FORWARD RECURRENCE AS DESCRIBED IN
!   J. RES. OF NAT. BUR. OF STANDARDS-B, MATHEMATICAL SCIENCES, VOL 77B,
!   P111-114, SEPTEMBER 1973, BESSEL FUNCTIONS I AND J OF COMPLEX ARGUMENT
!   AND INTEGER ORDER, BY D. J. SOOKNE.

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  CRATI

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: cy(n)
REAL (dp), INTENT(IN)      :: tol

COMPLEX (dp)  :: cdfnu, pt, p1, p2, rz, t1
REAL (dp)     :: ak, amagz, ap1, ap2, arg, az, dfnu, fdnu, flam, fnup,  &
                 rap1, rho, test, test1
INTEGER       :: i, id, idnu, inu, itime, k, kk, magz

REAL (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp)

az = ABS(z)
inu = fnu
idnu = inu + n - 1
fdnu = idnu
magz = az
amagz = magz + 1
fnup = MAX(amagz, fdnu)
id = idnu - magz - 1
itime = 1
k = 1
rz = (cone+cone) / z
t1 = fnup * rz
p2 = -t1
p1 = cone
t1 = t1 + rz
IF (id > 0) id = 0
ap2 = ABS(p2)
ap1 = ABS(p1)
!-----------------------------------------------------------------------
!     THE OVERFLOW TEST ON K(FNU+I-1,Z) BEFORE THE CALL TO CBKNX
!     GUARANTEES THAT P2 IS ON SCALE. SCALE TEST1 AND ALL SUBSEQUENT P2
!     VALUES BY AP1 TO ENSURE THAT AN OVERFLOW DOES NOT OCCUR PREMATURELY.
!-----------------------------------------------------------------------
arg = (ap2+ap2) / (ap1*tol)
test1 = SQRT(arg)
test = test1
rap1 = 1.0_dp / ap1
p1 = p1 * rap1
p2 = p2 * rap1
ap2 = ap2 * rap1

10 k = k + 1
ap1 = ap2
pt = p2
p2 = p1 - t1 * p2
p1 = pt
t1 = t1 + rz
ap2 = ABS(p2)
IF (ap1 <= test) GO TO 10
IF (itime /= 2) THEN
  ak = ABS(t1) * 0.5_dp
  flam = ak + SQRT(ak*ak - 1.0_dp)
  rho = MIN(ap2/ap1, flam)
  test = test1 * SQRT(rho/(rho*rho - 1.0_dp))
  itime = 2
  GO TO 10
END IF
kk = k + 1 - id
ak = kk
dfnu = fnu + (n-1)
cdfnu = dfnu
t1 = ak
p1 = 1.0_dp/ap2
p2 = czero
DO  i = 1, kk
  pt = p1
  p1 = rz * (cdfnu+t1) * p1 + p2
  p2 = pt
  t1 = t1 - cone
END DO
IF (REAL(p1, KIND=dp) == 0.0_dp .AND. AIMAG(p1) == 0.0_dp) THEN
  p1 = CMPLX(tol, tol, KIND=dp)
END IF
cy(n) = p2 / p1
IF (n == 1) RETURN
k = n - 1
ak = k
t1 = ak
cdfnu = fnu * rz
DO  i = 2, n
  pt = cdfnu + t1 * rz + cy(k+1)
  IF (REAL(pt, KIND=dp) == 0.0_dp .AND. AIMAG(pt) == 0.0_dp) THEN
    pt = CMPLX(tol, tol, KIND=dp)
  END IF
  cy(k) = cone / pt
  t1 = t1 - cone
  k = k - 1
END DO
RETURN
END SUBROUTINE crati



SUBROUTINE cbknu(z, fnu, kode, n, y, nz, tol, elim, alim)
!***BEGIN PROLOGUE  CBKNU
!***REFER TO  CBESI,CBESK,CAIRY,CBESH

!     CBKNU COMPUTES THE K BESSEL FUNCTION IN THE RIGHT HALF Z PLANE

!***ROUTINES CALLED  CKSCL,CSHCH,GAMLN,I1MACH,R1MACH,CUCHK
!***END PROLOGUE  CBKNU

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: cch, ck, coef, crsc, cs, cscl, csh, csr(3), css(3), cz,  &
                 f, fmu, p, pt, p1, p2, q, rz, smu, st, s1, s2, zd, celm, &
                 cy(2)
REAL (dp)     :: aa, ak, ascle, a1, a2, bb, bk, bry(3), caz, dnu, dnu2,  &
                 etest, fc, fhs, fk, fks, g1, g2, p2i, p2m, p2r, rk, s,  &
                 tm, t1, t2, xx, yy, helim, elm, xd, yd, alas, as
INTEGER       :: i, iflag, inu, k, kflag, kk, koded, nw, j, ic, inub

INTEGER, PARAMETER       :: kmax = 30
REAL (dp), PARAMETER     :: r1 = 2.0_dp
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp), cone = (1.0_dp,0.0_dp), &
                            ctwo = (2.0_dp,0.0_dp)

REAL (dp), PARAMETER  :: pi = 3.14159265358979324_dp,  &
    rthpi = 1.25331413731550025_dp, spi = 1.90985931710274403_dp,  &
    hpi = 1.57079632679489662_dp, fpi = 1.89769999331517738_dp,  &
    tth = 6.66666666666666666D-01

REAL (dp), PARAMETER  :: cc(8) = (/ 5.77215664901532861D-01,  &
 -4.20026350340952355D-02, -4.21977345555443367D-02, 7.21894324666309954D-03, &
 -2.15241674114950973D-04, -2.01348547807882387D-05, 1.13302723198169588D-06, &
  6.11609510448141582D-09 /)

xx = REAL(z, KIND=dp)
yy = AIMAG(z)
caz = ABS(z)
cscl = 1.0_dp/tol
crsc = tol
css(1) = cscl
css(2) = cone
css(3) = crsc
csr(1) = crsc
csr(2) = cone
csr(3) = cscl
bry(1) = 1.0E+3 * TINY(0.0_dp) / tol
bry(2) = 1.0_dp / bry(1)
bry(3) = HUGE(0.0_dp)
nz = 0
iflag = 0
koded = kode
rz = ctwo / z
inu = fnu + 0.5_dp
dnu = fnu - inu
IF (ABS(dnu) /= 0.5_dp) THEN
  dnu2 = 0.0_dp
  IF (ABS(dnu) > tol) dnu2 = dnu * dnu
  IF (caz <= r1) THEN
!-----------------------------------------------------------------------
!     SERIES FOR ABS(Z) <= R1
!-----------------------------------------------------------------------
    fc = 1.0_dp
    smu = LOG(rz)
    fmu = smu * dnu
    CALL cshch(fmu, csh, cch)
    IF (dnu /= 0.0_dp) THEN
      fc = dnu * pi
      fc = fc / SIN(fc)
      smu = csh / dnu
    END IF
    a2 = 1.0_dp + dnu
!-----------------------------------------------------------------------
!     GAM(1-Z)*GAM(1+Z)=PI*Z/SIN(PI*Z), T1=1/GAM(1-DNU), T2=1/GAM(1+DNU)
!-----------------------------------------------------------------------
    t2 = EXP(-gamln(a2))
    t1 = 1.0_dp / (t2*fc)
    IF (ABS(dnu) <= 0.1_dp) THEN
!-----------------------------------------------------------------------
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
!-----------------------------------------------------------------------
      ak = 1.0_dp
      s = cc(1)
      DO  k = 2, 8
        ak = ak * dnu2
        tm = cc(k) * ak
        s = s + tm
        IF (ABS(tm) < tol) EXIT
      END DO
      g1 = -s
    ELSE
      g1 = (t1-t2) / (dnu+dnu)
    END IF
    g2 = 0.5_dp * (t1+t2) * fc
    g1 = g1 * fc
    f = g1 * cch + smu * g2
    pt = EXP(fmu)
    p = CMPLX(0.5_dp/t2, 0.0_dp, KIND=dp) * pt
    q = CMPLX(0.5_dp/t1, 0.0_dp, KIND=dp) / pt
    s1 = f
    s2 = p
    ak = 1.0_dp
    a1 = 1.0_dp
    ck = cone
    bk = 1.0_dp - dnu2
    IF (inu <= 0 .AND. n <= 1) THEN
!-----------------------------------------------------------------------
!     GENERATE K(FNU,Z), 0.0D0  <=  FNU  <  0.5D0 AND N=1
!-----------------------------------------------------------------------
      IF (caz >= tol) THEN
        cz = z * z * 0.25_dp
        t1 = 0.25_dp * caz * caz

        30 f = (f*ak + p + q) / bk
        p = p / (ak-dnu)
        q = q / (ak+dnu)
        rk = 1.0_dp / ak
        ck = ck * cz * rk
        s1 = s1 + ck * f
        a1 = a1 * t1 * rk
        bk = bk + ak + ak + 1.0_dp
        ak = ak + 1.0_dp
        IF (a1 > tol) GO TO 30
      END IF
      y(1) = s1
      IF (koded == 1) RETURN
      y(1) = s1 * EXP(z)
      RETURN
    END IF
!-----------------------------------------------------------------------
!     GENERATE K(DNU,Z) AND K(DNU+1,Z) FOR FORWARD RECURRENCE
!-----------------------------------------------------------------------
    IF (caz >= tol) THEN
      cz = z * z * 0.25_dp
      t1 = 0.25_dp * caz * caz

      40 f = (f*ak + p + q) / bk
      p = p / (ak-dnu)
      q = q / (ak+dnu)
      rk = 1.0_dp / ak
      ck = ck * cz * rk
      s1 = s1 + ck * f
      s2 = s2 + ck * (p - f*ak)
      a1 = a1 * t1 * rk
      bk = bk + ak + ak + 1.0_dp
      ak = ak + 1.0_dp
      IF (a1 > tol) GO TO 40
    END IF
    kflag = 2
    bk = REAL(smu, KIND=dp)
    a1 = fnu + 1.0_dp
    ak = a1 * ABS(bk)
    IF (ak > alim) kflag = 3
    p2 = s2 * css(kflag)
    s2 = p2 * rz
    s1 = s1 * css(kflag)
    IF (koded == 1) GO TO 100
    f = EXP(z)
    s1 = s1 * f
    s2 = s2 * f
    GO TO 100
  END IF
END IF
!-----------------------------------------------------------------------
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
!-----------------------------------------------------------------------
coef = rthpi / SQRT(z)
kflag = 2
IF (koded /= 2) THEN
  IF (xx > alim) GO TO 200
!     BLANK LINE
  a1 = EXP(-xx) * REAL(css(kflag), KIND=dp)
  pt = a1 * CMPLX(COS(yy), -SIN(yy), KIND=dp)
  coef = coef * pt
END IF

50 IF (ABS(dnu) == 0.5_dp) GO TO 210
!-----------------------------------------------------------------------
!     MILLER ALGORITHM FOR ABS(Z) > R1
!-----------------------------------------------------------------------
ak = COS(pi*dnu)
ak = ABS(ak)
IF (ak == 0.0_dp) GO TO 210
fhs = ABS(0.25_dp - dnu2)
IF (fhs == 0.0_dp) GO TO 210
!-----------------------------------------------------------------------
!     COMPUTE R2=F(E). IF ABS(Z) >= R2, USE FORWARD RECURRENCE TO
!     DETERMINE THE BACKWARD INDEX K. R2=F(E) IS A STRAIGHT LINE ON
!     12 <= E <= 60. E IS COMPUTED FROM 2**(-E)=B**(1-DIGITS(0.0_dp))=
!     TOL WHERE B IS THE BASE OF THE ARITHMETIC.
!-----------------------------------------------------------------------
t1 = (DIGITS(0.0_dp) - 1) * LOG10( REAL( RADIX(0.0_dp), KIND=dp) ) * 3.321928094_dp
t1 = MAX(t1,12.0_dp)
t1 = MIN(t1,60.0_dp)
t2 = tth * t1 - 6.0_dp
IF (xx == 0.0_dp) THEN
  t1 = hpi
ELSE
  t1 = ATAN(yy/xx)
  t1 = ABS(t1)
END IF
IF (t2 <= caz) THEN
!-----------------------------------------------------------------------
!     FORWARD RECURRENCE LOOP WHEN ABS(Z) >= R2
!-----------------------------------------------------------------------
  etest = ak / (pi*caz*tol)
  fk = 1.0_dp
  IF (etest < 1.0_dp) GO TO 80
  fks = 2.0_dp
  rk = caz + caz + 2.0_dp
  a1 = 0.0_dp
  a2 = 1.0_dp
  DO  i = 1, kmax
    ak = fhs / fks
    bk = rk / (fk+1.0_dp)
    tm = a2
    a2 = bk * a2 - ak * a1
    a1 = tm
    rk = rk + 2.0_dp
    fks = fks + fk + fk + 2.0_dp
    fhs = fhs + fk + fk
    fk = fk + 1.0_dp
    tm = ABS(a2) * fk
    IF (etest < tm) GO TO 70
  END DO
  GO TO 220

  70 fk = fk + spi * t1 * SQRT(t2/caz)
  fhs = ABS(0.25_dp-dnu2)
ELSE
!-----------------------------------------------------------------------
!     COMPUTE BACKWARD INDEX K FOR ABS(Z) < R2
!-----------------------------------------------------------------------
  a2 = SQRT(caz)
  ak = fpi * ak / (tol*SQRT(a2))
  aa = 3.0_dp * t1 / (1.0_dp+caz)
  bb = 14.7_dp * t1 / (28.0_dp+caz)
  ak = (LOG(ak) + caz*COS(aa)/(1.0_dp + 0.008_dp*caz)) / COS(bb)
  fk = 0.12125_dp * ak * ak / caz + 1.5_dp
END IF

80 k = fk
!-----------------------------------------------------------------------
!     BACKWARD RECURRENCE LOOP FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
fk = k
fks = fk * fk
p1 = czero
p2 = tol
cs = p2
DO  i = 1, k
  a1 = fks - fk
  a2 = (fks+fk) / (a1+fhs)
  rk = 2.0_dp / (fk + 1.0_dp)
  t1 = (fk+xx) * rk
  t2 = yy * rk
  pt = p2
  p2 = (p2*CMPLX(t1, t2, KIND=dp) - p1) * a2
  p1 = pt
  cs = cs + p2
  fks = a1 - fk + 1.0_dp
  fk = fk - 1.0_dp
END DO
!-----------------------------------------------------------------------
!     COMPUTE (P2/CS)=(P2/ABS(CS))*(CONJG(CS)/ABS(CS)) FOR BETTER SCALING
!-----------------------------------------------------------------------
tm = ABS(cs)
pt = CMPLX(1.0_dp/tm, 0.0_dp, KIND=dp)
s1 = pt * p2
cs = CONJG(cs) * pt
s1 = coef * s1 * cs
IF (inu <= 0 .AND. n <= 1) THEN
  zd = z
  IF (iflag == 1) GO TO 190
  GO TO 130
END IF
!-----------------------------------------------------------------------
!     COMPUTE P1/P2=(P1/ABS(P2)*CONJG(P2)/ABS(P2) FOR SCALING
!-----------------------------------------------------------------------
tm = ABS(p2)
pt = CMPLX(1.0_dp/tm, 0.0_dp, KIND=dp)
p1 = pt * p1
p2 = CONJG(p2) * pt
pt = p1 * p2
s2 = s1 * (cone + (CMPLX(dnu+0.5_dp, 0.0_dp, KIND=dp) - pt)/z)
!-----------------------------------------------------------------------
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION WITH
!     SCALING NEAR EXPONENT EXTREMES ON KFLAG=1 OR KFLAG=3
!-----------------------------------------------------------------------
100 ck = CMPLX(dnu+1.0_dp, 0.0_dp, KIND=dp) * rz
IF (n == 1) inu = inu - 1
IF (inu <= 0) THEN
  IF (n == 1) s1 = s2
  zd = z
  IF (iflag == 1) GO TO 190
  GO TO 130
END IF
inub = 1
IF (iflag == 1) GO TO 160

110 p1 = csr(kflag)
ascle = bry(kflag)
DO  i = inub, inu
  st = s2
  s2 = ck * s2 + s1
  s1 = st
  ck = ck + rz
  IF (kflag < 3) THEN
    p2 = s2 * p1
    p2r = REAL(p2, KIND=dp)
    p2i = AIMAG(p2)
    p2r = ABS(p2r)
    p2i = ABS(p2i)
    p2m = MAX(p2r,p2i)
    IF (p2m > ascle) THEN
      kflag = kflag + 1
      ascle = bry(kflag)
      s1 = s1 * p1
      s2 = p2
      s1 = s1 * css(kflag)
      s2 = s2 * css(kflag)
      p1 = csr(kflag)
    END IF
  END IF
END DO
IF (n == 1) s1 = s2

130 y(1) = s1 * csr(kflag)
IF (n == 1) RETURN
y(2) = s2 * csr(kflag)
IF (n == 2) RETURN
kk = 2

140 kk = kk + 1
IF (kk > n) RETURN
p1 = csr(kflag)
ascle = bry(kflag)
DO  i = kk, n
  p2 = s2
  s2 = ck * s2 + s1
  s1 = p2
  ck = ck + rz
  p2 = s2 * p1
  y(i) = p2
  IF (kflag < 3) THEN
    p2r = REAL(p2, KIND=dp)
    p2i = AIMAG(p2)
    p2r = ABS(p2r)
    p2i = ABS(p2i)
    p2m = MAX(p2r,p2i)
    IF (p2m > ascle) THEN
      kflag = kflag + 1
      ascle = bry(kflag)
      s1 = s1 * p1
      s2 = p2
      s1 = s1 * css(kflag)
      s2 = s2 * css(kflag)
      p1 = csr(kflag)
    END IF
  END IF
END DO
RETURN
!-----------------------------------------------------------------------
!     IFLAG=1 CASES, FORWARD RECURRENCE ON SCALED VALUES ON UNDERFLOW
!-----------------------------------------------------------------------
160 helim = 0.5_dp * elim
elm = EXP(-elim)
celm = elm
ascle = bry(1)
zd = z
xd = xx
yd = yy
ic = -1
j = 2
DO  i = 1, inu
  st = s2
  s2 = ck * s2 + s1
  s1 = st
  ck = ck + rz
  as = ABS(s2)
  alas = LOG(as)
  p2r = -xd + alas
  IF (p2r >= -elim) THEN
    p2 = -zd + LOG(s2)
    p2r = REAL(p2, KIND=dp)
    p2i = AIMAG(p2)
    p2m = EXP(p2r) / tol
    p1 = p2m * CMPLX(COS(p2i), SIN(p2i), KIND=dp)
    CALL cuchk(p1, nw, ascle, tol)
    IF (nw == 0) THEN
      j = 3 - j
      cy(j) = p1
      IF (ic == i-1) GO TO 180
      ic = i
      CYCLE
    END IF
  END IF
  IF (alas >= helim) THEN
    xd = xd - elim
    s1 = s1 * celm
    s2 = s2 * celm
    zd = CMPLX(xd, yd, KIND=dp)
  END IF
END DO
IF (n == 1) s1 = s2
GO TO 190

180 kflag = 1
inub = i + 1
s2 = cy(j)
j = 3 - j
s1 = cy(j)
IF (inub <= inu) GO TO 110
IF (n == 1) s1 = s2
GO TO 130

190 y(1) = s1
IF (n /= 1) THEN
  y(2) = s2
END IF
ascle = bry(1)
CALL ckscl(zd, fnu, n, y, nz, rz, ascle, tol, elim)
inu = n - nz
IF (inu <= 0) RETURN
kk = nz + 1
s1 = y(kk)
y(kk) = s1 * csr(1)
IF (inu == 1) RETURN
kk = nz + 2
s2 = y(kk)
y(kk) = s2 * csr(1)
IF (inu == 2) RETURN
t2 = fnu + (kk-1)
ck = t2 * rz
kflag = 1
GO TO 140
!-----------------------------------------------------------------------
!     SCALE BY EXP(Z), IFLAG = 1 CASES
!-----------------------------------------------------------------------
200 koded = 2
iflag = 1
kflag = 2
GO TO 50
!-----------------------------------------------------------------------
!     FNU=HALF ODD INTEGER CASE, DNU=-0.5
!-----------------------------------------------------------------------
210 s1 = coef
s2 = coef
GO TO 100

220 nz = -2
RETURN
END SUBROUTINE cbknu



SUBROUTINE ckscl(zr, fnu, n, y, nz, rz, ascle, tol, elim)
!***BEGIN PROLOGUE  CKSCL
!***REFER TO  CBKNU,CUNK1,CUNK2

!     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
!     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
!     RETURN WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.

!***ROUTINES CALLED  CUCHK
!***END PROLOGUE  CKSCL

COMPLEX (dp), INTENT(IN)   :: zr
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
COMPLEX (dp), INTENT(IN)   :: rz
REAL (dp), INTENT(IN OUT)  :: ascle
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim

COMPLEX (dp)  :: ck, cs, cy(2), s1, s2, zd, celm
REAL (dp)     :: aa, acs, as, csi, csr, fn, xx, zri, elm, alas, helim
INTEGER       :: i, ic, kk, nn, nw
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp,0.0_dp)

nz = 0
ic = 0
xx = REAL(zr, KIND=dp)
nn = MIN(2,n)
DO  i = 1, nn
  s1 = y(i)
  cy(i) = s1
  as = ABS(s1)
  acs = -xx + LOG(as)
  nz = nz + 1
  y(i) = czero
  IF (acs >= -elim) THEN
    cs = -zr + LOG(s1)
    csr = REAL(cs, KIND=dp)
    csi = AIMAG(cs)
    aa = EXP(csr) / tol
    cs = aa * CMPLX(COS(csi), SIN(csi), KIND=dp)
    CALL cuchk(cs, nw, ascle, tol)
    IF (nw == 0) THEN
      y(i) = cs
      nz = nz - 1
      ic = i
    END IF
  END IF
END DO
IF (n == 1) RETURN
IF (ic <= 1) THEN
  y(1) = czero
  nz = 2
END IF
IF (n == 2) RETURN
IF (nz == 0) RETURN
fn = fnu + 1.0_dp
ck = fn * rz
s1 = cy(1)
s2 = cy(2)
helim = 0.5_dp * elim
elm = EXP(-elim)
celm = elm
zri = AIMAG(zr)
zd = zr

!     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
!     S2 GETS LARGER THAN EXP(ELIM/2)

DO  i = 3, n
  kk = i
  cs = s2
  s2 = ck * s2 + s1
  s1 = cs
  ck = ck + rz
  as = ABS(s2)
  alas = LOG(as)
  acs = -xx + alas
  nz = nz + 1
  y(i) = czero
  IF (acs >= -elim) THEN
    cs = -zd + LOG(s2)
    csr = REAL(cs, KIND=dp)
    csi = AIMAG(cs)
    aa = EXP(csr) / tol
    cs = aa * CMPLX(COS(csi), SIN(csi), KIND=dp)
    CALL cuchk(cs, nw, ascle, tol)
    IF (nw == 0) THEN
      y(i) = cs
      nz = nz - 1
      IF (ic == kk-1) GO TO 30
      ic = kk
      CYCLE
    END IF
  END IF
  IF (alas >= helim) THEN
    xx = xx - elim
    s1 = s1 * celm
    s2 = s2 * celm
    zd = CMPLX(xx, zri, KIND=dp)
  END IF
END DO
nz = n
IF (ic == n) nz = n - 1
GO TO 40

30 nz = kk - 2

40 y(1:nz) = czero
RETURN
END SUBROUTINE ckscl



SUBROUTINE cacon(z, fnu, kode, mr, n, y, nz, rl, fnul, tol, elim, alim)
!***BEGIN PROLOGUE  CACON
!***REFER TO  CBESK,CBESH

!     CACON APPLIES THE ANALYTIC CONTINUATION FORMULA

!         K(FNU,ZN*EXP(MP))=K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!                 MP=PI*MR*CMPLX(0.0,1.0)

!     TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT
!     HALF Z PLANE

!***ROUTINES CALLED  CBINU,CBKNU,CS1S2,R1MACH
!***END PROLOGUE  CACON

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: mr
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN OUT)  :: rl
REAL (dp), INTENT(IN OUT)  :: fnul
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN OUT)  :: elim
REAL (dp), INTENT(IN OUT)  :: alim

COMPLEX (dp)  :: ck, cs, cscl, cscr, csgn, cspn, css(3), csr(3), c1, c2,  &
                 rz, sc1, sc2, st, s1, s2, zn, cy(2)
REAL (dp)     :: arg, ascle, as2, bscle, bry(3), cpn, c1i, c1m, c1r,  &
                 fmr, sgn, spn, yy
INTEGER       :: i, inu, iuf, kflag, nn, nw
REAL (dp), PARAMETER     :: pi = 3.14159265358979324_dp
COMPLEX (dp), PARAMETER  :: cone = (1.0_dp, 0.0_dp)

nz = 0
zn = -z
nn = n
CALL cbinu(zn, fnu, kode, nn, y, nw, rl, fnul, tol, elim, alim)
IF (nw >= 0) THEN
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
  nn = MIN(2, n)
  CALL cbknu(zn, fnu, kode, nn, cy, nw, tol, elim, alim)
  IF (nw == 0) THEN
    s1 = cy(1)
    fmr = mr
    sgn = -SIGN(pi, fmr)
    csgn = CMPLX(0.0_dp, sgn, KIND=dp)
    IF (kode /= 1) THEN
      yy = -AIMAG(zn)
      cpn = COS(yy)
      spn = SIN(yy)
      csgn = csgn * CMPLX(cpn, spn, KIND=dp)
    END IF
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
    inu = fnu
    arg = (fnu - inu) * sgn
    cpn = COS(arg)
    spn = SIN(arg)
    cspn = CMPLX(cpn, spn, KIND=dp)
    IF (MOD(inu, 2) == 1) cspn = -cspn
    iuf = 0
    c1 = s1
    c2 = y(1)
    ascle = 1.0E+3 * TINY(0.0_dp) / tol
    IF (kode /= 1) THEN
      CALL cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      sc1 = c1
    END IF
    y(1) = cspn * c1 + csgn * c2
    IF (n == 1) RETURN
    cspn = -cspn
    s2 = cy(2)
    c1 = s2
    c2 = y(2)
    IF (kode /= 1) THEN
      CALL cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
      nz = nz + nw
      sc2 = c1
    END IF
    y(2) = cspn * c1 + csgn * c2
    IF (n == 2) RETURN
    cspn = -cspn
    rz = 2.0_dp / zn
    ck = CMPLX(fnu+1.0_dp, 0.0_dp, KIND=dp) * rz
!-----------------------------------------------------------------------
!     SCALE NEAR EXPONENT EXTREMES DURING RECURRENCE ON K FUNCTIONS
!-----------------------------------------------------------------------
    cscl = 1.0_dp/tol
    cscr = tol
    css(1) = cscl
    css(2) = cone
    css(3) = cscr
    csr(1) = cscr
    csr(2) = cone
    csr(3) = cscl
    bry(1) = ascle
    bry(2) = 1.0_dp / ascle
    bry(3) = HUGE(0.0_dp)
    as2 = ABS(s2)
    kflag = 2
    IF (as2 <= bry(1)) THEN
      kflag = 1
    ELSE
      IF (as2 >= bry(2)) THEN
        kflag = 3
      END IF
    END IF
    bscle = bry(kflag)
    s1 = s1 * css(kflag)
    s2 = s2 * css(kflag)
    cs = csr(kflag)
    DO  i = 3, n
      st = s2
      s2 = ck * s2 + s1
      s1 = st
      c1 = s2 * cs
      st = c1
      c2 = y(i)
      IF (kode /= 1) THEN
        IF (iuf >= 0) THEN
          CALL cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
          nz = nz + nw
          sc1 = sc2
          sc2 = c1
          IF (iuf == 3) THEN
            iuf = -4
            s1 = sc1 * css(kflag)
            s2 = sc2 * css(kflag)
            st = sc2
          END IF
        END IF
      END IF
      y(i) = cspn * c1 + csgn * c2
      ck = ck + rz
      cspn = -cspn
      IF (kflag < 3) THEN
        c1r = REAL(c1, KIND=dp)
        c1i = AIMAG(c1)
        c1r = ABS(c1r)
        c1i = ABS(c1i)
        c1m = MAX(c1r, c1i)
        IF (c1m > bscle) THEN
          kflag = kflag + 1
          bscle = bry(kflag)
          s1 = s1 * cs
          s2 = st
          s1 = s1 * css(kflag)
          s2 = s2 * css(kflag)
          cs = csr(kflag)
        END IF
      END IF
    END DO
    RETURN
  END IF
END IF
nz = -1
IF (nw == -2) nz = -2
RETURN
END SUBROUTINE cacon



SUBROUTINE cbinu(z, fnu, kode, n, cy, nz, rl, fnul, tol, elim, alim)
!***BEGIN PROLOGUE  CBINU
!***REFER TO  CBESH,CBESI,CBESJ,CBESK,CAIRY,CBIRY

!     CBINU COMPUTES THE I FUNCTION IN THE RIGHT HALF Z PLANE

!***ROUTINES CALLED  CASYI,CBUNI,CMLRI,CSERI,CUOIK,CWRSK
!***END PROLOGUE  CBINU

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: cy(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: rl
REAL (dp), INTENT(IN)      :: fnul
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: cw(2)
REAL (dp)     :: az, dfnu
INTEGER       :: inw, nlast, nn, nui, nw
COMPLEX (dp), PARAMETER  :: czero = (0.0_dp, 0.0_dp)

nz = 0
az = ABS(z)
nn = n
cy = czero
dfnu = fnu + (n-1)
IF (az > 2.0_dp) THEN
  IF (az*az*0.25_dp > dfnu+1.0_dp) GO TO 10
END IF
!-----------------------------------------------------------------------
!     POWER SERIES
!-----------------------------------------------------------------------
CALL cseri(z, fnu, kode, nn, cy, nw, tol, elim, alim)
inw = ABS(nw)
nz = nz + inw
nn = nn - inw
IF (nn == 0) RETURN
IF (nw >= 0) GO TO 80
dfnu = fnu + (nn-1)

10 IF (az >= rl) THEN
  IF (dfnu > 1.0_dp) THEN
    IF (az+az < dfnu*dfnu) GO TO 20
  END IF
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z
!-----------------------------------------------------------------------
  CALL casyi(z, fnu, kode, nn, cy, nw, rl, tol, elim, alim)
  IF (nw < 0) GO TO 90
  GO TO 80
END IF
IF (dfnu <= 1.0_dp) GO TO 40
!-----------------------------------------------------------------------
!     OVERFLOW AND UNDERFLOW TEST ON I SEQUENCE FOR MILLER ALGORITHM
!-----------------------------------------------------------------------
20 CALL cuoik(z, fnu, kode, 1, nn, cy, nw, tol, elim, alim)
IF (nw < 0) GO TO 90
nz = nz + nw
nn = nn - nw
IF (nn == 0) RETURN
dfnu = fnu + (nn-1)
IF (dfnu > fnul) GO TO 70
IF (az > fnul) GO TO 70

30 IF (az > rl) GO TO 50
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES
!-----------------------------------------------------------------------
40 CALL cmlri(z, fnu, kode, nn, cy, nw, tol)
IF (nw < 0) GO TO 90
GO TO 80
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE WRONSKIAN
!-----------------------------------------------------------------------
!     OVERFLOW TEST ON K FUNCTIONS USED IN WRONSKIAN
!-----------------------------------------------------------------------
50 CALL cuoik(z, fnu, kode, 2, 2, cw, nw, tol, elim, alim)
IF (nw < 0) THEN
  nz = nn
  cy(1:nn) = czero
  RETURN
END IF
IF (nw > 0) GO TO 90
CALL cwrsk(z, fnu, kode, nn, cy, nw, cw, tol, elim, alim)
IF (nw < 0) GO TO 90
GO TO 80
!-----------------------------------------------------------------------
!     INCREMENT FNU+NN-1 UP TO FNUL, COMPUTE AND RECUR BACKWARD
!-----------------------------------------------------------------------
70 nui = fnul - dfnu + 1
nui = MAX(nui, 0)
CALL cbuni(z, fnu, kode, nn, cy, nw, nui, nlast, fnul, tol, elim, alim)
IF (nw < 0) GO TO 90
nz = nz + nw
IF (nlast /= 0) THEN
  nn = nlast
  GO TO 30
END IF

80 RETURN

90 nz = -1
IF (nw == -2) nz = -2
RETURN
END SUBROUTINE cbinu



FUNCTION gamln(z) RESULT(fn_val)

! N.B. Argument IERR has been removed.

!***BEGIN PROLOGUE  GAMLN
!***DATE WRITTEN   830501   (YYMMDD)
!***REVISION DATE  830501   (YYMMDD)
!***CATEGORY NO.  B5F
!***KEYWORDS  GAMMA FUNCTION,LOGARITHM OF GAMMA FUNCTION
!***AUTHOR  AMOS, DONALD E., SANDIA NATIONAL LABORATORIES
!***PURPOSE  TO COMPUTE THE LOGARITHM OF THE GAMMA FUNCTION
!***DESCRIPTION

!   GAMLN COMPUTES THE NATURAL LOG OF THE GAMMA FUNCTION FOR Z > 0.
!   THE ASYMPTOTIC EXPANSION IS USED TO GENERATE VALUES GREATER THAN ZMIN
!   WHICH ARE ADJUSTED BY THE RECURSION G(Z+1)=Z*G(Z) FOR Z <= ZMIN.
!   THE FUNCTION WAS MADE AS PORTABLE AS POSSIBLE BY COMPUTIMG ZMIN FROM THE
!   NUMBER OF BASE 10 DIGITS IN A WORD,
!   RLN = MAX(-LOG10(EPSILON(0.0_dp)), 0.5E-18)
!   LIMITED TO 18 DIGITS OF (RELATIVE) ACCURACY.

!   SINCE INTEGER ARGUMENTS ARE COMMON, A TABLE LOOK UP ON 100
!   VALUES IS USED FOR SPEED OF EXECUTION.

!  DESCRIPTION OF ARGUMENTS

!      INPUT
!        Z      - REAL ARGUMENT, Z > 0.0_dp

!      OUTPUT
!        GAMLN  - NATURAL LOG OF THE GAMMA FUNCTION AT Z
!        IERR   - ERROR FLAG
!                 IERR=0, NORMAL RETURN, COMPUTATION COMPLETED
!                 IERR=1, Z <= 0.0_dp,    NO COMPUTATION

!***REFERENCES  COMPUTATION OF BESSEL FUNCTIONS OF COMPLEX ARGUMENT
!                 BY D. E. AMOS, SAND83-0083, MAY 1983.
!***ROUTINES CALLED  I1MACH,R1MACH
!***END PROLOGUE  GAMLN

REAL (dp), INTENT(IN)  :: z
REAL (dp)              :: fn_val

INTEGER    :: i, i1m, k, mz, nz
REAL (dp)  :: fln, fz, rln, s, tlg, trm, tst, t1, wdtol,  &
              zdmy, zinc, zm, zmin, zp, zsq

!           LNGAMMA(N), N=1,100
REAL (dp), PARAMETER  :: gln(100) = (/ 0.00000000000000000_dp,  &
  0.00000000000000000_dp, 6.93147180559945309D-01, 1.79175946922805500_dp,   &
  3.17805383034794562_dp, 4.78749174278204599_dp, 6.57925121201010100_dp,    &
  8.52516136106541430_dp, 1.06046029027452502D+01, 1.28018274800814696D+01,  &
  1.51044125730755153D+01, 1.75023078458738858D+01, 1.99872144956618861D+01, &
  2.25521638531234229D+01, 2.51912211827386815D+01, 2.78992713838408916D+01, &
  3.06718601060806728D+01, 3.35050734501368889D+01, 3.63954452080330536D+01, &
  3.93398841871994940D+01, 4.23356164607534850D+01, 4.53801388984769080D+01, &
  4.84711813518352239D+01, 5.16066755677643736D+01, 5.47847293981123192D+01, &
  5.80036052229805199D+01, 6.12617017610020020D+01, 6.45575386270063311D+01, &
  6.78897431371815350D+01, 7.12570389671680090D+01, 7.46582363488301644D+01, &
  7.80922235533153106D+01, 8.15579594561150372D+01, 8.50544670175815174D+01, &
  8.85808275421976788D+01, 9.21361756036870925D+01, 9.57196945421432025D+01, &
  9.93306124547874269D+01, 1.02968198614513813D+02, 1.06631760260643459D+02, &
  1.10320639714757395D+02, 1.14034211781461703D+02, 1.17771881399745072D+02, &
  1.21533081515438634D+02, 1.25317271149356895D+02, 1.29123933639127215D+02, &
  1.32952575035616310D+02, 1.36802722637326368D+02, 1.40673923648234259D+02, &
  1.44565743946344886D+02, 1.48477766951773032D+02, 1.52409592584497358D+02, &
  1.56360836303078785D+02, 1.60331128216630907D+02, 1.64320112263195181D+02, &
  1.68327445448427652D+02, 1.72352797139162802D+02, 1.76395848406997352D+02, &
  1.80456291417543771D+02, 1.84533828861449491D+02, 1.88628173423671591D+02, &
  1.92739047287844902D+02, 1.96866181672889994D+02, 2.01009316399281527D+02, &
  2.05168199482641199D+02, 2.09342586752536836D+02, 2.13532241494563261D+02, &
  2.17736934113954227D+02, 2.21956441819130334D+02, 2.26190548323727593D+02, &
  2.30439043565776952D+02, 2.34701723442818268D+02, 2.38978389561834323D+02, &
  2.43268849002982714D+02, 2.47572914096186884D+02, 2.51890402209723194D+02, &
  2.56221135550009525D+02, 2.60564940971863209D+02, 2.64921649798552801D+02, &
  2.69291097651019823D+02, 2.73673124285693704D+02, 2.78067573440366143D+02, &
  2.82474292687630396D+02, 2.86893133295426994D+02, 2.91323950094270308D+02, &
  2.95766601350760624D+02, 3.00220948647014132D+02, 3.04686856765668715D+02, &
  3.09164193580146922D+02, 3.13652829949879062D+02, 3.18152639620209327D+02, &
  3.22663499126726177D+02, 3.27185287703775217D+02, 3.31717887196928473D+02, &
  3.36261181979198477D+02, 3.40815058870799018D+02, 3.45379407062266854D+02, &
  3.49954118040770237D+02, 3.54539085519440809D+02, 3.59134205369575399D+02 /)

!             COEFFICIENTS OF ASYMPTOTIC EXPANSION
REAL (dp), PARAMETER  :: cf(22) = (/  &
    8.33333333333333333D-02, -2.77777777777777778D-03,  &
    7.93650793650793651D-04, -5.95238095238095238D-04,  &
    8.41750841750841751D-04, -1.91752691752691753D-03,  &
    6.41025641025641026D-03, -2.95506535947712418D-02,  &
    1.79644372368830573D-01, -1.39243221690590112_dp,   &
    1.34028640441683920D+01, -1.56848284626002017D+02,  &
    2.19310333333333333D+03, -3.61087712537249894D+04,  &
    6.91472268851313067D+05, -1.52382215394074162D+07,  &
    3.82900751391414141D+08, -1.08822660357843911D+10,  &
    3.47320283765002252D+11, -1.23696021422692745D+13,  &
    4.88788064793079335D+14, -2.13203339609193739D+16 /)

!             LN(2*PI)
REAL (dp), PARAMETER  :: con = 1.83787706640934548_dp

!***FIRST EXECUTABLE STATEMENT  GAMLN
IF (z > 0.0_dp) THEN
  IF (z <= 101.0_dp) THEN
    nz = z
    fz = z - nz
    IF (fz <= 0.0_dp) THEN
      IF (nz <= 100) THEN
        fn_val = gln(nz)
        RETURN
      END IF
    END IF
  END IF
  wdtol = EPSILON(0.0_dp)
  wdtol = MAX(wdtol, 0.5D-18)
  i1m = DIGITS(0.0_dp)
  rln = LOG10( REAL( RADIX(0.0_dp), KIND=dp) ) * i1m
  fln = MIN(rln,20.0_dp)
  fln = MAX(fln,3.0_dp)
  fln = fln - 3.0_dp
  zm = 1.8000_dp + 0.3875_dp * fln
  mz = zm + 1
  zmin = mz
  zdmy = z
  zinc = 0.0_dp
  IF (z < zmin) THEN
    zinc = zmin - nz
    zdmy = z + zinc
  END IF
  zp = 1.0_dp / zdmy
  t1 = cf(1) * zp
  s = t1
  IF (zp >= wdtol) THEN
    zsq = zp * zp
    tst = t1 * wdtol
    DO  k = 2, 22
      zp = zp * zsq
      trm = cf(k) * zp
      IF (ABS(trm) < tst) EXIT
      s = s + trm
    END DO
  END IF

  IF (zinc == 0.0_dp) THEN
    tlg = LOG(z)
    fn_val = z * (tlg-1.0_dp) + 0.5_dp * (con-tlg) + s
    RETURN
  END IF
  zp = 1.0_dp
  nz = zinc
  DO  i = 1, nz
    zp = zp * (z + (i-1))
  END DO
  tlg = LOG(zdmy)
  fn_val = zdmy * (tlg-1.0_dp) - LOG(zp) + 0.5_dp * (con-tlg) + s
  RETURN
END IF

WRITE(*, *) '** ERROR: Zero or -ve argument for function GAMLN **'
RETURN
END FUNCTION gamln



SUBROUTINE cuchk(y, nz, ascle, tol)
!***BEGIN PROLOGUE  CUCHK
!***REFER TO CSERI,CUOIK,CUNK1,CUNK2,CUNI1,CUNI2,CKSCL

!   Y ENTERS AS A SCALED QUANTITY WHOSE MAGNITUDE IS GREATER THAN
!   EXP(-ALIM) = ASCLE = 1.0E+3*TINY(0.0_dp)/TOL.  THE TEST IS MADE TO SEE
!   IF THE MAGNITUDE OF THE REAL OR IMAGINARY PART WOULD UNDERFLOW WHEN Y IS
!   SCALED (BY TOL) TO ITS PROPER VALUE.  Y IS ACCEPTED IF THE UNDERFLOW IS AT
!   LEAST ONE PRECISION BELOW THE MAGNITUDE OF THE LARGEST COMPONENT; OTHERWISE
!   THE PHASE ANGLE DOES NOT HAVE ABSOLUTE ACCURACY AND AN UNDERFLOW IS ASSUMED.

!***ROUTINES CALLED  (NONE)
!***END PROLOGUE  CUCHK

COMPLEX (dp), INTENT(IN)  :: y
INTEGER, INTENT(OUT)      :: nz
REAL (dp), INTENT(IN)     :: ascle
REAL (dp), INTENT(IN)     :: tol

REAL (dp)  :: ss, st, yr, yi

nz = 0
yr = REAL(y, KIND=dp)
yi = AIMAG(y)
yr = ABS(yr)
yi = ABS(yi)
st = MIN(yr, yi)
IF (st > ascle) RETURN
ss = MAX(yr, yi)
st = st / tol
IF (ss < st) nz = 1
RETURN
END SUBROUTINE cuchk



SUBROUTINE cacai(z, fnu, kode, mr, n, y, nz, rl, tol, elim, alim)
!***BEGIN PROLOGUE  CACAI
!***REFER TO  CAIRY

!  CACAI APPLIES THE ANALYTIC CONTINUATION FORMULA

!      K(FNU,ZN*EXP(MP)) = K(FNU,ZN)*EXP(-MP*FNU) - MP*I(FNU,ZN)
!              MP = PI*MR*CMPLX(0.0,1.0)

!  TO CONTINUE THE K FUNCTION FROM THE RIGHT HALF TO THE LEFT HALF Z PLANE
!  FOR USE WITH CAIRY WHERE FNU=1/3 OR 2/3 AND N=1.
!  CACAI IS THE SAME AS CACON WITH THE PARTS FOR LARGER ORDERS AND
!  RECURRENCE REMOVED.  A RECURSIVE CALL TO CACON CAN RESULT IF CACON
!  IS CALLED FROM CAIRY.

!***ROUTINES CALLED  CASYI,CBKNU,CMLRI,CSERI,CS1S2,R1MACH
!***END PROLOGUE  CACAI

COMPLEX (dp), INTENT(IN)   :: z
REAL (dp), INTENT(IN)      :: fnu
INTEGER, INTENT(IN)        :: kode
INTEGER, INTENT(IN OUT)    :: mr
INTEGER, INTENT(IN)        :: n
COMPLEX (dp), INTENT(OUT)  :: y(n)
INTEGER, INTENT(OUT)       :: nz
REAL (dp), INTENT(IN)      :: rl
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(IN)      :: elim
REAL (dp), INTENT(IN)      :: alim

COMPLEX (dp)  :: csgn, cspn, c1, c2, zn, cy(2)
REAL (dp)     :: arg, ascle, az, cpn, dfnu, fmr, sgn, spn, yy
INTEGER       :: inu, iuf, nn, nw
REAL (dp), PARAMETER  :: pi = 3.14159265358979324_dp

nz = 0
zn = -z
az = ABS(z)
nn = n
dfnu = fnu + (n-1)
IF (az > 2.0_dp) THEN
  IF (az*az*0.25_dp > dfnu+1.0_dp) GO TO 10
END IF
!-----------------------------------------------------------------------
!     POWER SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
CALL cseri(zn, fnu, kode, nn, y, nw, tol, elim, alim)
GO TO 20

10 IF (az >= rl) THEN
!-----------------------------------------------------------------------
!     ASYMPTOTIC EXPANSION FOR LARGE Z FOR THE I FUNCTION
!-----------------------------------------------------------------------
  CALL casyi(zn, fnu, kode, nn, y, nw, rl, tol, elim, alim)
  IF (nw < 0) GO TO 30
ELSE
!-----------------------------------------------------------------------
!     MILLER ALGORITHM NORMALIZED BY THE SERIES FOR THE I FUNCTION
!-----------------------------------------------------------------------
  CALL cmlri(zn, fnu, kode, nn, y, nw, tol)
  IF (nw < 0) GO TO 30
END IF
!-----------------------------------------------------------------------
!     ANALYTIC CONTINUATION TO THE LEFT HALF PLANE FOR THE K FUNCTION
!-----------------------------------------------------------------------
20 CALL cbknu(zn, fnu, kode, 1, cy, nw, tol, elim, alim)
IF (nw == 0) THEN
  fmr = mr
  sgn = -SIGN(pi, fmr)
  csgn = CMPLX(0.0_dp, sgn, KIND=dp)
  IF (kode /= 1) THEN
    yy = -AIMAG(zn)
    cpn = COS(yy)
    spn = SIN(yy)
    csgn = csgn * CMPLX(cpn, spn, KIND=dp)
  END IF
!-----------------------------------------------------------------------
!     CALCULATE CSPN=EXP(FNU*PI*I) TO MINIMIZE LOSSES OF SIGNIFICANCE
!     WHEN FNU IS LARGE
!-----------------------------------------------------------------------
  inu = fnu
  arg = (fnu - inu) * sgn
  cpn = COS(arg)
  spn = SIN(arg)
  cspn = CMPLX(cpn, spn, KIND=dp)
  IF (MOD(inu,2) == 1) cspn = -cspn
  c1 = cy(1)
  c2 = y(1)
  IF (kode /= 1) THEN
    iuf = 0
    ascle = 1.0E+3 * TINY(0.0_dp) / tol
    CALL cs1s2(zn, c1, c2, nw, ascle, alim, iuf)
    nz = nz + nw
  END IF
  y(1) = cspn * c1 + csgn * c2
  RETURN
END IF

30 nz = -1
IF (nw == -2) nz = -2
RETURN
END SUBROUTINE cacai

END Module Complex_Bessel
