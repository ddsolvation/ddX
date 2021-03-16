!***************************************************************
!*     Purpose: This program computes the modified spherical   *
!*              Bessel functions of the first kind in(x) and   *
!*              in'(x) using subroutine SPHI_bessel                   *
!* ----------------------------------------------------------- *
!*     Input :  x --- Argument of in(x)                        *
!*              n --- Order of in(x) ( 0 to 250 )              *
!*     Output:  SI(n) --- in(x)                                *
!*              DI(n) --- in'(x)                               *
!*     Example: x = 10.0                                       *
!*                n          in(x)               in'(x)        *
!*              --------------------------------------------   *
!*                0     .1101323287D+04     .9911909633D+03    *
!*                1     .9911909633D+03     .9030850948D+03    *
!*                2     .8039659985D+03     .7500011637D+03    *
!*                3     .5892079640D+03     .5682828129D+03    *
!*                4     .3915204237D+03     .3934477522D+03    *
!*                5     .2368395827D+03     .2494166741D+03    *
!* ----------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special     *
!*             Functions,                                      *
!*             jin.ece.uiuc.edu/routines/routines.html".       *
!*                                                             *
!*                           F90 Release By J-P Moreau, Paris. *
!*                                  (www.jpmoreau.fr)          *
!***************************************************************
module bessel
implicit none

real*8, parameter :: tol_zero = 1.0D-15
real*8, parameter :: tol_inf = 1.0D30
contains

        

!      ========================================================
!      Purpose: Compute modified spherical Bessel functions
!               of the first kind, in(x) and in'(x)
!      Input :  
!               x --- Argument of in(x)
!               n --- Order of in(x) ( n = 0,1,2,... )
!      Output:  
!               SI(n) --- in(x)
!               DI(n) --- in'(x)
!               NM --- Highest order computed
!      Routines called:
!               MSTA1 and MSTA2 for computing the starting
!               point for backward recurrence
!      ========================================================
      SUBROUTINE SPHI_bessel(N,X,NM,SI,DI)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, intent(in) :: N
        INTEGER, INTENT(INOUT) :: NM
        double precision,  intent(in) :: X
        double precision, DIMENSION(0:N), intent(inout) :: SI, DI
        
        double precision:: CS, F, F0, F1, SI0
        INTEGER :: K, M, MP        
        logical :: flag
                
        NM=N
        IF (DABS(X).LT.1.0D-100) THEN
           DO 10 K=0,N
              SI(K)=0.0D0
10            DI(K)=0.0D0
           SI(0)=1.0D0
           DI(1)=0.333333333333333D0
           RETURN
        ENDIF
        
        SI(0)=DSINH(X)/X
        if (N .eq. 0) then
            return
        end if
        SI(1)=-(DSINH(X)/X-DCOSH(X))/X
        
        if (SI(0) .gt. tol_inf ) then
            SI = tol_inf + 1
            DI = tol_inf+1
            !write(*,*) 'Bessel error: too large input value for the bessel function!' 
            return
        end if
        
        if (X**N .lt. 1.0D0-6) then
            flag = .false.
        else
            flag = .true.
            !flag = .false.
        end if
        
        if (flag) then
        
            F0 = SI(0)
            F1 = SI(1)
            if (N .ge. 2) then
                DO K= 2, N
                    F= -(2.0D0*K-1.0D0)*F1/X+F0
                    SI(K)=F
                    F0=F1
                    F1=F
                END DO
             end if
             
        else
            SI0=SI(0)
            IF (N.GE.2) THEN
               !M= 50 !
               M = MSTA1(X,200)
               IF (M.LT.N) THEN
                  NM=M
               ELSE
                  M=MSTA2(X,N,15)
               ENDIF
               F0=0.0D0
               F1=1.0D0-100
               DO 15 K=M,0,-1
                  F=(2.0D0*K+3.0D0)*F1/X+F0
                  IF (K.LE.NM) SI(K)=F
                  F0=F1
15              F1=F
               CS=SI0/F
               DO 20 K=0,NM
20             SI(K)=CS*SI(K)
            ENDIF
            
        end if

        DI(0)=SI(1)
        DO 25 K=1,NM
25         DI(K)=SI(K-1)-(K+1.0D0)/X*SI(K)
        
!          write(*,*)  M , X, DSINH(X)/X, F
!         write(*,*)  'SI0 = ', SI0
!         write(*,*)  'CS = ', CS
          !write(*,*) 'SI = ', SI
          !write(*,*) 'DI = ', DI
          !stop
        
        END subroutine SPHI_bessel

        
!***************************************************************
!
!*       Purpose: This program computes the modified spherical    *
!*                Bessel functions kn(x) and kn'(x) using         *
!*                subroutine SPHK_bessel                                 *
!*       Input :  x --- Argument of kn(x)  ( x > 0 )              *
!*                n --- Order of kn(x) ( n = 0 to 250 )           *
!*       Output:  SK(n) --- kn(x)                                 *
!*                DK(n) --- kn'(x)                                *
!*       Example: x= 10.0                                         *
!*                  n          kn(x)               kn'(x)         *
!*                --------------------------------------------    *
!*                  0     .7131404291D-05    -.7844544720D-05     * 
!*                  1     .7844544720D-05    -.8700313235D-05     *
!*                  2     .9484767707D-05    -.1068997503D-04     * 
!*                  3     .1258692857D-04    -.1451953914D-04     *
!*                  4     .1829561771D-04    -.2173473743D-04     *
!*                  5     .2905298451D-04    -.3572740841D-04     *
!* -------------------------------------------------------------- *
!* REFERENCE: "Fortran Routines for Computation of Special        *
!*             Functions,                                         *
!*             jin.ece.uiuc.edu/routines/routines.html".          *
!*                                                                *
!*                            F90 Release By J-P Moreau, Paris.   *
!*                                   (www.jpmoreau.fr)            *
!******************************************************************

        

!       =====================================================
!       Purpose: Compute modified spherical Bessel functions
!                of the second kind, kn(x) and kn'(x)
!       Input :  x --- Argument of kn(x)  ( x ò 0 )
!                n --- Order of kn(x) ( n = 0,1,2,... )
!       Output:  SK(n) --- kn(x)
!                DK(n) --- kn'(x)
!                NM --- Highest order computed
!       =====================================================
      SUBROUTINE SPHK_bessel(N,X,NM,SK,DK)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
         integer, intent(in) :: N
         INTEGER, INTENT(INOUT) :: NM
        DIMENSION SK(0:N),DK(0:N)
        INTEGER :: K
        
        PI=3.141592653589793D0
        NM=N
        
        IF (X.LT.1.0D-60) THEN
           DO 10 K=0,N
              SK(K)=1.0D+300
10            DK(K)=-1.0D+300
           RETURN
        ENDIF
        SK(0)= 1.0D0/X*DEXP(-X) !!! I HAVE CHANGED HERE
        SK(1)=SK(0)*(1.0D0+1.0D0/X)
        F0=SK(0)
        F1=SK(1)
        DO 15 K=2,N
           F=(2.0D0*K-1.0D0)*F1/X+F0
           SK(K)=F
           IF (DABS(F).GT.1.0D+300) GO TO 20
           F0=F1
15         F1=F
20      NM=K-1
        DK(0)=-SK(1)
        DO 25 K=1,NM
25         DK(K)=-SK(K-1)-(K+1.0D0)/X*SK(K)
        RETURN
        END subroutine SPHK_bessel

!      ===================================================
!      Purpose: Determine the starting point for backward  
!               recurrence such that the magnitude of    
!               Jn(x) at that point is about 10^(-MP)
!      Input :  x     --- Argument of Jn(x)
!               MP    --- Value of magnitude
!      Output:  MSTA1 --- Starting point   
!      ===================================================
        INTEGER FUNCTION MSTA1(X,MP)
        
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, intent(in) :: MP
        real*8, intent(in) :: X
        integer :: IT, N0, N1, NN
        
        A0=DABS(X)
        N0=INT(1.1*A0)+1
        F0=ENVJ(N0,A0)-MP
        N1=N0+5
        F1=ENVJ(N1,A0)-MP
        DO 10 IT=1,20             
           NN=N1-(N1-N0)/(1.0D0-F0/F1)                  
           F=ENVJ(NN,A0)-MP
           IF(ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
 10        F1=F
 20     MSTA1=NN
        RETURN
        
        END


        
!      ===================================================
!      Purpose: Determine the starting point for backward
!               recurrence such that all Jn(x) has MP
!               significant digits
!      Input :  x  --- Argument of Jn(x)
!               n  --- Order of Jn(x)
!               MP --- Significant digit
!      Output:  MSTA2 --- Starting point
!      ===================================================
      INTEGER FUNCTION MSTA2(X,N,MP)
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        integer, intent(in) :: N, MP
        real*8, intent(in) :: X
        integer :: IT, N0, N1, NN
        
        A0=DABS(X)
        HMP=0.5D0*MP
        EJN=ENVJ(N,A0)
        IF (EJN.LE.HMP) THEN
           OBJ=MP
           N0=INT(1.1*A0)
        ELSE
           OBJ=HMP+EJN
           N0=N
        ENDIF
        F0=ENVJ(N0,A0)-OBJ
        N1=N0+5
        F1=ENVJ(N1,A0)-OBJ
        DO 10 IT=1,20
           NN=N1-(N1-N0)/(1.0D0-F0/F1)
           F=ENVJ(NN,A0)-OBJ
           IF (ABS(NN-N1).LT.1) GO TO 20
           N0=N1
           F0=F1
           N1=NN
10         F1=F
20      MSTA2=NN+10
        RETURN
        
        END

        REAL*8 FUNCTION ENVJ(N,X)
        integer, intent(in) :: N
        DOUBLE PRECISION X
        ENVJ=0.5D0*DLOG10(6.28D0*N)-N*DLOG10(1.36D0*X/N)
        RETURN
        END
        

end module bessel

