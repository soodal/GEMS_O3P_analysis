MODULE bounded_nonlin_LS

  ! The systems are described in a user's guide UMINF-109/110.84
  ! by Per Lindstroem, Institute of Information Processing,
  ! University of Umea, Sweden.

  ! More details can be found in
  ! Wedin, Lindstroem: Methods and Software for Nonlinear Least
  ! Squares Problems, UMINF-133.87 (rev. July 1988, 1993) and in "Gauss-Newton
  ! based algorithms for constrained least squares problems" of which the last
  ! one can be downloaded from
  ! http://www.cs.umu.se/~perl/reports/alg.ps.gz

  ! Remember: Anything that comes free has no guarantee.

  ! Code converted using TO_F90 by Alan Miller
  ! Date: 2000-04-29  Time: 15:46:47

  ! TPK addition
  ! ============
  USE OMSAO_precision_module, ONLY : dp
  IMPLICIT NONE
  !INTEGER, PARAMETER    :: dp = KIND (1.0D0) !SELECTED_REAL_KIND(12, 60)

  REAL (dp), PARAMETER  :: zero = 0.0_dp, one = 1.0_dp

  !     COMMON VARIABLES CONTAINING INFORMATION CONCERNING THE PREVIOUS
  !     2 POINTS.THE SUFFICES KM2 AND KM1 IN THE NAMES OF THE VARIABLES
  !     REPRESENT TIME STEP K-2 AND K-1 RESPECTIVELY.
  !     THESE VARIABLES ARE UPDATED ONLY INSIDE THE ROUTINE EVREUC

  ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
  !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1
  INTEGER, SAVE    :: rngkm2,kodkm2,rngkm1,kodkm1
  REAL (dp), SAVE  :: betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,betkm1,d1km1,  &
       fsqkm1,dxnkm1,alfkm1,aupkm1

  !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
  !     SRELPR = SINGLE RELATIVE PRECISION

  ! COMMON /machin/ srelpr
  REAL (dp), SAVE  :: srelpr

  !     COMMON VARIABLES USED FOR RESTART INFORMATION

  ! COMMON /back/ icount,lattry,itotal,philat,bestpg,irank
  ! COMMON /INDEX/ imax,ival
  INTEGER, SAVE    :: icount,lattry,itotal,irank,imax,ival
  REAL (dp), SAVE  :: philat,bestpg

  ! COMMON /negdir/ ifree,indic
  INTEGER, SAVE    :: ifree,indic



CONTAINS


  !ELSUNC

  SUBROUTINE elsunc(x,n,mdc,m,ffunc,bnd,bl,bu, p,w, EXIT,f,c)

    ! N.B. Argument MDW has been removed.

    REAL (dp), INTENT(IN OUT)  :: x(:) !x(1:n)!
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: mdc
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: bnd
    REAL (dp), INTENT(IN OUT)  :: bl(:) !bl(1:n)!
    REAL (dp), INTENT(IN OUT)  :: bu(:) !bu(1:n) !
    !tpk INTEGER, INTENT(OUT)       :: p(11)
    INTEGER, INTENT(INOUT)       :: p(11) !tpk
    !tpk REAL (dp), INTENT(OUT)     :: w(6)
    REAL (dp), INTENT(INOUT)     :: w(6)  !tpk
    INTEGER, INTENT(OUT)       :: EXIT
    REAL (dp), INTENT(OUT)     :: f(:) !f(1:m)!
    REAL (dp), INTENT(OUT)     ::  c(:,:) !c(mdc,n)!

    EXTERNAL ffunc

    
    !   ***********************************************************
    !   * THIS IS AN EASY-TO-USE VERSION OF THE SUBROUTINE LSUNC. *
    !   * LSUNC IS DEVELOPED BY PER LINDSTR\M AND PER-]KE WEDIN   *
    !   * AT THE INSTITUTE OF INFORMATION PROCESSING UNIVERSITY   *
    !   * OF UME], S-90187 UME], SWEDEN                           *
    !   ***********************************************************

    !   THE FOLLOWING PARAMETERS ARE GIVEN DEFAULT VALUES
    !   PROVIDED THE USER HAS GIVEN A NEGATIVE VALUE TO THE
    !   CORRESPONDING LOCATIONS IN THE AREAS P AND W RESPECTIVELY.

    !          IPRINT=1  (WRITE EVERY STEP)
    !          NOUT=10  (WRITING IS DONE ON UNIT NO. 10)
    !          MAXIT=20*N  (20 TIMES THE NO. OF PARAMETERS)
    !          SEC=TRUE  (THE METHOD OF NEWTON IS ALLOWED AT THE
    !                     END OF THE ITERATION IN SOME SITUATIONS
    !                     (LARGE RESIDUALS)  )
    !          SCALE=0  (INTERNAL SCALING IS NOT USED TO COMPUTE
    !                    THE GAUSS-NEWTON SEARCH DIRECTION )
    !          TOL=SQRT(SINGLE RELATIVE PRECISION)
    !          EPSREL=SQRT(SINGLE RELATIVE PRECISION)
    !          EPSABS=SINGLE RELATIVE PRECISION
    !          EPSX=SQRT(SINGLE RELATIVE PRECISION)

    !   PURPOSE...
    !   SOLVE THE NONLINEAR LEAST SQUARES PROBLEM

    !           MINIMIZE  PHI(X) = 0.5*|| F(X) ||**2
    !           SUBJECT TO
    !                      BL(I) <= X(I) <= BU(I)  I=1,2,.....,N


    !   WHERE                                T
    !        F(X)= (F (X) F (X)...... F (X) )
    !                1     2           M

    !   AND X IS N-DIMENSIONAL.

    !   ON ENTRY

    !   X()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        A FIRST APPROXIMATION OF THE PARAMETERS (UNKNOWNS)
    !   N    INTEGER SCALAR CONTAINING THE NUMBER OF UNKNOWNS
    !   MDC  INTEGER SCALAR CONTAINING LEADING DIMENSION OF THE ARRAY C
    !        (MDC MUST BE >= M)
    !   M    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAY F
    !   THE FOLLOWING 9 PARAMETERS ARE GIVEN DEFAULT VALUES PROVIDED
    !   THE CORRESPONDING LOCATION IS < 0 ON ENTRY.  THE DEFAULT VALUES
    !   ARE GIVEN IN BRACKETS TOGETHER WITH THE ORIGINAL NAME OF
    !   THE PARAMETER.  WHERE APPLICABLE,THESE DEFAULT VALUES ARE RETURNED
    !   IN THE CORRESPONDING LOCATION WHEN THE LOCATION IS <0 ON ENTRY.

    !   P(1) (IPRINT=1) STEP BETWEEN WRITING
    !   P(2) (NOUT=10) FORTRAN UNIT FOR WRITING
    !   P(3) (MAXIT=20*N) MAXIMUM NO. OF ALLOWED ITERATIONS
    !   P(4) (NEWTON ALLOWED)INDICATOR FOR ALLOWING THE METHOD OF NEWTON
    !   P(5) (NO INTERNAL) INTERNAL SCALING IS NOT USED AS DEFAULT
    !   W(1) (TOL=SQRT(SRELPR)) PSEUDO RANK TOLERANCE CONSTANT
    !   W(2) (EPSREL=SQRT(SRELPR)) RELATIVE CONVERGENCE CONSTANT
    !   W(3) (EPSABS=SRELPR) ABSOLUTE CONVERGENCE CONSTANT
    !   W(4) (EPSX=SQRT(SRELPR)) PARAMETER CONVERGENCE CONSTANT
    !   FFUNC        SUBROUTINE NAME-THAT MUST BE DECLEARED EXTERNAL IN
    !                THE CALLING PROGRAM.  A ROUTINE NAMED FFUNC IS USED
    !                TO EVALUATE THE FUNCTION F(X) AND/OR THE JACOBIAN AT
    !                A CERTAIN POINT X AND SHOULD BE WRITTEN AS FOLLOWS

    !                SUBROUTINE FFUNC(X,N,F,M,CTRL,C,MDC)
    !                INTEGER N,M,CTRL,MDC
    !                REAL X(N),F(M),C(MDC,N)
    !                -----------------------
    !                CTRL CAN HAVE 3 DIFFERENT VALUES ON ENTRY
    !       CTRL= 1  MEANS EVALUATE THE FUNCTIONS AT THE POINT X AND
    !                RETURN THIS VECTOR IN THE ARRAY F IF THE FUNCTIONS
    !                ARE COMPUTABLE.
    !                ON RETURN THE USER CAN INDICATE UNCOMPUTABILITY BY
    !                SETTING CTRL=-1
    !                DO NOT ALTER ARRAY X.
    !       CTRL=-1  MEANS EVALUATE THE FUNCTIONS AT THE POINT X AND
    !                RETURN THIS VECTOR IN THE ARRAY F IF THE FUNCTIONS
    !                ARE COMPUTABLE. DO NOT ALTER ARRAY X.
    !                POSSIBLE UNCOMPUTABILITY OF THE FUNCTIONS MUST BE
    !                INDICATED BY SETTING CTRL TO A VALUE <-10 ON RETURN.
    !       CTRL= 2  MEANS CALCULATE THE JACOBIAN OF F(X) AT THE POINT X
    !                AND RETURN THIS MATRIX IN THE ARRAY C IF THE JACOBIAN
    !                IS SUPPLIED ANALYTICALLY.
    !                POSSIBLE UNCOMPUTABILITY OF THE JACOBIAN MUST BE
    !                INDICATED BY SETTING CTRL TO A VALUE <-10 ON RETURN.
    !                IF THE USER WANTS THE JACOBIAN BEING COMPUTED
    !                NUMERICALLY THAT SHOULD BE INDICATED BY SETTING
    !                CTRL=0 ON RETURN.
    !                DO NOT ALTER ARRAYS X AND F.
    !                ------------------------------
    !                RETURN
    !                END
    !   BND  INTEGER SCALAR CONTAINING A CODE FOR ASSESSING
    !        THE BOUNDS.
    !        BND = 0 MEANS AN UNCONSTRAINED PROBLEM
    !            = 1 MEANS THE SAME LOWER BOUND FOR ALL UNKNOWNS AND
    !                THE SAME UPPER BOUND FOR ALL UNKNOWNS. THE LOWER
    !                AND UPPER BOUND MUST BE STORED IN BL(1) AND BU(1)
    !                RESPECTIVELY.
    !            = ANYTHING BUT 0 AND 1 MEANS THAT THE BOUNDS MUST BE
    !              SUPPLIED BY THE USER IN BL(I) AND BU(I) I=1,2,....,N
    !   BL() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE LOWER BOUNDS OF THE UNKNOWNS
    !   BU() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE UPPER BOUNDS OF THE UNKNOWNS

    !   ON RETURN     AND EXIT.NE.-1 AND EXIT.GE.-10

    !   X()  CONTAINS THE LATEST (BEST) ESTIMATE OF THE SOLUTION POINT

    !   THE CONVERGENCE CRITERIA ARE

    !   1) RELATIVE PREDICTED REDUCTION IN THE OBJECTIVE FUNCTION
    !      IS LESS THAN EPSREL**2
    !   2) THE SUM OF SQUARES IS LESS THAN EPSABS**2
    !   3) THE RELATIVE CHANGE IN X IS LESS THAN EPSX
    !   4) WE ARE COMPUTING AT NOISE LEVEL
    !      THE LAST DIGIT IN THE CONVERGENCE CODE (SEE BELOW) INDICATES
    !      HOW THE LAST STEPS WERE COMPUTED
    !      = 0 NO TROUBLE (GAUSS-NEWTON THE LAST 3 STEPS)
    !      = 1 PRANK<N AT THE TERMINATION POINT
    !      = 2 THE METHOD OF NEWTON WAS USED (AT LEAST) IN THE LAST STEP
    !      = 3 THE 2:ND BUT LAST STEP WAS SUBSPACE MINIMIZATION BUT THE
    !          LAST TWO WERE GAUSS-NEWTON STEPS
    !      = 4 THE STEPLENGTH WAS NOT UNIT IN BOTH THE LAST TWO STEPS


    !    THE ABNORMAL TERMINATION CRITERIA ARE

    !   5) NO. OF ITERATIONS HAS EXCEEDED MAXIMUM ALLOWED ITERATIONS
    !   6) THE HESSIAN EMANATING FROM 2:ND ORDER METHOD IS NOT POS DEF
    !   7) THE ALGORITHM WOULD LIKE TO USE 2:ND DERIVATIVES BUT IS
    !      NOT ALLOWED TO DO THAT
    !   8) AN UNDAMPED STEP WITH NEWTONS METHOD IS A FAILURE
    !   9) THE LATEST SEARCH DIRECTION COMPUTED USING SUBSPACE
    !      MINIMIZATION WAS NOT A DESCENT DIRECTION (PROBABLY CAUSED
    !      BY WRONGLY COMPUTED JACOBIAN)

    !   EXIT INTEGER SCALAR THAT INDICATES WHY THE RETURN IS TAKEN
    !        =10000  CONVERGENCE DUE TO CRITERION NO. 1
    !        = 2000  CONVERGENCE DUE TO CRITERION NO. 2
    !        =  300  CONVERGENCE DUE TO CRITERION NO. 3
    !        =   40  CONVERGENCE DUE TO CRITERION NO. 4
    !        =    X   WHERE X EQUALS 0,1,2,3 OR 4

    !        <0   INDICATES THAT NO CONVERGENCE CRITERION IS FULFILLED
    !             BUT SOME ABNORMAL TERMINATION CRITERION IS SATISFIED
    !        =  -1  IF M<N OR N<=0 OR M<=0 OR MDC<M OR MDW<N*N+5*N+3*M+6
    !                   OR MAXIT<=0 OR EPSREL<0 OR EPSABS<0 OR EPSX<0
    !                   OR INVALID STARTING POINT   ON ENTRY
    !        =  -2   TERMINATION DUE TO CRITERION NO. 5
    !        =  -3   TERMINATION DUE TO CRITERION NO. 6
    !        =  -4   TERMINATION DUE TO CRITERION NO. 7
    !        =  -5   TERMINATION DUE TO CRITERION NO. 8
    !        =  -6   TERMINATION DUE TO CRITERION NO. 9
    !        =  -7   THERE IS ONLY ONE FEASIBLE POINT, NAMELY
    !                X(I)=BL(I)=BU(I)  ; I=1,2,.....,N
    !        <  -10  TERMINATION DUE TO USER STOP INDICATOR

    !   P(6) (ITER)INTEGER SCALAR CONTAINING THE NO. OF ITERATIONS UNTIL
    !        TERMINATION
    !   P(7) (FUNCEV)INTEGER SCALAR CONTAINING THE TOTAL NO. OF FUNCTION
    !        EVALUATIONS DONE INSIDE THIS ROUTINE
    !   P(8) (JACEV)INTEGER SCALAR CONTAINING THE NO. OF FUNCTION
    !        EVALUATIONS CAUSED BY COMPUTING JACOBIANS WITH DIFFERENCE
    !        METHODS
    !   P(9) (SECEV)INTEGER SCALAR CONTAINING THE NO. OF FUNCTION
    !        EVALUATIONS CAUSED BY USING THE METHOD OF NEWTON
    !   P(10) (LINEV)INTEGER SCALAR CONTAINING THE NO. OF FUNCTION
    !         EVALUATIONS CAUSED BY THE LINESEARCH ALGORITHM
    !   P(11) (PRANK)INTEGER SCALAR CONTAINING THE PSEUDO RANK OF THE MATRIX
    !         JHAT AT THE TERMINATION POINT X.  JHAT IS THE JACOBIAN MATRIX WITH
    !         COLUMNS CORRESPONDING TO ACTIVE CONSTRAINTS DELETED
    !   W(5) (PHI)REAL SCALAR CONTAINING THE VALUE OF THE OBJECTIVE
    !        FUNCTION AT THE TERMINATION POINT X
    !   W(6) (SPEED)REAL SCALAR CONTAINING AN ESTIMATE OF THE LINEAR
    !        CONVERGENCE FACTOR
    !   F()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTAINING THE
    !        VALUE OF THE RESIDUALS F(I) AT THE TERMINATION POINT X
    !   C(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDC*MAX(4,N)
    !        (MDC>=M) CONTAINING THE MAIN PART OF THE COVARINCE MATRIX
    !            T   -1
    !          (J *J)     WHERE J IS THE JACOBIAN MATRIX COMPUTED AT X

    !   WORKING AREAS

    !   W()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION MDW
    !        (MDW MUST BE >= N*N+5*N+3*M+6 WHEN NEWTON IS ALLOWED)
    !        (MDW MUST BE >= 6*N+3*M+6 WHEN NEWTON IS NOT ALLOWED)
    !   P()  INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION 11+2*N

    !   COMMON VARIABLES CONTAINING INFORMATION CONCERNING THE PREVIOUS
    !   2 POINTS.THE SUFFICES KM2 AND KM1 IN THE NAMES OF THE VARIABLES
    !   REPRESENT TIME STEP K-2 AND K-1 RESPECTIVELY.
    !   THESE VARIABLES ARE UPDATED ONLY INSIDE THE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    !     COMMON VARIABLES USED FOR RESTART INFORMATION

    ! COMMON /back/ icount,lattry,itotal,philat,bestpg,irank
    ! COMMON /INDEX/ imax,ival

    !     THIS PROGRAM PACKAGE USES THE FOLLOWING LINPACK AND BLAS ROUTINES

    !     SQRDC,SQRSL,STRSL,SCHDC,SPOSL
    !     SAXPY,SDOT,SSCAL,SSWAP,SNRM2,SCOPY

    !     THE FOLLOWING MINPACK-1 ROUTINE IS USED

    !     COVAR

    !     INTERNAL VARIABLES

    INTEGER   :: mdg,scale,i,j
    REAL (dp) :: rootsp,bl1,bu1
    LOGICAL   :: sec
    REAL (dp), PARAMETER  :: big = 1.0E32_dp

    !     THE WORKING AREA W() WAS STRUCTURED IN THE FOLLOWING WAY

    !     TOL=W(1)
    !     EPSREL=W(2)
    !     EPSABS=W(3)
    !     EPSX=W(4)
    !     PHI=W(5)
    !     SPEED=W(6)
    !     DX = W(7).....W(N+6)
    !      G = W(N+7)....W(2*N+6)
    !     W0 = W(2*N+7).....W(3*N+6)
    !     W1 = W(3*N+7)....W(3*N+M+6)
    !     W2 = W(3*N+M+7).....W(3*N+2*M+6)
    !     W3 = W(3*N+2*M+7).......W(4*N+2*M+6)
    !      V = W(4*N+2*M+7).......W(4*N+3*M+6)
    !     DIAG = W(4*N+3*M+7).....W(5*N+3*M+6)
    !     GMAT = W(5*N+3*M+7)......

    REAL (dp)  :: dx(n), g(n), w0(n), w1(m), w2(m), w3(n), v(m), diag(n),  &
         gmat(n,n)
    INTEGER    :: ip(n), jp(n)

    ! START of code

    !     VALIDATE SOME INPUT PARAMETERS

    EXIT=0
    ! IF(p(4) >= 0 .AND. mdw < 6*n + 3*m+6) EXIT=-1
    ! IF(p(4) < 0 .AND. mdw < n*n + 5*n + 3*m + 6) EXIT=-1
    ! IF(EXIT < 0) RETURN

    !     INITIATE THE BOUND ARRAYS

    bl1=-big
    bu1=big
    IF(bnd == 1) THEN
       bl1=bl(1)
       bu1=bu(1)
    END IF
    IF(bnd == 0 .OR. bnd == 1) THEN
       DO  i=1,n
          bl(i)=bl1
          bu(i)=bu1
       END DO
    END IF

    !     COMPUTE SINGLE RELATIVE PRECISION

    CALL releps(srelpr)
    rootsp=SQRT(srelpr)

    !     SET DEFAULT VALUES OF MISSING PARAMETERS

    mdg=n
    IF(p(1) < 0) p(1)=1
    IF(p(2) < 0) p(2)=10
    IF(p(3) < 0) p(3)=20*n
    sec=p(4) < 0
    scale=1
    IF(p(5) < 0) scale=0
    IF(w(1) < zero) w(1)=rootsp
    IF(w(2) < zero) w(2)=10.0_dp*rootsp
    IF(w(3) < zero) w(3)=srelpr
    IF(w(4) < zero) w(4)=10.0_dp*rootsp

    !     CALL LSUNC TO MINIMIZE

    CALL lsunc(x,n,mdc,mdg,m,w(1),w(2),w(3),w(4),sec,p(1),p(2),p(3),  &
         scale,ffunc,bl,bu,bnd,EXIT,w(5),p(6),p(7),p(8),p(9),p(10),p(11), &
         ip,f,c,diag,w(6),jp,dx,g,w0,w1,w2,w3,v,gmat)

    IF(EXIT < 0) RETURN
    IF(MOD(EXIT,10) == 1 .AND. p(1) > 0 .AND. p(2) > 0) WRITE(p(2),1000)
1000 FORMAT(//t21, 31('*')/ t21, '* pseudo-rank of the jacobian *'/  &
         t21, '* is NOT full ( <(n-nract) ). *'/  &
         t21, '* a run with internal scaling *'/  &
         t21, '* activated(p(5)>=0 on ENTRY) *'/  &
         t21, '* might give another solution *'/  &
         t21, 31('*'))
   
    CALL covar(n,c,ip,w(1),dx)
    IF(scale <= 0) RETURN

    !     COMPENSATE FOR THE INTERNAL SCALING BY FORMING
    !          C = D*C*D
    !     WHERE D IS A DIAGONAL MATRIX

    DO  j=1,n
       c(1:n,j)=c(1:n,j)*diag(j)
    END DO
    DO  i=1,n
       c(i,1:n)=c(i,1:n)*diag(i)
    END DO

    RETURN
  END SUBROUTINE elsunc


  !COVAR

  SUBROUTINE covar(n, r, ipvt, tol, wa)

    ! Argument LDR has been removed.

    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: r(:,:)   ! r(ldr,n)
    INTEGER, INTENT(IN)     :: ipvt(:)
    REAL (dp), INTENT(IN)   :: tol
    REAL (dp), INTENT(OUT)  :: wa(:)

    !     THIS ROUTINE IS STOLEN FROM MINPACK-1

    INTEGER   :: i,ii,j,jj,k,km1,l
    LOGICAL   :: sing
    REAL (dp) :: temp,tolr

    !     FORM THE INVERSE OF R IN THE FULL UPPER TRIANGLE OF R.

    tolr=tol*ABS(r(1,1))
    l=0
    DO  k=1,n
       IF(ABS(r(k,k)) <= tolr) EXIT
       r(k,k)=one/r(k,k)
       km1=k-1
       DO  j=1,km1
          temp=r(k,k)*r(j,k)
          r(j,k)=zero
          DO  i=1,j
             r(i,k)=r(i,k) - temp*r(i,j)
          END DO
       END DO

       l=k
    END DO

    !     FORM THE FULL UPPER TRIANGLE OF THE INVERSE OF (R TRANSPOSE)*R
    !     IN THE FULL UPPER TRIANLE OF R.

    DO  k=1,l
       km1=k-1
       DO  j=1,km1
          temp=r(j,k)
          DO  i=1,j
             r(i,j)=r(i,j) + temp*r(i,k)
          END DO
       END DO

       temp=r(k,k)
       DO  i=1,k
          r(i,k)=temp*r(i,k)
       END DO
    END DO

    !     FORM THE FULL LOWER TRIANGLE OF THE COVARIANCE MATRIX
    !     IN THE STRICT LOWER TRIANGLE OF R AND IN WA

    DO  j=1,n
       jj=ipvt(j)
       sing=j > l
       DO  i=1,j
          IF(sing) r(i,j)=zero
          ii=ipvt(i)
          IF(ii > jj) r(ii,jj)=r(i,j)
          IF(ii < jj) r(jj,ii)=r(i,j)
       END DO
       wa(jj)=r(j,j)
    END DO

    !     SYMMETRIZE THE COVARIANCE MATRIX IN R.

    DO  j=1,n
       DO  i=1,j
          r(i,j)=r(j,i)
       END DO
       r(j,j)=wa(j)
    END DO
    RETURN

    !     LAST CARD OF SUBROUTINE COVAR.

  END SUBROUTINE covar


  !LSUNC

  SUBROUTINE lsunc(x,n,mdc,mdg,m,tol,epsrel,epsabs,epsx,sec,  &
       iprint,nout,maxit,scale,ffunc,bl,bu,bnd,  &
       EXIT,phi,k,funcev,jacev,secev,linev,prank,p,f,c,diag,speed, &
       aset,dx,g,w0,w1,w2,w3,v,gmat)

    REAL (dp), INTENT(IN OUT)  :: x(:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: mdc
    INTEGER, INTENT(IN)        :: mdg
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(IN)      :: tol
    REAL (dp), INTENT(IN)      :: epsrel
    REAL (dp), INTENT(IN)      :: epsabs
    REAL (dp), INTENT(IN)      :: epsx
    LOGICAL, INTENT(IN OUT)    :: sec
    INTEGER, INTENT(IN OUT)    :: iprint
    INTEGER, INTENT(IN OUT)    :: nout
    INTEGER, INTENT(IN)        :: maxit
    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN)      :: bl(:)
    REAL (dp), INTENT(IN)      :: bu(:)
    INTEGER, INTENT(IN)        :: bnd
    INTEGER, INTENT(OUT)       :: EXIT
    REAL (dp), INTENT(OUT)     :: phi
    INTEGER, INTENT(OUT)       :: k
    INTEGER, INTENT(OUT)       :: funcev
    INTEGER, INTENT(OUT)       :: jacev
    INTEGER, INTENT(OUT)       :: secev
    INTEGER, INTENT(OUT)       :: linev
    INTEGER, INTENT(IN OUT)    :: prank
    INTEGER, INTENT(IN OUT)    :: p(:)
    REAL (dp), INTENT(OUT)     :: f(:)
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    REAL (dp), INTENT(IN OUT)  :: diag(:)
    REAL (dp), INTENT(OUT)     :: speed
    INTEGER, INTENT(OUT)       :: aset(:)
    REAL (dp), INTENT(IN OUT)  :: dx(:)
    REAL (dp), INTENT(IN OUT)  :: g(:)
    REAL (dp), INTENT(IN OUT)  :: w0(:)
    REAL (dp), INTENT(IN OUT)  :: w1(:)
    REAL (dp), INTENT(IN OUT)  :: w2(:)
    REAL (dp), INTENT(IN OUT)  :: w3(:)
    REAL (dp), INTENT(IN OUT)  :: v(:)
    REAL (dp), INTENT(IN OUT)  :: gmat(:,:)

    EXTERNAL ffunc

    
    !   ****************************************************************
    !   * LSUNC IS DEVELOPED BY PER LINDSTR\M AND PER-]KE WEDIN AT THE *
    !   * INSTITUTE OF INFORMATION PROCESSING UNIVERSITY OF UME],      *
    !   * S-90187 UME], SWEDEN                                         *
    !   ****************************************************************

    !     PURPOSE...
    !     SOLVE THE NONLINEAR LEAST SQUARES PROBLEM

    !             MINIMIZE  PHI(X) = 0.5*II F(X) II**2
    !     WHERE                                T
    !          F(X)= (F (X) F (X)...... F (X) )
    !                  1     2           M
    !     SUBJECT TO
    !                 BL(I)<=X(I)<=BU(I)  ; I=1,2,.....,N
    !     AND X IS N-DIMENSIONAL

    !     ON ENTRY

    !     X()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          A FIRST APPROXIMATION OF THE PARAMETERS (UNKNOWNS)
    !     N    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS X,P,DIAG,
    !          DX,G,W0,W3 AND THE ORDER OF THE SQUARE ARRAY GMAT
    !     MDC  INTEGER SCALAR CONTAINING LEADING DIMENSION OF ARRAY C
    !          (MDC MUST BE >= M)
    !     MDG  INTEGER SCALAR CONTAINING LEADING DIMENSION OF ARRAY GMAT
    !          (MDG MUST BE >= N IF SEC IS TRUE)
    !     M    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS F,W1,V,W2
    !     TOL  REAL SCALAR CONTAINING A SMALL VALUE >=0 USED TO DETERMINE
    !          PSEUDO RANK OF THE JACOBIAN J (SEE PRANK)
    !     EPSREL ARE ALL REAL SCALARS CONTAINING SMALL POSITIVE VALUES
    !     EPSABS USED TO TEST CONVERGENCE
    !     EPSX   EPSREL USED BY CRITERION NO. 1
    !            EPSABS USED BY CRITERION NO. 2
    !            EPSX USED BY CRITERION NO. 3
    !     SEC  LOGICAL SCALAR THAT IS TRUE IF THE USER PERMITS USE OF 2:ND
    !          DERIVATIVES AND THAT IS FALSE OTHERWISE
    !     IPRINT INTEGER SCALAR CONTAINING CODE TO CONTROL PRINTING
    !           IF <=0 THEN NO PRINTING IS DONE INSIDE THIS ROUTINE
    !                  ELSE CERTAIN INFORMATION CONSERNING THE FIRST AND
    !                       EVERY IPRINT ITERATION IS PRINTED ON UNIT NOUT
    !                       IN TABLE FORM WITH A HEADER
    !     NOUT INTEGER SCALAR CONTAINING LOGICAL NUMBER FOR OUTPUT UNIT
    !          WHERE PRINTING CAN BE DONE (DEPENDING ON PARAMETER IPRINT)
    !          RECORD LENGTH MUST NOT BE LESS THAN 120
    !     maxit  INTEGER SCALAR CONTAINING MAXIMUM NO. OF ALLOWED ITERATIONS
    !     SCALE INTEGER SCALAR CONTAINING CODE TO CONTROL INTERNAL SCALING
    !           IF =0 THEN NO SCALING IS DONE
    !                 ELSE THE COLUMNS OF THE JACOBIAN ARE SCALED TO
    !                      HAVE UNIT LENGTH
    !     FFUNC        SUBROUTINE NAME-THAT MUST BE DECLEARED EXTERNAL IN
    !                  THE CALLING PROGRAM. A ROUTINE NAMED FFUNC IS USED
    !                  TO EVALUATE THE FUNCTION F(X) AND/OR THE JACOBIAN AT
    !                  A CERTAIN POINT X AND SHOULD BE WRITTEN AS FOLLOWS

    !                  SUBROUTINE FFUNC(X,N,F,M,CTRL,C,MDC)
    !                  INTEGER N,M,CTRL,MDC
    !                  REAL X(N),F(M),C(MDC,N)
    !                  -----------------------
    !                  CTRL CAN HAVE 3 DIFFERENT VALUES ON ENTRY
    !         CTRL= 1  MEANS EVALUATE THE FUNCTIONS AT THE POINT X AND
    !                  RETURN THIS VECTOR IN THE ARRAY F IF THE FUNCTIONS
    !                  ARE COMPUTABLE.
    !                  ON RETURN THE USER CAN INDICATE UNCOMPUTABILITY BY
    !                  SETTING CTRL=-1
    !                  DO NOT ALTER ARRAY X.
    !         CTRL=-1  MEANS EVALUATE THE FUNCTIONS AT THE POINT X AND
    !                  RETURN THIS VECTOR IN THE ARRAY F IF THE FUNCTIONS
    !                  ARE COMPUTABLE. DO NOT ALTER ARRAY X.
    !                  POSSIBLE UNCOMPUTABILITY OF THE FUNCTIONS MUST BE
    !                  INDICATEDBY SETTING CTRL TO A VALUE <-10 ON RETURN
    !         CTRL= 2  MEANS CALCULATE THE JACOBIAN OF F(X) AT THE POINT X
    !                  AND RETURN THIS MATRIX IN THE ARRAY C IF THE JACOBIAN
    !                  IS SUPPLIED ANALYTICALLY.
    !                  POSSIBLE UNCOMPUTABILITY OF THE JACOBIAN MUST BE
    !                  INDICATEDBY SETTING CTRL TO A VALUE <-10 ON RETURN
    !                  IF THE USER WANTS THE JACOBIAN BEING COMPUTED
    !                  NUMERICALLY THAT SHOULD BE INDICATED BY SETTING
    !                  CTRL=0 ON RETURN.
    !                  DO NOT ALTER ARRAYS X AND F.
    !                  ------------------------------
    !                  RETURN
    !                  END
    !     BL() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          THE LOWER BOUNDS OF THE UNKNOWNS
    !     BU() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          THE UPPER BOUNDS OF THE UNKNOWNS

    !     ON RETURN     AND EXIT >=0

    !     X()  CONTAINS THE LATEST (BEST) ESTIMATE OF THE SOLUTION POINT

    !     THE CONVERGENCE CRITERIA ARE

    !     1) RELATIVE PREDICTED REDUCTION IN THE OBJECTIVE FUNCTION
    !        IS LESS THAN EPSREL**2
    !     2) THE SUM OF SQUARES IS LESS THAN EPSABS**2
    !     3) THE RELATIVE CHANGE IN X IS LESS THAN EPSX
    !     4) WE ARE COMPUTING AT NOISE LEVEL
    !        THE LAST DIGIT IN THE CONVERGENCE CODE (SEE BELOW) INDICATES
    !        HOW THE LAST STEPS WERE COMPUTED
    !        = 0 NO TROUBLE (GAUSS-NEWTON THE LAST 3 STEPS)
    !        = 1 PRANK<>N AT THE TERMINATION POINT
    !        = 2 THE METHOD OF NEWTON WAS USED (AT LEAST) IN THE LAST STEP
    !        = 3 THE 2:ND BUT LAST STEP WAS SUBSPACE MINIMIZATION BUT THE
    !            LAST TWO WERE GAUSS-NEWTON STEPS
    !        = 4 THE STEPLENGTH WAS NOT UNIT IN BOTH THE LAST TWO STEPS


    !      THE ABNORMAL TERMINATION CRITERIA ARE

    !     5) NO. OF ITERATIONS HAS EXCEEDED MAXIMUM ALLOWED ITERATIONS
    !     6) THE HESSIAN EMANATING FROM 2:ND ORDER METHOD IS NOT POS DEF
    !     7) THE ALGORITHM WOULD LIKE TO USE 2:ND DERIVATIVES BUT IS
    !        NOT ALLOWED TO DO THAT
    !     8) AN UNDAMPED STEP WITH NEWTONS METHOD IS A FAILURE
    !     9) THE LATEST SEARCH DIRECTION COMPUTED USING SUBSPACE
    !        MINIMIZATION WAS NOT A DESCENT DIRECTION (PROBABLY CAUSED
    !        BY WRONGLY COMPUTED JACOBIAN)

    !     EXIT INTEGER SCALAR THAT INDICATE WHY THE RETURN IS TAKEN
    !          =10000  CONVERGENCE DUE TO CRITERION NO. 1
    !          = 2000  CONVERGENCE DUE TO CRITERION NO. 2
    !          =  300  CONVERGENCE DUE TO CRITERION NO. 3
    !          =   40  CONVERGENCE DUE TO CRITERION NO. 4
    !          =    X   WHERE X EQUALS 0,1,2,3 OR 4

    !          <0   INDICATES THAT NO CONVERGENCE CRITERION IS FULFILLED
    !               BUT SOME ABNORMAL TERMINATION CRITERION IS SATISFIED
    !          =  -1  IF M<N OR N<=0 OR M<=0 OR MDC<M OR MDG<N OR MAX<=0
    !                     OR SCALE<0 OR TOL<0 OR ANY OF EPSILON VALUES <0
    !                     OR INVALID STARTING POINT   ON ENTRY
    !          =  -2   TERMINATION DUE TO CRITERION NO. 5
    !          =  -3   TERMINATION DUE TO CRITERION NO. 6
    !          =  -4   TERMINATION DUE TO CRITERION NO. 7
    !          =  -5   TERMINATION DUE TO CRITERION NO. 8
    !          =  -6   TERMINATION DUE TO CRITERION NO. 9
    !          =  -7   THERE IS ONLY ONE FEASIBLE POINT, NAMELY
    !                  X(I)=BL(I)=BU(I)  I=1,2,....,N
    !          <  -10  TERMINATION DUE TO USER STOP INDICATOR
    !     PHI  REAL SCALAR CONTAINING THE VALUE OF THE OBJECTIVE FUNCTION
    !          AT THE TERMINATING POINT HELD IN X()
    !     K    INTEGER SCALAR CONTAINING NO. OF ITERATIONS UNTIL
    !          TERMINATION
    !     FUNCEV INTEGER SCALAR CONTAINING THE TOTAL NO. OF FUNCTION
    !            EVALUATIONS
    !     JACEV  INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !           CAUSED BY COMPUTING JACOBIANS WITH DIFFERENCE METHODS
    !     SECEV  INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !           CAUSED BY USING 2:ND DERIVATIVE INFORMATION
    !     LINEV  INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !           CAUSED BY THE LINESEARCH ALGORITHM
    !     PRANK INTEGER SCALAR CONTAINING THE ESTIMATED PSEUDO RANK OF THE
    !           JACOBIAN MATRIX J AT THE TERMINATING POINT
    !     P()  INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          PERMUTATION MATRIX E IN THE DECOMPOSITION   T
    !                                                     Q *J*D*E = (R)
    !                                                                (0)
    !          P(I) CONTAINS THE INDEX OF THE COLUMN THAT WAS MOVED
    !          INTO COLUMN I
    !     F()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTAINING THE
    !          VALUE OF THE FUNCTIONS F (X) AT THE TERMINATING POINT
    !                                  I       I=1,2,.....,M
    !     C(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDC*MAX(4,N)
    !          (MDC >= M)
    !          IF EXIT <0 ON RETURN ARRAY C IS UNDEFINED
    !          OTHERWISE C CONTAINS THE MATRIX R FROM THE
    !          DECOMPOSITION    T
    !                          Q *J*D*E = (R)
    !                                     (0)
    !          IN THE UPPER TRIANGLE OF C
    !           WHERE
    !           Q IS ORTHOGONAL (M*M)
    !           J IS THE JACOBIAN (M*N) AT THE TERMINATING POINT
    !           D IS A DIAGONAL MATRIX (N*N)
    !           E IS A PERMUTATION MATRIX (N*N)
    !           R IS AN UPPER TRIANGULAR MATRIX (N*N)
    !     I.E.           T  -1
    !         J = Q*(R)*E *D
    !               (0)
    !     AND
    !          T      -1    T    T  -1     T   -1        -1   T -1  T
    !         J *J = D  *E*R *R*E *D     (J *J)   = D*E*R  *(R )  *E *D

    !     WHICH IS THE MAINPART OF THE COVARIANCE MATRIX

    !     DIAG() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !           THE DIAGONAL ELEMENTS OF THE DIAGONAL MATRIX D IF SCALE>0
    !           ON ENTRY. OTHERWISE UNDEFINED
    !     SPEED   REAL SCALAR CONTAINING AN ESTIMATE OF THE LINEAR
    !             CONVERGENCE FACTOR
    !     ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !           CONTAINING A CODE WHICH INDICATES WHETHER AN UNKNOWN IS
    !           ACTIVE OR NOT.
    !           ASET(I) = 0 WHEN X(I) IS FREE
    !                   =+1 WHEN X(I) EQUALS BL(I)
    !                   =-1 WHEN X(I) EQUALS BU(I)

    !     WORKING AREAS
    !     DX(),G(),W0(),W3() REAL SINGLY SUBSCRIPTED ARRAYS OF DIMENSION N
    !     W1(),W2(),V()  REAL SINGLY SUBSCRIPTED ARRAYS OF DIMENSION M
    !     GMAT(,)        REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDG*N
    !                    (MDG >= N). IF SEC=FALSE GMAT CAN BE A SINGLY
    !                    SUBSCRIPTED ARRAY OF DIMENSION N

    !     COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !     2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !     REPRESENT TIME STEP K-2 RESP. K-1
    !     THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    !     COMMON VARIABLES USED FOR RESTART INFORMATION

    ! COMMON /back/ icount,lattry,itotal,philat,bestpg,irank
    ! COMMON /INDEX/ imax,ival
    !*************
    ! COMMON /negdir/ ifree,indic
    !*************

    !     THIS PROGRAM PACKAGE USES THE FOLLOWING LINPACK AND BLAS ROUTINES

    !     SQRDC,SQRSL,STRSL,SCHDC,SPOSL
    !     SAXPY,SDOT,SSCAL,SSWAP,SNRM2,SCOPY


    !     INTERNAL VARIABLES

    INTEGER   :: i,kod,error,ctrl,eval,nract,nrcons
    REAL (dp) :: tau,fsum,xnorm,d1sqs,beta,alpha,gnorm,dxnorm,alphup
    REAL (dp) :: alfnoi,prekm1,xdiff
    LOGICAL   :: restar
    REAL (dp), SAVE :: c1 = 0.0_dp

    !     VALIDATE INPUT VALUES

    EXIT=0
    IF(m < n .OR. n <= 0 .OR. m <= 0 .OR. mdc < m ) EXIT=-1
    IF(mdg < n .OR. maxit <= 0 .OR. scale < 0 .OR. tol <= zero) exit=-1
    IF(epsrel < zero .OR. epsabs < zero .OR. epsx < zero) exit=-1
    IF(epsrel+epsabs+epsx <= zero) EXIT=-1

    !     INITIATE ASET(I) AND CHECK INPUT VALUES IN BL(I) AND BU(I)

    nract=0
    nrcons=0
    DO  i=1,n
       IF(bl(i) > bu(i)) EXIT=-1
       aset(i)=0
       IF(x(i) < (bl(i)+c1)) x(i)=bl(i)
       IF(x(i) > (bu(i)-c1)) x(i)=bu(i)
       IF(x(i) == bu(i)) aset(i)=-1
       IF(x(i) == bl(i) .OR. bl(i) == bu(i)) aset(i)=1
       IF(aset(i) /= 0) nract=nract+1
       IF(bl(i) == bu(i)) nrcons=nrcons+1
    END DO
    IF(EXIT < 0) RETURN

    !     COMPUTE SRELPR = SINGLE RELATIVE PRECISION

    CALL releps(srelpr)

    !     INITIATE VARIABLES

    k=0
    !************
    ifree=0
    indic=0
    !************
    error=0
    funcev=0
    jacev=0
    secev=0
    linev=0
    restar=.false.
    icount=0
    itotal=0
    tau=tol
    xdiff=dnrm2(n,x,1)*2.0_dp*epsx
    betkm1=zero
    alfkm1=one
    alpha=zero
    alfkm2=one
    kodkm2=1
    kod=1
    rngkm1=MAX(1,n-nract)
    lattry=rngkm1
    kodkm1=1
    bestpg=zero
    aupkm1=zero

    !     EVALUATE AT USER SUPPLIED POINT

    ctrl=1
 
    CALL ffunc(x,n,f,m,ctrl,c,mdc)
    IF(-1 == ctrl) EXIT=-1
    IF(EXIT < 0) RETURN
    funcev=funcev+1
    phi=0.5_dp*dnrm2(m,f,1)**2
    fsqkm1=2.0_dp*phi
    IF(nrcons == n) EXIT=-7 ! n=# of fitvar

    IF(EXIT < 0) RETURN

    !     MAIN LOOP OF ITERATION STARTS HERE

10  IF(restar.OR.(error < 0)) GO TO 20

    !     COMPUTE GAUSS-NEWTON SEARCH DIRECTION (DX) BY SOLVING THE
    !     LINEARIZED PROBLEM  IF NOT A RESTART STEP

    !*************

15  CALL soliuc(x,n,f,m,ffunc,tau,scale,mdc,bl,bu,aset,nract,c,funcev,jacev, &
         dx,v,p,diag,g,w0,d1sqs,beta,prekm1,gnorm,prank,error)
    IF(error < -10) GO TO 50

    !   ON RETURN BETA  IS THE NORM OF THE ORTHOGONAL PROJECTION OF F
    !                   ONTO THE SPACE SPANNED BY THE COLUMNS OF THE
    !                   JACOBIAN
    !             PREKM1 IS PREDICTED REDUCTION IF PSEUDORANK-1
    !                    FROM PREVIOUS STEP IS USED
    !             D1SQS IS PREDICTED REDUCTION IN THE OBJECTIVE
    !                   FUNCTION IF DX IS USED AS SEARCH DIRECTION
    !             GNORM IS THE NORM OF THE GRADIENT DIVIDED BY THE
    !                   LENGTH OF THE LONGEST COLUMN OF THE JACOBIAN
    !             PRANK IS PSEUDO RANK USED TO COMPUTE THE DIRECTION DX

    !   COMPUTE THE NORM OF CURRENT POINT
    !           THE NORM OF THE SEARCH DIRECTION
    !           THE SUM OF SQUARES AT THE CURRENT POINT
    !           A NOISE LEVEL INDICATOR

    xnorm=dnrm2(n,x,1)
    dxnorm=dnrm2(n,dx,1)
    fsum=dnrm2(m,f,1)**2
    !      ALFNOI=SQRT(SRELPR)/(DXNORM+SRELPR)
    IF(dxnorm < srelpr) THEN
       alfnoi=one
    ELSE
       alfnoi=SQRT(srelpr)/dxnorm
    END IF

    !     CHECK MAIN TERMINATION CRITERIA

20  CALL termuc(prank,error,restar,maxit,k,fsum,d1sqs,n,xnorm,  &
         alfnoi,xdiff,epsabs,epsrel,epsx,g,bl,bu,aset,nract,EXIT)
      
    !     IF EXIT<>0 THE ITERATION IS FINISHED

    IF(EXIT /= 0) GO TO 40

    !     ANALYSE THE PAST AND SOMETIMES RECOMPUTE THE SEARCH DIRECTION

    CALL analuc(k,restar,kod,fsum,d1sqs,beta,dxnorm,prekm1,ffunc,x,  &
         dx,diag,w0,p,n,prank,scale,f,v,c,mdc,m,mdg,sec,nract,aset,  &
         error,eval)

    !      write(10,*) 'POINT',(x(i),i=1,n)
    !      write(10,*) 'PRANK',prank,'DIRECTION',(dx(i),i=1,n)
    IF(error < -10) GO TO 50

    !     ACCUMULATE FUNCTION EVALUATIONS

    funcev=funcev + eval
    secev=secev + eval

    !     CHECK SOME ABNORMAL TERMINATION CRITERIA
    !     ERROR = -3 IF NOT POSITIVE DEFINITE HESSIAN
    !             -4 IF NOT ALLOWED TO USE 2:ND DERIVATIVES

    !      IF(ERROR.LE.(-3)) GOTO 20
    !*********************
    IF(-3 == error)THEN
       !       write(10,*) '********Hessian NOT POS-DEF: Try Gauss-Newton'
       ifree=5
       indic=-3
       GO TO 15
    END IF
    !*********************

    !     SAVE CURRENT POINT

    CALL saveuc(x,n,w3)

    !     COMPUTE STEP LENGTH AND TAKE A STEP

    CALL stepuc(x,g,dx,n,f,v,m,phi,ffunc,kod,prank,bl,bu,aset,nract,bnd,eval,  &
         alpha,alphup,xdiff,error)
    IF(error < -10) GO TO 50

    !     ACCUMULATE FUNCTION EVALUATIONS

    funcev=funcev + eval
    linev=linev + eval

    !     POSSIBLY A RESTART IS DONE
    !     IF NO RESTART IS DONE VARIABLES REPRESENTING THE PAST IS UPDATED

    CALL evreuc(w3,n,m,k,ffunc,funcev,alpha,alphup,d1sqs,beta,fsum,  &
         dxnorm,kod,prank,phi,error,x,f,restar,aset,bl,bu,nract)
    IF(error < -10) GO TO 50
    !***********************
    !      if(ifree.gt.0) goto 15
    !***********************

    !     PRINT SOME INFORMATION DEPENDING ON IPRINT

    CALL outuc(iprint,k,restar,nout,phi,gnorm,eval,aset,n,nract,speed)
    !**********************
    IF(ifree > 0) GO TO 15
    !**********************

    !     REPEAT FROM START OF ITERATION LOOP

    GO TO 10

40  speed=zero
    IF(betkm1 /= zero) speed=beta/betkm1
    RETURN

    !     INDICATE USER STOP AND RETURN

50  EXIT=error

    RETURN
  END SUBROUTINE lsunc


  !SOLIUC

  SUBROUTINE soliuc(x,n,f,m,ffunc,tau,scale,mdc,bl,bu,aset,nract,c,funcev,  &
       jacev,dx,v,p,diag,g,qraux,d1sqs,beta,prekm1,gnorm,  &
       prank,ERR)

    ! N.B. Arguments WORK, W1 & W2 have been removed.

    REAL (dp), INTENT(IN OUT)  :: x(:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: f(:)
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(IN)      :: tau
    INTEGER, INTENT(IN)        :: scale
    INTEGER, INTENT(IN)        :: mdc
    REAL (dp), INTENT(IN)      :: bl(:)
    REAL (dp), INTENT(IN)      :: bu(:)
    INTEGER, INTENT(IN OUT)    :: aset(:)
    INTEGER, INTENT(IN OUT)    :: nract
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN OUT)    :: funcev
    INTEGER, INTENT(IN OUT)    :: jacev
    REAL (dp), INTENT(IN OUT)  :: dx(:)
    REAL (dp), INTENT(IN OUT)  :: v(:)
    INTEGER, INTENT(IN OUT)    :: p(:)
    REAL (dp), INTENT(IN OUT)  :: diag(:)
    REAL (dp), INTENT(OUT)     :: g(:)
    REAL (dp), INTENT(IN OUT)  :: qraux(:)
    REAL (dp), INTENT(IN OUT)  :: d1sqs
    REAL (dp), INTENT(OUT)     :: beta
    REAL (dp), INTENT(OUT)     :: prekm1
    REAL (dp), INTENT(OUT)     :: gnorm
    INTEGER, INTENT(IN OUT)    :: prank
    INTEGER, INTENT(OUT)       :: ERR

    EXTERNAL ffunc

    
    !   COMPUTE THE PSEUDO INVERSE SOLUTION (DX) OF THE SYSTEM
    !                    J*DX = -F
    !   WHERE
    !        J IS THE JACOBIAN OF F(X) AT THE CURRENT POINT
    !        F IS THE VECTOR OF RESIDUALS
    !   ON ENTRY

    !   X()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE CURRENT POINT
    !   N    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS X,DX,P,
    !        DIAG,G AND QRAUX
    !   F()  REAL SINGLY SUBSCRIPTED ARRAY CONTAINING THE VECTOR OF RESIDUALS
    !   M    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS F,V,WORK AND W1
    !   FFUNC SUBROUTINE NAME USED TO EVALUATE THE VECTOR OF RESIDUALS
    !   TAU  REAL SCALAR CONTAINING A SMALL POSITIVE VALUE USED TO
    !        DETERMINE THE PSEUDO RANK OF THE JACOBIAN
    !   SCALE INTEGER SCALAR CONTAINING A CODE THAT CONTROLS INTERNAL
    !        SCALING. =0 IMPLIES NO INTERNAL SCALING
    !                 =1 IMPLIES COLUMN SCALING OF THE JACOBIAN
    !   MDC  INTEGER SCALAR CONTAINING LEADING DIMENSION OF ARRAY C
    !   BL() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE LOWER BOUNDS OF THE UNKNOWNS
    !   BU() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE UPPER BOUNDS OF THE UNKNOWNS
    !   NRACT INTEGER SCALAR CONTAINING THE NO. OF ACTIVE CONSTRAINTS
    !         EQUALS -NO. OF ACTIVE CONSTRAINTS WHEN ONE IS DELETED IN THIS STEP
    !   ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING A CODE
    !         WHICH INDICATES WHETHER AN UNKNOWN IS ACTIVE OR NOT.
    !         ASET(I) = 0 WHEN X(I) IS FREE
    !                 =+1 WHEN X(I) EQUALS BL(I)
    !                 =-1 WHEN X(I) EQUALS BU(I)

    !   ON RETURN

    !   C(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDC*N (MDC>=M)
    !        CONTAINING THE NON ZERO PART OF THE UPPER TRIANGULAR N*N
    !        MATRIX R FROM THE QR-
    !        DECOMPOSITION    T
    !                        Q *J*D*E = (R)
    !                                   (0)
    !        IN THE UPPER TRIANGLE.
    !        INFORMATION NEEDED TO FORM MATRIX Q IS STORED IN THE LOWER
    !        PART OF C AS THE OUTPUT FROM LINPACK ROUTINE SQRDC
    !   FUNCEV INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS DONE SO FAR
    !   JACEV INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !         CAUSED BY COMPUTING JACOBIANS WITH DIFFERENCE METHODS
    !   DX() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE PSEUDO INVERSE SOLUTION THE SYSTEM ABOVE
    !   V()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTAINING
    !        THE ORTHOGONAL PROJECTION OF F ONTO THE SPACE SPANNED BY
    !        THE FIRST PRANK LINEAR INDEPENDANT COLUMNS OF THE JACOBIAN
    !   P()  INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !        CONTAINING THE PERMUTATION MATRIX E (FROM ABOVE) STORED
    !        IN A SPECIAL WAY (SEE OUTPUT FROM LINPACK ROUTINE SQRDC)
    !   DIAG() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !          CONTAINING THE DIAGONAL ELEMENTS OF THE DIAGONAL
    !          MATRIX D  (FROM ABOVE) IF SCALE=1
    !   G()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE GRADIENT OF THE OBJECTIVE AT THE CURRENT POINT
    !   QRAUX() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !           CONTAINING ADDITIONAL INFORMATION NEEDED TO FORM
    !           MATRIX Q (FROM ABOVE) (SEE LINPACK ROUTINE SQRDC)
    !   D1SQS REAL SCALAR CONTAINING THE PREDICTED REDUCTION IN
    !         THE OBJECTIVE FUNCTION IF GN-DIRECTION IS USED
    !   BETA  REAL SCALAR CONTAINING THE NORM OF THE ORTHOGONAL
    !         PROJECTION OF F ONTO THE SPACE SPANNED BY THE COLUMNS
    !         OF THE JACOBIAN
    !   PREKM1 REAL SCALAR CONTAINING THE PREDICTED REDUCTION IN THE
    !         OBJECTIVE IF PSEUDORANK-1 FROM PREVIOUS STEP IS USED
    !   GNORM REAL SCALAR CONTAINING THE NORM OF THE GRADIENT
    !         DIVIDED BY THE NORM OF THE JACOBIAN
    !   PRANK INTEGER SCALAR CONTAINING THE PSEUDO RANK OF MATRIX
    !         J THAT WAS USED TO COMPUTE DX
    !   ERR   INTEGER SCALAR CONTAINING A USER STOP INDICATOR <-10 OR EQUALS ZERO

    !   WORKING AREAS

    !   WORK() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N USED
    !          IN CALLING LINPACK ROUTINE SQRDC
    !   W1() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M USED IF
    !        THE JACOBIAN IS COMPUTED NUMERICALLY
    !   W2() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M

    REAL (dp)  :: work(n), w1(m), w2(m), work2(n,1)

    !   COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !   2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !   REPRESENT TIME STEP K-2 RESP. K-1
    !   THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     COMMON VARIABLES CONSERNING RESTART STEPS

    ! COMMON /INDEX/ imax,ival

    !     INTERNAL VARIABLES

    INTEGER   :: i,j,k,inds,ip(1),kk,ii,lprank,mr,info,kmax
    REAL (dp) :: cnorm,tol,scaqr(1),dummy(1),dykk(1),d1max,d1new,ymax
    REAL (dp), SAVE  :: factor = 2.0_dp

    !      write(10,*) 'Current point in SOLIUC'
    !      write(10,*) (x(i),i=1,n)

    !     COMPUTE JACOBIAN AT THE CURRENT POINT AND STORE IN ARRAY C

    imax=0
    ERR=0
    CALL newunc(x,n,f,m,ffunc,funcev,jacev,c,mdc,ERR,w1)
    ! 101  format(10e8.1)
    IF(ERR /= 0) RETURN

    !     COMPUTE THE GRADIENT OF THE OBJECTIVE AT CURRENT POINT

    CALL grad(c,m,n,f,g)

    !     SCALE THE JACOBIAN IF SCALE INDICATOR SAYS SO
    !       CNORM:=LENGTH OF LONGEST COLUMN OF THE JACOBIAN

    CALL scaunc(scale,c,m,n,cnorm,diag)
    !      write(10,*) '****** NaN test ******'
    gnorm=dnrm2(n,g,1)/cnorm
    !      write(10,*) 'cnorm,gnorm= ', cnorm,gnorm

    !     MAKE A QR-DECOMPOSITION OF (POSSIBLY SCALED) JACOBIAN
    !     I.E.                T
    !                        Q *J*D*E = (R)
    !                                   (0)
    !     BY USING THE LINPACK ROUTINE         SQRDC
    !     ALSO DETERMINE PSEUDO RANK OF MATRIX J*D

    tol=SQRT(REAL(n))*tau
    CALL triunc(c,m,n,tol,aset,nract,p,qraux,prank)

    !     COMPUTE THE "PSEUDO INVERSE" SOLUTION (DX) BY USING PRANK AS
    !     PSEUDO RANK OF MATRIX J*D
    !     ALSO COMPUTE  V= J*DX
    !     AND       D1SQS= PREDICTED REDUCTION IN THE OBJECTIVE
    !               BETA = THE NORM OF ORTH. PROJ. OF F ONTO COLUMN
    !                      SPACE OF THE JACOBIAN J

    CALL gndunc(f,m,c,n,scale,diag,p,qraux,prank,dx,d1sqs,v, work,w1)
    !      write(10,*) 'The right hand side'
    !      write(10,101) (w1(i),i=1,n)
    beta=dnrm2(n-nract,w1,1)
    info=rngkm1-1
    IF(alfkm1 == aupkm1) info=rngkm1-2
    info=MAX(1,info)
    prekm1=beta
    IF(ABS(kodkm1) == 1) prekm1=dnrm2(info,w1,1)

    !     DELETE THE CONSTRAINT THAT CAUSES THE LARGEST PREDICTED
    !     REDUCTION IN THE OBJECTIVE

    !      write(10,*) 'NRACT,PRANK,IMAX=',NRACT,PRANK,IMAX
    !      IF((NRACT.EQ.0).OR.(PRANK.LT.(N-NRACT)).OR.(IMAX.NE.0)) RETURN
    IF((nract == 0).OR.(imax /= 0)) RETURN

    !     INITIATE

    ii=1
    kk=n-nract+1
    mr=n-kk+1
    d1max=d1sqs
    kmax=0
    w1(1:m)=-f(1:m)
    CALL sqrsl(c,m,n,qraux,w1,dummy,w1,w2,dummy,dummy,1000,info)

    !     TRY EACH "BOUNDED" COLUMN TO SEE WHAT PREDICTED REDUCTION IT GIVES

    inds=0
    !      do 35 i=1,n
    !       w2(i)=dx(i)
    !   35 continue
    DO  i=1,n
       IF(aset(i) == 0) CYCLE
       IF(prank < n-nract)THEN
          k=1
          kk=1
          !       write(10,*) 'Not full rank: prank= ',prank
          inds=aset(i)
          aset(i)=0
          CALL newunc(x,n,f,m,ffunc,funcev,jacev,c,mdc,ERR,w1)
          !tpk IF(ERR /= 0) STOP
          IF(ERR /= 0) RETURN
          CALL triunc(c,m,n,tol,aset,nract-1,p,qraux,lprank)
          !       write(10,*) 'lprank= ', lprank
          CALL gndunc(f,m,c,n,0,diag,p,qraux,lprank,dx,d1new,v,work,w1)
          !       write(10,*) 'New search direction'
          !       write(10,*) (dx(jj),jj=1,n)
          aset(i)=inds
          IF(dx(i)*inds*(bu(i)-bl(i)) <= zero) CYCLE
       ELSE
          k=n-nract+ii
          ii=ii+1
          !                       T
          !     PUT W2 EQUAL TO -Q *F

          w2(1:m)=w1(1:m)

          !     ZERO THE K.TH COLUMN BELOW ITS (N-NRACT)+1 ELEMENT

          DO  j=1,n
             work(j)=c(j,k)
             IF(j > k) work(j)=zero
          END DO
          IF(k /= kk) THEN
             ip(1)=1
             work2(1:mr,1) = work(kk:kk+mr-1)
             CALL sqrdc(work2,mr,1,scaqr,ip,0)

             !     TANSFORM W2 IN THE SAME WAY

             CALL sqrsl(work2,mr,1,scaqr,w2(kk:),dummy,w2(kk:),  &
                  dykk,dummy,dummy,100,info)
             work(kk:kk+mr-1) = work2(1:mr,1)
          END IF

          !     DETERMINE THE KK:TH ELEMENT IN DY

          !           write(10,*) 'SOLIUC:WORK(KK),TOL,C(1,1)=',WORK(KK),TOL,C(1,1)
          IF(ABS(work(kk)) <= tol*ABS(c(1,1))) CYCLE
          dykk(1)=w2(kk)/work(kk)
          !    write(10,*) 'In SOLIUC: DYKK=', DYKK
          !    write(10,*) 'ASET(I),BU(I),BL(I)=',ASET(I),BU(I),BL(I)
          IF(dykk(1)*aset(i)*(bu(i)-bl(i)) <= zero) CYCLE
          d1new=d1sqs + w2(kk)*w2(kk)
          !           IF(D1SQS.GT.FACTOR*D1NEW) GOTO 70
       END IF
       !           write(10,*) 'In SOLIUC: D1NEW,D1SQS=',D1NEW,D1SQS
       IF(d1new < factor*d1sqs) CYCLE
       !    write(10,*) 'Delete a fixed variable'
       IF(d1new < d1max) CYCLE
       d1max=d1new
       ymax=w2(kk)
       kmax=k
       imax=i
    END DO
    !      IF(KMAX.EQ.0) RETURN
    !      if(kmax.eq.0) then
    !       if(inds.ne.0) then
    ! do 75 i=1,n
    !dx(i)=w2(i)
    !   75   continue
    !       else
    ! return
    !       endif
    !      endif
    IF(prank < (n-nract))THEN
       IF(imax /= 0)THEN
          ival=aset(imax)
          aset(imax)=0
          nract=nract-1
          !       imax=0
       END IF
       CALL newunc(x,n,f,m,ffunc,funcev,jacev,c,mdc,ERR,w1)
       CALL triunc(c,m,n,tol,aset,nract,p,qraux,prank)
       CALL gndunc(f,m,c,n,0,diag,p,qraux,prank,dx,d1sqs,v,work,w1)
       beta=dnrm2(n-nract,w1,1)
       prekm1=beta
    ELSE IF(kmax /= 0)THEN

       !     FORM J*DX WHERE DX IS THE AUGMENTED GN-DIRECTION.
       !         J = Q*H*(D1)
       !                 ( 0)
       !                                T
       !     FIRST COMPUTE W2:= (D1)= -Q *F
       !                        (D2)

       w2(1:m)=w1(1:m)
       IF(kmax == kk) GO TO 100
       DO  j=1,n
          work(j)=c(j,kmax)
          IF(j > kmax) work(j)=zero
       END DO
       work2(1:mr,1) = work(kk:kk+mr-1)
       CALL sqrdc(work2,mr,1,scaqr,ip,0)
       w2(kk)=ymax

       !     W2:= H*(D1)
       !            ( 0)

100    ip(1)=kk+1
       w2(ip(1):m)=zero
       IF(kmax /= kk) CALL sqrsl(work2,mr,1,scaqr,w2(kk:),w2(kk:),dummy,  &
            dummy,dummy,dummy,10000,info)
       work(kk:kk+mr-1) = work2(1:mr,1)

       !     V:= Q*H*(D1)
       !             ( 0)

       CALL sqrsl(c,m,n,qraux,w2,v,dummy,dummy,dummy,dummy,10000,info)

       !     MOVE COLUMN KMAX TO COLUMN KK

       IF(kmax == kk) GO TO 140
       c(1:kk,kk)=work(1:kk)

       !     CHANGE IN PERMUTATION VECTOR

140    ip=p(kk)
       p(kk)=p(kmax)
       p(kmax)=ip(1)

       !     UPDATE ACTIVE SET VARIABLES

       ival=aset(imax)
       aset(imax)=0
       !      imax=0
       prank=prank+1
       nract=nract-1
       beta=dnrm2(n-nract,w1,1)
       prekm1=beta
       d1sqs=d1max

       !     COMPUTE THE AUGMENTED GN-DIRECTION

       dx(1:kk)=w1(1:kk)
       dx(kk)=ymax
       CALL strsl(c,prank,dx,1,info)

       !     BACK TRANSFORM

       CALL btrunc(dx,n,prank,diag,p,scale,w2)
       !      write(10,*) 'Full rank new search direction'
       !      write(10,*) (dx(i),i=1,n)
    END IF

    RETURN
  END SUBROUTINE soliuc


  !NEWUNC

  SUBROUTINE newunc(x,n,f,m,ffunc,funcev,jacev,c,mdc,ERR,w1)

    REAL (dp), INTENT(IN OUT)  :: x(:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: f(:)
    INTEGER, INTENT(IN)        :: m
    !tpk INTEGER, INTENT(OUT)       :: funcev
    INTEGER, INTENT(INOUT)       :: funcev  !tpk
    !tpk INTEGER, INTENT(OUT)       :: jacev
    INTEGER, INTENT(INOUT)       :: jacev  !tpk
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: mdc
    INTEGER, INTENT(OUT)       :: ERR
    REAL (dp), INTENT(IN OUT)  :: w1(:)

    EXTERNAL ffunc

    
    !   COMPUTE THE JACOBIAN OF F(X) AT THE CURRENT POINT AND STORE IN ARRAY C

    !   ON ENTRY

    !   X()  CONTAINS THE CURRENT POINT
    !   F()  CONTAINS THE VECTOR OF RESIDUALS AT THE CURRENT POINT
    !   N,M  THE LENGTH OF X AND F RESPECTIVELY
    !   MDC  THE LEADING DIMENSION OF ARRAY C
    !   ERR  EQUALS ZERO

    !   ON RETURN

    !   C(,) CONTAINS THE M*N JACOBIAN MATRIX
    !   FUNCEV CONTAINS THE TOTAL NO.OF FUNCTION EVALUATIONS DONE SO FAR
    !   JACEV  CONTAINS NO. OF FUNCTION EVALUATION CAUSED BY COMPUTING
    !          THE JACOBIAN WITH DIFFERENCE METHODS
    !   ERR    CONTAINS A USER STOP INDICATION <-10 OR UNCHANGED

    !   INTERNAL VARIABLES

    INTEGER :: ctrl

    ctrl=2

    !     IF CTRL=2 ON ENTRY TO FFUNC THE JACOBIAN MATRIX IS REQUESTED.
    !     HOWEVER, IF THE ANALYTICAL JACOBIAN IS NOT AVAILABLE THE USER
    !     SIGNALS THAT BY SETTING CTRL TO ZERO ON RETURN
    
    CALL ffunc(x,n,f,m,ctrl,c,mdc)
    IF(ctrl == 2) GO TO 10
    IF(ctrl < -10) GO TO 20

    !     COMPUTE THE JACOBIAN USING FORWARD DIFFERENCES
   
    IF(ctrl == 0) CALL jacdif(x,n,f,m,ffunc,c,mdc,ctrl,w1)
    IF(ctrl < -10) GO TO 20
    funcev=funcev+n
    jacev=jacev+n
    
10  RETURN

    !     USER STOOP INDICATION DETECTED

20  ERR=ctrl

    RETURN
  END SUBROUTINE newunc


  !SCAUNC

  SUBROUTINE scaunc(scale,c,m,n,cnorm,diag)

    ! Argument MDC has been removed.

    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(OUT)     :: cnorm
    REAL (dp), INTENT(OUT)     :: diag(:)

    !   IF SCALE=0
    !     THEN NO SCALING IS DONE
    !     ELSE SCALE THE M*N MATRIX C SO THAT EACH COLUMN HAS UNIT LENGTH

    !   ON RETURN

    !   CNORM SET TO THE MAXIMUM COLUMN LENGTH OF THE MATRIX C
    !   C(,)  IF SCALING IS DONE  C:=C*D WHERE D IS A DIAGONAL MATRIX
    !         WITH DIAGONAL ELEMENTS D(I)=1/LENGTH(I) WHERE LENGTH(I)
    !         IS THE LENGTH OF COLUMN NO. I UNLESS LENGTH(I)=0
    !         WHEN D(I) IS SET TO 1.0
    !   DIAG() IF SCALING IS DONE DIAG(I) HOLDS THE DIAGONAL ELEMENTS
    !          OF MATRIX D ABOVE  (I=1,2,......,N)


    !   INTERNAL VARIABLES

    INTEGER   :: j
    REAL (dp) :: colj

    cnorm=zero
    DO  j=1,n
       colj=dnrm2(m,c(:,j),1)
       IF(colj > cnorm) cnorm=colj
       IF(scale == 0) CYCLE
       IF(colj == zero) colj=one
       c(1:m,j)=c(1:m,j)/colj
       diag(j)=one/colj
    END DO
    IF(cnorm == zero) cnorm=one

    RETURN
  END SUBROUTINE scaunc


  !TRIUNC

  SUBROUTINE triunc(c,m,n,tol,aset,nract,p,qraux,prank)

    ! Arguments MDC & WORK have been removed.

    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN)      :: tol
    INTEGER, INTENT(IN)        :: aset(:)
    INTEGER, INTENT(IN)        :: nract
    INTEGER, INTENT(OUT)       :: p(:)
    REAL (dp), INTENT(IN OUT)  :: qraux(:)
    INTEGER, INTENT(OUT)       :: prank

    !   MAKE A QR-DECOMPOSITION OF THE M*N MATRIX C BY USING THE
    !   ROUTINE SQRDC FROM LINPACK
    !   I.E. DETEMINE MATRICES Q,E,R SO THAT      T
    !                                            Q *C*E = (R)
    !                                                     (0)
    !   WHERE  Q IS M*M ORTHOGONAL
    !          E IS N*N PERMUTATION MATRIX
    !          R IS N*N UPPER TRIANGULAR

    !   ON ENTRY

    !   F()  contains the right hand side in C*dx = f
    !   C(,) CONTAINS THE M*N MATRIX TO DECOMPOSE
    !   MDC  THE LEADING DIMENSION OF ARRAY C
    !   M    NO. OF ROWS IN MATRIX C
    !   N    NO. OF COLUMNS IN MATRIX C
    !   TOL  A SMALL POSITIVE VALUE USED TO DETERMINE PSEUDO RANK OF
    !        MATRIX C  (SEE PRANK)
    !   NRACT INTEGER SCALAR CONTAINING THE NO. OF ACTIVE CONSTRAINTS
    !   ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING A CODE
    !         WHICH INDICATES WHETHER AN UNKNOWN IS ACTIVE OR NOT.
    !         ASET(I) = 0 WHEN X(I) IS FREE
    !                 =+1 WHEN X(I) EQUALS BL(I)
    !                 =-1 WHEN X(I) EQUALS BU(I)

    !   ON RETURN

    !   THE ARRAYS C,P AND QRAUX CONTAIN THE CORRESPONDING OUTPUT FROM
    !   LINPACK ROUTINE SQRDC
    !   ABS(R(1,1)) >= ABS(R(2,2)) >=.......>=ABS(R(N,N))
    !   PRANK IS determined as the "pseudo rank" that is most suitable

    !   INTERNAL VARIABLES

    INTEGER   :: i,j,k,nn
    REAL (dp) :: r11

    prank=n
    IF(n == 0) RETURN

    !     INITIATE PIVOT VECTOR SO THAT ALL COLUMNS CORRESPONDING TO
    !     ACTIVE BOUNDS ARE CONSIDERED AS FINAL COLUMNS

    DO  i=1,n
       p(i)=0
       IF(aset(i) /= 0) p(i)=-1
    END DO

    !     DECOMPOSE MATRIX C

    CALL sqrdc(c,m,n,qraux,p,1)

    !     DETERMINE PSEUDO RANK

    k=0
    r11=ABS(c(1,1))
    nn=n-nract
    IF(nn == 0) GO TO 30
    DO  i=1,nn
       IF(ABS(c(i,i)) >= tol) k=i
    END DO
    IF(k == 0) THEN
       DO  j=1,nn
          IF(ABS(c(j,j)) <= tol*r11) EXIT
          k=j
       END DO
    END IF

30  prank=k

    RETURN
  END SUBROUTINE triunc


  !GNDUNC

  SUBROUTINE gndunc(f,m,c,n,scale,diag,p,qraux,prank,dx,d1sqs,v,work,y)

    ! N.B. Argument MDC has been removed.

    REAL (dp), INTENT(IN)      :: f(:)
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN OUT)  :: diag(:)
    INTEGER, INTENT(IN OUT)    :: p(:)
    REAL (dp), INTENT(IN OUT)  :: qraux(:)
    !tpk INTEGER, INTENT(OUT)       :: prank
    INTEGER, INTENT(INOUT)       :: prank  !tpk
    REAL (dp), INTENT(OUT)     :: dx(:)
    REAL (dp), INTENT(OUT)     :: d1sqs
    REAL (dp), INTENT(IN OUT)  :: v(:)
    REAL (dp), INTENT(IN OUT)  :: work(:)
    REAL (dp), INTENT(OUT)     :: y(:)

    !   COMPUTE THE SOLUTION (DX) OF THE SYSTEM
    !                           T
    !              R*D*E*DX = -Q *F
    !   WHERE
    !         R IS N*N UPPER TRIANGULAR
    !         D IS N*N DIAGONAL MATRIX
    !         E IS N*N PERMUTATION MATRIX
    !         Q IS M*M ORTHOGONAL MATRIX
    !         F IS AN M-VECTOR
    !   THE SOLUTION WILL BE COMPUTED BY ONLY USING THE FIRST PRANK <= N
    !   COLUMNS OF R

    !   ON ENTRY

    !   F()  CONTAINS THE VECTOR OF RESIDUALS AT CURRENT POINT
    !   M    IS THE LENGTH OF THE ARRAYS F,V AND Y
    !   C(,) CONTAINS THE MATRIX R IN THE UPPER TRIANGLE AND
    !        INFORMATION NEEDED TO FORM MATRIX Q IN THE LOWER PART
    !   MDC  THE LEADING DIMENSION OF ARRAY C
    !   N    IS THE LENGHT OF THE ARRAYS DIAG,P,QRAUX AND DX
    !   SCALE =0 IF DIAGONAL MATRIX D IS THE UNIT MATRIX
    !         =1 IF SCALING IS DONE
    !   DIAG() CONTAINS THE DIAGONAL ELEMENTS OF MATRIX D IF SCALE=1
    !          IF SCALE=0 THEN DIAG(K) K=1,2,....,N IS UNDEFINED
    !   P(),QRAUX() CONTAIN INFORMATION RETURNED FROM THE LINPACK
    !               ROUTINE SQRDC AS JPVT AND QRAUX
    !   PRANK IS THE SUGGESTED PSEUDO RANK OF MATRIX R

    !   ON RETURN

    !   DX() THE PSEUDO INVERSE SOLUTION (GAUSS-NEWTON SEARCH DIRECTION)
    !   D1SQS IS THE PREDICTED REDUCTION IN THE OBJECTIVE
    !   V()   IS THE OTHOGONAL PROJECTION OF F ONTO THE SPACE SPANNED
    !        BY THE COLUMNS OF THE JACOBIAN J FOR WHICH
    !                T
    !               Q *J*D*E = (R)
    !                          (0)           HOLDS.

    !   Y()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M  !!!!!
    !        CONTAINS THE FIRST N ELEMENTS OF    T
    !                                          -Q *F

    !   INTERNAL VARIABLES

    INTEGER   :: i,info
    REAL (dp) :: dummy(1),temp

    d1sqs=zero
    y(1:m)=-f(1:m)
    IF(n == 0 .OR. prank <= 0) RETURN
    !                                                         T
    !     COMPUTE THE SOLUTION DX,THE PROJECTION  -V  AND Y=-Q *F

    CALL sqrsl(c,m,prank,qraux,y,dummy,y,dx,dummy,v,1101,info)
    d1sqs=dnrm2(prank,y,1)**2
    IF(info == 0) GO TO 40

    !     SQRSL HAS DETECTED EXACT SINGULARITY OF MATRIX R
    !     INFO = INDEX OF THE FIRST ZERO DIAGONAL ELEMENT OF R
    !     SOLVE UPPER TRIANGULAR SYSTEM

    prank=info-1
    dx(1:n)=y(1:n)
    IF(prank <= 0) GO TO 40
    CALL strsl(c,prank,y,1,info)

    !     MOVE SOLUTION OF TRIANGULAR SYSTEM TO DX

    DO  i=1,n
       temp=dx(i)
       dx(i)=y(i)
       y(i)=temp
    END DO

    !     DO BACKTRANSFORMATIONS   DX:=D*E*DX

40  CALL btrunc(dx,n,prank,diag,p,scale,work)

    RETURN
  END SUBROUTINE gndunc


  !BTRUNC

  SUBROUTINE btrunc(dx,n,prank,diag,p,scale,work)

    REAL (dp), INTENT(OUT)     :: dx(:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: prank
    REAL (dp), INTENT(IN)      :: diag(:)
    INTEGER, INTENT(IN OUT)    :: p(:)
    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN OUT)  :: work(:)

    !   BACKTRANSFORM SO THAT   DY:=D*E*DX
    !   WHERE DX IS AN N-VECTOR
    !         E IS A PERMUTATION MATRIX N*N
    !         D IS A DIAGONAL MATRIX N*N

    !   ON ENTRY

    !   DX() CONTAINS THE VECTOR TO BE TRANSFORMED IN THE ELEMENTS
    !        DX(I) I=1,2...,PRANK
    !        IF PRANK<N  THEN DX(PRANK+J):=0 J=1,2,...(N-PRANK) BEFORE
    !                         TRANSFORMATION
    !   DIAG() CONTAINS THE DIAGONAL ELEMENTS OF DIAGONAL MATRIX D
    !          IF SCALE=1. OTHERWISE DIAG() IS UNDEFINED
    !   P()  CONTAINS THE PERMUTATION MATRIX E STORED IN A SPECIAL
    !        WAY (SEE LINPACK ROUTINE SQRDC)
    !   SCALE =0 IF MATRIX D IS THE UNIT MATRIX
    !         =1 IF SCALING WAS DONE BY MATRIX D

    !   ON RETURN

    !   DX() CONTAINS THE TRANSFORMED VECTOR DY ABOVE

    !   INTERNAL VARIABLES

    INTEGER :: is

    IF(n == 0) RETURN
    IF(prank /= n) THEN
       is=prank+1
       dx(is:n)=zero
    END IF

    !     DO THE PIVOTING

    CALL pivec(p,n,dx,work)

    !     EVENTUALLY RESCALE

    IF(scale == 0) GO TO 40
    dx(1:n)=diag(1:n)*dx(1:n)

40  RETURN
  END SUBROUTINE btrunc


  !JACDIF

  SUBROUTINE jacdif(x,n,f,m,ffunc,c,mdc,ERR,w1)

    REAL (dp), INTENT(IN OUT)  :: x(:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN)      :: f(:)
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(OUT)     :: c(:,:)
    INTEGER, INTENT(IN)        :: mdc
    INTEGER, INTENT(OUT)       :: ERR
    REAL (dp), INTENT(IN OUT)  :: w1(:)

    !   COMPUTE THE M*N JACOBIAN OF F(X) AT THE CURRENT POINT BY
    !   USING FORWARD DIFFERENCES

    !   ON ENTRY

    !   X()  CONTAINS THE CURRENT POINT
    !   N    IS THE LENGTH OF THE ARRAY X
    !   F()  CONTAINS THE VECTOR OF RESIDUALS AT CURRENT POINT
    !   M    IS THE LENGTH OF THE ARRAYS F AND W1
    !   MDC  IS THE LEADING DIMENSION OF THE ARRAY C

    !   ON RETURN

    !   C(,) CONTAINS THE APPROXIMATE JACOBIAN IN THE M*N UPPER PART
    !   ERR  CONTAINS A USER STOP INDICATOR OR ZERO

    !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    EXTERNAL ffunc

    
    !     INTERNAL VARIABLES

    INTEGER   :: i,j,ctrl
    REAL (dp) :: delta,xtemp,deltaj

    ERR=0
    delta=SQRT(srelpr)
     
    DO  j=1,n
       xtemp=x(j) ! for most cases, x< one except for sin
       deltaj=MAX(ABS(xtemp),one)*delta 
       x(j)=xtemp + deltaj
       ctrl=-1
       
       CALL ffunc(x,n,w1,m,ctrl,c,mdc)
       
       IF(ctrl < -10) GO TO 30
       DO  i=1,m
          c(i,j)=(w1(i)-f(i))/deltaj        
       END DO
       x(j)=xtemp
    END DO
      
    RETURN

    !     USER STOP INDICATION DETECTED

30  ERR=ctrl

    RETURN
  END SUBROUTINE jacdif


  !PIVEC

  SUBROUTINE pivec(p,n,v,work)

    INTEGER, INTENT(IN)        :: p(:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: v(:)
    REAL (dp), INTENT(OUT)     :: work(:)

    !     COMPUTE W:=P*V  WHERE P REPRESENTS AN N*N PERMUTATION MATRIX
    !     AND V IS AN N-VECTOR
    !     P(J) CONTAINS THE INDEX OF THE ELEMENT THAT WAS MOVED INTO POSITION J

    !     ON RETURN

    !     V()  CONTAINS THE TRANSFORMED VECTOR W ABOVE

    !     INTERNAL VARIABLES

    INTEGER :: i,k

    IF(n == 0) RETURN
    DO  i=1,n
       k=p(i)
       work(k)=v(i)
    END DO
    v(1:n)=work(1:n)

    RETURN
  END SUBROUTINE pivec


  !ANALUC

  SUBROUTINE analuc(k,restar,kod,fsum,d1sqs,beta,dxnorm,prekm1,  &
       ffunc,x,dx,diag,qraux,p,n,prank,scale,f,v,c,mdc,m,mdg,sec,  &
       nract,aset,error,eval)

    ! N.B. Arguments W1, W3 & GMAT have been removed.

    INTEGER, INTENT(IN OUT)    :: k
    LOGICAL, INTENT(OUT)       :: restar
    INTEGER, INTENT(OUT)       :: kod
    REAL (dp), INTENT(IN OUT)  :: fsum
    REAL (dp), INTENT(IN OUT)  :: d1sqs
    REAL (dp), INTENT(IN OUT)  :: beta
    REAL (dp), INTENT(OUT)     :: dxnorm
    REAL (dp), INTENT(IN OUT)  :: prekm1
    REAL (dp), INTENT(IN OUT)  :: x(:)
    REAL (dp), INTENT(IN OUT)  :: dx(:)
    REAL (dp), INTENT(IN OUT)  :: diag(:)
    REAL (dp), INTENT(IN OUT)  :: qraux(:)
    INTEGER, INTENT(IN OUT)    :: p(:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN OUT)    :: prank
    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN OUT)  :: f(:)
    REAL (dp), INTENT(IN OUT)  :: v(:)
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: mdc
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: mdg
    LOGICAL, INTENT(IN)        :: sec
    INTEGER, INTENT(IN OUT)    :: nract
    INTEGER, INTENT(IN OUT)    :: aset(:)
    INTEGER, INTENT(OUT)       :: error
    INTEGER, INTENT(OUT)       :: eval

    EXTERNAL ffunc

    
    !   CHECK IF PREVIOUS STEP WAS SUFFICIENTLY GOOD AND DECIDE IF THE SEARCH
    !   DIRECTION (DX=GAUSS-NEWTON DIRECTION) FOR CURRENT STEP SHALL BE RECOMPUTED

    !   ON ENTRY

    !   K       INTEGER SCALAR CONTAINING ITERATION COUNT
    !   RESTAR  LOGICAL SCALAR = FALSE IF THIS IS NOT A RESTART STEP
    !                          = TRUE IF THIS IS A RESTART STEP
    !   KOD     INTEGER SCALAR CONTAINING A CODE THAT SAYS HOW THE
    !           PREVIOUS SEARCH DIRECTION WAS COMPUTED
    !           = 1 IF GAUSS-NEWTON DIRECTION
    !           =-1 IF SUBSPACE DIRECTION
    !           =-2 IF NEWTON DIRECTION
    !   FSUM    REAL SCALAR CONTAINING SUM OF SQUARES AT CURRENT POINT
    !   D1SQS   REAL SCALAR CONTAINING PREDICTED REDUCTION IN OBJECTIVE
    !           FUNCTION IF GAUSS-NEWTON DIRECTION IS USED
    !   BETA    REAL SCALAR CONTAINING THE NORM OF THE ORTHOGONAL
    !           PROJECTION OF F ONTO THE COLUMN SPACE OF THE JACOBIAN J
    !   PREKM1  REAL SCALAR CONTAINING
    !           THE SAME AS BETA IF ABS(KOD)=2
    !           IF ABS(KOD)=1 IT CONTAINS THE PROJECTION MENTIONED
    !           IN THE EXPLANATION OF BETA BUT NOT NECESSARYLY ONTO THE
    !           FULL COLUMN SPACE OF THE JACOBIAN. THE DIMENSION OF THE
    !           SPACE IS EQUAL TO 1 LESS THE DIMENSION OF THE SPACE WHERE
    !           MINIMIZATION WAS DONE IN THE PREVIOUS STEP
    !   FFUNC   SUBROUTINE NAME USED TO EVALUATE THE VECTOR OF RESIDUALS
    !           AT A CERTAIN POINT
    !   X()     REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !           THE CURRENT POINT
    !   DX()    REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !           CURRENT GAUSS-NEWTON SEARCH DIRECTION
    !   DIAG()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING THE
    !           DIAGONAL ELEMENTS OF THE DIAGONAL MATRIX D FROM (1) BELOW
    !   QRAUX() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !           INFORMATION NEEDED TO FORM MATRIX Q IN THE QR-
    !           DECOMPOSITION OF THE JACOBIAN OF F(X)
    !   P()     INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !           THE PERMUTATION MATRIX E FROM (1) BELOW
    !   FOR MORE EXACT EXPLANATION OF THE ARRAYS QRAUX AND P SEE THE
    !   EXPLANATION OF ARRAYS QRAUX AND JPVT ON RETURN FROM THE
    !   LINPACK ROUTINE SQRDC WHICH HAS BEEN USED TO MAKE THE
    !   QR-DECOMPOSITION OF MATRIX J*D

    !   N       INTEGER SCALAR CONTAINING LENGTH OF ARRAYS X,DX,DIAG, DIAG AND P
    !   PRANK   INTEGER SCALAR CONTAINING THE PSEUDO RANK USED WHEN
    !           CURRENT GAUSS-NEWTON SEARCH DIRECTION WAS COMPUTED
    !   SCALE   INTEGER SCALAR CONTAINING
    !           = 0 IF SCALING OF JACOBIAN J WAS DONE BEFORE QR-DECOMP.
    !           = 1 IF SCALING WAS DONE
    !   F()     REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTAINING
    !           THE VECTOR OF RESIDUALS AT CURRENT POINT
    !   V()     REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTAINING
    !           THE ORTHOGONAL PROJECTION OF F ONTO THE SPACE SPANNED BY
    !           THE COLUMNS OF THE JACOBIAN J
    !   C(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDC*N  (MDC>=M)
    !        CONTAINING THE MATRIX R FROM THE
    !        DECOMPOSITION    T
    !                        Q *J*D*E = (R)
    !                                   (0)           (1)
    !        IN THE UPPER TRIANGLE OF C
    !         WHERE
    !         Q IS ORTHOGONAL (M*M)
    !         J IS THE JACOBIAN (M*N) AT THE TERMINATING POINT
    !         D IS A DIAGONAL MATRIX (N*N)
    !         E IS A PERMUTATION MATRIX (N*N)
    !         R IS AN UPPER TRIANGULAR MATRIX (N*N)
    !   MDC     INTEGER SCALAR CONTAINING LEADING DIMENSION OF ARRAY C
    !           MDC MUST BE >= M
    !   M       INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS F,V,W1
    !   MDG     INTEGER SCALAR CONTAINING LEADING DIMENSION OF ARRAY GMAT.
    !           MDG MUST BE >= N
    !   SEC    LOGICAL SCALAR = TRUE IF THE USER HAS ALLOWED USE OF
    !          2:ND DERIVATES. FALSE OTHERWISE
    !   NRACT INTEGER SCALAR CONTAINING THE NO. OF ACTIVE CONSTRAINTS
    !          (POSITIVE OR NEGATIVE)
    !   ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !           CONTAINING +1 OR -1 TO INDICATE AN ACTIVE BOUND (OTHERWISE ZERO)

    !   ON RETURN

    !   RESTAR  IS SET TO FALSE
    !   KOD     CONTAINS A CODE THAT INDICATES HOW THE SEARCH DIRECTION
    !           CONTAINED IN ARRAY DX ON RETURN IS COMPUTED
    !           = 1 IF THE DIRECTION IN DX ON ENTRY IS ACCEPTED
    !           =-1 IF MINIMIZATION IN SUBSPACE HAS BEEN USED
    !           TO COMPUTE DX
    !           =-2 IF THE METHOD OF NEWTON HAS BEEN USED TO COMPUTE DX
    !   NRACT IS PUT EQUAL TO -NRACT-1 IF NRACT<0 ON ENTRY
    !   ERROR   INTEGER SCALAR CONTAINING -3 IF COEFFICIENT MATRIX IN
    !           SYSTEM ARISING FROM SECOND DERIVATIVE METHOD IS NOT
    !           POSITIVE DEFINITE
    !           = -4 IF THE ALGORITHM WOULD LIKE TO USE SECOND DERIVATIVES
    !                BUT IS NOT ALLOWED TO DO THAT
    !           < -10 AS A USER STOP INDICATOR
    !           OTHERWISE ERROR = 0 ON RETURN
    !   EVAL    INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !           DONE INSIDE THIS ROUTINE

    !   IF KOD=1 ON RETURN ALL THE OTHER INPUT PARAMETERS ARE UNCHANGED
    !   IF KOD=-1 OR ABS(KOD)=2 RETURN THE FOLLWOING PARAMATERS ARE
    !   CHANGED ON RETURN

    !   DX()    CONTAINS THE SUGGESTED SEARCH DIRECTION
    !   D1SQS   CONTAINS PREDICTED REDUCTION IN OBJECTIVE IF SUGGESTED
    !           SEARCH DIRECTION IS USED
    !   DXNORM  CONTAINS THE NORM OF THE SUGGESTED SEARCH DIRECTION
    !   PRANK   CONTAINS THE DIMENSION OF THE SUBSPACE USED TO COMPUTE
    !           SEARCH DIRECTION OR (IF ABS(KOD)=2 ON RETURN) PRANK=-N
    !   V()     CONTAINS THE ORTHOGONAL PROJECTION OF F ONTO THE SPACE
    !           SPANNED BY THE PRANK FIRST LINEAR INDEPENDENT COLUMNS
    !           OF THE JACOBIAN J  (IF KOD=-1 ON RETURN)
    !           UNCHANGED IF KOD=2 ON RETURN
    !           CONTAINS THE VECTOR OF ESIDUALS IF KOD=-2 ON RETURN
    !   C(,)    UNCHANGED IF KOD=-1 ON RETURN
    !           COMPLETELY DESTROYED IF ABS(KOD)=2 ON RETURN
    !   QRAUX() UNCHANGED IF KOD=-1 ON RETURN
    !           COMPLETELY DESTROYED IF ABS(KOD)=2 ON RETURN

    !   INTERNAL VARIABLES

    INTEGER :: rank,secind
    LOGICAL :: gndok

    eval=0

    !   SET GNDOK=TRUE IF SEARCH DIRECTION CONTAINED IN DX ON ENTRY IS ACCEPTED

    !   SET SECIND=
    !              -2 IF NEWTON METHOD
    !              -1 IF SUBSPACE MINIMIZATION IS THE ALTERNATIVE
    !   FOR RECOMPUTATION OF DX

    CALL gnavuc(k,restar,kod,secind,fsum,beta,prekm1,nract,aset,gndok)
    IF(.NOT.gndok) GO TO 10

    !     SEARCH DIRECTION CONTAINED IN DX ON ENTRY IS ACCEPTED

    kod=1
    error=0
    GO TO 50

    !     RECOMPUTE SEARCH DIRECTION

10  IF(ABS(secind) == 2) GO TO 30

    !     USE MINIMIZATION IN SUBSPACE TO RECOMPUTE

    CALL subuc(restar,kod,fsum,d1sqs,dxnorm,f,m,c,mdc,n,scale,  &
         diag,p,qraux,prank,dx,v,rank,k)
    kod=-1
    IF(rank == prank) kod=1
    prank=rank
    error=0
    GO TO 50

    !     USE 2:ND DERIVATIVES TO RECOMPUTE IF WE ARE ALLOWED TO

30  IF(.NOT.sec) THEN
       error=-4
       kod=secind
    ELSE
       CALL secuc(ffunc,x,dx,diag,p,qraux,n,scale,f,v,c,mdc,m,mdg,nract,error,eval)
       kod=secind
       prank=-(n-nract)
       dxnorm=dnrm2(n,dx,1)
    END IF
50  restar=.false.

    RETURN
  END SUBROUTINE analuc


  !GNAVUC

  SUBROUTINE gnavuc(k,restar,kod,secind,fsum,beta,prekm1,nract,aset,gndok)

    INTEGER, INTENT(IN)      :: k
    LOGICAL, INTENT(IN)      :: restar
    INTEGER, INTENT(IN)      :: kod
    INTEGER, INTENT(OUT)     :: secind
    REAL (dp), INTENT(IN)    :: fsum
    REAL (dp), INTENT(IN)    :: beta
    REAL (dp), INTENT(IN)    :: prekm1
    !tpk INTEGER, INTENT(OUT)     :: nract
    INTEGER, INTENT(INOUT)     :: nract !tpk
    INTEGER, INTENT(IN OUT)  :: aset(:)
    LOGICAL, INTENT(OUT)     :: gndok

    !   THIS IS THE GAUSS-NEWTON DIRECTION ADVISOR ROUTINE
    !   IT ACCEPTS THE GN-DIRECTION AS SEARCH DIRECTION FOR CURRENT
    !   STEP IF CONDITION 1) AND ONE OF CONDITION 2)-3) ARE FULLFILLED

    !   CONDITIONS

    !   1) A PRIMARY DEMAND FOR ACCEPTING THE GN-DIRECTION IS THAT A GN-DIRECTION
    !      WAS USED IN THE PREVIOUS STEP AND NO RESTART WAS DONE IN THIS STEP
    !   2) THE DECREASE IN OBJECTIVE FUNCTION VALUE IN LATEST STEP WAS GOOD ENOUGH
    !      AND WE ARE OUTSIDE THE UMBRELLA
    !   3) A TOLLERABLE VALUE OF CONVERGENCE FACTOR IS OBSERVED

    !   ON ENTRY

    !   K     CONTAINS THE ITERATION COUNT
    !   RESTAR =TRUE IF THIS STEP IS A RESTART STEP
    !          =FALSE IF THIS STEP IS NOT A RESTART STEP
    !   KOD   CONTAINS A CODE THAT SAYS HOW THE PREVIOUS SEARCH
    !         DIRECTION WAS COMPUTED
    !         = 1 IF GN-DIRECTION
    !         =-1 IF SUBSPACE DIRECTION
    !         =-2 IF NEWTON DIRECTION
    !   FSUM  CONTAINS THE SUM OF SQUARES AT CURRENT POINT
    !   BETA  CONTAINS THE NORM OF THE ORTHOGONAL PROJECTION OF F
    !         ONTO THE COLUMN SPACE OF THE JACOBIAN
    !   PREKM1  SEE EXPLANATION IN SUBROUTINE ANALUC
    !   NRACT INTEGER SCALAR CONTAINING THE NO.OF ACTIVE CONSTRAINTS
    !          (POSITIVE OR NEGATIVE)
    !   ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !          CONTAINING +1 OR -1 TO INDICTATE AN ACTIVE BOUND (OTHERWISE ZERO)

    !   ON RETURN

    !   SECIND  INDICATES THE ALTERNATIVE IF GN-DIRECTION IS NOT ACCEPTED
    !         =-1 IF SUBSPACE MINIMIZATION
    !         =-2 IF NEWTONS METHOD
    !   NRACT IS PUT EQUAL TO -NRACT-1 IF NRACT IS <0 ON ENTRY
    !   GNDOK = TRUE IF CURRENT GN-DIRECTION IS ACCEPTED
    !         = FALSE IF RECOMPUTATION OF SEARCH DIRECTION IS NEEDED

    !   COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !   2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !   REPRESENT TIME STEP K-2 RESP. K-1
    !   THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     COMMON VARIABLES TO HOLD SOME INFORMATION CONCERNING
    !     RESTART STEPS

    ! COMMON /back/ icount,lattry,itotal,philat,bestpg,irank
    ! COMMON /INDEX/ imax,ival
    !**************************
    ! COMMON /negdir/ ifree,indic
    !**************************


    !     INTERNAL VARIABLES

    REAL (dp) :: pgress,fnorm,ac
    REAL (dp), PARAMETER  :: c1 = 0.5_dp, c2 = 0.1_dp, c3 = 4.0_dp,  &
         c4 = 10.0_dp, c5 = 0.016_dp

    gndok=.true.
    pgress=-fsum
    fnorm=SQRT(fsum)

    !     CONDITION 1

    IF(.NOT.restar .AND. (k == 0 .OR. imax /= 0)) RETURN
    !**************************
    IF(ifree > 0) THEN
       ifree=ifree-1
       IF(restar) THEN
          secind=-1
          gndok=.false.
       ELSE
          gndok=.true.
       END IF
       RETURN
    END IF
    !**************************
    gndok=.false.
    pgress=fsqkm1 - fsum
    IF(ABS(kod) == 2) GO TO 30
    IF(-1 == kod .OR. restar) GO TO 10
    gndok=.true.

    !     CONDITION 2

    ac=MIN(one,aupkm1)
    IF(pgress > c2*ac*(2.0_dp-ac)*d1km1 .AND. fnorm < c3*beta) return

    !     CONDITION 3

    IF(beta < c1*betkm1) RETURN
    gndok=.false.
    IF(fnorm <= c4*beta) secind=-1
    IF(fnorm > c4*beta) secind=-2
    RETURN

10  secind=-1
    IF(restar) GO TO 20
    IF(alfkm1 < c5*aupkm1 .AND. prekm1 < c2*beta) secind=-2
    RETURN

20  IF(prekm1 < c2*beta) secind=-2
    IF(icount /= 1) RETURN
    IF(imax == 0) RETURN
    aset(imax)=ival
    imax=0
    nract=nract+1
    RETURN

30  secind=kod

    RETURN
  END SUBROUTINE gnavuc


  !SUBUC

  SUBROUTINE subuc(restar,kod,fsum,d1sqs,dxnorm,f,m,c,mdc,n,  &
       scale,diag,p,qraux,prank,dx,v,rank,iter)

    ! N.B. Arguments WORK & W1 have been removed.

    LOGICAL, INTENT(IN OUT)    :: restar
    INTEGER, INTENT(IN OUT)    :: kod
    REAL (dp), INTENT(IN OUT)  :: fsum
    REAL (dp), INTENT(IN OUT)  :: d1sqs
    REAL (dp), INTENT(OUT)     :: dxnorm
    REAL (dp), INTENT(IN OUT)  :: f(:)
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: mdc
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN OUT)  :: diag(:)
    INTEGER, INTENT(IN OUT)    :: p(:)
    REAL (dp), INTENT(IN OUT)  :: qraux(:)
    INTEGER, INTENT(IN OUT)    :: prank
    REAL (dp), INTENT(IN OUT)  :: dx(:)
    REAL (dp), INTENT(IN OUT)  :: v(:)
    INTEGER, INTENT(IN OUT)    :: rank
    INTEGER, INTENT(IN OUT)    :: iter

    !   COMPUTE A SEARCH DIRECTION (DX) BY MINIMIZING IN A SUBSPACE
    !   I.E.  SOLVE THE UPPER TRIANGULAR SYSTEM            T
    !                                             R*DX = -Q *F

    !   BY ONLY USING THE RANK (WILL BE DETERMINED ) FIRST COLUMNS IN MATRIX R
    !   AND SETTING THE REST OF THE ELEMENTS IN DX TO ZERO

    !   ON ENTRY

    !   RESTAR = FALSE IF CURRENT STEP IS NOT A RESTART STEP
    !          = TRUE IF IT IS
    !   KOD  CONTAINS A CODE THAT SAYS HOW THE PREVIOUS SEARCH DIRECTION WAS
    !        COMPUTED
    !        = 1 IF GN-DIRECTION
    !        =-1 IF SUBSPACE DIRECTION
    !        =-2 IF NEWTON DIRECTION
    !   FSUM CONTAINS THE SUM OF SQUARES AT CURRENT POINT
    !   D1SQS CONTAINS PREDICTED REDUCTION IN THE OBJECTIVE FUNCTION
    !         IF GN-DIRECTION IS USED AS SEARCH DIRECTION
    !   DXNORM CONTAINS THE NORM OF THE GN-DIRECTION CONTAINED IN DX
    !   F()  CONTAINS THE VECTOR OF RESIDUALS AT CURRENT POINT
    !   M    CONTAINS LENGTH OF THE ARRAYS F,V AND W1
    !   C(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDC*N  (MDC>=M)
    !        CONTAINING THE MATRIX R FROM THE
    !        DECOMPOSITION    T
    !                        Q *J*D*E = (R)
    !                                   (0)
    !        IN THE UPPER TRIANGLE OF C
    !         WHERE
    !         Q IS ORTHOGONAL (M*M)
    !         J IS THE JACOBIAN (M*N) AT THE TERMINATING POINT
    !         D IS A DIAGONAL MATRIX (N*N)
    !         E IS A PERMUTATION MATRIX (N*N)
    !         R IS AN UPPER TRIANGULAR MATRIX (N*N)
    !   MDC  CONTAINS LEADING DIMENSION OF ARRAY C
    !   N    CONTAINS LENGTH OF THE ARRAYS DIAG,P,QRAUX,DX AND WORK
    !   SCALE CONTAINS =0 IF NO SCALING OF THE JACOBIAN J IS DONE
    !                  =1 IF SCALING IS DONE
    !   DIAG() CONTAINS THE DIAGONAL ELEMENTS OF THE DIAGONAL MATRIX
    !          D ABOVE (IF SCALING IS DONE OTHERWISE UNDEFINED)
    !   P()    CONTAINS THE PERMUTATION MATRIX E ABOVE
    !   QRAUX() CONTAINS INFO. NEEDED TO FORM MATRIX Q ABOVE
    !   PRANK CONTAINS PSEUDO RANK USED TO COMPUTE GN-DIRECTION IN DX
    !   DX()  CONTAINS GN-DIRECTION
    !   V()   CONTAINS THE ORTHOGONAL PROJECTION OF F ONTO THE SPACE
    !         SPANNED BY THE COLUMNS OF THE JACOBIAN OF F(X)

    !   ON RETURN

    !   D1SQS CONTAINS PREDICTED REDUCTION IN THE OBJECTIVE FUNCTION
    !         IF SUGGESTED SEARCH DIRECTION IS USED
    !   DXNORM CONTAINS THE NORM OF THE SUGGESTED SEARCH DIRECTION
    !   DX()  CONTAINS THE SUGGESTED SEARCH DIRECTION
    !   V()   CONTAINS THE ORTHOGONAL PROJECTION OF F ONTO THE SPACE
    !         SPANNED BY RANK LINEARLY INDEPENDENT COLUMNS OF THE
    !         JACOBIAN OF F(X)
    !   RANK  CONTAINS DIMENSION OF THE SUBSPACED USED TO COMPUTE
    !         THE SUGGESTED SEARCH DIRECTION

    !   WORKING AREAS

    !   WORK() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !   W1()   REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M

    REAL (dp) :: work(n), w1(m)

    !   FIND APPROPRIATE SUBSPACE

    CALL dimsub(restar,kod,fsum,f,m,c,n,  &
         scale,diag,p,qraux,prank,dx,v,rank,iter)

    !     COMPUTE SEARCH DIRECTION BY MINIMIZING IN A SUBSPACE OF DIMENSION RANK

    CALL gndunc(f,m,c,n,scale,diag,p,qraux,rank,dx,d1sqs,v,work,w1)
    dxnorm=dnrm2(n,dx,1)

    RETURN
  END SUBROUTINE subuc


  !DIMSUB

  SUBROUTINE dimsub(restar,kod,fsum,f,m,c,n,  &
       scale,diag,p,qraux,prank,dx,v,rank,iter)

    ! Arguments MDC, WORK & W1 have been removed.

    LOGICAL, INTENT(IN)        :: restar
    INTEGER, INTENT(IN)        :: kod
    REAL (dp), INTENT(IN)      :: fsum
    REAL (dp), INTENT(IN)      :: f(:)
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(IN OUT)  :: c(:,:)   ! c(mdc,n)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN)      :: diag(:)
    INTEGER, INTENT(IN)        :: p(:)
    REAL (dp), INTENT(IN OUT)  :: qraux(:)
    INTEGER, INTENT(IN)        :: prank
    REAL (dp), INTENT(IN OUT)  :: dx(:)
    REAL (dp), INTENT(IN OUT)  :: v(:)
    INTEGER, INTENT(OUT)       :: rank
    !tpk INTEGER, INTENT(OUT)       :: iter
    INTEGER, INTENT(INOUT)       :: iter !tpk

    ! COMMON /back/ icount,lattry,itotal,philat,bestpg,irank

    !   FIND APPROPRIATE DIMENSION OF SUBSPACE WHERE MINIMIZATION SHALL BE DONE

    !   ON ENTRY

    !   RESTAR = FALSE IF THIS STEP IS NOT A RESTART STEP
    !            TRUE IF IT IS
    !   KOD  CONTAINS A CODE THAT SAYS HOW THE PREVIOUS SEARCH
    !        DIRECTION WAS COMPUTED
    !        = 1 IF GN-DIRECTION
    !        =-1 IF SUBSPACE DIRECTION
    !        =-2 IF NEWTON DIRECTION
    !   FSUM CONTAINS THE SUM OF SQUARES AT CURRENT POINT
    !   F()  CONTAINS THE VECTOR OF RESIDUALS AT CURRENT POINT
    !   M    CONTAINS LENGTH OF THE ARRAYS F,V AND W1
    !   C(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDC*N  (MDC>=M)
    !        CONTAINING THE MATRIX R FROM THE
    !        DECOMPOSITION    T
    !                        Q *J*D*E = (R)
    !                                   (0)
    !        IN THE UPPER TRIANGLE OF C
    !         WHERE
    !         Q IS ORTHOGONAL (M*M)
    !         J IS THE JACOBIAN (M*N) AT THE TERMINATING POINT
    !         D IS A DIAGONAL MATRIX (N*N)
    !         E IS A PERMUTATION MATRIX (N*N)
    !         R IS AN UPPER TRIANGULAR MATRIX (N*N)
    !   MDC  CONTAINS LEADING DIMENSION OF ARRAY C
    !   N    CONTAINS LENGTH OF THE ARRAYS QRAUX,DX,P,DIAG AND WORK
    !   SCALE  =0 IF NO SCALING OF THE JACOBIAN J IS DONE
    !          =1 IF COLUMN SCALING IS DONE
    !   DIAG() CONTAINS THE DIAGONAL ELEMENTS OF THE DIAGONAL MATRIX
    !          D ABOVE (IF SCALING IS DONE, OTHERWISE UNDEFINED)
    !   P()    CONTAINS THE PERMUTATION MATRIX E FROM ABOVE
    !   QRAUX() CONTAINS INFO. NEEDED TO FORM MATRIX Q ABOVE
    !   PRANK CONTAINS PSEUDO RANK USED TO COMPUTE GN-DIRECTION IN DX
    !   DX()  CONTAINS GN-DIRECTION     T
    !   V()   CONTAINS   J*DX

    !   ON RETURN

    !   DX() COMPLETELY DESTROYED
    !   RANK CONTAINS SUGGESTED DIMENSION OF SUBSPACE WHERE THE
    !        MINIMIZATION SHALL BE DONE

    !   WORKING AREAS

    !   WORK() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !   W1()   REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M


    !   COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !   2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !   REPRESENT TIME STEP K-2 RESP. K-1
    !   THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     INTERNAL VARIABLES

    INTEGER   :: i,info,j,mindim
    REAL (dp) :: bn,sn,dummy(1),pgress,work(n),w1(m)
    REAL (dp), SAVE  :: rabs = 0.1_dp

    !     CHECK IF A RESTART STEP

    ! ===================================================================
    ! Modified by T.KUROSU to fix a strange Infinite-Loop bug encountered
    ! when compiling this routine with pgf90 and optimizaion.
    !
    ! Date of modification: 26 October 2004
    ! ===================================================================

    IF ( .NOT. restar ) THEN 

       !   IN THIS POSITION WE KNOW THAT PREVIOUS SEARCH DIRECTION WAS
       !   COMPUTED BY MINIMIZING IN A SUBSPACE OR WAS THE GN-DIRECTION
       !
       !   TO DETERMINE A SUITABLE SUBSPACE WE TRY TO ESTIMATE THE
       !   BEST REDUCTION IN THE OBJECTIVE FUNCTION BY INVESTIGATING
       !   THE RIGHT HAND SIDE OF THE UPPER TRIANGULAR
       !   SYSTEM                T
       !                R*DX = -Q *F          (1)
       !
       !   RNGKM1 = PSEUDO RANK USED TO COMPUTE PREVIOUS SEARCH DIRECTION


       !   FORM RIGHT HAND SIDE OF (1) AND STORE IN W1

       w1(1:m)=-f(1:m)
       CALL sqrsl(c,m,prank,qraux,w1,dummy,w1,dummy,dummy,dummy,1000,info)

       !     COMPUTE ESTIMATES OF STEPLENGTHS W1(I) AND PROGRESS WORK(I)

       work(1)=w1(1)
       j=p(1)
       IF(scale == 0) w1(1)=w1(1)/c(1,1)
       IF(scale /= 0) w1(1)=w1(1)*diag(j)/c(1,1)
       DO  i=2,prank
          work(i)=w1(i)
          IF(scale /= 0) THEN
             j=p(i)
             w1(i)=w1(i)*diag(j)/c(i,i)
          ELSE
             w1(i)=w1(i)/c(i,i)
          END IF
          w1(i)=dnrm2(2,w1(i-1:),1)
          work(i)=dnrm2(2,work(i-1:),1)
       END DO
       !tpk sn=w1(prank)
       !tpk bn=work(prank)
       sn = 0.0_dp ; bn = 0.0_dp
       IF ( prank >= 1 ) THEN
          sn=w1(prank) ; bn=work(prank)
       END IF

       !     DETERMINE THE LOWEST POSSIBLE DIMENSION
       
       DO  i=1,prank
          mindim=i
          IF(work(i) > rabs*bn) EXIT
       END DO

       IF(kod == 1) CALL pregn(w1,sn,work,bn,mindim,prank,rank)
       IF(-1 == kod) CALL presub(w1,work,bn,rabs,prank,fsum,rank)
       RETURN !GO TO 100
    END IF

    !     SUGGEST NEW PSEUDO RANK = LATEST PSEUDO RANK -1
    !     ALSO SAVE PSEUDO RANK THAT REDUCES THE OBJECTIVE BEST

!!tpk60  IF(icount > 1) GO TO 80
!!tpk70  bestpg=0.5_dp*fsum - philat
!!tpk    irank=lattry
!!tpk80  pgress=0.5_dp*fsum - philat
!!tpk    PRINT *, pgress, bestpg, philat
!!tpk    IF(pgress > bestpg) GO TO 70
!!tpk    rank=lattry-1
!!tpk    IF(lattry <= 1) THEN
!!tpk       rank=irank
!!tpk       iter=iter+1
!!tpk       lattry=0
!!tpk    END IF

    pgress=0.5_dp*fsum - philat
    IF ( icount <= 1 .OR. pgress > bestpg ) THEN
       bestpg=0.5_dp*fsum - philat
       irank=lattry
    END IF
    rank=lattry-1
    IF(lattry <= 1) THEN
       rank=irank
       iter=iter+1
       lattry=0
    END IF
!tpk    IF(restar) GO TO 60
!tpk
!tpk    !   IN THIS POSITION WE KNOW THAT PREVIOUS SEARCH DIRECTION WAS
!tpk    !   COMPUTED BY MINIMIZING IN A SUBSPACE OR WAS THE GN-DIRECTION
!tpk
!tpk    !   TO DETERMINE A SUITABLE SUBSPACE WE TRY TO ESTIMATE THE
!tpk    !   BEST REDUCTION IN THE OBJECTIVE FUNCTION BY INVESTIGATING
!tpk    !   THE RIGHT HAND SIDE OF THE UPPER TRIANGULAR
!tpk    !   SYSTEM                T
!tpk    !                R*DX = -Q *F          (1)
!tpk
!tpk    !   RNGKM1 = PSEUDO RANK USED TO COMPUTE PREVIOUS SEARCH DIRECTION
!tpk
!tpk
!tpk    !   FORM RIGHT HAND SIDE OF (1) AND STORE IN W1
!tpk
!tpk    w1(1:m)=-f(1:m)
!tpk    CALL sqrsl(c,m,prank,qraux,w1,dummy,w1,dummy,dummy,dummy,1000,info)
!tpk
!tpk    !     COMPUTE ESTIMATES OF STEPLENGTHS W1(I) AND PROGRESS WORK(I)
!tpk
!tpk    work(1)=w1(1)
!tpk    j=p(1)
!tpk    IF(scale == 0) w1(1)=w1(1)/c(1,1)
!tpk    IF(scale /= 0) w1(1)=w1(1)*diag(j)/c(1,1)
!tpk    DO  i=2,prank
!tpk       work(i)=w1(i)
!tpk       IF(scale /= 0) THEN
!tpk          j=p(i)
!tpk          w1(i)=w1(i)*diag(j)/c(i,i)
!tpk       ELSE
!tpk          w1(i)=w1(i)/c(i,i)
!tpk       END IF
!tpk       w1(i)=dnrm2(2,w1(i-1:),1)
!tpk       work(i)=dnrm2(2,work(i-1:),1)
!tpk    END DO
!tpk    !tpk sn=w1(prank)
!tpk    !tpk bn=work(prank)
!tpk    sn = 0.0_dp ; bn = 0.0_dp
!tpk    IF ( prank >= 1 ) THEN
!tpk       sn=w1(prank) ; bn=work(prank)
!tpk    END IF
!tpk
!tpk    !     DETERMINE THE LOWEST POSSIBLE DIMENSION
!tpk
!tpk    DO  i=1,prank
!tpk       mindim=i
!tpk       IF(work(i) > rabs*bn) EXIT
!tpk    END DO
!tpk
!tpk    IF(kod == 1) CALL pregn(w1,sn,work,bn,mindim,prank,rank)
!tpk    IF(-1 == kod) CALL presub(w1,work,bn,rabs,prank,fsum,rank)
!tpk    GO TO 100
!tpk
!tpk    !     SUGGEST NEW PSEUDO RANK = LATEST PSEUDO RANK -1
!tpk    !     ALSO SAVE PSEUDO RANK THAT REDUCES THE OBJECTIVE BEST
!tpk
!tpk60  IF(icount > 1) GO TO 80
!tpk70  bestpg=0.5_dp*fsum - philat
!tpk    irank=lattry
!tpk80  pgress=0.5_dp*fsum - philat
!tpk    IF(pgress > bestpg) GO TO 70
!tpk    rank=lattry-1
!tpk    IF(lattry <= 1) THEN
!tpk       rank=irank
!tpk       iter=iter+1
!tpk       lattry=0
!tpk    END IF

100 RETURN
  END SUBROUTINE dimsub


  !PREGN

  SUBROUTINE pregn(s,sn,b,bn,mindim,prank,dim)

    !   GN-STEP IN PREVIOUS STEP
    !   TAKE DIM AS THE LARGEST K (MINDIM<=K<=PRANK-1) FOR WHICH
    !   S(K) < SMAX*SN AND B(K) > RMIN*BN
    !   IF NO SUCH K EXISTS TAKE DIM=PRANK-1 PROVIDED (PRANK-1) >= MINDIM

    REAL (KIND=dp), INTENT(IN OUT)  :: s(:)
    REAL (KIND=dp), INTENT(IN OUT)  :: sn
    REAL (KIND=dp), INTENT(IN OUT)  :: b(:)
    REAL (KIND=dp), INTENT(IN OUT)  :: bn
    INTEGER, INTENT(IN)        :: mindim
    INTEGER, INTENT(IN)        :: prank
    INTEGER, INTENT(OUT)       :: dim

    REAL (dp), SAVE  :: smax = 0.2_dp, rmin = 0.5_dp
    INTEGER          :: i, k, m1

    m1=prank-1
    k=prank
    IF(mindim > m1) GO TO 20
    DO  i=mindim,m1
       k=m1-i+mindim
       IF(s(k) < smax*sn .AND. b(k) > rmin*bn) GO TO 20
    END DO
    dim=MAX(mindim,prank-1)
    RETURN

20  dim=k

    RETURN
  END SUBROUTINE pregn


  !PRESUB

  SUBROUTINE presub(s,b,bn,rabs,prank,fsum,dim)

    REAL (dp), INTENT(IN)  :: s(:)
    REAL (dp), INTENT(IN)  :: b(:)
    REAL (dp), INTENT(IN)  :: bn
    REAL (dp), INTENT(IN)  :: rabs
    INTEGER, INTENT(IN)    :: prank
    REAL (dp), INTENT(IN)  :: fsum
    INTEGER, INTENT(OUT)   :: dim

    !     SUBSPACE MINIMIZATION IN LATEST STEP

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    REAL (dp), SAVE  :: stepb = 0.2_dp, pgb1 = 0.3_dp, pgb2 = 0.1_dp,  &
         predb = 0.7_dp, rlenb = 2.0_dp
    INTEGER          :: i, i1
    REAL (dp)        :: pgress

    !     IF THE LATEST STEP WAS FAIRLY GOOD THE DIMENSION MUST NOT BE DECREASED

    pgress=SQRT(fsqkm1-fsum)
    IF(alfkm1 >= stepb .OR. pgress > pgb1*SQRT(d1km1) .OR.  &
         pgress > pgb2*betkm1 ) GO TO 10

    !     A BAD STEP

    dim=MAX(1,rngkm1-1)
    IF(rngkm1 > 1 .AND. b(dim) > rabs*bn) RETURN

10  dim=rngkm1
    IF(b(dim) > predb*bn .AND. rlenb*s(dim) < s(dim+1)) RETURN
    i1=rngkm1+1
    DO  i=i1,prank
       dim=i
       IF(b(i) > predb*bn) EXIT
    END DO

    RETURN
  END SUBROUTINE presub


  !SECUC

  SUBROUTINE secuc(ffunc,x,dx,diag,p,qraux,n,scale,f,v,c,mdc,m,mdg,nract,  &
       error,eval)

    ! N.B. Arguments W1, W3 & GMAT have been removed.

    REAL (dp), INTENT(IN OUT)  :: x(:)
    REAL (dp), INTENT(OUT)     :: dx(:)
    REAL (dp), INTENT(IN)      :: diag(:)
    INTEGER, INTENT(IN)        :: p(:)
    REAL (dp), INTENT(OUT)     :: qraux(:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: scale
    REAL (dp), INTENT(IN OUT)  :: f(:)
    REAL (dp), INTENT(OUT)     :: v(:)
    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: mdc
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: mdg
    INTEGER, INTENT(IN)        :: nract
    INTEGER, INTENT(OUT)       :: error
    INTEGER, INTENT(OUT)       :: eval

    EXTERNAL ffunc

    
    !   COMPUTE THE SOLUTION (DX) OF THE N*N SYSTEM
    !                         T
    !            GMAT*DX =  -J *F               (1)
    !   WHERE
    !            T       M
    !    GMAT = J *J + SIGMA(F *G )
    !                   K=1   K  K
    !   AND

    !   J IS THE M*N JACOBIAN OF F(X) AT THE CURRENT POINT
    !   G  (N*N) IS THE K:TH HESSIAN OF F(X) AT THE CURRENT POINT
    !    K
    !        WE KNOW THE QR-DECOMPOSITION OF THE JACOBIAN J ON ENTRY

    !   ON ENTRY

    !   FFUNC SUBROUTINE NAME USED TO EVALUATE THE FUNCTIONS AT A CERTAIN POINT
    !   X()   CONTAINS THE CURRENT POINT
    !   DIAG() CONTAINS THE DIAGONAL ELEMENTS OF THE DIAGONAL MATRIX
    !          D BELOW (IF SCALING IS DONE. OTHERWISE UNDEFINED)
    !   P()    CONTAINS THE PERMUTATION MATRIX E BELOW
    !   QRAUX() CONTAINS INFO. NEEDED TO FORM MATRIX Q BELOW
    !   N    CONTAINS LENGTH OF THE ARRAYS X,DX,DIAG,P,QRAUX,W3
    !   SCALE CONTAINS = 0 IF NO SCALING OF THE JACOBIAN IS DONE
    !                  = 1 IF SCALING IS DONE
    !   F()  CONTAINS THE VECTOR OF RESIDUALS AT CURRENT POINT
    !   V()  CONTAINS THE ORTHOGONAL PROJECTION OF -F ONTO THE
    !        SPACE SPANNED BY THE COLUMNS THE JACOBIAN AT CURRENT POINT
    !   C(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDC*N  (MDC>=M)
    !        CONTAINING THE MATRIX R FROM THE
    !        DECOMPOSITION    T
    !                        Q *J*D*E = (R)
    !                                   (0)           (2)
    !        IN THE UPPER TRIANGLE OF C
    !         WHERE
    !         Q IS ORTHOGONAL (M*M)
    !         J IS THE JACOBIAN (M*N) AT THE CURRENT POINT
    !         D IS A DIAGONAL MATRIX (N*N)
    !         E IS A PERMUTATION MATRIX (N*N)
    !         R IS AN UPPER TRIANGULAR MATRIX (N*N)
    !   MDC     INTEGER SCALAR CONTAINING LEADING DIMENSION OF ARRAY C
    !           MDC MUST BE >= M
    !   M       INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS F,V,W1
    !   MDG     INTEGER SCALAR CONTAINING LEADING DIMENSION OF ARRAY GMAT.
    !           MDG MUST BE >= N
    !   NRACT    INTEGER SCALAR CONTAINING NO. OF ACTIVE CONSTRAINTS

    !   ON RETURN

    !   DX() THE COMPUTED SOLUTION DX OF SYSTEM (1)
    !   QRAUX()  COMPLETELY DESTROYED
    !   V()      COMPLETELY DESTROYED
    !   C(,)     COMPLETELY DESTROYED
    !   ERROR    CONTAINS = -3 IF THE MATRIX GMAT IN SYSTEM (1)
    !                          IS NOT POSITIVE DEFINITE
    !                     < -10 IF A USER STOP IS INDICATED
    !                     =  0 IF GMAT IS POSITIVE DEFINITE
    !   EVAL     CONTAINS NO. OF FUNCTION EVALUATIONS DONE INSIDE THIS ROUTINE

    !   WORKING AREAS

    !   W1()    REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M
    !   W3()    REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !   GMAT(,) REAL DOUBLY SUBSCRIPTED ARRAY OF DIMENSION MDG*N

    REAL (dp)  :: w1(m), w3(n), gmat(n,n)

    !   INTERNAL VARIABLES

    INTEGER   :: i,j,k,l,nn,info,idummy(1)
    REAL (dp) :: dummy(n)

    !     INSTEAD OF SOLVING (1) WE SOLVE A TRANSFORMED SYSTEM
    !     WHICH HAS N-NRACT UNKNOWNS DY

    !     I.E. SOLVE  T    T      M                        T     T
    !               (R *R+E *D*(SIGMA(F *G ))*D*E)*DY = -(R :0)*Q *F  (2)
    !                            K=1   K  K

    !           DX = D*E*DY  (THE LAST NRACT ELEMENTS OF DY EQUAL TO ZERO)

    !     FORM RIGHT HAND SIDE OF SYSTEM (2)
    !     I.E.          T     T
    !          FORM  -(R :0)*Q *F   AND STORE IN DX
    !                 T
    !     FIRST FORM Q *F BY USING SQRSL AND STORE IN W1

    nn=n - nract
    CALL sqrsl(c,m,nn,qraux,f,dummy,w1,dummy,dummy,dummy,1000,info)

    !              T
    !     FORM  -(R :0)*W1 AND STORE IN DX

    CALL rtrw1(c,nn,w1,dx)
    !                T
    !     FORM  C = R *R

    CALL jtrj(c,nn,qraux)
    !      write(10,*) 'R(tr)*R'
    !      do 6 i=1,nn
    !       write(10,*) (c(i,j),j=1,nn)
    !    6 continue

    !     SAVE THE 4 FIRST COLUMNS OF THE MATRIX C

    DO  i=1,n
       qraux(i)=c(i,1)
       v(i)=c(i,2)
       w1(i)=c(i,3)
       w3(i)=c(i,4)
    END DO

    !     COMPUTE HESSIANS G   AND FORM FINAL MATRIX GMAT OF SYSTEM (2)
    !                       K

    !     THE 4 FIRST COLUMNS OF ARRAY C ARE USED AS WORKING STORAGE
    !     INSIDE THE ROUTINE HESS

    error=0
    CALL hess(ffunc,gmat,x,n,f,c(:,1),c(:,2),c(:,3),c(:,4),m,error)
    IF(error < -10) RETURN
    eval=2*n*(n+1)

    !     RESAVE THE 4 COLUMNS OF MATRIX C

    DO  i=1,n
       c(i,1)=qraux(i)
       c(i,2)=v(i)
       c(i,3)=w1(i)
       c(i,4)=w3(i)
    END DO
    !           T
    !     FORM E *D*GMAT*D*E

    IF(scale == 0) GO TO 50
    DO  i=1,n
       DO  j=1,n
          gmat(i,j)=diag(i)*gmat(i,j)
          gmat(j,i)=gmat(j,i)*diag(i)
       END DO
    END DO

    !     DO THE PIVOTING OF COLUMNS AND ROWS IN MATRIX GMAT
    !     AND ADD TO MATRIX C

    !      write(10,*) 'The nonlinear part. Only the upper triangular.'
    !      write(10,*) 'Row-wize written. I.e. the first 4 values make'
    !      write(10,*) 'up the 1st row. The following 3 values make up'
    !      write(10,*) 'the 2nd row starting at the diagonal elem.'
50  DO  i=1,nn
       l=p(i)
       DO  j=i,nn
          k=p(j)
          c(i,j)=c(i,j)+gmat(l,k)
          !      write(10,*) GMAT(L,K)
       END DO
    END DO

    !     PERFORM CHOLESKY DECOMPOSITION OF MATRIX C

    CALL schdc(c,nn,qraux,idummy,0,info)
    error=0
    IF(nn == info) GO TO 90

    !     MATRIX C IS NOT POSITIVE DEFINITE

    error=-3
    RETURN
    !                              T    T
    !     SOLVE SYSTEM   C*DY = -(R 0)*Q *F

90  CALL sposl(c,nn,dx)

    !     FORM  DX = D*E*DY

    IF(nract == 0) GO TO 110
    DO  i=1,nract
       j=nn+i
       dx(j)=zero
    END DO

110 CALL pivec(p,n,dx,qraux)
    IF(scale == 0) GO TO 130
    DO  i=1,n
       dx(i)=diag(i)*dx(i)
    END DO

130 RETURN
  END SUBROUTINE secuc


  !RTRW1

  SUBROUTINE rtrw1(c,n,w1,dx)

    ! N.B. Arguments MDC & M have been removed.

    REAL (dp), INTENT(IN)   :: c(:,:)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: w1(:)
    REAL (dp), INTENT(OUT)  :: dx(:)

    !                 T
    !   FORM  DX:= -(R :0)*W1
    !   THE N*N UPPER TRIANGULAR MATRIX R IS CONTAINED IN THE UPPER
    !   LEFT TRIANGLE OF ARRAY C

    !   INTERNAL VARIABLES

    INTEGER  :: j

    DO  j=1,n
       dx(j)=-DOT_PRODUCT( c(1:j,j), w1(1:j) )
    END DO

    RETURN
  END SUBROUTINE rtrw1


  !JTRJ

  SUBROUTINE jtrj(c,n,w)

    ! N.B. Argument MDC has been removed.

    !                                   T
    !   FORM THE N*N SYMMETRIC MATRIX  C *C   AND STORE IN C
    !
    !                                 T
    !   FIRST FORM THE LOWER PART OF C *C  AND STORE IN THE LOWER PART OF C

    REAL (dp), INTENT(IN OUT)  :: c(:,:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(OUT)     :: w(:)

    !   INTERNAL VARIABLES

    INTEGER :: i,j,k

    DO  j=1,n
       w(1:n)=c(1:n,j)
       DO  k=j,n
          c(k,j)=DOT_PRODUCT( c(1:j,k), w(1:j) )
       END DO
    END DO

    !     MOVE THE LOWER PART OF C TO THE UPPER PART OF C

    IF(n == 1) RETURN
    DO  i=2,n
       k=i-1
       DO  j=1,k
          c(j,i)=c(i,j)
       END DO
    END DO

    RETURN
  END SUBROUTINE jtrj


  !HESS

  SUBROUTINE hess(ffunc,b,x,n,v,f1,f2,f3,f4,m,error)

    ! Argument MDB has been removed.

    !   COMPUTE THE N*N MATRIX          M
    !                          B :=   SIGMA(V *G )
    !                                  K=1   K  K
    !   WHERE G   IS THE HESSIAN OF F (X)
    !          K                     K

    REAL (dp), INTENT(OUT)     :: b(:,:)
    REAL (dp), INTENT(IN OUT)  :: x(:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN)      :: v(:)
    REAL (dp), INTENT(OUT)     :: f1(:)
    REAL (dp), INTENT(OUT)     :: f2(:)
    REAL (dp), INTENT(OUT)     :: f3(:)
    REAL (dp), INTENT(OUT)     :: f4(:)
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(OUT)       :: error

    EXTERNAL  ffunc

    

    !   COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !   SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    !     INTERNAL VARIABLES

    INTEGER   :: j,k,l,ctrl
    REAL (dp) :: athird,eps1,eps2,epsk,epsj,xk,xj,term,dummy(1,1),sum

    ctrl=-1
    athird=one/3.0_dp
    eps2=srelpr**athird
    eps1=eps2
    DO  k=1,n
       xk=x(k)
       epsk=MAX(ABS(xk),one)*eps2
       DO  j=1,k
          xj=x(j)
          epsj=MAX(ABS(xj),one)*eps1
          x(k)=xk+epsk
          x(j)=x(j)+epsj
          CALL ffunc(x,n,f1,m,ctrl,dummy,1)
          IF(ctrl < -10) GO TO 40
          x(k)=xk
          x(j)=xj
          x(k)=xk+epsk
          x(j)=x(j)-epsj
          CALL ffunc(x,n,f2,m,ctrl,dummy,1)
          IF(ctrl < -10) GO TO 40
          x(k)=xk
          x(j)=xj
          x(k)=xk-epsk
          x(j)=x(j)+epsj
          CALL ffunc(x,n,f3,m,ctrl,dummy,1)
          IF(ctrl < -10) GO TO 40
          x(k)=xk
          x(j)=xj
          x(k)=xk - epsk
          x(j)=x(j) - epsj
          CALL ffunc(x,n,f4,m,ctrl,dummy,1)
          IF(ctrl < -10) GO TO 40
          x(k)=xk
          x(j)=xj
          sum=zero
          DO  l=1,m
             term=f1(l)-f2(l)-f3(l)+f4(l)
             sum=sum + term/(4.0_dp*epsk*epsj)*v(l)
          END DO
          b(k,j)=sum
          IF(k == j) CYCLE
          b(j,k)=sum
       END DO
    END DO
    RETURN

    !     USER STOP INDICATION IS DETECTED

40  error=ctrl

    RETURN
  END SUBROUTINE hess


  !TERMUC

  SUBROUTINE termuc(prank,error,restar,maxit,k,fsum,d1sqs,n,xnorm,  &
       alfnoi,xdiff,epsabs,epsrel,epsx,g,bl,bu,aset,nract,EXIT)

    INTEGER, INTENT(IN)        :: prank
    INTEGER, INTENT(IN)        :: error
    LOGICAL, INTENT(IN)        :: restar
    INTEGER, INTENT(IN)        :: maxit
    INTEGER, INTENT(IN)        :: k
    REAL (dp), INTENT(IN)      :: fsum
    REAL (dp), INTENT(IN)      :: d1sqs
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN)      :: xnorm
    REAL (dp), INTENT(IN OUT)  :: alfnoi
    REAL (dp), INTENT(IN)      :: xdiff
    REAL (dp), INTENT(IN)      :: epsabs
    REAL (dp), INTENT(IN)      :: epsrel
    REAL (dp), INTENT(IN)      :: epsx
    REAL (dp), INTENT(IN OUT)  :: g(:)
    REAL (dp), INTENT(IN)      :: bl(:)
    REAL (dp), INTENT(IN)      :: bu(:)
    INTEGER, INTENT(IN OUT)    :: aset(:)
    INTEGER, INTENT(IN)        :: nract
    INTEGER, INTENT(OUT)       :: EXIT

    !   CHECK IF ANY OF THE TERMINATION CRITERIA ARE SATISFIED
    !   THE CONVERGENCE CRITERIA ARE ONLY CHECKED IF THE LATEST SEARCH
    !   DIRECTION WAS COMPUTED USING THE GAUSS-NEWTON WITH FULL PSEUDO
    !   RANK OR IF THE METHOD OF NEWTON WAS USED AND PROVIDED THAT NO
    !   RESTART WAS DONE


    !   THE CONVERGENCE CRITERIA ARE

    !   1) RELATIVE PREDICTED REDUCTION IN THE OBJECTIVE FUNCTION
    !      IS LESS THAN EPSREL**2
    !   2) THE SUM OF SQUARES IS LESS THAN EPSABS**2
    !   3) THE RELATIVE CHANGE IN X IS LESS THAN EPSX
    !   4) WE ARE COMPUTING AT NOISE LEVEL
    !      THE LAST DIGIT IN THE CONVERGENCE CODE (SEE BELOW) INDICATES
    !      HOW THE LAST STEPS WERE COMPUTED
    !      = 0 NO TROBLE (GAUSS-NEWTON THE LAST 3 STEPS)
    !      = 1 PRANK<>N AT THE TERMINATION POINT
    !      = 2 THE METHOD OF NEWTON WAS USED (AT LEAST) IN THE LAST STEP
    !      = 3 THE 2:ND BUT LAST STEP WAS SUBSPACE MINIMIZATION BUT THE
    !          LAST TWO WERE GAUSS-NEWTON STEPS
    !      = 4 THE STEPLENGTH WAS NOT UNIT IN BOTH THE LAST TWO STEPS


    !    THE ABNORMAL TERMINATION CRITERIA ARE

    !   5) NO. OF ITERATIONS HAS EXCEEDED MAXIMUM ALLOWED ITERATIONS
    !   6) THE HESSIAN EMANATING FROM 2.ND ORDER METHOD IS NOT POS. DEF.
    !   7) THE ALGORITHM WOULD LIKE TO USE 2:ND DERIVATIVES BUT IS
    !      NOT ALLOWED TO DO THAT
    !   8) AN UNDAMPED STEP WITH NEWTONS METHOD IS A FAILURE
    !   9) THE LATEST SEARCH DIRECTION COMPUTED USING SUBSPACE
    !      MINIMIZATION WAS NOT A DESCENT DIRECTION (PROBABLY CAUSED
    !   BY WRONGLY COMPUTED JACOBIAN)

    !   ON ENTRY

    !   PRANK   INTEGER SCALAR CONTAINING PSEUDO RANK USED TO COMPUTE
    !           THE SEARCH DIRECTION IF SUBSPACE MINIMIZATION IS USED
    !           = -N IF SECOND ORDER METHOD IS USED
    !   ERROR   INTEGER SCALAR CONTAINING
    !          =-1 IF PREVIOUS DIRECTION WAS NOT A DESCENT DIRECTION
    !          =-2 IF THE CHANGE IN X < SQRT(SRELPR)
    !          =-3 IF NOT POS. DEF. HESSIAN FROM SECOND ORDER METHOD
    !          =-4 IF THE ALGORITHM WANTS TO USE SECOND DERIVATIVES
    !              BUT IS NOT ALLOWED TO DO THAT
    !          =-5 IF AN UNDAMPED NEWTON STEP FAILS
    !          =0 OTHERWISE
    !   RESTAR  LOGICAL SCALAR =TRUE IF A RESTART STEP
    !           =FALSE IF THIS STEP IS NOT A RESTART STEP
    !   maxit     INTEGER SCALAR CONTAINING MAXIMUM ALLOWED ITERATIONS
    !   K       INTEGER SCALAR CONTAINING ITERATION NO.
    !   FSUM    REAL SCALAR CONTAINING SUM OF SQUARES AT CURRENT POINT
    !   D1SQS   REAL SCALAR CONTAINING PREDICTED REDUCTION IN THE
    !           OBJECTIVE FUNCTION IF ONE MORE GAUSS-NEWTON STEP
    !           IS TAKEN
    !   N       INTEGER SCALAR CONTAINING NUMBER OF PARAMETERS
    !   XNORM   REAL SCALAR CONTAINING NORM OF CURRENT POINT
    !   ALFNOI  REAL SCALAR CONTAINING A NOISE LEVEL
    !   XDIFF   REAL SCALAR CONTAINING CHANGE IN LATEST SUCCESIVE POINTS
    !   EPSABS  REAL SCALARS CONTAINING SMALL POSITIVE VALUES
    !   EPSREL  USED IN CHECKING THE CONVERGENCE CRITERIA
    !   EPSX
    !   G()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE GRADIENT OF THE OBJECTIVE AT THE CURRENT POINT
    !   BL() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE LOWER BOUNDS OF THE UNKNOWNS
    !   BU() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !        THE UPPER BOUNDS OF THE UNKNOWNS
    !   NRACT INTEGER SCALAR CONTAINING THE NO. OF ACTIVE CONSTRAINTS
    !   ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !         CONTAINING A CODE WHICH INDICATES WHETHER AN UNKNOWN IS
    !         ACTIVE OR NOT.
    !         ASET(I) = 0 WHEN X(I) IS FREE
    !                 =+1 WHEN X(I) EQUALS BL(I)
    !                 =-1 WHEN X(I) EQUALS BU(I)

    !   ON RETURN

    !   EXIT    INTEGER SCALAR CONTAINING
    !         = 0 IF NO TERMINATION CRITERIUM IS SATISFIED
    !         = 10000 CRITERIUM 1) IS SATISFIED
    !         =  2000     #     2)      #
    !         =   300     #     3)      #
    !         =    40     #     4)      #
    !         =     X  WHERE X EQUALS 0,1,2,3 OR 4
    !         =    -2     #     5)      #
    !         =    -3     #     6)      #
    !         =    -4     #     7)      #
    !         =    -5     #     8)      #
    !         =    -6     #     9)      #


    !     COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !     2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !     REPRESENT TIME STEP K-2 RESP. K-1
    !     THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    !     COMMON VARIABLES CONTAINING INFORMATION FOR RESTART STEPS

    ! COMMON /INDEX/ imax,ival

    REAL (dp)        :: rlmax
    INTEGER          :: i
    REAL (dp), SAVE  :: noise = 0.1_dp, eps = 0.01_dp, epsl = 0.001_dp

    EXIT=0
    rlmax=zero

    !     THE CONVERGENCE CRITERIA ARE NOT CHECKED IF THE LATEST STEP
    !     WAS A RESTART STEP OR WAS NOT A FULL PSEUDO RANK STEP

    IF(imax /= 0) GO TO 60
    IF(restar) GO TO 60
    IF(-1 == kodkm1 .AND. alfnoi <= noise) GO TO 60
    IF(error < 0 .AND. -2 /= error) GO TO 60
    !      write(10,*) 'In TERMUC:'

    !     CRITERIUM NO. 1

    IF(d1sqs <= fsum*epsrel**2) EXIT=EXIT+10000

    !     CRITERIUM NO. 2

    IF(fsum <= epsabs**2) EXIT=EXIT+2000

    !     CRITERIUM NO. 3

    IF(xdiff < xnorm*epsx) EXIT=EXIT+300

    !     CRITERIUM NO. 4

    IF(alfnoi > noise .OR. -2 == error) EXIT=EXIT+40
    IF(EXIT == 0) GO TO 60

    !     CHECK LAGRANGE MULTIPLIERS

    !      write(10,*) 'Checking Lagrange mult. in TERMUC'
    !      write(10,*) (g(i),i=1,n)
    IF(nract == 0) GO TO 15
    DO  i=1,n
       IF(ABS(g(i)) < epsl) CYCLE
       IF(g(i)*aset(i)*(bu(i)-bl(i)) >= zero) CYCLE
       IF(rlmax > ABS(g(i))) CYCLE
       rlmax=ABS(g(i))
       EXIT=0
    END DO

15  IF(EXIT == 0) GO TO 60

    !     DETERMINE THE LAST DIGIT IN THE CONVERGENCE CODE

    IF(ABS(kodkm1) /= 2) GO TO 20
    EXIT=EXIT+2
    RETURN

20  IF(prank == (n-nract)) GO TO 30
    EXIT=EXIT+1
    RETURN

30  IF(kodkm2 == 1 .OR. k < 2) GO TO 40
    EXIT=EXIT+3
    RETURN

40  IF(ABS(alfkm2-one) <= eps .AND. ABS(alfkm1-one) <= eps) RETURN
    EXIT=EXIT+4
    RETURN

    !     CRITERIUM NO. 5

60  IF(k > maxit) EXIT=-2

    !     CRITERIUM NO. 9

    IF(-1 == error) EXIT=-6

    !     CRITERIUM NO. 6-8

    IF(error >= -5 .AND. error <= -3) EXIT=error

    RETURN
  END SUBROUTINE termuc


  !STEPUC

  SUBROUTINE stepuc(x,g,dx,n,f,v1,m,phi,ffunc,kod,prank,bl,bu,  &
       aset,nract,bnd,eval,alpha,alphup,xdiff,error)

    ! N.B. Arguments FNEW, V2 & GMOD have been removed.

    REAL (dp), INTENT(IN OUT)  :: x(:)
    REAL (dp), INTENT(IN)      :: g(:)
    REAL (dp), INTENT(IN)      :: dx(:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(OUT)     :: f(:)
    REAL (dp), INTENT(IN OUT)  :: v1(:)
    INTEGER, INTENT(IN)        :: m
    !tpk REAL (dp), INTENT(OUT)     :: phi
    REAL (dp), INTENT(INOUT)     :: phi !tpk
    INTEGER, INTENT(IN OUT)    :: kod
    INTEGER, INTENT(IN)        :: prank
    REAL (dp), INTENT(IN)      :: bl(:)
    REAL (dp), INTENT(IN)      :: bu(:)
    INTEGER, INTENT(IN)        :: aset(:)
    INTEGER, INTENT(IN OUT)    :: nract
    INTEGER, INTENT(IN)        :: bnd
    INTEGER, INTENT(OUT)       :: eval
    REAL (dp), INTENT(OUT)     :: alpha
    REAL (dp), INTENT(OUT)     :: alphup
    REAL (dp), INTENT(OUT)     :: xdiff
    INTEGER, INTENT(IN OUT)    :: error

    EXTERNAL ffunc

    
    !     COMPUTE THE STEPLENGTH ALPHA IN THE ITERATION
    !            XNEW:=XOLD+ALPHA*DX           (1)

    !     ON ENTRY
    !     PARAMETERS FLAGGED WITH (*) ARE ONLY DEFINED IF ABS(KOD)=1

    !     X()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          THE CURRENT POINT
    !  *  G()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          GRADIENT OF THE OBJECTIVE FUNCTION AT CURRENT POINT
    !     DX() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          THE SEARCH DIRECTION
    !     N    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS X,G,DX
    !     F()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTAINING
    !          THE VECTOR OF RESIDUALS AT CURRENT POINT
    !  *  V1() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTAINING
    !          THE DERIVATIVE OF F(XOLD+ALPHA*DX) WITH RESPECT TO ALPHA
    !          AT ALPHA=0
    !     M    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAYS F,V1,
    !          FNEW AND V2
    !     PHI  REAL SCALAR CONTAINING THE VALUE OF THE OBJECTIVE
    !          FUNCTION AT THE CURRENT POINT
    !     FFUNC SUBROUTINE NAME USED TO EVALUATE THE VECTOR OF
    !           RESIDUALS AT A CERTAIN POINT
    !     KOD  INTEGER SCALAR CONTAINING A CODE THAT SAYS HOW THE
    !           SEARCH DIRECTION DX HAS BEEN COMPUTED
    !           = 1 IF GAUSS-NEWTON DIRECTION
    !           =-1 IF SUBSPACE DIRECTION
    !           =-2 IF NEWTON DIRECTION
    !     PRANK  INTEGER SCALAR CONTAINING PSEUDO RANK OF MATRIX C
    !     BL() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          THE LOWER BOUNDS OF THE UNKNOWNS
    !     BU() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N CONTAINING
    !          THE UPPER BOUNDS OF THE UNKNOWNS
    !     NRACT INTEGER SCALAR CONTAINING THE NO. OF ACTIVE CONSTRAINTS
    !     ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !           CONTAINING A CODE WHICH INDICATES WHETHER AN UNKNOWN IS
    !           ACTIVE OR NOT.
    !           ASET(I) = 0 WHEN X(I) IS FREE
    !                   =+1 WHEN X(I) EQUALS BL(I)
    !                   =-1 WHEN X(I) EQUALS BU(I)

    !     ON RETURN

    !     X()  THE NEW POINT XNEW IN  (1)
    !     F()  THE VALUE OF THE RESIDUALS AT THE NEW POINT
    !     PHI  THE VALUE OF THE OBJECTIVE AT THE NEW POINT
    !     ASET() MODIFIED IN ONE ELEMENT IF A BOUND IS CROSSED
    !     NRACT  THE CURRENT NO. OF ACTIVE CONSTRAINTS
    !     EVAL INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !          DONE INSIDE THE LINESEARCH ROUTINE
    !     ALPHA REAL SCALAR CONTAINING THE STEPSIZE
    !     ALPHUP REAL SCALAR CONTAINING THE UPPER BOUND OF STEPZISE
    !     XDIFF REAL SCALAR CONTAINING THE EUCLIDEAN DISTANCE BETWEEN
    !          THE TWO POINTS XOLD AND XNEW  FROM  (1)
    !     ERROR INTEGER SCALAR CONTAINING
    !           =-1 IF DX IS NOT A DESCENT DIRECTION
    !           =-2 IF XDIFF < SQRT(SRELPR) OR IF ALPHA<= A LOWER BOUND
    !           < -10 AS A USER STOP INDICATOR
    !           = 0 OTHERWISE

    !     WORKING AREAS

    !     FNEW(),V2()  REAL SINGLY SUBSCRIPTED ARRAYS OF DIMENSION M
    !     GMOD()       REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N

    REAL (dp)  :: gmod(n)

    !     COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !     2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !     REPRESENT TIME STEP K-2 RESP. K-1
    !     THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr
    !*****************
    ! COMMON /negdir/ ifree,indic
    !******************

    !     INTERNAL VARIABLES

    INTEGER   :: i,ctrl,iev
    REAL (dp) :: alplow,phikp1,xel,dummy(1,1),magfy,dphize

    !   IF ABS(KOD)=2
    !     THEN TAKE AN UNDAMPED STEP
    !     ELSE USE LINEUC TO COMPUTE STEPSIZE ALPHA

    !   DETERMINE LOWER AND UPPER BOUND OF STEP LENGTH AND
    !   SUGGEST ALPHA AS STEPSIZE

    eval=0
    alphup=3.0_dp
    !************************
    IF(ifree > 0 .AND. indic < 0) alphup=10.0_dp*alphup
    !************************
    IF(bnd /= 0)THEN
       !      write(10,*) 'ASET: ',(aset(i),i=1,n)
       !      write(10,*) (dx(i),i=1,n)
       DO  i=1,n
          IF(aset(i) /= 0) CYCLE
          IF(dx(i) <= zero) GO TO 10
          IF(alphup*dx(i) < bu(i)-x(i)) CYCLE
          alphup=(bu(i)-x(i))/dx(i)
          CYCLE

10        IF(-alphup*dx(i) <= x(i)-bl(i)) CYCLE
          alphup=-(x(i)-bl(i))/dx(i)
       END DO
    END IF
    alplow=alphup/3000.0_dp
    IF(ABS(kod) == 2) GO TO 40
    magfy=3.0_dp
    IF(prank < rngkm1) magfy=6.0_dp
    alpha=MIN(one,magfy*alfkm1,alphup)
    !      write(10,*) 'In STEPUC: ALPHA,ALFKM1,ALPHUP='
    !      write(10,*) alpha,alfkm1,alphup
    !**************************
    IF(ifree > 0 .AND. indic < 0) alpha=MIN(10.0_dp,alphup)
    !*************************
    alpha=MAX(alpha,alplow)*2.0_dp

    !     COMPUTE THE DERIVATIVE OF PHI(X+ALPHA*DX) AT ALPHA=0

    dphize=DOT_PRODUCT( g(1:n), dx(1:n) )

30  alpha=alpha*0.5_dp

    !     COMPUTE A STEP LENGTH

    CALL lineuc(x,gmod,dx,f,v1,m,n,alpha,phi,dphize,alplow,ffunc,  &
         alphup,phikp1,xdiff,iev,error)
    
    IF(error < -10) RETURN
    eval=eval + iev
    IF(-3 == error) GO TO 30
    phi=phikp1
    GO TO 60

    !     TAKE AN UNDAMPED STEP

40  alpha=MIN(one,alphup)
    error=0
    xdiff=zero
    DO  i=1,n
       xel=x(i)
       x(i)=xel + alpha*dx(i)
       xdiff=xdiff + (x(i)-xel)**2
    END DO
    xdiff=SQRT(xdiff)
    

    !     COMPUTE THE VECTOR OF RESIDUALS AT THE NEW POINT

    ctrl=1
    CALL ffunc(x,n,f,m,ctrl,dummy,1)
    eval=1

    !     INDICATE A LARGER VALUE OF THE OBJECTIVE IF THE RESIDUALS ARE
    !     UNCOMPUTABLE  (I.E. CTRL=-1 )

    IF(-1 == ctrl) phi=phi*2.0_dp
    IF(ctrl == 1)  phi=0.5_dp*dnrm2(m,f,1)**2
    IF(ctrl < -10) error=ctrl

    !     CHECK IF ANY BOUND IS CROSSED. IF SO, MAKE THE
    !     CORRESPONDING CONSTRAINT ACTIVE

    !      DO 80 I=1,N
    !           IF(ASET(I).NE.0) GOTO 80
    !           IF(X(I).GE.(BL(I)+C1)) GOTO 70
    !           ASET(I)=1
    !           X(I)=BL(I)
    !           NRACT=NRACT+1
    !           GOTO 80
    !   70      CONTINUE
    !           IF(X(I).LE.(BU(I)-C1)) GOTO 80
    !           ASET(I)=-1
    !           X(I)=BU(I)
    !           NRACT=NRACT+1
    !   80 CONTINUE
      
60  RETURN
  END SUBROUTINE stepuc


  !EVREUC

  SUBROUTINE evreuc(xkm1,n,m,k,ffunc,funcev,alpha,alphup,d1sqs,beta,fsum,  &
       dxnorm,kod,prank,phi,error,x,f,restar,aset,bl,bu,nract)

    REAL (dp), INTENT(IN)      :: xkm1(:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: m
    !tpk INTEGER, INTENT(OUT)       :: k
    INTEGER, INTENT(INOUT)       :: k  !tpk
    !tpk INTEGER, INTENT(OUT)       :: funcev
    INTEGER, INTENT(INOUT)       :: funcev  !tpk
    REAL (dp), INTENT(IN OUT)  :: alpha
    REAL (dp), INTENT(IN OUT)  :: alphup
    REAL (dp), INTENT(IN)      :: d1sqs
    REAL (dp), INTENT(IN)      :: beta
    REAL (dp), INTENT(IN)      :: fsum
    REAL (dp), INTENT(IN)      :: dxnorm
    INTEGER, INTENT(IN)        :: kod
    INTEGER, INTENT(IN)        :: prank
    REAL (dp), INTENT(IN OUT)  :: phi
    INTEGER, INTENT(IN OUT)    :: error
    REAL (dp), INTENT(IN OUT)  :: x(:)
    REAL (dp), INTENT(OUT)     :: f(:)
    LOGICAL, INTENT(OUT)       :: restar
    INTEGER, INTENT(IN OUT)    :: aset(:)
    REAL (dp), INTENT(IN)      :: bl(:)
    REAL (dp), INTENT(IN)      :: bu(:)
    !tpk INTEGER, INTENT(OUT)       :: nract
    INTEGER, INTENT(INOUT)       :: nract !tpk

    EXTERNAL  ffunc

    
    ! COMMON /back/ icount,lattry,itotal,philat,bestpg,irank
    !*******************
    ! COMMON /negdir/ ifree,indic
    !*******************

    !   CHECK IF PREVIOUS STEP WAS A COMPLETE FAILURE
    !   I.E.  IF THE VALUE OF THE OBJECTIVE FUNCTION INCREASED
    !         IN CURRENT STEP WE RESTART AT PREVIOUS POINT
    !   IF NO RESTART IS DONE WE UPDATE CERTAIN COMMON-VARIABLES
    !   THAT HOLD INFORMATION OF THE TWO LATEST POINTS IN THE
    !   ITERATION      K+1    K
    !                 X   := X +ALPHA*DX

    !   ON ENTRY

    !   XKM1() REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION CONTAINING
    !          THE PREVIOUS POINT (TIME STEP K)
    !   N    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAY XKM1,X
    !   M    INTEGER SCALAR CONTAINING THE LENGTH OF THE ARRAY F
    !   K    INTEGER SCALAR CONTAINING ITERATION NUMBER
    !   FFUNC SUBROUTINE NAME USED TO EVALUATE THE VECTOR OF
    !         RESIDUALS AT A CERTAIN POINT
    !   FUNCEV INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !          DONE SO FAR
    !   ALPHA REAL SCALAR CONTAINING THE STEPSIZE ALPHA USED IN
    !         PREVIOUS STEP
    !   ALPHUP REAL SCALAR CONTAINING THE UPPER BOUND OF STEP LENGTH
    !          USED IN PREVIOUS STEP
    !   D1SQS REAL SCALAR CONTAINING PREDICTED REDUCTION IN THE
    !         OBJECTIVE FUNCTION FOR LATEST STEP
    !   BETA REAL SCALAR CONTAINING THE NORM OF ORTHOGONAL PROJECTION
    !        OF F ONTO COLUMN SPACE OF THE JACOBIAN
    !   FSUM REAL SCALAR CONTAINIG THE SUM OF SQUARS AT PREVIOUS POINT
    !        (TIME STEP K)
    !   DXNORM REAL SCALAR CONTAINING THE NORM OF THE SEARCH
    !          DIRECTION (DX)
    !   KOD  INTEGER SCALAR CONTAINING A CODE THAT SAYS HOW THE SEARCH
    !        DIRECTION DX WAS COMPUTED
    !        = 1 GAUSS-NEWTON DIRECTION
    !        =-1 SUBSPACE DIRECTION
    !        =-2 IF NEWTON DIRECTION
    !   PRANK INTEGER SCALAR CONTAINING A CODE DEPENDING OF KOD
    !         IF ABS(KOD)=1
    !           THEN PRANK EQUALS PSEUDO RANK USED TO COMPUTE DX
    !           ELSE PRANK EQUALS -N
    !   PHI  REAL SCALAR CONTAINING THE VALUE OF OBJECTIVE FUNCTION
    !          AT THE CURRENT POINT (TIME STEP K+1)
    !   ERROR INTEGER SCALAR CONTAINING
    !         -1 IF DX WAS NOT A DESCENT DIRECTION
    !         -2 IF ALPHA*II DX II <SQRT(SRELPR)
    !         -3 IF NOT POSITIVE DEFINITE MATRIX FROM SECOND DERIVATIVE
    !            SYSTEM
    !         -4 IF THE ALGORITHM WANTS TO USE SECOND DERIVATIVES BUT
    !            IS NOT ALLOWED TO DO THAT
    !          < -10 AS A USER STOP INDICATOR
    !          0 OTHERWISE

    !   ON RETURN

    !   RESTAR LOGICAL SCALAR SET TO TRUE IF RESTART IS RECOMMENDABLE
    !          SET TO FALSE OTHERWISE
    !   IF RESTART IS NOT RECOMMENDABLE THE FOLLOWING ARE CHANGED
    !   AND VARIABLES HELD IN COMMON-BLOCK /PREUNC/ REPRESENTING THE
    !   PAST ARE UPDATED

    !   K     IS INCREASED BY 1
    !   ERROR IS SET TO -5 IF THE OBJECTIVE IS INCREASED IN THE
    !         LATEST STEP

    !   IF RESTART IS RECOMMENDABLE THE FOLLOWING ARE CHANGED
    !   FUNCEV IS INCREASED BY 1
    !   X()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION N TAKES
    !        THE VALUE OF THE PREVIOUS POINT (XKM1)
    !   F()  REAL SINGLY SUBSCRIPTED ARRAY OF DIMENSION M CONTIANS
    !        THE VECTOR OF RESIDUALS AT THE PREVIOUS POINT


    !   COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !   2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !   REPRESENT TIME STEP K-2 RESP. K-1
    !   THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC

    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1

    !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    !     INTERNAL VARIABLES

    INTEGER   :: i,ctrl
    REAL (dp) :: c2,dummy(1,1),dsqrel

    restar=.false.
    dsqrel=SQRT(srelpr)

    !     CALL RESTART ADVISOR ROUTINE

    CALL reavuc(error,alpha,alphup,fsum,phi,restar)
    IF(restar) GO TO 20

    !     CHECK IF A NEWTON STEP HAS FAILED

    !      IF(-5.NE.ERROR) GOTO 10
    !      DO 5 I=1,N
    !           X(I)=XKM1(I)
    !    5 CONTINUE
    !      CTRL=-1
    !      CALL FFUNC(X,N,F,M,CTRL,DUMMY,1)
    !      IF(CTRL.LT.-10) GOTO 40
    !      FUNCEV=FUNCEV+1
    !**************************
    IF(-5 /= error) GO TO 10
    !       write(10,*) '!!!!!!!!!Undamped Newton  does not work: Try GN'
    ifree=5
    indic=3
    k=k+1
    GO TO 25
    !**************************

    !     NO RESTART. UPDATE HISTORY VARIABLES

10  icount=0
    lattry=prank
    betkm2=betkm1
    d1km2=d1km1
    fsqkm2=fsqkm1
    dxnkm2=dxnkm1
    alfkm2=alfkm1
    aupkm2=aupkm1
    rngkm2=rngkm1
    kodkm2=kodkm1
    betkm1=beta
    d1km1=d1sqs
    fsqkm1=fsum
    dxnkm1=dxnorm
    alfkm1=alpha
    aupkm1=alphup
    rngkm1=prank
    kodkm1=kod

    !     CHECK IF ANY BOUND IS CROSSED. IF SO, MAKE THE
    !     CORRESPONDING CONSTRAINT ACTIVE

    DO  i=1,n
       IF(aset(i) /= 0) CYCLE
       !   c2=min(abs(BL(i)),one)
       c2=0.0
       IF(x(i) >= (bl(i)+c2+dsqrel)) GO TO 70
       !           IF(X(I).GE.(BL(I)+C1)) GOTO 70
       aset(i)=1
       x(i)=bl(i)
       nract=nract+1
       CYCLE

       !   c2=min(abs(BU(i)),one)

70     c2=0.0
       IF(x(i) <= (bu(i)-c2-dsqrel)) CYCLE
       !           IF(X(I).LE.(BU(I)-C1)) GOTO 80
       aset(i)=-1
       x(i)=bu(i)
       nract=nract+1
    END DO

    !     INCREASE ITERATION COUNT

    k=k+1
    RETURN

    !      write(10,*) 'In EVREUC: Restart!'
20  itotal=itotal+1
    icount=icount+1
    lattry=prank
    philat=phi

    !     RESTART FROM PREVIOUS POINT

    !****************************
25  error=0
    !****************************
    x(1:n)=xkm1(1:n)
    ctrl=-1

    !     EVALUATE THE FUNCTIONS AT PREVIOUS POINT

    CALL ffunc(x,n,f,m,ctrl,dummy,1)
    IF(ctrl < -10) GO TO 40
    alpha=alfkm1
    alphup=aupkm1
    phi=0.5_dp*dnrm2(m,f,1)**2
    funcev=funcev+1
    RETURN

    !     A USER STOP INDICATION IS DETECTED

40  error=ctrl

    RETURN
  END SUBROUTINE evreuc


  !REAVUC

  SUBROUTINE reavuc(error,alpha,alphup,fsum,phi,restar)

    INTEGER, INTENT(IN OUT)  :: error
    REAL (dp), INTENT(IN)    :: alpha
    REAL (dp), INTENT(IN)    :: alphup
    REAL (dp), INTENT(IN)    :: fsum
    REAL (dp), INTENT(IN)    :: phi
    LOGICAL, INTENT(IN OUT)  :: restar

    ! COMMON /back/ icount,lattry,itotal,philat,bestpg,irank

    !   THIS IS THE RESTART ADVISOR ROUTINE FOR UNCONSTRAINED PROBLEMS

    !   ON ENTRY

    !   ERROR INTEGER SCALAR CONTAINING
    !         -1 IF NO DESCENT DIRECTION
    !         -2 IF LATEST ALPHA*II DX II < SQRT(SRELPR) OR
    !            IF ALPHA <= A LOWER BOUND ON THE STEPLENGTH
    !         -3 IF NOT POS. DEF. MATRIX FROM 2:ND ORDER METHOD
    !         -4 IF NOT ALLOWED TO USE 2:ND DERIVATIVES
    !          0 OTHERWISE
    !   ALPHA REAL SCALAR CONTAINING THE LATEST STEP LENGTH
    !   ALPHUP REAL SCALAR CONTAINING THE LATEST UPPER BOUND
    !   FOR THE STEPLENGTH
    !   FSUM  REAL SCALAR CONTAINING THE SUM OF SQUARES
    !         AT THE PREVIOUS POINT
    !   PHI   REAL SCALAR CONTAINING THE VALUE OF THE OBJECTIVE
    !         FUNCTION AT THE CURRENT POINT
    !   RESTAR LOGICAL SCALAR EQUAL TO FALSE
    !   LATTRY  INTEGER SCALAR STORED IN COMMON /BACK/ EQUAL TO ZERO
    !           IF NO CHECK SHALL BE DONE

    !   ON RETURN

    !   ERROR   IS SET =-5 IF THE OBJECTIVE INCREASES
    !   RESTAR IS SET = TRUE IF THE LATEST STEPLENGTH IS TOO SHORT

    !      write(10,*) 'In REAVUC: alpha,alphup=',alpha,alphup
    IF(lattry == 0 .AND. bestpg > zero) RETURN
    IF(lattry == 0 .AND. bestpg <= zero) GO TO 5
    IF(error == -1 .OR. error <= -3) RETURN
    IF(-2 == error) GO TO 10
    IF(phi < 0.5_dp*fsum) GO TO 10

5   error=-5
    RETURN

10  IF(alpha <= alphup/3000.0_dp) restar=.true.

    RETURN
  END SUBROUTINE reavuc


  !OUTUC

  SUBROUTINE outuc (iprint, k, restar, UNIT, phi, gnorm, eval, aset,  &
       n, nract, speed)

    INTEGER, INTENT(IN)     :: iprint
    INTEGER, INTENT(IN)     :: k
    LOGICAL, INTENT(IN)     :: restar
    INTEGER, INTENT(IN)     :: UNIT
    REAL (dp), INTENT(IN)   :: phi
    REAL (dp), INTENT(IN)   :: gnorm
    INTEGER, INTENT(IN)     :: eval
    INTEGER, INTENT(IN)     :: aset(:)
    INTEGER, INTENT(IN)     :: n
    INTEGER, INTENT(IN)     :: nract
    REAL (dp), INTENT(OUT)  :: speed

    !
    !     IF IPRINT > 0
    !       THEN WRITE SOME INFO. EVERY IPRINT STEP ON UNIT
    !       ELSE NO PRINTING IS DONE
    !
    !     ON ENTRY
    !
    !     IPRINT INTEGER SCALAR CONTAINING ZERO OR A POSITIVE VALUE
    !     K    INTEGER SCALAR CONTAINING ITERATION COUNT
    !     RESTAR LOGICAL SCALAR =TRUE IF RESTART IS DONE
    !                           =FALSE IF NO RESTART IS DONE
    !     UNIT INTEGER SCALAR CONTAINING LOGICAL DEVICE NUMBER WHERE
    !          PRINTING SHALL BE DONE
    !     EVAL INTEGER SCALAR CONTAINING NO. OF FUNCTION EVALUATIONS
    !          DONE INSIDE THE LINESEARCH ROUTINE FOR LATEST STEP
    !     PHI  REAL SCALAR CONTAINING THE VALUE OF THE OBJECTIVE
    !          FUNCTION AT CURRENT POINT
    !     GNORM REAL SCALAR CONTAINING THE NORM OF THE GRADIENT OF
    !           THE OBJECTIVE FUNCTION AT THE PREVIOUS POINT DIVIDED BY
    !           THE LONGEST COLUMN OF THE JACOBIAN
    !
    !     N      INTEGER SCALAR CONTAINING THE NO. OF UNKNOWNS
    !     NRACT INTEGER SCALAR CONTAINING THE NO. OF ACTIVE CONSTRAINTS
    !     ASET() INTEGER SINGLY SUBSCRIPTED ARRAY OF DIMENSION N
    !           CONTAINING A CODE WHICH INDICATES WHETHER AN UNKNOWN IS
    !           ACTIVE OR NOT.
    !           ASET(I) = 0 WHEN X(I) IS FREE
    !                   =+1 WHEN X(I) EQUALS BL(I)
    !                   =-1 WHEN X(I) EQUALS BU(I)
    !     ON RETURN:
    !
    !     SPEED   REAL SCALAR CONTAINING AN ESTIMATE OF THE LINEAR
    !             CONVERGENCE FACTOR
    !
    !     COMMON VARIABLES CONTAINING INFORMATION CONCERNINING PREVIOUS
    !     2 POINTS.THE SUFFICES KM2 RESP. KM1 IN THE NAMES OF THE VARIABLES
    !     REPRESENT TIME STEP K-2 RESP. K-1
    !     THESE VARIABLES ARE UPDATED ONLY INSIDE ROUTINE EVREUC
    !
    ! COMMON /preunc/ betkm2,d1km2,fsqkm2,dxnkm2,alfkm2,aupkm2,rngkm2,  &
    !    kodkm2,betkm1,d1km1,fsqkm1,dxnkm1,alfkm1,aupkm1,rngkm1,kodkm1
    ! COMMON /negdir/ ifree,indic
    !
    !     INTERNAL VARIABLES
    !
    INTEGER   :: itno
    REAL (dp) :: reduc,ac,pred

    IF (restar) RETURN
    IF (ifree == 5) RETURN
    speed=zero
    itno=k-1
    IF (itno /= 0 .AND. betkm2 /= zero) speed=betkm1/betkm2
    IF (iprint <= 0) RETURN
    IF (itno/iprint*iprint /= itno) RETURN
    reduc=fsqkm1 - 2.0_dp*phi
    ac=MIN(one,aupkm1)
    pred=ac*(2.0_dp-ac)*d1km1
    IF (itno == 0) WRITE (UNIT,10)
    WRITE (UNIT,20) itno,fsqkm1,gnorm,dxnkm1,kodkm1,rngkm1,alfkm1,  &
         eval,pred,speed,phi,reduc
    IF (nract > 0) WRITE (UNIT,30) aset(1:n)
    RETURN

10  FORMAT (////t11, 'COLLECTED INFORMATION FOR ITERATION STEPS'//  &
         ' K    FSUM(K)     GNORM(K)    DXNORM    CODE  ALPHA ',  &
         '  EVAL    PREDICTED  CONV.SPEED   PHI(K+1)   REDUCTION',  &
         '  ACTIVE SET'/)
20  FORMAT (' ', i4, 3g12.3, 2I3, g12.3, i3, g13.3, 3g12.3)
30  FORMAT ('', (t60, 5I3))

  END SUBROUTINE outuc


  !SAVEUC

  SUBROUTINE saveuc(x,n,xkm1)

    REAL (dp), INTENT(IN)   :: x(:)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: xkm1(:)

    !     STORE THE ARRAY X OF DIMENSION N INTO THE ARRAY XKM1

    xkm1(1:n)=x(1:n)

    RETURN
  END SUBROUTINE saveuc


  !GRAD

  SUBROUTINE grad(a,m,n,f,g)

    ! Argument MDA has been removed.

    !   COMPUTE THE GRADIENT G(X) TO  0.5* II F(X) II**2       T
    !          WHERE  F(X) = (F1(X),F2(X),...............FM(X))
    !   THE MATRIX A  (M*N)  IS THE JACOBIAN OF  F(X)

    REAL (dp), INTENT(IN)   :: a(:,:)   ! a(mda,n)
    INTEGER, INTENT(IN)     :: m
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(IN)   :: f(:)
    REAL (dp), INTENT(OUT)  :: g(:)

    ! Internal variable
    INTEGER  :: i

    DO  i=1,n
       g(i)=DOT_PRODUCT( a(1:m,i), f(1:m) )
    END DO

    RETURN
  END SUBROUTINE grad


  !RELEPS

  SUBROUTINE releps(srelpr)

    REAL (dp), INTENT(OUT)  :: srelpr

    !   COMPUTE SRELPR = DOUBLE RELATIVE PRECISION FOR A BINARY
    !   MACHINE   I.E.
    !   DETERMINE THE SMALLEST POSITIVE NUMBER 0.5**K FOR WHICH
    !      (1.0+0.5**K) > 1.0  AND  (1.0+0.5**(K+1)) = 1.0
    !   WHERE K IS A POSITIVE INTEGER

    srelpr = EPSILON(one)

    RETURN
  END SUBROUTINE releps


  !LINEUC

  SUBROUTINE lineuc(xold,g,p,f,v1,m,n,alpha,phizer,dphize,alflow,  &
       ffunc,alfupp,phialf,xdiff,eval,EXIT)

    ! Arguments FNEW & V2 have been removed.

    REAL (dp), INTENT(IN OUT)  :: xold(:)
    REAL (dp), INTENT(IN OUT)  :: g(:)
    REAL (dp), INTENT(IN)      :: p(:)
    REAL (dp), INTENT(IN OUT)  :: f(:)
    REAL (dp), INTENT(IN OUT)  :: v1(:)
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: n
    !tpk REAL (dp), INTENT(OUT)     :: alpha
    REAL (dp), INTENT(INOUT)     :: alpha !tpk
    REAL (dp), INTENT(IN)      :: phizer
    REAL (dp), INTENT(IN)      :: dphize
    REAL (dp), INTENT(IN)      :: alflow
    REAL (dp), INTENT(IN)      :: alfupp
    REAL (dp), INTENT(OUT)     :: phialf
    REAL (dp), INTENT(OUT)     :: xdiff
    INTEGER, INTENT(OUT)       :: eval
    INTEGER, INTENT(OUT)       :: EXIT

    EXTERNAL ffunc

    
    !     COMMON VARIABLES CONTAINING MACHINE DEPENDENT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    !   THIS IS A LINE SEARCH ROUTINE FOR UNCONSTRAINED LEAST SQUARES PROBLEMS

    !   COMPUTE THE STEPLENGTH ALPHA FOR THE ITERATION
    !   XNEW:=XOLD+ALPHA*P
    !   WHERE   XOLD IS THE CURRENT POINT
    !           P    IS THE SEARCH DIRECTION
    !   ALPHA IS CLOSE TO THE SOLUTION OF THE PROBLEM
    !           MINIMIZE  PHI(ALPHA)
    !   WITH THE RESTRICTION  0<ALFLOW<=ALPHA<=ALFUPP
    !   HOWEVER, IF WE ARE FORCED TO TAKE A PURE GOLDSTEIN-ARMIJO STEP
    !   THE STEPLENGTH CAN BE SLIGHTLY LOWER THAN ALFLOW

    !   PHI(ALPHA)=0.5*(IIF(XOLD+ALPHA*P)II**2)
    !   F(X)= (F1(X),F2(X),.........,FM(X)) TRANSPOSE

    !   ON ENTRY:

    !   XOLD(N)      REAL ARRAY OF LENGTH N-THE CURRENT POINT
    !                (CHANGED ON RETURN)
    !   G(N)         REAL ARRAY OF LENGTH N-GRADIENT OF PHI AT ALPHA=0
    !                (CHANGED ON RETURN)
    !   P(N)         REAL ARRAY OF LENGTH N-SEARCH DIRECTION
    !   F(M)         REAL ARRAY OF LENGTH M -F(1)....F(M) CONTAIN THE
    !                VALUE OF F(X) AT THE CURRENT POINT.
    !                (CHANGED ON RETURN)
    !   V1(M)        REAL ARRAY 0F LENGTH M: THE VECTOR C*P WHERE
    !                C IS THE JACOBIAN AND P THE SEARCH DIRECTION
    !   M            INTEGER-NO. OF FUNCTIONS IN F(X)=F1(X).......,FM(X)
    !   N            INTEGER-NO. OF UNKNOWNS
    !   ALPHA        REAL-A FIRST GUESS OF STEPLENGTH
    !   PHIZER       REAL-PHI(ALPHA) AT ALPHA=0
    !   DPHIZE       REAL-THE DERIVATIVE OF PHI(ALPHA) AT ALPHA=0
    !   ALFLOW       REAL-FIX LOWER BOUND OF THE STEPLENGTH
    !   FFUNC        SUBROUTINE NAME-THAT MUST BE DECLEARED EXTERNAL IN
    !                THE CALLING PROGRAM. A ROUTINE NAMED FFUNC IS USED
    !                TO EVALUATE THE FUNCTION F(X) AND/OR THE JACOBIAN AT
    !                A CERTAIN POINT X AND SHOULD BE WRITTEN AS FOLLOWS

    !                SUBROUTINE FFUNC(X,N,F,M,CTRL,C,MDC)
    !                INTEGER N,M,CTRL,MDC
    !                REAL X(N),F(M),C(MDC,N)
    !                -----------------------
    !                CTRL CAN HAVE 3 DIFFERENT VALUES ON ENTRY
    !       CTRL= 1  MEANS EVALUATE THE FUNCTIONS AT THE POINT X AND
    !                RETURN THIS VECTOR IN THE ARRAY F IF THE FUNCTIONS
    !                ARE COMPUTABLE.
    !                ON RETURN THE USER CAN INDICATE UNCOMPUTABILITY BY
    !                SETTING CTRL=-1
    !                DO NOT ALTER ARRAY X.
    !       CTRL=-1  MEANS EVALUATE THE FUNCTIONS AT THE POINT X AND
    !                RETURN THIS VECTOR IN THE ARRAY F IF THE FUNCTIONS
    !                ARE COMPUTABLE. DO NOT ALTER ARRAY X.
    !                POSSIBLE UNCOMPUTABILITY OF THE FUNCTIONS MUST BE
    !                INDICATED BY SETTING CTRL TO A VALUE <-10 ON RETURN
    !       CTRL= 2  MEANS CALCULATE THE JACOBIAN OF F(X) AT THE POINT X
    !                AND RETURN THIS MATRIX IN THE ARRAY C IF THE JACOBIAN
    !                IS SUPPLIED ANALYTICALLY.
    !                POSSIBLE UNCOMPUTABILITY OF THE JACOBIAN MUST BE
    !                INDICATED BY SETTING CTRL TO A VALUE <-10 ON RETURN.
    !                IF THE USER WANTS THE JACOBIAN BEING COMPUTED
    !                NUMERICALLY THAT SHOULD BE INDICATED BY SETTING
    !                CTRL=0 ON RETURN.
    !                DO NOT ALTER ARRAYS X AND F.
    !                ------------------------------
    !                RETURN
    !                END
    !   ALFUPP       REAL-FIX UPPER BOUND OF THE STEPLENGTH

    !   ON RETURN:     AND EXIT = 0

    !   XOLD(N)      THE NEW POINT
    !   F(M)         THE VALUE OF F(X) AT THE NEW POINT
    !   ALPHA        THE COMPUTED STEPLENGTH
    !   PHIALF       REAL-PHI(ALPHA) AT THE NEW POINT
    !   XDIFF        REAL-THE 2-NORM  II XOLD-XNEW II
    !   EXIT         INTEGER-  =-1 IF THE SEARCH DIRECTION IS NOT A
    !                              DESCENT DIRECTION.
    !                          =-2 IF ALPHA*2NORM(P) < SQRT(SRELPR)
    !                              OR IF ALPHA<=ALFLOW
    !                          =-3 IF THE FIRST GUESS OF STEPLENGTH
    !                              MAKES THE FUNCTIONS UNCOMPUTABLE AT X
    !                          < -10 AS A USER STOP INDICATOR
    !                          = 0 OTHERWISE
    !   EVAL         INTEGER- =NO. OF FUNCTION EVALUATIONS DONE
    !                INSIDE THIS ROUTINE

    !   WORKING AREAS:
    !   FNEW(M)      REAL ARRAY OF LENGTH M
    !   V2(M)        REAL ARRAY OF LENGTH M

    !   THE LINEUC-PACKAGE USES THE BLAS ROUTINE SNRM2

    !   INTERNAL VARIABLES

    INTEGER   :: ctrl,i,k
    REAL (dp) :: eta,tau,gamma,alfmax,alfmin,alfk,alfkm1,alfkm2,alfkp1,  &
         pmax,beta,diff,pbeta,phik,phikm1,phikm2,pk,xmin
    REAL (dp) :: fnew(m),v2(m)
    LOGICAL   :: reduce

    k=0
    xdiff=zero
    alfkm1=zero
    phikm1=phizer

    !     SET VALUES OF THE CONSTANTS ETA,TAU AND GAMMA.
    !     COMPUTE ALFMIN,ALFMAX,ALFK AND PMAX= THE EUCLIDEAN NORM OF P

    eta=0.3_dp
    tau=0.25_dp
    gamma=0.4_dp
    alfmax=alfupp
    alfmin=alflow
    alfk=MIN(alpha,alfmax)
    alfk=MAX(alfk,alfmin)
    pmax=dnrm2(n,p,1)

    !     P MUST BE A DESCENT DIRECTION

    EXIT=-1
    IF(dphize >= zero) GO TO 1020
    EXIT=0

    !     COMPUTE PHIK=PHI(ALF0) AND TEST UNCOMPUTABILITY AT XOLD

    ctrl=1
    CALL fsumsq(xold,p,n,alfk,g,fnew,m,ffunc,ctrl,phik)
    !      write(10,*) 'In LINEUC no.1: PHIK,alf0= ',phik,alfk
    k=k+1
    eval=k
    IF(-1 == ctrl) EXIT=-3
    IF(alfk <= zero) EXIT=-2
    IF(ctrl < -10) EXIT=ctrl
    IF(EXIT < 0) RETURN

    !     COMPUTE THE VECTOR V2 SO THAT A ONE DIMENSIONAL
    !     MINIMIZATION IN R(M) CAN BE DONE

    CALL linuc2(m,v1,fnew,f,alfk,v2)
    diff=phizer-phik

    !     SET XMIN = THE BEST OF THE POINTS 0 AND ALF0

    IF(diff >= zero) xmin=alfk
    IF(diff < zero) xmin=zero

    !     MINIMIZE IN R(M). USE TWO POINTS : 0 AND ALF0
    !     NEW SUGGESTION OF STEPLENGTH IS ALFKP1
    !     PK IS THE VALUE OF THE APPROXIMATING FUNCTION AT ALFKP1

    CALL minrm(f,v1,v2,m,alfmin,alfmax,xmin,alfkp1,pk,beta,pbeta)

    !     POSSIBLY THE OTHER ROOT IS CHOSEN

    IF(alfkp1 == beta) GO TO 20
    IF(pk <= pbeta) GO TO 20
    IF(beta > alfk) GO TO 20
    alfkp1=beta
    pk=pbeta

20  alfkm1=zero
    phikm1=phizer
    CALL update(alfkm2,phikm2,alfkm1,phikm1,alfk,phik,alfkp1)

    !     TEST TERMINATION CONDITION AT ALPHA = ALF0

    IF(.NOT.(-diff <= tau*dphize*alfkm1 .OR. phikm1 < gamma*phizer )) GO TO 100

    !     TERMINATION CONDITION SATISFIED AT ALPHA = ALF0

30  diff=phizer-phik

    !     CHECK IF ESSENTIAL REDUCTION IS LIKELY

    CALL reduc(alfkm1,phikm1,alfk,pk,diff,eta,xold,p,f,  &
         g,fnew,m,n,ffunc,k,phik,reduce)
    IF(k < -10) GO TO 1030
    IF(.NOT.reduce) GO TO 1000

    !     THE VALUE OF THE OBJECTIVE FUNCTION CAN MOST LIKELY BE REDUCED

    !     MINIMIZE IN R(N). USE THREE POINTS: ALFKM2,ALFKM1,ALFK
    !     NEW SUGGESTION OF THE STEPLENGTH IS ALFKP1
    !     PK IS THE VALUE OF THE APPROXIMATING FUNCTION AT ALFKP1

    CALL minrn(alfk,phik,alfkm1,phikm1,alfkm2,phikm2,alfmin,  &
         alfmax,pmax,srelpr,alfkp1,pk)
    CALL update(alfkm2,phikm2,alfkm1,phikm1,alfk,phik,alfkp1)
    GO TO 30

    !     TERMINATION CONDITION NOT SATISFIED AT ALPHA = ALF0

    !     COMPUTE PHIK=PHI(ALF1)

100 ctrl=-1
    CALL fsumsq(xold,p,n,alfk,g,fnew,m,ffunc,ctrl,phik)
    !      write(10,*) 'In LINEUC no.2: PHIK,alf1= ',phik,alfk
    IF(ctrl < -10) k=ctrl
    IF(k < 0) GO TO 1030

    !     TEST TERMINATION CONDITION AT ALPHA = ALF1

    diff=phizer-phik
    IF(.NOT.(-diff <= tau*dphize*alfk .OR. phik < gamma*phizer)) GO TO 200


    !     TERMINATION CONDITION SATISFIED AT ALPHA = ALF1

    !     CHECK IF ALF0 IS SOMEWHAT GOOD

    IF(phizer > phikm1) GO TO 120

    !     SINCE PHI(0) <= PHI(ALF0), ALF0 IS NO GOOD GUESS AND WE TRY
    !     AN OTHER MINIMIZATION IN R(M). USE TWO POINTS: 0 AND ALF1.
    !     THE NEW SUGGESTION OF STEPLENGTH IS ALFKP1.
    !     PK IS THE VALUE OF THE APPROXIMATING FUNCTION AT ALFKP1

    xmin=alfk
    CALL linuc2(m,v1,fnew,f,alfk,v2)
    CALL minrm(f,v1,v2,m,alfmin,alfmax,xmin,alfkp1,pk,beta,pbeta)

    !     POSSIBLY THE OTHER ROOT IS CHOSEN

    IF(alfkp1 == beta) GO TO 115
    IF(pk <= pbeta) GO TO 115
    IF(beta > alfk) GO TO 115
    alfkp1=beta
    pk=pbeta

115 alfkm1=zero
    phikm1=phizer
    GO TO 130

    !     MINIMIZE IN R(N). USE THREE POINTS: 0,ALF0 AND ALF1.
    !     THE NEW SUGGESTION OF THE STEPLENGTH IS ALFKP1 .
    !     PK IS THE VALUE OF THE APPROXIMATING FUNCTION AT ALFKP1

120 CALL minrn(alfk,phik,alfkm1,phikm1,alfkm2,phikm2,alfmin,  &
         alfmax,pmax,srelpr,alfkp1,pk)
130 k=k+1
140 diff=phizer - phik
    CALL update(alfkm2,phikm2,alfkm1,phikm1,alfk,phik,alfkp1)

    !     CHECK IF ESSENTIAL REDUCTION IS LIKELY

    CALL reduc(alfkm1,phikm1,alfk,pk,diff,eta,xold,p,f,  &
         g,fnew,m,n,ffunc,k,phik,reduce)
    IF(k < -10) GO TO 1030
    IF(.NOT.reduce) GO TO 1000

    !     MINIMIZE IN R(N). USE THREE POINTS: ALFKM2,ALFKM1 AND ALFK.
    !     THE NEW SUGGESTION OF STEPLENGTH IS ALFKP1.
    !     PK IS THE VALUE OF THE APPROXIMATING FUNCTION AT ALFKP1

    CALL minrn(alfk,phik,alfkm1,phikm1,alfkm2,phikm2,alfmin,  &
         alfmax,pmax,srelpr,alfkp1,pk)
    GO TO 140

    !     TAKE A PURE GOLDSTEIN-ARMIJO STEP

200 k=k+1
    CALL gauc(xold,p,f,m,n,ffunc,k,alfmin,EXIT,g,fnew,phizer,dphize,  &
         alfk,phik,tau,pmax)!tpk ,srelpr)
    IF(k < -10) GO TO 1030

    !     COMPARE TWO EXPRESSIONS FOR THE DERIVATIVE OF PHI(XOLD+ALPHA*P)
    !     AT ALPHA=0. IF INCONSISTENCY IS DETECTED K IS SET=-1 ON RETURN.

    IF(-2 == EXIT) CALL chder(dphize,phizer,xold,p,m,n,ffunc,k,EXIT,alfk,phik)
    IF(k < -10) GO TO 1030
    alfkm1=alfk
    phikm1=phik

    !     COMPUTE THE NEW POINT AND THE DIFFERENCE II XOLD-XNEW II

1000 DO  i=1,n
       g(i)=xold(i)
       xold(i)=xold(i) + alfkm1*p(i)
       g(i)=g(i)-xold(i)
    END DO
    xdiff=dnrm2(n,g,1)

1020 alpha=alfkm1
    phialf=phikm1
    eval=k
    RETURN

    !     A USER STOP INDICATION IS DETECTED

1030 EXIT=k

    RETURN
  END SUBROUTINE lineuc


  !LINUC2

  SUBROUTINE linuc2(m,v1,fnew,f,alfk,v2)

    !     FORM THE VECTOR V2 AS THE DIVIDED DIFFERENCE VECTOR OF SECOND ORDER.

    INTEGER, INTENT(IN)     :: m
    REAL (dp), INTENT(IN)   :: v1(:)
    REAL (dp), INTENT(IN)   :: fnew(:)
    REAL (dp), INTENT(IN)   :: f(:)
    REAL (dp), INTENT(IN)   :: alfk
    REAL (dp), INTENT(OUT)  :: v2(:)

    REAL (dp) :: f1,f2,f3,alf
    INTEGER   :: i

    alf=alfk
    DO  i=1,m
       f1=fnew(i)
       f2=f(i)
       f3=v1(i)
       v2(i)=((f1-f2)/alf - f3)/alf
    END DO

    RETURN
  END SUBROUTINE linuc2


  !REDUC

  SUBROUTINE reduc(alf,phialf,alfk,pk,diff,eta,xold,  &
       p,f,xnew,fnew,m,n,ffunc,k,phik,reduce)

    !   REDUCE IS SET TO TRUE IF ESSENTIAL REDUCTION OF THE
    !   OBJECTIVE FUNCTION IS LIKELY.
    !   OTHERWISE REDUCE IS SET TO FALSE

    REAL (dp), INTENT(OUT)     :: alf
    REAL (dp), INTENT(IN OUT)  :: phialf
    REAL (dp), INTENT(IN)      :: alfk
    REAL (dp), INTENT(IN)      :: pk
    REAL (dp), INTENT(IN OUT)  :: diff
    REAL (dp), INTENT(IN)      :: eta
    REAL (dp), INTENT(IN OUT)  :: xold(:)
    REAL (dp), INTENT(IN)      :: p(:)
    REAL (dp), INTENT(OUT)     :: f(:)
    REAL (dp), INTENT(OUT)     :: xnew(:)
    REAL (dp), INTENT(OUT)     :: fnew(:)
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(OUT)       :: k
    REAL (dp), INTENT(OUT)     :: phik
    LOGICAL, INTENT(OUT)       :: reduce

    EXTERNAL ffunc

    
    !     COMMON VARIABLES CONTAINING MACHINE DEPENDANT CONSTANTS
    !     SRELPR = SINGLE RELATIVE PRECISION

    ! COMMON /machin/ srelpr

    INTEGER    :: ctrl
    REAL (dp)  :: c1
    REAL (dp), SAVE  :: delta = 0.2_dp

    c1=m*SQRT(srelpr)
    IF(.NOT.(phialf-pk > eta*diff .OR. pk+c1 < delta*phialf)) GO TO 20
    f(1:m)=fnew(1:m)
    ctrl=-1
    CALL fsumsq(xold,p,n,alfk,xnew,fnew,m,ffunc,ctrl,phik)
    !      write(10,*) 'In REDUC: phik,alfk= ',phik,alfk
    IF(ctrl < -10) k=ctrl
    IF(k < 0) RETURN
    k=k+1
    reduce=.true.
    IF(phialf-phik > eta*diff .OR. phik < delta*phialf) RETURN

    !     TERMINATE BUT CHOOSE THE BEST POINT OUT OF ALF AND ALFK

    IF(phialf <= phik) GO TO 40
    alf=alfk
    phialf=phik

20  f(1:m)=fnew(1:m)

40  reduce=.false.

    RETURN
  END SUBROUTINE reduc


  !GAUC

  SUBROUTINE gauc(xold,p,f,m,n,ffunc,k,alfmin,EXIT,xnew,  &
       fnew,phi0,dphi0,u,phiu,tau,pmax)!tpk ,srelpr)

    !   THIS IS A ROUTINE FOR UNCONSTRAINED LEAST SQUARES PROBLEMS THAT HALVES
    !   THE VALUE OF U UNTIL A GOLDSTEIN-ARMIJO CONDITION IS SATISFIED OR UNTIL
    !   THE NORM OF THE SEARCH DIRECTION TIMES THE THE STEPLENGTH IS REDUCED BELOW
    !   SQRT(RELATIVE PRECISION) OR UNTIL THE STEPLENGTH IS SMALLER THAN ALFMIN

    !   PHI(ALPHA)=0.5*(IIF(XOLD+ALPHA*P)II**2)
    !   CHOOSE ALPHA=X SO THAT
    !                PHI(X)<=PHI(0)+TAU*X*DPHI(0)   (1)
    !   WE KNOW THAT PHI(U)>PHI(0)+TAU*U*DPHI(0)
    !   THE SIMPLEST WE CAN DO IS TO SET   U=U*0.5 AND
    !   TEST IF CONDITION  (1) IS SATISFIED FOR X=U

    REAL (dp), INTENT(IN OUT)  :: xold(:)
    REAL (dp), INTENT(IN)      :: p(:)
    REAL (dp), INTENT(OUT)     :: f(:)
    INTEGER, INTENT(IN)        :: m
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(OUT)       :: k
    REAL (dp), INTENT(IN OUT)  :: alfmin
    INTEGER, INTENT(OUT)       :: EXIT
    REAL (dp), INTENT(OUT)     :: xnew(:)
    REAL (dp), INTENT(IN OUT)  :: fnew(:)
    REAL (dp), INTENT(IN)      :: phi0
    REAL (dp), INTENT(IN)      :: dphi0
    REAL (dp), INTENT(IN OUT)  :: u
    REAL (dp), INTENT(OUT)     :: phiu
    REAL (dp), INTENT(IN)      :: tau
    REAL (dp), INTENT(IN OUT)  :: pmax
    !tpk REAL (dp), INTENT(IN OUT)  :: srelpr

    EXTERNAL ffunc

    
    REAL (dp)  :: phix, sqreps, x
    INTEGER    :: ctrl

    sqreps=SQRT(srelpr)
    phix=phi0
    x=u

10  IF(pmax*x < sqreps .OR. x <= alfmin) EXIT=-2
    IF(-2 == EXIT) GO TO 20
    x=x*0.5_dp
    ctrl=-1
    CALL fsumsq(xold,p,n,x,xnew,f,m,ffunc,ctrl,phix)
    ! write(10,*) 'In GAUC: PHIK,alfk= ',phix,x

    !     POSSIBLY THE USER CAN TERMINATE

    IF(ctrl < -10) k=ctrl
    IF(k < 0) RETURN
    k=k+1
    IF(phix > phi0 + tau*x*dphi0) GO TO 10

20  u=x
    phiu=phix

    RETURN
  END SUBROUTINE gauc


  !CHDER

  SUBROUTINE chder(dphize,phizer,xold,p,m,n,ffunc,k,EXIT,alfk,phik)

    ! Arguments XNEW & FNEW have been removed.

    REAL (dp), INTENT(IN)  :: dphize
    REAL (dp), INTENT(IN)  :: phizer
    REAL (dp), INTENT(IN)  :: xold(:)
    REAL (dp), INTENT(IN)  :: p(:)
    INTEGER, INTENT(IN)    :: m
    INTEGER, INTENT(IN)    :: n
    INTEGER, INTENT(OUT)   :: k
    INTEGER, INTENT(OUT)   :: EXIT
    REAL (dp), INTENT(IN)  :: alfk
    REAL (dp), INTENT(IN)  :: phik

    EXTERNAL ffunc

    
    !   MAKE A CONSISTENCY CHECK OF THE DERIVATIVE APPROXIMATION
    !   BASED ON THE JACOBIAN MATRIX
    !   THE FUNCTION UNDER CONSERN IS
    !        PHI(ALPHA) = F(XOLD+ALPHA*P)

    !   ON ENTRY

    !   DPHIZE THE DERIVATIVE OF PHI(ALPHA) AT ALPHA=0 COMPUTED
    !          BY P(TR)*G(XOLD)    WHERE G(X) IS THE GRADIENT OF F(X)
    !   PHIZER PHI(0)
    !   XOLD() CONTAINS THE STARTING POINT OF THE LINESEARCH
    !   P()    CONTAINS THE SEARCH DIRECTION
    !   M      NO. OF RESIDUALS IN THE VECTOR VALUED FUNCTION F(X)
    !   N      NO. OF UNKNOWNS
    !   FFUNC  SUBROUTINE NAME
    !   K      NO. OF FUNCTION EVALUATIONS DONE SO FAR IN THE LINESEARCH
    !   EXIT   =-2
    !   ALFK   THE ALPHA VALUE FOR WHICH THE DIFFERENCES ARE COMPUTED
    !   PHIK   PHI(ALFK)

    !   ON RETURN

    !   K      INCREASED BY 1 OR SET TO < -10 AS A USER STOP INDICATOR
    !   EXIT   SET = -1 IF INCONSISTENCY IS DETECTED

    !   WORKING AREAS

    !   XNEW() OF DIMENSION N
    !   FNEW() OF DIMENSION M

    INTEGER   :: ctrl
    REAL (dp) :: phimk,dphifo,dphiba,dphice,maxdif,xnew(n),fnew(m)

    !     COMPUTE PHI(-ALFK)

    ctrl=-1
    CALL fsumsq(xold,p,n,-alfk,xnew,fnew,m,ffunc,ctrl,phimk)
    IF(ctrl < -10) k=ctrl
    IF(k < 0) RETURN
    k=k+1

    !     COMPUTE APPROXIMATIONS OF THE DERIVATIVE BY USING FORWARD,
    !     BACKWARD AND CENTRAL DIFFERENCES

    IF(alfk <= zero) GO TO 20
    dphifo=(phik-phizer)/alfk
    dphiba=(phizer-phimk)/alfk
    dphice=(phik-phimk)/2.0_dp/alfk
    maxdif=ABS(dphifo-dphiba)
    maxdif=MAX(maxdif,ABS(dphifo-dphice))
    maxdif=MAX(maxdif,ABS(dphiba-dphice))
    IF((ABS(dphifo-dphize) > maxdif).AND. (ABS(dphice-dphize) > maxdif)) EXIT=-1
    RETURN

20  EXIT=-1

    RETURN
  END SUBROUTINE chder


  !QUAMIN

  SUBROUTINE quamin(x,fx,w,fw,v,fv,u)

    !     COMPUTE THE MINIMUM POINT U OF A QUADRATIC POLYNOMIAL PASSING
    !     THROUGH (V,F(V)), (W,F(W)) AND (X,F(X))

    REAL (dp), INTENT(IN)   :: x
    REAL (dp), INTENT(IN)   :: fx
    REAL (dp), INTENT(IN)   :: w
    REAL (dp), INTENT(IN)   :: fw
    REAL (dp), INTENT(IN)   :: v
    REAL (dp), INTENT(IN)   :: fv
    REAL (dp), INTENT(OUT)  :: u

    ! Internal variables
    REAL (dp)  :: d1, d2, s, q

    d1=fv-fx
    d2=fw-fx
    s=(w-x)*(w-x)*d1-(v-x)*(v-x)*d2
    q=2.d0*((v-x)*d2-(w-x)*d1)
    u=x-s/q

    RETURN
  END SUBROUTINE quamin


  !MINRN

  SUBROUTINE minrn(x,fx,w,fw,v,fv,alfmin,alfmax,pmax,srelpr,u,pu)

    !     PROVIDED THE POINTS X,V AND W ARE NOT TOO CLOSE, THE QUADRATIC
    !     PASSING THROUGH (X,FX), (V,FV) AND (W,FW) IS DETERMINED.
    !     U = THE MINIMUM POINT OF THIS QUADRATIC.
    !     PU = THE VALUE OF THE QUADRATIC AT U

    REAL (dp), INTENT(IN)      :: x
    REAL (dp), INTENT(IN)      :: fx
    REAL (dp), INTENT(IN OUT)  :: w
    REAL (dp), INTENT(IN)      :: fw
    REAL (dp), INTENT(IN OUT)  :: v
    REAL (dp), INTENT(IN)      :: fv
    REAL (dp), INTENT(IN OUT)  :: alfmin
    REAL (dp), INTENT(IN OUT)  :: alfmax
    REAL (dp), INTENT(IN)      :: pmax
    REAL (dp), INTENT(IN OUT)  :: srelpr
    REAL (dp), INTENT(OUT)     :: u
    REAL (dp), INTENT(OUT)     :: pu

    ! Internal variables
    REAL (dp)  :: eps, t1, t2, t3

    eps=SQRT(srelpr)/pmax
    u=x
    pu=fx
    IF(ABS(v-x) < eps .OR. ABS(w-x) < eps .OR. ABS(w-v) < eps) RETURN
    CALL quamin(x,fx,w,fw,v,fv,u)
    u=MIN(u,alfmax)
    u=MAX(u,alfmin)
    t1=(u-x)*(u-v)*fw/(w-x)/(w-v)
    t2=(u-w)*(u-v)*fx/(x-w)/(x-v)
    t3=(u-w)*(u-x)*fv/(v-x)/(v-w)
    pu=t1+t2+t3

    RETURN
  END SUBROUTINE minrn


  !UPDATE

  SUBROUTINE update(x,fx,w,fw,v,fv,u)

    REAL (dp), INTENT(OUT)     :: x
    REAL (dp), INTENT(OUT)     :: fx
    REAL (dp), INTENT(IN OUT)  :: w
    REAL (dp), INTENT(IN OUT)  :: fw
    REAL (dp), INTENT(IN OUT)  :: v
    REAL (dp), INTENT(IN)      :: fv
    REAL (dp), INTENT(IN)      :: u

    x=w
    fx=fw
    w=v
    fw=fv
    v=u

    RETURN
  END SUBROUTINE update


  !MINRM

  SUBROUTINE minrm(v0,v1,v2,m,alfmin,alfmax,xmin,x,px,y,py)

    REAL (dp), INTENT(IN OUT)  :: v0(:)
    REAL (dp), INTENT(IN OUT)  :: v1(:)
    REAL (dp), INTENT(IN OUT)  :: v2(:)
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(IN)      :: alfmin
    REAL (dp), INTENT(IN)      :: alfmax
    REAL (dp), INTENT(IN)      :: xmin
    REAL (dp), INTENT(IN OUT)  :: x
    REAL (dp), INTENT(OUT)     :: px
    REAL (dp), INTENT(OUT)     :: y
    REAL (dp), INTENT(OUT)     :: py

    !   A SUBROUTINE WHICH FINDS THE POINT X WHERE THE FUNCTION
    !   P(X)= 0.5*II V0+V1*X+V2*X**2 II**2  IS MINIMIZED.
    !   V0,V1 AND V2 BELONG TO R(M) AND X IS A SCALAR.
    !   THE VALUES OF V0,V1 AND V2 ARE BASED ON TW0 FUNCTION
    !   VALUES :  F(0) AND F(ALPHA).
    !   THE FUNCTION P(X) IS ALWAYS >=0 AND IT IS A POLYNOMIAL OF 4:TH DEGREE.
    !   THE MINIMUM VALUE OF P(X) IS ATTAINED AT A POINT X WHERE
    !   THE FIRST DERIVATIVE OF P(X)=DP(X) IS ZERO.
    !   DP(X) IS A POLYNOMIAL OF 3:RD DEGREE.
    !   IN CASE OF THREE ROOTS (X1<X2<X3), X2 CORRESPONDS TO A LOCAL
    !   MAXIMUM. CHOOSE THE ONE (OF X1 AND X3) THAT IS AT THE SAME
    !   SIDE OF THE MAXIMUM AS XMIN IS PROVIDED NO EXTRAPOLATION IS DONE.
    !   HOWEVER,THE MINIMUM POINT X MUST LIE IN THE INTERVALL (ALFMIN,ALFMAX) .
    !   PX IS SET TO P(X) AT THE MINIMUM POINT
    !   WHERE P(X) IS THE POLYNOMIAL ABOVE.
    !   Y IS SET TO THE OTHER MINIMIZER OF P(X)  (IF TWO ARE DETERMINED)
    !   OTHERWISE Y = X.
    !   PY = P(Y)

    !     INTERNAL VARIABLES

    INTEGER   :: k
    REAL (dp) :: v0norm,v1norm,v2norm,scv0v1,scv0v2,scv1v2,pprim,pbiss,  &
         d,error,h0,dm,hm,p,q,delta,a1div3,x0,x1,x2,x3
    REAL (dp), SAVE  :: eps = 1.E-04_dp

    !     COMPUTE NORMS AND SCALAR PRODUCTS

    CALL minrm1(v0,v1,v2,m,v0norm,v1norm,v2norm,scv0v1,scv0v2,scv1v2)
    pprim=pol3(scv0v1,scv0v2,v1norm,scv1v2,v2norm,xmin)
    pbiss=xmin**2*v2norm*v2norm*6.0_dp + xmin*scv1v2*6.0_dp + 2.0_dp*scv0v2 + &
         v1norm**2
    h0=ABS(pprim/pbiss)
    dm=ABS(6.d0*scv1v2 + xmin*v2norm*v2norm*12.0_dp) + h0*v2norm*v2norm*24.0d0

    !     DETERMINE IF DP(X)=0 SHOULD BE SOLVED BY USING NEWTONS METHOD

    hm=MAX(h0,one)
    IF(pbiss > 20.d0*hm*dm) GO TO 40

    !     COMPUTE QUANTITIES P,Q,DELTA AND A1DIV3 SO THAT THE SOLUTION OF
    !     X**3 + A1*X**2 + A2*X + A3 = 0    IS   X=T-A1/3
    !     WHERE T SOLVES     T**3 + P*T + Q = 0

    CALL minrm2(v1norm,v2norm,scv0v1,scv0v2,scv1v2,p,q,delta,a1div3)

    !     MATHEMATICALLY WE DESTINGWISH THREE DIFFERENT CASES:
    !     IF DELTA>0 THEN DF(X)=0 HAS ONE REAL ROOT
    !     IF DELTA=0 THEN DF(X)=0 HAS ONE SINGLE AND ONE DOUBLE REAL ROOT
    !     IF DELTA<0 THEN DF(X)=0 HAS THREE DIFFERENT REAL ROOTS

    !     IF DELTA=0 THE ONLY ROOT OF INTEREST IS THE SINGLE ONE,SO
    !     NUMERICALLY WE DISTINGWISH TWO CASES: DELTA>=0 AND DELTA<0

    IF(delta < zero) GO TO 30

    !     DELTA>=0 ONE INTERESTING ROOT. X

    CALL oner(q,delta,a1div3,x)
    y=x
    GO TO 100

    !     DELTA<0  TWO INTERESTING ROOTS.  Y AND Z, Y<Z

30  CALL twor(p,q,delta,a1div3,x1,x2,x3)

    !     CHOOSE X= X1 OR X2 OR X3

    CALL choose(x1,x2,x3,xmin,v0,v1,v2,m,x,y,py)
    GO TO 100

40  delta=one

    !     ITERATE USING NEWTONS METHOD

    k=0
    x0=xmin

50  pprim=pol3(scv0v1,scv0v2,v1norm,scv1v2,v2norm,x0)
    pbiss=x0**2*v2norm*v2norm*6.0_dp + x0*scv1v2*6.0_dp + 2.0_dp*scv0v2 + v1norm**2
    d=-pprim/pbiss
    x=x0 + d
    error=2.d0*dm*d*d/ABS(pbiss)
    x0=x
    k=k+1
    IF(error > eps .AND. k < 3) GO TO 50
    y=x

    !     MAKE THE MINIMUM POINT X LIE IN THE INTERVALL
    !     (ALFMIN,ALFMAX) AND EVALUATE F(X) AT THE MINIMUM POINT

100 x=MIN(x,alfmax)
    x=MAX(x,alfmin)
    px=pol4(v0,v1,v2,m,x)
    y=MIN(y,alfmax)
    y=MAX(y,alfmin)
    IF(delta >= zero) THEN
       y=x
       py=px
    END IF

    RETURN
  END SUBROUTINE minrm


  !MINRM1

  SUBROUTINE minrm1(v0,v1,v2,m,v0norm,v1norm,v2norm,scv0v1, scv0v2,scv1v2)

    !     COMPUTE THE EUCLIDEAN NORM OF THE M-DIMENSIONAL VECTORS V0, V1 AND V2

    REAL (dp), INTENT(IN OUT)  :: v0(:)
    REAL (dp), INTENT(IN OUT)  :: v1(:)
    REAL (dp), INTENT(IN OUT)  :: v2(:)
    INTEGER, INTENT(IN)        :: m
    REAL (dp), INTENT(OUT)     :: v0norm
    REAL (dp), INTENT(OUT)     :: v1norm
    REAL (dp), INTENT(OUT)     :: v2norm
    REAL (dp), INTENT(OUT)     :: scv0v1
    REAL (dp), INTENT(OUT)     :: scv0v2
    REAL (dp), INTENT(OUT)     :: scv1v2

    ! Internal variables
    REAL (dp)  :: sc1, sc2, sc3
    INTEGER    :: i

    v0norm=dnrm2(m,v0,1)
    v1norm=dnrm2(m,v1,1)
    v2norm=dnrm2(m,v2,1)

    !     SCALE THE VECTORS

    IF(v0norm /= zero) CALL scalv(v0,v0norm,m)
    IF(v1norm /= zero) CALL scalv(v1,v1norm,m)
    IF(v2norm /= zero) CALL scalv(v2,v2norm,m)

    !     COMPUTE THE SCALAR PRODUCTS V0(T)*V1, V0(T)*V2, V1(T)*V2

    sc1=zero
    sc2=zero
    sc3=zero
    DO  i=1,m
       sc1=sc1 + v0(i)*v1(i)
       sc2=sc2 + v0(i)*v2(i)
       sc3=sc3 + v1(i)*v2(i)
    END DO
    scv0v1=sc1*v0norm*v1norm
    scv0v2=sc2*v0norm*v2norm
    scv1v2=sc3*v1norm*v2norm

    !     RESCALE THE VECTORS

    IF(v0norm /= zero) CALL scalv(v0,one/v0norm,m)
    IF(v1norm /= zero) CALL scalv(v1,one/v1norm,m)
    IF(v2norm /= zero) CALL scalv(v2,one/v2norm,m)

    RETURN
  END SUBROUTINE minrm1


  !SCALV

  SUBROUTINE scalv(v,factor,n)

    REAL (dp), INTENT(IN OUT)  :: v(:)
    REAL (dp), INTENT(IN)      :: factor
    INTEGER, INTENT(IN)        :: n

    !     FORM THE NEW CONTENTS OF THE VECTOR V AS
    !            V(I)=1/FACTOR*V(I)  I=1,2,.....,N

    v(1:n)=v(1:n) / factor

    RETURN
  END SUBROUTINE scalv


  !MINRM2

  SUBROUTINE minrm2(v1nrm,v2nrm,scv0v1,scv0v2,scv1v2,p,q,delta,a1div3)

    REAL (dp), INTENT(IN)   :: v1nrm
    REAL (dp), INTENT(IN)   :: v2nrm
    REAL (dp), INTENT(IN)   :: scv0v1
    REAL (dp), INTENT(IN)   :: scv0v2
    REAL (dp), INTENT(IN)   :: scv1v2
    REAL (dp), INTENT(OUT)  :: p
    REAL (dp), INTENT(OUT)  :: q
    REAL (dp), INTENT(OUT)  :: delta
    REAL (dp), INTENT(OUT)  :: a1div3

    !     SUPPOSE WE HAVE A FOURTH DEGREE POLYNOMIAL II V0+V1*X+V2*X**2 II**2.
    !     THE DERIVATIVE IS A THIRD DEGREE POLYNOMIAL WHICH CAN BE
    !     WRITTEN    T**3+P*T+Q
    !     WHERE T=X+A1 AND P,Q DEFINED BELOW.

    ! Internal variables
    REAL (dp)  :: a1, a2, a3

    a1=1.5_dp*scv1v2/v2nrm/v2nrm
    a2=0.5_dp*((v1nrm/v2nrm)**2 + scv0v2/v2nrm/v2nrm*2.0_dp)
    a3=0.5_dp*scv0v1/v2nrm/v2nrm
    p=a2 - a1*a1/3._dp
    q=a3-a1*a2/3.d0 + 2.d0*a1*a1*a1/27.0_dp
    delta=(q/2.d0)**2 + (p/3.d0)**3
    a1div3=a1/3.d0

    RETURN
  END SUBROUTINE minrm2


  !ONER

  SUBROUTINE oner(q,delta,a1div3,x)

    !     COMPUTE THE ROOT OF A THIRD DEGREE POLYNOMIAL WHEN THERE
    !     IS ONLY ONE REAL ROOT

    REAL (dp), INTENT(IN)   :: q
    REAL (dp), INTENT(IN)   :: delta
    REAL (dp), INTENT(IN)   :: a1div3
    REAL (dp), INTENT(OUT)  :: x

    ! Internal variables
    REAL (dp)  :: arg1, arg2, a3rd, sqd, s1, s2, t

    sqd=SQRT(delta)
    arg1=-q/2.d0 + sqd
    s1=SIGN(1.d0,arg1)
    arg2=-q/2.d0 - sqd
    s2=SIGN(1.d0,arg2)
    a3rd=1.d0/3.d0
    t=s1*ABS(arg1)**a3rd + s2*ABS(arg2)**a3rd
    x=t - a1div3

    RETURN
  END SUBROUTINE oner


  !TWOR

  SUBROUTINE twor(p,q,delta,a1div3,x1,x2,x3)

    REAL (dp), INTENT(IN)   :: p
    REAL (dp), INTENT(IN)   :: q
    REAL (dp), INTENT(IN)   :: delta
    REAL (dp), INTENT(IN)   :: a1div3
    REAL (dp), INTENT(OUT)  :: x1
    REAL (dp), INTENT(OUT)  :: x2
    REAL (dp), INTENT(OUT)  :: x3

    !     COMPUTE THE THREE ROOTS OF A THIRD DEGREE POLYNOMIAL WHEN
    !     THERE ARE 3 REAL ROOTS

    ! Internal variables
    REAL (dp), PARAMETER  :: eps = 1.0E-8_dp
    REAL (dp)             :: fi, pi, sqd, t, tanfi

    sqd=SQRT(-delta)
    IF(ABS(q) <= 2.0_dp*eps*sqd) THEN
       fi=ATAN(one)*2.0_dp
    ELSE
       tanfi=ABS(2.0_dp*sqd/q)
       fi=ATAN(tanfi)
    END IF

    t=2.d0*SQRT(-p/3.d0)
    IF(q > zero) t=-t
    x1=t*COS(fi/3.d0) - a1div3
    pi=4.d0*ATAN(one)
    x2=t*COS((fi+2.d0*pi)/3.d0) - a1div3
    x3=t*COS((fi+4.d0*pi)/3.d0) - a1div3

    RETURN
  END SUBROUTINE twor


  !CHOOSE

  SUBROUTINE choose(x1,x2,x3,xmin,v0,v1,v2,m,root1,root2,proot2)

    REAL (dp), INTENT(IN)   :: x1
    REAL (dp), INTENT(IN)   :: x2
    REAL (dp), INTENT(IN)   :: x3
    REAL (dp), INTENT(IN)   :: xmin
    REAL (dp), INTENT(IN)   :: v0(:)
    REAL (dp), INTENT(IN)   :: v1(:)
    REAL (dp), INTENT(IN)   :: v2(:)
    INTEGER, INTENT(IN)     :: m
    REAL (dp), INTENT(OUT)  :: root1
    REAL (dp), INTENT(OUT)  :: root2
    REAL (dp), INTENT(OUT)  :: proot2

    !   X1, X2 AND X3 ARE THREE REAL ROOTS OF A 3:RD DEGREE POLYNOMIAL
    !   WHICH TENDS TO MINUS INFINITY WHEN X TENDS TO MINUS INFINITY.
    !   CHOOSE ONE OF THE OUTER ROOTS FOR ROOT1.
    !   ROOT2 = THE OTHER OUTER ONE. PROOT2 = THE VALUE OF THE CORRESPONDING
    !   4:TH DEGREE POLYNOMIAL AT ROOT2.

    ! Internal variables
    REAL (dp)  :: x, y, z

    x=MIN(x1,x2,x3)
    z=MAX(x1,x2,x3)
    IF(x1 <= x2 .AND. x1 <= x3) y=MIN(x2,x3)
    IF(x2 <= x1 .AND. x2 <= x3) y=MIN(x1,x3)
    IF(x3 <= x1 .AND. x3 <= x2) y=MIN(x1,x2)
    IF(xmin <= y) root1=x
    IF(xmin <= y) root2=z
    IF(xmin > y) root1=z
    IF(xmin > y) root2=x
    proot2=pol4(v0,v1,v2,m,root2)

    RETURN
  END SUBROUTINE choose


  !POL4

  FUNCTION pol4(v0,v1,v2,m,x) RESULT(fn_val)

    !     EVALUATE THE 4:TH DEGREE POLYNOMIAL || V2*X**2+V1*X+V0 ||**2 AT X

    REAL (dp), INTENT(IN)  :: v0(:)
    REAL (dp), INTENT(IN)  :: v1(:)
    REAL (dp), INTENT(IN)  :: v2(:)
    INTEGER, INTENT(IN)    :: m
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val

    ! Internal variables
    REAL (dp)  :: p, s
    INTEGER    :: i

    s=zero
    DO  i=1,m
       p=v0(i) + x*v1(i) + x**2*v2(i)
       s=s + p*p
    END DO
    fn_val=0.5_dp*s

    RETURN
  END FUNCTION pol4


  !POL3

  FUNCTION pol3(scv0v1,scv0v2,v1nrm,scv1v2,v2nrm,x) RESULT(fn_val)

    REAL (dp), INTENT(IN)  :: scv0v1
    REAL (dp), INTENT(IN)  :: scv0v2
    REAL (dp), INTENT(IN)  :: v1nrm
    REAL (dp), INTENT(IN)  :: scv1v2
    REAL (dp), INTENT(IN)  :: v2nrm
    REAL (dp), INTENT(IN)  :: x
    REAL (dp)              :: fn_val

    !   EVALUATE A 3:RD DEGREE POLYNOMIAL
    !   A3*X**3+A2*X**2+A1*X+A0   AT X
    !   WHERE  A0 = SCV0V1
    !         A1 = 2*SCV0V2+V1NRM**2
    !         A2 = 3*SCV1V2
    !         A3 = 2*V2NRM**2

    fn_val=scv0v1 + x*v1nrm*v1nrm + 2.0_dp*x*scv0v2 + x**2*3.0_dp*scv1v2 +  &
         x**3*v2nrm*v2nrm*2.0_dp

    RETURN
  END FUNCTION pol3


  !FSUMSQ

  SUBROUTINE fsumsq(xold,p,n,alfk,xnew,fnew,m,ffunc,ctrl,fn_val)

    !   EVALUATE FUNCTION VALUES AT THE POINT XOLD+ALFK*P .
    !   POSSIBLY THE USER CAN SIGNAL UNCOMPUTABILTY ON RETURN FROM
    !   THE USER WRITTEN ROUTINE FFUNC

    REAL (dp), INTENT(IN)    :: xold(:)
    REAL (dp), INTENT(IN)    :: p(:)
    INTEGER, INTENT(IN)      :: n
    REAL (dp), INTENT(IN)    :: alfk
    REAL (dp), INTENT(OUT)   :: xnew(:)
    REAL (dp), INTENT(OUT)   :: fnew(:)
    INTEGER, INTENT(IN)      :: m
    INTEGER, INTENT(IN OUT)  :: ctrl
    REAL (dp), INTENT(OUT)   :: fn_val

    
    REAL (dp) :: dummy(1,1)
    INTEGER   :: lctrl

    EXTERNAL ffunc


    fn_val=zero
    xnew(1:n)=xold(1:n) + alfk*p(1:n)
    lctrl=ctrl
    CALL ffunc(xnew,n,fnew,m,lctrl,dummy,1)
    !      write(10,*) 'XNEW and FNEW inside FSUMSQ'
    !      write(10,*) (xnew(i),i=1,n)
    !      write(10,*) (fnew(i),i=1,m)
    IF(ctrl == 1) GO TO 20

    !     CTRL=-1 ON ENTRY

    IF(-1 == lctrl) fn_val=0.5_dp*dnrm2(m,fnew,1)**2
    !      write(10,*) 'FSUMSQ when ctrl=-1', fsumsq
    IF(lctrl < -10) ctrl=lctrl
    RETURN

20  IF(lctrl == 1) fn_val=0.5_dp*dnrm2(m,fnew,1)**2
    !      write(10,*) 'FSUMSQ when ctrl= 1', fsumsq
    ctrl=lctrl

    RETURN
  END SUBROUTINE fsumsq


  !LINPAC

  SUBROUTINE sposl(a,n,b)

    ! Argument LDA has been removed.

    REAL (dp), INTENT(IN)      :: a(:,:)
    INTEGER, INTENT(IN)        :: n
    REAL (dp), INTENT(IN OUT)  :: b(:)

    !   SPOSL SOLVES THE REAL SYMMETRIC POSITIVE DEFINITE SYSTEM
    !   A * X = B
    !   USING THE FACTORS COMPUTED BY SPOCO OR SPOFA.

    !   ON ENTRY

    !      A       REAL(LDA, N)
    !              THE OUTPUT FROM SPOCO OR SPOFA.

    !      LDA     INTEGER
    !              THE LEADING DIMENSION OF THE ARRAY  A .

    !      N       INTEGER
    !              THE ORDER OF THE MATRIX  A .

    !      B       REAL(N)
    !              THE RIGHT HAND SIDE VECTOR.

    !   ON RETURN

    !      B       THE SOLUTION VECTOR  X .

    !   ERROR CONDITION

    !      A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
    !      A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
    !      SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
    !      ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
    !      CORRECTLY AND  INFO .EQ. 0 .

    !   TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
    !   WITH  P  COLUMNS
    !         CALL SPOCO(A,LDA,N,RCOND,Z,INFO)
    !         IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
    !         DO 10 J = 1, P
    !            CALL SPOSL(A,LDA,N,C(1,J))
    !      10 CONTINUE

    !   LINPACK.  THIS VERSION DATED 08/14/78 .
    !   CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.

    !   SUBROUTINES AND FUNCTIONS

    !   BLAS SAXPY,SDOT

    !   INTERNAL VARIABLES

    REAL (dp) :: t
    INTEGER   :: k,kb

    ! =======================================================
    ! Changes by Thomas P. Kurosu (tkurosu@cfa.harvard.edu).
    ! -------------------------------------------------------
    ! Date: 30 April, 2002
    ! Change: Avoid A(1:K-1) for K=1, since this will lead to
    !         core dumps if array checking is turned on.
    !         2 lines commented out, 3 lines added ("!tpk")
    ! =======================================================

    !     SOLVE TRANS(R)*Y = B

    DO  k = 1, n
       !tpk t = DOT_PRODUCT( a(1:k-1,k), b(1:k-1) )
       t = 0.0_dp  !tpk
       IF ( k > 1 ) t = DOT_PRODUCT( a(1:k-1,k), b(1:k-1) ) !tpk
       b(k) = (b(k) - t)/a(k,k)
    END DO

    !     SOLVE R*X = Y

    DO  kb = 1, n
       k = n + 1 - kb
       b(k) = b(k)/a(k,k)
       t = -b(k)
       !tpk b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
       IF ( k > 1 ) b(1:k-1) = b(1:k-1) + t * a(1:k-1,k) !tpk
    END DO

    RETURN
  END SUBROUTINE sposl


  !STRSL

  SUBROUTINE strsl(t,n,b,job,info)

    ! N.B. Argument LDT has been removed.

    REAL (dp), INTENT(IN)   :: t(:,:)
    INTEGER, INTENT(IN)     :: n
    REAL (dp), INTENT(OUT)  :: b(:)
    INTEGER, INTENT(IN)     :: job
    INTEGER, INTENT(OUT)    :: info

    !   STRSL SOLVES SYSTEMS OF THE FORM

    !                 T * X = B
    !   OR
    !                 TRANS(T) * X = B

    !   WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE TRANS(T)
    !   DENOTES THE TRANSPOSE OF THE MATRIX T.

    !   ON ENTRY

    !       T         REAL(LDT,N)
    !                 T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO
    !                 ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND
    !                 THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE
    !                 USED TO STORE OTHER INFORMATION.

    !       LDT       INTEGER
    !                 LDT IS THE LEADING DIMENSION OF THE ARRAY T.

    !       N         INTEGER
    !                 N IS THE ORDER OF THE SYSTEM.

    !       B         REAL(N).
    !                 B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM.

    !       JOB       INTEGER
    !                 JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED.
    !                 IF JOB IS

    !                      00   SOLVE T*X=B, T LOWER TRIANGULAR,
    !                      01   SOLVE T*X=B, T UPPER TRIANGULAR,
    !                      10   SOLVE TRANS(T)*X=B, T LOWER TRIANGULAR,
    !                      11   SOLVE TRANS(T)*X=B, T UPPER TRIANGULAR.

    !   ON RETURN

    !       B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0.
    !                 OTHERWISE B IS UNALTERED.

    !       INFO      INTEGER
    !                 INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR.
    !                 OTHERWISE INFO CONTAINS THE INDEX OF
    !                 THE FIRST ZERO DIAGONAL ELEMENT OF T.

    !   LINPACK. THIS VERSION DATED 08/14/78 .
    !   G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.

    !   SUBROUTINES AND FUNCTIONS

    !   BLAS SAXPY,SDOT
    !   FORTRAN MOD

    !   INTERNAL VARIABLES

    REAL (dp) :: temp
    INTEGER   :: case,j,jj

    !     BEGIN BLOCK PERMITTING ...EXITS TO 150

    !        CHECK FOR ZERO DIAGONAL ELEMENTS.

    DO  info = 1, n
       !     ......EXIT
       IF (t(info,info) == zero) GO TO 150
    END DO
    info = 0

    !        DETERMINE THE TASK AND GO TO IT.

    case = 1
    IF (MOD(job,10) /= 0) case = 2
    IF (MOD(job,100)/10 /= 0) case = case + 2
    SELECT CASE ( case )
    CASE (    1)
       !        SOLVE T*X=B FOR T LOWER TRIANGULAR

       b(1) = b(1)/t(1,1)
       DO  j = 2, n
          temp = -b(j-1)
          b(j:n) = b(j:n) + temp * t(j:n,j-1)
          b(j) = b(j)/t(j,j)
       END DO

    CASE (    2)
       !        SOLVE T*X=B FOR T UPPER TRIANGULAR.

       b(n) = b(n)/t(n,n)
       DO  jj = 2, n
          j = n - jj + 1
          temp = -b(j+1)
          b(1:j) = b(1:j) + temp * t(1:j,j+1)
          b(j) = b(j)/t(j,j)
       END DO

    CASE (    3)
       !        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.

       b(n) = b(n)/t(n,n)
       DO  jj = 2, n
          j = n - jj + 1
          b(j) = b(j) - DOT_PRODUCT( t(j+1:n,j), b(j+1:n) )
          b(j) = b(j)/t(j,j)
       END DO

    CASE (    4)
       !        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.

       b(1) = b(1)/t(1,1)
       DO  j = 2, n
          b(j) = b(j) - DOT_PRODUCT( t(1:j-1,j), b(1:j-1) )
          b(j) = b(j)/t(j,j)
       END DO

    END SELECT

150 RETURN
  END SUBROUTINE strsl


  !SCHDC

  SUBROUTINE schdc(a,p,work,jpvt,job,info)

    ! N.B. Argument LDA has been removed.

    REAL (dp), INTENT(IN OUT)  :: a(:,:)
    INTEGER, INTENT(IN)        :: p
    REAL (dp), INTENT(OUT)     :: work(:)
    INTEGER, INTENT(IN OUT)    :: jpvt(:)
    INTEGER, INTENT(IN)        :: job
    INTEGER, INTENT(OUT)       :: info

    !     SCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE
    !     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE
    !     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK
    !     OF A POSITIVE SEMIDEFINITE MATRIX.

    !     ON ENTRY

    !         A      REAL(LDA,P).
    !                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO
    !                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED.
    !                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED.

    !         LDA    INTEGER.
    !                LDA IS THE LEADING DIMENSION OF THE ARRAY A.

    !         P      INTEGER.
    !                P IS THE ORDER OF THE MATRIX.

    !         WORK   REAL.
    !                WORK IS A WORK ARRAY.

    !         JPVT   INTEGER(P).
    !                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION
    !                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED.
    !                EACH DIAGONAL ELEMENT A(K,K)
    !                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE
    !                VALUE OF JPVT(K).

    !                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL
    !                                      ELEMENT.

    !                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT.

    !                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT.

    !                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS
    !                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO
    !                THE BEGINNING OF THE ARRAY A AND FINAL
    !                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS
    !                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
    !                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE
    !                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT
    !                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT
    !                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF
    !                JOB .EQ. 0.

    !        JOB     INTEGER.
    !                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
    !                IF JOB .EQ. 0, NO PIVOTING IS DONE.
    !                IF JOB .NE. 0, PIVOTING IS DONE.

    !     ON RETURN

    !         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR
    !                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING.

    !         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT
    !                OF A THAT WAS MOVED INTO THE J-TH POSITION,
    !                PROVIDED PIVOTING WAS REQUESTED.

    !         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL
    !                ELEMENT OF THE CHOLESKY FACTOR.

    !     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN.
    !     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL
    !     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN
    !     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO
    !     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE
    !     INFO TO BE LESS THAN P.

    !     LINPACK. THIS VERSION DATED 03/19/79 .
    !     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND
    !     UNIVERSITY OF MARYLAND.


    !     BLAS SAXPY,SSWAP
    !     FORTRAN SQRT

    !     INTERNAL VARIABLES

    INTEGER   :: pu,pl,plp1,j,jp,jt,k,kb,km1,kp1,l,maxl
    REAL (dp) :: temp
    REAL (dp) :: maxdia
    LOGICAL   :: swapk,negk

    pl = 1
    pu = 0
    info = p
    IF (job == 0) GO TO 160

    !        PIVOTING HAS BEEN REQUESTED.
    !        REARRANGE THE THE ELEMENTS ACCORDING TO JPVT.

    DO  k = 1, p
       swapk = jpvt(k) > 0
       negk = jpvt(k) < 0
       jpvt(k) = k
       IF (negk) jpvt(k) = -jpvt(k)
       IF (.NOT.swapk) CYCLE
       IF (k == pl) GO TO 50
       CALL dswap(pl-1, a(:,k), 1, a(:,pl), 1)
       temp = a(k,k)
       a(k,k) = a(pl,pl)
       a(pl,pl) = temp
       plp1 = pl + 1
       DO  j = plp1, p
          IF (j < k) THEN
             temp = a(pl,j)
             a(pl,j) = a(j,k)
             a(j,k) = temp
          ELSE
             IF (j == k) CYCLE
             temp = a(k,j)
             a(k,j) = a(pl,j)
             a(pl,j) = temp
          END IF
       END DO

       jpvt(k) = jpvt(pl)
       jpvt(pl) = k

50     pl = pl + 1
    END DO

    pu = p
    DO  kb = pl, p
       k = p - kb + pl
       IF (jpvt(k) >= 0) CYCLE
       jpvt(k) = -jpvt(k)
       IF (pu == k) GO TO 120
       CALL dswap(k-1,a(:,k),1,a(:,pu),1)
       temp = a(k,k)
       a(k,k) = a(pu,pu)
       a(pu,pu) = temp
       kp1 = k + 1
       DO  j = kp1, p
          IF (j < pu) THEN
             temp = a(k,j)
             a(k,j) = a(j,pu)
             a(j,pu) = temp
          ELSE
             IF (j == pu) CYCLE
             temp = a(k,j)
             a(k,j) = a(pu,j)
             a(pu,j) = temp
          END IF
       END DO

       jt = jpvt(k)
       jpvt(k) = jpvt(pu)
       jpvt(pu) = jt

120    pu = pu - 1
    END DO

160 DO  k = 1, p

       !        REDUCTION LOOP.

       maxdia = a(k,k)
       kp1 = k + 1
       maxl = k

       !        DETERMINE THE PIVOT ELEMENT.

       IF (k < pl .OR. k >= pu) GO TO 190
       DO  l = kp1, pu
          IF (a(l,l) > maxdia) THEN
             maxdia = a(l,l)
             maxl = l
          END IF
       END DO

       !        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE.

190    IF (maxdia <= zero) THEN
          info = k - 1
          !     ......EXIT
          EXIT
       END IF
       IF (k == maxl) GO TO 210

       !           START THE PIVOTING AND UPDATE JPVT.

       km1 = k - 1
       CALL dswap(km1,a(:,k),1,a(:,maxl),1)
       a(maxl,maxl) = a(k,k)
       a(k,k) = maxdia
       jp = jpvt(maxl)
       jpvt(maxl) = jpvt(k)
       jpvt(k) = jp

       !        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.

210    work(k) = SQRT(a(k,k))
       a(k,k) = work(k)
       DO  j = kp1, p
          IF (k == maxl) GO TO 240
          IF (j < maxl) THEN
             temp = a(k,j)
             a(k,j) = a(j,maxl)
             a(j,maxl) = temp
          ELSE
             IF (j == maxl) GO TO 240
             temp = a(k,j)
             a(k,j) = a(maxl,j)
             a(maxl,j) = temp
          END IF

240       a(k,j) = a(k,j)/work(k)
          work(j) = a(k,j)
          temp = -a(k,j)
          a(kp1:j,j) = a(kp1:j,j) + temp * work(kp1:j)
       END DO
    END DO

    RETURN
  END SUBROUTINE schdc


  !SQRDC

  SUBROUTINE sqrdc(x,n,p,qraux,jpvt,job)

    ! N.B. Arguments LDX & WORK have been removed.

    REAL (dp), INTENT(IN OUT)  :: x(:,:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: p
    REAL (dp), INTENT(OUT)     :: qraux(:)
    INTEGER, INTENT(IN OUT)    :: jpvt(:)
    INTEGER, INTENT(IN)        :: job

    !   SQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR
    !   FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING
    !   BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE
    !   PERFORMED AT THE USERS OPTION.

    !   ON ENTRY

    !      X       REAL(LDX,P), WHERE LDX >= N.
    !              X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE COMPUTED.

    !      LDX     INTEGER.
    !              LDX IS THE LEADING DIMENSION OF THE ARRAY X.

    !      N       INTEGER.
    !              N IS THE NUMBER OF ROWS OF THE MATRIX X.

    !      P       INTEGER.
    !              P IS THE NUMBER OF COLUMNS OF THE MATRIX X.

    !      JPVT    INTEGER(P).
    !              JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION OF THE
    !              PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X IS PLACED IN ONE OF
    !              THREE CLASSES ACCORDING TO THE VALUE OF JPVT(K).

    !                 IF JPVT(K) > 0, THEN X(K) IS AN INITIAL COLUMN.

    !                 IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN.

    !                 IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN.

    !              BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS
    !              ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL
    !              COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS
    !              ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY
    !              FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE
    !              REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN
    !              IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST
    !              REDUCED NORM.  JPVT IS NOT REFERENCED IF
    !              JOB .EQ. 0.

    !      WORK    REAL(P).
    !              WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF
    !              JOB .EQ. 0.

    !      JOB     INTEGER.
    !              JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING.
    !              IF JOB .EQ. 0, NO PIVOTING IS DONE.
    !              IF JOB .NE. 0, PIVOTING IS DONE.

    !   ON RETURN

    !      X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER TRIANGULAR MATRIX R
    !              OF THE QR FACTORIZATION.   BELOW ITS DIAGONAL X CONTAINS
    !              INFORMATION FROM WHICH THE ORTHOGONAL PART OF THE
    !              DECOMPOSITION CAN BE RECOVERED.
    !              NOTE THAT IF PIVOTING HAS BEEN REQUESTED, THE DECOMPOSITION IS
    !              NOT THAT OF THE ORIGINAL MATRIX X BUT THAT OF X
    !              WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT.

    !      QRAUX   REAL(P).
    !              QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER
    !              THE ORTHOGONAL PART OF THE DECOMPOSITION.

    !      JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE
    !              ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO
    !              THE K-TH COLUMN, IF PIVOTING WAS REQUESTED.

    !   LINPACK. THIS VERSION DATED 08/14/78 .
    !   G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.

    !   SQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.

    !   BLAS SAXPY,SDOT,SSCAL,SSWAP,SNRM2
    !   FORTRAN ABS,MAX,MIN,SQRT

    !   INTERNAL VARIABLES

    INTEGER   :: j,jj,jp,l,lp1,lup,maxj,pl,pu
    REAL (dp) :: maxnrm,tt
    REAL (dp) :: nrmxl,t,work(p)
    LOGICAL   :: negj,swapj

    pl = 1
    pu = 0
    IF (job == 0) GO TO 60

    !     PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS ACCORDING TO JPVT.

    DO  j = 1, p
       swapj = jpvt(j) > 0
       negj = jpvt(j) < 0
       jpvt(j) = j
       IF (negj) jpvt(j) = -j
       IF (.NOT.swapj) CYCLE
       IF (j /= pl) CALL dswap(n,x(:,pl),1,x(:,j),1)
       jpvt(j) = jpvt(pl)
       jpvt(pl) = j
       pl = pl + 1
    END DO
    pu = p
    DO  jj = 1, p
       j = p - jj + 1
       IF (jpvt(j) >= 0) CYCLE
       jpvt(j) = -jpvt(j)
       IF (j /= pu) THEN
          CALL dswap(n,x(:,pu),1,x(:,j),1)
          jp = jpvt(pu)
          jpvt(pu) = jpvt(j)
          jpvt(j) = jp
       END IF
       pu = pu - 1
    END DO

    !     COMPUTE THE NORMS OF THE FREE COLUMNS.

60  DO  j = pl, pu
       qraux(j) = dnrm2(n,x(:,j),1)
       work(j) = qraux(j)
    END DO

    !     PERFORM THE HOUSEHOLDER REDUCTION OF X.

    lup = MIN(n,p)
    DO  l = 1, lup
       IF (l < pl .OR. l >= pu) GO TO 120

       !     LOCATE THE COLUMN OF LARGEST NORM AND BRING IT INTO THE PIVOT POSITION.

       maxnrm = zero
       maxj = l
       DO  j = l, pu
          IF (qraux(j) <= maxnrm) CYCLE
          maxnrm = qraux(j)
          maxj = j
       END DO
       IF (maxj == l) GO TO 120
       CALL dswap(n,x(:,l),1,x(:,maxj),1)
       qraux(maxj) = qraux(l)
       work(maxj) = work(l)
       jp = jpvt(maxj)
       jpvt(maxj) = jpvt(l)
       jpvt(l) = jp

120    qraux(l) = zero
       IF (l == n) CYCLE

       !           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.

       nrmxl = dnrm2(n-l+1,x(l:,l),1)
       IF (nrmxl == zero) CYCLE
       IF (x(l,l) /= zero) nrmxl = SIGN(nrmxl,x(l,l))
       x(l:n,l) = x(l:n,l) / nrmxl
       x(l,l) = one + x(l,l)

       !              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
       !              UPDATING THE NORMS.

       lp1 = l + 1
       DO  j = lp1, p
          t = -DOT_PRODUCT( x(l:n,l), x(l:n,j) ) / x(l,l)
          x(l:n,j) = x(l:n,j) + t * x(l:n,l)
          IF (j < pl .OR. j > pu) CYCLE
          IF (qraux(j) == zero) CYCLE
          tt = one - (ABS(x(l,j))/qraux(j))**2
          tt = MAX(tt,zero)
          t = tt
          tt = one + 0.05_dp*tt*(qraux(j)/work(j))**2
          IF (tt /= one) THEN
             qraux(j) = qraux(j)*SQRT(t)
          ELSE
             qraux(j) = dnrm2(n-l,x(l+1:,j),1)
             work(j) = qraux(j)
          END IF
       END DO

       !              SAVE THE TRANSFORMATION.

       qraux(l) = x(l,l)
       x(l,l) = -nrmxl
    END DO

    RETURN
  END SUBROUTINE sqrdc


  !SQRSL

  SUBROUTINE sqrsl(x,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)

    ! N.B. Argument LDX has been removed.

    REAL (dp), INTENT(IN OUT)  :: x(:,:)
    INTEGER, INTENT(IN)        :: n
    INTEGER, INTENT(IN)        :: k
    REAL (dp), INTENT(IN)      :: qraux(:)
    REAL (dp), INTENT(IN)      :: y(:)
    REAL (dp), INTENT(OUT)     :: qy(:)
    REAL (dp), INTENT(OUT)     :: qty(:)
    REAL (dp), INTENT(OUT)     :: b(:)
    REAL (dp), INTENT(OUT)     :: rsd(:)
    REAL (dp), INTENT(OUT)     :: xb(:)
    INTEGER, INTENT(IN)        :: job
    INTEGER, INTENT(OUT)       :: info

    !   SQRSL APPLIES THE OUTPUT OF SQRDC TO COMPUTE COORDINATE
    !   TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS.
    !   FOR K <= MIN(N,P), LET XK BE THE MATRIX

    !          XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K)))

    !   FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL
    !   N X P MATRIX X THAT WAS INPUT TO SQRDC (IF NO PIVOTING WAS
    !   DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR
    !   ORIGINAL ORDER).  SQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q
    !   AND AN UPPER TRIANGULAR MATRIX R SUCH THAT

    !            XK = Q * (R)
    !                     (0)

    !   THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS X AND QRAUX.

    !   ON ENTRY

    !      X      REAL(LDX,P).
    !             X CONTAINS THE OUTPUT OF SQRDC.

    !      LDX    INTEGER.
    !             LDX IS THE LEADING DIMENSION OF THE ARRAY X.

    !      N      INTEGER.
    !             N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST
    !             HAVE THE SAME VALUE AS N IN SQRDC.

    !      K      INTEGER.
    !             K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K
    !             MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE
    !             SAME AS IN THE CALLING SEQUENCE TO SQRDC.

    !      QRAUX  REAL(P).
    !             QRAUX CONTAINS THE AUXILIARY OUTPUT FROM SQRDC.

    !      Y      REAL(N)
    !             Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED
    !             BY SQRSL.

    !      JOB    INTEGER.
    !             JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS
    !             THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING MEANING.

    !                  IF A.NE.0, COMPUTE QY.
    !                  IF B,C,D, OR E .NE. 0, COMPUTE QTY.
    !                  IF C.NE.0, COMPUTE B.
    !                  IF D.NE.0, COMPUTE RSD.
    !                  IF E.NE.0, COMPUTE XB.

    !             NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB
    !             AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR
    !             WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING SEQUENCE.

    !   ON RETURN

    !      QY     REAL(N).
    !             QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN REQUESTED.

    !      QTY    REAL(N).
    !             QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS BEEN REQUESTED.
    !             HERE TRANS(Q) IS THE TRANSPOSE OF THE MATRIX Q.

    !      B      REAL(K)
    !             B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM

    !                  MINIMIZE NORM2(Y - XK*B),

    !             IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT
    !             IF PIVOTING WAS REQUESTED IN SQRDC, THE J-TH
    !             COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J)
    !             OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO SQRDC.)

    !      RSD    REAL(N).
    !             RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B,
    !             IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS
    !             ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE
    !             ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK.

    !      XB     REAL(N).
    !             XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B,
    !             IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO
    !             THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE OF X.

    !      INFO   INTEGER.
    !             INFO IS ZERO UNLESS THE COMPUTATION OF B HAS BEEN REQUESTED AND
    !             R IS EXACTLY SINGULAR.  IN THIS CASE, INFO IS THE INDEX OF THE
    !             THE FIRST ZERO DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED.

    !   THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED IF THEIR
    !   COMPUTATION IS NOT REQUESTED AND IN THIS CASE CAN BE REPLACED BY DUMMY
    !   VARIABLES IN THE CALLING PROGRAM.
    !   TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME ARRAY FOR
    !   DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A FREQUENTLY OCCURRING
    !   EXAMPLE IS WHEN ONE WISHES TO COMPUTE ANY OF B, RSD, OR XB AND DOES NOT
    !   NEED Y OR QTY.  IN THIS CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD,
    !   OR XB, WHILE PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE
    !   COMPUTED.  THUS THE CALLING SEQUENCE

    !        CALL SQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)

    !   WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD OVERWRITING Y.
    !   MORE GENERALLY, EACH ITEM IN THE FOLLOWING LIST CONTAINS GROUPS OF
    !   PERMISSIBLE IDENTIFICATIONS FOR A SINGLE CALLINNG SEQUENCE.

    !        1. (Y,QTY,B) (RSD) (XB) (QY)

    !        2. (Y,QTY,RSD) (B) (XB) (QY)

    !        3. (Y,QTY,XB) (B) (RSD) (QY)

    !        4. (Y,QY) (QTY,B) (RSD) (XB)

    !        5. (Y,QY) (QTY,RSD) (B) (XB)

    !        6. (Y,QY) (QTY,XB) (B) (RSD)

    !   IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO
    !   THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP.

    !   LINPACK. THIS VERSION DATED 08/14/78 .
    !   G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.

    !   SQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.

    !   BLAS SAXPY,SCOPY,SDOT
    !   FORTRAN ABS,MIN,MOD

    !   INTERNAL VARIABLES

    INTEGER   :: j,jj,ju,kp1
    REAL (dp) :: t,temp
    LOGICAL   :: cb,cqy,cqty,cr,cxb

    !     SET INFO FLAG.

    info = 0

    !     DETERMINE WHAT IS TO BE COMPUTED.

    cqy = job/10000 /= 0
    cqty = MOD(job,10000) /= 0
    cb = MOD(job,1000)/100 /= 0
    cr = MOD(job,100)/10 /= 0
    cxb = MOD(job,10) /= 0
    ju = MIN(k,n-1)

    !     SPECIAL ACTION WHEN N=1.

    IF (ju == 0) THEN
       IF (cqy) qy(1) = y(1)
       IF (cqty) qty(1) = y(1)
       IF (cxb) xb(1) = y(1)
       IF (.NOT.cb) GO TO 30
       IF (x(1,1) == zero) THEN
          info = 1
       ELSE
          b(1) = y(1)/x(1,1)
       END IF

30     IF (cr) rsd(1) = zero
       GO TO 240
    END IF

    !        SET UP TO COMPUTE QY OR QTY.

    IF (cqy) qy(1:n) = y(1:n)
    IF (cqty) qty(1:n) = y(1:n)
    IF (cqy) THEN

       !           COMPUTE QY.

       DO  jj = 1, ju
          j = ju - jj + 1
          IF (qraux(j) == zero) CYCLE
          temp = x(j,j)
          x(j,j) = qraux(j)
          t = -DOT_PRODUCT( x(j:n,j), qy(j:n) ) / x(j,j)
          qy(j:n) = qy(j:n) + t * x(j:n,j)
          x(j,j) = temp
       END DO
    END IF

    IF (cqty) THEN

       !           COMPUTE TRANS(Q)*Y.

       DO  j = 1, ju
          IF (qraux(j) == zero) CYCLE
          temp = x(j,j)
          x(j,j) = qraux(j)
          t = -DOT_PRODUCT( x(j:n,j), qty(j:n) ) / x(j,j)
          qty(j:n) = qty(j:n) + t * x(j:n,j)
          x(j,j) = temp
       END DO
    END IF

    !        SET UP TO COMPUTE B, RSD, OR XB.

    IF (cb) b(1:k) = qty(1:k)
    kp1 = k + 1
    IF (cxb) xb(1:k) = qty(1:k)
    IF (cr .AND. k < n) rsd(kp1:n) = qty(kp1:n)
    IF (.NOT.cxb .OR. kp1 > n) GO TO 120
    xb(kp1:n) = zero

120 IF (.NOT.cr) GO TO 140
    rsd(1:k) = zero

140 IF (.NOT.cb) GO TO 190

    !           COMPUTE B.

    DO  jj = 1, k
       j = k - jj + 1
       IF (x(j,j) == zero) THEN
          info = j
          !           ......EXIT
          EXIT
       END IF
       b(j) = b(j)/x(j,j)
       IF (j == 1) CYCLE
       t = -b(j)
       b(1:j-1) = b(1:j-1) + t * x(1:j-1,j)
    END DO

190 IF (.NOT.cr .AND. .NOT.cxb) GO TO 240

    !           COMPUTE RSD OR XB AS REQUIRED.

    DO  jj = 1, ju
       j = ju - jj + 1
       IF (qraux(j) == zero) CYCLE
       temp = x(j,j)
       x(j,j) = qraux(j)
       IF (cr) THEN
          t = -DOT_PRODUCT( x(j:n,j), rsd(j:n) ) / x(j,j)
          rsd(j:n) = rsd(j:n) + t * x(j:n,j)
       END IF
       IF (cxb) THEN
          t = -DOT_PRODUCT( x(j:n,j), xb(j:n) ) / x(j,j)
          xb(j:n) = xb(j:n) + t * x(j:n,j)
       END IF
       x(j,j) = temp
    END DO

240 RETURN
  END SUBROUTINE sqrsl


  FUNCTION dnrm2 ( n, x, incx) RESULT(fn_val)

    !  Euclidean norm of the n-vector stored in x() with storage increment incx .
    !  if n <= 0 return with result = 0.
    !  if n >= 1 then incx must be >= 1

    !  c.l.lawson, 1978 jan 08
    !  modified to correct failure to update ix, 1/25/92.
    !  modified 3/93 to return if incx <= 0.
    !  This version by Alan.Miller @ vic.cmis.csiro.au
    !  Latest revision - 22 January 1999

    !  four phase method using two built-in constants that are
    !  hopefully applicable to all machines.
    !      cutlo = maximum of  SQRT(u/eps)  over all known machines.
    !      cuthi = minimum of  SQRT(v)      over all known machines.
    !  where
    !      eps = smallest no. such that eps + 1. > 1.
    !      u   = smallest positive no.   (underflow limit)
    !      v   = largest  no.            (overflow  limit)

    !  brief outline of algorithm..

    !  phase 1    scans zero components.
    !  move to phase 2 when a component is nonzero and <= cutlo
    !  move to phase 3 when a component is > cutlo
    !  move to phase 4 when a component is >= cuthi/m
    !  where m = n for x() real and m = 2*n for complex.

    INTEGER, INTENT(IN)   :: n, incx
    REAL (dp), INTENT(IN) :: x(:)
    REAL (dp)             :: fn_val

    ! Local variables
    INTEGER    :: i, ix, j, next
    REAL (dp)  :: cuthi, cutlo, hitest, sum, xmax

    IF(n <= 0 .OR. incx <= 0) THEN
       fn_val = zero
       RETURN
    END IF

    ! Set machine-dependent constants

    cutlo = SQRT( TINY(one) / EPSILON(one) )
    cuthi = SQRT( HUGE(one) )

    next = 1
    sum = zero
    i = 1
    ix = 1
 
    !                                                 begin main loop
20  SELECT CASE (next)
    CASE (1)
       IF( ABS(x(i)) > cutlo) GO TO 85
       next = 2
       xmax = zero
       GO TO 20

    CASE (2)
       !                   phase 1.  sum is zero

       IF( x(i) == zero) GO TO 200
       IF( ABS(x(i)) > cutlo) GO TO 85

       !                                prepare for phase 2.   x(i) is very small.
       next = 3
       GO TO 105

    CASE (3)
       !                   phase 2.  sum is small.
       !                             scale to avoid destructive underflow.

       IF( ABS(x(i)) > cutlo ) THEN
          !                  prepare for phase 3.

          sum = (sum * xmax) * xmax
          GO TO 85
       END IF

    CASE (4)
       GO TO 110
    END SELECT

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !                     common code for phases 2 and 4.
    !                     in phase 4 sum is large.  scale to avoid overflow.

110 IF( ABS(x(i)) <= xmax ) GO TO 115
    sum = one + sum * (xmax / x(i))**2
    xmax = ABS(x(i))
    GO TO 200

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !                   phase 3.  sum is mid-range.  no scaling.

    !     for real or d.p. set hitest = cuthi/n
    !     for complex      set hitest = cuthi/(2*n)

85  hitest = cuthi / n

    DO j = ix, n
       IF(ABS(x(i)) >= hitest) GO TO 100
       sum = sum + x(i)**2
       i = i + incx

    END DO
    fn_val = SQRT( sum )
    
    RETURN

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !                                prepare for phase 4.
    !                                ABS(x(i)) is very large
100 ix = j
    next = 4
    sum = (sum / x(i)) / x(i)
    !                                Set xmax; large if next = 4, small if next = 3
105 xmax = ABS(x(i))

115 sum = sum + (x(i)/xmax)**2
    
200 ix = ix + 1
    i = i + incx
    IF( ix <= n ) GO TO 20

    !              end of main loop.

    !              compute square root and adjust for scaling.

    fn_val = xmax * SQRT(sum)

    RETURN
  END FUNCTION dnrm2


  ! ============= dswap.f ==============
  SUBROUTINE dswap (n, x, incx, y, incy)

    !     interchanges two vectors.

    INTEGER, INTENT(IN)       :: n, incx, incy
    REAL (dp), INTENT(IN OUT) :: x(:), y(:)

    ! Local variables
    REAL (dp) :: temp(n)

    IF(n <= 0) RETURN
    IF(incx == 1 .AND. incy == 1) THEN
       temp = x(:n)
       x(:n) = y(:n)
       y(:n) = temp
       RETURN
    END IF

    temp = x(:n*incx:incx)
    x(:n*incx:incx) = y(:n*incy:incy)
    y(:n*incy:incy) = temp

    RETURN
  END SUBROUTINE dswap

END MODULE bounded_nonlin_LS
