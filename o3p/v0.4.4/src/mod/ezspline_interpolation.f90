SUBROUTINE ezspline_interpolation ( n_in, x_in, y_in, n_out, x_out, y_out, errstat )

  USE EZspline_obj
  USE EZspline  
  USE OMSAO_precision_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER,                           INTENT (IN) :: n_in, n_out
  REAL (KIND=dp), DIMENSION (n_in),  INTENT (IN) :: x_in, y_in
  REAL (KIND=dp), DIMENSION (n_out), INTENT (IN) :: x_out

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER,                           INTENT (OUT) :: errstat
  REAL (KIND=dp), DIMENSION (n_out), INTENT (OUT) :: y_out

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER, PARAMETER :: ezs_r8 = SELECTED_REAL_KIND(12,100) ! real*8 kind

  REAL (KIND=ezs_r8), DIMENSION (n_in)  :: x1, f ! independent variable and function
  REAL (KIND=ezs_r8), DIMENSION (n_out) :: z1, fz
  INTEGER,            DIMENSION (2)     :: BCS1(2)
  TYPE (EZspline1_r8)                   :: spline_o ! 1-d EZspline object

  INTEGER :: i, errstat1

  !! ---------------
  !! Allocate memory
  !! ---------------
  !IF ( ALLOCATED(x1) ) DEALLOCATE(x1) ;  IF ( ALLOCATED(f) ) DEALLOCATE(f)
  !ALLOCATE ( x1(n_in), f(n_in), STAT=errstat )
  BCS1 = (/  0,  0 /) ! not a knot

  ! --------------------------
  ! Initialize/allocate memory
  ! --------------------------

  CALL EZspline_init (spline_o, n_in, BCS1, errstat)
  IF ( errstat /= pge_errstat_ok ) THEN
       CALL Ezspline_free (spline_o, errstat1);  RETURN
  ENDIF

  ! ---------------------------
  ! Set explicitely spline_o%x1
  ! ---------------------------
  spline_o%x1 = REAL(x_in(1:n_in), KIND=ezs_r8)

  ! ---------------------
  ! Assign function value
  ! ---------------------
  f(1:n_in) = REAL ( y_in(1:n_in), KIND=ezs_r8 )

  ! ----------------------------
  ! Setting up interpolation ...
  ! ----------------------------
  CALL EZspline_setup(spline_o, f, errstat)
  IF ( errstat /= pge_errstat_ok ) THEN
       CALL Ezspline_free (spline_o, errstat1); RETURN
  ENDIF

  ! -------------------
  ! Array interpolation
  ! -------------------
  !IF ( ALLOCATED(z1) ) DEALLOCATE(z1) ;  IF ( ALLOCATED(fz) ) DEALLOCATE(fz)
  !ALLOCATE ( z1(n_out), fz(n_out), STAT=errstat )

  z1 = REAL ( x_out, KIND=ezs_r8 )

  CALL EZspline_interp (spline_o, n_out, z1, fz, errstat)
  IF ( errstat /= pge_errstat_ok ) THEN
       CALL Ezspline_free (spline_o, errstat1); RETURN
  ENDIF

  ! -----------------------------------------------
  ! Assign interpolated results to output variables
  ! -----------------------------------------------
  y_out(1:n_out) = fz(1:n_out)

  ! ---------------------------
  ! Clean up and free up memory
  ! ---------------------------
  CALL Ezspline_free (spline_o, errstat)
  IF ( errstat /= pge_errstat_ok ) RETURN

  RETURN
END SUBROUTINE ezspline_interpolation


! This subroutine combine spline and splint fnuction
SUBROUTINE BSPLINE(xa, ya, n, x, y, np, errstat)

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  
  INTEGER, INTENT (IN)                      :: n, np
  INTEGER, INTENT (OUT)                     :: errstat
  REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
  REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y

  
  REAL (KIND=dp), DIMENSION(n)              :: y2a, xacpy
  REAL (KIND=dp), DIMENSION(np)             :: xcpy
  REAL (KIND=dp), DIMENSION(n-1)            :: diff
  REAL (KIND=dp)                            :: xmin, xmax, xamin, xamax
  
  errstat = 0
  IF (n < 3) THEN
     errstat = - 1; RETURN
  ENDIF
  
  diff = xa(2:n) - xa(1:n-1)
  IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
     errstat =  -2; RETURN
  ENDIF

  xmin = MINVAL(x); xmax = MAXVAL(x)
  xamin = MINVAL(xa); xamax = MAXVAL(xa)
  IF (xmin < xamin .OR. xmax > xamax) THEN
     errstat =  -3; RETURN
  ENDIF
  
  IF (xa(1) < xa(n)) THEN
     CALL spline(xa, ya, n, y2a)
     CALL splint(xa, ya, y2a, n, x, y, np)
  ELSE
     xacpy = -xa; xcpy = -x
     CALL spline(xacpy, ya, n, y2a)
     CALL splint(xacpy, ya, y2a, n, xcpy, y, np)
  ENDIF

  RETURN
END SUBROUTINE BSPLINE


! This subroutine combine spline and splint function
SUBROUTINE BSPLINE1(xa, ya, n, x, y, dydx, np, errstat)

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  
  INTEGER, INTENT (IN)                      :: n, np
  INTEGER, INTENT (OUT)                     :: errstat
  REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
  REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y, dydx
  
  REAL (KIND=dp), DIMENSION(n)              :: y2a, xacpy
  REAL (KIND=dp), DIMENSION(np)             :: xcpy
  REAL (KIND=dp), DIMENSION(n-1)            :: diff
  REAL (KIND=dp)                            :: xmin, xmax, xamin, xamax
  
  errstat = 0
  IF (n < 3) THEN
     errstat = - 1; RETURN
  ENDIF
  
  diff = xa(2:n) - xa(1:n-1)
  IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
     errstat =  -2; RETURN
  ENDIF

  xmin = MINVAL(x); xmax = MAXVAL(x)
  xamin = MINVAL(xa); xamax = MAXVAL(xa)
  IF (xmin < xamin .OR. xmax > xamax) THEN
     errstat =  -3; RETURN
  ENDIF

  IF (xa(1) < xa(n)) THEN
     CALL spline(xa, ya, n, y2a)
     CALL splint1(xa, ya, y2a, n, x, y, dydx, np)
  ELSE
     xacpy = -xa; xcpy = -x
     CALL spline(xacpy, ya, n, y2a)
     CALL splint1(xacpy, ya, y2a, n, xcpy, y, dydx, np)
  ENDIF
  
  RETURN
END SUBROUTINE BSPLINE1

! This subroutine provide a switch for calculating shifting wf
SUBROUTINE BSPLINE2(xa, ya, n, cal_shiwf, x, y, dydx, np, errstat)

  IMPLICIT NONE  
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  
  INTEGER, INTENT (IN)                      :: n, np
  INTEGER, INTENT (OUT)                     :: errstat
  LOGICAL, INTENT (IN)                      :: cal_shiwf
  REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
  REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y, dydx
  
  REAL (KIND=dp), DIMENSION(n)              :: y2a, xacpy
  REAL (KIND=dp), DIMENSION(np)             :: xcpy
  REAL (KIND=dp), DIMENSION(n-1)            :: diff
  REAL (KIND=dp)                            :: xmin, xmax, xamin, xamax
  
  errstat = 0
  IF (n < 3) THEN
     errstat = - 1; RETURN
  ENDIF
  
  diff = xa(2:n) - xa(1:n-1)
  IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
     errstat =  -2; RETURN
  ENDIF

  xmin = MINVAL(x); xmax = MAXVAL(x)
  xamin = MINVAL(xa); xamax = MAXVAL(xa)
  IF (xmin < xamin .OR. xmax > xamax) THEN
     errstat =  -3; RETURN
  ENDIF

  IF (xa(1) < xa(n)) THEN
     CALL spline(xa, ya, n, y2a)
     IF (cal_shiwf) THEN
        CALL splint1(xa, ya, y2a, n, x, y, dydx, np)
     ELSE
        CALL splint(xa, ya, y2a, n, x, y, np)
     END IF
  ELSE
     xacpy = -xa; xcpy = -x
     CALL spline(xacpy, ya, n, y2a)
     IF (cal_shiwf) THEN
        CALL splint1(xacpy, ya, y2a, n, xcpy, y, dydx, np)
     ELSE
        CALL splint(xacpy, ya, y2a, n, xcpy, y, np)
     END IF
  ENDIF
  
  RETURN
END SUBROUTINE BSPLINE2


! modified to always use "natural" boundary conditions
SUBROUTINE SPLINE (x, y, n, y2)
  
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  INTEGER, INTENT(IN) :: n
  REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: x, y  
  REAL (KIND=dp), DIMENSION(n), INTENT(OUT) :: y2
  
  REAL (KIND=dp), DIMENSION(n)     		 :: u
  INTEGER       :: i, k
  REAL(KIND=dp) :: sig, p, qn, un
  
  y2 (1) = 0.0
  u (1) = 0.0
  
  DO i = 2, n - 1
     sig = (x (i) - x (i - 1)) / (x (i + 1) -x (i - 1))
     p = sig * y2 (i - 1) + 2.D0
     y2 (i) = (sig - 1.) / p
     u (i) = (6._dp * ((y (i + 1) - y (i)) / (x (i + 1) - x (i)) -  & 
          (y (i) - y (i - 1)) / (x (i) - x (i - 1))) / (x (i + 1) - &
          x (i - 1)) - sig * u (i - 1)) / p
  ENDDO
  
  qn = 0.0
  un = 0.0
  y2 (n) = (un - qn * u (n - 1)) / (qn * y2 (n - 1) + 1.0)
  DO k = n - 1, 1, -1
     y2 (k) = y2 (k) * y2 (k + 1) + u (k)
  ENDDO
  
  RETURN
END SUBROUTINE SPLINE


SUBROUTINE SPLINT1 (xa, ya, y2a, n, x, y, dy1, m)
  
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  INTEGER, INTENT(IN) :: n, m 
  REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: xa, ya, y2a
  REAL (KIND=dp), DIMENSION(m), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(m), INTENT(OUT):: y, dy1
  
  INTEGER        :: ii, klo, khi, k 
  REAL (KIND=dp) :: h, a, b
    
  DO ii = 1, m 
     klo = 1
     khi = n
     DO WHILE (khi - klo > 1)
        k = (khi + klo) / 2
        IF (xa (k) > x(ii)) THEN
           khi = k
        ELSE
           klo = k
        ENDIF
     ENDDO
     
     h = xa (khi) - xa (klo)
     IF (h == 0.0) STOP 'Bad xa input in: splint!!!'
     a = (xa (khi) - x(ii)) / h
     b = (x(ii) - xa (klo)) / h
     
     y(ii) = a * ya (klo) + b * ya (khi) + ((a**3 - a) * y2a (klo) + &
          (b**3 - b) * y2a (khi)) * (h**2) / 6.0
     
     dy1(ii)=(-1.0D0/h) * ya(klo) + (1.0D0/h) * ya(khi) + &
          (-(3.D0*a**2 - 1.0D0)  * y2a(klo) / h + (3.0D0*b**2 - 1.0D0) * &
          y2a(khi)/h) * (h**2) / 6.D0     
  ENDDO
  
  RETURN
END SUBROUTINE SPLINT1

! This code could be optimized if x is in increasing/descreasing order
SUBROUTINE SPLINT (xa, ya, y2a, n, x, y, m)
  
  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  INTEGER, INTENT(IN) :: n, m
  REAL (KIND=dp), DIMENSION(n), INTENT(IN) :: xa, ya, y2a
  REAL (KIND=dp), DIMENSION(m), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(m), INTENT(OUT):: y
  
  INTEGER        :: ii, klo, khi, k 
  REAL (KIND=dp) :: h, a, b

  !klo = 1; khi = n
 
  DO ii = 1, m 
     klo = 1; khi = n
    
     !IF ( khi - klo == 1) THEN
     !   IF (x(ii) > xa(khi) )   THEN
     !       khi = n
     !   ENDIF
     !ENDIF

     DO WHILE (khi - klo > 1)
        k = (khi + klo) / 2
        IF (xa (k) > x(ii)) THEN
           khi = k
        ELSE
           klo = k
        ENDIF
     ENDDO
     
     h = xa (khi) - xa (klo)
     IF (h == 0.0) STOP 'Bad xa input in: splint!!!'
     a = (xa (khi) - x(ii)) / h
     b = (x(ii) - xa (klo)) / h
     
     y(ii) = a * ya (klo) + b * ya (khi) + ((a**3 - a) * y2a (klo) + &
          (b**3 - b) * y2a (khi)) * (h**2) / 6.0
  ENDDO
  
  RETURN
END SUBROUTINE SPLINT

! This is not optimized (could save the time to find the index by using previous indices)
SUBROUTINE INTERPOL(xa, ya, n, x, y, np, errstat)

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  
  INTEGER, INTENT (IN)                      :: n, np
  INTEGER, INTENT (OUT)                     :: errstat
  REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
  REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y

  REAL (KIND=dp), DIMENSION(n-1)            :: diff
  REAL (KIND=dp), DIMENSION(np)             :: frac
  INTEGER,        DIMENSION(np)             :: fidx
  INTEGER                                   :: i
 
  errstat = 0
  IF (n <= 1) THEN
     errstat = - 1; RETURN
  ENDIF
  
  diff = xa(2:n) - xa(1:n-1)
  IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
     errstat =  -2; RETURN
  ENDIF
  
  ! find the node
  IF (xa(1) < xa(n)) THEN
     DO i = 1, np
        IF (x(i) <= xa(1)) THEN
           fidx(i) = 1  ! extrapolation
        ELSE 
           fidx(i) = MIN(n-1, MINVAL(MAXLOC(xa, MASK=(xa <= x(i)))))
        ENDIF
     ENDDO
  ELSE
     DO i = 1, np
        IF (x(i) >= xa(1)) THEN
           fidx(i) = 1     ! extrapolation
        ELSE
           fidx(i) = MIN(n-1, MINVAL(MINLOC(xa, MASK=(xa >= x(i)))))
        ENDIF
    ENDDO
  ENDIF
  
  frac = (x - xa(fidx)) / (xa(fidx+1) - xa(fidx))

  y = (1.0D0 - frac) * ya(fidx) + frac * ya(fidx + 1)
  
  RETURN

END SUBROUTINE INTERPOL

! Also return the first derivative of dy/dx
SUBROUTINE INTERPOL2(xa, ya, n, calc_shiwf, x, y, dydx, np, errstat)

  IMPLICIT NONE
  INTEGER, PARAMETER  :: dp = KIND(1.0D0)
  
  INTEGER, INTENT (IN)                      :: n, np
  LOGICAL, INTENT (IN)                      :: calc_shiwf
  INTEGER, INTENT (OUT)                     :: errstat
  REAL (KIND=dp), DIMENSION(n),  INTENT(IN) :: xa, ya
  REAL (KIND=dp), DIMENSION(np), INTENT(IN) :: x
  REAL (KIND=dp), DIMENSION(np), INTENT(OUT):: y, dydx

  REAL (KIND=dp), DIMENSION(n-1)            :: diff
  REAL (KIND=dp), DIMENSION(np)             :: frac
  INTEGER,        DIMENSION(np)             :: fidx
  INTEGER                                   :: i
 
  errstat = 0
  IF (n <= 1) THEN
     errstat = - 1; RETURN
  ENDIF
  
  diff = xa(2:n) - xa(1:n-1)
  IF (.NOT. (ALL(diff > 0) .OR. ALL(diff < 0))) THEN
     errstat =  -2; RETURN
  ENDIF
  
  ! find the node
  IF (xa(1) < xa(n)) THEN
     DO i = 1, np
        IF (x(i) <= xa(1)) THEN
           fidx(i) = 1  ! extrapolation
        ELSE 
           fidx(i) = MIN(n-1, MINVAL(MAXLOC(xa, MASK=(xa <= x(i)))))
        ENDIF
     ENDDO
  ELSE
     DO i = 1, np
        IF (x(i) >= xa(1)) THEN
           fidx(i) = 1     ! extrapolation
        ELSE
           fidx(i) = MIN(n-1, MINVAL(MINLOC(xa, MASK=(xa >= x(i)))))
        ENDIF
    ENDDO
  ENDIF
  
  frac = (x - xa(fidx)) / (xa(fidx+1) - xa(fidx))

  y = (1.0D0 - frac) * ya(fidx) + frac * ya(fidx + 1)

  IF (calc_shiwf) THEN
     dydx = (ya(fidx+1) - ya(fidx)) / (xa(fidx+1) - xa(fidx))
  ENDIF
  
  RETURN

END SUBROUTINE INTERPOL2



SUBROUTINE blend_101 ( r, x0, x1, x )
!
!*******************************************************************************
!
!! BLEND_101 extends scalar endpoint data to a line.
!
!
!  Diagram:
!
!    0-----r-----1
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    14 December 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, the coordinate where an interpolated value is desired.  
!
!    Input, real X0, X1, the data values at the ends of the line.
!
!    Output, real X, the interpolated data value at (R).
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)
  
  REAL (KIND=dp) :: r
  REAL (KIND=dp) :: x
  REAL (KIND=dp) :: x0
  REAL (KIND=dp) :: x1
  
  x = ( 1.0E+00 - r ) * x0 + r * x1

  RETURN
END SUBROUTINE blend_101


SUBROUTINE blend_102 ( r, s, x00, x01, x10, x11, x )
!
!*******************************************************************************
!
!! BLEND_102 extends scalar point data into a square.
!
!
!  Diagram:
!
!    01------------11
!     |      .      |
!     |      .      |
!     |.....rs......|
!     |      .      |
!     |      .      |
!    00------------10
!
!  Formula:
!
!    Written in terms of R and S, the map has the form:
!
!      X(R,S) =
!               1     * ( + x00 )
!             + r     * ( - x00 + x10 )
!             + s     * ( - x00       + x01 )
!             + r * s * ( + x00 - x10 - x01 + x11 )
!
!    Written in terms of the coefficients, the map has the form:
!
!      X(R,S) = x00 * ( 1 - r - s + r * s )
!             + x01 * (         s - r * s )
!             + x10 * (     r     - r * s )
!             + x11 * (             r * s )
!
!             = x00 * ( 1 - r ) * ( 1 - s )
!             + x01 * ( 1 - r ) *       s
!             + x10 *       r   * ( 1 - s )
!             + x11 *       r           s
!
!    The nonlinear term ( r * s ) has an important role:
!
!      If ( x01 + x10 - x00 - x11 ) is zero, then the input data lies in
!      a plane, and the mapping is affine.  All the interpolated data 
!      will lie on the plane defined by the four corner values.  In 
!      particular, on any line through the square, data values at 
!      intermediate points will lie between the values at the endpoints.  
!
!      If ( x01 + x10 - x00 - x11 ) is not zero, then the input data does
!      not lie in a plane, and the interpolation map is nonlinear.  On
!      any line through the square, data values at intermediate points
!      may lie above or below the data values at the endpoints.  The
!      size of the coefficient of r * s will determine how severe this
!      effect is.
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, the coordinates where an interpolated value is 
!    desired.  
!
!    Input, real X00, X01, X10, X11, the data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S).
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  REAL (KIND=dp) :: r
  REAL (KIND=dp) :: s
  REAL (KIND=dp) :: x
  REAL (KIND=dp) :: x00
  REAL (KIND=dp) :: x01
  REAL (KIND=dp) :: x10
  REAL (KIND=dp) :: x11

  x =             + x00 &
       + r *     ( - x00 + x10 ) & 
       + s *     ( - x00       + x01 ) &
       + r * s * ( + x00 - x10 - x01 + x11 )
  
  RETURN
END SUBROUTINE blend_102


SUBROUTINE blend_103(r, s, t, x000, x001, x010, x011, x100, x101, x110, x111, x )
!
!*******************************************************************************
!
!! BLEND_103 extends scalar point data into a cube.
!
!
!  Diagram:
!
!    011--------------111 
!      |               |
!      |               | 
!      |               |
!      |               |
!      |               |
!    001--------------101
!
!
!      *---------------*
!      |               |
!      |               |
!      |      rst      |
!      |               |
!      |               |
!      *---------------*
!
!
!    010--------------110
!      |               |
!      |               |
!      |               |
!      |               | 
!      |               |
!    000--------------100 
!
!
!  Formula:
!
!    Written as a polynomial in R, S and T, the interpolation map has the 
!    form:
!
!      X(R,S,T) =
!        1         * ( + x000 )
!      + r         * ( - x000 + x100 )
!      +     s     * ( - x000        + x010 )
!      +         t * ( - x000               + x001 )
!      + r * s     * ( + x000 - x100 - x010                       + x110 )
!      + r     * t * ( + x000 - x100        - x001        + x101 )
!      +     s * t * ( + x000        - x010 - x001 + x011 )
!      + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
!
!  Reference:
!
!    William Gordon,
!    Blending-Function Methods of Bivariate and Multivariate Interpolation
!      and Approximation,
!    SIAM Journal on Numerical Analysis,
!    Volume 8, Number 1, March 1971, pages 158-177.
!
!    William Gordon and Charles Hall,
!    Transfinite Element Methods: Blending-Function Interpolation over
!      Arbitrary Curved Element Domains,
!    Numerische Mathematik,
!    Volume 21, Number 1, 1973, pages 109-129.
!
!    William Gordon and Charles Hall,
!    Construction of Curvilinear Coordinate Systems and Application to
!      Mesh Generation,
!    International Journal of Numerical Methods in Engineering,
!    Volume 7, 1973, pages 461-477.
!
!    Joe Thompson, Bharat Soni, Nigel Weatherill,
!    Handbook of Grid Generation,
!    CRC Press, 1999.
!
!  Modified:
!
!    11 October 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, S, T, the coordinates where an interpolated value
!    is desired.
!
!    Input, real X000, X001, X010, X011, X100, X101, X110, X111, the
!    data values at the corners.
!
!    Output, real X, the interpolated data value at (R,S,T).
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)

  REAL (KIND=dp) :: r
  REAL (KIND=dp) :: s
  REAL (KIND=dp) :: t
  REAL (KIND=dp) :: x
  REAL (KIND=dp) :: x000
  REAL (KIND=dp) :: x001
  REAL (KIND=dp) :: x010
  REAL (KIND=dp) :: x011
  REAL (KIND=dp) :: x100
  REAL (KIND=dp) :: x101
  REAL (KIND=dp) :: x110
  REAL (KIND=dp) :: x111

!  Interpolate the interior point.

  x =  1.0E+00     * ( + x000 ) &
       + r         * ( - x000 + x100 ) &
       +     s     * ( - x000        + x010 ) &
       +         t * ( - x000               + x001 ) &
       + r * s     * ( + x000 - x100 - x010                      + x110 ) &
       + r     * t * ( + x000 - x100        - x001        + x101 ) &
       +     s * t * ( + x000        - x010 - x001 + x011 ) &
       + r * s * t * ( - x000 + x100 + x010 + x001 - x011 - x101 - x110 + x111 )
  
  RETURN
END SUBROUTINE blend_103


SUBROUTINE reverse ( inarr, num )
  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = KIND(1.0D0)
  
  INTEGER, INTENT(IN) :: num
  INTEGER             :: i
  REAL (KIND=dp), DIMENSION(1: num), INTENT(INOUT) :: inarr
  REAL (KIND=dp), DIMENSION(1: num)                :: temp

  DO i = 1, num
     temp(i) = inarr(num - i + 1)
  ENDDO
  inarr = temp

  RETURN
END SUBROUTINE reverse


