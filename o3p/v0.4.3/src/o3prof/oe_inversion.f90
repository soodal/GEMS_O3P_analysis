! The original code, developed by Roeland van Oss, was converted to FORTRAN 90
! and modified for by X. Liu

! Note from X. Liu
! do_sa_diagonal: assume no covariance in the a priori covariance matrix
! do_oe_output:   write intermediate results
! file_unit:      file unitto write intermediate results
! delta_chi_min:  convergence criteria
! last_iter:      is this the last iteration (will calculate averaging kernels, 
!                 sx, sxn, contri, otherwise unnecessary)
! num_iter:       current iteration
! y:              differences between observed and simulated radiances
! dy:             measurement error
! rk:             weighting function matrix
! xap:            a priori state vector
! xold:           previus state vector
! sa:             a priori covariance matrix
! sidx, eidx:     start and end positions of key variables (e.g., ozone) in the state vector
! x:              updated retrieval vector
! sx:             covariance matrix (both smoothing and random-noise)
! sxn:            random-noise covariance matrix
! conv:           Does it converge?
! rKernel:        averaging kernels
! contri:         contribution function
! ozdfs :         DFS for ozone variables
! ozinfo:         information content for ozone variables
! chi_new:        updated cost function
! ynew:           differences between observed and simulated radiances after updating state vector                    
                                                                        
! Notes from Roeland's code                                            
! =========================                                            
                                                                        
! Optimal estimation retrieval:                                         
! the optimal estimation for state vector x obeying:                    
!       y=Kx                                                            
! is calculated.                                                        
! Sy is the measurement error covariance matrix; only the sqrt diagonal 
! is passed as vector dy(ny).                                           
!                                                                       
! xa is the apriori for the state vector, Sa the aprior covariance matri
! and x the retrieved state vector.                                     
!                                                                       
! Method described in: C.D. Rogers "Inverse Methods for Atmospheric Soun
! Theory and Practice" (Draft version 06-02-1998) Chapter 2.6.          
!                                                                       
! Firstly, all quantities are transformed such that the error covariance
! matrices are unit matrices.                                           
! The new matrix Ktil = Sy^(-1/2)*K*Sa^(1/2) is then decomposed         
! with SVD: Ktil = U*W*V^T, using NumRec DSVDCMP.                       
!                                                                       
! The retrieved state vector x is calculated by:                        
! x_prime = (W^t*W+I)^(-1)(W*y_prime + xa_prime)                        
! with:                                                                 
!       x(_a)_prime=VT*Sa^(-1/2)*x(_a),                                 
!       x(_a)=Sa^(1/2)*V*x(_a)_prime,                                   
!       y_prime=UT*Sy^(-1/2)*y.                                         
! The error covariance of the state vector x is calculated by:          
!       Sx=Sa^(1/2)*V*(WT*W+I)^(-1)*VT*Sa^(1/2)                         
!                                                                       
!               Author: Roeland van Oss, KNMI, 2000    
                                                                 
SUBROUTINE oe_inversion (do_sa_diagonal, do_oe_output, file_unit, delta_chi_min, last_iter, &
     num_iter, ny, nx, y, dy, rk, xap, xold, sa, xname, sidx, eidx, x, sx, sxn, conv, &
     rkernel, contr, ozdfs, ozinfoh, chi_new, ynew)  

 USE OMSAO_precision_module
 IMPLICIT NONE
     
 ! ====================================
 ! Input / Output variaibles
 ! ====================================
 LOGICAL, INTENT(IN)         :: do_sa_diagonal, do_oe_output, last_iter
 INTEGER, INTENT(IN)         :: file_unit, num_iter, nx, ny, sidx, eidx
 REAL (KIND=dp), INTENT (IN) :: delta_chi_min
 REAL (KIND=dp), DIMENSION(ny, nx), INTENT (IN)  :: rK 
 REAL (KIND=dp), DIMENSION(ny), INTENT (IN)      :: y, dy
 REAL (KIND=dp), DIMENSION(nx, nx), INTENT(IN)   :: Sa 
 REAL (KIND=dp), DIMENSION(nx), INTENT (IN)      :: xap, xold
 CHARACTER(len=6), DIMENSION(nx), INTENT(IN)     :: xname
 
 LOGICAL, INTENT(OUT)                            :: conv                                   
 REAL (KIND=dp), INTENT(OUT)                     :: ozdfs, ozinfoh, chi_new
 REAL (KIND=dp), DIMENSION(nx), INTENT(OUT)      :: x
 REAL (KIND=dp), DIMENSION(ny), INTENT (OUT)     :: ynew
 REAL (KIND=dp), DIMENSION(nx, nx), INTENT (OUT) :: Sx, Sxn, rkernel
 REAL (KIND=dp), DIMENSION(nx, ny), INTENT (OUT) :: contr
                                                                        
 ! =======================                                                              
 ! Local variables                                                       
 ! =======================         
 INTEGER                          :: i, j, l, k 
 REAL(KIND=dp)                    :: dfn, chi_old, wsa_min, dfs, h,  wsa_max, delta_chi 
 REAL(KIND=dp), DIMENSION(nx)     :: y_prime, xa_prime, x_prime, wsa, wsa_inv, w, xastd, xa, &
      tmpw2, tmpw2p1, tmpratio
 REAL (KIND=dp), DIMENSION(ny)    :: ytmp
 REAL(KIND=dp), DIMENSION(ny, nx) :: rK_tilde, u, tmp2
 REAL(KIND=dp), DIMENSION(nx, nx) :: usa, vsa, sasqp, sasqn, tmp, tmp1, v

 conv = .FALSE. 
 xa = xap - xold
 DO i = 1, nx
    xastd(i) = SQRT(Sa (i, i))
 ENDDO

 IF (do_sa_diagonal) THEN      ! Diagonal Apriori    
    ! Construct rK_tilde = Sy^(-1/2)*K*Sa^(1/2) (u=Ktil for dsvdcmp) 
    DO j = 1, nx           
       rK_tilde (:, j) = xastd(j) * rK (:, j) / dy 
    ENDDO
    u = rK_tilde
    
    ! SVD of Ktil:                                                          
    CALL dsvdcmp (u, ny, nx, ny, nx, w, v) 
    tmpw2 = w ** 2;    tmpw2p1 = tmpw2 + 1.0d0; tmpratio = tmpw2 / tmpw2p1
    
    ! Construct y_prime=UT*Sy^(-1/2)*y:                                     
    DO i = 1, nx 
       y_prime (i) = SUM ( u(:, i) * y / dy )
    ENDDO
                                                                        
    ! Construct xa_prime=VT*Sa^(-1/2)*xa:                                   
    DO i = 1, nx 
       xa_prime (i) = SUM ( v(:, i) * xa  / xastd ) 
    ENDDO
                                                                        
    ! ACTUAL RETRIEVAL: Calculate x_prime = (WT*W+I)^(-1)(W*y_prime + xa_prime) 
    x_prime = (w * y_prime + xa_prime ) / tmpw2p1
    
    ! Construct x=Sa^(1/2)*V*x_prime:                                       
    DO i = 1, nx 
       x (i) = SUM ( v (i, :) * x_prime * xastd(i) )
    ENDDO
    
    ! Construct Sx=Sa^(1/2)*V* (WT*W+I)^(-1) *VT*Sa^(1/2)                   
    DO i = 1, nx 
       DO j = 1, nx 
          Sx (i, j) = SUM( v (i, :) * v (j, :) / tmpw2p1 )  * xastd(i) * xastd(j)
       ENDDO
    ENDDO
    
    IF (last_iter .OR. do_oe_output) THEN
       ! Construct rkernel=Sa^(1/2)*V* (W^t*W+I)^(-1)*WT*W *V^t*Sa^(-1/2)      
       DO i = 1, nx 
          DO j = 1, nx 
             rkernel (i, j) = SUM (v (i, :) * v (j, :) * tmpratio ) * xastd(i) * xastd(j)
          ENDDO
       ENDDO
    ENDIF     
                                                                                                   
 ELSE ! Sa is not diagonal      
    
    ! SVD of Sa to calculate Sa^(1/2) and Sa^(-1/2)                         
    ! NOTE: Sa is square and symmetric: u=v                                 
    usa = sa
    CALL dsvdcmp (usa, nx, nx, nx, nx, wsa, vsa) 
   
    ! Make reciprokal of very small eigenvalues zero: ! replaced with follow
    wsa_min = 1.d-16
    wsa_inv = 0.0d0
    DO i = 1, nx 
       IF (wsa (i) > wsa_min) wsa_inv (i) = 1.0d0 / wsa (i) 
    ENDDO
                                                                        
    ! Sa^(1/2) = U*W^(1/2)*UT and Sa^(-1/2) = U*W^(-1/2)*UT                 
    DO i = 1, nx 
       DO j = 1, nx 
          sasqp (i, j) = SUM ( usa (i, :) * usa (j, :) * SQRT (wsa) )
          sasqn (i, j) = SUM ( usa (i, :) * usa (j, :) * SQRT (wsa_inv) )
       ENDDO
    ENDDO
    
    ! Construct rK_tilde = Sy^(-1/2)*K*Sa^(1/2) (u=Ktil for dsvdcmp)        
    ! NOTE: Sy is assumed to be diagonal!!!:                                
    DO i = 1, ny 
       DO j = 1, nx 
          rK_tilde (i, j) = SUM (rK (i, :) * sasqp (:, j))
          u (i, j) = rK_tilde (i, j) / dy (i) 
       ENDDO
    ENDDO
                                                                          
    ! SVD of Ktil:                                                          
    CALL dsvdcmp (u, ny, nx, ny, nx, w, v) 
    tmpw2 = w ** 2;    tmpw2p1 = tmpw2 + 1.0d0; tmpratio = tmpw2 / tmpw2p1
    
    ! Construct y_prime=UT*Sy^(-1/2)*y:                                     
    DO i = 1, nx 
       y_prime (i) = SUM (u (:, i) * y / dy) 
    ENDDO
    
    ! Construct xa_prime=VT*Sa^(-1/2)*xa: 
    DO i = 1, nx 
       DO j = 1, nx 
          tmp (i, j) = SUM (v (:, i) * sasqn (:, j))
       ENDDO
       xa_prime (i) = SUM (tmp (i, :) * xa)
    ENDDO
    
    ! ACTUAL RETRIEVAL: Calculate x_prime = (WT*W+I)^(-1)(W*y_prime + xa_prime)
    x_prime = (w * y_prime + xa_prime ) / tmpw2p1   
    
    ! Construct x=Sa^(1/2)*V*x_prime:                                                                                   
    DO i = 1, nx 
       DO j = 1, nx 
          tmp (i, j) = SUM(sasqp(i, :) * v (:, j)) 
       ENDDO
       x (i) =  SUM(tmp (i, :) * x_prime)
    ENDDO
    
    ! Construct Sx=Sa^(1/2)*V* (WT*W+I)^(-1) *VT*Sa^(1/2)                   
    ! First: tmp = Sa^(1/2)*V* (WT*W+I)^(-1)                                
    DO i = 1, nx 
       DO j = 1, nx 
          tmp (i, j) =  SUM (sasqp (i, :) * v (:, j) / tmpw2p1(j)) 
       ENDDO
    ENDDO

    ! then: tmp1 = tmp*VT:                                                  
    DO i = 1, nx 
       DO j = 1, nx 
          tmp1 (i, j) = SUM( tmp (i, :) * v (j, :))
       ENDDO
    ENDDO

    ! then: sx = tmp1*Sa^(1/2):                                             
    DO i = 1, nx 
       DO j = 1, nx 
          sx (i, j) = SUM (tmp1 (i, :) * sasqp (:, j) )
       ENDDO
    ENDDO
    
    ! IF (last_iter .OR. do_oe_output) THEN                        
    
    ! Construct rkernel=Sa^(1/2)*V* (W^t*W+I)^(-1)*WT*W *V^t*Sa^(-1/2)      
    ! First: tmp = Sa^(1/2)*V* (WT*W+I)^(-1)*WT*W                           
    DO i = 1, nx 
       DO j = 1, nx 
          tmp (i, j) = SUM (sasqp (i, :) * v (:, j) * tmpratio(j) )
       ENDDO
    ENDDO

    ! then: tmp1 = tmp*VT:                                                  
    DO i = 1, nx 
       DO j = 1, nx 
          tmp1 (i, j) = SUM(tmp (i, :) * v (j, :))  
       ENDDO
    ENDDO

    ! then: rkernel = tmp1*Sa^(-1/2):                                       
    DO i = 1, nx 
       DO j = 1, nx 
          rkernel (i, j) = SUM (tmp1 (i, :) * sasqn (:, j) )
       ENDDO
    ENDDO
                                                                        
    ! endif          
 ENDIF  ! Sa diagonal or not diagonal  
                                                                        
! Improvement in Chi:                                                   
 chi_old = SQRT(SUM( (y / dy) ** 2) / ny)
 DO i = 1, ny 
    ytmp(i) = SUM (rK (i, :) * x )
 ENDDO
 chi_new   =  SQRT (SUM(((y - ytmp) / dy ) **2) / ny)
 ynew      = y - ytmp
 delta_chi = ABS ( (chi_new - chi_old) / chi_old) 
                                                                        
 ! Check convergence                                                     
 IF (delta_chi < delta_chi_min) conv = .TRUE. 
 
 IF (last_iter .OR. do_oe_output) THEN   
    ! Construct contribution function: Sx K^T Sy^(-1):                  
    DO i = 1, nx 
       DO j = 1, ny 
          contr(i, j) = SUM(sx(i, :) * rK(j, :)) / (dy(j) ** 2) 
       ENDDO
    ENDDO  
 ENDIF
                                                                        
! Construct retrieval noise covariance matrix: Sx K^T Sy^(-1) K Sx      
 DO i = 1, ny 
    DO j = 1, nx 
       tmp2(i, j) = SUM(rK (i, :) * sx (:, j)) / (dy (i) ** 2)    
    ENDDO
 ENDDO

 DO i = 1, nx 
    DO j = 1, nx 
       tmp(i, j) =  SUM(rK (:, i) * tmp2 (:, j))
    ENDDO
 ENDDO

 DO i = 1, nx 
    DO j = 1, nx 
       Sxn(i, j) = SUM(Sx(:, i) * tmp (:, j))
    ENDDO
 ENDDO
 
 ! Degrees of Freedom Noise and Signal, Information content (dfn,dfs,h): 
 IF (last_iter .OR. do_oe_output) THEN  
    
    dfn = SUM(1.0d0 / tmpw2p1)
    h = SUM(0.5d0 * LOG (tmpw2p1) )
    dfs = nx - dfn 
    
    ozdfs = 0.0 
    DO i = sidx, eidx 
       ozdfs = ozdfs + rkernel (i, i) 
    ENDDO

    ! need to check for this later        
    ozinfoh = ozdfs / dfs * h    
 ENDIF
                                                                        
 !  Level 2 output debug                                                 
 !  --------------------                                                 
 IF (do_oe_output) THEN
    ! chi-square                                                            
    WRITE (file_unit, '(A, I5)')    'Iteration = ', num_iter 
    WRITE (file_unit, '(A, D14.6)') 'Old Chi   = ', chi_old 
    WRITE (file_unit, '(A, D14.6)') 'New Chi   = ', chi_new 
    WRITE (file_unit, '(A, D14.6, A1,D14.6)') 'Delchi / limit value = ', &
         delta_chi, '/', delta_chi_min                                  
    
    ! Degrees of Freedom Noise and Signal, Information content (dfn,dfs,h): 
    WRITE (file_unit, '(2(A, D14.6))') 'DFN =       ', dfn, ' OZDFN  = ',  (eidx - sidx + 1)  - ozdfs                                   
    WRITE (file_unit, '(2(A, D14.6))') 'DFS =       ', dfs, ' OZDFS  = ', ozdfs                                                         
    WRITE (file_unit, '(2(A, D14.6))') 'Information content = ', h,  ' OZINFO = ', ozinfoh                                          
    
    !! Eigenvalues:                                                          
    !WRITE(file_unit,'(A)') ' Eigenvalues:'                       
    !WRITE(file_unit,'(A)') ' columns, points'                    
    !WRITE(file_unit,'(I3,I5)') 1, nx                               
    !DO i = 1, nx                                                 
    !   WRITE(file_unit,'(I3,5D15.7)') i, w(i)                      
    !ENDDO
    
    ! A priori and its error and retrieved state and error                                       
    WRITE (file_unit, '(A6,6A14)') '  Var  ', ' retrieved ', ' noise error', 'smooth error', &
         'Previous ', ' apriori  ', ' apr. error '  

    DO i = 1, nx 
       WRITE (file_unit, '(A6, 6D14.6)') xname (i), x (i) + xold (i), SQRT (sxn (i, i) ), &
            SQRT (Sx (i, i) ), xold (i), xap (i), xastd(i)                                             
    ENDDO
    
 ENDIF

 RETURN 
END SUBROUTINE oe_inversion


! --------------------------------------------------------------------
SUBROUTINE inverse (a, n, b) ! Invert matrix by Gauss method
! --------------------------------------------------------------------
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL(8), DIMENSION(n, n), INTENT(IN)  :: a
REAL(8), DIMENSION(n, n), INTENT(OUT) :: b

! - - - Local Variables - - -
REAL(8) :: c, d, temp(n)
INTEGER :: i, j, k, m, imax(1), ipvt(n)
! - - - - - - - - - - - - - -

b = a
ipvt = (/ (i, i = 1, n) /)

DO k = 1,n
imax = MAXLOC(ABS(b(k:n,k)))
m = k-1+imax(1)

IF (m /= k) THEN
ipvt( (/m,k/) ) = ipvt( (/k,m/) )
b((/m,k/),:) = b((/k,m/),:)
END IF
d = 1/b(k,k)

temp = b(:,k)
DO j = 1, n
c = b(k,j)*d
b(:,j) = b(:,j)-temp*c
b(k,j) = c
END DO
b(:,k) = temp*(-d)
b(k,k) = d
END DO

b(:,ipvt) = b

END SUBROUTINE inverse
