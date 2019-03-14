  ! ***************************************************************
  ! Author:  xiong liu (xliu)
  ! Date  :  July 23, 2003
  ! Purpose: Generate Tikhonov regularization matrix
  !   
  ! Jan 14, 2004: Correct the Sobolev form, get L^T not L
  ! Jan 14, 2004: Use a priori covariance for smoothing (option 6)
  ! *** Note ***: For option 6, need to add diagonal elements for 
  !               other non-ozone state variables
  ! ***************************************************************

SUBROUTINE get_reg_matrix(sa, ptr_nump, ptr_b, pge_error_status)

  USE OMSAO_precision_module
  USE ozprof_data_module,        ONLY: nlay_fit, ozfit_start_index, &
       ozfit_end_index, ptr_order, ptr_w0, ptr_w1, ptr_w2, &
       tf_fidx, tf_lidx, nlay
  USE OMSAO_parameters_module,   ONLY: maxchlen                 
  USE OMSAO_variables_module,    ONLY: n_fitvar_rad
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! -----------------------
  ! Iinput/Output variable
  ! -----------------------
  REAL (KIND=dp), DIMENSION(nlay_fit, nlay_fit), INTENT (IN)          :: sa
  REAL (KIND=dp), DIMENSION(n_fitvar_rad, n_fitvar_rad), INTENT (OUT) :: ptr_b  
  INTEGER, INTENT(OUT) :: ptr_nump, pge_error_status
  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER                                       :: i, j
  REAL (KIND=dp), DIMENSION(nlay_fit)           :: parr
  REAL (KIND=dp), DIMENSION(nlay_fit, nlay_fit) :: bb0, bb1, bb2, bb_all, sainv

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=14), PARAMETER :: modulename = 'get_reg_matrix'

  
  pge_error_status = pge_errstat_ok  ! Initialize exit status

  IF (ptr_order == 0) THEN           ! zero order
     ptr_nump = nlay_fit             ! only apply to ozone profile variables
  ELSE IF (ptr_order == 1) THEN      ! first order
     ptr_nump = nlay_fit - 1      
  ELSE IF (ptr_order == 2) THEN
     ptr_nump = nlay_fit - 2    
  ELSE IF (ptr_order == 3) THEN      ! x - xavg
     ptr_nump = nlay_fit   
  ELSE IF (ptr_order == 4) THEN      ! Sobolev form, need further check, does sound correct
     ptr_nump = nlay_fit   
  ELSE IF (ptr_order == 5) THEN      ! include temperature variables, zero order
     ptr_nump = nlay_fit + tf_lidx - tf_fidx + 1
  ELSE IF (ptr_order == 6) THEN      ! use a priori covariance matrix L^T * L = Sa^(-1) 
     ptr_nump = nlay_fit   
  ELSE
     WRITE(*, *) modulename, ' : This order of PTR is not implemented!!!'
     pge_error_status = pge_errstat_error; RETURN
  END IF

  IF (ptr_nump <= 3) THEN
     WRITE(*, *) modulename, ' : Please increase # regularized parameters!!!'
     pge_error_status = pge_errstat_error; RETURN
  END IF

  ! initialize regularization matrix (ptr_nump, n_fitvar_rad)
  ptr_b(1:ptr_nump, 1:n_fitvar_rad) = 0.0_dp

  IF (ptr_order == 0) THEN
     i = 0
     DO j = ozfit_start_index, ozfit_end_index
        i = i + 1
        ptr_b(i, j) = 1.0_dp
     END DO

  ELSE IF (ptr_order == 5) THEN
     i = 0
     DO j = ozfit_start_index, tf_lidx
        i = i + 1
        ptr_b(i, j) = 1.0_dp
     END DO

  ELSE IF (ptr_order == 1) THEN
     i = 0
     DO j = ozfit_start_index, ozfit_end_index - 1
        i = i + 1
        ptr_b(i, j) = -1.0_dp
        ptr_b(i, j + 1) = 1.0_dp
     END DO

  ELSE IF (ptr_order == 2) THEN
     i = 0
     DO j = ozfit_start_index, ozfit_end_index - 2
        i = i + 1
        ptr_b(i, j) = -1.0_dp
        ptr_b(i, j + 1) = 2.0_dp
        ptr_b(i, j + 2) = -1.0_dp
     END DO

  ELSE IF (ptr_order == 3) THEN
     ptr_b(1:ptr_nump, ozfit_start_index:ozfit_end_index) = -1.0_dp / &
          REAL (nlay_fit, KIND=dp)
     i = 0
     DO j = ozfit_start_index, ozfit_end_index
        i = i + 1
        ptr_b(i, j) = (REAL(nlay_fit, KIND=dp) - 1.0_dp) / &
             REAL(nlay_fit, KIND=dp)
     END DO

  ELSE IF (ptr_order == 4) THEN
     bb0(1:ptr_nump, 1:ptr_nump) = 0.0_dp
     bb1(1:ptr_nump, 1:ptr_nump) = 0.0_dp
     bb2(1:ptr_nump, 1:ptr_nump) = 0.0_dp

     ! bb0 = In
     DO j = 1, ptr_nump
        bb0(j, j) = 1.0_dp
     END DO

     ! bb1
     bb1 (1, 1) = 1.0_dp 
     bb1 (1, 2) = -1.0_dp 
     bb1 (ptr_nump, ptr_nump) = 1.0_dp
     bb1 (ptr_nump, ptr_nump - 1) = -1.0_dp

     DO j = 2, ptr_nump - 1
        bb1 (j, j) = 2.0_dp
        bb1 (j, j - 1) = -1.0_dp
        bb1 (j, j + 1) = -1.0_dp
     END DO

     ! bb2
     bb2 (1, 1) = 1.0_dp
     bb2 (1, 2) = -2.0_dp    
     bb2 (1, 3) = 1.0_dp

     bb2 (ptr_nump, ptr_nump) = 1.0_dp
     bb2 (ptr_nump, ptr_nump - 1) = -2.0_dp    
     bb2 (ptr_nump, ptr_nump - 2) = 1.0_dp

     bb2 (2, 1) = -2.0_dp
     bb2 (2, 2) = 5.0_dp    
     bb2 (2, 3) = -4.0_dp
     bb2 (2, 4) = 1.0_dp

     bb2 (ptr_nump - 1, ptr_nump) = -2.0_dp
     bb2 (ptr_nump - 1, ptr_nump - 1) = 5.0_dp    
     bb2 (ptr_nump - 1, ptr_nump - 2) = -4.0_dp
     bb2 (ptr_nump - 1, ptr_nump - 3) = 1.0_dp

     DO j = 3, ptr_nump - 2
        bb2 (j, j) = 6.0_dp
        bb2 (j, j - 1) = -4.0_dp
        bb2 (j, j - 2) = 1.0_dp
        bb2 (j, j + 1) = -4.0_dp
        bb2 (j, j + 2) = 1.0_dp
     END DO

     ! get weighted sum of the 0, 1st, 2nd regularization
     bb_all(1:ptr_nump, 1:ptr_nump) = (ptr_w0 * bb0(1:ptr_nump, 1:ptr_nump) &
          + ptr_w1 * bb1(1:ptr_nump, 1:ptr_nump) + ptr_w2 * &
          bb2(1:ptr_nump, 1:ptr_nump)) / (ptr_w0 + ptr_w1 + ptr_w2)

     CALL choldc(bb_all, ptr_nump, ptr_nump, parr, pge_error_status)
     IF ( pge_error_status == pge_errstat_error) THEN
        WRITE(*, *) modulename, ' : CHOLDC failed!!!'
        RETURN
     END IF

     ptr_b(1:ptr_nump, 1:n_fitvar_rad) = 0.0_dp

	 ! get L^T not L
     DO i = 1, ptr_nump
        DO j = 1, i - 1
           ptr_b(j, i + ozfit_start_index - 1) = bb_all (i, j)
        END DO
        ptr_b(i, i + ozfit_start_index - 1) = parr(i)
     END DO
  ELSE IF (ptr_order == 6) THEN
     
     CALL inverse(sa, sainv, nlay_fit)
     !IF ( pge_error_status == pge_errstat_error) THEN
     !   WRITE(*, *) modulename, ' : Singular Matrix!!!'
     !   RETURN
     !END IF

     CALL choldc(sainv, nlay_fit, nlay_fit, parr, pge_error_status)
     IF ( pge_error_status == pge_errstat_error) THEN
        WRITE(*, *) modulename, ' : CHOLDC failed!!!'
        RETURN
     END IF
     
     DO i =1, nlay_fit
        ptr_b(i, i + ozfit_start_index - 1) = parr(i)
        DO j = 1, i - 1
           ptr_b(j, i + ozfit_start_index - 1) = sainv(i, j)
        ENDDO
     ENDDO
     
  END IF
  
  RETURN
END SUBROUTINE get_reg_matrix
