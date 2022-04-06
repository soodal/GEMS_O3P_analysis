! Purpose : GEMS SNR Calculation, converted from read_gemssnr_evol.pro
! Author: Juseon Bak
! date  : 2020-11-10
! update: Dae Sung Choi
! update date : 2020-11-27


MODULE m_get_snr_gems
   character(len=100), parameter :: data_dir = '/data1/gems/o3p/ds/GEMS_O3P_analysis/SNR/'
   integer, parameter :: cr_lines=1001, ir_lines=1001, sr_lines =1001,sl_lines=1001
   real (kind=8), allocatable, dimension (:), save :: cr_wave, deg_bol, deg_eol  
   real (kind=8), allocatable, dimension (:), save :: ir_wave, ir_nom, ir_max, ir_num 
   real (kind=8), allocatable, dimension (:), save :: sr_wave, sr_cal_pri, sr_resp_pri, &
                                                      sr_cal_red, sr_resp_red
   real (kind=8), allocatable, dimension (:), save :: sl_wave, sl_per
   public :: get_snr_gems
  
   private
   CONTAINS
   ! wave(1) should be 300, wave should be 0.2 nm intervals, wave(nw) is free.
   ! rad is in the unit of photon
   ! 1/snr *100 = relative measurement errors (%).
   SUBROUTINE get_snr_gems (nw, wave, rad, snr)
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: nw
       REAL (KIND=8), INTENT(IN), DIMENSION(nw) :: rad, wave
       REAL (KIND=8), INTENT(OUT), DIMENSION(nw) :: snr
       ! local
       LOGICAL :: error
       REAL (KIND=8), PARAMETER :: pi    = 3.14159265358979d0
       REAL (KIND=8), PARAMETER :: n_bin = 1.0 !4.0
       ! imager system parameters
       REAL (KIND=8), PARAMETER :: T_i=70.50, EW_dis=7.722, N_f=34.0, N_bl=2.0*n_bin, N_bp=1.0, B=14.0
       ! Optical system parameter
       REAL (KIND=8), PARAMETER :: lambl=300.0, lambn=500.0, &
                                   res_samp=0.1985, Spec_SW=51.3, F_n=3.00, IFOV=67.21, nr=7.0, nt=3.0, lgr=325.01
       ! Detector subsystem parameter
       REAL (KIND=8), PARAMETER :: N_l      = 2048d0    !# of pixels along spatial dimension [#]
       REAL (KIND=8), PARAMETER :: N_p      = 1032d0    !# of pixels along wavelength dimension [#]
       REAL (KIND=8), PARAMETER :: d_x      = 18.0d0    !Pixel length (Spatial dimension) [um]
       REAL (KIND=8), PARAMETER :: d_l      = 18.0d0    !Pixel width (wavelength dimension) [um]
       REAL (KIND=8), PARAMETER :: TdarkRef = -20.0d0   !Dark Current Reference Temperature [degreeC]
       REAL (KIND=8), PARAMETER :: Ddark    = 12627.0d0 !EOL dark current density at referencetemperature [e-/pix/s]
       REAL (KIND=8), PARAMETER :: T_r      = 2.13d0    !Frame transfer time [ms]
       REAL (KIND=8), PARAMETER :: conv     = 1.0d0     !charge to voltage conversion factor[uV/e-]
       REAL (KIND=8), PARAMETER :: N_sat    = 263.0d0   !Saturation signal capacity [ke-]
       REAL (KIND=8), PARAMETER :: T        = -20.0d0   !CCD operating temperature [dC]
       REAL (KIND=8), PARAMETER :: m        = 4.0d0     !Number of output taps
       REAL (KIND=8), PARAMETER :: n_r      = 51.0d0    !read noise [e-]
       REAL (KIND=8), PARAMETER :: n_a      = 24d0*0.0d0 !CDS & pre-amp noise [e-]
       REAL (KIND=8), PARAMETER :: N_cds    = 0.0d0     !Number of CDS resets
       REAL (KIND=8), PARAMETER :: J_d = Ddark*2.0d0**((T-TdarkRef)/6.0d0) !EOL Dark current density at operational temperature [e-/pix/s]
       REAL (KIND=8), PARAMETER :: n_q = 1.0d0/((12.0d0)**(1d0/2d0))*(1000.0d0*N_sat)/(2.0d0**(B)*0.95d0) !Unbinned quantization noise [e-]
       !------------- CCD parameters --------------------------------------------
       REAL (KIND=8), PARAMETER ::perc_pix_illu = 100.0d0  !% pixels illuminated[%]
       REAL (KIND=8), PARAMETER ::CTE           = 0.99993d0 !CTE
       REAL (KIND=8), PARAMETER ::MaxN_PixTrans = 2024.00d0 !Max number of Pixel transfers [split frame]  
       REAL (KIND=8), PARAMETER ::r_t           = 51.05d0   !readout time [ms]
       REAL (KIND=8), PARAMETER ::RTSR          = 74.0d0!16.37*exp(0.0747*T)*6240000*0.00000324 Random Telegraph Signal Rate [e-/pix/s]
       !--------------- Satellite Viewing Gemetry ------------------------------
       REAL (KIND=8), PARAMETER ::Re    = 6378.0d0          !Radius of Earth [km]
       REAL (KIND=8), PARAMETER ::Ro    = 35786d0           !Geostationary Altitude of Satellite [km]
       REAL (KIND=8), PARAMETER ::Set_time = 0.1d0          !Step mirror settle time [s]
       REAL (KIND=8), PARAMETER ::f0    = 37.3d0            !Image point in Latitude [deg]
       REAL (KIND=8), PARAMETER ::f_image  = 1000.0d0*d_x/IFOV       !ocal length of imaging system [mm]
       ! calculated variables
       REAL (KIND=8) :: a, phi0, SL0, theta0, GSD_ew, GSD_ns, EGSD_ew, EGSD_ns
       REAL (KIND=8) :: nl, np, Spec_dis, Spec_bp, Spec_sp, theta_p, SFOV_l, SFPV_w, &
                        OMIS, PS, Sph_abb_cor, D_ta, D_area, TA_area, Omega, A_omega
       REAL (KIND=8), DIMENSION(nw) :: TotSig_pix, Pixel_sat, Pix_dkph_noise, Out_sig,smear_noise,sat
       REAL (KIND=8) :: phFl_pix, TotDk_pix, CTE_noise, Out_noise, SignleF_SNR
       INTEGER :: idx_400, ntmp,i, nw_data

       error = .false.
       ! read constants -----------------------------------------------------------
       CALL read_gems_snr_data (nw,wave(1), wave(nw), error)
       IF (error) then 
         WRITE(*,*) 'get_snr_gems: error in read_gems_snr_data'; stop
       endif
       !--------------- Calculated Geometric Parameters -------------------------
       a       = Re+Ro !Re+h [km]
       phi0    = f0    !Earth angle to image point [deg]
       SL0     = dsqrt(a**(2.0d0) + Re**(2.0d0) - 2.0d0*a*Re*dcos(pi*phi0/180.0d0))  !Slant range to image center [km]
       theta0  = (180.0d0/pi)*datan((dsin(pi*phi0/180.0d0))/(a/Re-dcos(pi*phi0/180.0d0)))
                 !Nominal sensor viewing angle (northward from nadir) [deg]
       GSD_ew  = SL0*(0.001d0*dmin1(Spec_SW, d_l)/f_image) !GSD at the image point (nominally along lat. line) [km]
       GSD_ns  = SL0*IFOV*0.000001d0/(dcos(pi*(phi0+theta0)/180.0d0))!GSD at the image point (nominally along lon. line) [km]
       EGSD_ew = GSD_ew                                           !Effective GSD [km]
       EGSD_ns = GSD_ns * N_bl                                    !Effective GSD [km]

      !--------------- Calculated sensor parameters ----------------------------
       nl       = (lambn-lambl)/res_samp           !number of wavelength bands [#]
       np       = nl                               !pixels per readout line (lambda dimension) [#]
       Spec_dis = res_samp/(0.001d0*d_l)             !spectral dispersion [nm/mm]
       Spec_bp  = 0.599d0                            !spectral bandpass [nm]
       Spec_sp  = Spec_bp/res_samp                 !spectral sampling [no/FWHM]
       theta_p  = 0.5d0*(180.0d0/pi)*(IFOV*0.000001d0)  !IFOV for single pixel - 1/2 angle
       SFOV_l   = 2.0d0*theta_p*N_l                  !Slit FOV - in slit length direction [deg]
       SFPV_w   = (180.0d0/pi)*(Spec_sw/f_image/IFOV)     !Slit FOV - in slit width direction [deg] 
       OMIS     = (0.000001d0*d_x) / (1000.0d0*GSD_ns) !Optical Magnification of Imaging System [m/m]
       PS       = 1000.0d0 * GSD_ew / d_x            !Plate Scale [km/mm]
       Sph_abb_cor = dcos(datan(0.5d0*OMIS))           !Spherical abberation correction
       D_ta     = (f_image/10.0d0)/F_n               !Telescope aperture diameter [cm]
       D_area   = d_x * d_l * 0.00000001d0           !detector area [cm^2]
       TA_area  = pi/4.0d0*(D_ta**(2.0d0))               !Telescope aperture area [cm^2]
       Omega    = 0.000001d0*(d_x/f_image)*dmin1(d_l,Spec_sw)/f_image !Solid angle view for single pixel [st]
       A_Omega  = TA_area * Omega                  !Optical throughput (etendue) [cm^2sr]
       !-------------- Calculate radiometric model --------------------------
       !op_trans0  = om_OST * 100.0 * (1.0 - degrad / 100.0)
       !grat_eff0  = om_GEC * 100.0
       !sys_trans0 = op_trans0 * grat_eff0 / 100.0
       !op_trans  = interpol(op_trans0, om_wave,    wave)
       !grat_eff  = spline  (om_wave,   grat_eff0,  wave)
       !sys_trans = spline  (om_wave,   sys_trans0, wave)

       !int_rad = interpol(rad, wave, sr_wave, /spline)
       DO i = 1, nw
         phFl_pix   = rad(i) / (sr_cal_red(i) * (1.0 + deg_eol(i)/100.0))
         !phFl_pix  = ir_nom / (sr_cal_red * (1 + deg_eol/100.0d))
         TotSig_pix(i) = phFl_pix * (T_i / 1000.0d0)
         TotDk_pix  = (J_d * (T_i + T_r)/1000.0d0)  + & 
                      (RTSR * (T_i + T_r * ir_num(i) /1024.0d0) / 1000.0d0) + &
                      (J_d * R_t / 1000.0d0 * ir_num(i) / 1024.0d0) + & 
                      (RTSR * R_t / 1000.0d0 *ir_num(i) / 1024.0d0)
         Pixel_sat(i)  = 100.0d0 * (TotSig_pix(i) * (1.0d0+sl_per(i)/100.0d0) + TotDk_pix) / (1000.0d0*N_sat)
         Pix_dkph_noise(i) = dsqrt(TotDk_pix + TotSig_pix(i) * (1.0d0 + sl_per(i)/100.0d0))
         Out_sig(i)        = TotSig_pix(i) * N_bp * N_bl        !Output signal [e-]
       ENDDO

       !idx_400 = where(wave ge 400., n_idx_400) !IDL
       idx_400 = COUNT(mask=(wave < 400) )
       !if (n_idx_400 /= 0) then  !IDL
       if (idx_400 /= nw) then 
         ! @ < 400
         ntmp = idx_400
         smear_noise(1:idx_400) = &
          dsqrt((sum(TotSig_pix(1:idx_400))/ntmp + &
                 sum(TotSig_pix(1:idx_400)* sl_per(1:idx_400)/100.0d0)/ntmp) * (perc_pix_illu/100.0d0) * T_r / (T_r + T_i))
         ! @ > 400
         ntmp = nw - idx_400 
         smear_noise(idx_400+1:nw) = &
            dsqrt((sum(TotSig_pix(idx_400+1:nw))/ntmp + & 
                   sum(TotSig_pix(idx_400+1:nw)*sl_per(idx_400+1:nw)/100.0d0)/ntmp) * (perc_pix_illu/100.0d0) * T_r /(T_r + T_i))
       else
         smear_noise(1:nw) = &
            dsqrt((sum(TotSig_pix)/nw + sum(TotSig_pix * sl_per/100.0d0)/nw) * (perc_pix_illu/100.0d0) * T_r / (T_r + T_i))
       endif
          
       DO i = 1, nw 
         CTE_noise = dsqrt(2.0d0 * (1.0d0 - CTE) * &
                  (TotSig_pix(i) + TotSig_pix(i) * sl_per(i)/100.0d0 + smear_noise(i))**(0.715d0) * MaxN_PixTrans)
         Out_noise = dsqrt(N_bp*N_bl*N_cds*n_a + N_bp*N_bl*n_r**(2.0d0) + &
                  N_bp*N_bl*n_q**(2.0d0) + N_bp*N_bl*Pix_dkph_noise(i)**(2.0d0) + &
                  N_bp*N_bl*smear_noise(i)**(2.0d0) + N_bp*N_bl*CTE_noise**(2.0d0))
         !Out_noise   = dsqrt(n_a + n_r^2.0 + n_q^2.0 + N_bl*N_bp*(Pix_dkph_noise^2.0) +&
         !N_bl*N_bp*(smear_noise^2.0) + N_bl*N_bp*(CTE_noise^2.0))

         SignleF_SNR = Out_sig(i) / Out_noise
         snr(i)         = SignleF_SNR *dsqrt(N_f)
         sat(i)         = Pixel_sat(i) ![%]
       ENDDO

   END SUBROUTINE

   SUBROUTINE read_gems_snr_data(nw, fwav, lwav, error)
       IMPLICIT NONE
       integer, intent(in) :: nw
       real (kind=8), intent(in) :: fwav, lwav
       logical, intent(inout) :: error
       ! local
       character(len=100)::  in_file
       integer       :: i ,c
       real (kind=8) :: wav, tmp1, tmp2, tmp3
       logical, save :: first = .true.
       logical :: file_exist
       if (first) then
        !allocate (cr_wave(cr_lines), deg_bol(cr_lines), deg_eol(cr_lines))
        !allocate (ir_wave(ir_lines), ir_nom(ir_lines), ir_max(ir_lines),ir_num(ir_lines))
        !allocate (sr_wave(sr_lines), sr_cal_pri(sr_lines), sr_resp_pri(sr_lines), &
        !          sr_cal_red(sr_lines), sr_resp_red(sr_lines))
        !allocate (sl_wave(sl_lines), sl_per(sl_lines))
        allocate (cr_wave(nw), deg_bol(nw), deg_eol(nw))
        allocate (ir_wave(nw), ir_nom(nw),  ir_max(nw), ir_num(nw))
        allocate (sr_wave(nw), sr_cal_pri(nw), sr_resp_pri(nw), sr_cal_red(nw), sr_resp_red(nw))
        allocate (sl_wave(nw), sl_per(nw))
        in_file = adjustl(trim(data_dir))//'CONTAMINATION_RADIATION.txt'
        inquire(file=in_file, exist = file_exist)
        if (.not. file_exist) then 
          write(*,*) 'read_gems_snr_data: not exsit for'//adjustl(trim(in_file)) ;stop
        endif 
        open(unit=1, file=adjustl(trim(in_file)), status='old')
        c = 0
        do i = 1, cr_lines 
          read(1,*) wav, tmp1
          if ( wav >= fwav .and. wav <= lwav) then
           c = c + 1
           cr_wave(c)=wav; deg_eol(c)=tmp1
          endif
        enddo
        close (1)
        IF (c /= nw) then 
           write(*,*) 'read_gems_snr_data: inconsistent cr_wave'
           error = .true. ; return
        endif
        ! input radiance
        in_file = adjustl(trim(data_dir))//'INPUT_RADIANCE.txt'
        inquire(file=in_file, exist = file_exist)
        if (.not. file_exist) then 
          write(*,*) 'read_gems_snr_data: not exsit for'//adjustl(trim(in_file)) ;stop
        endif 
        open(unit=1, file=adjustl(trim(in_file)), status='old')
        c = 0
        do i = 1, ir_lines
          read(1,*) wav, tmp1, tmp2, tmp3
          if (wav >= fwav .and. wav <= lwav) then
           c = c + 1
           ir_wave(c) = wav;  ir_nom(c)=tmp1 ;  ir_max(c)=tmp2; ir_num(c)=tmp3
          endif
        enddo
        close (1)
        IF (c /= nw) then 
           write(*,*) 'read_gems_snr_data: inconsistent ir_wave'
           error = .true. ; return
        endif
        ! systematic responsivity
        in_file = adjustl(trim(data_dir))//'SYSTEM_RESPONSIVITY.txt'
        inquire(file=in_file, exist = file_exist)
        if (.not. file_exist) then 
         write(*,*) 'read_gems_snr_data: not exsit for'//adjustl(trim(in_file)) ;stop
        endif 
        open(unit=1, file=adjustl(trim(in_file)), status='old')
        c= 0
        do i = 1, sr_lines
          read(1,*) wav, tmp1, tmp2
          if (wav >= fwav .and. wav <= lwav) then
          c = c + 1
          sr_wave(c)=wav; sr_cal_pri(c)=tmp1; sr_cal_red(c)=tmp2
          endif
        enddo
        close (1)
        IF (c /= nw) then 
           write(*,*) 'read_gems_snr_data: inconsistent sr_wave'
           error = .true. ; return
        endif
        sr_resp_pri = 1.0/sr_cal_pri
        sr_resp_red = 1.0/sr_cal_red
        ! straylight
        in_file = adjustl(trim(data_dir))//'STRAY_LIGHT.txt'
        inquire(file=in_file, exist = file_exist)
        if (.not. file_exist) then 
         write(*,*) 'read_gems_snr_data: not exsit for'//adjustl(trim(in_file)) ;stop
        endif 
        open(unit=1, file=adjustl(trim(in_file)), status='old')
        c = 0
        do i = 1, sl_lines
          read(1,*) wav, tmp1
          if (wav >= fwav .and. wav <= lwav) then
          c = c + 1
          sl_wave(c)=wav; sl_per(c)=tmp1    
          endif
        enddo
        close (1)
        IF (c /= nw) then 
           write(*,*) 'read_gems_snr_data: inconsistent sl_wave'
           error = .true. ; return
        endif
        ! straylight
        first = .false.
       endif
 END SUBROUTINE
END MODULE m_get_snr_gems
