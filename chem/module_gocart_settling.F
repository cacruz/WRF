MODULE MODULE_GOCART_SETTLING
!
! Original GOCART
! KYKang, MKlose, CLWu various changes 
! YShao [2011/09/30] complete update
!
! 11/29/18 - A. Ukhov, bug fix: dust (and sea salt) mass balance was
! violated in the "settling" routine 
!
! NUWRF: Dec 2021 - ZTao
!   Add in subroutines to consider settling of sulf, NH4 and NO3 aerosols
!   The code was adapted from the relevant components of GEOS-5  

CONTAINS

SUBROUTINE gocart_settling_driver(dt,config_flags,t_phy,moist,      &
           chem,rho_phy,dz8w,p8w,p_phy,                             &
           dustin,seasin,dx,g,                                      &
           dustgraset_1,dustgraset_2,dustgraset_3,                  & 
           dustgraset_4,dustgraset_5,                               & 
           setvel_1,setvel_2,setvel_3,setvel_4,setvel_5,imod,       & 
           ids,ide,jds,jde,kds,kde,                                 &
           ims,ime,jms,jme,kms,kme,                                 &
           its,ite,jts,jte,kts,kte)

  USE module_configure
  USE module_state_description
  USE module_data_gocart_dust
  USE module_data_gocart_seas
  USE module_model_constants, ONLY: mwdry
  IMPLICIT NONE

  TYPE(grid_config_rec_type), INTENT(IN) :: config_flags

  INTEGER, INTENT(IN) :: ids,ide, jds,jde, kds,kde,                  &
                         ims,ime, jms,jme, kms,kme,                  &
                         its,ite, jts,jte, kts,kte

  REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_moist), INTENT(IN) :: moist
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme,num_chem), INTENT(INOUT) :: chem
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: t_phy,p_phy,dz8w,p8w,rho_phy
  REAL,  DIMENSION( ims:ime , jms:jme, 5 ),                          &
          INTENT(IN   ) ::  dustin,seasin  

  INTEGER, INTENT(IN) :: imod          
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(INOUT)::                  &
        dustgraset_1,dustgraset_2,dustgraset_3,                      &
        dustgraset_4,dustgraset_5,                                   &
        setvel_1,setvel_2,setvel_3,setvel_4,setvel_5 
                       
  REAL*8,  DIMENSION (1,1,5)          :: graset_dust,grasetvel_dust 
  REAL*8,  DIMENSION (1,1,4)          :: graset_ss,grasetvel_ss

  REAL, INTENT(IN) :: dt,dx,g
  INTEGER          :: kkk,nmx,i,j,k,kk,lmx,iseas,idust
  REAL*8, DIMENSION (1,1,kte-kts+1)   :: tmp,airden,p_mid,delz,rh
  REAL*8, DIMENSION (1,1,kte-kts+1,5) :: ddust
  REAL*8, DIMENSION (1,1,kte-kts+1,4) :: sea_salt

!NUWRF added variables
  INTEGER         :: iflag                                          !iflag used for control particle swelling
  REAL, PARAMETER :: reff_sulf = 0.2659e-6, reff_nh4 = 0.2659e-6    !effective dry particle radium (m)
  REAL, PARAMETER :: reff_no3an1 = 0.2659e-6, reff_no3an2 = 2.1e-6
  REAL, PARAMETER :: reff_no3an3 = 7.57e-6

  REAL, PARAMETER :: rho_sulf = 1769., rho_nh4 = 1769.        !dry particle density (kg/m^3)
  REAL, PARAMETER :: rho_no3an1 = 1725., rho_no3an2 = 2200.
  REAL, PARAMETER :: rho_no3an3 = 2650.

  REAL, PARAMETER :: mwso4 = 96.06, mwnh4 = 18.04, mwair = 28.97  !molecular weight (g/mol)
  REAL*8, DIMENSION(1,1,kte-kts+1) :: tc_species, delzz
  
  INTEGER :: uoc_flag  ! flag for UoC dust schemes 
!
!   REAL*8, DIMENSION (5), PARAMETER :: den_dust(5)=(/2500.,2650.,2650.,2650.,2650./)
!   REAL*8, DIMENSION (5), PARAMETER :: reff_dust(5)=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/)
!   REAL*8, DIMENSION (4), PARAMETER :: den_seas(4)=(/2200.,2200.,2200.,2290./)
!   REAL*8, DIMENSION (4), PARAMETER :: reff_seas(4)=(/0.30D-6,1.00D-6,3.25D-6,7.50D-6/)
!  
  REAL*8 conver, converi
  conver=1.e-9
  converi=1.e9
    
  uoc_flag = 0
  if (config_flags%dust_opt .eq. 4) uoc_flag = 1

  lmx=kte-kts+1

  do j=jts,jte
    do i=its,ite

       IF(config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11) THEN
         do kkk=1,5
            ddust(1,1,kts,kkk)=dustin(i,j,kkk)*conver
         enddo
         kk=0
         do k=kts,kte
            kk=kk+1
            
            p_mid(1,1,kk) =.01*p_phy(i,kte-k+kts,j)
            delz(1,1,kk)  =dz8w(i,kte-k+kts,j)
            airden(1,1,kk)=rho_phy(i,k,j)
            tmp(1,1,kk)= t_phy(i,k,j)
            rh(1,1,kk) = .95
            rh(1,1,kk) = MIN( .95, moist(i,k,j,p_qv) / &
                 (3.80*exp(17.27*(t_phy(i,k,j)-273.)/ &
                 (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
            rh(1,1,kk) = max(1.0D-1,rh(1,1,kk))
         enddo
       ELSE  
         kk=0       
         DO k=kts,kte
            kk=kk+1
!
            ddust(1,1,kk,1)=chem(i,k,j,p_dust_1)                         ! chem and dust in [ug/kg]
            ddust(1,1,kk,2)=chem(i,k,j,p_dust_2)
            ddust(1,1,kk,3)=chem(i,k,j,p_dust_3)
            ddust(1,1,kk,4)=chem(i,k,j,p_dust_4)
            ddust(1,1,kk,5)=chem(i,k,j,p_dust_5)
            
            p_mid(1,1,kk)=.01*p_phy(i,kte-k+kts,j)
            delz(1,1,kk)=dz8w(i,kte-k+kts,j)                            ! delz(1) = dz8w(kte), delz(lmx)=dz8w(kts)
            airden(1,1,kk)=rho_phy(i,k,j)
            tmp(1,1,kk)   =t_phy(i,k,j)
            rh(1,1,kk)    = .95
            rh(1,1,kk)    = MIN( .95, moist(i,k,j,p_qv) /               &
                 (3.80*exp(17.27*(t_phy(i,k,j)-273.)/                   &
                 (t_phy(i,k,j)-36.))/(.01*p_phy(i,k,j))))
            rh(1,1,kk)=max(1.0D-1,rh(1,1,kk))
         ENDDO
       ENDIF
       graset_dust(1,1,:)=0.
       graset_ss(1,1,:)=0.                                           

       iseas=0
       idust=1
       CALL settling(1,1,lmx,5,g,dyn_visc,ddust,tmp,p_mid,delz,       &
                     imod,graset_dust,grasetvel_dust, uoc_flag,       & 
                     den_dust,reff_dust,dt,rh,idust,iseas,airden)
          
       IF (config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) THEN
          kk=1
          do kkk=1,5
             if (kkk .le. 4) sea_salt(1,1,kts,kkk)=seasin(i,j,kkk)*conver
             if(ddust(1,1,kk,kkk) .ge. dustin(i,j,kkk)) ddust(1,1,kk,kkk)=dustin(i,j,kkk)
          enddo
          chem(i,kts,j,p_p25i)=chem(i,kts,j,p_p25i) &
                       +.25*(ddust(1,1,kk,1)+.286*ddust(1,1,kk,2))*converi
          chem(i,kts,j,p_p25i)=max(chem(i,kts,j,p_p25i),1.e-16)
          chem(i,kts,j,p_p25j)=chem(i,kts,j,p_p25j) &
                       +.75*(ddust(1,1,kk,1)+.286*ddust(1,1,kk,2))*converi
          chem(i,kts,j,p_p25j)=max(chem(i,kts,j,p_p25j),1.e-16)
          chem(i,kts,j,p_soila)=chem(i,kts,j,p_soila) &
                       +(.714*ddust(1,1,kk,2)+ddust(1,1,kk,3))*converi
          chem(i,kts,j,p_soila)=max(chem(i,kts,j,p_soila),1.e-16)
       ELSE                     
          kk = 0
          DO k = kts,kte
             kk = kk+1
             chem(i,k,j,p_dust_1)=ddust(1,1,kk,1)                         ! dust for size bin 1 [ug/kg]
             chem(i,k,j,p_dust_2)=ddust(1,1,kk,2)                         ! ...
             chem(i,k,j,p_dust_3)=ddust(1,1,kk,3)                         ! ...
             chem(i,k,j,p_dust_4)=ddust(1,1,kk,4)                         ! ...
             chem(i,k,j,p_dust_5)=ddust(1,1,kk,5)                         ! dust for size bin 5 (dust_opt 3: for all size bins) [ug/kg]
             
             sea_salt(1,1,kk,1)=chem(i,k,j,p_seas_1)                    ! salt [ug/kg]
             sea_salt(1,1,kk,2)=chem(i,k,j,p_seas_2)
             sea_salt(1,1,kk,3)=chem(i,k,j,p_seas_3)
             sea_salt(1,1,kk,4)=chem(i,k,j,p_seas_4)              
          ENDDO
       ENDIF
!
! gravitional settling in [ug/m2/s]; from settling, graset_dust in [ug/kg][m/s]
!
       dustgraset_1(i,j)=dustgraset_1(i,j)+conver*graset_dust(1,1,1)*airden(1,1,1)*dt   !dustgraset [kg/m2]
       dustgraset_2(i,j)=dustgraset_2(i,j)+conver*graset_dust(1,1,2)*airden(1,1,1)*dt
       dustgraset_3(i,j)=dustgraset_3(i,j)+conver*graset_dust(1,1,3)*airden(1,1,1)*dt
       dustgraset_4(i,j)=dustgraset_4(i,j)+conver*graset_dust(1,1,4)*airden(1,1,1)*dt
       dustgraset_5(i,j)=dustgraset_5(i,j)+conver*graset_dust(1,1,5)*airden(1,1,1)*dt

       setvel_1(i,j)=grasetvel_dust(1,1,1)                            ! settling velocity [m/s]
       setvel_2(i,j)=grasetvel_dust(1,1,2)
       setvel_3(i,j)=grasetvel_dust(1,1,3)
       setvel_4(i,j)=grasetvel_dust(1,1,4)
       setvel_5(i,j)=grasetvel_dust(1,1,5)

       iseas=1
       idust=0

       CALL settling(1,1,lmx,4,g,dyn_visc,sea_salt,tmp,p_mid,delz,        &
                     imod,graset_ss,grasetvel_ss, uoc_flag,               &
                     den_seas,reff_seas,dt,rh,idust,iseas,airden)
       IF (config_flags%chem_opt == 2 .or. config_flags%chem_opt == 11 ) THEN
          kk=1
          do kkk=1,4
             if(sea_salt(1,1,kk,kkk) .ge. seasin(i,j,kkk))sea_salt(1,1,kk,kkk)=seasin(i,j,kkk)
          enddo
          chem(i,kts,j,p_naai)=chem(i,kts,j,p_naai) &
                       +.25*(sea_salt(1,1,kk,1)+.942*sea_salt(1,1,kk,2))*converi
          chem(i,kts,j,p_naai)=max(1.e-16,chem(i,kts,j,p_naai))
          chem(i,kts,j,p_naaj)=chem(i,kts,j,p_naaj) &
                       +.75*(sea_salt(1,1,kk,1)+.942*sea_salt(1,1,kk,2))*converi
          chem(i,kts,j,p_naaj)=max(1.e-16,chem(i,kts,j,p_naaj))
          chem(i,kts,j,p_seas)=chem(i,kts,j,p_seas) &
                       +(.058*sea_salt(1,1,kk,2)+sea_salt(1,1,kk,3))*converi
          chem(i,kts,j,p_seas)=max(1.e-16,chem(i,kts,j,p_seas))
       ELSE
          kk=0
          DO k=kts,kte
             kk=kk+1
             chem(i,k,j,p_seas_1)=sea_salt(1,1,kk,1)
             chem(i,k,j,p_seas_2)=sea_salt(1,1,kk,2)
             chem(i,k,j,p_seas_3)=sea_salt(1,1,kk,3)
             chem(i,k,j,p_seas_4)=sea_salt(1,1,kk,4)
          ENDDO
       ENDIF 
       
! NUWRF: begin process of NH4/NO3 settling (ZTao)
!! 1) sulf (ppmv -> ug/kg)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        tc_species(1,1,kk) = 1.e3 * chem(i,k,j,p_sulf) * mwso4/mwair
        delzz(1,1,kk) = dz8w(i,k,j)
      end do
      iflag = 3                                                         !ammonium sulfate

      call settling_nitrate(1,1,1,1,lmx, iflag, reff_sulf, rho_sulf,  &
                             dt, tc_species, tmp, airden, rh, delzz)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        chem(i,k,j,p_sulf) = 1.e-3 * tc_species(1,1,kk) * mwair/mwso4   !ug/kg -> ppmv
      end do

!! 2) NH4 (ppmv -> ug/kg)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        tc_species(1,1,kk) = 1.e3 * chem(i,k,j,p_nh4) * mwnh4/mwair
      end do
      iflag = 3

      call settling_nitrate(1,1,1,1,lmx, iflag, reff_nh4, rho_nh4,    &
                            dt, tc_species, tmp, airden, rh, delzz)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        chem(i,k,j,p_nh4) = 1.e-3 * tc_species(1,1,kk) * mwair/mwnh4    !ug/kg -> ppmv
      end do

!! 3) NO3an1 (ug/kg)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        tc_species(1,1,kk) = chem(i,k,j,p_no3an1) 
      end do
      iflag = 3                                                         !like (NH4)SO4

      call settling_nitrate(1,1,1,1,lmx, iflag, reff_no3an1, rho_no3an1,  &
                               dt, tc_species, tmp, airden, rh, delzz)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        chem(i,k,j,p_no3an1) = tc_species(1,1,kk)
      end do

!! 4) NO3an2 (ug/kg)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        tc_species(1,1,kk) = chem(i,k,j,p_no3an2) 
      end do
      iflag = 2                                                          !like sea salt

      call settling_nitrate(1,1,1,1,lmx, iflag, reff_no3an2, rho_no3an2,  &
                               dt, tc_species, tmp, airden, rh, delzz)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        chem(i,k,j,p_no3an2) = tc_species(1,1,kk)
      end do

!! 5) NO3an3 (ug/kg)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        tc_species(1,1,kk) = chem(i,k,j,p_no3an3) 
      end do
      iflag = 0                                                          !like dust

      call settling_nitrate(1,1,1,1,lmx, iflag, reff_no3an3, rho_no3an3,  &
                               dt, tc_species, tmp, airden, rh, delzz)
      kk = 0
      do k = kts, kte
        kk = kk + 1
        chem(i,k,j,p_no3an3) = tc_species(1,1,kk)
      end do

    enddo  ! i
  enddo  ! j 


 END SUBROUTINE gocart_settling_driver


 subroutine settling(imx,jmx,lmx,nmx,g0,dyn_visc,tc,tmp,p_mid,delz, &
                    imod,graset, grasetvel, uoc,                    & 
                    den_in,reff_in,dt,rh,idust,iseas,airden)

! ****************************************************************************
! *                                                                          *
! *  Calculate the loss by settling, using an implicit method                *
! *                                                                          *
! *  Input variables:                                                        *
! *    SIGE(k)         - sigma coordinate of the vertical edges              *
! *    PS(i,j)         - Surface pressure (mb)                               *
! *    TMP(i,j,k)      - Air temperature  (K)                                *
! *    CT(i,j)         - Surface exchange coeff for moisture
! *                                                                          *
! **************************************************************************** 

  IMPLICIT  NONE

  INTEGER, INTENT(IN) :: imx, jmx, lmx, nmx,iseas,idust
  INTEGER             :: ntdt
  REAL,    INTENT(IN) :: dt,g0,dyn_visc
  REAL*8,  INTENT(IN) :: tmp(imx,jmx,lmx), delz(imx,jmx,lmx),  &
                         rh(imx,jmx,lmx), p_mid(imx,jmx,lmx),airden(imx,jmx,lmx)
! 
  REAL*8, INTENT(IN)  :: den_in(nmx), reff_in(nmx)

  REAL*8, INTENT(INOUT) :: tc(imx,jmx,lmx,nmx)

  INTEGER, INTENT(IN)   :: imod, uoc
  REAL*8, INTENT(INOUT) :: graset(imx,jmx,nmx)
  REAL*8, INTENT(OUT)   :: grasetvel(imx,jmx,nmx)
  
  REAL*8    :: den(nmx), reff(nmx)                    ! local variables here
  REAL*8    :: dt_settl(nmx), rcm(nmx), rho(nmx)
  INTEGER   :: ndt_settl(nmx)
  REAL*8    :: dzmin, vsettl, dtmax, pres, rhb, rwet(nmx), ratio_r(nmx)
  REAL*8    :: c_stokes, free_path, c_cun, viscosity, growth_fac
  REAL*8    :: vd_cor(lmx), vd_wk1
  INTEGER   :: k, n, i, j, l, l2

  REAL*8    :: transfer_to_below_level,temp_tc

! for sea-salt:
  REAL*8, PARAMETER :: c1=0.7674, c2=3.079, c3=2.573E-11, c4=-1.424 

! for OMP:
  REAL*8 :: rwet_priv(nmx), rho_priv(nmx)

! IF (type) /= 'dust' .AND. TRIM(aero_type) /= 'sea_salt') RETURN
  IF ( idust.ne.1 .and. iseas.ne.1 ) RETURN
  WHERE ( tc(:,:,:,:) < 0.0 ) tc(:,:,:,:) = 1.0D-32

  den = den_in  
  reff = reff_in  
  dzmin = MINVAL(delz(:,:,:))
  IF (idust == 1) growth_fac = 1.0
  IF (iseas == 1) growth_fac = 3.0
!
  IF (idust == 1 .and. uoc == 1) then
      den(1) = 2650.            ! constant density is use in UoC dust emission schemes; 
!                               dust radii are now consistent with GOCART dust schemes.  [mklose, 03082015]
  ENDIF

  DO k = 1, nmx                                                   ! k for different size bins
!
! Settling velocity (m/s) for each tracer (Stokes Law)
! DEN         density                        (kg/m3)
! REFF        effective radius               (m)
! dyn_visc    dynamic viscosity              (kg/m/s)
! g0          gravity                        (m/s2)
! 3.0         corresponds to a growth of a factor 3 of radius with 100% RH
! 0.5         upper limit with temp correction
!
     vsettl = 4.0/9.0 * g0 * den(k) * (growth_fac*reff(k))**2 / dyn_visc
!
! Determine the maximum time-step satisying CFL: dt <= (dz)_min / v_settl
!
     ntdt = INT(dt)
     dtmax = dzmin / vsettl
     ndt_settl(k) = MAX( 1,INT(ntdt/dtmax) )

! Limit maximum number of iterations
     IF (ndt_settl(k) > 12) ndt_settl(k) = 12
     dt_settl(k) = REAL(ntdt) / REAL(ndt_settl(k))

! Particles radius in centimeters
     IF (iseas.eq.1)rcm(k) = reff(k)*100.0
     IF (idust.eq.1) THEN
      rwet(k) = reff(k)
      ratio_r(k) = 1.0
      rho(k) = den(k)
     ENDIF
  ENDDO

! Solve the bidiagonal matrix (l,l)

!$OMP PARALLEL DO &
!$OMP DEFAULT( SHARED ) &
!$OMP PRIVATE( i,   j,   l,   l2, n,   k,   rhb, rwet_priv, ratio_r, c_stokes)&
!$OMP PRIVATE( free_path, c_cun, viscosity, rho_priv, vd_cor )

! Loop over latitudes

  DO j = 1,jmx                      ! lat loop
    DO i = 1,imx                    ! lon loop
      DO k = 1,nmx                  ! bin loop
        graset(i,j,k)=0.
        grasetvel(i,j,k)=0.
    
        IF (idust.eq.1) THEN
           rwet_priv(k) = rwet(k)
           rho_priv(k)  = rho(k)
        END IF

        DO n = 1,ndt_settl(k)        ! time loop 

          transfer_to_below_level=0

          DO l = lmx,1,-1            ! height loop, from top
             l2 = lmx - l + 1

             IF (iseas.eq.1) THEN
                rhb = MIN(9.9D-1, rh(i,j,l))                             ! Aerosol growth with relative humidity (Gerber, 1985)
                rwet_priv(k) = 0.01*(c1*rcm(k)**c2/(c3*rcm(k)**c4 -  &
                LOG10(rhb)) + rcm(k)**3)**0.33                           ! td changed to LOG10
                ratio_r(k) = (reff(k)/rwet_priv(k))**3.0
             END IF

             c_stokes = 1.458E-6*tmp(i,j,l)**1.5/(tmp(i,j,l) + 110.4)    ! Dynamic viscosity
             free_path = 1.1E-3/p_mid(i,j,l2)/SQRT(tmp(i,j,l))           ! Free path as func of pres(mb) and temp(K); order of p_mid: top->sfc

             c_cun = 1.0+free_path/rwet_priv(k)*                     &
             (1.257 + 0.4*EXP(-1.1*rwet_priv(k)/free_path))              ! Slip correction
             viscosity = c_stokes / c_cun                                ! Corrected dynamic viscosity (kg/m/s)

             IF (iseas.eq.1) THEN
                rho_priv(k) = ratio_r(k)*den(k) + (1.-ratio_r(k))*1000.
             END IF

             vd_cor(l) = 2./9.*g0*rho_priv(k)*rwet_priv(k)**2/viscosity  ! Settling velocity, depends on temp

             ! Update mixing ratio; order of delz: top->sfc
             temp_tc=tc(i,j,l,k)                                          ! temp_tc - for temporal storage [ug/kg]            
             vd_wk1 = dt_settl(k)*vd_cor(l)/delz(i,j,l2)                  ! fraction to leave level

             tc(i,j,l,k)   =  tc(i,j,l,k)*(1.- vd_wk1)+transfer_to_below_level          ! [ug/kg]
            
             IF (l==1) THEN
                graset(i,j,k)=graset(i,j,k)+vd_cor(l)*temp_tc/ndt_settl(k) ! [ug/kg][m/s]
                grasetvel(i,j,k)=vd_cor(l)                                 ! [m/s]
             ELSE
                transfer_to_below_level = (temp_tc*vd_wk1)*((delz(i,j,l2)*airden(i,j,l))/(delz(i,j,l2+1)*airden(i,j,l-1))) ! [ug/kg]
             ENDIF
          ENDDO                 !l, height
        ENDDO                   !n, time
      ENDDO                     !k, bin
!
    ENDDO                     !i
  END DO                      !j
!$OMP END PARALLEL DO

END SUBROUTINE settling

subroutine settling_nitrate(i1,i2,j1,j2,km, flag, radiusInp, rhopInp, &
                                cdt, tc, tmpu, rhoa, rh, dz)

! This subroutine was adapted from subroutine Chem_SettlingSimple
!   in GEOS-5 for NU-WRF application.
!
! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop) 
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the 
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! REVISION HISTORY:
!   12/27/2021   Z. Tao

   implicit none

   integer, intent(in) :: i1, i2, j1, j2, km
   integer, intent(in) :: flag                 ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt                  ! main dt (time_step in namelist) 
   real, intent(in)    :: radiusInp, rhopInp   ! particle effective radius (m) and density (kg/m^3)
   real*8, dimension(i1:i2,j1:j2,1:km), intent(in) :: tmpu, rhoa   ! air temperature (K) and density (kg/m^3) 
   real*8, dimension(i1:i2,j1:j2,1:km), intent(in) :: rh, dz       ! relative humidity and layer thickness (m)
   real*8, dimension(i1:i2,j1:j2,1:km), intent(inout) :: tc        ! mass mixing ratio (ug/kg)

! Local Variables
   integer  ::  i, j, k, iit
   integer  ::  nsubstep             ! # of sub steps in each time step
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]
   real :: vsettle(i1:i2,j1:j2,km)   ! fall speed [m s-1]
   real :: diff_coef                 ! Brownian diffusion coefficient [m2 s-1]
   real :: qdel(i1:i2,j1:j2), qsrc(i1:i2,j1:j2)

! The following parameters relate to the swelling of seasalt like particles
!    following Fitzgerald, Journal of Applied Meteorology, 1975.
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2

! parameter from Gerber 1985 (units require radius in cm, see rcm)
   real :: rcm
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424
! parameters for ammonium sulfate
   real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428

   real :: sat, rrat
   real :: radius, rhop   ! particle radius and density passed to fall velocity calculation
   real :: dtmax          ! max time step to satisfy CFL
   real :: dt_settle      ! sub time step of settling

   radius = radiusInp
   rhop = rhopInp
   if(radius .le. 0.) return    !get out of loop if particle radius <= 0

   do k = 1, km
     do j = j1, j2
       do i = i1, i2

! Adjust the particle size for relative humidity effects

       sat = max(rh(i,j,k),tiny(1.0)) !        to avoid zero FPE
       if(flag .eq. 1 .and. sat .ge. 0.80) then                                 ! Fitzgerald
! parameterization blows up for RH > 0.995, so set that as max rh needs to be scaled 0 - 1
        sat = min(0.995,sat)
! Calculate the alpha and beta parameters for the wet particle relative to amonium sulfate
        beta = exp( (0.00077*sat) / (1.009-sat) )
        if(sat .le. 0.97) then
         theta = 1.058
        else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
        endif
        alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
        f1 = 10.2 - 23.7*sat + 14.5*sat**2.
        f2 = -6.7 + 15.5*sat - 9.2*sat**2.
        alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
        alpha = alphaNaCl * (alpha1*alpharat)
        radius = alpha * radiusInp**beta      ! radius of the wet particle
        rrat = (radiusInp/radius)**3.
        rhop = rrat*rhopInp + (1.-rrat)*rhow

       else if(flag .eq. 2) then         ! Gerber
        sat = min(0.995,sat)
        rcm = radiusInp*100.
        radius = 0.01 * (c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                          + rcm**3.)**(1./3.)
        rrat = (radiusInp/radius)**3.
        rhop = rrat*rhopInp + (1.-rrat)*rhow

       else if(flag .eq. 3) then        ! Gerber parameterization for Ammonium Sulfate
        sat = min(0.995,sat)
        rcm = radiusInp*100.
        radius = 0.01 * (SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                      + rcm**3.)**(1./3.)
        rrat = (radiusInp/radius)**3.
        rhop = rrat*rhopInp + (1.-rrat)*rhow

       else if(flag .eq. 4) then        ! Petters and Kreidenweis (ACP2007) parameterization
        sat = min(0.99,sat)
        radius = (radiusInp**3 * (1+1.19*sat/(1-sat)))**(1./3.)
        rrat = (radiusInp/radius)**3
        rhop = rrat*rhopInp + (1.-rrat)*rhow
       endif

! Calculate the settling velocity
       call Chem_CalcVsettle(radius, rhop, rhoa(i,j,k),   &
                 tmpu(i,j,k), diff_coef, vsettle(i,j,k))
       end do
     end do
   end do 

! Determine the maximum time-step satisying CFL: dt <= min(dz/vsettle)
   dtmax = MINVAL((dz(:,:,:)/vsettle(:,:,:)))
   if(dtmax .ge. cdt) then
     nsubstep = 1
     dt_settle = cdt
   else
     nsubstep = cdt/dtmax + 1
     dt_settle = cdt/nsubstep
   end if

! Loop over sub-timeste
   do iit = 1, nsubstep
! Update mixing ratio with order of top->sfc
     qdel(i1:i2,j1:j2) = tc(i1:i2,j1:j2,km) * dt_settle * vsettle(i1:i2,j1:j2,km)/dz(i1:i2,j1:j2,km)  
     tc(i1:i2,j1:j2,km) = tc(i1:i2,j1:j2,km) - qdel(i1:i2,j1:j2)           !ug/kg

     do k = km-1, 1, -1
       qsrc(i1:i2,j1:j2) = qdel(i1:i2,j1:j2) * (dz(i1:i2,j1:j2,k+1)*rhoa(i1:i2,j1:j2,k+1))    &
                                             /(dz(i1:i2,j1:j2,k)*rhoa(i1:i2,j1:j2,k))
       qdel(i1:i2,j1:j2) = tc(i1:i2,j1:j2,k) * dt_settle * vsettle(i1:i2,j1:j2,k)/dz(i1:i2,j1:j2,k)
       tc(i1:i2,j1:j2,k) = tc(i1:i2,j1:j2,k) - qdel(i1:i2,j1:j2) + qsrc(i1:i2,j1:j2)
     end do
   end do  ! iit=1,nsubstep

end subroutine settling_nitrate


subroutine Chem_CalcVsettle(radius, rhop, rhoa, tmpu, diff_coef, vsettle)
!
! DESCRIPTION: Calculates the aerosol settling velocity and Brownian diffusion
!              coefficient
!              Follows discussions in Seinfeld and Pandis, Pruppacher and
!              Klett, and the coding in CARMA (Toon et al., 1988)
!              Should work satisfactorily for al reasonable sized aerosols
!              (up to Reynolds number 300)
   implicit none

   real, intent(in)    :: radius              ! Particle radius [m]
   real, intent(in)    :: rhop                ! Particle density [kg m-3]
   real, intent(out)   :: diff_coef           ! Brownian diffusion coefficient [m2 s-1]
   real, intent(out)   :: vsettle             ! Layer fall speed [m s-1]
   real*8, intent(in)    :: rhoa              ! Layer air density [kg m-3]
   real*8, intent(in)    :: tmpu              ! Layer temperature [K]

! Local Variables
   real*8 rmu                       ! Dynamic viscosity [kg m-1 s-1]
   real*8 vt                        ! Thermal velocity of air molecule [m s-1]
   real*8 rmfp                      ! Air molecule mean free path [m]
   real*8 bpm                       ! Cunningham slip correction factor
   real*8 rkn                       ! Knudsen number
   real*8 re, x, y                  ! reynold's number and parameters
   real, parameter :: grav = 9.81   ! Earth gravity constant [m s-2], added by ZTao
   real, parameter :: kb = 1.3807e-23    ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]
   real, parameter :: pi = 3.141529265

!  Dynamic viscosity from corrected Sutherland's Equation
   rmu = 1.8325e-5*(416.16/(tmpu+120.))*(tmpu/296.16)**1.5

!  Thermal velocity of air molecule
   vt = sqrt(8.*kb*tmpu/pi/m_air)

!  Air molecule mean free path
   rmfp = 2.*rmu/rhoa/vt

!  Knudsen number
   rkn = rmfp/radius

!  Cunningham slip correction factor
   bpm = 1. + 1.246*rkn + 0.42*rkn*exp(-0.87/rkn)

!  Brownian diffusion coefficient
   diff_coef = kb*tmpu*bpm/3./pi/rmu/(2.*radius)

!  Fall speed (assumes Reynolds # < 0.01)
   vsettle = 2./9.*rhop*radius**2.*grav*bpm/rmu

!  Check the Reynold's number to see if we need a drag correction
!  First guess at Reynold's number using Stoke's calculation
   re = 2.*rhoa*radius*vsettle/rmu

!  If Re > 0.01 then apply drag correction following Pruppacher and
!  Klett regime 2 (eq. 10-142).  Assuming reasonable aerosols we
!  do not consider that particle Re may exceed 300.
   if(re .gt. 0.01) then
    x = log(24.*re/bpm)
    y = -3.18657 + 0.992696   *x     - .00153193   *x**2. &
                 - 0.000987059*x**3. - .000578878  *x**4. &
                 + 8.55176E-05*x**5. -  3.27815E-06*x**6.
    re = exp(y)*bpm
    vsettle = rmu*re/2./rhoa/radius
   endif
end subroutine Chem_CalcVsettle


END MODULE MODULE_GOCART_SETTLING
