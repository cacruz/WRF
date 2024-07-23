MODULE MODULE_GOCART_WETSCAV

!-------------------------------------------------------------------------
!  NUWRF:
!       Calculate aerosol and gas (NH3) loss due to large-scale (LS)
!	   wet scavange in GOCART related aerosol schemes
!       Code is adapted from NASA/GSFC GMAO by Z. Tao
!    
!   Revision history
!    1. Feb. 2022  Z. Tao
!-------------------------------------------------------------------------

  implicit none
  real, parameter :: grav = 9.81

!-------------------------------------------------------------------------
CONTAINS
!
!------------------------------------------------------------------------------
! subroutines wetscav_gocart_driver: calculate aerosol wet removal due 
!                                    to large scale processes.
!------------------------------------------------------------------------------

      subroutine wetscav_gocart_driver(cdt, aero_type, kin, fwet, qa, ple, tmpu, rhoa,   &
                                          precr, preci, precs, precg, raincv, rainncv,   &
                                          ids,ide, jds,jde, kds,kde,                     &
                                          ims,ime, jms,jme, kms,kme,                     &
                                          its,ite, jts,jte, kts,kte)

!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
!-- cdt		  chemistry time step (s)
!-- aero_type     aerosol names (character)
!-- kin		  flag (.true. for aersols)
!-- fwet          Wet removal efficiency
!-- qa		  species mixing ratio
!-- ple		  pressure at full levels (Pa)
!-- tmpu	  temperature (K)
!-- rhoa	  air density (kg/m^3)
!-- precr         rain precipitation rate at all levels (kg/m2/s)
!-- preci         ice precipitation rate at all levels (kg/m2/s)
!-- precs         snow precipitation rate at all levels (kg/m2/s)
!-- precg         graupel precipitation rate at all levels (kg/m2/s)                             &
!-- raincv        time-step cumulus precipitation (mm)
!-- rainncv       time-step nonconvectiv precipitation (mm)
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!======================================================================

        implicit none

        integer, intent(in) :: ids,ide, jds,jde, kds,kde,  &
                               ims,ime, jms,jme, kms,kme,  &
                               its,ite, jts,jte, kts,kte  
	real, intent(in) :: cdt, fwet
	character(len=*) :: aero_type
	logical, intent(in) :: kin
	real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: qa
	real, dimension(ims:ime, kms:kme, jms:jme), intent(in) ::   &
                 ple, tmpu, rhoa, precr, preci, precs, precg
	real, dimension(ims:ime , jms:jme), intent(in) :: raincv, rainncv

! !Local Variables
	integer  ::  i, j, k, iit, LH, kk
	real :: pac                       ! tot precip = large-scale (ls) + convective (cv) [mm/time-step]
	real :: qls(kms:kme)              ! ls portion of moisture tendency [kg m-3 s-1]
	real :: fd(kms:kme)               ! flux across layers (kg/m^2)
	real :: qmx, qd, A                ! temporary variables on moisture
	real :: F, B, BT                  ! temporary variables on cloud, freq.
	real :: DC                        ! scavenge change in mass mixing ratio
	real :: WASHFRAC, WASHFRAC_F_14
        real, dimension(ims:ime, kms:kme, jms:jme) :: c_h2o, cldliq, cldice 
        real, dimension(ims:ime, kms:kme, jms:jme) :: pdog   !air mass factor = dp/g (kg/m^2)
        real, dimension(ims:ime, kms:kme, jms:jme) :: delz   !layer thickness = dp/g/rhoa (m)
        real, dimension(ims:ime, kms:kme, jms:jme) :: dpfli  !vertical gradient of LS ice+rain precip flux
        character (len=1000) :: msg

!  Rain parameters from Liu et al.
	real, parameter :: B0_ls = 1.0e-4
	real, parameter :: F0_ls = 1.0
	real, parameter :: XL_ls = 5.0e-4
	real, parameter :: B0_cv = 1.5e-3
	real, parameter :: F0_cv = 0.3
	real, parameter :: XL_cv = 2.0e-3
	real, parameter :: k_wash = 1.d0  ! first order washout rate, constant, [cm^-1]
!  Duration of rain: ls = model timestep, cv = 1800 s (<= cdt)
	real            :: Td_ls
	real, parameter :: Td_cv = 1800.
	REAL*8,  PARAMETER   :: R = 8.2057d-2  ! universal gas constant [L*atm/moles/K]
	REAL*8,  PARAMETER   :: INV_T0 = 1d0 / 298d0
	REAL*8,  PARAMETER   :: conv_NH3 = 5.69209978831d-1 ! 0.6*SQRT(0.9) for ice to gas ratio
	REAL*8  :: k_rain, Kstar298, H298_R, I2G, L2G, C_TOT, F_L, F_I
	REAL*8  :: PP, LP
	logical :: snow_scavenging
	logical :: rain_scavenging
!  Efficiency of dust wet removal (since dust is really not too hygroscopic)
!  Applied only to in-cloud scavenging
	real :: effRemoval

! Initialize local arrays to make sure no bad number exists
        c_h2o = 0.
        pdog = 0.
        dpfli = 0.

!  Initialize local variables
!  --------------------------
!  c_h2o, cldliq, and cldice are respectively the 
!  water mixing ratio, cloud liquid water mixing ratio,
!  and cloud ice water mixing ratio
        c_h2o  = (10**(-2663.5/tmpu(:,:,:) + 12.537)) /       &
                 (ple(:,kms:kme-1,:)+ple(:,kms+1:kme,:))/2.
        cldliq = 0.
        where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.) / 20. )
        where(tmpu >= 268.) cldliq = 1.d-6
        cldice = 1.d-6 - cldliq
	Td_ls = cdt

!  Accumulate the 3-dimensional arrays of rhoa and pdog
	pdog = (ple(:,kms:kme-1,:)-ple(:,kms+1:kme,:)) / grav
	delz = pdog / rhoa
	dpfli = precr(:,kms:kme-1,:) - precr(:,kms+1:kme,:) +   &
	        preci(:,kms:kme-1,:) - preci(:,kms+1:kme,:) +   &
	        precs(:,kms:kme-1,:) - precs(:,kms+1:kme,:) +   &
	        precg(:,kms:kme-1,:) - precg(:,kms+1:kme,:) 

	if (.not. kin) then              ! Gases
	  if (aero_type == 'NH3' .or. aero_type == 'nh3') then   ! Only for NH3 at present
! values adopted in Umich/IMPACT and GMI, effective Henry's law coefficient at pH=5
	    Kstar298 = 1.05d6
	    H298_R = -4.2d3
	  else
	    write(msg, '("GOCART gaseous wet scaveging only works for NH3")') 
	    CALL wrf_error_fatal( trim(msg) )
	  endif
	endif

! Rain scavenging flag
	rain_scavenging = .true.
	if ( aero_type == 'OC1' ) rain_scavenging = .false.

!  Snow scavenging flag
	snow_scavenging = .true.
	if (aero_type == 'BC2' .or. aero_type == 'OC2'      .or.     &
	    aero_type == 'SEAS_1' .or. aero_type == 'SEAS_2' .or.    &
	    aero_type == 'SEAS_3' .or. aero_type == 'SEAS_4') then
	  snow_scavenging = .false.
	endif

	do j = jts, jte
	  do i = its, ite
	    pac = raincv(i,j) + rainncv(i,j)
	    if(pac .le. 0.) goto 100              !no wet scavenge if there is no precip.

	    qls = 0.
	    fd = 0.

!    Find the highest model layer experiencing rainout. 
!	Assumes no scavenging if T < 258 K
	    LH = kte + 1
	    do k = kte, kts, -1
	      if(dpfli(i,k,j) .gt. 0. ) then
	        LH = k
	        exit
              endif
            end do
            if(LH .gt. kte) goto 100

!    change to moisture tendency (kg/m3/s)
	    do k = LH, kts, -1
              qls(k) = dpfli(i,k,j)/pdog(i,k,j)*rhoa(i,k,j)
            end do

! Loop over vertical from top down to do the scavenging!
            do k = LH, kts, -1

!-----------------------------------------------------------------------------
!   (1) LARGE-SCALE RAINOUT:             
!       Tracer loss by rainout = TC0 * F * [1.-exp(-B*dt)]
!         where B = precipitation frequency,
!               F = fraction of grid box covered by precipitating clouds.
!       We assume that tracer scavenged by rain is falling down to the
!       next level, where a fraction could be re-evaporated to gas phase
!       if Qls is less then 0 in that level.
!-----------------------------------------------------------------------------
              if (qls(k) .gt. 0.) then
                F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
                k_rain  = B0_ls/F0_ls +1./(F0_ls*XL_ls/qls(k))
                if ( kin ) then     ! Aerosols
                  B = k_rain
                else                ! Gases
   ! ice to gas ratio
                  if ( c_h2o(i,k,j) > 0.) then
                    I2G = (cldice(i,k,j) / c_h2o(i,k,j)) * conv_NH3
                  else
                    I2G = 0.d0
                  endif

                  L2G = cldliq(i,k,j) * R * tmpu(i,k,j) * Kstar298 *      &
		      EXP( -H298_R * ( ( 1d0 / tmpu(i,k,j) ) - INV_T0 ) )
   ! fraction of NH3 in liquid & ice phases
                  C_TOT = 1d0 + L2G + I2G
                  F_L = L2G / C_TOT
                  F_I = I2G / C_TOT
   ! compute kg, the retention factor for liquid NH3 is 0 at T < 248K 
   !    and 0.05 at 248K < T < 268K
                  if (tmpu(i,k,j) >=268d0) then
                    B = k_rain * ( F_L+F_I )
                  else if ( (248d0 < tmpu(i,k,j)) .and. (tmpu(i,k,j) < 268d0) ) then
                    B = k_rain * ( (0.05*F_L)+F_I )
                  else
                    B = k_rain * F_I
                  end if
                end if ! kin

                BT = B * Td_ls
                if (BT.gt.10.) BT = 10.               !< Avoid overflow >
   ! Adjust du level:
   ! supress scavenging at cold T except for HNO3
                if (tmpu(i,k,j) >= 258d0 .and. .not.rain_scavenging) F = 0.d0
                if (tmpu(i,k,j) < 258d0 .and. .not.snow_scavenging) F = 0.d0
                effRemoval = fwet
                if (tmpu(i,k,j) < 258d0 .and. snow_scavenging) effRemoval = 1.

                DC = qa(i,k,j) * F * effRemoval * (1.-exp(-BT))
                if (DC.lt.0.) DC = 0.
                qa(i,k,j) = qa(i,k,j) - DC
                if (qa(i,k,j) .lt. 1.0E-32) qa(i,k,j) = 1.0E-32

   ! Flux down:  unit is kg m-2
   ! Formulated in terms of production in the layer.  In the revaporation step
   ! we consider possibly adding flux from above...
                fd(k) = DC * pdog(i,k,j)
              end if     !qls > 0

!-----------------------------------------------------------------------------
! * (2) LARGE-SCALE WASHOUT:
! *     Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------
              if(k .lt. LH .and. qls(k) .ge. 0.) then
                if(qls(k) .lt. qls(k+1)) then
   ! Find a maximum F overhead until the level where Qls<0.
                  qmx = 0.
                  do kk = k+1, LH
                    if (qls(kk) .gt. 0.) then
                      qmx = max(qmx, qls(kk))
                    else
                      goto 333
                    end if
                  end do
333               continue
                  F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qmx*cdt/Td_ls))
                  if (F.lt.0.01) F = 0.01

!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa. Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------
!
   ! Aerosols
                  qd = qmx / rhoa(i,k,j) * pdog(i,k,j)     !back to flux (kg/m2/s)
                  if (qd .ge. 50.) then
                    B = 0.
                  else
                    B = qd * 0.1
                  end if
                  BT = B * cdt
                  if (BT .gt. 10.) BT = 10.
   ! Gases
                  if ( .not. kin ) then
                    IF ( tmpu(i,k,j) >= 268d0 ) THEN
   !------------------------
   ! T >= 268K: Do washout
   !------------------------
   ! Rainwater content in the grid box (Eq. 17, Jacob et al, 2000)
   ! densities(kg/m3): water = 1000; ice = 917; snow = 100; graupel = 200
                      PP = (precr(i,k,j)/1000d0 + preci(i,k,j)/917d0  +        &
                            precs(i,k,j)/100d0 + precg(i,k,j)/200d0) * 100d0 ! from kg H2O/m2/s to cm3 H2O/cm2/s
                      LP = ( PP * cdt ) / ( F * delz(i,k,j) * 100.d0 )      ! DZ*100.d0 in cm
   ! Compute liquid to gas ratio for H2O2, using the appropriate 
   ! parameters for Henry's law -- also use rainwater content Lp
   ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
   !CALL COMPUTE_L2G( Kstar298, H298_R, tmpu(i,k,j), LP, L2G )
                      L2G = Kstar298 * EXP( -H298_R*((1d0/tmpu(i,k,j))-INV_T0) ) &
                                * LP * R * tmpu(i,k,j)
   ! Washout fraction from Henry's law (Eq. 16, Jacob et al, 2000)
                      WASHFRAC = L2G / ( 1d0 + L2G )
   ! Washout fraction / F from Eq. 14, Jacob et al, 2000
   ! Note: WASHFRAC_F_14 should match what's used for HNO3 (hma, 13aug2011)
                      WASHFRAC_F_14 = 1d0 - EXP( -K_WASH * ( PP / F ) * cdt )
   ! Do not let the Henry's law washout fraction exceed
                      IF ( WASHFRAC > WASHFRAC_F_14 ) WASHFRAC = WASHFRAC_F_14
                    ELSE
   !------------------------
   ! T < 268K: No washout
   !------------------------
                      WASHFRAC = 0d0
                    ENDIF
                  endif    !not kin

   ! Adjust du level:
                  if ( kin ) then
                    DC = qa(i,k,j) * F * (1.-exp(-BT))
                  else
                    DC = qa(i,k,j) * F * WASHFRAC
                  end if

                  if (DC .lt. 0.) DC = 0.
                  qa(i,k,j) = qa(i,k,j) - DC
                  if (qa(i,k,j) .lt. 1.0E-32) qa(i,k,j) = 1.0E-32
	        end if      !qls(k) < qls(k+1)
	      end if

!-----------------------------------------------------------------------------
!  (3) RE-EVAPORATION.  Assume that SO2 is re-evaporated as SO4 since it
!      has been oxidized by H2O2 at the level above. 
!-----------------------------------------------------------------------------
   ! Add in the flux from above, which will be subtracted if reevaporation occurs
              if(k .lt. LH) then
                fd(k) = fd(k) + fd(k+1)
   ! Is there evaporation in the currect layer?
                if (dpfli(i,k,j) .lt. 0.) then
   ! Fraction evaporated = H2O(k)evap / H2O(next condensation level).
                  if (dpfli(i,k+1,j) .gt. 0.) then
                    A =  abs(dpfli(i,k,j) / dpfli(i,k+1,j))
                    if (A .gt. 1.) A = 1.
   ! Adjust tracer in the level
                    DC =  fd(k+1) / pdog(i,k,j) * A
                    qa(i,k,j) = qa(i,k,j) + DC
                    qa(i,k,j) = max(qa(i,k,j), 1.e-32)
   ! Adjust the flux out of the bottom of the layer
                    fd(k) = fd(k) - DC * pdog(i,k,j)
	          end if
	        end if
	      end if
	    end do   ! k = LH, kts, -1

100	    continue
	  end do
	end do

	return
      end subroutine wetscav_gocart_driver

END MODULE MODULE_GOCART_WETSCAV
