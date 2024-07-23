MODULE module_dust_load

! This module for calculation of dust loading
 CONTAINS

   SUBROUTINE dust_load_driver ( config_flags,                         &
            alt, chem, dz8w, dustload_1, dustload_2, dustload_3,       &
            dustload_4, dustload_5,                                    &
            dustscon_1, dustscon_2, dustscon_3,                        &  !dkim6
            dustscon_4, dustscon_5,                                    &  !dkim6
            dod,                                                       &  !dkim6
            ids,ide, jds,jde, kds,kde,                                 &
            ims,ime, jms,jme, kms,kme,                                 &
            its,ite, jts,jte, kts,kte                                  )

   USE module_configure
   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::                                   &
                                      ids,ide, jds,jde, kds,kde,       &
                                      ims,ime, jms,jme, kms,kme,       &
                                      its,ite, jts,jte, kts,kte

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),             &
         INTENT(IN ) :: chem

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                       &
         INTENT(IN ) :: alt, dz8w

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: dustload_1,   &
                dustload_2, dustload_3, dustload_4, dustload_5

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: dustscon_1,   &  !dkim6
                dustscon_2, dustscon_3, dustscon_4, dustscon_5            !dkim6

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: dod              !dkim6
   !dkim6, m2/g, 550nm, 5 bins, table from GEOS
   REAL, DIMENSION (5), PARAMETER :: mee_du(5)=(/2.016,0.644,0.334,0.167,0.084/)

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags

   INTEGER :: i, j, k

!      do j=jts,jte
!       do i=its,ite
       dustload_1(its:ite,jts:jte) = 0.
       dustload_2(its:ite,jts:jte) = 0.
       dustload_3(its:ite,jts:jte) = 0.
       dustload_4(its:ite,jts:jte) = 0.
       dustload_5(its:ite,jts:jte) = 0.
!       enddo
!      enddo

       ! dkim6
       dustscon_1(its:ite,jts:jte) = 0.
       dustscon_2(its:ite,jts:jte) = 0.
       dustscon_3(its:ite,jts:jte) = 0.
       dustscon_4(its:ite,jts:jte) = 0.
       dustscon_5(its:ite,jts:jte) = 0.
       dod(its:ite,jts:jte) = 0.

     do j=jts,jte
      do i=its,ite
       do k=kts,kte
! chem(p_dust) : [ug/kg_dryair], alt : [m3/kg], dz8w : [m] -> dustload : [ug/m2]
        dustload_1(i,j)= dustload_1(i,j) + chem(i,k,j,p_dust_1)/alt(i,k,j) * dz8w(i,k,j)
        dustload_2(i,j)= dustload_2(i,j) + chem(i,k,j,p_dust_2)/alt(i,k,j) * dz8w(i,k,j)
        dustload_3(i,j)= dustload_3(i,j) + chem(i,k,j,p_dust_3)/alt(i,k,j) * dz8w(i,k,j)
        dustload_4(i,j)= dustload_4(i,j) + chem(i,k,j,p_dust_4)/alt(i,k,j) * dz8w(i,k,j)
        dustload_5(i,j)= dustload_5(i,j) + chem(i,k,j,p_dust_5)/alt(i,k,j) * dz8w(i,k,j)
!         if (j.eq.int(0.5*(ite-its)).and.i.eq.int(0.5*(jte-jts))) write(6,*) 'dload', chem(i,k,j,p_dust_5)
        !dkim6, calculate dod = (ug/m2)*(m2/g)*10^-6(g/ug)
        dod(i,j)=dustload_1(i,j)*mee_du(1)+dustload_2(i,j)*mee_du(2) &
                +dustload_3(i,j)*mee_du(3)+dustload_4(i,j)*mee_du(4) &
                +dustload_5(i,j)*mee_du(5)
        dod(i,j)=dod(i,j)*(1.e-6)
       enddo
      enddo
     enddo

     !dkim6, dust scon
     do j=jts,jte
      do i=its,ite
       k=kts        !surface 
       !do k=kts,kte
       ! chem(p_dust) : [ug/kg_dryair], alt : [m3/kg] -> dustscon : [ug/m3]
        dustscon_1(i,j)= chem(i,k,j,p_dust_1)/alt(i,k,j)
        dustscon_2(i,j)= chem(i,k,j,p_dust_2)/alt(i,k,j)
        dustscon_3(i,j)= chem(i,k,j,p_dust_3)/alt(i,k,j)
        dustscon_4(i,j)= chem(i,k,j,p_dust_4)/alt(i,k,j)
        dustscon_5(i,j)= chem(i,k,j,p_dust_5)/alt(i,k,j)
       !enddo
      enddo
     enddo

   END SUBROUTINE dust_load_driver

END MODULE module_dust_load

