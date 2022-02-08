MODULE p4zmicro
   !!======================================================================
   !!                         ***  MODULE p4zmicro  ***
   !! TOP :   PISCES Compute the sources/sinks for microzooplankton
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_micro      : Compute the sources/sinks for microzooplankton
   !!   p4z_micro_init : Initialize and read the appropriate namelist
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zlim          ! Co-limitations
   USE p4zprod         ! production
   USE iom             ! I/O manager
   USE prtctl_trc      ! print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_micro         ! called in p4zbio.F90
   PUBLIC   p4z_micro_init    ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::   part        !: part of calcite not dissolved in microzoo guts
   REAL(wp), PUBLIC ::   xprefc      !: microzoo preference for POC 
   REAL(wp), PUBLIC ::   xprefn      !: microzoo preference for nanophyto
   REAL(wp), PUBLIC ::   xprefd      !: microzoo preference for diatoms
   REAL(wp), PUBLIC ::   xthreshdia  !: diatoms feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshphy  !: nanophyto threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthreshpoc  !: poc threshold for microzooplankton 
   REAL(wp), PUBLIC ::   xthresh     !: feeding threshold for microzooplankton 
   REAL(wp), PUBLIC ::   resrat      !: exsudation rate of microzooplankton
   REAL(wp), PUBLIC ::   mzrat       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::   grazrat     !: maximal microzoo grazing rate
   REAL(wp), PUBLIC ::   xkgraz      !: Half-saturation constant of assimilation
   REAL(wp), PUBLIC ::   unass       !: Non-assimilated part of food
   REAL(wp), PUBLIC ::   sigma1      !: Fraction of microzoo excretion as DOM 
   REAL(wp), PUBLIC ::   epsher      !: growth efficiency for grazing 1 
   REAL(wp), PUBLIC ::   epshermin   !: minimum growth efficiency for grazing 1
   REAL(wp), PUBLIC ::   e13c_calz   !: C13 microzoo calification fractionation
   REAL(wp), PUBLIC ::   e15n_ex     !: N15 microzoo excretion fractionation
   REAL(wp), PUBLIC ::   e15n_in     !: N15 microzoo ingestion fractionation
   REAL(wp), PUBLIC ::   e18oxy_zoo  !: O18 microzoo respiration fractionation

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmicro.F90 10374 2018-12-06 09:49:35Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_micro( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_micro  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for microzooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time step
      INTEGER, INTENT(in) ::   knt   ! ??? 
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaz , zcompaph, zcompapoc
      REAL(wp) :: zgraze  , zdenom, zdenom2
      REAL(wp) :: zfact   , zfood, zfoodlim, zbeta
      REAL(wp) :: zepsherf, zepshert, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem, zgrafer, zgrapoc, zprcaca, zmortz
      REAL(wp) :: zrespz, ztortz, zgrasrat, zgrasratn
      REAL(wp) :: zgrazp, zgrazm, zgrazsd
      REAL(wp) :: zgrazmf, zgrazsf, zgrazpf
      REAL(wp) :: zgrazp13, zgrazm13, zgrazsd13, zgraztotc13
      REAL(wp) :: zgrarem_13, zgrapoc_13, zgrasig_13, zmortz_13, zr13_dic
      REAL(wp) :: zgrazp15, zgrazm15, zgrazsd15, zgraztotc15
      REAL(wp) :: zgrarem_15, zgrapoc_15, zgrasig_15, zgrasigex_15, zmortz_15
      REAL(wp) :: zr18_oxy
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrazing1, zfezoo
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: excretion1, excretion1_13, excretion1_15
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zw3d, zzligprod
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_micro')
      !
      IF (ln_ligand) THEN
         ALLOCATE( zzligprod(jpi,jpj,jpk) )
         zzligprod(:,:,:) = 0._wp
      ENDIF
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompaz = MAX( ( trb(ji,jj,jk,jpzoo) - 1.e-9 ), 0.e0 )
               zfact   = xstep * tgfunc2(ji,jj,jk) * zcompaz

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz = resrat * zfact * trb(ji,jj,jk,jpzoo) / ( xkmort + trb(ji,jj,jk,jpzoo) )  &
                  &   + resrat * zfact * 3. * nitrfac(ji,jj,jk)

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation.
               !  ---------------------------------------------------------------
               ztortz = mzrat * 1.e6 * zfact * trb(ji,jj,jk,jpzoo) * (1. - nitrfac(ji,jj,jk))

               zcompadi  = MIN( MAX( ( trb(ji,jj,jk,jpdia) - xthreshdia ), 0.e0 ), xsizedia )
               zcompaph  = MAX( ( trb(ji,jj,jk,jpphy) - xthreshphy ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,jk,jppoc) - xthreshpoc ), 0.e0 )
               
               !     Microzooplankton grazing
               !     ------------------------
               zfood     = xprefn * zcompaph + xprefc * zcompapoc + xprefd * zcompadi
               zfoodlim  = MAX( 0. , zfood - min(xthresh,0.5*zfood) )
               zdenom    = zfoodlim / ( xkgraz + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )
               zgraze    = grazrat * xstep * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpzoo) * (1. - nitrfac(ji,jj,jk))

               zgrazp    = zgraze  * xprefn * zcompaph  * zdenom2 
               zgrazm    = zgraze  * xprefc * zcompapoc * zdenom2 
               zgrazsd   = zgraze  * xprefd * zcompadi  * zdenom2 

               IF( ln_c13 ) THEN
                  zgrazp13    = zgrazp  * ( (trb(ji,jj,jk,jp13phy)+rtrn) / (trb(ji,jj,jk,jpphy)+rtrn) )
                  zgrazm13    = zgrazm  * ( (trb(ji,jj,jk,jp13poc)+rtrn) / (trb(ji,jj,jk,jppoc)+rtrn) )
                  zgrazsd13   = zgrazsd * ( (trb(ji,jj,jk,jp13dia)+rtrn) / (trb(ji,jj,jk,jpdia)+rtrn) )
                  zgraztotc13 = zgrazp13 + zgrazm13 + zgrazsd13
               ENDIF
               IF( ln_n15 ) THEN
                  zgrazp15    = zgrazp  * ( (trb(ji,jj,jk,jp15phy)+rtrn) / (trb(ji,jj,jk,jpphy)+rtrn) )
                  zgrazm15    = zgrazm  * ( (trb(ji,jj,jk,jp15poc)+rtrn) / (trb(ji,jj,jk,jppoc)+rtrn) )
                  zgrazsd15   = zgrazsd * ( (trb(ji,jj,jk,jp15dia)+rtrn) / (trb(ji,jj,jk,jpdia)+rtrn) )
                  zgraztotc15 = zgrazp15 + zgrazm15 + zgrazsd15
               ENDIF

               zgrazpf   = zgrazp  * trb(ji,jj,jk,jpnfe) / (trb(ji,jj,jk,jpphy) + rtrn)
               zgrazmf   = zgrazm  * trb(ji,jj,jk,jpsfe) / (trb(ji,jj,jk,jppoc) + rtrn)
               zgrazsf   = zgrazsd * trb(ji,jj,jk,jpdfe) / (trb(ji,jj,jk,jpdia) + rtrn)
               !
               zgraztotc = zgrazp  + zgrazm  + zgrazsd 
               zgraztotf = zgrazpf + zgrazsf + zgrazmf 
               zgraztotn = zgrazp * quotan(ji,jj,jk) + zgrazm + zgrazsd * quotad(ji,jj,jk)

               ! Grazing by microzooplankton
               zgrazing1(ji,jj,jk) = zgraztotc

               !    Various remineralization and excretion terms
               !    --------------------------------------------
               zgrasrat  = ( zgraztotf + rtrn ) / ( zgraztotc + rtrn )
               zgrasratn = ( zgraztotn + rtrn ) / ( zgraztotc + rtrn )
               zepshert  =  MIN( 1., zgrasratn, zgrasrat / ferat3)
               zbeta     = MAX(0., (epsher - epshermin) )
               zepsherf  = epshermin + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta )
               zepsherv  = zepsherf * zepshert 

               zgrafer   = zgraztotc * MAX( 0. , ( 1. - unass ) * zgrasrat - ferat3 * zepsherv ) 
               zgrarem   = zgraztotc * ( 1. - zepsherv - unass )
               zgrapoc   = zgraztotc * unass

               IF ( ln_c13 ) THEN
                  zgrarem_13 = zgraztotc13 * ( 1. - zepsherv - unass )
                  zgrapoc_13 = zgraztotc13 * unass
                  zgrasig_13 = zgrarem_13 * sigma1
                  zmortz_13  = (ztortz + zrespz) * ( (trb(ji,jj,jk,jp13zoo)+rtrn) / (trb(ji,jj,jk,jpzoo)+rtrn) )
               ENDIF
               IF ( ln_n15 ) THEN
                  zgrarem_15 = zgraztotc15 * ( 1. - zepsherv - unass )
                  zgrapoc_15 = zgraztotc15 * unass
                  zgrasig_15 = zgrarem_15 * sigma1
                  zgrasigex_15 = ( 1. - epsher - unass ) * zgraztotc15 * sigma1 * zgrasratn
                    ! zgrasigex_15 = amount of NH4 excreted by zooplankton
                    ! according to measure of food quality [0,1] (zgrasratn) multiplied by
                    ! the minimum possible excretion (( 1. - epsher - unass ) * zgraztotc15 * sigma1)
                  zmortz_15  = (ztortz + zrespz) * ( (trb(ji,jj,jk,jp15zoo)+rtrn) / (trb(ji,jj,jk,jpzoo)+rtrn) )
               ENDIF

               !  Update of the TRA arrays
               !  ------------------------
               zgrarsig  = zgrarem * sigma1
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               excretion1(ji,jj,jk) = zgrarsig
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem - zgrarsig
               !
               IF( ln_ligand ) THEN
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + (zgrarem - zgrarsig) * ldocz
                  zzligprod(ji,jj,jk) = (zgrarem - zgrarsig) * ldocz
               ENDIF
               !
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer
               zfezoo(ji,jj,jk)    = zgrafer
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zgrapoc
               prodpoc(ji,jj,jk)   = prodpoc(ji,jj,jk) + zgrapoc
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zgraztotf * unass
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig

               IF( ln_c13 ) THEN
                  tra(ji,jj,jk,jp13doc) = tra(ji,jj,jk,jp13doc) + zgrarem_13 - zgrasig_13
                  tra(ji,jj,jk,jp13poc) = tra(ji,jj,jk,jp13poc) + zgrapoc_13
                  tra(ji,jj,jk,jp13dic) = tra(ji,jj,jk,jp13dic) + zgrasig_13
                  excretion1_13(ji,jj,jk) = zgrasig_13
               ENDIF
               IF( ln_n15 ) THEN
                  excretion1_15(ji,jj,jk) = zgrasigex_15 * ( 1. - e15n_ex/1000.0 ) + ( zgrasig_15 - zgrasigex_15 )
                  tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4) + excretion1_15(ji,jj,jk)
                  tra(ji,jj,jk,jp15doc) = tra(ji,jj,jk,jp15doc) + zgrarem_15 - zgrasig_15
                  tra(ji,jj,jk,jp15poc) = tra(ji,jj,jk,jp15poc) + zgrapoc_15 * ( 1. - e15n_in/1000.0 )
               ENDIF
               IF( ln_o18 ) THEN
                  zr18_oxy = ( trb(ji,jj,jk,jp18oxy) + rtrn ) / ( trb(ji,jj,jk,jpoxy) + rtrn )
                  tra(ji,jj,jk,jp18oxy) = tra(ji,jj,jk,jp18oxy) - o2ut * zgrarsig * zr18_oxy * (1. - e18oxy_zoo/1000.)
               ENDIF
               !   Update the arrays TRA which contain the biological sources and sinks
               !   --------------------------------------------------------------------
               zmortz = ztortz + zrespz
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zmortz + zepsherv * zgraztotc 
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazp
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazsd
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazp  * trb(ji,jj,jk,jpnch)/(trb(ji,jj,jk,jpphy)+rtrn)
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazsd * trb(ji,jj,jk,jpdch)/(trb(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazsd * trb(ji,jj,jk,jpdsi)/(trb(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazsd * trb(ji,jj,jk,jpdsi)/(trb(ji,jj,jk,jpdia)+rtrn)
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgrazpf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazsf
               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zmortz - zgrazm
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zmortz
               conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazm
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + ferat3 * zmortz - zgrazmf
               !
               ! calcite production
               zprcaca = xfracal(ji,jj,jk) * zgrazp
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)
               !
               zprcaca = part * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * zprcaca
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) + zprcaca
               IF ( ln_c13 ) THEN
                  zr13_dic = ( (trb(ji,jj,jk,jp13dic)+rtrn) / (trb(ji,jj,jk,jpdic)+rtrn) )
                  tra(ji,jj,jk,jp13zoo) = tra(ji,jj,jk,jp13zoo) - zmortz_13 + zepsherv * zgraztotc13
                  tra(ji,jj,jk,jp13phy) = tra(ji,jj,jk,jp13phy) - zgrazp13
                  tra(ji,jj,jk,jp13dia) = tra(ji,jj,jk,jp13dia) - zgrazsd13
                  tra(ji,jj,jk,jp13poc) = tra(ji,jj,jk,jp13poc) + zmortz_13 - zgrazm13
                  tra(ji,jj,jk,jp13dic) = tra(ji,jj,jk,jp13dic) - zprcaca * zr13_dic * (1. - e13c_calz/1000.)
                  tra(ji,jj,jk,jp13cal) = tra(ji,jj,jk,jp13cal) + zprcaca * zr13_dic * (1. - e13c_calz/1000.)
                  prodcal13(ji,jj,jk) = prodcal13(ji,jj,jk) + zprcaca * zr13_dic * (1. - e13c_calz/1000.)
               ENDIF
               IF ( ln_n15 ) THEN
                  tra(ji,jj,jk,jp15zoo) = tra(ji,jj,jk,jp15zoo) - zmortz_15 + zepsherv * zgraztotc15  &
                  &                       + zgrasigex_15 * (e15n_ex/1000.0) + zgrapoc_15 * (e15n_in/1000.0)
                  tra(ji,jj,jk,jp15phy) = tra(ji,jj,jk,jp15phy) - zgrazp15
                  tra(ji,jj,jk,jp15dia) = tra(ji,jj,jk,jp15dia) - zgrazsd15
                  tra(ji,jj,jk,jp15poc) = tra(ji,jj,jk,jp15poc) + zmortz_15 - zgrazm15
               ENDIF
            END DO
         END DO
      END DO
      !
      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
           ALLOCATE( zw3d(jpi,jpj,jpk) )
           IF( iom_use( "GRAZ1" ) ) THEN
              zw3d(:,:,:) = zgrazing1(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !  Total grazing of phyto by zooplankton
              CALL iom_put( "GRAZ1", zw3d )
           ENDIF
           IF( iom_use( "EXCR1" ) ) THEN
              zw3d(:,:,:) = excretion1(:,:,:) * rno3 * 1.e+3 * rfact2r * tmask(:,:,:) !  Total excretion of NH4 by zooplankton
              CALL iom_put( "EXCR1", zw3d )
           ENDIF
           IF( iom_use( "FEZOO" ) ) THEN
              zw3d(:,:,:) = zfezoo(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
              CALL iom_put( "FEZOO", zw3d )
           ENDIF
           IF( iom_use( "LPRODZ" ) .AND. ln_ligand )  THEN
              zw3d(:,:,:) = zzligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)
              CALL iom_put( "LPRODZ"  , zw3d )
           ENDIF
           DEALLOCATE( zw3d )
         ENDIF
      ENDIF
      !
      IF (ln_ligand)  DEALLOCATE( zzligprod )
      !
      IF(ln_ctl) THEN      ! print mean trends (used for debugging)
         WRITE(charout, FMT="('micro')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_micro')
      !
   END SUBROUTINE p4z_micro


   SUBROUTINE p4z_micro_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_micro_init  ***
      !!
      !! ** Purpose :   Initialization of microzooplankton parameters
      !!
      !! ** Method  :   Read the nampiszoo namelist and check the parameters
      !!                called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampiszoo
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zzoo/ part, grazrat, resrat, mzrat, xprefn, xprefc, &
         &                xprefd,  xthreshdia,  xthreshphy,  xthreshpoc, &
         &                xthresh, xkgraz, epsher, epshermin, sigma1, unass, &
         &                e13c_calz, e15n_ex, e15n_in, e18oxy_zoo
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p4z_micro_init : Initialization of microzooplankton parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampiszoo in reference namelist : Pisces microzooplankton
      READ  ( numnatp_ref, namp4zzoo, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zzoo in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampiszoo in configuration namelist : Pisces microzooplankton
      READ  ( numnatp_cfg, namp4zzoo, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zzoo in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zzoo )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zzoo'
         WRITE(numout,*) '      part of calcite not dissolved in microzoo guts  part        =', part
         WRITE(numout,*) '      microzoo preference for POC                     xprefc      =', xprefc
         WRITE(numout,*) '      microzoo preference for nano                    xprefn      =', xprefn
         WRITE(numout,*) '      microzoo preference for diatoms                 xprefd      =', xprefd
         WRITE(numout,*) '      diatoms feeding threshold  for microzoo         xthreshdia  =', xthreshdia
         WRITE(numout,*) '      nanophyto feeding threshold for microzoo        xthreshphy  =', xthreshphy
         WRITE(numout,*) '      poc feeding threshold for microzoo              xthreshpoc  =', xthreshpoc
         WRITE(numout,*) '      feeding threshold for microzooplankton          xthresh     =', xthresh
         WRITE(numout,*) '      exsudation rate of microzooplankton             resrat      =', resrat
         WRITE(numout,*) '      microzooplankton mortality rate                 mzrat       =', mzrat
         WRITE(numout,*) '      maximal microzoo grazing rate                   grazrat     =', grazrat
         WRITE(numout,*) '      non assimilated fraction of P by microzoo       unass       =', unass
         WRITE(numout,*) '      Efficicency of microzoo growth                  epsher      =', epsher
         WRITE(numout,*) '      Minimum efficicency of microzoo growth          epshermin   =', epshermin
         WRITE(numout,*) '      Fraction of microzoo excretion as DOM           sigma1      =', sigma1
         WRITE(numout,*) '      half sturation constant for grazing 1           xkgraz      =', xkgraz
         WRITE(numout,*) '      C13 microzoo calcification fractionation        e13c_calz   =', e13c_calz
         WRITE(numout,*) '      N15 microzoo excretion fractionation            e15n_ex     =', e15n_ex
         WRITE(numout,*) '      N15 microzoo ingestion fractionation            e15n_in     =', e15n_in
         WRITE(numout,*) '      O18 microzoo respiration fractionation          e18oxy_zoo  =', e18oxy_zoo
      ENDIF
      !
   END SUBROUTINE p4z_micro_init

   !!======================================================================
END MODULE p4zmicro
