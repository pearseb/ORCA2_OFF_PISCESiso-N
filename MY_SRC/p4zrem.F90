MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zlim
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   REAL(wp), PUBLIC ::   xremikc    !: remineralisation rate of DOC 
   REAL(wp), PUBLIC ::   xremikn    !: remineralisation rate of DON 
   REAL(wp), PUBLIC ::   xremikp    !: remineralisation rate of DOP 
   REAL(wp), PUBLIC ::   xremik     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   nitrif     !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::   xsirem     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsiremlab  !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsilab     !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::   feratb     !: Fe/C quota in bacteria
   REAL(wp), PUBLIC ::   xkferb     !: Half-saturation constant for bacteria Fe/C
   REAL(wp), PUBLIC ::   mu_aoa     !: Growth rate of ammonia-oxidising archaea
   REAL(wp), PUBLIC ::   kaoanh4    !: NH4 half-saturation constant for ammonia-oxidising archaea 
   REAL(wp), PUBLIC ::   kaoafer    !: Fer half-saturation constant for ammonia-oxidising archaea 
   REAL(wp), PUBLIC ::   mu_nob     !: Growth rate of nitrite-oxidising bacteria
   REAL(wp), PUBLIC ::   knobno2    !: NO2 half-saturation constant for nitrite-oxidising bacteria
   REAL(wp), PUBLIC ::   knobfer    !: Fer half-saturation constant for nitrite-oxidising bacteria
   REAL(wp), PUBLIC ::   ranammox   !: Maximum rate of anammox (per day)
   REAL(wp), PUBLIC ::   e15n_nar   !: N15 fractionation due to nitrate reduction
   REAL(wp), PUBLIC ::   e15n_nir   !: N15 fractionation due to nitrite reduction
   REAL(wp), PUBLIC ::   e15n_amo   !: N15 fractionation due to ammonium oxidation
   REAL(wp), PUBLIC ::   e15n_nio   !: N15 fractionation due to nitrite oxidation
   REAL(wp), PUBLIC ::   e15n_amm   !: N15 fractionation due to ammonification [DOC-->NH4]
   REAL(wp), PUBLIC ::   e15n_xamo  !: N15 fractionation due to ammonium oxidiation (anammox)
   REAL(wp), PUBLIC ::   e15n_xnir  !: N15 fractionation due to nitrite reduction (anammox)
   REAL(wp), PUBLIC ::   e15n_xnio  !: N15 fractionation due to nitrite oxidation (anammox)
   REAL(wp), PUBLIC ::   e18o_nar   !: O18 fractionation due to nitrate reduction
   REAL(wp), PUBLIC ::   e18o_nir   !: O18 fractionation due to nitrite reduction
   REAL(wp), PUBLIC ::   e18o_o2_h2o!: O18 fractionation due to oxygen incorportation from O2 and H2O
   REAL(wp), PUBLIC ::   e18o_h2o_2 !: O18 fractionation due to oxygen incorportation from H2O (nitrite oxidation)
   REAL(wp), PUBLIC ::   e18o_no2   !: O18 fractionation due to oxygen selection of N2O (nitrite oxidation)
   REAL(wp), PUBLIC ::   e18o_xnir  !: O18 fractionation due to nitrite reduction (anammox)
   REAL(wp), PUBLIC ::   e18o_xnio  !: O18 fractionation due to nitrite oxidation (anammox)
   REAL(wp), PUBLIC ::   e18oxy_res !: O18 fractionation (O2) due to aerobic respiration
   REAL(wp), PUBLIC ::   e18oxy_amo !: O18 fractionation (O2) due to ammonium oxidation
   REAL(wp), PUBLIC ::   e18oxy_nio !: O18 fractionation (O2) due to nitrite oxidation
   REAL(wp), PUBLIC ::   fb_h2ono2  !: fraction of biotic oxygen atom exchange between nitrite and water
  
   LOGICAL , PUBLIC ::   ln_newnitr !: New nitrification parameterisation

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitrnh4    !: ammonia oxidiation array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitrno2    !: nitrite oxidiation array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitr15nh4  !: nitrification array (NH4-->NO2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitr15no2  !: nitrification array (NO2-->NO3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitr18nh4  !: nitrification array (NH4-->NO2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitr18no2  !: nitrification array (NO2-->NO3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zabiox18no2  !: abiotic exchange of oxygen atoms (NO2<-->H2O)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr       !: total denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr15doc  !: denitrification array (DOC-->NH4)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitrno2    !: NO3-->NO2 denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitrno3    !: NO2-->N2 denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr15no3  !: denitrification array (NO3-->NO2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr15no2  !: denitrification array (NO2-->N2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr18no3  !: denitrification array (NO3-->NO2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr18no2  !: denitrification array (NO2-->N2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zaltrem      !: anaerobic remin without O2, NO3 or NO2
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zaltrem15    !: non O2, NO3 or NO2 anoxic remin (DOC--NH4) array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zanammox     !: anammox array (NH4+NO2-->N2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zanammox15no2!: anammox array (NH4 + 1.3NO2 --> N2 + 0.3NO3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zanammox15nh4!: anammox array (NH4 + 1.3NO2 --> N2 + 0.3NO3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zanammox15no3!: anammox array (NH4 + 1.3NO2 --> N2 + 0.3NO3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zanammox18no2!: anammox array (NH4 + 1.3NO2 --> N2 + 0.3NO3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zanammox18no3!: anammox array (NH4 + 1.3NO2 --> N2 + 0.3NO3)

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zrem.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremik, zremikc, zremikn, zremikp, zsiremin, zfact 
      REAL(wp) ::   zsatur, zsatur2, znusil, znusil2, zdep, zdepmin, zfactdep
      REAL(wp) ::   zbactfer, zolimit, zrfact2
      REAL(wp) ::   zammonic, zoxyremn, zoxyremp, znitrate2ton
      REAL(wp) ::   zosil, ztem, zolimic, zolimin, zolimip, zdenitrn, zdenitrp
      REAL(wp) ::   znh3, zlimaoan, zlimaoaf, zlimnobn, zlimnobf
      REAL(wp) ::   zr15_doc, zr15_no3, zr15_no2, zr15_nh4
      REAL(wp) ::   zr18_no3, zr18_no2, zr18_oxy
      REAL(wp) ::   d18Oh2o, tk, zph, k_h2ono2, e18o_eq
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj    ) :: ztempbac
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdepbac, zolimi, zdepprod, zfacsi, zfacsib, zdepeff, zfebact
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zolimi15
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_rem')
      !
      ! Initialisation of arrys
      zdepprod(:,:,:) = 1._wp
      zdepeff (:,:,:) = 0.3_wp
      ztempbac(:,:)   = 0._wp
      zfacsib(:,:,:)  = xsilab / ( 1.0 - xsilab )
      zfebact(:,:,:)  = 0._wp
      zfacsi(:,:,:)   = xsilab
      zolimi(:,:,:) = 0._wp

      IF ( ln_n15 ) THEN
         zolimi15(:,:,:)  = 0._wp
      ENDIF

      ! Computation of the mean phytoplankton concentration as
      ! a crude estimate of the bacterial biomass
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria
      ! -------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep = MAX( hmld(ji,jj), heup(ji,jj) )
               IF( gdept_n(ji,jj,jk) < zdep ) THEN
                  zdepbac(ji,jj,jk) = MIN( 0.7 * ( trb(ji,jj,jk,jpzoo) + 2.* trb(ji,jj,jk,jpmes) ), 4.e-6 )
                  ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
               ELSE
                  zdepmin = MIN( 1., zdep / gdept_n(ji,jj,jk) )
                  zdepbac (ji,jj,jk) = zdepmin**0.683 * ztempbac(ji,jj)
                  zdepprod(ji,jj,jk) = zdepmin**0.273
                  zdepeff (ji,jj,jk) = zdepeff(ji,jj,jk) * zdepmin**0.3
               ENDIF
            END DO
         END DO
      END DO

      IF( ln_p4z ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  zremik = xremik * xstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk) 
                  zremik = MAX( zremik, 2.74e-4 * xstep )

                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremik * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,jk,jpdoc) 
                  zolimi(ji,jj,jk) = MIN( ( trb(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimit ) 

                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremik * nitrfac(ji,jj,jk) * trb(ji,jj,jk,jpdoc) ! DOC remineralised by anaerobically
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) ) ! anaerobic remin due to denitrification

                  ! Do two-step denitrification
                  ! -------------------------------------------------------
                  znitrate2ton = ( trb(ji,jj,jk,jpno3)+rtrn ) / ( trb(ji,jj,jk,jpno3)+trb(ji,jj,jk,jpno2)+rtrn )
                  denitrno3(ji,jj,jk) = denitr(ji,jj,jk) * znitrate2ton ! NO3 --> NO2 denitrification
                  denitrno2(ji,jj,jk) = denitr(ji,jj,jk) * (1.0 - znitrate2ton) ! NO2 --> N2 denitrification
                  denitrno3(ji,jj,jk)  = MIN((trb(ji,jj,jk,jpno3)-rtrn) / rdenit, denitrno3(ji,jj,jk)) ! only remove available NO3
                  denitrno2(ji,jj,jk)  = MIN((trb(ji,jj,jk,jpno2)-rtrn) / rdenit, denitrno2(ji,jj,jk)) ! only remove available NO2

                  ! anaerobic remineralisation without NO3, NO2 or O2
                  ! -------------------------------------------------------
                  zaltrem(ji,jj,jk) = zammonic - denitr(ji,jj,jk)

                  ! make sure all arrays are positive
                  ! -------------------------------------------------------
                  zolimi(ji,jj,jk) = MAX( 0.e0, zolimi (ji,jj,jk) )
                  denitrno3(ji,jj,jk) = MAX( 0.e0, denitrno3(ji,jj,jk) )
                  denitrno2(ji,jj,jk) = MAX( 0.e0, denitrno2(ji,jj,jk) )
                  zaltrem(ji,jj,jk) = MAX( 0.e0, zaltrem(ji,jj,jk) )
                  
                  ! update the full denitrification array
                  ! -------------------------------------------------------
                  denitr(ji,jj,jk) = denitrno3(ji,jj,jk) + denitrno2(ji,jj,jk)

                  ! update the tracer arrays
                  ! -------------------------------------------------------
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jpno2) = tra(ji,jj,jk,jpno2) - denitrno2(ji,jj,jk) * rdenit + denitrno3(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitrno3(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimi(ji,jj,jk) - denitr(ji,jj,jk) - zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimi(ji,jj,jk) * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) +                                                &
                  &                     rno3 * ( zolimi(ji,jj,jk) + zaltrem(ji,jj,jk + denitr(ji,jj,jk) )    &
                  &                     + rdenit * denitrno2(ji,jj,jk) )  ! removal of NO2-->N2 increases alkalinity
                        ! Wolf-Gladrow et al. (2007) 
                        ! Alkalinity increases by 1 mol for every 1 mol NO3/NO2 removed
                        ! Alkalinity decreases by 1 mol for every 1 mol NH4 removed
                  IF( ln_n15) THEN
                     ! get isotopic signatures of major tracers
                     zr15_doc = ( (trb(ji,jj,jk,jp15doc)+rtrn) / (trb(ji,jj,jk,jpdoc)+rtrn) )
                     zr15_no3 = ( (trb(ji,jj,jk,jp15no3)+rtrn) / (trb(ji,jj,jk,jpno3)+rtrn) )
                     zr15_no2 = ( (trb(ji,jj,jk,jp15no2)+rtrn) / (trb(ji,jj,jk,jpno2)+rtrn) )
                     ! save flux arrays with isotopic effects
                     zolimi15(ji,jj,jk) = zolimi(ji,jj,jk) * ( 1.0 - e15n_amm/1000.0 ) * zr15_doc
                     denitr15no3(ji,jj,jk) = denitrno3(ji,jj,jk) * ( 1.0 - e15n_nar/1000.0 ) * zr15_no3
                     denitr15no2(ji,jj,jk) = denitrno2(ji,jj,jk) * ( 1.0 - e15n_nir/1000.0 ) * zr15_no2
                     denitr15doc(ji,jj,jk) = denitr(ji,jj,jk) * ( 1.0 - e15n_amm/1000.0 ) * zr15_doc
                     zaltrem15(ji,jj,jk) = zaltrem(ji,jj,jk) * ( 1.0 - e15n_amm/1000.0 ) * zr15_doc
                     ! update tracer arrays with these fluxes
                     tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4) + zolimi15(ji,jj,jk) + denitr15doc(ji,jj,jk) + zaltrem15(ji,jj,jk)
                     tra(ji,jj,jk,jp15no3) = tra(ji,jj,jk,jp15no3) - denitr15no3(ji,jj,jk) * rdenit
                     tra(ji,jj,jk,jp15no2) = tra(ji,jj,jk,jp15no2) + denitr15no3(ji,jj,jk) * rdenit - denitr15no2(ji,jj,jk) * rdenit
                     tra(ji,jj,jk,jp15doc) = tra(ji,jj,jk,jp15doc) - zolimi15(ji,jj,jk) - denitr15doc(ji,jj,jk) - zaltrem15(ji,jj,jk)
                  ENDIF
                  IF( ln_o18) THEN
                     ! get isotopic signatures of major tracers
                     zr18_oxy = ( (trb(ji,jj,jk,jp18oxy)+rtrn) / (trb(ji,jj,jk,jpoxy)+rtrn) )
                     zr18_no3 = ( (trb(ji,jj,jk,jp18no3)+rtrn) / (trb(ji,jj,jk,jpno3)+rtrn) )
                     zr18_no2 = ( (trb(ji,jj,jk,jp18no2)+rtrn) / (trb(ji,jj,jk,jpno2)+rtrn) )
                     ! save flux arrays with isotopic effects
                     denitr18no3(ji,jj,jk) = denitrno3(ji,jj,jk) * ( 1.0 - e18o_nar/1000.0 ) * zr18_no3
                     denitr18no2(ji,jj,jk) = denitrno2(ji,jj,jk) * ( 1.0 - e18o_nir/1000.0 ) * zr18_no2
                     ! update tracer arrays with these fluxes
                     tra(ji,jj,jk,jp18oxy) = tra(ji,jj,jk,jp18oxy) - zolimi(ji,jj,jk) * o2ut * ( 1. - e18oxy_res/1000. ) * zr18_oxy
                     tra(ji,jj,jk,jp18no3) = tra(ji,jj,jk,jp18no3) - denitr18no3(ji,jj,jk) * rdenit
                     tra(ji,jj,jk,jp18no2) = tra(ji,jj,jk,jp18no2) + denitr18no3(ji,jj,jk) * rdenit - denitr18no2(ji,jj,jk) * rdenit
                  ENDIF
               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  ! -----------------------------------------------------------------
                  zremik = xstep / 1.e-6 * MAX(0.01, xlimbac(ji,jj,jk)) * zdepbac(ji,jj,jk) 
                  zremik = MAX( zremik, 2.74e-4 * xstep / xremikc )

                  zremikc = xremikc * zremik
                  zremikn = xremikn / xremikc
                  zremikp = xremikp / xremikc

                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremikc * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,jk,jpdoc) 
                  zolimic = MAX( 0.e0, MIN( ( trb(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimit ) ) 
                  zolimi(ji,jj,jk) = zolimic
                  zolimin = zremikn * zolimic * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zolimip = zremikp * zolimic * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn ) 

                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremikc * nitrfac(ji,jj,jk) * trb(ji,jj,jk,jpdoc)
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
                  denitr(ji,jj,jk)  = MAX(0., MIN(  ( trb(ji,jj,jk,jpno3) - rtrn ) / rdenit, denitr(ji,jj,jk) ) )
                  zaltrem(ji,jj,jk) = MAX(0., zammonic - denitr(ji,jj,jk))
                  zdenitrn  = zremikn * denitr(ji,jj,jk) * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zdenitrp  = zremikp * denitr(ji,jj,jk) * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zoxyremn  = zremikn * zaltrem(ji,jj,jk) * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zoxyremp  = zremikp * zaltrem(ji,jj,jk) * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn )

                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimip + zdenitrp + zoxyremp
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimin + zdenitrn + zoxyremn
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitr(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimic - denitr(ji,jj,jk) - zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) - zolimin - zdenitrn - zoxyremn
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) - zolimip - zdenitrp - zoxyremp
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimic * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimic + denitr(ji,jj,jk) + zaltrem(ji,jj,jk)
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zolimin + zoxyremn + ( rdenit + 1.) * zdenitrn )
               END DO
            END DO
         END DO
         !
      ENDIF


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! NH4 nitrification to NO3. Ceased for oxygen concentrations
               ! below 2 umol/L. Inhibited at strong light 
               ! ----------------------------------------------------------
               zonitrnh4(ji,jj,jk)  = nitrif * xstep * trb(ji,jj,jk,jpnh4) * ( 1.- nitrfac(ji,jj,jk) )  &
               &                      / ( 1.+ emoy(ji,jj,jk) ) * ( 1. + fr_i(ji,jj) * emoy(ji,jj,jk) ) 
               zonitrno2(ji,jj,jk)  = nitrif * xstep * trb(ji,jj,jk,jpno2) * ( 1.- nitrfac(ji,jj,jk) )  &
               &                      / ( 1.+ emoy(ji,jj,jk) ) * ( 1. + fr_i(ji,jj) * emoy(ji,jj,jk) ) 

               ! New nitrification parameterisations
               if ( ln_newnitr ) then
                 ! Ammonia oxidation
                 !zph = (-1)*log10(hi(ji,jj,jk))
                 !znh3 = min(1.0, 10**(zph - 9.3) / 10**(8.0 - 9.3) ) 
                 zlimaoan = (trb(ji,jj,jk,jpnh4)+rtrn) / ( trb(ji,jj,jk,jpnh4) + kaoanh4 + rtrn)
                 zlimaoaf = (trb(ji,jj,jk,jpfer)+rtrn) / ( trb(ji,jj,jk,jpfer) + kaoafer + rtrn)
                 zonitrnh4(ji,jj,jk) = mu_aoa * xstep * trb(ji,jj,jk,jpnh4) * min(zlimaoan, zlimaoaf)   &
                 &                     * (1. - nitrfac(ji,jj,jk)) / ( 1. + emoy(ji,jj,jk) )
                 ! Nitrite oxidation
                 zlimnobn = (trb(ji,jj,jk,jpno2)+rtrn) / ( trb(ji,jj,jk,jpno2) + knobno2 + rtrn)
                 zlimnobf = (trb(ji,jj,jk,jpfer)+rtrn) / ( trb(ji,jj,jk,jpfer) + knobfer + rtrn)
                 zonitrno2(ji,jj,jk) = mu_nob * xstep * trb(ji,jj,jk,jpno2) * min(zlimnobn,zlimnobf)    &
                 &                     / ( 1. + emoy(ji,jj,jk) )
               endif

               ! Loss of NH4 and NO2 due to anammox
               ! ----------------------------------------------------------
               zanammox(ji,jj,jk) = ranammox * xstep * trb(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk)
               zanammox(ji,jj,jk) = max(0.0, min(  ( trb(ji,jj,jk,jpno2) - rtrn ) / 1.3, zanammox(ji,jj,jk) ) ) 
                        ! 1.3 mol of NO2 required per mol of NH4 oxidised by anammox (Brunner et al., 2013 PNAS)

               ! Make sure that the multiple sources and sinks of nitrite are not removing more NO2 than is available
               ! and if so, reduce nitrification and anammox by equal proportions (ignores denitrification above)
               ! ----------------------------------------------------------
               zfact = max(0.0, min(1.0, (trb(ji,jj,jk,jpno2) + tra(ji,jj,jk,jpno2) - rtrn) /  &
               &                         (zonitrno2(ji,jj,jk) + zanammox(ji,jj,jk)*1.3 + rtrn) ) )
               zonitrno2(ji,jj,jk) = zonitrno2(ji,jj,jk) * zfact
               zanammox(ji,jj,jk) = zanammox(ji,jj,jk) * zfact

               ! Update of the tracers trends
               ! ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zonitrnh4(ji,jj,jk) - zanammox(ji,jj,jk)
               tra(ji,jj,jk,jpno2) = tra(ji,jj,jk,jpno2) + zonitrnh4(ji,jj,jk) - zonitrno2(ji,jj,jk) - zanammox(ji,jj,jk)*1.3
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zonitrno2(ji,jj,jk) + zanammox(ji,jj,jk)*0.3
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zonitrnh4(ji,jj,jk)*o2nit*(3./4) - zonitrno2(ji,jj,jk)*o2nit*(1./4)
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2*rno3*zonitrnh4(ji,jj,jk)
                        ! Wolf-Gladrow et al. (2007) 
                        ! Alkalinity increases by 1 mol for every 1 mol NO3/NO2 removed
                        ! Alkalinity decreases by 1 mol for every 1 mol NH4 removed
                        ! ammonia-oxidation removes two mol Alk because -1 NH4, +1 NO2 = -1 Alk, -1 Alk
                        ! Anammox has no net effect because -1 NH4, -1.3 NO3 and +0.3 NO3 = -1 Alk, +1.3 Alk, -0.3 Alk
               IF( ln_n15 ) THEN
                  ! isotopic signatures of major tracers
                  zr15_nh4 = ( (trb(ji,jj,jk,jp15nh4)+rtrn) / (trb(ji,jj,jk,jpnh4)+rtrn) )
                  zr15_no2 = ( (trb(ji,jj,jk,jp15no2)+rtrn) / (trb(ji,jj,jk,jpno2)+rtrn) )
                  zr15_no3 = ( (trb(ji,jj,jk,jp15no3)+rtrn) / (trb(ji,jj,jk,jpno3)+rtrn) )
                  ! save fluxes
                  zonitr15nh4(ji,jj,jk) = zonitrnh4(ji,jj,jk) * zr15_nh4 * ( 1.0 - e15n_amo/1000.0 )
                  zonitr15no2(ji,jj,jk) = zonitrno2(ji,jj,jk) * zr15_no2 * ( 1.0 - e15n_nio/1000.0 )
                  zanammox15nh4(ji,jj,jk) = zanammox(ji,jj,jk) * zr15_nh4 * ( 1.0 - e15n_xamo/1000.0 )
                  zanammox15no2(ji,jj,jk) = zanammox(ji,jj,jk) * 1.3 * zr15_no2 * ( 1.0 - e15n_xnir/1000.0 )
                  zanammox15no3(ji,jj,jk) = zanammox(ji,jj,jk) * 0.3 * zr15_no2 * ( 1.0 - e15n_xnio/1000.0 )
                  ! update the arrays of major tracers with fluxes
                  tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4) - zonitr15nh4(ji,jj,jk) - zanammox15nh4(ji,jj,jk)
                  tra(ji,jj,jk,jp15no2) = tra(ji,jj,jk,jp15no2) + zonitr15nh4(ji,jj,jk) - zonitr15no2(ji,jj,jk)  &
                  &                       - zanammox15no2(ji,jj,jk)
                  tra(ji,jj,jk,jp15no3) = tra(ji,jj,jk,jp15no3) + zonitr15no2(ji,jj,jk) + zanammox15no3(ji,jj,jk)
               ENDIF
               IF( ln_o18 ) THEN
                  ! estimate d18O of seawater (coefficients from regression analysis of Legrande & Schmidt 2006 database)
                  d18Oh2o = (0.0333*tsn(ji,jj,jk,jp_tem) + 0.424*tsn(ji,jj,jk,jp_sal) - 14.8255)
                  ! isotopic signatures of major tracers
                  zr18_oxy = ( (trb(ji,jj,jk,jp18oxy)+rtrn) / (trb(ji,jj,jk,jpoxy)+rtrn) )
                  zr18_no2 = ( (trb(ji,jj,jk,jp18no2)+rtrn) / (trb(ji,jj,jk,jpno2)+rtrn) )

                  !! OLD EQUATIONS
                  !! save fluxes
                  !zonitr18nh4(ji,jj,jk) = zonitrnh4(ji,jj,jk) * zr18_no2 * ( 1.0 - e18o_amo/1000.0 ) ! need to change
                  !zonitr18no2(ji,jj,jk) = zonitrno2(ji,jj,jk) * zr18_no2 * ( 1.0 - e18o_nio/1000.0 )
                  !zanammox18no2(ji,jj,jk) = zanammox(ji,jj,jk) * 1.3 * zr18_no2 * ( 1.0 - e18o_xnir/1000.0 )
                  !zanammox18no3(ji,jj,jk) = zanammox(ji,jj,jk) * 0.3 * zr18_no2 * ( 1.0 - e18o_xnio/1000.0 )
                  !! update the arrays of major tracers with fluxes
                  !tra(ji,jj,jk,jp18oxy) = tra(ji,jj,jk,jp18oxy) - zr18_oxy                               &
                  !&                       * ( zonitrnh4(ji,jj,jk)*o2nit*(3./4)*(1.0 - e18oxy_amo/1000.)  &
                  !&                       + zonitrno2(ji,jj,jk)*o2nit*(1./4)*(1.0 - e18oxy_nio/1000.) ) 
                  !tra(ji,jj,jk,jp18no2) = tra(ji,jj,jk,jp18no2) + zonitr18nh4(ji,jj,jk) - zonitr18no2(ji,jj,jk)  &
                  !&                       - zanammox18no2(ji,jj,jk)
                  !tra(ji,jj,jk,jp18no3) = tra(ji,jj,jk,jp18no3) + zonitr18no2(ji,jj,jk) + zanammox18no3(ji,jj,jk)
 
                  ! NEW EQUATIONS
                  tk = tsn(ji,jj,jk,jp_tem)+273.15
                  zph = (-1)*log10(hi(ji,jj,jk))

                  ! Calculate the abiotic exchange rate (day-1) of nitrite and water 
                  ! and the associated equilibrium fractionation (Buchwald & Casciotti, 2013, Nat. Geo.)
                  k_h2ono2 = max(0.0, ( -0.0166 - 1.4123*exp( -(log(tk/317.71)/0.0462)**2.0 ) )            &
                  &          * zph + 0.1467 + 10.193*exp( -(log(tk/316.49)/0.0447)**2.0 ) ) 
                  e18o_eq = -0.12 * tk + 48.79 

                  ! Calculate the isotopic fluxes                   
                  zonitr18nh4(ji,jj,jk) = zonitrnh4(ji,jj,jk) * (                                          &
                  &                       0.5 * ( zr18_oxy + (1. + (d18Oh2o - e18o_o2_h2o)/1000.0 ) )      &
                  &                       * (1. - fb_h2ono2) + (1.0 + (d18Oh2o + e18o_eq)/1000.0) * fb_h2ono2 )
                  zabiox18no2(ji,jj,jk) = trb(ji,jj,jk,jpno2) * k_h2ono2 * xstep * (                       & 
                  &                       ( 1.0 + (d18Oh2o + e18o_eq)/1000. ) - zr18_no2 )
                  zonitr18no2(ji,jj,jk) = zonitrno2(ji,jj,jk) * ( 2/3. * (zr18_no2 - e18o_no2/1000.0)      &
                  &                       + 1/3. * (1.0 + (d18Oh2o - e18o_h2o_2)/1000.) )   ! Eq 3 Casciotti et al. 2010 L&O
                  zanammox18no2(ji,jj,jk) = zanammox(ji,jj,jk) * 1.3 * zr18_no2 * ( 1.0 - e18o_xnir/1000.0 )
                  zanammox18no3(ji,jj,jk) = zanammox(ji,jj,jk) * 0.3 * zr18_no2 * ( 1.0 - e18o_xnio/1000.0 )

                  ! update the arrays of major tracers with fluxes
                  tra(ji,jj,jk,jp18oxy) = tra(ji,jj,jk,jp18oxy) - zr18_oxy                               &
                  &                       * ( zonitrnh4(ji,jj,jk)*o2nit*(3./4)*(1.0 - e18oxy_amo/1000.)  &
                  &                       + zonitrno2(ji,jj,jk)*o2nit*(1./4)*(1.0 - e18oxy_nio/1000.) ) 
                  tra(ji,jj,jk,jp18no2) = tra(ji,jj,jk,jp18no2) + zonitr18nh4(ji,jj,jk) + zabiox18no2(ji,jj,jk)  &
                  &                       - zonitr18no2(ji,jj,jk) - zanammox18no2(ji,jj,jk)
                  tra(ji,jj,jk,jp18no3) = tra(ji,jj,jk,jp18no3) + zonitr18no2(ji,jj,jk) + zanammox18no3(ji,jj,jk)
               ENDIF
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! Bacterial uptake of iron. No iron is available in DOC. So
               ! Bacteries are obliged to take up iron from the water. Some
               ! studies (especially at Papa) have shown this uptake to be significant
               ! ----------------------------------------------------------
               zbactfer = feratb *  rfact2 * 0.6_wp / rday * tgfunc(ji,jj,jk) * xlimbacl(ji,jj,jk)     &
                  &              * trb(ji,jj,jk,jpfer) / ( xkferb + trb(ji,jj,jk,jpfer) )    &
                  &              * zdepprod(ji,jj,jk) * zdepeff(ji,jj,jk) * zdepbac(ji,jj,jk)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zbactfer*0.33
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zbactfer*0.25
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zbactfer*0.08
               zfebact(ji,jj,jk)   = zbactfer * 0.33
               blim(ji,jj,jk)      = xlimbacl(ji,jj,jk)  * zdepbac(ji,jj,jk) / 1.e-6 * zdepprod(ji,jj,jk)
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      ! Initialization of the array which contains the labile fraction
      ! of bSi. Set to a constant in the upper ocean
      ! ---------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep     = MAX( hmld(ji,jj), heup_01(ji,jj) )
               zsatur   = MAX( rtrn, ( sio3eq(ji,jj,jk) - trb(ji,jj,jk,jpsil) ) / ( sio3eq(ji,jj,jk) + rtrn ) )
               zsatur2  = ( 1. + tsn(ji,jj,jk,jp_tem) / 400.)**37
               znusil   = 0.225  * ( 1. + tsn(ji,jj,jk,jp_tem) / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25
               ! Remineralization rate of BSi depedant on T and saturation
               ! ---------------------------------------------------------
               IF ( gdept_n(ji,jj,jk) > zdep ) THEN
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk-1) * EXP( -0.5 * ( xsiremlab - xsirem )  &
                  &                   * znusil * e3t_n(ji,jj,jk) / wsbio4(ji,jj,jk) )
                  zfacsi(ji,jj,jk)  = zfacsib(ji,jj,jk) / ( 1.0 + zfacsib(ji,jj,jk) )
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk) * EXP( -0.5 * ( xsiremlab - xsirem )    &
                  &                   * znusil * e3t_n(ji,jj,jk) / wsbio4(ji,jj,jk) )
               ENDIF
               zsiremin = ( xsiremlab * zfacsi(ji,jj,jk) + xsirem * ( 1. - zfacsi(ji,jj,jk) ) ) * xstep * znusil
               zosil    = zsiremin * trb(ji,jj,jk,jpgsi)
               !
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) - zosil
               tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zosil
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF( knt == nrdttrc ) THEN
          zrfact2 = 1.e3 * rfact2r
          ALLOCATE( zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "REMIN" ) )  THEN
              zw3d(:,:,:) = zolimi(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMIN"  , zw3d )
          ENDIF
          IF( iom_use( "NITRNH4" ) )  THEN
              zw3d(:,:,:) = zonitrnh4(:,:,:) * rno3 * tmask(:,:,:) * zfact ! 1st step of nitrification
              CALL iom_put( "NITRNH4"  , zw3d )
          ENDIF
          IF( iom_use( "NITRNO2" ) )  THEN
              zw3d(:,:,:) = zonitrno2(:,:,:) * rno3 * tmask(:,:,:) * zfact ! 2nd step of nitrification
              CALL iom_put( "NITRNO2"  , zw3d )
          ENDIF
          IF( iom_use( "DENITNO3" ) )  THEN
              zw3d(:,:,:) = denitrno3(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENITNO3"  , zw3d )
          ENDIF
          IF( iom_use( "DENITNO2" ) )  THEN
              zw3d(:,:,:) = denitrno2(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENITNO2"  , zw3d )
          ENDIF
          IF( iom_use( "ANAMMOX" ) )  THEN
              zw3d(:,:,:) = zanammox(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "ANAMMOX"  , zw3d )
          ENDIF
          IF( iom_use( "ALTREM" ) )  THEN
              zw3d(:,:,:) = zaltrem(:,:,:) * tmask(:,:,:) * zfact ! anaerobic remin (DOC-->NH4) without O2,NO3,NO2
              CALL iom_put( "ALTREM"  , zw3d )
          ENDIF
          IF( iom_use( "BACT" ) )  THEN
               zw3d(:,:,:) = zdepbac(:,:,:) * 1.E6 * tmask(:,:,:)  ! Bacterial biomass
               CALL iom_put( "BACT", zw3d )
          ENDIF
          IF( iom_use( "FEBACT" ) )  THEN
               zw3d(:,:,:) = zfebact(:,:,:) * 1E9 * tmask(:,:,:) * zrfact2   ! Bacterial iron consumption
               CALL iom_put( "FEBACT" , zw3d )
          ENDIF
          !
          DEALLOCATE( zw3d )
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, nitrif, xsirem, xsiremlab, xsilab, feratb, xkferb,   & 
         &                xremikc, xremikn, xremikp, mu_aoa, kaoanh4, kaoafer, mu_nob, &
         &                knobno2, knobfer, ranammox, e15n_nar, e15n_nir, e15n_amo,    &
         &                e15n_nio, e15n_amm, e15n_xamo, e15n_xnir, e15n_xnio, e18o_nar&
         &                , e18o_nir, e18o_o2_h2o, e18o_h2o_2, e18o_no2, e18o_xnir,    & 
         &                e18o_xnio, e18oxy_res, e18oxy_amo, e18oxy_nio, fb_h2ono2, ln_newnitr
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of remineralization parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisrem in reference namelist : Pisces remineralization
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisrem in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisrem in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisrem in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist parameters for remineralization, nampisrem'
         IF( ln_p4z ) THEN
            WRITE(numout,*) '      remineralization rate of DOC              xremik    =', xremik
         ELSE
            WRITE(numout,*) '      remineralization rate of DOC              xremikc   =', xremikc
            WRITE(numout,*) '      remineralization rate of DON              xremikn   =', xremikn
            WRITE(numout,*) '      remineralization rate of DOP              xremikp   =', xremikp
         ENDIF
         WRITE(numout,*) '      remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '      fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '      fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '      NH4 nitrification rate                    nitrif    =', nitrif
         WRITE(numout,*) '      Bacterial Fe/C ratio                      feratb    =', feratb
         WRITE(numout,*) '      Half-saturation constant for bact. Fe/C   xkferb    =', xkferb
         WRITE(numout,*) '      Growth rate of ammonia-oxidising archaea  mu_aoa    =', mu_aoa
         WRITE(numout,*) '      NH4 half-saturation constant for AOA      kaoanh4   =', kaoanh4
         WRITE(numout,*) '      Fer half-saturation constant for AOA      kaoafer   =', kaoafer
         WRITE(numout,*) '      Growth rate of nitrite-oxidising bacteria mu_nob    =', mu_nob
         WRITE(numout,*) '      NO2 half-saturation constant for NOB      knobno2   =', knobno2
         WRITE(numout,*) '      Fer half-saturation constant for NOB      knobfer   =', knobfer
         WRITE(numout,*) '      Max rate of anammox (per day)             ranammox  =', ranammox
         WRITE(numout,*) '      N15 fractionation - nitrate reductase     e15n_nar  =', e15n_nar
         WRITE(numout,*) '      N15 fractionation - nitrite reductase     e15n_nir  =', e15n_nir
         WRITE(numout,*) '      N15 fractionation - ammonium oxidation    e15n_amo  =', e15n_amo
         WRITE(numout,*) '      N15 fractionation - nitrite oxidation     e15n_nio  =', e15n_nio
         WRITE(numout,*) '      N15 fractionation - ammonification        e15n_amm  =', e15n_amm
         WRITE(numout,*) '      N15 fractionation - anammox NH4 oxidation e15n_xamo =', e15n_xamo
         WRITE(numout,*) '      N15 fractionation - anammox NO2 reduction e15n_xnir =', e15n_xnir
         WRITE(numout,*) '      N15 fractionation - anammox NO2 oxidation e15n_xnio =', e15n_xnio
         WRITE(numout,*) '      O18 fractionation - nitrate reductase     e18o_nar  =', e18o_nar
         WRITE(numout,*) '      O18 fractionation - nitrite reductase     e18o_nir  =', e18o_nir
         WRITE(numout,*) '      O18 fractionation - O incorp (O2+H2O)     e18o_o2_h2o =', e18o_o2_h2o
         WRITE(numout,*) '      O18 fractionation - O incorp (H2O)        e18o_h2o_2=', e18o_h2o_2
         WRITE(numout,*) '      O18 fractionation - O selection of N2O    e18o_no2  =', e18o_no2
         WRITE(numout,*) '      O18 fractionation - anammox NO2 reduction e18o_xnir =', e18o_xnir
         WRITE(numout,*) '      O18 fractionation - anammox NO2 oxidation e18o_xnio =', e18o_xnio
         WRITE(numout,*) '      Logical for new nitrification param       ln_newnitr=', ln_newnitr 
         WRITE(numout,*) '      O18 fractionation (O2)- respiration       e18oxy_res=', e18oxy_res
         WRITE(numout,*) '      O18 fractionation (O2)- ammonia oxidation e18oxy_amo=', e18oxy_amo
         WRITE(numout,*) '      O18 fractionation (O2)- nitrite oxidation e18oxy_nio=', e18oxy_nio
         WRITE(numout,*) '      fraction biotic O-atom exchange NO2-H2O   fb_h2ono2 =', fb_h2ono2
      ENDIF
      !
      zonitrnh4(:,:,:) = 0._wp
      zonitrno2(:,:,:) = 0._wp
      zonitr15nh4(:,:,:) = 0._wp
      zonitr15no2(:,:,:) = 0._wp
      zonitr18nh4(:,:,:) = 0._wp
      zonitr18no2(:,:,:) = 0._wp
      zabiox18no2(:,:,:) = 0._wp
      denitr(:,:,:) = 0._wp
      denitr15doc(:,:,:) = 0._wp
      denitrno3(:,:,:) = 0._wp
      denitrno2(:,:,:) = 0._wp
      denitr15no2(:,:,:) = 0._wp
      denitr15no3(:,:,:) = 0._wp
      denitr18no2(:,:,:) = 0._wp
      denitr18no3(:,:,:) = 0._wp
      zaltrem(:,:,:) = 0._wp
      zaltrem15(:,:,:) = 0._wp
      zanammox(:,:,:) = 0._wp
      zanammox15no2(:,:,:) = 0._wp
      zanammox15nh4(:,:,:) = 0._wp
      zanammox15no3(:,:,:) = 0._wp
      zanammox18no2(:,:,:) = 0._wp
      zanammox18no3(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( zonitrnh4(jpi,jpj,jpk), zonitrno2(jpi,jpj,jpk), denitr(jpi,jpj,jpk),   &
      &         denitrno3(jpi,jpj,jpk), denitrno2(jpi,jpj,jpk), zaltrem(jpi,jpj,jpk),  &
      &         zanammox(jpi,jpj,jpk),                                                 &
      &         denitr15doc(jpi,jpj,jpk), denitr15no3(jpi,jpj,jpk),                    &
      &         denitr15no2(jpi,jpj,jpk), zonitr15nh4(jpi,jpj,jpk),                    &
      &         denitr18no2(jpi,jpj,jpk), denitr18no3(jpi,jpj,jpk),                    &
      &         zonitr15no2(jpi,jpj,jpk), zaltrem15(jpi,jpj,jpk),                      &
      &         zonitr18no2(jpi,jpj,jpk), zonitr18nh4(jpi,jpj,jpk),                    &
      &         zanammox15nh4(jpi,jpj,jpk), zanammox15no2(jpi,jpj,jpk),                &
      &         zanammox15no3(jpi,jpj,jpk), zanammox18no2(jpi,jpj,jpk),                &
      &         zanammox18no3(jpi,jpj,jpk), zabiox18no2(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_rem_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_rem_alloc

   !!======================================================================
END MODULE p4zrem
