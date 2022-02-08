MODULE p4zmeso
   !!======================================================================
   !!                         ***  MODULE p4zmeso  ***
   !! TOP :   PISCES Compute the sources/sinks for mesozooplankton
   !!======================================================================
   !! History :   1.0  !  2002     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_meso       : Compute the sources/sinks for mesozooplankton
   !!   p4z_meso_init  : Initialization of the parameters for mesozooplankton
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zprod         ! production
   USE prtctl_trc      ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_meso              ! called in p4zbio.F90
   PUBLIC   p4z_meso_init         ! called in trcsms_pisces.F90

   REAL(wp), PUBLIC ::  part2        !: part of calcite not dissolved in mesozoo guts
   REAL(wp), PUBLIC ::  xpref2d      !: mesozoo preference for diatoms
   REAL(wp), PUBLIC ::  xpref2n      !: mesozoo preference for nanophyto
   REAL(wp), PUBLIC ::  xpref2z      !: mesozoo preference for microzooplankton
   REAL(wp), PUBLIC ::  xpref2c      !: mesozoo preference for POC 
   REAL(wp), PUBLIC ::  xthresh2zoo  !: zoo feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2dia  !: diatoms feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2phy  !: nanophyto feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2poc  !: poc feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  xthresh2     !: feeding threshold for mesozooplankton 
   REAL(wp), PUBLIC ::  resrat2      !: exsudation rate of mesozooplankton
   REAL(wp), PUBLIC ::  mzrat2       !: microzooplankton mortality rate 
   REAL(wp), PUBLIC ::  grazrat2     !: maximal mesozoo grazing rate
   REAL(wp), PUBLIC ::  xkgraz2      !: non assimilated fraction of P by mesozoo 
   REAL(wp), PUBLIC ::  unass2       !: Efficicency of mesozoo growth 
   REAL(wp), PUBLIC ::  sigma2       !: Fraction of mesozoo excretion as DOM 
   REAL(wp), PUBLIC ::  epsher2      !: growth efficiency
   REAL(wp), PUBLIC ::  epsher2min   !: minimum growth efficiency at high food for grazing 2
   REAL(wp), PUBLIC ::  grazflux     !: mesozoo flux feeding rate
   REAL(wp), PUBLIC ::  e13c_cal2    !: C13 mesozoo calcification fractionation
   REAL(wp), PUBLIC ::  e15n_ex2     !: N15 mesozoo excretion fractionation
   REAL(wp), PUBLIC ::  e15n_in2     !: N15 mesozoo ingestion fractionation
   REAL(wp), PUBLIC ::  e18oxy_mes   !: O18 mesozoo respiration fractionation

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zmeso.F90 10367 2018-12-03 11:35:19Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_meso( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_meso  ***
      !!
      !! ** Purpose :   Compute the sources/sinks for mesozooplankton
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step and ???
      !
      INTEGER  :: ji, jj, jk
      REAL(wp) :: zcompadi, zcompaph, zcompapoc, zcompaz, zcompam
      REAL(wp) :: zgraze2 , zdenom, zdenom2
      REAL(wp) :: zfact   , zfood, zfoodlim, zproport, zbeta
      REAL(wp) :: zmortzgoc, zfrac, zfracfe, zratio, zratio2, zfracal, zgrazcal
      REAL(wp) :: zepsherf, zepshert, zepsherv, zgrarsig, zgraztotc, zgraztotn, zgraztotf
      REAL(wp) :: zgrarem2, zgrafer2, zgrapoc2, zprcaca, zmortz, zgrasrat, zgrasratn
      REAL(wp) :: zrespz, ztortz, zgrazd, zgrazz, zgrazpof
      REAL(wp) :: zgrazn, zgrazpoc, zgraznf, zgrazf
      REAL(wp) :: zgrazfffp, zgrazfffg, zgrazffep, zgrazffeg
      REAL(wp) :: zgrazn13, zgrazd13, zgrazz13, zgrazpoc13, zgrazffp13, zgrazffg13
      REAL(wp) :: zmortz13, zfrac13, zgraztotc13
      REAL(wp) :: zgrarem2_13, zgrapoc2_13, zgrasig2_13, zmortzgoc_13, zr13_dic, zr13_cal
      REAL(wp) :: zgrazn15, zgrazd15, zgrazz15, zgrazpoc15, zgrazffp15, zgrazffg15
      REAL(wp) :: zmortz15, zfrac15, zgraztotc15
      REAL(wp) :: zgrarem2_15, zgrapoc2_15, zgrasig2_15, zgrasigex2_15, zmortzgoc_15
      REAL(wp) :: zr18_oxy
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zgrazing2, zfezoo2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: excretion2, excretion2_13, excretion2_15
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zw3d, zz2ligprod
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_meso')
      !
      zgrazing2(:,:,:) = 0._wp 
      zfezoo2 (:,:,:) = 0._wp
      excretion2(:,:,:) = 0._wp 
      IF( ln_n15 ) THEN
        excretion2_15(:,:,:) = 0._wp
      ENDIF
      !
      IF (ln_ligand) THEN
         ALLOCATE( zz2ligprod(jpi,jpj,jpk) )
         zz2ligprod(:,:,:) = 0._wp
      ENDIF
      !
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zcompam   = MAX( ( trb(ji,jj,jk,jpmes) - 1.e-9 ), 0.e0 )
               zfact     = xstep * tgfunc2(ji,jj,jk) * zcompam

               !  Respiration rates of both zooplankton
               !  -------------------------------------
               zrespz    = resrat2 * zfact * ( trb(ji,jj,jk,jpmes) / ( xkmort + trb(ji,jj,jk,jpmes) )  &
               &           + 3. * nitrfac(ji,jj,jk) )

               !  Zooplankton mortality. A square function has been selected with
               !  no real reason except that it seems to be more stable and may mimic predation
               !  ---------------------------------------------------------------
               ztortz    = mzrat2 * 1.e6 * zfact * trb(ji,jj,jk,jpmes)  * (1. - nitrfac(ji,jj,jk) )
               !
               zcompadi  = MAX( ( trb(ji,jj,jk,jpdia) - xthresh2dia ), 0.e0 )
               zcompaz   = MAX( ( trb(ji,jj,jk,jpzoo) - xthresh2zoo ), 0.e0 )
               zcompapoc = MAX( ( trb(ji,jj,jk,jppoc) - xthresh2poc ), 0.e0 )
               ! Size effect of nanophytoplankton on grazing : the smaller it is, the less prone
               ! it is to predation by mesozooplankton
               ! -------------------------------------------------------------------------------
               zcompaph  = MAX( ( trb(ji,jj,jk,jpphy) - xthresh2phy ), 0.e0 ) &
                  &      * MIN(1., MAX( 0., ( quotan(ji,jj,jk) - 0.2) / 0.3 ) )

               !   Mesozooplankton grazing
               !   ------------------------
               zfood     = xpref2d * zcompadi + xpref2z * zcompaz + xpref2n * zcompaph + xpref2c * zcompapoc 
               zfoodlim  = MAX( 0., zfood - MIN( 0.5 * zfood, xthresh2 ) )
               zdenom    = zfoodlim / ( xkgraz2 + zfoodlim )
               zdenom2   = zdenom / ( zfood + rtrn )

               ! get grazing rate (increase with temperature and decrease with low O2)
               zgraze2   = grazrat2 * xstep * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpmes) * (1. - nitrfac(ji,jj,jk)) 

               zgrazd    = zgraze2  * xpref2d  * zcompadi  * zdenom2 
               zgrazz    = zgraze2  * xpref2z  * zcompaz   * zdenom2 
               zgrazn    = zgraze2  * xpref2n  * zcompaph  * zdenom2 
               zgrazpoc  = zgraze2  * xpref2c  * zcompapoc * zdenom2 

               zgraznf   = zgrazn   * trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) + rtrn)
               zgrazf    = zgrazd   * trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) + rtrn)
               zgrazpof  = zgrazpoc * trb(ji,jj,jk,jpsfe) / ( trb(ji,jj,jk,jppoc) + rtrn)

               !  Mesozooplankton flux feeding on GOC
               !  ----------------------------------
               zgrazffeg = grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpgoc) * trb(ji,jj,jk,jpmes) &
               &           * (1. - nitrfac(ji,jj,jk))
               zgrazffep = grazflux  * xstep *  wsbio3(ji,jj,jk)     &
               &           * tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpmes) &
               &           * (1. - nitrfac(ji,jj,jk))
               ! grazflux = flux feeding rate (per m molC/L)
               ! wsbio4(ji,jj,jk) * trb(ji,jj,jk,jpgoc) * xstep = rate of GOC flux (m molC/L)
               ! tgfunc2(ji,jj,jk) * trb(ji,jj,jk,jpmes) = growth rate of zooplankton (molC/L)
               ! (1. - nitrfac(ji,jj,jk)) = no grazing without oxygen
               ! zgrazffeg = assimilation of big particles by flux feeding (molC/L)
               ! zgrazffep = assimilation of small particles by flux feeding (molC/L)

               ! Iron flux feeding
               zgrazfffp = zgrazffep * trb(ji,jj,jk,jpsfe) / (trb(ji,jj,jk,jppoc) + rtrn)
               zgrazfffg = zgrazffeg * trb(ji,jj,jk,jpbfe) / (trb(ji,jj,jk,jpgoc) + rtrn)

               ! Total feeding of Mesozooplankton [POC,GOC,PHY,PHY2,ZOO] --> [ZOO2]
               zgraztotc = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg ! total feeding

               ! get the d13C and d15N signature of the food sources
               IF ( ln_c13 ) THEN
                  zgrazn13 = zgrazn * ( (trb(ji,jj,jk,jp13phy)+rtrn) / (trb(ji,jj,jk,jpphy)+rtrn) )
                  zgrazd13 = zgrazd * ( (trb(ji,jj,jk,jp13dia)+rtrn) / (trb(ji,jj,jk,jpdia)+rtrn) )
                  zgrazz13 = zgrazz * ( (trb(ji,jj,jk,jp13zoo)+rtrn) / (trb(ji,jj,jk,jpzoo)+rtrn) )
                  zgrazpoc13 = zgrazpoc * ( (trb(ji,jj,jk,jp13poc)+rtrn) / (trb(ji,jj,jk,jppoc)+rtrn) )
                  zgrazffg13 = zgrazffeg * ( (trb(ji,jj,jk,jp13goc)+rtrn) / (trb(ji,jj,jk,jpgoc)+rtrn) )
                  zgrazffp13 = zgrazffep * ( (trb(ji,jj,jk,jp13poc)+rtrn) / (trb(ji,jj,jk,jppoc)+rtrn) )
               ENDIF
               IF ( ln_n15 ) THEN
                  zgrazn15 = zgrazn * ( (trb(ji,jj,jk,jp15phy)+rtrn) / (trb(ji,jj,jk,jpphy)+rtrn) )
                  zgrazd15 = zgrazd * ( (trb(ji,jj,jk,jp15dia)+rtrn) / (trb(ji,jj,jk,jpdia)+rtrn) )
                  zgrazz15 = zgrazz * ( (trb(ji,jj,jk,jp15zoo)+rtrn) / (trb(ji,jj,jk,jpzoo)+rtrn) )
                  zgrazpoc15 = zgrazpoc * ( (trb(ji,jj,jk,jp15poc)+rtrn) / (trb(ji,jj,jk,jppoc)+rtrn) )
                  zgrazffg15 = zgrazffeg * ( (trb(ji,jj,jk,jp15goc)+rtrn) / (trb(ji,jj,jk,jpgoc)+rtrn) )
                  zgrazffp15 = zgrazffep * ( (trb(ji,jj,jk,jp15poc)+rtrn) / (trb(ji,jj,jk,jppoc)+rtrn) )
               ENDIF

               !!! COMPUTE THE FRACTIONATION OF POC INTO GOC BY FILTER FEEDING [GOC] --> [POC]
               ! compute the proportion of filter feeding of total feeding rate
               zproport  = (zgrazffep + zgrazffeg + rtrn)/(zgraztotc + rtrn)
               ! Compute fractionation of aggregates. It is assumed that 
               ! diatoms based aggregates are more prone to fractionation
               ! since they are more porous (marine snow instead of fecal pellets)
               zratio    = trb(ji,jj,jk,jpgsi) / ( trb(ji,jj,jk,jpgoc) + rtrn )
               zratio2   = zratio * zratio
               zfrac     = zproport * grazflux  * xstep * wsbio4(ji,jj,jk)      &
               &          * trb(ji,jj,jk,jpgoc) * trb(ji,jj,jk,jpmes)          &
               &          * ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 ) )
               ! zproport = proportion of flux/filter feeder within prey food bank
               ! grazflux * trb(ji,jj,jk,jpmes) = bulk flux feeding rate of meso
               ! wsbio4(ji,jj,jk) * trb(ji,jj,jk,jpgoc) * xstep = rate of GOC flux ariving at depth level (m molC/L)
               ! ( 0.2 + 3.8 * zratio2 / ( 1.**2 + zratio2 )) = random parameterisation between 0.2 and 2
               ! zfrac = amount of GOC carbon converted to POC carbon

               ! Iron
               zfracfe   = zfrac * trb(ji,jj,jk,jpbfe) / (trb(ji,jj,jk,jpgoc) + rtrn)

               ! Multiply proportion of filter feeders by the grazing rate?
               zgrazffep = zproport * zgrazffep
               zgrazffeg = zproport * zgrazffeg
               zgrazfffp = zproport * zgrazfffp
               zgrazfffg = zproport * zgrazfffg
               zgraztotc = zgrazd + zgrazz + zgrazn + zgrazpoc + zgrazffep + zgrazffeg
               zgraztotn = zgrazd * quotad(ji,jj,jk) + zgrazz + zgrazn * quotan(ji,jj,jk)   &
               &   + zgrazpoc + zgrazffep + zgrazffeg
               zgraztotf = zgrazf + zgraznf + zgrazz * ferat3 + zgrazpof + zgrazfffp + zgrazfffg

               ! calculate the fractionation of particulates by filter feeders
               IF ( ln_c13 ) THEN
                  zfrac13     = zfrac * ( (trb(ji,jj,jk,jp13goc)+rtrn) / (trb(ji,jj,jk,jpgoc)+rtrn) )
                  zgrazffg13  = zgrazffg13 * zproport
                  zgrazffp13  = zgrazffp13 * zproport
                  zgraztotc13 = zgrazn13 + zgrazd13 + zgrazz13 + zgrazpoc13 + zgrazffp13 + zgrazffg13
               ENDIF
               IF ( ln_n15 ) THEN
                  zfrac15     = zfrac * ( (trb(ji,jj,jk,jp15goc)+rtrn) / (trb(ji,jj,jk,jpgoc)+rtrn) )
                  zgrazffg15  = zgrazffg15 * zproport
                  zgrazffp15  = zgrazffp15 * zproport
                  zgraztotc15 = zgrazn15 + zgrazd15 + zgrazz15 + zgrazpoc15 + zgrazffp15 + zgrazffg15
               ENDIF

               ! Total grazing ( grazing by microzoo is already computed in p4zmicro )
               zgrazing2(ji,jj,jk) = zgraztotc

               !    Mesozooplankton efficiency
               !    --------------------------
               zgrasrat  =  ( zgraztotf + rtrn )/ ( zgraztotc + rtrn )
               zgrasratn =  ( zgraztotn + rtrn )/ ( zgraztotc + rtrn )
               zepshert  = MIN( 1., zgrasratn, zgrasrat / ferat3)
               zbeta     = MAX(0., (epsher2 - epsher2min) )
               zepsherf  = epsher2min + zbeta / ( 1.0 + 0.04E6 * 12. * zfood * zbeta ) 
               zepsherv  = zepsherf * zepshert ! [0,0.5], 0.0 --> 0.35 in p4z
                 !  zepshert = quality of food accounting for Fe and N quotas [0,1]
                 !  zepsherf = efficiency of meso growth [0.2,0.5], 0.35 in p4z
                 !  zepsherv = proportion of food assimilated into MES
                 !             [POC,GOC,PHY,PHY2,ZOO] --> [ZOO2] (max 0.35 in p4z)
                 !  unass2   = proportion of food not assimilated into MES
                 !             [POC,GOC,PHY,PHY2,ZOO] --> [POC] (0.3)
                 !  epsher2  = proportion of food assimilated into MES if  maximum efficiency achieved
                 !             [POC,GOC,PHY,PHY2,ZOO] --> [ZOO2] (0.35)

               zgrarem2  = zgraztotc * ( 1. - zepsherv - unass2 ) &  ![food] --> [NH4/DOC]
               &         + ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 ) * ztortz ![mes] --> [NH4/DOC]
               zgrafer2  = zgraztotc * MAX( 0. , ( 1. - unass2 ) * zgrasrat - ferat3 * zepsherv )    &
               &         + ferat3 * ( ( 1. - epsher2 - unass2 ) /( 1. - epsher2 ) * ztortz )
               zgrapoc2  = zgraztotc * unass2 ![food] --> [GOC]

                 ! ( 1. - zepsherv - unass2 ) = proportion of food excreted [0.35 --> 0.7 in P4Z]
                 ! ( 1. - epsher2 - unass2 ) = proportion of food excreted assuming maximum efficiency [0.35]
                 ! ( 1. - epsher2 ) = proportion of food not assimilated (excreted + not consumed) max efficiency [0.65]
                 ! ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 ) = ratio of food excreted to not assimilated
                 !                                                (excreted + not consumed) assuming max efficiency [0.538]

                 ! zgraztotc * ( 1. - zepsherv - unass2 )                = total food carbon excreted
                 ! ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 ) * ztortz = mortality that goes to excretion
                 !                                                         (because not assimilating everything)
                 ! ( unass2 / ( 1. - epsher2 ) = unassimilated C : assimilated + excreted C (0.462)
                 ! ( unass2 / ( 1. - epsher2 ) * ztortz = mortality that goes to large particles (POC)

                 !  zgrarem2 = carbon excreted during feeding |+| mortality
                 !             [POC,GOC,PHY,PHY2,ZOO] --> [DOC,NH4] |+| [ZOO2] --> [DOC,NH4]
                 !  zgrapoc2 = carbon unassimilated during mesozooplankton eating [POC,GOC,PHY,PHY2,ZOO] --> [GOC]
                 !  zgrasig2 = carbon remineralised directly into NH4 (bulk term, includes direct remin + excretion )
                 !  zgrasigex2 = food excreted to  directly into NH4 via excretion
                 !               (accounts for food quality and achieved growth efficiency)
                 !               (excretion is greater with poor food quality and unrealised growth efficiency)
                 !  zmotrzgoc = carbon formed into large particles from zooplankton mortality and respiration
                 !              [ZOO2] --> [GOC]

               ! calculate the remin (-->NH4/DOC), remin[excretion] (-->NH4), unassim (-->GOC), mort (-->GOC)
               IF ( ln_c13 ) THEN
                  zgrarem2_13   = zgraztotc13 * ( 1. - zepsherv - unass2 ) &  ! [food] --> [NH4/DOC]
                  &               + ztortz * ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 )  & ![mes] --> [NH4/DOC]
                  &               * ( (trb(ji,jj,jk,jp13mes)+rtrn) / (trb(ji,jj,jk,jpmes)+rtrn) )
                  zgrapoc2_13   = zgraztotc13 * unass2 ! [food] --> [GOC]
                  zgrasig2_13   = zgrarem2_13 * sigma2 ! [food/mes] --> [NH4]
                  zmortzgoc_13  = ( unass2 / ( 1. - epsher2 ) * ztortz + zrespz )  & ! [mes] --> [GOC]
                  &               * ( (trb(ji,jj,jk,jp13mes)+rtrn) / (trb(ji,jj,jk,jpmes)+rtrn) )
               ENDIF
               IF ( ln_n15 ) THEN
                  zgrarem2_15   = zgraztotc15 * ( 1. - zepsherv - unass2 ) &  ! [food] --> [NH4/DOC]
                  &               + ztortz * ( 1. - epsher2 - unass2 ) / ( 1. - epsher2 )  & ![mes] --> [NH4/DOC]
                  &               * ( (trb(ji,jj,jk,jp15mes)+rtrn) / (trb(ji,jj,jk,jpmes)+rtrn) )
                  zgrapoc2_15   = zgraztotc15 * unass2 ! [food] --> [GOC]
                  zgrasig2_15   = zgrarem2_15 * sigma2 ! [food/mes] --> [NH4]
                  zgrasigex2_15 = zgraztotc15 * ( 1. - epsher2 - unass2 ) * sigma2 * zgrasratn ![food] --> [NH4]
                  zmortzgoc_15  = ( unass2 / ( 1. - epsher2 ) * ztortz + zrespz )  & ! [mes] --> [GOC]
                  &               * ( (trb(ji,jj,jk,jp15mes)+rtrn) / (trb(ji,jj,jk,jpmes)+rtrn) )
               ENDIF

               !   Update the arrays TRA which contain the biological sources and sinks
               zgrarsig  = zgrarem2 * sigma2 ! 60% of carbon excreted to NH4 [food/mes] --> [NH4]
               tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zgrarsig
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zgrarsig
               excretion2(ji,jj,jk) = zgrarsig
               tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zgrarem2 - zgrarsig
               !
               IF( ln_ligand ) THEN 
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + (zgrarem2 - zgrarsig) * ldocz
                  zz2ligprod(ji,jj,jk) = (zgrarem2 - zgrarsig) * ldocz
               ENDIF
               !
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2ut * zgrarsig
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + zgrafer2
               zfezoo2(ji,jj,jk)   = zgrafer2
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrarsig
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zgrarsig              

               IF ( ln_c13 ) THEN
                  tra(ji,jj,jk,jp13doc) = tra(ji,jj,jk,jp13doc) + zgrarem2_13 - zgrasig2_13
                  tra(ji,jj,jk,jp13dic) = tra(ji,jj,jk,jp13dic) + zgrasig2_13
                  excretion2_13(ji,jj,jk) = zgrasig2_13
               ENDIF
               IF ( ln_n15 ) THEN
                  excretion2_15(ji,jj,jk) = zgrasigex2_15 * ( 1. - e15n_ex2/1000.0)  + ( zgrasig2_15 - zgrasigex2_15 )
                  tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4) + excretion2_15(ji,jj,jk)
                  ! zgrasigex2_15 is always less than zgrasig2_15, because:
                  !     (1) it does not include the excretion of NH4 during mortality
                  !     (2) although it assumes maximum efficiency of assimilation (epsher2),
                  !         it multiples this by zgrazratn, which is [0,1] and is a measure of the food quality via N content
                  ! So, zgrasigex is the amount of carbon excreted to NH4 due to N content of food
                  tra(ji,jj,jk,jp15doc) = tra(ji,jj,jk,jp15doc) + zgrarem2_15 - zgrasig2_15
               ENDIF
               IF ( ln_o18 ) THEN
                  zr18_oxy = ( trb(ji,jj,jk,jp18oxy) + rtrn ) / ( trb(ji,jj,jk,jpoxy) + rtrn )
                  tra(ji,jj,jk,jp18oxy) = tra(ji,jj,jk,jp18oxy) - o2ut * zgrarsig * zr18_oxy * (1. - e18oxy_mes/1000.)
               ENDIF

               zmortz = ztortz + zrespz
               zmortzgoc = unass2 / ( 1. - epsher2 ) * ztortz + zrespz
               tra(ji,jj,jk,jpmes) = tra(ji,jj,jk,jpmes) - zmortz + zepsherv * zgraztotc 
               tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) - zgrazd
               tra(ji,jj,jk,jpzoo) = tra(ji,jj,jk,jpzoo) - zgrazz
               tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) - zgrazn
               tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) - zgrazn * trb(ji,jj,jk,jpnch) / ( trb(ji,jj,jk,jpphy) + rtrn )
               tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) - zgrazd * trb(ji,jj,jk,jpdch) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) - zgrazd * trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) + zgrazd * trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn )
               tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) - zgraznf
               tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) - zgrazf

               tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zgrazpoc - zgrazffep + zfrac
               prodpoc(ji,jj,jk) = prodpoc(ji,jj,jk) + zfrac
               conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zgrazpoc - zgrazffep
               tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zmortzgoc - zgrazffeg + zgrapoc2 - zfrac
               prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zmortzgoc + zgrapoc2
               consgoc(ji,jj,jk) = consgoc(ji,jj,jk) - zgrazffeg - zfrac
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zgrazpof - zgrazfffp + zfracfe
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + ferat3 * zmortzgoc - zgrazfffg     &
                 &                + zgraztotf * unass2 - zfracfe
               zfracal = trb(ji,jj,jk,jpcal) / (trb(ji,jj,jk,jppoc) + trb(ji,jj,jk,jpgoc) + rtrn )
               zgrazcal = (zgrazffeg + zgrazpoc) * (1. - part2) * zfracal !fraction calcite from GOC/POC dissolved in zooplankton guts

               ! calcite production
               zprcaca = xfracal(ji,jj,jk) * zgrazn
               prodcal(ji,jj,jk) = prodcal(ji,jj,jk) + zprcaca  ! prodcal=prodcal(nanophy)+prodcal(microzoo)+prodcal(mesozoo)

               ! zgrazcal [GOC/POC] --> [MES] --> [DIC]
               ! zprcaca [PHY] --> [MES] --> [CAL]
               zprcaca = part2 * zprcaca
               tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zgrazcal - zprcaca
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2. * ( zgrazcal + zprcaca )
               tra(ji,jj,jk,jpcal) = tra(ji,jj,jk,jpcal) - zgrazcal + zprcaca

               IF ( ln_c13 ) THEN
                  zmortz13 = zmortz * ( (trb(ji,jj,jk,jp13mes)+rtrn) / (trb(ji,jj,jk,jpmes)+rtrn) )
                  tra(ji,jj,jk,jp13mes) = tra(ji,jj,jk,jp13mes) - zmortz13 + zepsherv * zgraztotc13
                  tra(ji,jj,jk,jp13dia) = tra(ji,jj,jk,jp13dia) - zgrazd13
                  tra(ji,jj,jk,jp13zoo) = tra(ji,jj,jk,jp13zoo) - zgrazz13
                  tra(ji,jj,jk,jp13phy) = tra(ji,jj,jk,jp13phy) - zgrazn13
                  tra(ji,jj,jk,jp13poc) = tra(ji,jj,jk,jp13poc) - zgrazpoc13 - zgrazffp13 + zfrac13
                  tra(ji,jj,jk,jp13goc) = tra(ji,jj,jk,jp13goc) + zmortzgoc_13 - zgrazffg13 - zfrac13 + zgrapoc2_13

                  zr13_dic = ( (trb(ji,jj,jk,jp13dic)+rtrn) / (trb(ji,jj,jk,jpdic)+rtrn) )
                  zr13_cal = ( (trb(ji,jj,jk,jp13cal)+rtrn) / (trb(ji,jj,jk,jpcal)+rtrn) )
                  tra(ji,jj,jk,jp13dic) = tra(ji,jj,jk,jp13dic) + zgrazcal * zr13_cal                  &
                  &                       - zprcaca * zr13_dic * (1. - e13c_cal2/1000.)
                  tra(ji,jj,jk,jp13cal) = tra(ji,jj,jk,jp13cal) - zgrazcal * zr13_cal                  &
                  &                       + zprcaca * zr13_dic * (1. - e13c_cal2/1000.)
                  prodcal13(ji,jj,jk) = prodcal13(ji,jj,jk) + zprcaca * zr13_dic * (1. - e13c_cal2/1000.)
               ENDIF
               IF ( ln_n15 ) THEN
                  zmortz15 = zmortz * ( (trb(ji,jj,jk,jp15mes)+rtrn) / (trb(ji,jj,jk,jpmes)+rtrn) )
                  tra(ji,jj,jk,jp15mes) = tra(ji,jj,jk,jp15mes) - zmortz15 + zepsherv * zgraztotc15    &
                  &                       + zgrasigex2_15 * (e15n_ex2/1000.0)                          &
                  &                       + zgrapoc2_15 * (e15n_in2/1000.0)
                  tra(ji,jj,jk,jp15dia) = tra(ji,jj,jk,jp15dia) - zgrazd15
                  tra(ji,jj,jk,jp15zoo) = tra(ji,jj,jk,jp15zoo) - zgrazz15
                  tra(ji,jj,jk,jp15phy) = tra(ji,jj,jk,jp15phy) - zgrazn15
                  tra(ji,jj,jk,jp15poc) = tra(ji,jj,jk,jp15poc) - zgrazpoc15 - zgrazffp15 + zfrac15
                  tra(ji,jj,jk,jp15goc) = tra(ji,jj,jk,jp15goc) + zmortzgoc_15 - zgrazffg15 - zfrac15  &
                  &                       + zgrapoc2_15 * ( 1. - e15n_in2/1000.0 )
                  ! zgrasigex2_15 therefore enriches meso in N15 during excretion of NH4
               ENDIF
            END DO
         END DO
      END DO
      !
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         ALLOCATE( zw3d(jpi,jpj,jpk) )
         IF( iom_use( "GRAZ2" ) ) THEN
            zw3d(:,:,:) = zgrazing2(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)  !   Total grazing of phyto by zooplankton
            CALL iom_put( "GRAZ2", zw3d )
         ENDIF
         IF( iom_use( "EXCR2" ) ) THEN
            zw3d(:,:,:) = excretion2(:,:,:) * rno3 * 1.e+3 * rfact2r * tmask(:,:,:)  ! Total excretion of NH4 by mesozooplankton
            CALL iom_put( "EXCR2", zw3d )
         ENDIF
         IF( iom_use( "PCAL" ) ) THEN
            zw3d(:,:,:) = prodcal(:,:,:) * 1.e+3 * rfact2r * tmask(:,:,:)   !  Calcite production
            CALL iom_put( "PCAL", zw3d )  
         ENDIF
         IF( iom_use( "FEZOO2" ) ) THEN
            zw3d(:,:,:) = zfezoo2(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)   !
            CALL iom_put( "FEZOO2", zw3d )
         ENDIF
         IF( iom_use( "LPRODZ2" ) .AND. ln_ligand )  THEN
            zw3d(:,:,:) = zz2ligprod(:,:,:) * 1e9 * 1.e+3 * rfact2r * tmask(:,:,:)
            CALL iom_put( "LPRODZ2"  , zw3d )
         ENDIF
         DEALLOCATE( zw3d )
      ENDIF
      !
      IF (ln_ligand)  DEALLOCATE( zz2ligprod )
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
        WRITE(charout, FMT="('meso')")
        CALL prt_ctl_trc_info(charout)
        CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_meso')
      !
   END SUBROUTINE p4z_meso


   SUBROUTINE p4z_meso_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_meso_init  ***
      !!
      !! ** Purpose :   Initialization of mesozooplankton parameters
      !!
      !! ** Method  :   Read the nampismes namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampismes
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zmes/ part2, grazrat2, resrat2, mzrat2, xpref2n, xpref2d, xpref2z,   &
         &                xpref2c, xthresh2dia, xthresh2phy, xthresh2zoo, xthresh2poc, &
         &                xthresh2, xkgraz2, epsher2, epsher2min, sigma2, unass2, grazflux, &
         &                e13c_cal2, e15n_ex2, e15n_in2, e18oxy_mes
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*) 
         WRITE(numout,*) 'p4z_meso_init : Initialization of mesozooplankton parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampismes in reference namelist : Pisces mesozooplankton
      READ  ( numnatp_ref, namp4zmes, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zmes in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampismes in configuration namelist : Pisces mesozooplankton
      READ  ( numnatp_cfg, namp4zmes, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zmes in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zmes )
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zmes'
         WRITE(numout,*) '      part of calcite not dissolved in mesozoo guts  part2        =', part2
         WRITE(numout,*) '      mesozoo preference for phyto                   xpref2n      =', xpref2n
         WRITE(numout,*) '      mesozoo preference for diatoms                 xpref2d      =', xpref2d
         WRITE(numout,*) '      mesozoo preference for zoo                     xpref2z      =', xpref2z
         WRITE(numout,*) '      mesozoo preference for poc                     xpref2c      =', xpref2c
         WRITE(numout,*) '      microzoo feeding threshold  for mesozoo        xthresh2zoo  =', xthresh2zoo
         WRITE(numout,*) '      diatoms feeding threshold  for mesozoo         xthresh2dia  =', xthresh2dia
         WRITE(numout,*) '      nanophyto feeding threshold for mesozoo        xthresh2phy  =', xthresh2phy
         WRITE(numout,*) '      poc feeding threshold for mesozoo              xthresh2poc  =', xthresh2poc
         WRITE(numout,*) '      feeding threshold for mesozooplankton          xthresh2     =', xthresh2
         WRITE(numout,*) '      exsudation rate of mesozooplankton             resrat2      =', resrat2
         WRITE(numout,*) '      mesozooplankton mortality rate                 mzrat2       =', mzrat2
         WRITE(numout,*) '      maximal mesozoo grazing rate                   grazrat2     =', grazrat2
         WRITE(numout,*) '      mesozoo flux feeding rate                      grazflux     =', grazflux
         WRITE(numout,*) '      non assimilated fraction of P by mesozoo       unass2       =', unass2
         WRITE(numout,*) '      Efficiency of Mesozoo growth                   epsher2      =', epsher2
         WRITE(numout,*) '      Minimum Efficiency of Mesozoo growth           epsher2min  =', epsher2min
         WRITE(numout,*) '      Fraction of mesozoo excretion as DOM           sigma2       =', sigma2
         WRITE(numout,*) '      half sturation constant for grazing 2          xkgraz2      =', xkgraz2
         WRITE(numout,*) '      C13 mesozoo calcification fractionation        e13c_cal2    =', e13c_cal2
         WRITE(numout,*) '      N15 mesozoo excretion fractionation            e15n_ex2     =', e15n_ex2
         WRITE(numout,*) '      N15 mesozoo ingestion fractionation            e15n_in2     =', e15n_in2
         WRITE(numout,*) '      O15 mesozoo respiration fractionation          e18oxy_mes   =', e18oxy_mes
      ENDIF
      !
   END SUBROUTINE p4z_meso_init

   !!======================================================================
END MODULE p4zmeso
