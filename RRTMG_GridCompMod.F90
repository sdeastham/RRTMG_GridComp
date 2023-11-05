#include "MAPL_Generic.h"
!-------------------------------------------------------------------------
!         GEOS-Chem High Performance Global Chemical Transport Model
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: RRTMG_GridCompMod
!
! !DESCRIPTION:Radiative transfer component.
!\\
!\\
! !INTERFACE:
!
module RRTMG_GridComp
!
! !USES:
!
   use ESMF
   use MAPL_Mod
   use pFlogger, only: logging, Logger
   USE PARRRTM,      ONLY : NBNDLW
   USE PARRRSW,      ONLY : NBNDSW
   
   implicit none
   private
!
! !PUBLIC MEMBER FUNCTIONS:
!
   public  :: SetServices
!
! !PRIVATE MEMBER FUNCTIONS:
!
   private :: Initialize
   private :: Run
   private :: Finalize
!
! !PUBLIC DATA MEMBERS:
!
   ! NOTE: Changed to pointers to get inputs from HEMCO (bmy, 10/30/18)
   ! NOTE: These should eventually go into fields of State_Chm
   REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDFNIR(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDFVIS(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDRNIR(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_ALBDRVIS(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_01(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_02(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_03(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_04(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_05(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_06(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_07(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_08(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_09(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_10(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_11(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_12(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_13(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_14(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_15(:,:  )
   REAL*4,  POINTER,     PUBLIC         :: MODIS_EMISS_16(:,:  )

   !MCICA cloud variables now stored for reuse
   !NOTE: These should eventually go into fields of State_Chm
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLDFMCL_LW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CIWPMCL_LW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLWPMCL_LW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: TAUCMCL_LW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLDFMCL_SW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CIWPMCL_SW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: CLWPMCL_SW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: TAUCMCL_SW(:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: SSACMCL   (:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: ASMCMCL   (:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: FSFCMCL   (:,:,:,:)
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: REICMCL   (:,:,:  )
   REAL*8,  ALLOCATABLE, PUBLIC, TARGET :: RELQMCL   (:,:,:  )
!
! !PRIVATE DATA MEMBERS:
!
   integer, parameter :: r8     = 8
   integer, parameter :: r4     = 4
   
!   logical :: meteorology_vertical_index_is_top_down
!   integer :: use_total_air_pressure_in_advection
!   integer :: correct_mass_flux_for_humidity
   
   class(Logger), pointer :: lgr => null()
!
! !REVISION HISTORY:
!  18 Jun 2013 - D.A. Ridley - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL*8,  ALLOCATABLE  :: LW_UFLUX (:,:,:)
   REAL*8,  ALLOCATABLE  :: LW_DFLUX (:,:,:)
   REAL*8,  ALLOCATABLE  :: SW_UFLUX (:,:,:)
   REAL*8,  ALLOCATABLE  :: SW_DFLUX (:,:,:)
   REAL*8,  ALLOCATABLE  :: LW_UFLUXC(:,:,:)
   REAL*8,  ALLOCATABLE  :: LW_DFLUXC(:,:,:)
   REAL*8,  ALLOCATABLE  :: SW_UFLUXC(:,:,:)
   REAL*8,  ALLOCATABLE  :: SW_DFLUXC(:,:,:)
 
   REAL*8  :: RRTMG_LMB(NBNDLW+NBNDSW)
 
   INTEGER :: ID_AER_LMB0 (NBNDLW+NBNDSW)
   INTEGER :: ID_AER_LMB1 (NBNDLW+NBNDSW)
!
! !REMARKS:
!  This file was adapted from a GEOS file developed at NASA GMAO.
!                                                                             .
!  NOTES:
!  - The abbreviation "PET" stands for "Persistent Execution Thread".
!    It is a synomym for CPU.
!
! !REVISION HISTORY:
!  06 Dec 2009 - A. da Silva - Initial version this file was adapted from
!
!EOP
!------------------------------------------------------------------------------
!BOC
   contains
!EOC

!-------------------------------------------------------------------------
!         GEOS-Chem High Performance Global Chemical Transport Model
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices -- Sets ESMF services for this component
!
! !INTERFACE:
!
   subroutine SetServices(GC, RC)
!
! !INPUT/OUTPUT PARAMETERS
!
      type(ESMF_GridComp), intent(INOUT) :: GC  ! gridded component
!
! !OUTPUT PARAMETERS
!
      integer, intent(OUT)               :: RC  ! return code
!
! !DESCRIPTION:
!   The SetServices needs to register its Initialize, Run, and Finalize.
!   It uses the MAPL_Generic construct for defining state specs.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer           :: STATUS, I
      type(ESMF_Config) :: CF

      character(len=ESMF_MAXSTR) :: COMP_NAME, msg, short_name, long_name
      character(len=ESMF_MAXSTR) :: IAm = 'SetServices'

      !================================
      ! SetServices starts here
      !================================
      
      ! Get gridded component name and set-up traceback handle
      ! -----------------------------------------------------------------
      call ESMF_GridCompGet(GC, NAME=COMP_NAME, CONFIG=CF, RC=STATUS)
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)
      lgr => logging%get_logger('RRTMG')
      
      ! Get whether to import mass fluxes from ExtData or derive from winds
      ! -----------------------------------------------------------------
!      call ESMF_ConfigGetAttribute(CF,                                     &
!                                   value=import_mass_flux_from_extdata,    &
!                                   label='IMPORT_MASS_FLUX_FROM_EXTDATA:', &
!                                   Default=.false.,                        &
!                                   __RC__)
!      if (import_mass_flux_from_extdata) then
!         msg = 'Configured to import mass fluxes from ''ExtData'''
!      else
!         msg = 'Configured to derive and export mass flux and courant numbers'
!      end if
!      call lgr%info(msg)
      
      ! Register services for this component
      ! -----------------------------------------------------------------
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, &
                                      RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_FINALIZE, Finalize, RC=STATUS)
      _VERIFY(STATUS)
      
      ! Define Import state
      ! -----------------------------------------------------------------
      call lgr%debug('Adding import specs')
      call MAPL_AddImportSpec(gc, &
                              SHORT_NAME='MODIS_ALBDFNIR', &
                              LONG_NAME='diffuse_near-ir_albedo',&
                              UNITS='-', &
                              DIMS=MAPL_DimsHorzOnly, &
                              VLOCATION=MAPL_VLocationEdge, &
                              RC=STATUS)
      _VERIFY(STATUS)
      
      call MAPL_AddImportSpec(gc, &
                              SHORT_NAME='MODIS_ALBDFVIS', &
                              LONG_NAME='diffuse_visible_albedo',&
                              UNITS='-', &
                              DIMS=MAPL_DimsHorzOnly, &
                              VLOCATION=MAPL_VLocationEdge, &
                              RC=STATUS)
      _VERIFY(STATUS)
      
      call MAPL_AddImportSpec(gc, &
                              SHORT_NAME='MODIS_ALBDRNIR', &
                              LONG_NAME='direct_near-ir_albedo',&
                              UNITS='-', &
                              DIMS=MAPL_DimsHorzOnly, &
                              VLOCATION=MAPL_VLocationEdge, &
                              RC=STATUS)
      _VERIFY(STATUS)
      
      call MAPL_AddImportSpec(gc, &
                              SHORT_NAME='MODIS_ALBDRVIS', &
                              LONG_NAME='direct_visible_albedo',&
                              UNITS='-', &
                              DIMS=MAPL_DimsHorzOnly, &
                              VLOCATION=MAPL_VLocationEdge, &
                              RC=STATUS)
      _VERIFY(STATUS)
     
      Do I=1, 16
          write(short_name,'(a,I0.2)') 'MODIS_EMISSIVITY_', I
          write(long_name,'(a,I0.2)') 'emissivity_in_modis_band_', I
          call MAPL_AddImportSpec(gc, &
                                  SHORT_NAME=trim(short_name), &
                                  LONG_NAME=trim(long_name),&
                                  UNITS='-', &
                                  DIMS=MAPL_DimsHorzOnly, &
                                  VLOCATION=MAPL_VLocationEdge, &
                                  RC=STATUS)
          _VERIFY(STATUS)
      End Do 

      ! Define Export State
      ! -----------------------------------------------------------------
      call lgr%debug('Adding export specs')
!      call MAPL_AddExportSpec(gc, &
!                              SHORT_NAME='SPHU0', &
!                              LONG_NAME='specific_humidity_before_advection', &
!                              UNITS='kg kg-1', &
!                              PRECISION=ESMF_KIND_R8, &
!                              DIMS=MAPL_DimsHorzVert, &
!                              VLOCATION=MAPL_VLocationCenter, &
!                              RC=STATUS)
      _VERIFY(STATUS)

      ! Set profiling timers
      !-------------------------
      call lgr%debug('Adding timers')
      call MAPL_TimerAdd(gc, name="INITIALIZE", RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_TimerAdd(gc, name="RUN", RC=STATUS)
      _VERIFY(STATUS)
      call MAPL_TimerAdd(gc, name="FINALIZE", RC=STATUS)
      _VERIFY(STATUS)
      
      call lgr%debug('Calling MAPL_GenericSetServices')

      ! Create children's gridded components and invoke their SetServices
      ! -----------------------------------------------------------------
      call MAPL_GenericSetServices(gc, RC=STATUS)
      _VERIFY(STATUS)
      
      _RETURN(ESMF_SUCCESS)
      
   end subroutine SetServices
!EOC
!-------------------------------------------------------------------------
!         GEOS-Chem High Performance Global Chemical Transport Model
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialize -- Initialized method for composite the CTMder
!
! !INTERFACE:
!
   subroutine Initialize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
!
      integer, optional,   intent(out)   :: RC     ! Error code
!
! !DESCRIPTION:
!  The Initialize method of the radiative transfer Gridded Component.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer                    :: comm
      character(len=ESMF_MAXSTR) :: COMP_NAME
      character(len=ESMF_MAXSTR) :: msg
      type(ESMF_Config)          :: CF
      type(ESMF_Grid)            :: esmfGrid
      type(ESMF_VM)              :: VM

      type(MAPL_MetaComp), pointer  :: ggState ! MAPL Generic State
      REAL, POINTER, DIMENSION(:,:) :: cellArea

      !================================
      ! Initialize starts here
      !================================

      __Iam__('Initialize')
      
      !  Get this gridded component name and set-up traceback handle
      ! -----------------------------------------------------------------
      call ESMF_GridCompGet(GC, NAME=COMP_NAME, CONFIG=CF, VM=VM, RC=STATUS)
      _VERIFY(STATUS)
      Iam = TRIM(COMP_NAME)//"::Initialize"
      
      !  Initialize MAPL_Generic
      ! -----------------------------------------------------------------
      call MAPL_GenericInitialize(gc, IMPORT, EXPORT, clock, RC=STATUS)
      _VERIFY(STATUS)
      
      !  Get internal MAPL_Generic state
      ! -----------------------------------------------------------------
      call MAPL_GetObjectFromGC(GC, ggState, RC=STATUS)
      _VERIFY(STATUS)

      ! Turn on timers
      ! -----------------------------------------------------------------
      call MAPL_TimerOn(ggSTATE, "TOTAL")
      call MAPL_TimerOn(ggSTATE, "INITIALIZE")
      
      ! Get grid-related information
      ! -----------------------------------------------------------------
      call ESMF_GridCompGet(GC, GRID=esmfGrid, rc=STATUS)
      _VERIFY(STATUS)

      ! Turn off timers
      ! -----------------------------------------------------------------
      call MAPL_TimerOff(ggSTATE,"INITIALIZE")
      call MAPL_TimerOff(ggSTATE,"TOTAL")
      
      _RETURN(ESMF_SUCCESS)
      
   end subroutine Initialize
!EOC
!-------------------------------------------------------------------------
!         GEOS-Chem High Performance Global Chemical Transport Model
!-------------------------------------------------------------------------
!BOP
!
! !INTERFACE:
!
   subroutine Run(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
!
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
!
      integer, optional,   intent(out) :: RC       ! Error code
!
! !DESCRIPTION:
! The Run method of the radiative transfer Gridded Component.
!
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer                      :: STATUS
      integer                      :: ndt
      character(len=ESMF_MAXSTR)   :: IAm = "Run"
      character(len=ESMF_MAXSTR)   :: COMP_NAME
      type(MAPL_MetaComp), pointer :: ggState
      type(ESMF_Grid)              :: esmfGrid
      real(r8)                     :: dt
      real(r8), pointer            :: PLE(:,:,:) ! Edge pressures

      ! Saved variables
      logical, save :: firstRun = .true.
      
#ifdef ADJOINT
      integer :: reverseTime
#endif

      !================================
      ! Run starts here
      !================================
      
      ! Get this component's name and set-up traceback handle.
      call ESMF_GridCompGet(GC, name=COMP_NAME, Grid=esmfGrid, RC=STATUS)
      _VERIFY(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)
      
      ! Get internal MAPL_Generic state
      call MAPL_GetObjectFromGC(GC, ggState, RC=STATUS)
      _VERIFY(STATUS)

      ! Turn on timers
      call MAPL_TimerOn(ggState,"TOTAL")
      call MAPL_TimerOn(ggState,"RUN")
      
      ! Retrieve timestep [s] and store as real
      call MAPL_GetResource( ggState,   &
                             ndt,       &
                             'RUN_DT:', &
                             default=0, &
                             RC=STATUS )
      _VERIFY(STATUS)
      dt = ndt
      call lgr%debug('RRTMG run beginning')
      
      ! Turn off timers
      call MAPL_TimerOff(ggState,"RUN")
      call MAPL_TimerOff(ggState,"TOTAL")
      
      _RETURN(ESMF_SUCCESS)

   end subroutine Run
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize - user supplied finalize routine
!
! !INTERFACE:
!
  subroutine Finalize(GC, IMPORT, EXPORT, CLOCK, RC)
!
! !INPUT/OUTPUT PARAMETERS:
      type(ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
      type(ESMF_State),    intent(inout) :: IMPORT ! Import state
      type(ESMF_State),    intent(inout) :: EXPORT ! Export state
      type(ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
!
! !OUTPUT PARAMETERS:
      integer, optional,   intent(  out) :: RC     ! Error code
!
! !DESCRIPTION:
!    Finalize merely destroys the FVadv object that was created in Initialize
!    and releases the space for the persistent data .
!
!EOP
!=============================================================================
!BOC
! !LOCAL VARIABLES:

      character(len=ESMF_MAXSTR)    :: IAm
      integer                       :: STATUS
      character(len=ESMF_MAXSTR)    :: COMP_NAME

      ! Get my name and set-up traceback handle
      ! ---------------------------------------

      Iam = 'Finalize'
      call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
      VERIFY_(STATUS)
      Iam = trim(COMP_NAME) // TRIM(Iam)
      call lgr%debug('RRTMG finalizing')

      !=================================================================
      ! Nullify pointers
      !=================================================================

      ! Albedoes
      MODIS_ALBDFNIR => NULL()
      MODIS_ALBDFVIS => NULL()
      MODIS_ALBDRNIR => NULL()
      MODIS_ALBDRVIS => NULL()

      ! Emissivity
      MODIS_EMISS_01 => NULL()
      MODIS_EMISS_02 => NULL()
      MODIS_EMISS_03 => NULL()
      MODIS_EMISS_04 => NULL()
      MODIS_EMISS_05 => NULL()
      MODIS_EMISS_06 => NULL()
      MODIS_EMISS_07 => NULL()
      MODIS_EMISS_08 => NULL()
      MODIS_EMISS_09 => NULL()
      MODIS_EMISS_10 => NULL()
      MODIS_EMISS_11 => NULL()
      MODIS_EMISS_12 => NULL()
      MODIS_EMISS_13 => NULL()
      MODIS_EMISS_14 => NULL()
      MODIS_EMISS_15 => NULL()
      MODIS_EMISS_16 => NULL()

      !=================================================================
      ! Deallocate surface radiation arrays
      !=================================================================
      IF ( ALLOCATED( LW_UFLUX        ) ) DEALLOCATE( LW_UFLUX        )
      IF ( ALLOCATED( LW_DFLUX        ) ) DEALLOCATE( LW_DFLUX        )
      IF ( ALLOCATED( SW_UFLUX        ) ) DEALLOCATE( SW_UFLUX        )
      IF ( ALLOCATED( SW_DFLUX        ) ) DEALLOCATE( SW_DFLUX        )
      IF ( ALLOCATED( LW_UFLUXC       ) ) DEALLOCATE( LW_UFLUXC       )
      IF ( ALLOCATED( LW_DFLUXC       ) ) DEALLOCATE( LW_DFLUXC       )
      IF ( ALLOCATED( SW_UFLUXC       ) ) DEALLOCATE( SW_UFLUXC       )
      IF ( ALLOCATED( SW_DFLUXC       ) ) DEALLOCATE( SW_DFLUXC       )

      !=================================================================
      ! Deallocate MCICA cloud arrays
      !=================================================================
      IF ( ALLOCATED( CLDFMCL_LW     ) ) DEALLOCATE( CLDFMCL_LW       )
      IF ( ALLOCATED( CIWPMCL_LW     ) ) DEALLOCATE( CIWPMCL_LW       )
      IF ( ALLOCATED( CLWPMCL_LW     ) ) DEALLOCATE( CLWPMCL_LW       )
      IF ( ALLOCATED( TAUCMCL_LW     ) ) DEALLOCATE( TAUCMCL_LW       )
      IF ( ALLOCATED( CLDFMCL_SW     ) ) DEALLOCATE( CLDFMCL_SW       )
      IF ( ALLOCATED( CIWPMCL_SW     ) ) DEALLOCATE( CIWPMCL_SW       )
      IF ( ALLOCATED( CLWPMCL_SW     ) ) DEALLOCATE( CLWPMCL_SW       )
      IF ( ALLOCATED( TAUCMCL_SW     ) ) DEALLOCATE( TAUCMCL_SW       )
      IF ( ALLOCATED( SSACMCL        ) ) DEALLOCATE( SSACMCL          )
      IF ( ALLOCATED( ASMCMCL        ) ) DEALLOCATE( ASMCMCL          )
      IF ( ALLOCATED( FSFCMCL        ) ) DEALLOCATE( FSFCMCL          )
      IF ( ALLOCATED( REICMCL        ) ) DEALLOCATE( REICMCL          )
      IF ( ALLOCATED( RELQMCL        ) ) DEALLOCATE( RELQMCL          )

      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, RC)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
      end subroutine Finalize
   
end module RRTMG_GridComp
