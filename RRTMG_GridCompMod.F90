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
!   logical, public :: import_mass_flux_from_extdata = .false.
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
      integer           :: STATUS
      type(ESMF_Config) :: CF

      character(len=ESMF_MAXSTR) :: COMP_NAME, msg
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
!      call MAPL_AddImportSpec(gc, &
!                              SHORT_NAME='PS1', &
!                              LONG_NAME='pressure_at_surface_before_advection',&
!                              UNITS='hPa', &
!                              DIMS=MAPL_DimsHorzOnly, &
!                              VLOCATION=MAPL_VLocationEdge, &
!                              RC=STATUS)
!      _VERIFY(STATUS)
      
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

      !if (.NOT. FV3_DynCoreIsRunning) then
      !   call fv_end(FV_Atm, grids_on_my_pe, .false.)
      !endif

      call MAPL_GenericFinalize(GC, IMPORT, EXPORT, CLOCK, RC)
      VERIFY_(STATUS)

      RETURN_(ESMF_SUCCESS)
      end subroutine Finalize
   
end module RRTMG_GridComp
