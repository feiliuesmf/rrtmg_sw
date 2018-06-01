!--------------- RRTMG NUOPC CAP -----------------
! This is the RRTMG model cap component that's NUOPC compiant.
!
! Author:  Fei.Liu@NOAA.GOV
!
! 3/22/18
! This is now acting as a cap/connector between NUOPC driver and RRTMG code.
! As a radiation cap, Fields are created on the Grid transferred over from ATM.
!

#define USE_MESH

module rrtmg_cap

  use rrtmg_sw_init
  use rrtmg_sw_rad

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS         => SetServices, &
    model_label_CheckImport  => label_CheckImport, &
    model_label_Advance      => label_Advance, &
    model_label_DataInitialize  => label_DataInitialize, &
    model_label_Finalize     => label_Finalize

  implicit none
  private
  public SetServices

  character(len=45)     :: exportFieldList(6) = (/ &
    "Total Sky Shortward Downard Flux            ", &
    "Total Sky Shortward Radiative Heating Rate  ", &
    "Total Sky Shortwave Upward Flux             ", &
    "Clear Sky Shortwave Downward Flux           ", &
    "Clear Sky Shortwave Radiative Heating Rate  ", &
    "Clear Sky Shortwave Upward Flux             " /)
  character(len=8)      :: exportFieldSN(6) = (/ &
    "tssdf   ", &
    "tssrhr  ", &
    "tssuf   ", &
    "cssdf   ", &
    "cssrhr  ", &
    "cssuf   " /)

  character(len=55)      :: importFieldList(39) = (/ &
  "Layer pressures (hPa, mb)                             ",&
  "Interface pressures (hPa, mb)                         ",&
  "Layer temperatures (K)                                ",&
  "Interface temperatures (K)                            ",&
  "Surface temperature (K)                               ",&
  "H2O volume mixing ratio                               ",&
  "O3 volume mixing ratio                                ",&
  "CO2 volume mixing ratio                               ",&
  "Methane volume mixing ratio                           ",&
  "Nitrous oxide volume mixing ratio                     ",&
  "Oxygen volume mixing ratio                            ",&
  "UV/vis surface albedo direct rad                      ",&
  "Near-IR surface albedo direct rad                     ",&
  "UV/vis surface albedo: diffuse rad                    ",&
  "Near-IR surface albedo: diffuse rad                   ",&
  "Day of the year (used to get Earth/Sun                ",&
  "Flux adjustment for Earth/Sun distance                ",&
  "Cosine of solar zenith angle                          ",&
  "Solar constant (W/m2)                                 ",&
  "Flag for solar variability method                     ",&
  "Facular and sunspot amplitude                         ",&
  "Solar variability scale factors                       ",&
  "Fraction of averaged 11-year solar cycle (0-1)        ",&
  "Flag for cloud optical properties                     ",&
  "Flag for ice particle specification                   ",&
  "Flag for liquid droplet specification                 ",&
  "Cloud fraction                                        ",&
  "In-cloud optical depth                                ",&
  "In-cloud single scattering albedo                     ",&
  "In-cloud asymmetry parameter                          ",&
  "In-cloud forward scattering fraction                  ",&
  "In-cloud ice water path (g/m2)                        ",&
  "In-cloud liquid water path (g/m2)                     ",&
  "Cloud ice effective radius (microns)                  ",&
  "Cloud water drop effective radius (microns)           ",&
  "Aerosol optical depth (iaer=10 only)                  ",&
  "Aerosol single scattering albedo (iaer=10 only)       ",&
  "Aerosol asymmetry parameter (iaer=10 only)            ",&
  "Aerosol optical depth at 0.55 micron (iaer=6 only)    " /)
  character(len=10)       :: importFieldSN(39) = (/ &
  "play      ", &
  "plev      ", &
  "tlay      ", &
  "tlev      ", &
  "tsfc      ", &
  "h2ovmr    ", &
  "o3vmr     ", &
  "co2vmr    ", &
  "ch4vmr    ", &
  "n2ovmr    ", &
  "o2vmr     ", &
  "asdir     ", &
  "aldir     ", &
  "asdif     ", &
  "aldif     ", &
  "dyofyr    ", &
  "adjes     ", &
  "coszen    ", &
  "scon      ", &
  "isolvar   ", &
  "indsolvar ", &
  "bndsolvar ", &
  "solcycfrac", &
  "inflgsw   ", &
  "iceflgsw  ", &
  "liqflgsw  ", &
  "cldfmcl   ", &
  "taucmcl   ", &
  "ssacmcl   ", &
  "asmcmcl   ", &
  "fsfcmcl   ", &
  "ciwpmcl   ", &
  "clwpmcl   ", &
  "reicmcl   ", &
  "relqmcl   ", &
  "tauaer    ", &
  "ssaaer    ", &
  "asmaer    ", &
  "ecaer     " /)

  integer   :: import_slice = 0
  integer   :: export_slice = 0

  type FieldListType
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: transferOfferGeom
    character(len=64) :: transferOfferField
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:),   pointer :: farrayPtr2D => null()
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr => null()
  end type FieldListType

  integer :: fldsToRRTMG_num = 0
  integer :: fldsFrRRTMG_num = 6

  !https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
  real(kind=rb)          :: cpdair =  1004      ! Cp = 1004 J/kg K at 298K

  integer(kind=im)       :: ncol = 1        ! Number of horizontal columns     
  integer(kind=im)       :: nlay = 50       ! Number of model layers
  !integer(kind=im), parameter :: ngptsw = 112

  integer(kind=im)       :: icld            ! Cloud overlap method
  integer(kind=im)       :: iaer            ! Aerosol option flag
  real                   :: p0 = 1013.25    ! 1 atm
  real(kind=rb), pointer :: play(:,:)          ! Layer pressures (hPa, mb)
  real(kind=rb), pointer :: plev(:,:)          ! Interface pressures (hPa, mb)
  real(kind=rb), pointer :: tlay(:,:)          ! Layer temperatures (K)
  real(kind=rb), pointer :: tlev(:,:)          ! Interface temperatures (K)
  real(kind=rb), pointer :: tsfc(:)            ! Surface temperature (K)
  real(kind=rb), pointer :: h2ovmr(:,:)        ! H2O volume mixing ratio
  real(kind=rb), pointer :: o3vmr(:,:)         ! O3 volume mixing ratio
  real(kind=rb), pointer :: co2vmr(:,:)        ! CO2 volume mixing ratio
  real(kind=rb), pointer :: ch4vmr(:,:)        ! Methane volume mixing ratio
  real(kind=rb), pointer :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
  real(kind=rb), pointer :: o2vmr(:,:)         ! Oxygen volume mixing ratio
  real(kind=rb), pointer :: asdir(:)           ! UV/vis surface albedo direct rad
  real(kind=rb), pointer :: aldir(:)           ! Near-IR surface albedo direct rad
  real(kind=rb), pointer :: asdif(:)           ! UV/vis surface albedo: diffuse rad
  real(kind=rb), pointer :: aldif(:)           ! Near-IR surface albedo: diffuse rad
  integer(kind=im)       :: dyofyr             ! Day of the year (used to get Earth/Sun
  real(kind=rb)          :: adjes              ! Flux adjustment for Earth/Sun distance
  real(kind=rb), pointer :: coszen(:)          ! Cosine of solar zenith angle
  real(kind=rb)          :: scon               ! Solar constant (W/m2)
  integer(kind=im)       :: isolvar            ! Flag for solar variability method
  real(kind=rb), pointer :: indsolvar(:)       ! Facular and sunspot amplitude 
  real(kind=rb), pointer :: bndsolvar(:)       ! Solar variability scale factors 
  real(kind=rb)          :: solcycfrac         ! Fraction of averaged 11-year solar cycle (0-1)
  integer(kind=im)       :: inflgsw            ! Flag for cloud optical properties
  integer(kind=im)       :: iceflgsw           ! Flag for ice particle specification
  integer(kind=im)       :: liqflgsw           ! Flag for liquid droplet specification
  real(kind=rb), pointer :: cldfmcl(:,:,:)     ! Cloud fraction
  real(kind=rb), pointer :: taucmcl(:,:,:)     ! In-cloud optical depth
  real(kind=rb), pointer :: ssacmcl(:,:,:)     ! In-cloud single scattering albedo
  real(kind=rb), pointer :: asmcmcl(:,:,:)     ! In-cloud asymmetry parameter
  real(kind=rb), pointer :: fsfcmcl(:,:,:)     ! In-cloud forward scattering fraction
  real(kind=rb), pointer :: ciwpmcl(:,:,:)     ! In-cloud ice water path (g/m2)
  real(kind=rb), pointer :: clwpmcl(:,:,:)     ! In-cloud liquid water path (g/m2)
  real(kind=rb), pointer :: reicmcl(:,:)       ! Cloud ice effective radius (microns)
  real(kind=rb), pointer :: relqmcl(:,:)       ! Cloud water drop effective radius (microns)
  real(kind=rb), pointer :: tauaer(:,:,:)      ! Aerosol optical depth (iaer=10 only)
  real(kind=rb), pointer :: ssaaer(:,:,:)      ! Aerosol single scattering albedo (iaer=10 only)
  real(kind=rb), pointer :: asmaer(:,:,:)      ! Aerosol asymmetry parameter (iaer=10 only)
  real(kind=rb), pointer :: ecaer(:,:,:)       ! Aerosol optical depth at 0.55 micron (iaer=6 only)

! --!--- Output -----

  real(kind=rb), pointer :: swuflx(:,:)       ! Total sky shortwave upward flux (W/m2)
                                                  !    Dimensions: (ncol,nlay+1)
  real(kind=rb), pointer :: swdflx(:,:)       ! Total sky shortwave downward flux (W/m2)
                                                  !    Dimensions: (ncol,nlay+1)
  real(kind=rb), pointer :: swhr(:,:)         ! Total sky shortwave radiative heating rate (K/d)
                                                  !    Dimensions: (ncol,nlay)
  real(kind=rb), pointer :: swuflxc(:,:)      ! Clear sky shortwave upward flux (W/m2)
                                                  !    Dimensions: (ncol,nlay+1)
  real(kind=rb), pointer :: swdflxc(:,:)      ! Clear sky shortwave downward flux (W/m2)
                                                  !    Dimensions: (ncol,nlay+1)
  real(kind=rb), pointer :: swhrc(:,:)        ! Clear sky shortwave radiative heating rate (K/d)
                                                  !    Dimensions: (ncol,nlay)

  contains
  !-----------------------------------------------------------------------
  !------------------- RRTMG code starts here -----------------------
  !-----------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(gcomp, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! switching to IPD versions
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p3"/), userRoutine=InitializeRealize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p4"/), userRoutine=InitializeAcceptGrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p5"/), userRoutine=InitializeCompleteField, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    call ESMF_MethodAdd(gcomp, label=model_label_Advance, &
      userRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_Finalize, &
      specRoutine=rrtmg_model_finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_MethodRemove(gcomp, model_label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_CheckImport, &
      specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(gcomp, specLabel=model_label_DataInitialize, &
      specRoutine=DataInitialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: model
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    rc = ESMF_SUCCESS

    ! Switch to IPDv03 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(model, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeAdvertise(gcomp, importState, exportState, clock, rc)

    type(ESMF_GridComp)                    :: gcomp
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Clock)                       :: clock
    integer, intent(out)                   :: rc

    ! Local Variables
    type(ESMF_VM)                          :: vm
    integer                                :: mpi_comm

    rc = ESMF_SUCCESS

    call ESMF_VMGetCurrent(vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call RRTMG_FieldsSetup(rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !call AdvertiseImpFields(importState, fldsToRRTMG_num, fldsToRRTMG, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    call AdvertiseExpFields(exportState, fldsFrRRTMG_num, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    contains

    subroutine AdvertiseExpFields(state, nfields, rc)

      type(ESMF_State), intent(inout)             :: state
      integer,intent(in)                          :: nfields
      integer, intent(inout)                      :: rc

      integer                                     :: i

      rc = ESMF_SUCCESS

      do i = 1, nfields

        call NUOPC_Advertise(state, &
          standardName=exportFieldList(i), &
          !name=exportFieldSN(i), &
          TransferOfferGeomObject="cannot provide", &
          rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

      enddo

    end subroutine AdvertiseExpFields

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(gcomp, rc)
    type(ESMF_GridComp)                    :: gcomp
    integer, intent(out)                   :: rc
    
    ! local variables
    integer                                :: i, j 
    type(ESMF_Clock)                       :: clock
    type(ESMF_State)                       :: importState, exportState
    type(ESMF_Time)                        :: currTime
    type(ESMF_TimeInterval)                :: timeStep

    rc = ESMF_SUCCESS

    import_slice = import_slice + 1
    export_slice = export_slice + 1

    !  integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
    !                                                  !    0: Clear only
    !                                                  !    1: Random
    !                                                  !    2: Maximum/random
    !                                                  !    3: Maximum
    icld = 0
    !  integer(kind=im), intent(inout) :: iaer         ! Aerosol option flag
    !                                                  !    0: No aerosol
    !                                                  !    6: ECMWF method
    !                                                  !    10:Input aerosol optical 
    !                                                  !       properties
    iaer = 0
    !  real(kind=rb), intent(in) :: play(:,:)          ! Layer pressures (hPa, mb)
    !                                                  !    Dimensions: (ncol,nlay)
    p0 = 1013.25
    do i = 1, ncol
      do j = 1, nlay
        play(i,j) = p0*exp(-DBLE(j)/nlay)
      enddo
    enddo
    !  real(kind=rb), intent(in) :: plev(:,:)          ! Interface pressures (hPa, mb)
    !                                                  !    Dimensions: (ncol,nlay+1)
    do i = 1, ncol
      do j = 1, nlay+1
        plev(i,j) = p0*exp(-DBLE(j-1)/nlay)
      enddo
    enddo
   
    !  real(kind=rb), intent(in) :: tlay(:,:)          ! Layer temperatures (K)
    !                                                  !    Dimensions: (ncol,nlay)
    do i = 1, ncol
      do j = 1, nlay
        tlay(i,j) = 298
      enddo
    enddo
    !  real(kind=rb), intent(in) :: tlev(:,:)          ! Interface temperatures (K)
    !                                                  !    Dimensions: (ncol,nlay+1)
    do i = 1, ncol
      do j = 1, nlay+1
        tlev(i,j) = 298
      enddo
    enddo
    !  real(kind=rb), intent(in) :: tsfc(:)            ! Surface temperature (K)
    !                                                  !    Dimensions: (ncol)
    do i = 1, ncol
      tsfc(i) = 298
    enddo
    !  real(kind=rb), intent(in) :: h2ovmr(:,:)        ! H2O volume mixing ratio
    !                                                  !    Dimensions: (ncol,nlay)
    ! https://www.researchgate.net/figure/Models-of-the-vertical-profile-of-the-water-vapor-volume-mixing-ratio-from-Irwin_fig5_254364473
    do i = 1, ncol
      do j = 1, nlay
        h2ovmr(i,j) = 0.001  ! 10^-8 -> 10^-2
      enddo
    enddo
  
    !  real(kind=rb), intent(in) :: o3vmr(:,:)         ! O3 volume mixing ratio
    !                                                  !    Dimensions: (ncol,nlay)
    ! https://www.briangwilliams.us/atmospheric-chemistry/atmospheric-ozone.html
    do i = 1, ncol
      do j = 1, nlay
        o3vmr(i,j) = 1.E-6  ! around 0.01-10 10^-6 ppmv
      enddo
    enddo

    !  real(kind=rb), intent(in) :: co2vmr(:,:)        ! CO2 volume mixing ratio
    !                                                  !    Dimensions: (ncol,nlay)
    ! http://slideplayer.com/slide/5092962/
    do i = 1, ncol
      do j = 1, nlay
        co2vmr(i,j) = 380*1E-6 
      enddo
    enddo
    !  real(kind=rb), intent(in) :: ch4vmr(:,:)        ! Methane volume mixing ratio
    !                                                  !    Dimensions: (ncol,nlay)
    do i = 1, ncol
      do j = 1, nlay
        ch4vmr(i,j) = 1.7*1E-6 
      enddo
    enddo
    !  real(kind=rb), intent(in) :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
    !                                                  !    Dimensions: (ncol,nlay)
    ! https://www.researchgate.net/figure/Volume-mixing-ratio-for-selected-gaseous-species_tbl6_233792155
    do i = 1, ncol
      do j = 1, nlay
        n2ovmr(i,j) = 0.318*1E-6 
      enddo
    enddo
    !  real(kind=rb), intent(in) :: o2vmr(:,:)         ! Oxygen volume mixing ratio
    !                                                  !    Dimensions: (ncol,nlay)
    do i = 1, ncol
      do j = 1, nlay
        o2vmr(i,j) = 0.21
      enddo
    enddo
    !  real(kind=rb), intent(in) :: asdir(:)           ! UV/vis surface albedo direct rad
    !                                                  !    Dimensions: (ncol)
    do i = 1, ncol
      asdir(i) = 0.21
    enddo
    !  real(kind=rb), intent(in) :: aldir(:)           ! Near-IR surface albedo direct rad
    !                                                  !    Dimensions: (ncol)
    do i = 1, ncol
      aldir(i) = 0.21
    enddo
    !  real(kind=rb), intent(in) :: asdif(:)           ! UV/vis surface albedo: diffuse rad
    !                                                  !    Dimensions: (ncol)
    ! 
    do i = 1, ncol
      asdif(i) = 0.21
    enddo
    !  real(kind=rb), intent(in) :: aldif(:)           ! Near-IR surface albedo: diffuse rad
    !                                                  !    Dimensions: (ncol)
    do i = 1, ncol
      aldif(i) = 0.21
    enddo

    !  integer(kind=im), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
    !                                                  !  distance if adjflx not provided)
    dyofyr = 180
    !  real(kind=rb), intent(in) :: adjes              ! Flux adjustment for Earth/Sun distance
    adjes = 1.
    !  real(kind=rb), intent(in) :: coszen(:)          ! Cosine of solar zenith angle
    !                                                  !    Dimensions: (ncol)
    do i = 1, ncol
      coszen(i) = 0.8
    enddo
    !  real(kind=rb), intent(in) :: scon               ! Solar constant (W/m2)
    !                                                  !    Total solar irradiance averaged 
    !                                                  !    over the solar cycle.
    !                                                  !    If scon = 0.0, the internal solar 
    !                                                  !    constant, which depends on the  
    !                                                  !    value of isolvar, will be used. 
    !                                                  !    For isolvar=-1, scon=1368.22 Wm-2,
    !                                                  !    For isolvar=0,1,3, scon=1360.85 Wm-2,
    !                                                  !    If scon > 0.0, the internal solar
    !                                                  !    constant will be scaled to the 
    !                                                  !    provided value of scon.
    scon = 0.0

    !  integer(kind=im), intent(in) :: isolvar         ! Flag for solar variability method
    !                                                  !   -1 = (when scon .eq. 0.0): No solar variability
    !                                                  !        and no solar cycle (Kurucz solar irradiance
    !                                                  !        of 1368.22 Wm-2 only);
    !                                                  !        (when scon .ne. 0.0): Kurucz solar irradiance
    !                                                  !        scaled to scon and solar variability defined
    !                                                  !        (optional) by setting non-zero scale factors
    !                                                  !        for each band in bndsolvar
    !                                                  !    0 = (when SCON .eq. 0.0): No solar variability 
    !                                                  !        and no solar cycle (NRLSSI2 solar constant of 
    !                                                  !        1360.85 Wm-2 for the 100-50000 cm-1 spectral 
    !                                                  !        range only), with facular and sunspot effects 
    !                                                  !        fixed to the mean of Solar Cycles 13-24;
    !                                                  !        (when SCON .ne. 0.0): No solar variability 
    !                                                  !        and no solar cycle (NRLSSI2 solar constant of 
    !                                                  !        1360.85 Wm-2 for the 100-50000 cm-1 spectral 
    !                                                  !        range only), is scaled to SCON
    !                                                  !    1 = Solar variability (using NRLSSI2  solar
    !                                                  !        model) with solar cycle contribution
    !                                                  !        determined by fraction of solar cycle
    !                                                  !        with facular and sunspot variations
    !                                                  !        fixed to their mean variations over the
    !                                                  !        average of Solar Cycles 13-24;
    !                                                  !        two amplitude scale factors allow
    !                                                  !        facular and sunspot adjustments from
    !                                                  !        mean solar cycle as defined by indsolvar 
    !                                                  !    2 = Solar variability (using NRLSSI2 solar
    !                                                  !        model) over solar cycle determined by 
    !                                                  !        direct specification of Mg (facular)
    !                                                  !        and SB (sunspot) indices provided
    !                                                  !        in indsolvar (scon = 0.0 only)
    !                                                  !    3 = (when scon .eq. 0.0): No solar variability
    !                                                  !        and no solar cycle (NRLSSI2 solar irradiance
    !                                                  !        of 1360.85 Wm-2 only);
    !                                                  !        (when scon .ne. 0.0): NRLSSI2 solar irradiance
    !                                                  !        scaled to scon and solar variability defined
    !                                                  !        (optional) by setting non-zero scale factors
    !                                                  !        for each band in bndsolvar
    isolvar = 0

    !  real(kind=rb), intent(in), optional :: indsolvar(:) ! Facular and sunspot amplitude 
    !                                                      ! scale factors (isolvar=1), or
    !                                                      ! Mg and SB indices (isolvar=2)
    !                                                      !    Dimensions: (2)
    indsolvar = 1

    !  real(kind=rb), intent(in), optional :: bndsolvar(:) ! Solar variability scale factors 
    !                                                      ! for each shortwave band
    !                                                      !    Dimensions: (nbndsw=14)
    bndsolvar = 1       ! no variability

    !  real(kind=rb), intent(in), optional :: solcycfrac   ! Fraction of averaged 11-year solar cycle (0-1)
    !                                                      !    at current time (isolvar=1)
    !                                                      !    0.0 represents the first day of year 1
    !                                                      !    1.0 represents the last day of year 11
    solcycfrac = 0.0

    !  integer(kind=im), intent(in) :: inflgsw         ! Flag for cloud optical properties
    inflgsw = 1
    !  integer(kind=im), intent(in) :: iceflgsw        ! Flag for ice particle specification
    iceflgsw = 1
    !  integer(kind=im), intent(in) :: liqflgsw        ! Flag for liquid droplet specification
    liqflgsw = 1

    !  real(kind=rb), intent(in) :: cldfmcl(:,:,:)     ! Cloud fraction
    !                                                  !    Dimensions: (ngptsw,ncol,nlay)
    cldfmcl = 0

    !  real(kind=rb), intent(in) :: taucmcl(:,:,:)     ! In-cloud optical depth
    !                                                  !    Dimensions: (ngptsw,ncol,nlay)
    taucmcl = 0

    !  real(kind=rb), intent(in) :: ssacmcl(:,:,:)     ! In-cloud single scattering albedo
    !                                                  !    Dimensions: (ngptsw,ncol,nlay)
    ssacmcl = 0

    !  real(kind=rb), intent(in) :: asmcmcl(:,:,:)     ! In-cloud asymmetry parameter
    !                                                  !    Dimensions: (ngptsw,ncol,nlay)
    asmcmcl = 0

    !  real(kind=rb), intent(in) :: fsfcmcl(:,:,:)     ! In-cloud forward scattering fraction
    !                                                  !    Dimensions: (ngptsw,ncol,nlay)
    fsfcmcl = 1

    !  real(kind=rb), intent(in) :: ciwpmcl(:,:,:)     ! In-cloud ice water path (g/m2)
    !                                                  !    Dimensions: (ngptsw,ncol,nlay)
    ciwpmcl = 0

    !  real(kind=rb), intent(in) :: clwpmcl(:,:,:)     ! In-cloud liquid water path (g/m2)
    !                                                  !    Dimensions: (ngptsw,ncol,nlay)
    clwpmcl = 0

    !  real(kind=rb), intent(in) :: reicmcl(:,:)       ! Cloud ice effective radius (microns)
    !                                                  !    Dimensions: (ncol,nlay)
    !                                                  ! specific definition of reicmcl depends on setting of iceflgsw:
    !                                                  ! iceflgsw = 0: (inactive)
    !                                                  ! 
    !                                                  ! iceflgsw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
    !                                                  !               r_ec range is limited to 13.0 to 130.0 microns
    !                                                  ! iceflgsw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
    !                                                  !               r_k range is limited to 5.0 to 131.0 microns
    !                                                  ! iceflgsw = 3: generalized effective size, dge, (Fu, 1996),
    !                                                  !               dge range is limited to 5.0 to 140.0 microns
    !                                                  !               [dge = 1.0315 * r_ec]
    reicmcl = 0

    !  real(kind=rb), intent(in) :: relqmcl(:,:)       ! Cloud water drop effective radius (microns)
    !                                                  !    Dimensions: (ncol,nlay)
    relqmcl = 0

    !  real(kind=rb), intent(in) :: tauaer(:,:,:)      ! Aerosol optical depth (iaer=10 only)
    !                                                  !    Dimensions: (ncol,nlay,nbndsw)
    !                                                  ! (non-delta scaled)      
    tauaer = 0

    !  real(kind=rb), intent(in) :: ssaaer(:,:,:)      ! Aerosol single scattering albedo (iaer=10 only)
    !                                                  !    Dimensions: (ncol,nlay,nbndsw)
    !                                                  ! (non-delta scaled)      
    ssaaer = 0

    !  real(kind=rb), intent(in) :: asmaer(:,:,:)      ! Aerosol asymmetry parameter (iaer=10 only)
    !                                                  !    Dimensions: (ncol,nlay,nbndsw)
    !                                                  ! (non-delta scaled)      
    asmaer = 0

    !  real(kind=rb), intent(in) :: ecaer(:,:,:)       ! Aerosol optical depth at 0.55 micron (iaer=6 only)
    !                                                  !    Dimensions: (ncol,nlay,naerec)
    !                                                  ! (non-delta scaled)      
    ecaer = 0


! --!--- Output -----

    !  real(kind=rb), intent(out) :: swuflx(:,:)       ! Total sky shortwave upward flux (W/m2)
    !                                                  !    Dimensions: (ncol,nlay+1)
    !  real(kind=rb), intent(out) :: swdflx(:,:)       ! Total sky shortwave downward flux (W/m2)
    !                                                  !    Dimensions: (ncol,nlay+1)
    !  real(kind=rb), intent(out) :: swhr(:,:)         ! Total sky shortwave radiative heating rate (K/d)
    !                                                  !    Dimensions: (ncol,nlay)
    !  real(kind=rb), intent(out) :: swuflxc(:,:)      ! Clear sky shortwave upward flux (W/m2)
    !                                                  !    Dimensions: (ncol,nlay+1)
    !  real(kind=rb), intent(out) :: swdflxc(:,:)      ! Clear sky shortwave downward flux (W/m2)
    !                                                  !    Dimensions: (ncol,nlay+1)
    !  real(kind=rb), intent(out) :: swhrc(:,:)        ! Clear sky shortwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)

    call  rrtmg_sw &
            (ncol    ,nlay    ,icld    ,iaer    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc   , &
             h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,scon    ,isolvar, &
             inflgsw ,iceflgsw,liqflgsw,cldfmcl , &
             taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl , &
             ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
! optional I/O
             bndsolvar,indsolvar,solcycfrac)

  end subroutine 

  subroutine rrtmg_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)     :: clock
    type(ESMF_Time)                        :: currTime

    rc = ESMF_SUCCESS

    deallocate(play)
    deallocate(plev)
    deallocate(tlay)
    deallocate(tlev)
    deallocate(tsfc)
    deallocate(h2ovmr)
    deallocate(o3vmr)
    deallocate(co2vmr)
    deallocate(ch4vmr)
    deallocate(n2ovmr)
    deallocate(o2vmr)
    deallocate(asdir)
    deallocate(aldir)
    deallocate(asdif)
    deallocate(aldif)
    deallocate(coszen)
    deallocate(cldfmcl)
    deallocate(taucmcl)
    deallocate(ssacmcl)
    deallocate(asmcmcl)
    deallocate(fsfcmcl)
    deallocate(ciwpmcl)
    deallocate(clwpmcl)
    deallocate(reicmcl)
    deallocate(relqmcl)
    deallocate(tauaer)
    deallocate(ssaaer)
    deallocate(asmaer)
    deallocate(ecaer)

    !write(*,*) swuflx
    !write(*,*) swdflx
    !write(*,*) swhr
    !write(*,*) swuflx
    !write(*,*) swdflx
    !write(*,*) swhrc

    deallocate(swuflx)
    deallocate(swdflx)
    deallocate(swhr)
    deallocate(swuflxc)
    deallocate(swdflxc)
    deallocate(swhrc)

    !call RRTMG_Finalize

  end subroutine rrtmg_model_finalize


  subroutine InitializeRealize(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("InitializeRealize -- ", ESMF_LOGMSG_INFO)

    ! As a radiation component, check no grid is provided by radiation component
    ! First examine import State
    call RealizeFields(importState, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call RealizeFields(exportState, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    contains

    subroutine RealizeFields(State, rc)
    type(ESMF_State)     :: State
    integer, intent(out) :: rc
    ! local variables    
    type(ESMF_Field)                       :: field
    character(ESMF_MAXSTR)                 :: transferAction
    integer                                :: i, icount
    character(64), allocatable             :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: typeList(:)

    call ESMF_StateGet(State, itemCount=icount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(typeList(icount), itemNameList(icount))
    call ESMF_StateGet(State, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i = 1, icount
      if(typeList(i) == ESMF_STATEITEM_FIELD) then
        call ESMF_LogWrite("Realize Field Name Initiated: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)
        ! This Field was marked with TransferOfferGeomObject="can provide", so here
        ! we need to see what TransferActionGeomObject the Connector determined for
        ! this Field:
        call ESMF_StateGet(State, itemName=itemNameList(i), field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
    
        call NUOPC_GetAttribute(field, name="TransferActionGeomObject", &
          value=transferAction, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (trim(transferAction)=="provide") then
          ! the Connector instructed the RRTMG to provide the Grid object for "player"
          call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="RRTMG cannot provide any Grid", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        else  ! transferAction=="accept"
          ! the Connector instructed the RRTMG to accept the Grid from RTM for "player"
          call ESMF_LogWrite("RRTMG is accepting Grid for Field", &
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        endif
      endif
    enddo
    deallocate(typeList, itemNameList)

    end subroutine
  end subroutine

  subroutine InitializeAcceptGrid(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(ESMF_Grid)      :: Grid
    type(ESMF_Mesh)      :: Mesh
    type(ESMF_DistGrid)  :: distGrid
    integer              :: lsize
    integer, allocatable              :: seqIndexList(:)
    
    rc = ESMF_SUCCESS
    call ESMF_LogWrite("InitializeAcceptGrid -- ", ESMF_LOGMSG_INFO)

#ifdef USE_GRID
    call AcceptGrid(importState, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call AcceptGrid(exportState, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

#ifdef USE_MESH
    !AcceptMesh(importState, rc)
    !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    !  line=__LINE__, &
    !  file=__FILE__)) &
    !  return  ! bail out
    Mesh=AcceptMesh(exportState, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_MeshGet(mesh, elementDistGrid=distgrid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_DistgridGet(distgrid, localDe=0, elementCount=lsize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    allocate(seqIndexList(lsize))
    call ESMF_DistgridGet(distgrid, localDe=0, seqIndexList=seqIndexList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    !write (*, *) seqIndexList
    deallocate(seqIndexList)

    ncol = lsize
    write (*, *) 'ncol = ', ncol
#endif

    call RRTMG_SW_INI(cpdair)

    allocate(play(ncol, nlay))
    allocate(plev(ncol, nlay+1))
    allocate(tlay(ncol, nlay))
    allocate(tlev(ncol, nlay+1))
    allocate(tsfc(ncol))
    allocate(h2ovmr(ncol, nlay))
    allocate(o3vmr(ncol, nlay))
    allocate(co2vmr(ncol, nlay))
    allocate(ch4vmr(ncol, nlay))
    allocate(n2ovmr(ncol, nlay))
    allocate(o2vmr(ncol, nlay))
    allocate(asdir(ncol))
    allocate(aldir(ncol))
    allocate(asdif(ncol))
    allocate(aldif(ncol))
    allocate(coszen(ncol))
    allocate(indsolvar(2))
    allocate(bndsolvar(14))         ! 14 shortwave bands
    allocate(cldfmcl(ngptsw, ncol, nlay))
    allocate(taucmcl(ngptsw, ncol, nlay))
    allocate(ssacmcl(ngptsw, ncol, nlay))
    allocate(asmcmcl(ngptsw, ncol, nlay))
    allocate(fsfcmcl(ngptsw, ncol, nlay))
    allocate(ciwpmcl(ngptsw, ncol, nlay))
    allocate(clwpmcl(ngptsw, ncol, nlay))
    allocate(reicmcl(ncol, nlay))
    allocate(relqmcl(ncol, nlay))
    allocate(tauaer(ncol, nlay, nbndsw))
    allocate(ssaaer(ncol, nlay, nbndsw))
    allocate(asmaer(ncol, nlay, nbndsw))
    allocate(ecaer(ncol, nlay, nbndsw))

    allocate(swuflx(ncol, nlay+1))
    allocate(swdflx(ncol, nlay+1))
    allocate(swhr(ncol, nlay))
    allocate(swuflxc(ncol, nlay+1))
    allocate(swdflxc(ncol, nlay+1))
    allocate(swhrc(ncol, nlay))

    call ESMF_LogWrite("RRTMG - InitializeP4: DONE", &
       ESMF_LOGMSG_INFO, rc=rc)

    contains 
  
      subroutine AcceptGrid(State, rc)

      type(ESMF_State)     :: State
      integer, intent(out) :: rc
      
      ! local variables
      type(ESMF_Field)              :: field
      type(ESMF_Grid)               :: grid
      integer                       :: localDeCount
      character(80)                 :: name
      character(160)                :: msgString

      type(ESMF_DistGrid)           :: distgrid
      integer                       :: dimCount, tileCount, arbDimCount
      integer, allocatable          :: minIndexPTile(:,:), maxIndexPTile(:,:)
      integer                       :: connectionCount
      type(ESMF_DistGridConnection), allocatable :: connectionList(:)
      logical                       :: regDecompFlag

      ! local variables    
      character(ESMF_MAXSTR)                 :: transferAction
      integer                                :: i, icount
      character(64), allocatable             :: itemNameList(:)
      type(ESMF_StateItem_Flag), allocatable :: typeList(:)

      call ESMF_StateGet(state, itemCount=icount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(typeList(icount), itemNameList(icount))
      call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      do i = 1, icount
        if(typeList(i) == ESMF_STATEITEM_FIELD) then
          call ESMF_LogWrite("Accept Grid Initiated: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)
      
        ! access the field in the State
        call ESMF_StateGet(State, field=field, itemName=itemNameList(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! construct a local Grid according to the transferred grid
        call ESMF_FieldGet(field, grid=grid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_GridGet(grid, distgrid=distgrid, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_DistGridGet(distgrid, dimCount=dimCount, tileCount=tileCount, &
          connectionCount=connectionCount, regDecompFlag=regDecompFlag, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        allocate(minIndexPTile(dimCount, tileCount), &
          maxIndexPTile(dimCount, tileCount))
        allocate(connectionList(connectionCount))
        call ESMF_DistGridGet(distgrid, minIndexPTile=minIndexPTile, &
          maxIndexPTile=maxIndexPTile, connectionList=connectionList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (regDecompFlag) then
          ! The provider used a regDecomp scheme for DistGrid creation:
          ! This means that the entire index space is covered (no holes), and
          ! it the easieast is just to use a regDecomp scheme on the acceptor
          ! side as well.
          distgrid = ESMF_DistGridCreate(minIndexPTile=minIndexPTile, &
            maxIndexPTile=maxIndexPTile, connectionList=connectionList, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          grid = ESMF_GridCreate(distgrid, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          ! swap out the transferred grid for the newly created one
          call ESMF_FieldEmptySet(field, grid=grid, rc=rc)    
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call ESMF_LogWrite("RRTMG - Just set Grid for Field"//trim(itemNameList(i)), &
            ESMF_LOGMSG_INFO, rc=rc)
        endif
        deallocate(minIndexPTile, maxIndexPTile, connectionList)
      endif

      enddo
      deallocate(typeList, itemNameList)

      end subroutine
      !-----------------------------------------------------------------------------
      function AcceptMesh(State, rc)

      type(ESMF_Mesh)      :: AcceptMesh

      type(ESMF_State)     :: State
      integer, intent(out) :: rc
      
      ! local variables
      type(ESMF_Field)              :: field
      type(ESMF_Mesh)               :: mesh
      character(80)                 :: name
      character(160)                :: msgString

      ! local variables    
      character(ESMF_MAXSTR)                 :: transferAction
      integer                                :: i, icount
      character(64), allocatable             :: itemNameList(:)
      type(ESMF_StateItem_Flag), allocatable :: typeList(:)
      logical                                :: foundMesh = .false.
      integer                                :: noe, mipde(1,4), ecpde(4)
      type(ESMF_DistGrid)                    :: ndg, edg

      call ESMF_StateGet(state, itemCount=icount, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      allocate(typeList(icount), itemNameList(icount))
      call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      do i = 1, icount
        if(typeList(i) == ESMF_STATEITEM_FIELD) then
          call ESMF_LogWrite("Accept Mesh Initiated: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)
      
          ! access the field in the State
          call ESMF_StateGet(State, field=field, itemName=itemNameList(i), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          ! construct a local Grid according to the transferred grid
          call ESMF_FieldGet(field, mesh=mesh, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          !call ESMF_FieldEmptySet(field, mesh=mesh, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          !  line=__LINE__, &
          !  file=__FILE__)) &
          !  return  ! bail out
          !call ESMF_MeshGet(mesh, nodalDistGrid=ndg, elementDistGrid=edg, numOwnedElements=noe, rc=rc)
          !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          !  line=__LINE__, &
          !  file=__FILE__)) &
          !  return  ! bail out
          !call ESMF_DistGridGet(ndg, maxIndexPDe=mipde, elementCountPDe=ecpde, rc=rc)
          !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          !  line=__LINE__, &
          !  file=__FILE__)) &
          !  return  ! bail out
          !call ESMF_DistGridGet(edg, maxIndexPDe=mipde, elementCountPDe=ecpde, rc=rc)
          !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          !  line=__LINE__, &
          !  file=__FILE__)) &
          !  return  ! bail out
          !  ! swap out the transferred grid for the newly created one
          call ESMF_LogWrite("RRTMG - Just set Mesh for Field"//trim(itemNameList(i)), &
            ESMF_LOGMSG_INFO, rc=rc)
          foundMesh = .true.
        endif
      enddo
      deallocate(typeList, itemNameList)

      if(foundMesh) then
        AcceptMesh = mesh
      else
        call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
            msg="RRTMG cannot provide any Mesh", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
      endif

    end function

  end subroutine
    
  !-----------------------------------------------------------------------------

  subroutine InitializeCompleteField(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    call ESMF_LogWrite("InitializeCompleteField -- ", ESMF_LOGMSG_INFO)
#ifdef USE_GRID
    call CompleteFieldGrid(importState, writeGrid=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call CompleteFieldGrid(exportState, writeGrid=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif
#ifdef USE_MESH
    call CompleteFieldMesh(exportState, writeMesh=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    contains

    subroutine CompleteFieldGrid(State, writeGrid, rc)
    type(ESMF_State)              :: State
    logical, intent(in), optional :: writeGrid
    integer, intent(out)          :: rc
    
    ! local variables
    type(ESMF_Field)              :: field
    type(ESMF_Grid)               :: grid
    type(ESMF_Array)              :: array
    character(80)                 :: name
    character(160)                :: msgString
    type(ESMF_FieldStatus_Flag)   :: fieldStatus
    integer                       :: staggerEdgeLWidth(2)
    integer                       :: staggerEdgeUWidth(2)
    integer                       :: staggerAlign(2)

    integer                                :: i, icount
    character(64), allocatable             :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: typeList(:)
    logical                                :: l_writeGrid = .false.
    type(FieldListType)                    :: fieldList(6)

    fieldList(1)%farrayPtr2D => swuflx
    fieldList(2)%farrayPtr2D => swdflx
    fieldList(3)%farrayPtr2D => swhr
    fieldList(4)%farrayPtr2D => swuflxc
    fieldList(5)%farrayPtr2D => swdflxc
    fieldList(6)%farrayPtr2D => swhrc

    call ESMF_StateGet(state, itemCount=icount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(typeList(icount), itemNameList(icount))
    call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i = 1, icount
      if(typeList(i) == ESMF_STATEITEM_FIELD) then
        call ESMF_LogWrite("Complete Field Name Initiated: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)

        ! access the field in the State
        call ESMF_StateGet(State, field=field, itemName=itemNameList(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! check status of field and decide on action
        call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (fieldStatus==ESMF_FIELDSTATUS_COMPLETE) then
          ! log a message
          call ESMF_LogWrite("RRTMG - The Field was already complete", &
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        else
          ! the transferred Grid is already set, allocate memory for data by complete
          print *, ubound(fieldList(i)%farrayPtr2D, 1), ubound(fieldList(i)%farrayPtr2D, 2)
          call ESMF_FieldEmptyComplete(field, farrayPtr=fieldList(i)%farrayPtr2D, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          ! log a message
          call ESMF_LogWrite("RRTMG - Just completed the Field", &
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        endif
      endif
    enddo
    deallocate(typeList, itemNameList)

    if(present(writeGrid)) l_writeGrid = writeGrid

    if(l_writeGrid) then
    ! Use the last Field to dump the Grid information
#ifdef TEST_MULTI_TILE_GRID    
    ! write cubed sphere grid out to VTK
    call ESMF_FieldGet(field, grid=grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridWriteVTK(grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
      filename="RRTMG-accepted-Grid-ssh_centers", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    ! inspect the Grid name
    call ESMF_FieldGet(field, grid=grid, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGet(grid, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    write (msgString,*) "RRTMG - InitializeP5: transferred Grid name = ", name
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! check the staggerEdgeWidth of the transferred grid 
#ifdef TEST_GRID_EDGE_WIDTHS
    ! center stagger
    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
      staggerEdgeLWidth=staggerEdgeLWidth, &
      staggerEdgeUWidth=staggerEdgeUWidth, &
      staggerAlign=staggerAlign, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    print *, "staggerEdgeLWidth", staggerEdgeLWidth
    print *, "staggerEdgeUWidth", staggerEdgeUWidth
    print *, "staggerAlign", staggerAlign
    if (any(staggerEdgeLWidth /= (/0,0/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeLWidth for ESMF_STAGGERLOC_CENTER", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
    if (any(staggerEdgeUWidth /= (/0,0/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeUWidth for ESMF_STAGGERLOC_CENTER", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
    ! corner stagger
    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_CORNER, &
      staggerEdgeLWidth=staggerEdgeLWidth, &
      staggerEdgeUWidth=staggerEdgeUWidth, &
      staggerAlign=staggerAlign, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    print *, "staggerEdgeLWidth", staggerEdgeLWidth
    print *, "staggerEdgeUWidth", staggerEdgeUWidth
    print *, "staggerAlign", staggerAlign
    if (any(staggerEdgeLWidth /= (/1,1/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeLWidth for ESMF_STAGGERLOC_CORNER", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
    if (any(staggerEdgeUWidth /= (/0,0/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeuWidth for ESMF_STAGGERLOC_CORNER", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
    ! edge1 stagger
    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_EDGE1, &
      staggerEdgeLWidth=staggerEdgeLWidth, &
      staggerEdgeUWidth=staggerEdgeUWidth, &
      staggerAlign=staggerAlign, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    print *, "staggerEdgeLWidth", staggerEdgeLWidth
    print *, "staggerEdgeUWidth", staggerEdgeUWidth
    print *, "staggerAlign", staggerAlign
    if (any(staggerEdgeLWidth /= (/0,1/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeLWidth for ESMF_STAGGERLOC_EDGE1", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
    if (any(staggerEdgeUWidth /= (/1,1/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeUWidth for ESMF_STAGGERLOC_EDGE1", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
    ! edge2 stagger
    call ESMF_GridGet(grid, staggerloc=ESMF_STAGGERLOC_EDGE2, &
      staggerEdgeLWidth=staggerEdgeLWidth, &
      staggerEdgeUWidth=staggerEdgeUWidth, &
      staggerAlign=staggerAlign, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    print *, "staggerEdgeLWidth", staggerEdgeLWidth
    print *, "staggerEdgeUWidth", staggerEdgeUWidth
    print *, "staggerAlign", staggerAlign
    if (any(staggerEdgeLWidth /= (/1,0/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeLWidth for ESMF_STAGGERLOC_EDGE2", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
    if (any(staggerEdgeUWidth /= (/0,1/))) then
      call ESMF_LogSetError(ESMF_RC_VAL_WRONG, &
        msg="Wrong staggerEdgeUWidth for ESMF_STAGGERLOC_EDGE2", &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)
      return  ! bail out
    endif
#endif

#if 1
    ! testing the output of coord arrays
    !TODO:
    ! Coords are currently written in 2D index space even if there is coordinate
    ! factorization used, e.g. in the Ufrm() GridCreate. Therefore the coord
    ! arrays have replicated dims, and underlying allocation is only 1D. This 
    ! should be changed in the ArrayWrite() where Arrays with replicated dims
    ! should write out only the non-degenerate data, i.e. according to the 
    ! actual data allocation. 
    ! -> here that would be a 1D array for each coordiante dim.
    ! center:
    call ESMF_GridGetCoord(grid, coordDim=1, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_center_coord1.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(grid, coordDim=2, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_center_coord2.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#ifdef TEST_GRID_EDGE_WIDTHS
    ! corner:
    call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_CORNER, &
      coordDim=1, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_corner_coord1.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_CORNER, &
      coordDim=2, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_corner_coord2.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! edge1:
    call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE1, &
      coordDim=1, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_edge1_coord1.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE1, &
      coordDim=2, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_edge1_coord2.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! edge2:
    call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE2, &
      coordDim=1, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_edge2_coord1.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(grid, staggerloc=ESMF_STAGGERLOC_EDGE2, &
      coordDim=2, array=array, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_ArrayWrite(array, "array_RRTMG-grid_edge2_coord2.nc", overwrite=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif
#endif

#if 1
    ! write out the Grid into VTK file for inspection
    call ESMF_GridWriteVTK(grid, staggerloc=ESMF_STAGGERLOC_CENTER, &
      filename="RRTMG-accepted-Grid-centers", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite("Done writing RRTMG-accepted-Grid-centers VTK", &
      ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    endif  ! if(l_writeGrid)
    end subroutine  ! internal subroutine CompleteField

    subroutine CompleteFieldMesh(State, writeMesh, rc)
    type(ESMF_State)              :: State
    logical, intent(in), optional :: writeMesh
    integer, intent(out)          :: rc
    
    ! local variables
    type(ESMF_Field)              :: field
    type(ESMF_Mesh)               :: mesh
    type(ESMF_FieldStatus_Flag)   :: fieldStatus

    integer                                :: i, icount
    character(64), allocatable             :: itemNameList(:)
    type(ESMF_StateItem_Flag), allocatable :: typeList(:)
    logical                                :: l_writeMesh = .false.
    type(FieldListType)                    :: fieldList(6)

    ! Names are in alphabetical order
    fieldList(1)%farrayPtr2D => swdflx
    fieldList(2)%farrayPtr2D => swhr
    fieldList(3)%farrayPtr2D => swuflx
    fieldList(4)%farrayPtr2D => swdflxc
    fieldList(5)%farrayPtr2D => swhrc
    fieldList(6)%farrayPtr2D => swuflxc

    call ESMF_StateGet(state, itemCount=icount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    allocate(typeList(icount), itemNameList(icount))
    call ESMF_StateGet(state, itemTypeList=typeList, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i = 1, icount
      if(typeList(i) == ESMF_STATEITEM_FIELD) then
        call ESMF_LogWrite("Complete Field Name Initiated: "//trim(itemNameList(i)), ESMF_LOGMSG_INFO)

        ! access the field in the State
        call ESMF_StateGet(State, field=field, itemName=itemNameList(i), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! check status of field and decide on action
        call ESMF_FieldGet(field, status=fieldStatus, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        if (fieldStatus==ESMF_FIELDSTATUS_COMPLETE) then
          ! log a message
          call ESMF_LogWrite("RRTMG - The Field was already complete", &
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        else
          ! the transferred Mesh is already set, allocate memory for data by complete
          print *, ubound(fieldList(i)%farrayPtr2D, 1), ubound(fieldList(i)%farrayPtr2D, 2)
          call ESMF_FieldGet(field, mesh=mesh, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call ESMF_FieldEmptySet(field, mesh=mesh, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
          call ESMF_FieldEmptyComplete(field, farrayPtr=fieldList(i)%farrayPtr2D, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out

          !call ESMF_FieldGet(field, mesh=mesh, rc=rc)
          !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          !  line=__LINE__, &
          !  file=__FILE__)) &
          !  return  ! bail out
          !field = ESMF_FieldCreate(mesh, farrayPtr=swuflx, meshloc=ESMF_MESHLOC_ELEMENT, name=exportFieldList(i), rc=rc)
          !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          !  line=__LINE__, &
          !  file=__FILE__)) &
          !  return  ! bail out
          !call NUOPC_Realize(State, field, rc=rc)
          !if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          !  line=__LINE__, &
          !  file=__FILE__)) &
          !  return  ! bail out
          ! log a message
          call ESMF_LogWrite("RRTMG - Just completed the Field", &
            ESMF_LOGMSG_INFO, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, &
            file=__FILE__)) &
            return  ! bail out
        endif
      endif
    enddo
    deallocate(typeList, itemNameList)

    if(present(writeMesh)) l_writeMesh = writeMesh

    if(l_writeMesh) then
    call ESMF_FieldGet(field, Mesh=mesh, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! write out the Grid into VTK file for inspection
    call ESMF_MeshWrite(mesh, &
      filename="RRTMG-accepted-Mesh", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite("Done writing RRTMG-accepted-Mesh VTK", &
      ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    endif  ! if(l_writeMesh)
  end subroutine  ! internal subroutine CompleteField

  end subroutine

  subroutine DataInitialize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! indicate that data initialization is complete (breaking out of init-loop)
    call NUOPC_CompAttributeSet(model, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  subroutine RRTMG_FieldsSetup(rc)

    integer, intent(out)                   :: rc
    rc = ESMF_SUCCESS

  end subroutine
  !-----------------------------------------------------------------------------
end module
