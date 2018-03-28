!--------------- RRTMG NUOPC CAP -----------------
! This is the RRTMG model cap component that's NUOPC compiant.
!
! Author:  Fei.Liu@NOAA.GOV
!
! 3/22/18
! This is now acting as a cap/connector between NUOPC driver and RRTMG code.
!

module rrtmg_sw_cap_mod

  use rrtmg_sw_init
  use rrtmg_sw_rad

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS      => SetServices, &
    model_label_CheckImport  => label_CheckImport, &
    model_label_Advance   => label_Advance, &
    model_label_Finalize  => label_Finalize

  implicit none
  private
  public SetServices

  !https://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
  real(kind=rb)      :: cpdair =  1004      ! Cp = 1004 J/kg K at 298K

  integer(kind=im), parameter :: ncol = 1        ! Number of horizontal columns     
  integer(kind=im), parameter :: nlay = 50       ! Number of model layers
  !integer(kind=im), parameter :: ngptsw = 112

  integer(kind=im)  :: icld            ! Cloud overlap method
  integer(kind=im)  :: iaer            ! Aerosol option flag
  real              :: p0 = 1013.25    ! 1 atm
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
  integer(kind=im)       :: dyofyr          ! Day of the year (used to get Earth/Sun
  real(kind=rb)          :: adjes              ! Flux adjustment for Earth/Sun distance
  real(kind=rb), pointer :: coszen(:)          ! Cosine of solar zenith angle
  real(kind=rb)          :: scon               ! Solar constant (W/m2)
  integer(kind=im)       :: isolvar         ! Flag for solar variability method
  real(kind=rb), pointer :: indsolvar(:) ! Facular and sunspot amplitude 
  real(kind=rb), pointer :: bndsolvar(:) ! Solar variability scale factors 
  real(kind=rb)          :: solcycfrac   ! Fraction of averaged 11-year solar cycle (0-1)
  integer(kind=im)       :: inflgsw         ! Flag for cloud optical properties
  integer(kind=im)       :: iceflgsw        ! Flag for ice particle specification
  integer(kind=im)       :: liqflgsw        ! Flag for liquid droplet specification
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

  integer   :: import_slice = 0
  integer   :: export_slice = 0

  type fld_list_type
    character(len=64) :: stdname
    character(len=64) :: shortname
    character(len=64) :: transferOffer
    logical           :: assoc    ! is the farrayPtr associated with internal data
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: farrayPtr
  end type fld_list_type

  integer,parameter :: fldsMax = 100
  integer :: fldsToRRTMG_num = 0
  type (fld_list_type) :: fldsToRRTMG(fldsMax)
  integer :: fldsFrRRTMG_num = 0
  type (fld_list_type) :: fldsFrRRTMG(fldsMax)

  character(len=256) :: tmpstr
  character(len=2048):: info
  logical :: isPresent
  integer :: dbrc     ! temporary debug rc value

  contains
  !-----------------------------------------------------------------------
  !------------------- RRTMG code starts here -----------------------
  !-----------------------------------------------------------------------

  subroutine SetServices(gcomp, rc)

    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    character(len=*),parameter  :: subname='(rrtmg_cap:SetServices)'

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
      phaseLabelList=(/"IPDv01p1"/), userRoutine=InitializeAdvertise, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv01p3"/), userRoutine=InitializeRealize, rc=rc)
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


    !call RRTMG_FieldsSetup()

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    rc = ESMF_SUCCESS

    ! Switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(gcomp, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv01p"/), rc=rc)
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
    character(len=*),parameter  :: subname='(rrtmg_cap:InitializeAdvertise)'

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
    call RRTMG_AdvertiseFields(importState, fldsToRRTMG_num, fldsToRRTMG, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call RRTMG_AdvertiseFields(exportState, fldsFrRRTMG_num, fldsFrRRTMG, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- initialization phase 1 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeRealize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! Local Variables
    type(ESMF_VM)                          :: vm
    type(ESMF_Grid)                        :: gridIn
    type(ESMF_Grid)                        :: gridOut
    character(len=*),parameter  :: subname='(rrtmg_cap:InitializeRealize)'

    rc = ESMF_SUCCESS

    ! We can check if npet is 4 or some other value to make sure
    ! RRTMG is configured to run on the correct number of processors.

    ! create a Grid object for Fields

    gridOut = gridIn ! for now out same as in

    call RRTMG_RealizeFields(importState, gridIn , fldsToRRTMG_num, fldsToRRTMG, "RRTMG import", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call RRTMG_RealizeFields(exportState, gridOut, fldsFrRRTMG_num, fldsFrRRTMG, "RRTMG export", rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write(info,*) subname,' --- initialization phase 2 completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, line=__LINE__, file=__FILE__, rc=dbrc)

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
    character(len=*),parameter  :: subname='(rrtmg_cap:rrtmg_model_advance)'

    rc = ESMF_SUCCESS
    write(info,*) subname,' --- run phase 1 called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
    

    import_slice = import_slice + 1
    export_slice = export_slice + 1

    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
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

    write(info,*) subname,' --- run phase 1 completed --- '

  end subroutine 

  subroutine rrtmg_model_finalize(gcomp, rc)

    ! input arguments
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)     :: clock
    type(ESMF_Time)                        :: currTime
    character(len=*),parameter  :: subname='(rrtmg_cap:rrtmg_model_finalize)'

    rc = ESMF_SUCCESS

    write(info,*) subname,' --- finalize called --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

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

    write(*,*) swuflx
    write(*,*) swdflx
    write(*,*) swhr
    write(*,*) swuflx
    write(*,*) swdflx
    write(*,*) swhrc

    deallocate(swuflx)
    deallocate(swdflx)
    deallocate(swhr)
    deallocate(swuflxc)
    deallocate(swdflxc)
    deallocate(swhrc)

    !call RRTMG_Finalize

    write(info,*) subname,' --- finalize completed --- '
    call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)

  end subroutine rrtmg_model_finalize

  subroutine RRTMG_AdvertiseFields(state, nfields, field_defs, rc)

    type(ESMF_State), intent(inout)             :: state
    integer,intent(in)                          :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    integer, intent(inout)                      :: rc

    integer                                     :: i
    character(len=*),parameter  :: subname='(rrtmg_cap:RRTMG_AdvertiseFields)'

    rc = ESMF_SUCCESS

    do i = 1, nfields

      call NUOPC_Advertise(state, &
        standardName=field_defs(i)%stdname, &
        name=field_defs(i)%shortname, &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo

  end subroutine RRTMG_AdvertiseFields

  subroutine RRTMG_RealizeFields(state, grid, nfields, field_defs, tag, rc)

    type(ESMF_State), intent(inout)             :: state
    type(ESMF_Grid), intent(in)                 :: grid
    integer, intent(in)                         :: nfields
    type(fld_list_type), intent(inout)          :: field_defs(:)
    character(len=*), intent(in)                :: tag
    integer, intent(inout)                      :: rc

    integer                                     :: i
    type(ESMF_Field)                            :: field
    integer                                     :: npet, nx, ny, pet, elb(2), eub(2), clb(2), cub(2), tlb(2), tub(2)
    type(ESMF_VM)                               :: vm
    character(len=*),parameter  :: subname='(rrtmg_cap:RRTMG_RealizeFields)'
 
    rc = ESMF_SUCCESS

      !call ESMF_VMGetCurrent(vm, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !call ESMF_VMGet(vm, petcount=npet, localPet=pet, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !call ESMF_GridGet(grid, exclusiveLBound=elb, exclusiveUBound=eub, &
      !                        computationalLBound=clb, computationalUBound=cub, &
      !                        totalLBound=tlb, totalUBound=tub, rc=rc)
      !if (rc /= ESMF_SUCCESS) call ESMF_Finalize()

      !write(info, *) pet, 'exc', elb, eub, 'comp', clb, cub, 'total', tlb, tub
      !call ESMF_LogWrite(subname // tag // " Grid "// info, &
      !  ESMF_LOGMSG_INFO, &
      !  line=__LINE__, &
      !  file=__FILE__, &
      !  rc=dbrc)

    do i = 1, nfields

      if (field_defs(i)%assoc) then
        write(info, *) subname, tag, ' Field ', field_defs(i)%shortname, ':', &
          lbound(field_defs(i)%farrayPtr,1), ubound(field_defs(i)%farrayPtr,1), &
          lbound(field_defs(i)%farrayPtr,2), ubound(field_defs(i)%farrayPtr,2), &
          lbound(field_defs(i)%farrayPtr,3), ubound(field_defs(i)%farrayPtr,3)
        call ESMF_LogWrite(info, ESMF_LOGMSG_INFO, rc=dbrc)
        field = ESMF_FieldCreate(grid=grid, &
          farray=field_defs(i)%farrayPtr, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, indexflag=ESMF_INDEX_DELOCAL, &
          name=field_defs(i)%shortname, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

      if (NUOPC_IsConnected(state, fieldName=field_defs(i)%shortname)) then
        call NUOPC_Realize(state, field=field, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
!        call ESMF_FieldPrint(field=field, rc=rc)
!        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!          line=__LINE__, &
!          file=__FILE__)) &
!          return  ! bail out
      else
        call ESMF_LogWrite(subname // tag // " Field "// field_defs(i)%stdname // " is not connected.", &
          ESMF_LOGMSG_INFO, &
          line=__LINE__, &
          file=__FILE__, &
          rc=dbrc)
        ! TODO: Initialize the value in the pointer to 0 after proper restart is setup
        !if(associated(field_defs(i)%farrayPtr) ) field_defs(i)%farrayPtr = 0.0
        ! remove a not connected Field from State
        call ESMF_StateRemove(state, (/field_defs(i)%shortname/), rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif

    enddo


  end subroutine RRTMG_RealizeFields

  !-----------------------------------------------------------------------------

  subroutine dumpRRTMGInternal(grid, slice, stdname, nop, farray)

    type(ESMF_Grid)          :: grid
    integer, intent(in)      :: slice
    character(len=*)         :: stdname
    character(len=*)         :: nop
    real(ESMF_KIND_R8), dimension(:,:,:), target :: farray

    type(ESMF_Field)         :: field
    real(ESMF_KIND_R8), dimension(:,:), pointer  :: f2d
    integer                  :: rc

    return ! remove this line to debug field connection

    field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, &
      indexflag=ESMF_INDEX_DELOCAL, &
      name=stdname, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldGet(field, farrayPtr=f2d, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    f2d(:,:) = farray(:,:,1)

    call ESMF_FieldWrite(field, filename='field_rrtmg_internal_'//trim(stdname)//'.nc', &
      timeslice=slice, rc=rc) 
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_FieldDestroy(field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------
end module
