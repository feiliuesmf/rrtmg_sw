!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2018, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module ATM

  !-----------------------------------------------------------------------------
  ! ATM Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    model_routine_SS            => SetServices, &
    model_label_DataInitialize  => label_DataInitialize

  implicit none
  
  private
  
  public SetServices

  TYPE(ESMF_VM) :: vm
  INTEGER :: lPet
  integer :: importSlice = 0

  ! Field Entries are in alphabetical order
  character(len=45)     :: importFieldList(6) = (/ &
    "Total Sky Shortward Downard Flux            ", &
    "Total Sky Shortward Radiative Heating Rate  ", &
    "Total Sky Shortwave Upward Flux             ", &
    "Clear Sky Shortwave Downward Flux           ", &
    "Clear Sky Shortwave Radiative Heating Rate  ", &
    "Clear Sky Shortwave Upward Flux             " /)
  character(len=8)      :: importFieldSN(6) = (/ &
    "tssdf   ", &
    "tssrhr  ", &
    "tssuf   ", &
    "cssdf   ", &
    "cssrhr  ", &
    "cssuf   " /)
  character(len=55), parameter     :: exportFieldList(6) = (/ &
  "Layer pressures                                       ",&
  "Interface pressures                                   ",&
  "Layer temperatures                                    ",&
  "Interface temperatures                                ",&
  "Surface temperature                                   ",&
  "H2O volume mixing ratio                               " /)
  character(len=8), parameter     :: exportFieldSN(6) = (/ &
    "Pa      ", &
    "Pa      ", &
    "K       ", &
    "K       ", &
    "K       ", &
    "1       " /)
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    call ESMF_MethodRemove(model, label_SetRunClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=label_SetRunClock, &
      specRoutine=SetRunClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    call ESMF_MethodRemove(model, label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_CheckImport, &
      specRoutine=NUOPC_NOOP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_DataInitialize, &
      specRoutine=DataInitialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ATM calling InitializeP1 -- ", ESMF_LOGMSG_INFO)
    
    ! ATM can and will provide Mesh for the RTM component
    ! Import
    call NUOPC_Advertise(importState, &
      standardNames=importFieldList, &
      TransferOfferGeomObject="can provide", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! Export
    call NUOPC_Advertise(exportState, &
      standardNames=exportFieldList, &
      TransferOfferGeomObject="can provide", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP2(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field)        :: field
    type(ESMF_Mesh)         :: mesh
    integer                 :: i, atm_layers
    character(len=10)       :: value
    ! Field Entries are in alphabetical order
    integer                 :: nz(6)=(/1,0,1,1,0,1/)
    integer                 :: nz1(6)=(/0,1,0,1,-1000,0/)
    
    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ATM calling InitializeP2 -- ", ESMF_LOGMSG_INFO)
    
    ! create a Mesh object for Fields
    mesh = ESMF_MeshCreate("data/ll2.5deg_grid.nc", ESMF_FILEFORMAT_SCRIP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_AttributeGet(model, name="atm_layers", value=value, defaultValue="50", &
      convention="NUOPC", purpose="Instance", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    atm_layers = ESMF_UtilString2Int(value, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    nz = nz + atm_layers
    nz1 = nz1 + atm_layers

    do i = 1, 6
      ! Vertical levels are hardcoded on the ATM side
      field = ESMF_FieldCreate(name=importFieldList(i), mesh=mesh, &
        meshloc=ESMF_MESHLOC_ELEMENT, typekind=ESMF_TYPEKIND_R8, &
        ungriddedLBound=(/1/), ungriddedUBound=(/nz(i)/), & 
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_Realize(importState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo
    do i = 1, 6
      if(nz1(i) .gt. 0) then
        ! Fields with vertical dimension
        field = ESMF_FieldCreate(name=exportFieldList(i), mesh=mesh, &
          meshloc=ESMF_MESHLOC_ELEMENT, typekind=ESMF_TYPEKIND_R8, &
          ungriddedLBound=(/1/), ungriddedUBound=(/nz1(i)/), & 
          rc=rc)
      else
        ! Fields without vertical dimension
        field = ESMF_FieldCreate(name=exportFieldList(i), mesh=mesh, &
          meshloc=ESMF_MESHLOC_ELEMENT, typekind=ESMF_TYPEKIND_R8, &
          rc=rc)
      endif
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_Realize(exportState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    enddo

  end subroutine

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_State)              :: importState, exportState

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.

    importSlice = importSlice + 1
    call NUOPC_Write(importState, filenamePrefix="atm_import_", &
      overwrite=.true., timeslice=importSlice, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    CALL ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    CALL ESMF_VMGet(vm, localPet=lPet,rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if(lPet.eq.0) then
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing ATM from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="--------------------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    endif

  end subroutine


  !-----------------------------------------------------------------------------

  subroutine SetRunClock(model, rc)

    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    ! local variables
    type(ESMF_Clock)          :: clock, driverClock
    type(ESMF_Time)           :: checkCurrTime, currTime, stopTime
    type(ESMF_TimeInterval)   :: checkTimeStep, timeStep
    type(ESMF_Direction_Flag) :: direction

    rc = ESMF_SUCCESS


    ! query component for clock and driver clock
    call NUOPC_ModelGet(model, modelClock=clock, driverClock=driverClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(driverClock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockSet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! check and set the component clock against the driver clock
    call NUOPC_CompCheckSetClock(model, driverClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="NUOPC INCOMPATIBILITY DETECTED: between model and driver clocks", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetRunClock

  subroutine DataInitialize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS
    call ESMF_LogWrite("ATM calling DataInitialize -- ", ESMF_LOGMSG_INFO)

    ! indicate that data initialization is complete (breaking out of init-loop)
    call NUOPC_CompAttributeSet(model, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

end module
