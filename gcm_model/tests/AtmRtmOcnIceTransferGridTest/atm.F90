!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2018, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================
!#define xUSE_ICE

module ATM

  !-----------------------------------------------------------------------------
  ! ATM Component.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
!  use NUOPC_Model, &
!    model_routine_SS      => SetServices, &
!    model_label_SetClock  => label_SetClock, &
!    model_label_Advance   => label_Advance
!==
  use NUOPC_Model, inheritModel => SetServices
!--
  
  implicit none
  
  private
  
  public SetServices

  TYPE(ESMF_VM) :: vm
  INTEGER :: lPet

  character(len=45)      :: importFieldList(6) = (/ &
    "Total Sky Shortwave Upward Flux             ", &
    "Total Sky Shortward Downard Flux            ", &
    "Total Sky Shortward Radiative Heating Rate  ", &
    "Clear Sky Shortwave Upward Flux             ", &
    "Clear Sky Shortwave Downward Flux           ", &
    "Clear Sky Shortwave Radiative Heating Rate  " /) 
  character(len=8)      :: importFieldSN(6) = (/ &
    "tssuf   ", &
    "tssdf   ", &
    "tssrhr  ", &
    "cssuf   ", &
    "cssdf   ", &
    "cssrhr  " /) 
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
!    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
!==
    call NUOPC_CompDerive(model, inheritModel, rc=rc)
!--
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
!===
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

!---

!    call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
!      specRoutine=SetClock, rc=rc)
!    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!      line=__LINE__, &
!      file=__FILE__)) &
!      return  ! bail out

!    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
!      specRoutine=ModelAdvance, rc=rc)
!==
    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
!--
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!==
    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    call ESMF_MethodRemove(model, label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_CheckImport, &
      specRoutine=CheckImport, rc=rc)
      !specRoutine=NUOPC_NOOP, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!--
    CALL ESMF_GridCompGet(model, vm=vm, rc=rc)
    CALL ESMF_VMGet(vm, localPet=lPet,rc=rc)

    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! importable field: sea_surface_temperature
    call NUOPC_Advertise(importState, &
      StandardName="sea_surface_temperature", name="sst", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!#ifdef USE_ICE
    ! importable field: sea_surface_salinity
    call NUOPC_Advertise(importState, &
      StandardName="sea_surface_salinity", name="sss", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!#endif   
 
    ! exportable field: air_pressure_at_sea_level
    call NUOPC_Advertise(exportState, &
      StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! exportable field: surface_net_downward_shortwave_flux
    call NUOPC_Advertise(exportState, &
      StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(importState, &
      standardNames=importFieldList, &
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
    type(ESMF_Grid)         :: gridIn
    type(ESMF_Grid)         :: gridOut
    integer                 :: i
    
    rc = ESMF_SUCCESS
    
    ! create a Grid object for Fields
    gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/20, 100/), &
      minCornerCoord=(/10._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
      maxCornerCoord=(/100._ESMF_KIND_R8, 200._ESMF_KIND_R8/), &
      coordSys=ESMF_COORDSYS_CART, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    gridOut = gridIn ! for now out same as in

    ! importable field: sea_surface_temperature
    field = ESMF_FieldCreate(name="sst", grid=gridIn, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!#ifdef USE_ICE
    ! importable field: sea_surface_salinity
    field = ESMF_FieldCreate(name="sss", grid=gridIn, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
!#endif

    ! exportable field: air_pressure_at_sea_level
    field = ESMF_FieldCreate(name="pmsl", grid=gridOut, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable field: surface_net_downward_shortwave_flux
    field = ESMF_FieldCreate(name="rsns", grid=gridOut, &
      typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i = 1, 6
      field = ESMF_FieldCreate(name=importFieldList(i), grid=gridOut, &
        typekind=ESMF_TYPEKIND_R8, rc=rc)
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

  subroutine xSetClock(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    !TODO: stabilityTimeStep should be read in from configuation
    !TODO: or computed from internal Grid information
!    call ESMF_TimeIntervalSet(stabilityTimeStep, m=5, rc=rc) ! 5 minute steps
    call ESMF_TimeIntervalSet(stabilityTimeStep, m=30, rc=rc) ! 30 minute steps
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetClock(model, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
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

  subroutine CheckImport1(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    ! local variables
    type(ESMF_Clock)                :: clock
    type(ESMF_Time)                 :: currTime, invalidTime
    type(ESMF_State)                :: importState
    logical                         :: timeCheck
    type(ESMF_Field),       pointer :: fieldList(:)
    integer                         :: i

    rc = ESMF_SUCCESS

    ! query the Component for its clock and importState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up invalid time (by convention)
    call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, &
      h=00, m=00, s=00, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Loop through all the field in the importState, and test whether they
    ! are at invalidTime (ignore them for now), or at currTime. Any other
    ! time coming in would flag an incompatibility.

    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i=1, size(fieldList)
      timeCheck = NUOPC_IsAtTime(fieldList(i), invalidTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (timeCheck) then
        ! The field is at invalidTime

        ! -> In a real application mark the field with a flag as invalid
        !    so the actual model code can act accordingly.

        ! Here for purpose of demonstration just log a message and continue on.

        call ESMF_LogWrite("ATM: detected import field at invalidTime", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        ! The field is NOT at invalidTime -> it must then be at currTime or it
        ! is incompatible!

        ! check that Fields in the importState show correct timestamp
        timeCheck = NUOPC_IsAtTime(fieldList(i), currTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        if (.not.timeCheck) then
          !TODO: introduce and use INCOMPATIBILITY return codes!!!!
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="NUOPC INCOMPATIBILITY DETECTED: "//&
            "Import Field not at current time", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        endif

      endif
    enddo

  end subroutine


  subroutine CheckImport(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    ! local variables
    type(ESMF_Clock)                :: clock
    type(ESMF_Time)                 :: currTime, invalidTime
    type(ESMF_State)                :: importState
    logical                         :: timeCheck
    type(ESMF_Field),       pointer :: fieldList(:)
    character(len=256)              :: fname
    integer                         :: i

    rc = ESMF_SUCCESS

    ! query the Component for its clock and importState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>currTime: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="--------------------------------> stopTime: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up invalid time (by convention)
    call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, &
      h=00, m=00, s=00, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Loop through all the field in the importState, and test whether they
    ! are at invalidTime (ignore them for now), or at currTime. Any other
    ! time coming in would flag an incompatibility.

    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i=1, size(fieldList)
      timeCheck = NUOPC_IsAtTime(fieldList(i), invalidTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (timeCheck) then
        ! The field is at invalidTime

        ! -> In a real application mark the field with a flag as invalid
        !    so the actual model code can act accordingly.

        ! Here for purpose of demonstration just log a message and continue on.

        call ESMF_LogWrite("RTM: detected import field at invalidTime", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        ! The field is NOT at invalidTime -> it must then be at currTime or it
        ! is incompatible!

        call ESMF_FieldGet(fieldList(i), name=fname, rc=rc)
        print *, 'Field Name = ', fname
        !call ESMF_TimePrint(currTime)

        ! check that Fields in the importState show correct timestamp
        timeCheck = NUOPC_IsAtTime(fieldList(i), currTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        if (.not.timeCheck) then
          !TODO: introduce and use INCOMPATIBILITY return codes!!!!
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="NUOPC INCOMPATIBILITY DETECTED: "//&
            "Import Field not at current time", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        endif

      endif
    enddo

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

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>SRC currTime: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockPrint(clock, options="stopTime", &
      preString="--------------------------------> SRC stopTime: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    if(lPet.eq.0) then
!      print *, 'atm driver timeStep: '
!      call ESMF_TimeIntervalPrint(timeStep,rc=rc)
!    endif

    call ESMF_ClockSet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!    if(lPet.eq.0) then
!      print *, 'atm model timeStep: '
!      call ESMF_TimeIntervalPrint(timeStep,rc=rc)
!    endif

    ! check and set the component clock against the driver clock
    call NUOPC_CompCheckSetClock(model, driverClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="NUOPC INCOMPATIBILITY DETECTED: between model and driver clocks", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out


  end subroutine SetRunClock

end module
